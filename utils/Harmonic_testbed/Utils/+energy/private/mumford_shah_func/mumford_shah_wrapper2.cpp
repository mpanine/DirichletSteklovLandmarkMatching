/**
 * This code accompanies the paper:
 *
 * "Partial Functional Correspondence"
 * Rodola, Cosmo, Bronstein, Torsello, Cremers
 * Computer Graphics Forum 2016
 *
 * Please cite the paper above if you use this code in your research.
 *
 * Written by Luca Cosmo
 * Mar 2015
 *
 * To compile:
 * mex mumford_shah_wrapper2.cpp COMPFLAGS="/openmp $COMPFLAGS"
 */
#include "mex.h"
#include <cmath>
#include <vector>
#include <queue>
#include <iostream>
#include "mesh.h"
#include <climits>

#include "class_handle.hpp"

#define ALMOST_ZERO 1e-10


void initialize(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    if (nrhs != 2 || nlhs > 1)
        mexErrMsgTxt("Usage: [mesh] = init_mesh(vertices, triangles).");
    
    const double* const pts = mxGetPr(prhs[0]);
    const int np = int( mxGetN(prhs[0]) );
    const double* const tri = mxGetPr(prhs[1]);
    const int nt = int( mxGetN(prhs[1]) );
    
    if (np == 3)
        mexErrMsgTxt("It seems like you only have 3 vertices. Please try to transpose the input matrix.");
    
    if (nt == 3)
        mexErrMsgTxt("It seems like you only have 3 triangles. Please try to transpose the input matrix.");
    
    // Load the mesh
    mesh_t* mesh = new mesh_t();
    
    std::vector< vec3d<double> > vertices(np);
    
    for (int i=0; i<np; ++i)
    {
        vec3d<double>& pt = vertices[i];
        pt.x = pts[i*3];
        pt.y = pts[i*3+1];
        pt.z = pts[i*3+2];
        //std::cout << pt.x << " " << pt.y << " " << pt.z << std::endl;
    }
    
    mesh->put_vertices(vertices);
    
    for (int i=0; i<nt; ++i)
    {
        int a, b, c;
        a = (int)tri[i*3] - 1; // 1-based to 0-based
        b = (int)tri[i*3+1] - 1;
        c = (int)tri[i*3+2] - 1;
        mesh->add_triangle(a,b,c);
        //std::cout << a << " " << b << " " << c << std::endl;
    }
    
    int dims[2]; dims[0] = 1; dims[1] = 1;
    plhs[0] = convertPtr2Mat<mesh_t>(mesh);
    
    for (int j=0; j<nt; ++j)
    {
        mesh_t::triangle_data& tri = mesh->triangles[j];
        
        const vec3d<double>& xj1 = mesh->vertices[tri.p0].p;
        const vec3d<double>& xj2 = mesh->vertices[tri.p1].p;
        const vec3d<double>& xj3 = mesh->vertices[tri.p2].p;
        
        tri.E = dot_product(xj2-xj1, xj2-xj1);
        tri.F = dot_product(xj2-xj1, xj3-xj1);
        tri.G = dot_product(xj3-xj1, xj3-xj1);
    }
//     mexLock();
//     mexPrintf("Initialized mesh %d -\n",mesh);
}

void compute_cost(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    if (nrhs != 6 || nlhs > 1)
        mexErrMsgTxt("Usage: [cost_per_tri] = calc_cost_mumford_shah(mesh_ptr, scalar_fun, sigma, area_weighted, mean, perturb).");
    
    const mesh_t* mesh = convertMat2Ptr<mesh_t>(prhs[0]);
    
    const double* const v = mxGetPr(prhs[1]);
    const int nvert = int( mxGetM(prhs[1]) );
    
    if (int( mxGetN(prhs[1]) ) != 1 && int( mxGetM(prhs[1]) ) != 1)
        mexErrMsgTxt("The scalar function must be either a row or a column vector.");
    
    const double sigma = *mxGetPr(prhs[2]);
    const bool area_weighted = *mxGetPr(prhs[3]);
    const double mean = *mxGetPr(prhs[4]);
    
    const double* const perturb = mxGetPr(prhs[5]);
    
    int nt = mesh->triangles.size();
    int np = mesh->vertices.size();
    
    if (nvert!=np)
        mexErrMsgTxt("COST: wrong number of vertices");
    
    plhs[0] = mxCreateDoubleMatrix(nt, 1, mxREAL);
    double* cost_ms = mxGetPr(plhs[0]);
    
    const double var = 2*sigma*sigma;
    
    double* h = new double[np];
    
#pragma omp parallel for
    for (int i=0; i<np; ++i)
    {
        h[i] = std::exp( - (v[i]-mean)*(v[i]-mean) / var );
    }
    
#pragma omp parallel for
    for (int j=0; j<nt; ++j)
    {
        const mesh_t::triangle_data& tri = mesh->triangles[j];
        
        const double E = tri.E;
        const double F = tri.F;
        const double G = tri.G;
        
        const double v1 = v[tri.p0];
        const double v2 = v[tri.p1];
        const double v3 = v[tri.p2];
        
        const double va = v2 - v1;
        const double vb = v3 - v1;
        
        double norm_grad_v = std::sqrt(va*va*G + vb*vb*E - 2*va*vb*F);
        
        const double hsum = (1.0/6.0)*(h[tri.p0] + h[tri.p1] + h[tri.p2]);
        
// 		if (hsum > ALMOST_ZERO)
//         cost_ms[j] = hsum * norm_grad_v;
        cost_ms[j] = hsum;
    }
    
    delete[] h;
}

void compute_gradient(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 6 || nlhs > 2)
        mexErrMsgTxt("Usage: [grad_ms] = calc_grad_mumford_shah(mesh_ptr, v, Phi=Phi'*diag(S), Psi , sigma, mean).");
    
    const mesh_t* mesh = convertMat2Ptr<mesh_t>(prhs[0]);
    
    const double* v = mxGetPr(prhs[1]);
    const double* PhitS = mxGetPr(prhs[2]);
    const double* Psi = mxGetPr(prhs[3]);
    
    const int k = mxGetM(prhs[2]);
    const int n = mxGetM(prhs[3]);
    
//     if (int( mxGetM(prhs[2]) ) != 1 && int( mxGetN(prhs[3]) ) != 1)
// 		mexErrMsgTxt("something bad!");
    
    const double sigma = *mxGetPr(prhs[4]);
    const double mean = *mxGetPr(prhs[5]);
    
    int nt = mesh->triangles.size();
    int np = mesh->vertices.size();
    
    if (n!=np)
        mexErrMsgTxt("wohoh!");
    
//     mexPrintf("k: %d, n: %d, nt: %d,  sigma=%f, mean=%f\n",k,n,nt,sigma,mean);
    
    
    plhs[0] = mxCreateDoubleMatrix( (mwSize)k, (mwSize)1, mxREAL);
    double* dC = mxGetPr(plhs[0]);
    
//     plhs[1] = mxCreateDoubleMatrix(nt, 1, mxREAL);
// 	double* grad_norm = mxGetPr(plhs[1]);
    
    const double var = sigma*sigma;
    
    double* h = new double[np];
    
#pragma omp parallel for
    for (int i=0; i<np; ++i)
        h[i] = std::exp( - ((v[i]-mean)*(v[i]-mean)) / (2*var) );
    
    double* A = new double[nt];
    double* Q = new double[nt];
    
#pragma omp parallel for
    for(int j=0; j<nt; j++)
    {
        const mesh_t::triangle_data& tri = mesh->triangles[j];
        
        const double v0 = v[tri.p0];
        const double v1 = v[tri.p1];
        const double v2 = v[tri.p2];
        
        const double va = v1 - v0;
        const double vb = v2 - v0;
        
        A[j] = (1.0/6.0)*(h[tri.p0]+h[tri.p1]+h[tri.p2]);
        Q[j] = sqrt(va*va*tri.G + vb*vb*tri.E - 2*va*vb*tri.F);
    }
    
    
//     for (int p=0; p<k; ++p)
//         for (int q=0; q<k; ++q)
//         {
//             dC[p + q*k]=0;
//         }
    
    {
//
#pragma omp parallel for
        for (int p=0; p<k; ++p)
        {
            double sumj=0;
            for(int j=0; j<nt; j++)
            {
                const mesh_t::triangle_data& tri = mesh->triangles[j];
                
                const double E = tri.E;
                const double F = tri.F;
                const double G = tri.G;
                
                const double v0 = v[tri.p0];
                const double v1 = v[tri.p1];
                const double v2 = v[tri.p2];
                
                const int p0 = tri.p0;
                const int p1 = tri.p1;
                const int p2 = tri.p2;
                
                const double va = v1 - v0;
                const double vb = v2 - v0;
                
                const double dA = -(1.0/(6.0*var))*(
                        (v0-mean)*Psi[p0 + n*p]*h[p0] +
                        (v1-mean)*Psi[p1 + n*p]*h[p1] +
                        (v2-mean)*Psi[p2 + n*p]*h[p2]);
                
                const double dQ = 2*G*va*(Psi[p1 + n*p]-Psi[p0 + n*p]) +
                        2*E*vb*(Psi[p2 + n*p]-Psi[p0 + n*p]) -
                        2*F*(vb*(Psi[p1 + n*p]-Psi[p0 + n*p]) + va*(Psi[p2 + n*p]-Psi[p0 + n*p]));
                
                const double dB = (0.5/Q[j])*dQ;
//                 sumj += (dA*Q[j]+dB*A[j]);
                sumj += (dA);
            }
            dC[p] = sumj;
        }
    }
    
    delete[] h;
    delete[] A;
    delete[] Q;
}

void destroy(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    if (nrhs != 1 || nlhs > 0)
        mexErrMsgTxt("Usage: destroy_mesh(mesh_ptr).");
    
    destroyObject<mesh_t>(prhs[0]);
//     mexUnlock();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{
    if (nrhs < 1)
        mexErrMsgTxt("Usage: init_mesh(action, ...) where actions={0:initialize, -1:destroy, 1:compute_cost, 2:compute_gradeint} .");
    
    if(mxGetClassID(prhs[0]) != mxDOUBLE_CLASS)
        mexErrMsgTxt("action should be a double value (double)");
    int action = *mxGetPr(prhs[0]);
    
    switch(action)
    {
        case 0: initialize( nlhs, plhs, nrhs-1, prhs+1); break;
        case -1: destroy( nlhs, plhs, nrhs-1, prhs+1); break;
        case 1: compute_cost( nlhs, plhs, nrhs-1, prhs+1); break;
        case 2: compute_gradient( nlhs, plhs, nrhs-1, prhs+1); break;
        default:
            mexPrintf("--- %d\n",action);
            mexErrMsgTxt("checazzo: init_mesh(action, ...) where actions={0:initialize, -1:destroy, 1:compute_cost, 2:compute_gradeint} .");
    }
}
