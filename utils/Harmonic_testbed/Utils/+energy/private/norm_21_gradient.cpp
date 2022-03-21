/**
 * This code accompanies the paper:
 *
 * "Partial Functional Correspondence"
 * Rodola, Cosmo, Bronstein, Torsello, Cremers
 * Computer Graphics Forum 2016
 *
 * Please cite the paper above if you use this code in your research.
 *
 * Written by Emanuele Rodola and Luca Cosmo
 * Mar 2015
 *
 * To compile:
 * mex norm_21_gradient.cpp COMPFLAGS="/openmp $COMPFLAGS"
 */  
#include <limits>
#include <iostream>
#include "mex.h"
#include <vector>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{ 
    if (nrhs != 3)
		mexErrMsgTxt("Usage: [S] = norm_21_gradient(C,A,B) where ||CA-B||_{2,1}.");
    if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments."); 

	// input
	const double *C = (double*)mxGetPr(prhs[0]);	
    const double *A = (double*)mxGetPr(prhs[1]);	
    const double *B = (double*)mxGetPr(prhs[2]);	
    
    const int m = int( mxGetM(prhs[1])); // number of harmonics
    const int n = int( mxGetN(prhs[1])); // number of probe functions

    // output
	plhs[0] = mxCreateDoubleMatrix( (mwSize)m, (mwSize)m, mxREAL);
	double* gC = mxGetPr(plhs[0]);
  
    std::vector<double> sums_i(n);
    std::vector<double> sums_ij(m*n);
    #pragma omp parallel for
    for(int j=0; j<n; ++j)
    {
        double sum_i = 0;
        for(int i=0; i<m; ++i)
        {
            double sum_r = 0;
            for(int r=0; r<m; ++r)
                sum_r += C[i + r*m]*A[r + j*m];
            sum_r -= B[i + j*m];
            sums_ij[i+m*j] = sum_r;
            sum_i += sum_r*sum_r;
        }
        if (sum_i <= 1e-10)
			sums_i[j] = 0.;
		else
			sums_i[j] = 1/std::sqrt(sum_i);
    }
	
	// run	
    #pragma omp parallel for
    for(int p=0; p<m; ++p)
    for(int q=0; q<m; ++q)
    {             
        double sum_j=0;
        for(int j=0; j<n; ++j)
        {
            sum_j += A[q + j*m]*sums_i[j]*sums_ij[p+m*j];
        }                            
         gC[p + q*m] = sum_j;
    }
    
    //delete[] sums1;
    return;
}
