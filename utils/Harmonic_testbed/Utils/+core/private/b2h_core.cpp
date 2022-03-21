#include "mex.h"
#include <limits>
#include <iostream>
#include "mex.h"
#include <vector>
#include <cmath>
#include <string>
#include <sstream>


void bin2hexChar(double bit1, double bit2, double bit3, double bit4, char *hexChar) {
    if (bit1==0 && bit2==0 && bit3==0 && bit4==0) {*hexChar = '0';}
    else if (bit1==0 && bit2==0 && bit3==0 && bit4==1) {*hexChar = '1';}
    else if (bit1==0 && bit2==0 && bit3==1 && bit4==0) {*hexChar = '2';}
    else if (bit1==0 && bit2==0 && bit3==1 && bit4==1) {*hexChar = '3';}
    else if (bit1==0 && bit2==1 && bit3==0 && bit4==0) {*hexChar = '4';}
    else if (bit1==0 && bit2==1 && bit3==0 && bit4==1) {*hexChar = '5';}
    else if (bit1==0 && bit2==1 && bit3==1 && bit4==0) {*hexChar = '6';}
    else if (bit1==0 && bit2==1 && bit3==1 && bit4==1) {*hexChar = '7';}
    else if (bit1==1 && bit2==0 && bit3==0 && bit4==0) {*hexChar = '8';}
    else if (bit1==1 && bit2==0 && bit3==0 && bit4==1) {*hexChar = '9';}
    else if (bit1==1 && bit2==0 && bit3==1 && bit4==0) {*hexChar = 'A';}
    else if (bit1==1 && bit2==0 && bit3==1 && bit4==1) {*hexChar = 'B';}
    else if (bit1==1 && bit2==1 && bit3==0 && bit4==0) {*hexChar = 'C';}
    else if (bit1==1 && bit2==1 && bit3==0 && bit4==1) {*hexChar = 'D';}
    else if (bit1==1 && bit2==1 && bit3==1 && bit4==0) {*hexChar = 'E';}
    else if (bit1==1 && bit2==1 && bit3==1 && bit4==1) {*hexChar = 'F';}
//     if (!bin_str.compare("0000"))
//     {
//         *hexChar = '0';
//     }
//     else if (!bin_str.compare("0001"))
//     {
//         *hexChar = '1';
//     }
//     else if (!bin_str.compare("0010"))
//     {
//         *hexChar = '2';
//     }
//     else if (!bin_str.compare("0011"))
//     {
//         *hexChar = '3';
//     }
//     else if (!bin_str.compare("0100"))
//     {
//         *hexChar = '4';
//     }
//     else if (!bin_str.compare("0101"))
//     {
//         *hexChar = '5';
//     }
//     else if (!bin_str.compare("0110"))
//     {
//         *hexChar = '6';
//     }
//     else if (!bin_str.compare("0111"))
//     {
//         *hexChar = '7';
//     }
//     else if (!bin_str.compare("1000"))
//     {
//         *hexChar = '8';
//     }
//     else if (!bin_str.compare("1001"))
//     {
//         *hexChar = '9';
//     }
//     else if (!bin_str.compare("1010"))
//     {
//         *hexChar = 'A';
//     }
//     else if (!bin_str.compare("1011"))
//     {
//         *hexChar = 'B';
//     }
//     else if (!bin_str.compare("1100"))
//     {
//         *hexChar = 'C';
//     }
//     else if (!bin_str.compare("1101"))
//     {
//         *hexChar = 'D';
//     }
//     else if (!bin_str.compare("1110"))
//     {
//         *hexChar = 'E';
//     }
//     else if (!bin_str.compare("1111"))
//     {
//         *hexChar = 'F';
//     }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{ 
    if (nrhs > 1)
		mexErrMsgTxt("Usage: [hex_string] = h2b(bin_vec).");
    if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments."); 

	// input
	double *bin_vec = (double*)mxGetData(prhs[0]);
    
    int size = int( mxGetNumberOfElements(prhs[0]));
    
    if(size%4 != 0)
        mexErrMsgTxt("The input of binaryVectorToHex has to be a multiple of 4.");
    
    int n_char = size/4;

    // output
// 	output_buf = mxCalloc(n_char, sizeof(char));
//     plhs[0] = mxCreateString(output_buf);
    
    std::string out_hex_str(n_char, '0');
    
//     std::ostringstream strs;
    for(int i=0;i<n_char;i++) {
//         strs << bin_vec[i*4] << bin_vec[i*4+1] << bin_vec[i*4+2] << bin_vec[i*4+3];
        char tmp_char = '0';
//         bin2hexChar(strs.str(),&tmp_char);
        bin2hexChar(bin_vec[i*4], bin_vec[i*4+1], bin_vec[i*4+2], bin_vec[i*4+3], &tmp_char);
        out_hex_str[i] = tmp_char;
//         strs.clear();
    }
    
    plhs[0] = mxCreateString(out_hex_str.c_str());
    return;
}

