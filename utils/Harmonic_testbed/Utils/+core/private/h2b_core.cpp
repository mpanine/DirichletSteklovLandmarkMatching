#include "mex.h"
#include <limits>
#include <iostream>
#include "mex.h"
#include <vector>
#include <cmath>
#include <string>


void hexChar2bin(char *hexChar, double *outBin) {
    
     switch(*hexChar){
         case '0': outBin[0] = 0;outBin[1] = 0;outBin[2] = 0;outBin[3] = 0;break;//0000
         case '1': outBin[0] = 0;outBin[1] = 0;outBin[2] = 0;outBin[3] = 1;break;//0001
         case '2': outBin[0] = 0;outBin[1] = 0;outBin[2] = 1;outBin[3] = 0;break;//0010
         case '3': outBin[0] = 0;outBin[1] = 0;outBin[2] = 1;outBin[3] = 1;break;//0011
         case '4': outBin[0] = 0;outBin[1] = 1;outBin[2] = 0;outBin[3] = 0;break;//0100
         case '5': outBin[0] = 0;outBin[1] = 1;outBin[2] = 0;outBin[3] = 1; break;//0101
         case '6': outBin[0] = 0;outBin[1] = 1;outBin[2] = 1;outBin[3] = 0; break;//0110
         case '7': outBin[0] = 0;outBin[1] = 1;outBin[2] = 1;outBin[3] = 1;break;//0111
         case '8': outBin[0] = 1;outBin[1] = 0;outBin[2] = 0;outBin[3] = 0;break;//1000
         case '9': outBin[0] = 1;outBin[1] = 0;outBin[2] = 0;outBin[3] = 1;break;//1001
         case 'A': outBin[0] = 1;outBin[1] = 0;outBin[2] = 1;outBin[3] = 0;break;//1010
         case 'B': outBin[0] = 1;outBin[1] = 0;outBin[2] = 1;outBin[3] = 1;break;//1011
         case 'C': outBin[0] = 1;outBin[1] = 1;outBin[2] = 0;outBin[3] = 0; break;//1100
         case 'D': outBin[0] = 1;outBin[1] = 1;outBin[2] = 0;outBin[3] = 1;break;//1101
         case 'E': outBin[0] = 1;outBin[1] = 1;outBin[2] = 1;outBin[3] = 0;break;//1110
         case 'F': outBin[0] = 1;outBin[1] = 1;outBin[2] = 1;outBin[3] = 1;break;//1111
         case 'a': outBin[0] = 1;outBin[1] = 0;outBin[2] = 1;outBin[3] = 0;break;//1010
         case 'b': outBin[0] = 1;outBin[1] = 0;outBin[2] = 1;outBin[3] = 1;break;//1011
         case 'c': outBin[0] = 1;outBin[1] = 1;outBin[2] = 0;outBin[3] = 0; break;//1100
         case 'd': outBin[0] = 1;outBin[1] = 1;outBin[2] = 0;outBin[3] = 1;break;//1101
         case 'e': outBin[0] = 1;outBin[1] = 1;outBin[2] = 1;outBin[3] = 0;break;//1110
         case 'f': outBin[0] = 1;outBin[1] = 1;outBin[2] = 1;outBin[3] = 1;break;//1111
     }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[])
{ 
    if (nrhs != 2)
		mexErrMsgTxt("Usage: [bit_array] = h2b(hex_string, n_bits).");
    if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments."); 

	// input
    const mxArray *hex_data = prhs[0];
    mwSize len = mxGetNumberOfElements(hex_data);
    std::string hex_str;
    hex_str.resize(len + 1);
    mxGetString(hex_data, &hex_str.front(), hex_str.size());
    hex_str.resize(len);
    int size = mxGetNumberOfElements(hex_data);

    int n_bits = mxGetScalar(prhs[1]);	

    // output
	plhs[0] = mxCreateDoubleMatrix( (mwSize)1, (mwSize)n_bits, mxREAL);
	double* outBin = mxGetPr(plhs[0]);
    
    std::string zeros_str;
    if(size > n_bits/4) {
        mexErrMsgTxt("The input of hexToBinaryVector has to be a hex string with n_bits/4 characters.");
    }
    else if(size < n_bits/4){
        zeros_str = std::string(n_bits/4-size, '0');
    }
    else
        zeros_str = std::string("");

    std::string full_hex_str;
    full_hex_str.append(zeros_str);
    full_hex_str.append(hex_str);

    double outBin_tmp[4];
    for(int i=0;i<n_bits/4;i++) {
        hexChar2bin(&full_hex_str[i], outBin_tmp);
        for(int j=0;j<4;j++) {
            outBin[4*i+j] = outBin_tmp[j];//4*(i-1)+1:4*(i-1)+4 4*i-3:4*i
        }
    }
    return;
}

