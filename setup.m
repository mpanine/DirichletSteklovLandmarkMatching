clear all;
close all;

addpath(genpath('./utils'));
%%
mex -setup C++;
%%
run('./utils/Harmonic_testbed/ann_mwrapper/ann_compile_mex.m');
%%
cd('./utils/Harmonic_testbed/Utils/+core/private');
%
mex -O -D_LITTLE_ENDIAN -DHAS_HG2 -largeArrayDims -DMATLABVER=901 b2h_core.cpp -v
%
mex -O -D_LITTLE_ENDIAN -DHAS_HG2 -largeArrayDims -DMATLABVER=901 h2b_core.cpp -v
%%
mex -O -D_LITTLE_ENDIAN -DHAS_HG2 -largeArrayDims -DMATLABVER=901 GetMD5.cpp -v
%%
cd('../../../../../');
movefile('./utils/Harmonic_testbed/Utils/+core/private/GetMD5.mex*','./utils/Harmonic_testbed/Utils/+core/');
%%
cd('./utils/Harmonic_testbed/geodesics/');
mex -output fastmarchmex my_heap.cpp unfold.cpp;
%%
cd('../../../');