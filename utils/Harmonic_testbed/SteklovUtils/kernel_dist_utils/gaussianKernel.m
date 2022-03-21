function GK = gaussianKernel(E1, E2, sigma)
%GAUSSIANKERNEL computes the Gaussian kernel of two embeddings
%E1: embedding 1
%E2: embedding 2 -> must already be deformed (if needed)
%sigma: "standard deviation" of gaussian kernel (kernel is not normalized as probability)


GK = exp( - pdist2(E1, E2, 'squaredeuclidean') /(sigma^2) );



end

