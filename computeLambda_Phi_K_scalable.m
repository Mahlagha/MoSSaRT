%--------------------------------------------------------------------------
% This function takes a DxN matrix of N data points in a D-dimensional 
% space and returns the regularization constant of the Gaussain Kernel
function [lambda, Phi,K] = computeLambda_Phi_K_scalable(X,kernel,sigmaSquared,feature_dim)

[~,N] = size(X);

X = X - repmat(mean(X,2),1,N);
Phi = findFeatures(X,kernel, sigmaSquared, feature_dim);
K = Phi'*Phi;
lambda = max(rowNorm(K));
    
