%% Parameter Setup
add_corruption = true;
feature_dim=20;
l_norm = 2;
rem_sim = 0.05;
rho = 10;
maxIter = 5000;
eps_abs = 1e-5;
eps_rel = 1e-5;

function [time, Theta, idx_nonZero, idx_sample] = Mossart_sample_demo (data, num_sample, add_corruption, feature_dim, l_norm, rem_sim, rho, maxIter, eps_abs, eps_rel)

%%%%%%%%%%%%%%%%%%
%% This function runs a demo for sampling representatives from a data matrix
%% Input arguments:
% data: N d-dimensional data points in the form of a matrix \in R^{d*N}
% num_sample: Upperbound on number of representatives to be selected (If you don't have a value for this, just set a cap.)
% add_corruption: If true, adds sparse corruption matrix to the data.
% feature_dim: Dimension of the features in the Explicit transformation via Fourier features
% l_norm: Regularization norm for the row-sparsity. l_norm \in {1,2}
% rho: Reconstruction Coefficient
% maxIter, eps_abs, eps_rel: Optimization parameteres
%
%% Output:
% time: Running time of the algorithm
% Theta: The obtained Reproduction profile Theta
% idx_nonZero: Indices for the chosen represntatives samples
% idx_sample: Indices for the chosen representative samples capped by the num_sample.
%%%%%%%%%%%%%%%%%%%

N = size(data,2);
n1=N;
d = size(data,1);


% Adding Sparse corruptions
S = zeros(size(data))
if add_corruption:
    maxim = 5;
    S = maxim*rand(size(data)).*data.*double(rand(size(data))<0.05);
D = data+S;

sigma = (mean(var(D,0,2)))^2;
sigmaSquared = [sigma*0.1,sigma,sigma*10];

for k=1: length(sigmaSquared)
 
        %%  Gaussian Kernel
        [~, time, Theta, Kernel,~] = train_scalable( D,rho,...
            @gaussianKernelMatrix, sigmaSquared(k),...
            feature_dim,maxIter,eps_abs, eps_rel,l_norm);
        
        p = rowNorm(C);
        [idx_nonZero, idx_sample] = findSample(Kernel,p, num_sample, rem_sim);
    end
end


