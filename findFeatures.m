function Phi = findFeatures(X,kernel, sigmaSquared, feature_dim)

    N =size(X,2);
    d = size(X,1);
    if strcmp('gaussianKernelMatrix',func2str(kernel))
        % find Fourier Transform of the Gaussian kernel
        % feature_dim iid samples from multivariate normal with Mu=0, Sigma = (2/(sigma^2))I
        Mu = zeros(1,d);
        Sigma = (2/sigmaSquared)*ones(1,d);
        W = mvnrnd(Mu,Sigma,feature_dim);
    end
    
    b = 2*pi*rand(feature_dim,1);
    
    Phi = sqrt(2/feature_dim)*cos(W*X+b);
end

    
    
