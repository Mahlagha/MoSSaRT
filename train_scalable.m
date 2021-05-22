function [timeElapsedSeconds, Theta, K, Phi] = ...
    train_scalable(X,rho,kernel,sigmaSquared,feature_dim, maxIter,eps_abs, eps_rel,l_norm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Optimizationa algorithm for the method proposed in the paper Sketches by MoSSaRT: Representative Selection from Manifolds with Gross Sparse Corruptions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<9
    l_norm=2;
end

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lambda = 1;
thr = 1*10^-6;

%% Compute Lambda and K
[computed_lambda, Phi,K] = computeLambda_Phi_K_scalable(X,kernel,sigmaSquared,feature_dim);
mu1 = rho* 1/computed_lambda; % This is in practice the coefficient for recosntruction error


N = size(K,1);
[~,p] = chol(K);
if p
    K = K + 1.1*abs(min(eig(K)))*eye(N);
end


% Start timer
startTime = cputime;

%% Initialization
Theta_init = zeros(N,N);
B_init = zeros(N,N);
Q_init = zeros(feature_dim,N);
U1_init = zeros(N,N);
U2_init = zeros(feature_dim,N);



Theta_prime = Theta_init;
B_prime = B_init;
Q_prime = Q_init;
U1_prime = U1_init;
U2_prime = U2_init;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Variable initialization
iter = 0;
notConverged = true;

while notConverged
    iter = iter+1;
    fprintf(1,'#%4d\titer=%d\n',0, iter);

    %% Update Q %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q = proxRegRow_l1(Phi-Phi*Theta_prime-U2_prime,1.0/computed_lambda);
    % the second argument was rho/computed_lambda*mu, since mu is set t
    % rho, so it gets simplified to 1.0/computed_lambda
    %% Update B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if l_norm==1
        B = proxRegRow_l1(Theta_prime-U1_prime,double(lambda)/rho);
    else
        B = proxRegRow_l1_l2(Theta_prime-U1_prime,double(lambda)/rho);
    end
    
    %% Update Theta %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Theta = (K+eye(N))\(B+ K -Phi'*(Q+U2_prime)+U1_prime);

    %% Update U %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U1 = U1_prime+B-Theta;
    U2 = U2_prime+Q-Phi+Phi*Theta;
    
    
    %% Check for convergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check if the algorithm has converged
    r_1_norm = norm(Q-Phi+Phi*Theta,'fro');
    r_2_norm = norm(B-Theta,'fro');
   
    s1_norm = norm(rho*(Phi*(Theta-Theta_prime)),'fro');
    s2_norm = norm(rho*(Theta_prime-Theta),'fro');
    
    eps_pri_1 = N*eps_abs + eps_rel*max(norm(Q,'fro'),norm(Phi-Phi*Theta,'fro'));
    eps_pri_2 = N*eps_abs + eps_rel*max(norm(B,'fro'),norm(Theta,'fro'));
    eps_dual_1 = N*eps_abs + eps_rel*rho* norm( U2, 'fro');
    eps_dual_2 = N*eps_abs + eps_rel*rho* norm( U1, 'fro');
    % Stopping Criteria for ADMM
    if r_1_norm <=eps_pri_1 && r_2_norm <=eps_pri_2 && s1_norm <=eps_dual_1 && s2_norm <=eps_dual_2
        notConverged = false;
        
        
        fprintf(1, 'CONVERGED!\n\n');
        
        break
    end
  

    % Check if the algorithm has exhausted the maximum number of iterations
    if iter >= maxIter && notConverged
        
        
        fprintf(1, 'STOPPED (exceeded maximum itearions)...\n\n');
        
        break
    end
    
    % if we did not need to halt the algorithm, update variables
    Theta_prime = Theta;
    Q_prime = Q;
    B_prime = B;
    U1_prime = U1;
    U2_prime = U2;
    
    
end

% Populate return arguments
timeElapsedSeconds = cputime - startTime;

end % main while-loop

function value = proxRegRow_l1(X,rho)
%value = (max((arrayfun(@abs,X)-rho),0)).*arrayfun(@sign,X);
value = max(X-rho,0) - max(-X-rho,0);
end

function value = proxRegRow_l1_l2(X,rho)
rNorm= rowNorm(X);
value = (max((1-((rho* ones(size(X)))./ repmat(rNorm,1,size(X,2)))),0)) .*X;
end

