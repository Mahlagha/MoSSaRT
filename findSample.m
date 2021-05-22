function [idx, sub_idx] = findSample(Kernel,p, S,r)
% r is the 1-threshold for removing similarity - if 0 means not removing anything 
if (nargin < 4)
    r=true;
end

% fin nonzero rows of C by Integral criteria, sigma(p_i)/sigma(p)>0.95

idx_nonZero  = find_nonZero_rows(p);

if r
    
    %% This removes based on Hilbert space distance
    % remove the similar ones
    if length(idx_nonZero)>1
        idx = rmv_sim_Hilbert(idx_nonZero, Kernel,p, 1-r );
    end
else
        idx = idx_nonZero;
        
end
% now idx is sorted according to p values 'descend'
if S<length(idx)
    sub_idx = idx(1:S);
else
    sub_idx = idx;
end

    
        
