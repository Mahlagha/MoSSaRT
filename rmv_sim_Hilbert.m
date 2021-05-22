%--------------------------------------------------------------------------
% This function takes the data matrix and the indices of the
% representatives and removes the representatives that are too close to
% each other in the Hilbert space
% Kernel: N*N Kernel
% index: indices of the representatives
% thr: threshold for pruning the representatives, typically in [0.9,0.99]


function index = rmv_sim_Hilbert(init_idx,Kernel,p, thr)

if (nargin < 4)
    thr = 0.95;
end
% in case the init_idx is not sorted! 
[~,pSorted_idx] = sort(p(init_idx),'descend');
init_idx_sorted = init_idx(pSorted_idx);


Ns = length(init_idx);
d = zeros(Ns);
for i =1:Ns
    for j = 1:Ns
        d(i,j) = sqrt(Kernel(init_idx_sorted(i),init_idx_sorted(i)) - ...
            2*Kernel(init_idx_sorted(i),init_idx_sorted(j)) + ...
            Kernel(init_idx_sorted(j),init_idx_sorted(j)));
    end
end

[dSorted,dSorted_idx] = sort(d,'descend');

idx = 1:Ns;
for i = 1:Ns
    if (ismember(i,idx));
        summ = 0;
        t = 0;
        while ( (summ/sum(dSorted(:,i))) <= thr )
            t = t + 1;
            summ = summ + dSorted(t,i);
        end
        idx = setdiff(idx,setdiff(dSorted_idx(t:end,i),1:i));
    end
end
index = init_idx_sorted(idx);