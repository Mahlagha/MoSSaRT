% This functiontakes the row-norms of the matrix C, and retruns the
% "non-zero" indices, w.r.t. integral criteria, i.e. 
% sigma(p_i)/sigma(p) >thr

function idx_nonZero = find_nonZero_rows(rNorm,thr)
if (nargin < 2)
    thr = 0.99;
end

[rNorm,idx] = sort(rNorm,'descend');
summ = 0;
i = 0; 
while (summ / sum(rNorm) <= thr)
    i = i+1;
    summ = summ + rNorm(i);
end
idx_nonZero = idx(1:i);