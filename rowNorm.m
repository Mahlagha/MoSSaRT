
function rNorm = rowNorm(X)
temp = X .* X;
%vColumnNorms = sqrt(sum(temp, 1))';
rNorm = sqrt(sum(temp, 2));
end