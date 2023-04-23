function [c] = computeCoeff(mesh,pattern,L,m)
%   pattern may be = [0, 1, 2]
%   [-1, 0, 1] or [-2, -1, 0]
%   L is diff operator
%   for example [0, 1, 0]
%   m - current node
    A = ones(length(L),length(pattern));
    for i = 1:(length(L))
        for j = 1:(length(pattern))
            A(i,j) = ((mesh(m+pattern(j))-mesh(m))^(i-1))/factorial(i-1);
        end
    end
    c = linsolve(A,transpose(L));
end