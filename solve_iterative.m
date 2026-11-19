function [T, hist] = solve_iterative(A, b, T0, method, tol, maxIter, omega)
% Solve A*T = b iteratively.
% INPUT : A, b, T0 (initial guess), method ('jacobi' | 'gs' | 'sor'),
%         tol (stopping tolerance), maxIter, omega (SOR factor)
% OUTPUT: T (solution), hist (struct: hist.res = residual per iteration,
%         hist.dT = change between iterations)
    error('Not implemented');
end
