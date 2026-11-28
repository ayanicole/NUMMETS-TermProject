function [T, hist] = solve_iterative(A, b, T0, method, tol, maxIter, omega)
% Solve A*T = b using Jacobi, Gauss-Seidel, or SOR.
% hist.res = residual per iteration
% hist.dT  = change between iterations

N = length(b);
T = T0(:);
b = b(:);

hist.res = zeros(maxIter, 1);
hist.dT = zeros(maxIter, 1);

if any(diag(A) == 0)
    error('A contains a zero diagonal entry.');
end

if ~any(strcmpi(method, {'jacobi', 'gs', 'sor'}))
    error('Unknown method.');
end

for k = 1:maxIter
    Told = T;

    for i = 1:N
        if strcmpi(method, 'jacobi')
            % Use only values from the previous iteration
            sigma = A(i,:) * Told - A(i,i) * Told(i);
            T(i) = (b(i) - sigma) / A(i,i);
        else
            % Use the latest available values
            sigma = A(i,:) * T - A(i,i) * T(i);
            Tgs = (b(i) - sigma) / A(i,i);

            if strcmpi(method, 'sor')
                T(i) = Told(i) + omega * (Tgs - Told(i));
            else
                T(i) = Tgs;
            end
        end
    end

    hist.res(k) = norm(A*T - b, inf);
    hist.dT(k) = norm(T - Told, inf);

    if hist.dT(k) < tol && hist.res(k) < tol
        break;
    end
end

hist.res = hist.res(1:k);
hist.dT = hist.dT(1:k);

if k == maxIter && (hist.dT(k) >= tol || hist.res(k) >= tol)
    warning('Solver reached maxIter before convergence.');
end

end
