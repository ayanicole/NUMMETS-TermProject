function T = solve_gauss_jordan(A, b)
% Solve A*T = b by Gauss-Jordan elimination.
% Intended for small grids only.

N = length(b);
M = [full(A), b(:)];

for i = 1:N

    % Find the largest available pivot
    [~, p] = max(abs(M(i:N, i)));
    p = p + i - 1;

    if M(p, i) == 0
        error('Matrix is singular.');
    end

    % Swap rows
    M([i p], :) = M([p i], :);

    % Normalize pivot row
    M(i, :) = M(i, :) / M(i, i);

    % Eliminate other rows
    for j = 1:N
        if j ~= i
            factor = M(j, i);
            M(j, :) = M(j, :) - factor * M(i, :);
        end
    end
end

T = M(:, end);

end
