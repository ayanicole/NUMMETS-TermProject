function Tq = interp_field(T, grid, xq, yq, deg)
% INTERP_FIELD  Temperature at off-node points by Lagrange interpolation.
%
%   Tq = interp_field(T, grid, xq, yq)        % cubic (deg = 3)
%   Tq = interp_field(T, grid, xq, yq, deg)   % deg = 1..5
%
%   Tensor-product Lagrange polynomial through the (deg+1) x (deg+1) nodes
%   nearest the query point. Exact for polynomial fields up to degree deg.
%
% INPUT : T (Ny x Nx map or N x 1 vector), grid, xq/yq (same size, in meters), deg
% OUTPUT: Tq (same size as xq)

    if nargin < 5 || isempty(deg), deg = 3; end
    Nx = grid.Nx;  Ny = grid.Ny;
    if isvector(T) && numel(T) == grid.N
        Tm = reshape(T, Nx, Ny).';
    else
        Tm = T;
    end
    if deg + 1 > min(Nx, Ny)
        error('interp_field:deg', 'deg too large for this grid.');
    end
    if any(xq(:) < grid.x(1) | xq(:) > grid.x(end) | yq(:) < grid.y(1) | yq(:) > grid.y(end))
        error('interp_field:range', 'Query point outside the plate.');
    end

    Tq = zeros(size(xq));
    for m = 1:numel(xq)
        i0 = window(xq(m), grid.x(1), grid.dx, deg, Nx);
        j0 = window(yq(m), grid.y(1), grid.dy, deg, Ny);
        xn = grid.x(i0:i0+deg);
        yn = grid.y(j0:j0+deg);
        wx = lagrange_weights(xn, xq(m));            % 1 x (deg+1)
        wy = lagrange_weights(yn, yq(m));
        Tq(m) = wy * Tm(j0:j0+deg, i0:i0+deg) * wx.';
    end
end

function s = window(xq, x1, h, deg, n)
% First index of the (deg+1)-node window centered on the query point.
    fi = (xq - x1)/h + 1;                            % fractional 1-based index
    s  = round(fi - deg/2);
    s  = max(1, min(n - deg, s));
end

function w = lagrange_weights(xn, xq)
% w(a) = prod_{c ~= a} (xq - xn(c)) / (xn(a) - xn(c))
    m = numel(xn);
    w = ones(1, m);
    for a = 1:m
        for c = [1:a-1, a+1:m]
            w(a) = w(a) * (xq - xn(c)) / (xn(a) - xn(c));
        end
    end
end
