function q = make_source(grid, kind, Ptot, t, varargin)
% MAKE_SOURCE  Heat-generation field q [W/m^3] with an EXACT total power.
%
%   q = make_source(grid, 'uniform', Ptot, t)
%   q = make_source(grid, 'rect',  Ptot, t, xc, yc, w)        % square, half-width w
%   q = make_source(grid, 'rect',  Ptot, t, xc, yc, wx, wy)   % rectangle
%   q = make_source(grid, 'gauss', Ptot, t, xc, yc, s)        % Gaussian, std dev s
%
%   Ptot : total power [W] generated over the plate thickness t [m]
%   Use:  P.q = make_source(...);   assemble_system then gives  nodes.Psrc*t == Ptot.
%
% WHY: a sharp-edged source sampled on nodes changes its area with roundoff at
% the edges (we measured 10 W requested -> 12.7 W generated). The shape is
% therefore normalized with the same trapezoid weights assemble_system uses, so
% "hotspot vs uniform" comparisons are at exactly equal power.

    Nx = grid.Nx;  Ny = grid.Ny;  dx = grid.dx;  dy = grid.dy;

    switch lower(kind)
        case 'uniform'
            S = ones(Ny, Nx);
        case 'rect'
            xc = varargin{1};  yc = varargin{2};  wx = varargin{3};
            if numel(varargin) >= 4, wy = varargin{4}; else, wy = wx; end
            tol = 1e-9*min(dx, dy);                       % symmetric inclusion at the edges
            S = double(abs(grid.X - xc) <= wx + tol & abs(grid.Y - yc) <= wy + tol);
            if nnz(S) < 4
                warning('make_source:coarse', 'Source covers only %d nodes - refine the grid.', nnz(S));
            end
        case 'gauss'
            xc = varargin{1};  yc = varargin{2};  s = varargin{3};
            S = exp(-((grid.X - xc).^2 + (grid.Y - yc).^2)/(2*s^2));
        otherwise
            error('make_source:kind', 'kind must be ''uniform'', ''rect'' or ''gauss''.');
    end

    wxv = dx*ones(1, Nx);  wxv([1 Nx]) = dx/2;            % trapezoid weights = control volumes
    wyv = dy*ones(Ny, 1);  wyv([1 Ny]) = dy/2;
    V = wyv*wxv;                                          % Ny x Nx
    A = sum(S(:).*V(:));
    if A <= 0
        error('make_source:empty', 'Source contains no nodes.');
    end
    q = (Ptot/t)*S/A;
end
