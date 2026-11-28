function [Q, qline] = heat_flow(T, grid, P, rule, edge, mode)
% HEAT_FLOW  Total heat flow across one edge by integrating the heat flux.
%
%   Default ('gradient'): flux from a second-order one-sided difference
%   (numerical differentiation) at the edge nodes, then integrated along the
%   edge with the composite Trapezoidal or Simpson's 1/3 rule (numerical
%   integration). Works on any edge.
%
%   mode = 'radiation': flux from the radiation law eps*sigma*(T^4 - Tamb^4)
%   at the edge nodes (use ONLY on radiating edges). With the 'trap' rule this
%   closes the discrete energy balance exactly: sum of radiating-edge Q = nodes.Psrc.
%
%   SIGN CONVENTION: Q > 0 means heat LEAVES the plate through that edge,
%   Q < 0 means heat ENTERS (outward-normal flux q'' = -k dT/dn).
%
% INPUT : T     temperature field, Ny x Nx map or N x 1 vector (see build_grid)
%         grid  struct from build_grid
%         P     parameters struct (uses P.k; eps, sigma, Tamb for 'radiation')
%         rule  'trap' | 'simpson'                              (default 'trap')
%         edge  'S' (y=0) | 'N' (y=Ly) | 'W' (x=0) | 'E' (x=Lx)   (default 'S')
%         mode  'gradient' | 'radiation'                       (default 'gradient')
% OUTPUT: Q     heat flow [W per unit depth, W/m]
%         qline outward flux along the edge [W/m^2], 1 x (nodes on that edge)
% Note: Simpson's rule needs an even number of intervals along the edge.

    if nargin < 4 || isempty(rule), rule = 'trap'; end
    if nargin < 5 || isempty(edge), edge = 'S'; end
    if nargin < 6 || isempty(mode), mode = 'gradient'; end

    Nx = grid.Nx;  Ny = grid.Ny;
    if isvector(T) && numel(T) == grid.N
        Tm = reshape(T, Nx, Ny).';
    else
        Tm = T;
    end
    if ~isequal(size(Tm), [Ny, Nx])
        error('heat_flow:size', 'T must be %d x %d or %d x 1.', Ny, Nx, grid.N);
    end

    k = P.k;
    switch upper(edge)
        case 'S'      % y = 0, outward normal -y  ->  q = +k dT/dy
            Te    = Tm(1, :);
            qline = k*(-3*Tm(1, :) + 4*Tm(2, :) - Tm(3, :)) / (2*grid.dy);
            h = grid.dx;
        case 'N'      % y = Ly, outward normal +y ->  q = -k dT/dy
            Te    = Tm(Ny, :);
            qline = -k*(3*Tm(Ny, :) - 4*Tm(Ny-1, :) + Tm(Ny-2, :)) / (2*grid.dy);
            h = grid.dx;
        case 'W'      % x = 0, outward normal -x  ->  q = +k dT/dx
            Te    = Tm(:, 1).';
            qline = k*(-3*Tm(:, 1) + 4*Tm(:, 2) - Tm(:, 3)).' / (2*grid.dx);
            h = grid.dy;
        case 'E'      % x = Lx, outward normal +x ->  q = -k dT/dx
            Te    = Tm(:, Nx).';
            qline = -k*(3*Tm(:, Nx) - 4*Tm(:, Nx-1) + Tm(:, Nx-2)).' / (2*grid.dx);
            h = grid.dy;
        otherwise
            error('heat_flow:edge', 'edge must be ''S'', ''N'', ''W'' or ''E''.');
    end

    switch lower(mode)
        case 'gradient'
            % qline already holds the conduction flux
        case 'radiation'
            qline = P.eps*P.sigma*(Te.^4 - P.Tamb^4);
        otherwise
            error('heat_flow:mode', 'mode must be ''gradient'' or ''radiation''.');
    end

    m = numel(qline);                     % nodes on the edge
    switch lower(rule)
        case 'trap'
            Q = h*(sum(qline) - 0.5*(qline(1) + qline(end)));
        case 'simpson'
            if mod(m - 1, 2) ~= 0
                error('heat_flow:simpson', ...
                    'Simpson''s rule needs an even number of intervals (got %d).', m - 1);
            end
            Q = h/3*(qline(1) + qline(end) ...
                     + 4*sum(qline(2:2:m-1)) + 2*sum(qline(3:2:m-2)));
        otherwise
            error('heat_flow:rule', 'rule must be ''trap'' or ''simpson''.');
    end
end
