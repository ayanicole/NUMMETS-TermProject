function [A, b, nodes] = assemble_system(P, grid)
% ASSEMBLE_SYSTEM  5-point finite-difference system for steady 2D conduction.
%
%   Governing equation (per unit depth, constant k):
%       k*(d2T/dx2 + d2T/dy2) + q = 0        (q = heat generation, W/m^3)
%
%   Every non-Dirichlet node uses an energy balance on its control volume
%   (a half cell on an edge, a quarter cell at a corner), divided by k*Volume:
%
%       diag*T_p - sum_nb a_nb*T_nb  +  rho_p*(T_p^4 - Tamb^4)  =  q_p/k + hl*Tamb
%
%   with  a_E = a_W = 1/(dx*wx),  a_N = a_S = 1/(dy*wy),
%         wx = dx (dx/2 on a W/E edge),  wy = dy (dy/2 on an S/N edge),
%         diag = sum of the a_nb that exist (+ hl),
%         rho_p = (eps*sigma/k)*(1/wx per radiating W/E face + 1/wy per radiating S/N face).
%   Interior nodes reduce to the classic 5-point stencil.
%
%   The matrix A holds the LINEAR part only (rows are diagonally dominant, so
%   Jacobi / Gauss-Seidel / SOR apply). The radiative term is returned
%   separately in nodes.rho, so the full nonlinear system is
%
%       F(T) = A*T - b + nodes.rho .* (T.^4 - Tamb^4) = 0.
%
%   For a radiating node p the scalar equation is (for the radiative root-find)
%       T - Tnb + beta*(T^4 - Tamb^4) = 0,   Tnb = (b_p - A(p,:)*T + A(p,p)*T_p)/A(p,p),
%                                            beta = rho_p/A(p,p).
%
% INPUT : P    parameters struct (needs Lx..Ny, k, Tamb, eps, sigma).
%         grid struct from build_grid.
%   Optional fields of P (all default to "off" so parameters.m still works):
%     P.bc.S/N/W/E : struct('type',T,'value',V), T = 'dirichlet'|'radiative'|'insulated'.
%                    value (Dirichlet only) = scalar, vector along the edge, or f(s).
%                    Default: south edge Dirichlet at P.Thot, other three edges radiative.
%                    Corners: a Dirichlet edge wins (order W, E, S, N, last wins).
%     P.q          heat generation [W/m^3]: scalar, Ny x Nx matrix, or f(X,Y) (vectorized).
%     P.h, P.t     top-face loss: h [W/m^2-K] over plate thickness t [m]; adds
%                  hl = h/(k*t) [1/m^2] to the diagonal. Default off.
% OUTPUT: A (N x N sparse), b (N x 1), nodes (struct):
%     nodes.type  N x 1: 0 interior, 1 Dirichlet, 2 radiating boundary, 3 other boundary
%     nodes.rho   N x 1 radiation factor [1/(m^2 K^3)], zero if not radiating
%     nodes.vol   N x 1 control-volume area (wx*wy) [m^2]
%     nodes.hl    top-face loss coefficient [1/m^2]
%     nodes.Psrc  total generated power per unit depth [W/m]
%     nodes.q     Ny x Nx generation field actually used [W/m^3]

    Nx = grid.Nx;  Ny = grid.Ny;  N = grid.N;
    dx = grid.dx;  dy = grid.dy;
    k  = P.k;

    bc = get_bc(P);

    hl = 0;
    if isfield(P, 'h') && ~isempty(P.h) && P.h > 0
        if ~isfield(P, 't') || isempty(P.t) || P.t <= 0
            error('assemble_system:face', 'P.h needs a positive plate thickness P.t.');
        end
        hl = P.h/(k*P.t);
    end

    q = get_source(P, grid);

    % Dirichlet values on the edges (NaN = not Dirichlet). Order: W, E, S, N.
    Tdir = nan(Ny, Nx);
    if strcmp(bc.W.type, 'dirichlet'), Tdir(:, 1)  = edge_values(bc.W, grid.y); end
    if strcmp(bc.E.type, 'dirichlet'), Tdir(:, Nx) = edge_values(bc.E, grid.y); end
    if strcmp(bc.S.type, 'dirichlet'), Tdir(1, :)  = edge_values(bc.S, grid.x).'; end
    if strcmp(bc.N.type, 'dirichlet'), Tdir(Ny, :) = edge_values(bc.N, grid.x).'; end

    radW = strcmp(bc.W.type, 'radiative');
    radE = strcmp(bc.E.type, 'radiative');
    radS = strcmp(bc.S.type, 'radiative');
    radN = strcmp(bc.N.type, 'radiative');

    I = zeros(5*N, 1);  J = zeros(5*N, 1);  V = zeros(5*N, 1);
    c = 0;
    b = zeros(N, 1);
    nodes.type = zeros(N, 1);
    nodes.rho  = zeros(N, 1);
    nodes.vol  = zeros(N, 1);

    for j = 1:Ny
        for i = 1:Nx
            p = grid.idx(j, i);

            wx = dx;  if i == 1 || i == Nx, wx = dx/2; end
            wy = dy;  if j == 1 || j == Ny, wy = dy/2; end
            nodes.vol(p) = wx*wy;

            if ~isnan(Tdir(j, i))                    % Dirichlet node
                c = c + 1;  I(c) = p;  J(c) = p;  V(c) = 1;
                b(p) = Tdir(j, i);
                nodes.type(p) = 1;
                continue;
            end

            ax = 1/(dx*wx);
            ay = 1/(dy*wy);
            dg = hl;
            if i > 1,  c = c+1; I(c) = p; J(c) = grid.idx(j, i-1); V(c) = -ax; dg = dg + ax; end
            if i < Nx, c = c+1; I(c) = p; J(c) = grid.idx(j, i+1); V(c) = -ax; dg = dg + ax; end
            if j > 1,  c = c+1; I(c) = p; J(c) = grid.idx(j-1, i); V(c) = -ay; dg = dg + ay; end
            if j < Ny, c = c+1; I(c) = p; J(c) = grid.idx(j+1, i); V(c) = -ay; dg = dg + ay; end
            c = c + 1;  I(c) = p;  J(c) = p;  V(c) = dg;

            r = 0;                                   % radiating faces / volume
            if i == 1  && radW, r = r + 1/wx; end
            if i == Nx && radE, r = r + 1/wx; end
            if j == 1  && radS, r = r + 1/wy; end
            if j == Ny && radN, r = r + 1/wy; end
            nodes.rho(p) = (P.eps*P.sigma/k)*r;

            b(p) = q(j, i)/k + hl*P.Tamb;

            onEdge = (i == 1 || i == Nx || j == 1 || j == Ny);
            if r > 0
                nodes.type(p) = 2;
            elseif onEdge
                nodes.type(p) = 3;
            end
        end
    end

    A = sparse(I(1:c), J(1:c), V(1:c), N, N);

    nodes.hl   = hl;
    nodes.q    = q;
    qv = reshape(q.', [], 1);
    nodes.Psrc = sum(qv .* nodes.vol);               % trapezoid-rule integral of q
end

% ------------------------------------------------------------------------
function bc = get_bc(P)
    if isfield(P, 'bc')
        bc = P.bc;
    else
        bc.S = struct('type', 'dirichlet', 'value', P.Thot);
        bc.N = struct('type', 'radiative');
        bc.W = struct('type', 'radiative');
        bc.E = struct('type', 'radiative');
    end
    names = {'S', 'N', 'W', 'E'};
    for m = 1:4
        if ~isfield(bc, names{m}) || ~isfield(bc.(names{m}), 'type')
            error('assemble_system:bc', 'P.bc.%s.type is missing.', names{m});
        end
        t = bc.(names{m}).type;
        if ~any(strcmp(t, {'dirichlet', 'radiative', 'insulated'}))
            error('assemble_system:bc', 'Unknown BC type "%s" on edge %s.', t, names{m});
        end
        if strcmp(t, 'dirichlet') && ~isfield(bc.(names{m}), 'value')
            error('assemble_system:bc', 'Dirichlet edge %s needs a value.', names{m});
        end
    end
end

function v = edge_values(spec, s)
% Dirichlet values along an edge. s = coordinate vector along that edge.
    val = spec.value;
    if isa(val, 'function_handle')
        v = val(s);
    elseif isscalar(val)
        v = val*ones(size(s));
    else
        v = val;
    end
    if numel(v) ~= numel(s)
        error('assemble_system:bc', 'Dirichlet value has %d entries, edge has %d nodes.', numel(v), numel(s));
    end
    v = v(:);                                        % column vector
end

function q = get_source(P, grid)
    sz = [grid.Ny, grid.Nx];
    if ~isfield(P, 'q') || isempty(P.q)
        q = zeros(sz);
    elseif isa(P.q, 'function_handle')
        q = P.q(grid.X, grid.Y);
        if isscalar(q), q = q*ones(sz); end
    elseif isscalar(P.q)
        q = P.q*ones(sz);
    else
        q = P.q;
    end
    if ~isequal(size(q), sz)
        error('assemble_system:source', 'P.q must be a scalar, %d x %d matrix, or f(X,Y).', sz(1), sz(2));
    end
end
