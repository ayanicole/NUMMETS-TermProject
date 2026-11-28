function test_model()
% TEST_MODEL  Self-checks for the Model Lead functions.
% Run from the repo folder (needs parameters.m, build_grid.m, assemble_system.m,
% heat_flow.m, material_k.m, interp_field.m on the path):   >> test_model
% Prints PASS/FAIL per check and errors out at the end if anything failed.

    ok = true;

    % ---- 1. grid and indexing contract --------------------------------
    P = base_params();  P.Nx = 5;  P.Ny = 4;  P.Lx = 1;  P.Ly = 0.6;
    g = build_grid(P);
    Tvec = (1:g.N).';
    Tmap = reshape(Tvec, g.Nx, g.Ny).';
    ok = rep('grid: size, spacing', g.N == 20 && abs(g.dx - 0.25) < 1e-14 && abs(g.dy - 0.2) < 1e-14, '') && ok;
    ok = rep('grid: idx(j,i) = i + (j-1)*Nx', g.idx(3, 2) == 12 && Tmap(3, 2) == 12, '') && ok;

    % ---- 2. quadratic patch test (interior stencil + source term) ------
    P = base_params();  P.Nx = 13;  P.Ny = 9;  P.Lx = 0.03;  P.Ly = 0.02;
    Tex = @(x, y) 300 + 1e4*x.^2 + 2e4*y.^2;            % laplacian = 6e4
    P.q = -P.k*6e4;
    P.bc.S = struct('type', 'dirichlet', 'value', @(s) Tex(s, 0));
    P.bc.N = struct('type', 'dirichlet', 'value', @(s) Tex(s, P.Ly));
    P.bc.W = struct('type', 'dirichlet', 'value', @(s) Tex(0, s));
    P.bc.E = struct('type', 'dirichlet', 'value', @(s) Tex(P.Lx, s));
    g = build_grid(P);
    [A, b, nd] = assemble_system(P, g);
    T = A\b;
    err = max(abs(T - reshape(Tex(g.X, g.Y).', [], 1)));
    ok = rep('patch: quadratic reproduced exactly', err < 1e-9, sprintf('(max err %.2e K)', err)) && ok;
    ok = rep('patch: Psrc = q*Lx*Ly', abs(nd.Psrc - P.q*P.Lx*P.Ly) < 1e-9*abs(P.q*P.Lx*P.Ly), '') && ok;
    ok = rep('patch: A diagonally dominant', all(2*abs(diag(A)) >= sum(abs(A), 2) - 1e-9), '') && ok;

    % ---- 3. analytical Laplace case: second-order convergence ----------
    %  T = 0 on W, E, S;  T = T0*sin(pi*x/L) on N  ->  T = T0 sin(pi x/L) sinh(pi y/L)/sinh(pi)
    T0 = 100;  L = 0.02;
    Ns = [11 21 41 81];
    e  = zeros(size(Ns));
    for m = 1:numel(Ns)
        [T, g, P] = laplace_case(Ns(m), T0, L);
        Tx = T0*sin(pi*g.X/L).*sinh(pi*g.Y/L)/sinh(pi);
        e(m) = max(abs(T - reshape(Tx.', [], 1)));
    end
    p = log(e(1:end-1)./e(2:end))/log(2);
    ok = rep('laplace: observed order ~ 2', all(p > 1.9 & p < 2.1), sprintf('(orders %s)', num2str(p, '%.3f '))) && ok;

    % ---- 4. heat flow: trapezoid / Simpson vs exact, energy balance -----
    Qex = -2*400*T0/tanh(pi);                             % heat enters through the hot N edge
    eTr = zeros(1, 2);  eSi = zeros(1, 2);
    for m = 1:2
        [T, g, P] = laplace_case(Ns(m + 2), T0, L);       % N = 41, 81
        eTr(m) = abs(heat_flow(T, g, P, 'trap', 'N')    - Qex)/abs(Qex);
        eSi(m) = abs(heat_flow(T, g, P, 'simpson', 'N') - Qex)/abs(Qex);
    end
    Qsum = zeros(1, 2);
    for m = 1:2
        [T, g, P] = laplace_case(Ns(m + 2), T0, L);       % N = 41, 81
        for ed = 'SNWE', Qsum(m) = Qsum(m) + heat_flow(T, g, P, 'simpson', ed); end
    end
    ok = rep('heat_flow: trap rel. error < 5e-3 (N=41)',    eTr(1) < 5e-3, sprintf('(%.2e)', eTr(1))) && ok;
    ok = rep('heat_flow: simpson rel. error < 5e-3 (N=41)', eSi(1) < 5e-3, sprintf('(%.2e)', eSi(1))) && ok;
    ok = rep('heat_flow: error falls ~4x per refinement',   eSi(1)/eSi(2) > 3 && eSi(1)/eSi(2) < 5, sprintf('(ratio %.2f)', eSi(1)/eSi(2))) && ok;
    ok = rep('heat_flow: edge flows sum to ~0, shrinking as h^2', abs(Qsum(1)) < 5e-3*abs(Qex) && ...
             abs(Qsum(1)/Qsum(2)) > 3 && abs(Qsum(1)/Qsum(2)) < 5, ...
             sprintf('(sum %.2e -> %.2e W/m, ratio %.2f)', Qsum(1), Qsum(2), abs(Qsum(1)/Qsum(2)))) && ok;
    [To, go, Po] = laplace_case(40, T0, L);               % 40 nodes = 39 intervals (odd)
    msg = '';
    try
        heat_flow(To, go, Po, 'simpson', 'N');
    catch ex
        msg = ex.identifier;
    end
    ok = rep('heat_flow: Simpson rejects odd interval count', strcmp(msg, 'heat_flow:simpson'), '') && ok;

    % ---- 5. radiative half-cell nodes (manufactured solution) -----------
    %  T = A0 + B x^2 with E edge radiating, W/N/S insulated, q = -2 k B.
    P = base_params();  P.Nx = 9;  P.Ny = 7;  P.Lx = 0.02;  P.Ly = 0.02;     % dx ~= dy on purpose
    qE = P.eps*P.sigma*(1000^4 - P.Tamb^4);
    B  = -qE/(2*P.k*P.Lx);   A0 = 1000 - B*P.Lx^2;
    P.q = -2*P.k*B;
    P.bc.W = struct('type', 'insulated');  P.bc.N = struct('type', 'insulated');
    P.bc.S = struct('type', 'insulated');  P.bc.E = struct('type', 'radiative');
    g = build_grid(P);
    [A, b, nd] = assemble_system(P, g);
    T = 1000*ones(g.N, 1);
    for it = 1:40
        F  = A*T - b + nd.rho.*(T.^4 - P.Tamb^4);
        Jm = A + spdiags(4*nd.rho.*T.^3, 0, g.N, g.N);
        dT = -(Jm\F);   T = T + dT;
        if max(abs(dT)) < 1e-9, break; end
    end
    err = max(abs(T - reshape((A0 + B*g.X.^2).', [], 1)));
    ok = rep('radiative: manufactured solution recovered', err < 1e-8, sprintf('(max err %.2e K, %d Newton its)', err, it)) && ok;
    Qs = 0;  for ed = 'SNWE', Qs = Qs + heat_flow(T, g, P, 'trap', ed); end
    ok = rep('radiative: sum of edge Q = generated power', abs(Qs - nd.Psrc) < 1e-6*abs(nd.Psrc), sprintf('(%.4f vs %.4f W/m)', Qs, nd.Psrc)) && ok;

    % ---- 5b. all four edges radiate: exact discrete energy balance (corners!) ----
    P = base_params();  P.Nx = 21;  P.Ny = 13;  P.Lx = 0.02;  P.Ly = 0.012;
    P.bc.W = struct('type', 'radiative');  P.bc.E = struct('type', 'radiative');
    P.bc.S = struct('type', 'radiative');  P.bc.N = struct('type', 'radiative');
    P.q = @(x, y) 4e7*(abs(x - 0.01) <= 0.004 & abs(y - 0.006) <= 0.002);   % hot spot
    g = build_grid(P);
    [A, b, nd] = assemble_system(P, g);
    T = 400*ones(g.N, 1);
    for it = 1:60
        F  = A*T - b + nd.rho.*(T.^4 - P.Tamb^4);
        Jm = A + spdiags(4*nd.rho.*T.^3, 0, g.N, g.N);
        dT = -(Jm\F);   T = T + dT;
        if max(abs(dT)) < 1e-9, break; end
    end
    Qr = 0;  Qg = 0;
    for ed = 'SNWE'
        Qr = Qr + heat_flow(T, g, P, 'trap', ed, 'radiation');
        Qg = Qg + heat_flow(T, g, P, 'simpson', ed);
    end
    ok = rep('4 radiating edges: sum Q_rad = Psrc exactly', abs(Qr - nd.Psrc) < 1e-8*nd.Psrc, ...
             sprintf('(%.6f vs %.6f W/m, %d Newton its, Tmax %.1f K)', Qr, nd.Psrc, it, max(T))) && ok;
    ok = rep('4 radiating edges: gradient fluxes agree within 5%', abs(Qg - nd.Psrc) < 5e-2*nd.Psrc, ...
             sprintf('(%.3f vs %.3f W/m)', Qg, nd.Psrc)) && ok;

    % ---- 5c. make_source: exact total power; radiative_bc inside a nodal sweep ----
    P = base_params();  P.Nx = 41;  P.Ny = 41;  P.t = 0.0025;
    g = build_grid(P);
    okS = true;
    for kind = {'uniform', 'rect', 'gauss'}
        switch kind{1}
            case 'uniform', q = make_source(g, 'uniform', 10, P.t);
            case 'rect',    q = make_source(g, 'rect', 10, P.t, 0.01, 0.01, 0.002);
            case 'gauss',   q = make_source(g, 'gauss', 10, P.t, 0.01, 0.01, 0.0015);
        end
        P.q = q;
        [~, ~, nd] = assemble_system(P, g);
        okS = okS && abs(nd.Psrc*P.t - 10) < 1e-12*10;
    end
    ok = rep('make_source: total power exactly 10 W (uniform/rect/gauss)', okS, '') && ok;
    q = make_source(g, 'rect', 10, P.t, 0.01, 0.01, 0.002);
    ok = rep('make_source: rect covers 9 x 9 nodes symmetrically', nnz(q) == 81, '') && ok;

    Ts = 350;  bt = 3e-10;  Tn = Ts + bt*(Ts^4 - 298^4);   P.Tamb = 298;
    [T1, n1] = radiative_bc(300,  P, 1e-12, Tn, bt);
    [T2, n2] = radiative_bc(1500, P, 1e-12, Tn, bt);
    ok = rep('radiative_bc: root found from left and right guesses', abs(T1 - Ts) < 1e-9 && abs(T2 - Ts) < 1e-9 && n1 < 15 && n2 < 15, ...
             sprintf('(%d and %d steps)', n1, n2)) && ok;

    % nodal Gauss-Seidel sweep with radiative_bc, on the manufactured case of check 5
    P = base_params();  P.Nx = 9;  P.Ny = 7;  P.Lx = 0.02;  P.Ly = 0.02;
    qE = P.eps*P.sigma*(1000^4 - P.Tamb^4);
    B  = -qE/(2*P.k*P.Lx);   A0 = 1000 - B*P.Lx^2;
    P.q = -2*P.k*B;
    P.bc.W = struct('type', 'dirichlet', 'value', A0);   % anchored: well-posed for a sweep
    P.bc.N = struct('type', 'insulated');
    P.bc.S = struct('type', 'insulated');  P.bc.E = struct('type', 'radiative');
    g = build_grid(P);
    [A, b, nd] = assemble_system(P, g);
    T = 1000*ones(g.N, 1);
    for sweep = 1:5000
        Told = T;
        for p = 1:g.N
            d   = A(p, p);
            Tnb = (b(p) - A(p, :)*T + d*T(p))/d;
            if nd.rho(p) > 0
                T(p) = radiative_bc(T(p), P, 1e-12, Tnb, nd.rho(p)/d);
            else
                T(p) = Tnb;
            end
        end
        if max(abs(T - Told)) < 1e-10, break; end
    end
    err = max(abs(T - reshape((A0 + B*g.X.^2).', [], 1)));
    ok = rep('nodal GS sweep + radiative_bc recovers manufactured solution', err < 1e-6, ...
             sprintf('(max err %.2e K after %d sweeps)', err, sweep)) && ok;

    % ---- 6. top-face loss -----------------------------------------------
    P = base_params();  P.Nx = 7;  P.Ny = 7;
    P.bc.W = struct('type', 'insulated');  P.bc.E = struct('type', 'insulated');
    P.bc.S = struct('type', 'insulated');  P.bc.N = struct('type', 'insulated');
    P.q = 5e6;  P.h = 20;  P.t = 0.002;
    g = build_grid(P);
    [A, b] = assemble_system(P, g);
    T = A\b;
    ok = rep('face loss: uniform T = Tamb + q*t/h', max(abs(T - (P.Tamb + P.q*P.t/P.h))) < 1e-8, '') && ok;

    % ---- 7. interpolation -------------------------------------------------
    P = base_params();  P.Nx = 15;  P.Ny = 12;  P.Lx = 1;  P.Ly = 1;
    g = build_grid(P);
    f = @(x, y) 1 + 2*x + 3*y + x.^2.*y + y.^3;           % degree 3
    xq = [0.013 0.5 0.987 0.31 0.77];   yq = [0.02 0.41 0.98 0.6 0.05];
    errI = max(abs(interp_field(f(g.X, g.Y), g, xq, yq, 3) - f(xq, yq)));
    ok = rep('Lagrange: cubic field reproduced', errI < 1e-10, sprintf('(max err %.2e)', errI)) && ok;
    kc = material_k('copper', [273 400 330]);   ka = material_k('aluminum', [273 400 330]);
    ok = rep('Newton DD: table points reproduced', abs(kc(1) - 401) < 1e-9 && abs(kc(2) - 392) < 1e-9 ...
             && abs(ka(1) - 236) < 1e-9 && abs(ka(2) - 240) < 1e-9, '') && ok;
    ok = rep('Newton DD: k(330 K) between neighbours', kc(3) > 392 && kc(3) < 401 && ka(3) > 235 && ka(3) < 241, ...
             sprintf('(Cu %.2f, Al %.2f)', kc(3), ka(3))) && ok;
    caught = false;
    try
        material_k('copper', 250);
    catch
        caught = true;
    end
    ok = rep('Newton DD: refuses extrapolation', caught, '') && ok;

    % ---- 8. the repo's own parameters.m ------------------------------------
    P = parameters();
    g = build_grid(P);
    [A, b, nd] = assemble_system(P, g);
    ok = rep('defaults: sizes', size(A, 1) == 441 && numel(b) == 441, '') && ok;
    ok = rep('defaults: south edge Dirichlet at Thot', all(nd.type(g.idx(1, :)) == 1) && all(b(g.idx(1, :)) == P.Thot), '') && ok;
    ok = rep('defaults: W/E/N edge nodes radiate', all(nd.type(g.idx(2:end, 1)) == 2) && all(nd.type(g.idx(2:end, end)) == 2) ...
             && all(nd.type(g.idx(end, 2:end-1)) == 2), '') && ok;
    rEdge = nd.rho(g.idx(end, 5));  rCorner = nd.rho(g.idx(end, 1));
    ok = rep('defaults: corner radiates over two half-faces', abs(rCorner/rEdge - 2) < 1e-12, '') && ok;
    T = P.Thot*ones(g.N, 1);
    for it = 1:40
        F  = A*T - b + nd.rho.*(T.^4 - P.Tamb^4);
        Jm = A + spdiags(4*nd.rho.*T.^3, 0, g.N, g.N);
        dT = -(Jm\F);   T = T + dT;
        if max(abs(dT)) < 1e-9, break; end
    end
    fprintf('       defaults: T range %.4f .. %.4f K (spread %.4f K)\n', min(T), max(T), max(T) - min(T));
    QS  = heat_flow(T, g, P, 'simpson', 'S');
    Qrd = 0;  for ed = 'WEN', Qrd = Qrd + heat_flow(T, g, P, 'simpson', ed, 'radiation'); end
    fprintf('       defaults: heat in via S edge %.3f W/m, radiated from W/E/N %.3f W/m (mismatch %.2f%%)\n', ...
            -QS, Qrd, 100*abs(-QS - Qrd)/Qrd);

    if ok
        fprintf('\nALL CHECKS PASSED\n');
    else
        error('test_model:fail', 'Some checks failed - see FAIL lines above.');
    end
end

% ------------------------------------------------------------------------
function P = base_params()
    P.Lx = 0.02;  P.Ly = 0.02;  P.Nx = 21;  P.Ny = 21;
    P.k = 400;  P.Thot = 358;  P.Tamb = 298;  P.eps = 0.9;  P.sigma = 5.670e-8;
end

function [T, g, P] = laplace_case(N, T0, L)
    P = base_params();  P.Nx = N;  P.Ny = N;  P.Lx = L;  P.Ly = L;  P.q = 0;
    P.bc.S = struct('type', 'dirichlet', 'value', 0);
    P.bc.W = struct('type', 'dirichlet', 'value', 0);
    P.bc.E = struct('type', 'dirichlet', 'value', 0);
    P.bc.N = struct('type', 'dirichlet', 'value', @(s) T0*sin(pi*s/L));
    g = build_grid(P);
    [A, b] = assemble_system(P, g);
    T = A\b;
end

function r = rep(name, cond, msg)
    r = all(logical(cond(:)));
    if r, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('[%s] %s %s\n', tag, name, msg);
end
