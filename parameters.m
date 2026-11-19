function P = parameters()
% Shared parameters. All constants live here and nowhere else.
% OUTPUT: P (struct of constants)

% Geometry
P.Lx = 0.02;  P.Ly = 0.02;      % plate size [m]
P.Nx = 21;    P.Ny = 21;        % grid nodes (20 intervals, even for Simpson)

% Material
P.k = 400;                      % thermal conductivity [W/m-K] (copper)

% Boundary conditions
P.Thot = 358;                   % hot boundary temperature [K]
P.Tamb = 298;                   % ambient temperature [K]
P.eps = 0.9;                    % emissivity [-]
P.sigma = 5.670e-8;             % Stefan-Boltzmann constant [W/m^2-K^4]

% Solver settings
P.tol = 1e-6;                   % stopping tolerance (change between iterations)
P.tolBC = 1e-10;                % radiative BC tolerance (tighter than P.tol)
P.maxIter = 1e5;                % max iterations
P.omega = 1.8;                  % SOR relaxation factor
end
