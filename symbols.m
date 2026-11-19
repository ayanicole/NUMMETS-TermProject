
function S = symbols()
% SYMBOLS
% Centralized documentation of all mathematical symbols.
% OUTPUT: S (struct containing symbols, meanings, units, and MATLAB names)

% Temperature
S.T = struct('symbol', 'T(x,y)', 'meaning', 'Temperature', 'unit', 'K');

S.Thot = struct('symbol', 'T_hot', 'meaning', 'Heated boundary temperature', 'unit', 'K');

S.Tamb = struct('symbol', 'T_amb', 'meaning', 'Ambient temperature', 'unit', 'K');

% Geometry
S.Lx = struct('symbol', 'L_x', 'meaning', 'Plate length', 'unit', 'm');

S.Ly = struct('symbol', 'L_y', 'meaning', 'Plate width', 'unit', 'm');

S.Nx = struct('symbol', 'N_x', 'meaning', 'Number of grid nodes in x', 'unit', 'none');

S.Ny = struct('symbol', 'N_y', 'meaning', 'Number of grid nodes in y', 'unit', 'none');

S.dx = struct('symbol', 'Delta x', 'meaning', 'Grid spacing in x', 'unit', 'm');

S.dy = struct('symbol', 'Delta y', 'meaning', 'Grid spacing in y', 'unit', 'm');

% Material properties
S.k = struct('symbol', 'k', 'meaning', 'Thermal conductivity', 'unit', 'W/(m*K)');

S.eps = struct('symbol', 'epsilon', 'meaning', 'Surface emissivity', 'unit', 'none');

S.sigma = struct('symbol', 'sigma', 'meaning', 'Stefan-Boltzmann constant', 'unit', 'W/(m^2*K^4)');

% Heat transfer
S.q = struct('symbol', 'q''''', 'meaning', 'Heat flux', 'unit', 'W/m^2');

S.Q = struct('symbol', 'Q', 'meaning', 'Total heat flow across an edge', 'unit', 'W');

% System equations
S.A = struct('symbol', 'A', 'meaning', 'System matrix', 'unit', 'none');

S.b = struct('symbol', 'b', 'meaning', 'Right-hand side vector', 'unit', 'none');

% Solver settings
S.omega = struct('symbol', 'omega', 'meaning', 'SOR relaxation factor', 'unit', 'none');

S.tol = struct('symbol', 'tol', 'meaning', 'Outer iteration stopping tolerance', 'unit', 'K');

S.tolBC = struct('symbol', 'tol_BC', 'meaning', 'Radiative boundary solve tolerance', 'unit', 'K');

% Error analysis
S.res = struct('symbol', 'r^(k)', 'meaning', 'Residual at iteration k', 'unit', 'K or W/m^2');

S.maxabs = struct('symbol', 'e_max', 'meaning', 'Maximum absolute error', 'unit', 'K');

S.rel = struct('symbol', 'e_rel', 'meaning', 'Relative error', 'unit', 'none');

S.L2 = struct('symbol', 'e_2', 'meaning', 'L2 norm of the error', 'unit', 'K');

S.p = struct('symbol', 'p', 'meaning', 'Observed order of accuracy', 'unit', 'none');

end
