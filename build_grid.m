function grid = build_grid(P)
% BUILD_GRID  Uniform 2D grid for the chip heat spreader.
% INPUT : P (parameters struct: Lx, Ly, Nx, Ny)
% OUTPUT: grid (struct)
%   grid.Nx, grid.Ny   nodes in x and y
%   grid.N             total number of nodes (Nx*Ny)
%   grid.x, grid.y     node coordinates, 1 x Nx and 1 x Ny [m]
%   grid.dx, grid.dy   node spacing [m]
%   grid.X, grid.Y     meshgrid arrays, Ny x Nx (use with surf/contourf)
%   grid.idx           Ny x Nx matrix of linear indices
%
% INDEXING CONTRACT (everyone uses this):
%   node (i,j), i = 1..Nx along x, j = 1..Ny along y, has linear index
%   p = i + (j-1)*Nx = grid.idx(j,i).
%   Vector -> map :  Tmap = reshape(Tvec, grid.Nx, grid.Ny).';   % Ny x Nx
%   Map -> vector :  Tvec = reshape(Tmap.', [], 1);              % N x 1
%   Row j of Tmap is at y = grid.y(j); column i is at x = grid.x(i).

    Nx = P.Nx;  Ny = P.Ny;
    if Nx < 3 || Ny < 3
        error('build_grid:size', 'Need at least 3 nodes in each direction.');
    end

    grid.Nx = Nx;
    grid.Ny = Ny;
    grid.N  = Nx*Ny;
    grid.x  = linspace(0, P.Lx, Nx);
    grid.y  = linspace(0, P.Ly, Ny);
    grid.dx = P.Lx/(Nx - 1);
    grid.dy = P.Ly/(Ny - 1);
    [grid.X, grid.Y] = meshgrid(grid.x, grid.y);
    grid.idx = reshape(1:grid.N, Nx, Ny).';
end
