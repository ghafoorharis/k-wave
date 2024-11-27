function scattering_map = makeBezierScatteringRegion(Nx, Ny, Nz, dx, dy, dz, P0, P1, P2, radius, plot_ball, binary)
% MAKEBEZIERSCATTERINGREGION Create a binary map of a scattering region around a 3D Bézier curve.
%
% DESCRIPTION:
%     makeBezierScatteringRegion creates a binary map of a scattering region 
%     within a 3D grid, where the region is defined by a 3D Bézier curve.
%     The scattering region will be defined by all points within a specified 
%     radius of the curve, and those points will be marked as `1` (or `true` 
%     if using logicals).
%
% USAGE:
%     scattering_map = makeBezierScatteringRegion(Nx, Ny, Nz, dx, dy, dz, P0, P1, P2, radius)
%     scattering_map = makeBezierScatteringRegion(Nx, Ny, Nz, dx, dy, dz, P0, P1, P2, radius, plot_ball)
%     scattering_map = makeBezierScatteringRegion(Nx, Ny, Nz, dx, dy, dz, P0, P1, P2, radius, [], binary)
%
% INPUTS:
%     Nx, Ny, Nz      - size of the 3D grid [grid points]
%     dx, dy, dz      - grid spacing in x, y, and z directions [m]
%     P0, P1, P2      - control points of the Bézier curve [grid points]
%     radius          - scattering region radius around the Bézier curve [grid points]
%
% OPTIONAL INPUTS:
%     plot_ball       - Boolean controlling whether the scattering region is
%                       plotted using voxelPlot (default = false)
%     binary          - Boolean controlling whether the map is returned as a
%                       double precision matrix (false) or a logical matrix (true) 
%                       (default = false)
%
% OUTPUTS:
%     scattering_map  - 3D binary map of the scattering region
%
% ABOUT:
%     author          - Bradley Treeby (adapted from makeBall)
%     date            - 1st September 2024
%
% This file is part of the k-Wave Toolbox (http://www.k-wave.org)

% Define literal
MAGNITUDE = 1;

% Check for plot_ball input
if nargin < 11 || isempty(plot_ball)
    plot_ball = false;
end

% Check for binary input
if nargin < 12 || isempty(binary)
    binary = false;
end

% Force integer values
Nx = round(Nx);
Ny = round(Ny);
Nz = round(Nz);
P0 = round(P0);
P1 = round(P1);
P2 = round(P2);

% Create empty matrix
if binary
    scattering_map = false(Nx, Ny, Nz);
else
    scattering_map = zeros(Nx, Ny, Nz);
end

% Define the t parameter for the Bézier curve (from 0 to 1)
t = linspace(0, 1, 100);  % Number of points along the Bézier curve

% Generate the Bézier curve in 3D
bezier_curve_x = (1 - t).^2 * P0(1) + 2 * (1 - t) .* t * P1(1) + t.^2 * P2(1);
bezier_curve_y = (1 - t).^2 * P0(2) + 2 * (1 - t) .* t * P1(2) + t.^2 * P2(2);
bezier_curve_z = (1 - t).^2 * P0(3) + 2 * (1 - t) .* t * P1(3) + t.^2 * P2(3);

% Combine x, y, and z coordinates into the Bézier curve array
bezier_curve = [bezier_curve_x', bezier_curve_y', bezier_curve_z'];

% Create a distance matrix for the grid points
[x_grid, y_grid, z_grid] = meshgrid(1:Nx, 1:Ny, 1:Nz);
distance_matrix = sqrt((x_grid(:) - bezier_curve(:,1)').^2 + (y_grid(:) - bezier_curve(:,2)').^2 + (z_grid(:) - bezier_curve(:,3)').^2);

% Mark points within the scattering region
for i = 1:length(t)
    % Find grid points within the defined radius of the Bézier curve
    scattering_map(distance_matrix(:, i) <= radius) = MAGNITUDE;
end

% Reshape the scattering map back to the 3D grid
scattering_map = reshape(scattering_map, Nx, Ny, Nz);

% Shift centre if required (optional)
if plot_ball
    voxelPlot(double(scattering_map));
end

end
