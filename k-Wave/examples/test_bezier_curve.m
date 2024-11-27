
% Grid size and spacing
Nx_tot = 108; 
Ny_tot = 150; 
Nz_tot = 54;  
dx = 1e-3;  % Grid spacing in x-direction (1 mm)
dy = 1e-3;  % Grid spacing in y-direction (1 mm)
dz = 1e-3;  % Grid spacing in z-direction (1 mm)

% Control points for Bézier curve (relative to grid center)
P0 = [54, 25, 27];  % Start point in px
P1 = [60, 90, 27];  % Control point px
P2 = [72, 75, 27];  % End point px

% Desired physical radius around the curve (in meters)
physical_radius = 0.01;  % 1 cm
scattering_map = randn([Nx_tot, Ny_tot, Nz_tot]);
% Generate scattering region
[scattering_map,bezier_curve] = makeBezierScatteringRegion(Nx_tot, Ny_tot, Nz_tot, dx, dy, dz, P0, P1, P2, physical_radius,scattering_map);

% Visualization of the scattering map (cross-section)
figure;
imagesc(scattering_map(:, :, round(Nz_tot / 2))); 
axis equal;
colormap(gray);
title('Scattering Shape from Bézier Curve');
xlabel('X Position [grid points]');
ylabel('Y Position [grid points]');

% Visualization of the Bézier curve in 3D to check its shape
figure;
plot3(bezier_curve(:, 1), bezier_curve(:, 2), bezier_curve(:, 3), '-o');
axis equal;
xlabel('X Position [grid points]');
ylabel('Y Position [grid points]');
zlabel('Z Position [grid points]');
title('3D Bézier Curve');
grid on;