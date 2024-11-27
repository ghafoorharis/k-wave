function [sound_speed_map, density_map, scattering_bezier_map] = define_phantom_bezier(c0, ...
            rho0, Nx_tot, Ny_tot,Nz_tot, dx, ...
            dy,dz,...
            background_map_mean,background_map_std, ...
            radius, ...
            x_shift,y_shift, ...
            scattering_mean_c1, ...
            scattering_std_c2, ...
            scattering_divider_c3)
    % Get the parameters properly
    constant_scattering_region_add = scattering_mean_c1;
    constant_scattering_region_multiplier = scattering_std_c2;
    constant_scattering_region_divider = scattering_divider_c3;
    % define a random distribution of scatterers for the medium
    background_map = background_map_mean + background_map_std * randn([Nx_tot, ...
                        Ny_tot, Nz_tot]);
    % define a random distribution of scatterers for the highly scattering
    scattering_map = randn([Nx_tot, Ny_tot, Nz_tot]);
    scattering_c0 = c0 + constant_scattering_region_add + constant_scattering_region_multiplier * scattering_map;
    scattering_rho0 = scattering_c0 / constant_scattering_region_divider;
   
    % Scattering region: initialize uniform medium
    sound_speed_map = c0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
    density_map = rho0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
    % Define Bezier Curve Params
    % Ensure the radius is within the max allowed limit
    max_radius = 13;  % Maximum radius allowed (e.g., 11 mm)
    phantom_radius = min(radius, max_radius);
    % Define the control points for the Bézier curve (circular arc)
    % Position of the Bézier curve in the center of the grid
    center_X = Nx_tot / 2;  % Center x-coordinate
    center_Y = Ny_tot / 2;  % Center y-coordinate
    center_Z = Nz_tot / 2;  % Center z-coordinate
    P0 = [center_X - x_shift, center_Y - y_shift, center_Z];  % Start point on the circle
    P2 = [center_X, center_Y, center_Z];  % End point on the circle
    P1 = [center_X + x_shift / 2, center_Y + y_shift, center_Z];  % Control point for the arc
    % ---- Define Shapes ----
    % scattering_circle = makeBall(Nx_tot, Ny_tot, Nz_tot, ...
    %                     round(x_pos/dx), round(y_pos/dx), ...
    %                     Nz_tot/2, round(radius/dx));
    % Generate the phantom using Bézier curve (function you defined)
    [scattering_bezier_map, ~] = makeBezierScatteringRegion(Nx_tot, Ny_tot, Nz_tot, dx, dy, dz, P0, P1, P2, ...
        phantom_radius, scattering_map);
    % Modify sound speed and density for each shape
    sound_speed_map(scattering_bezier_map == 1) = scattering_c0(scattering_bezier_map == 1);
    density_map(scattering_bezier_map == 1) = scattering_rho0(scattering_bezier_map == 1);
end
