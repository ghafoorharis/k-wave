function [scattering_map,bezier_curve] = makeBezierScatteringRegion(Nx_tot, Ny_tot, Nz_tot, dx, dy, dz, P0, P1, P2, physical_radius,scattering_map)
    % This function generates a scattering region around a 3D Bézier curve
    % Nx_tot, Ny_tot, Nz_tot: Number of grid points in x, y, and z directions
    % dx, dy, dz: Grid spacing in x, y, and z directions (meters)
    % P0, P1, P2: Control points of the 3D Bézier curve (in grid points)
    % physical_radius: Physical radius (in meters) of the scattering region around the curve
    
    % % Position of the Bézier curve in the center of the grid
    % center_x = round(Nx_tot / 2);
    % center_y = round(Ny_tot / 2);
    % center_z = round(Nz_tot / 2);

    % Define the parameter t for the Bézier curve (from 0 to 1)
    t = linspace(0, 1, 100);  % 100 points along the Bézier curve

    % Generate the Bézier curve in 3D using the formula
    bezier_curve_x = (1 - t).^2 * P0(1) + 2 * (1 - t) .* t * P1(1) + t.^2 * P2(1);
    bezier_curve_y = (1 - t).^2 * P0(2) + 2 * (1 - t) .* t * P1(2) + t.^2 * P2(2);
    bezier_curve_z = (1 - t).^2 * P0(3) + 2 * (1 - t) .* t * P1(3) + t.^2 * P2(3);

    % Combine x, y, and z coordinates into the Bézier curve array
    bezier_curve = [bezier_curve_x', bezier_curve_y', bezier_curve_z'];

    % % Visualization of the Bézier curve in 3D to check its shape
    % figure;
    % plot3(bezier_curve(:, 1), bezier_curve(:, 2), bezier_curve(:, 3), '-o');
    % axis equal;
    % xlabel('X Position [grid points]');
    % ylabel('Y Position [grid points]');
    % zlabel('Z Position [grid points]');
    % title('3D Bézier Curve');
    % grid on;

    % % Initialize the scattering map
    % scattering_map = zeros(Nx_tot, Ny_tot, Nz_tot); 

    % Convert physical radius to grid points
    radius = round(physical_radius / dx);  % Convert to grid points based on grid spacing

    % Iterate over the grid points and check if the point lies within the radius of the Bézier curve
    for i = 1:Nx_tot
        for j = 1:Ny_tot
            for k = 1:Nz_tot
                % Calculate the distance from each point on the Bézier curve
                min_distance = inf;  % Initialize with a very large distance

                % Check the distance from each point on the Bézier curve
                for n = 1:length(t)
                    distance = sqrt((bezier_curve(n, 1) - i)^2 + ...
                                    (bezier_curve(n, 2) - j)^2 + ...
                                    (bezier_curve(n, 3) - k)^2);
                    % Update min_distance
                    if distance < min_distance
                        min_distance = distance;
                    end
                end

                % If the point is within the defined radius of the Bézier curve, mark it as part of the scattering region
                if min_distance <= radius
                    scattering_map(i, j, k) = 1;
                end
            end
        end
    end
end
