function bezier_points = generateBezier(control_points, t_values)
    % Generate points on a Bézier curve
    % control_points: Nx2 array of (x, y) control points
    % t_values: parameter values for Bézier curve [0, 1]
    
    n = size(control_points, 1) - 1; % Degree of the curve
    bezier_points = zeros(length(t_values), 2); % Initialize points on the curve
    
    for t_idx = 1:length(t_values)
        t = t_values(t_idx);
        bezier_points(t_idx, :) = [0, 0];
        for i = 0:n
            binomial_coeff = nchoosek(n, i);
            bezier_points(t_idx, :) = bezier_points(t_idx, :) + ...
                binomial_coeff * (1 - t)^(n - i) * t^i * control_points(i + 1, :);
        end
    end
end
% Define control points for the Bézier curve
control_points = [15, 30; 20, 35; 25, 25; 30, 30]; % Example control points (x, y)

% Generate curve points
t_values = linspace(0, 1, 100); % Parameter t [0, 1]
bezier_curve = generateBezier(control_points, t_values);
% Initialize a blank grid
bezier_mask = zeros(Nx_tot, Ny_tot, Nz_tot);

% Rasterize Bézier curve (map curve points to grid)
for i = 1:size(bezier_curve, 1)
    x_idx = round(bezier_curve(i, 1) / dx);
    y_idx = round(bezier_curve(i, 2) / dy);
    if x_idx > 0 && x_idx <= Nx_tot && y_idx > 0 && y_idx <= Ny_tot
        bezier_mask(x_idx, y_idx, round(Nz_tot / 2)) = 1;
    end
end

% Assign properties to the Bézier-defined region
sound_speed_map(bezier_mask == 1) = scattering_c0(bezier_mask == 1);
density_map(bezier_mask == 1) = scattering_rho0(bezier_mask == 1);
% Plot the Bézier curve
figure;
plot(bezier_curve(:, 1), bezier_curve(:, 2), 'r-', 'LineWidth', 2);
hold on;
plot(control_points(:, 1), control_points(:, 2), 'bo-', 'LineWidth', 1);
title('Bézier Curve and Control Points');
xlabel('x [m]');
ylabel('y [m]');
grid on;
axis equal;

% Visualize rasterized mask
figure;
imagesc(squeeze(bezier_mask(:, :, round(Nz_tot / 2))));
colormap('gray');
title('Rasterized Bézier Curve Region');
xlabel('x');
ylabel('y');
axis equal;
