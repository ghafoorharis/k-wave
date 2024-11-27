clearvars;
clc;
close all;

% Define the base folder for saving the phantoms
base_folder = 'Generated_Bezier_Phantoms_Final_Pakka';
% base_folder = 'test_bezier';

if ~exist(base_folder, 'dir')
    mkdir(base_folder);  % Create the base folder if it doesn't exist
end
% =========================================================================
% DEFINE THE K-WAVE GRID
% =========================================================================

% set the size of the perfectly matched layer (PML)
sc = 2;
pml_x_size = 20/sc;                % [grid points]
pml_y_size = 10/sc;                % [grid points]
pml_z_size = 10/sc;                % [grid points]

% set total number of grid points not including the PML
Nx = 256/sc - 2 * pml_x_size;      % [grid points]
Ny = 128/sc - 2 * pml_y_size;      % [grid points]
Nz = 128/sc - 2 * pml_z_size;      % [grid points]

% set desired grid size in the x-direction not including the PML
x = 40e-3;                      % [m]

% calculate the spacing between the grid points
dx = x / Nx;                    % [m]
dy = dx;                        % [m]
dz = dx;                        % [m]

% create the k-space grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% =========================================================================
% DEFINE THE MEDIUM PARAMETERS
% =========================================================================

% define the properties of the propagation medium
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;

% create the time array
t_end = (Nx * dx) * 2.2 / c0;   % [s]
kgrid.makeTime(c0, [], t_end);

% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================

% define properties of the input signal
source_strength = 1e6;          % [Pa]
tone_burst_freq = 1.5e6/sc;        % [Hz]
tone_burst_cycles = 4;

% create the input signal using toneBurst 
input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);

% scale the source magnitude by the source_strength divided by the
% impedance (the source is assigned to the particle velocity)
input_signal = (source_strength ./ (c0 * rho0)) .* input_signal;

% =========================================================================
% DEFINE THE ULTRASOUND TRANSDUCER
% =========================================================================

% physical properties of the transducer
transducer.number_elements = 32/sc;  	% total number of transducer elements
transducer.element_width = 2;       % width of each element [grid points]
transducer.element_length = 24/sc;  	% length of each element [grid points]
transducer.element_spacing = 0;  	% spacing (kerf  width) between the elements [grid points]
transducer.radius = inf;            % radius of curvature of the transducer [m]

% calculate the width of the transducer in grid points
transducer_width = transducer.number_elements * transducer.element_width ...
    + (transducer.number_elements - 1) * transducer.element_spacing;

% use this to position the transducer in the middle of the computational grid
transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);
% print(transducer.position)
% properties used to derive the beamforming delays
transducer.sound_speed = c0;                    % sound speed [m/s]
transducer.focus_distance = 20e-3;              % focus distance [m]
transducer.elevation_focus_distance = 19e-3;    % focus distance in the elevation plane [m]
transducer.steering_angle = 0;                  % steering angle [degrees]

% apodization
transducer.transmit_apodization = 'Hanning';    
transducer.receive_apodization = 'Rectangular';

% define the transducer elements that are currently active
number_active_elements = 32;
transducer.active_elements = ones(transducer.number_elements, 1);

% append input signal used to drive the transducer
transducer.input_signal = input_signal;

% create the transducer using the defined settings
transducer = kWaveTransducer(kgrid, transducer);

% print out transducer properties
transducer.properties;

% =========================================================================
% DEFINE THE MEDIUM PROPERTIES
% =========================================================================

% define a large image size to move across
number_scan_lines = 96/sc;
Nx_tot = Nx;
Ny_tot = Ny + number_scan_lines * transducer.element_width;
Nz_tot = Nz;
% Phantom center (fixed)
center_X = Nx_tot / 2;  % Center x-coordinate
center_Y = Ny_tot / 2;  % Center y-coordinate
center_Z = Nz_tot / 2;  % Center z-coordinate
% define a random distribution of scatterers for the medium
background_map_mean = 1;
background_map_std = 0.012;
background_map = background_map_mean + background_map_std * randn([Nx_tot, Ny_tot, Nz_tot]);

% define a random distribution of scatterers for the highly scattering
% region
scattering_map = randn([Nx_tot, Ny_tot, Nz_tot]);
scattering_c0 = c0 + 25 + 75 * scattering_map;
scattering_c0(scattering_c0 > 1600) = 1600;
scattering_c0(scattering_c0 < 1400) = 1400;
scattering_rho0 = scattering_c0 / 1.5;

% Scattering region: initialize uniform medium
sound_speed_map = c0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
density_map = rho0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;

% Define ranges for the parameter grid
% x_shift = [-30,-10,-5, 0,5,9,20];  % X shift for phantom
% y_shift = [-30,-10,-5, 0,6,7,20];  % Y shift for phantom
x_shift = linspace(-20,20,20);
y_shift = linspace(-20,20,20);
bg_mean_range = [1];           % Background mean (set to constant for simplicity)
bg_std_range = [0.001];  % Noise levels
scattering_mean_c1_range = [25];  % Scattering mean
scattering_std_c2_range = [2.5];  % Scattering standard deviation
scattering_divider_c3_range = [1.5];  % Scattering divider
range_radii = [1,3,5,7,9,11] * 1e-03;  % Different radii of the phantom

% Generate all combinations of parameters
[x_grid, y_grid, bg_mean_grid, bg_std_grid, ...
 scattering_mean_grid, scattering_std_grid, scattering_div_grid, radii_grid] = ndgrid(... 
    x_shift, y_shift, bg_mean_range, bg_std_range, ...
    scattering_mean_c1_range, scattering_std_c2_range, ...
    scattering_divider_c3_range, range_radii);

% Reshape the grids into column vectors and combine them into one matrix
param_combinations = [x_grid(:), y_grid(:), bg_mean_grid(:), bg_std_grid(:), ...
    scattering_mean_grid(:), scattering_std_grid(:), scattering_div_grid(:), radii_grid(:)];

% Experiment iteration
exp_num = 1;
% Loop through each combination of parameters
for idx = 1:size(param_combinations, 1)
    % Extract parameters from the combination
    x_shift = param_combinations(idx, 1);
    y_shift = param_combinations(idx, 2);
    background_map_mean = param_combinations(idx, 3);
    background_map_std = param_combinations(idx, 4);
    scattering_mean_c1 = param_combinations(idx, 5);
    scattering_std_c2 = param_combinations(idx, 6);
    scattering_divider_c3 = param_combinations(idx, 7);
    phantom_radius = param_combinations(idx, 8);
    % Random bg_map
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

    % Ensure the radius is within the max allowed limit
    max_radius = 100;  % Maximum radius allowed (e.g., 11 mm)
    phantom_radius = min(phantom_radius, max_radius);
    physical_radius = phantom_radius;
    % Define the control points for the Bézier curve (circular arc)
    P0 = [center_X - x_shift, center_Y - y_shift, center_Z];  % Start point on the circle
    P2 = [center_X, center_Y, center_Z];  % End point on the circle
    P1 = [center_X + x_shift / 2, center_Y + y_shift, center_Z];  % Control point for the arc
    
    % Generate the phantom using Bézier curve (function you defined)
    [scattering_map, bezier_curve] = makeBezierScatteringRegion(Nx_tot, Ny_tot, Nz_tot, dx, dy, dz, P0, P1, P2, ...
        physical_radius, scattering_map);
    % Modify sound speed and density in the scattering region
    sound_speed_map(scattering_map == 1) = scattering_c0(scattering_map==1);  % Increase sound speed in scatterers
    density_map(scattering_map == 1) = scattering_rho0(scattering_map==1);    % Increase density in scatterers

    % Visualization of the phantom
    figure;
    % 2D Cross-section (horizontal slice)
    % subplot(1, 2, 1);
    imagesc(sound_speed_map(:, :, round(Nz_tot / 2)));  % Middle slice
    axis equal;
    colormap(gray);
    title(sprintf('Phantom Slice Exp %d, Radius %.3f mm, Noise %.3f', exp_num, phantom_radius * 1e3, background_map_std));
    xlabel('X [mm]');
    ylabel('Y [mm]');

    % % 3D Visualization (Voxel plot)
    % subplot(1, 2, 2);
    % voxelPlot(scattering_map);
    % title('3D Visualization');

    % Save images and data for each phantom
    results_path = fullfile(base_folder, sprintf('experiment_%d', exp_num));
    if ~exist(results_path, 'dir')
        mkdir(results_path);  % Create folder for this experiment
    end

    % Save the generated phantom image
    saveas(gcf, fullfile(results_path, sprintf('phantom_image_exp_%d_radius_%.3fmm_noise_%.3f.png', exp_num, phantom_radius * 1e3, background_map_std)));
    
    % Save the phantom data
    save(fullfile(results_path, sprintf('phantom_data_exp_%d_radius_%.3fmm_noise_%.3f.mat', exp_num, phantom_radius * 1e3, background_map_std)), 'sound_speed_map', 'bezier_curve');
    
    % Increment experiment number
    exp_num = exp_num + 1;
    close(gcf);  % Close figure to save memory
    close all;
    % break;
end

disp('Phantom generation and visualization complete!');
close all;
clc;
clear;