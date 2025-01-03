clear;

base_folder = 'US-Bezier-Simulation-V0';
% base_folder = "test";
if ~exist(base_folder, 'dir')
    mkdir(base_folder);  % Create the base folder if it doesn't exist
end
% Simulation settings
DATA_CAST = 'gpuArray-single'; % Set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION = true; % Set to true to run the simulation

% Define shifting values for x and y positions of the phantom
% x_shift = [0,2,-2,3,-3,4,-4,6,-6].* 1e-3; extreme cases
% % y_shift = [0,2,-2,3,-3,-4,4,6,-6] .* 1e-3;extreme cases
% x_shift = [0,3,-3,6,-6].* 1e-3;
% y_shift = [0,3,-3,-6,-6] .* 1e-3;
x_shift = linspace(-20,20,10);
y_shift = linspace(-20,20,10);
bg_mean_range = [1];
% bg_std_range = [0,0.002,0.006,0.012]; extreme cases 
bg_std_range = [0.002,0.006,0.012]; 
scattering_mean_c1_range = [25,45,125];
% scattering_std_c2_range = [2.5,75,255]; % Uncomment for more samples
scattering_std_c2_range = [2.5,255];
scattering_divider_c3_range = [1.5];
% generate the possible radii of phantoms
% range_radii = [3,5,7,9,10,11] .* 1e-3; extreme cases
range_radii = [3,6,11] .* 1e-3;
% range_radii = [1,3,5,7,9,11] * 1e-03;  % Different radii of the phantom
num_pairings = length(x_shift) * length(y_shift) * length(bg_mean_range) * length(bg_std_range) * length(scattering_mean_c1_range) * length(scattering_std_c2_range) *  length(scattering_divider_c3_range) * length(range_radii);
% Generate all possible combinations using ndgrid
[x_grid, y_grid, bg_mean_grid, bg_std_grid, ...
 scattering_mean_grid, scattering_std_grid, scattering_div_grid, radii_grid] = ndgrid(...
    x_shift, y_shift, bg_mean_range, bg_std_range, ...
    scattering_mean_c1_range, scattering_std_c2_range, ...
    scattering_divider_c3_range, range_radii);

% Reshape the grids into column vectors and combine them into one matrix
pairings = [x_grid(:), y_grid(:), bg_mean_grid(:), bg_std_grid(:), ...
    scattering_mean_grid(:), scattering_std_grid(:), scattering_div_grid(:), radii_grid(:)];

% % Generate all unique pairings and shuffle them
% [all_x_shifts, all_y_shifts] = ndgrid(unique(x_shift), unique(y_shift));
% pairings = [all_x_shifts(:), all_y_shifts(:)];  % Flatten the grids and combine

% save('pairings.mat', 'pairings');

% Define constants
% Grid Size Reduction Constant
sc = 2;
% Medium Properties
c0 = 1540;                      % [m/s]
rho0 = 1000;                    % [kg/m^3]
% Input Signal Constants
source_strength = 1e6; % [Pa]
tone_burst_freq = 1.5e6 / sc; % [Hz]
tone_burst_cycles = 4;
% Domain Properties
pml_x_size = 20/sc;
pml_y_size = 10/sc;
pml_z_size = 10/sc;
Nx = 256/sc - 2 * pml_x_size; % [grid points]
Ny = 128/sc - 2 * pml_y_size; % [grid points]
Nz = 128/sc - 2 * pml_z_size; % [grid points]
x = 40e-3; % [m]
dx = x / Nx; % [m]
dy = dx;
dz = dx;
% Create kwave grid
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);
% Define sampling frequency of the ultrasound wavei.e 1/kgrid.dt
t_end = (Nx * dx) * 2.2 / c0; % [s] with sound speed c0 = 1540 m/s
kgrid.makeTime(c0, [], t_end); % Assuming c0 = 1540 m/s for medium
% define medium properties
c0 = 1540; % [m/s]
rho0 = 1000; % [kg/m^3]
medium.alpha_coeff = 0.75; % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;
medium.BonA = 6;
% =========================================================================
% DEFINE THE INPUT SIGNAL
% =========================================================================
input_signal = toneBurst(1 / kgrid.dt, tone_burst_freq, tone_burst_cycles);
input_signal = (source_strength / (c0 * rho0)) * input_signal; % Scale signal

% Define the transducer here
number_scan_lines = 96/sc;

% Define it in the loop
% Define phantom properties
% background_map_mean = 1;
% background_map_std = 0.008;
center_X = 24e-3;
center_Y = 27e-3;
% constant_scattering_region_add = 25;
% constant_scattering_region_multiplier = 75;
% constant_scattering_region_divider = 1.5;
% Experimentation with variable parametets for generating different
% phantoms
% set the input settings
input_args = {...
    'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
    'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};
exp_num = 1;
starting_experiment_index = 1;
% pairings = pairings(exp_num,:);
for params = starting_experiment_index:length(pairings)
    % params = 2111;
    % Get the values from the iterator
    row_params = pairings(params,:);
    x_shift = row_params(1);
    y_shift = row_params(2);
    background_map_mean = row_params(3);
    background_map_std = row_params(4);
    scattering_mean_c1 = row_params(5);
    scattering_std_c2 = row_params(6);
    scattering_divider_c3 = row_params(7);
    phantom_radius = row_params(8);
    % transform the parameters related to the simulation
    % x_pos = center_X + x_shift;
    % y_pos = center_Y + y_shift;
    radius = phantom_radius;
    % SAVE PHANTOM IMAGE AND DATA
    results_path = fullfile(base_folder,sprintf('experiment_%d_results', exp_num));
    if ~exist(results_path, 'dir')
        mkdir(results_path);
    end
    % =========================================================================
    % DEFINE THE PHANTOM WITH SHIFTING CENTER
    % =========================================================================    
    clear transducer;
    transducer = define_transducer(kgrid, c0, input_signal,sc,Ny,Nz);
    % Define image size to scan across
    Nx_tot = Nx;
    Ny_tot = Ny + number_scan_lines * transducer.element_width;
    Nz_tot = Nz;
    % [sound_speed_map, density_map, scattering_circle] = define_phantom( ...
    %                         c0, rho0, ...
    %                         Nx_tot, Ny_tot, Nz_tot, dx, ...
    %                         background_map_mean, ...
    %                         background_map_std, ...
    %                         radius, ...
    %                         x_pos, ...
    %                         y_pos, ...
    %                         scattering_mean_c1, ...
    %                         scattering_std_c2, ...
    %                         scattering_divider_c3);
    [sound_speed_map, density_map, scattering_bezier_map] = define_phantom_bezier( ...
                        c0, rho0, ...
                        Nx_tot, Ny_tot, Nz_tot, dx, ...
                        dy,dz,...
                        background_map_mean, ...
                        background_map_std, ...
                        radius, ...
                        x_shift, ...
                        y_shift, ...
                        scattering_mean_c1, ...
                        scattering_std_c2, ...
                        scattering_divider_c3);
    %=======================================================================
    % RUN THE SIMULATION
    % ===================================================================
    if RUN_SIMULATION
        % Preallocate space for scan lines
        scan_lines = zeros(number_scan_lines, kgrid.Nt);
        
        % Initialize the medium position
        medium_position = 1;
        
        % loop through the scan lines
        for scan_line_index = 1:number_scan_lines
            disp(['Computing scan line ' num2str(scan_line_index) ' of ' num2str(number_scan_lines)]);
            % load the current section of the medium
            medium.sound_speed = sound_speed_map(:, medium_position:medium_position + Ny - 1, :);
            medium.density = density_map(:, medium_position:medium_position + Ny - 1, :);
            
            % run the simulation
            sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});
    
            % extract the scan line from the sensor data
            scan_lines(scan_line_index, :) = transducer.scan_line(sensor_data);
            
            % update medium position
            medium_position = medium_position + transducer.element_width;
    
        end
        % Save the scan lines and phantom data
        save(fullfile(results_path, sprintf('scan_lines_experiment_%d.mat', ...
            exp_num)), 'scan_lines', 'sound_speed_map', ...
            'density_map', 'scattering_bezier_map');
    else
        % disp('s')
        load_results_path = fullfile(base_folder, sprintf('experiment_%d_results', exp_num));
        % Load scan lines from file if the simulation is turned off
        loaded_data = load(fullfile(load_results_path, sprintf('scan_lines_experiment_%d.mat', exp_num)), 'scan_lines');
        % break;
        scan_lines = loaded_data.scan_lines;
    end
    
    % =========================================================================
    % PROCESS THE RESULTS
    % =========================================================================
    [scan_line_example,t0,r,scan_lines_fund,scan_lines_harm] = process_results(kgrid, input_signal, ...
            medium,scan_lines,tone_burst_freq, ...
            c0,number_scan_lines);

    % =========================================================================
    % VISUALIZE AND SAVE RESULTS
    % =========================================================================
    B_mode_image = scan_lines_fund;
    harmonic_image = scan_lines_harm;
    % Save the b_mode and harmonic image
    save(fullfile(results_path, sprintf('images_data_%d.mat', ...
        exp_num)), 'scan_lines_fund','scan_lines_harm');
    scale_factor = 2;
    % Calculate horizontal axis and depth axis for the images
    % r = c0 * ((1:length(kgrid.t_array)) * kgrid.dt - t0) / 2; % Depth axis for B-mode and Harmonic images
    
    % % Visualize and save the B-mode and Harmonic images
    % visualize_and_save_results(B_mode_image.', harmonic_image.', horz_axis, ...
    %     results_path, exp_num, r, ...
    %     sound_speed_map, dx, dy, Ny, Nz, ...
    %     Nx_tot, transducer, number_scan_lines);
    % plot the medium, truncated to the field of view
    figure;
    horz_axis_phantom = (0:number_scan_lines * transducer.element_width - 1) * dy * 1e3;
    imagesc(horz_axis_phantom, (0:Nx_tot-1) * dx * 1e3, sound_speed_map(:, 1 + Ny/2:end - Ny/2, Nz/2));
    axis image;
    colormap(gray);
    set(gca, 'YLim', [5, 40]);
    title('Scattering Phantom');
    xlabel('Horizontal Position [mm]');
    ylabel('Depth [mm]');
    % Save the phantom image
    saveas(gcf, fullfile(results_path, sprintf('phantom_image_exp_%d.png', exp_num)));
    
    % close(gcf); % Close figure to save memory
    
    % % plot the processing steps
    % figure;
    % stackedPlot(kgrid.t_array * 1e6, {'1. Beamformed Signal', '2. Time Gain Compensation', '3. Frequency Filtering', '4. Envelope Detection', '5. Log Compression'}, scan_line_example);
    % xlabel('Time [\mus]');
    % set(gca, 'XLim', [5, t_end * 1e6]);
    % 
    % plot the processed b-mode ultrasound image
    figure;
    horz_axis = (0:length(scan_lines_fund(:, 1)) - 1) * transducer.element_width * dy / scale_factor * 1e3;
    imagesc(horz_axis, r * 1e3, scan_lines_fund.');
    axis image;
    colormap(gray);
    set(gca, 'YLim', [5, 40]);
    title('B-mode Image');
    xlabel('Horizontal Position [mm]');
    ylabel('Depth [mm]');
    % Save the B-mode image
    saveas(gcf, fullfile(results_path, sprintf('B_mode_image_exp_%d.png', exp_num)));
    % Save the frame into the rendered image
    % Capture the rendered image data from the figure
    frame = getframe(gca);           % Capture the displayed axes with aspect adjustments
    rendered_bmode_image = frame.cdata;  % Extract the image data as an RGB matrix
    
    % Save the rendered image data as a .mat file
    save(fullfile(results_path, sprintf('B_mode_image_exp_%d.mat', exp_num)), 'rendered_bmode_image');
    
    % close(gcf); % Close figure to save memory
    
    % plot the processed harmonic ultrasound image
    figure;
    imagesc(horz_axis, r * 1e3, scan_lines_harm.');
    axis image;
    colormap(gray);
    set(gca, 'YLim', [5, 40]);
    title('Harmonic Image');
    xlabel('Horizontal Position [mm]');
    ylabel('Depth [mm]');
    % % Save the Harmonic image
    saveas(gcf, fullfile(results_path, sprintf('harmonic_image_exp_%d.png', exp_num)));
    
    % % Plot the processed harmonic ultrasound image
    % figure;
    % imagesc(horz_axis, r * 1e3, scan_lines_harm.');
    % axis image;
    % colormap(gray);
    % set(gca, 'YLim', [5, 40]);
    % title('Harmonic Image');
    % xlabel('Horizontal Position [mm]');
    % ylabel('Depth [mm]');
    % % Capture the rendered image data from the figure
    frame = getframe(gca);           % Capture the displayed axes with aspect adjustments
    rendered_harmonic_image = frame.cdata;  % Extract the image data as an RGB matrix
    
    % Save the rendered image data as a .mat file
    save(fullfile(results_path, sprintf('harmonic_image_exp_%d.mat', exp_num)), 'rendered_harmonic_image');

    % Save the image as mat
    % close(gcf); % Close figure to save memory
    % Clear variables to save memory for the next experiment
    clear transducer;
    % clear kgrid;
    % clear scan_line_example;
    % clear scan_lines;
    % close;
    exp_num = exp_num + 1;
    if RUN_SIMULATION
        close all;
    else
        break
    end
    % break;
end
