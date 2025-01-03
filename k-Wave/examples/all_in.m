clearvars;

% Simulation settings
DATA_CAST       = 'gpuArray-single'; % set to 'single' or 'gpuArray-single' to speed up computations
RUN_SIMULATION  = false; % set to true to run simulation

% Number of experiments
num_experiments = 5;

% Define shifting values for x and y positions of the phantom
x_shift_range = linspace(0, 10e-3, num_experiments);  % [m]
y_shift_range = linspace(0, 5e-3, num_experiments);   % [m]

for exp_num = 1:num_experiments
    
    % =========================================================================
    % DEFINE THE K-WAVE GRID
    % =========================================================================

    % PML settings
    sc = 2;
    pml_x_size = 20/sc;                % [grid points]
    pml_y_size = 10/sc;                % [grid points]
    pml_z_size = 10/sc;                % [grid points]

    % Total grid points excluding PML
    Nx = 256/sc - 2 * pml_x_size;      % [grid points]
    Ny = 128/sc - 2 * pml_y_size;      % [grid points]
    Nz = 128/sc - 2 * pml_z_size;      % [grid points]

    % Grid spacing and dimensions
    x = 40e-3;                      % [m]
    dx = x / Nx;                    % [m]
    dy = dx;
    dz = dx;

    % Create k-wave grid
    kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

    % =========================================================================
    % DEFINE THE MEDIUM PARAMETERS
    % =========================================================================
    c0 = 1540;                      % [m/s]
    rho0 = 1000;                    % [kg/m^3]
    medium.alpha_coeff = 0.75;      % [dB/(MHz^y cm)]
    medium.alpha_power = 1.5;
    medium.BonA = 6;

    % Create time array
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
    % DEFINE THE PHANTOM WITH SHIFTING CENTER
    % =========================================================================

    % Define image size to scan across
    number_scan_lines = 96/sc;
    Nx_tot = Nx;
    Ny_tot = Ny + number_scan_lines * transducer.element_width;
    Nz_tot = Nz;

    % Background scatterer
    background_map_mean = 1;
    background_map_std = 0.008;
    background_map = background_map_mean + background_map_std * randn([Nx_tot, Ny_tot, Nz_tot]);

    % Scattering map for different region
    scattering_map = randn([Nx_tot, Ny_tot, Nz_tot]);
    scattering_c0 = c0 + 25 + 75 * scattering_map;
    scattering_c0(scattering_c0 > 1600) = 1600;
    scattering_c0(scattering_c0 < 1400) = 1400;
    scattering_rho0 = scattering_c0 / 1.5;

    % Uniform medium
    sound_speed_map = c0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
    density_map = rho0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;

    % Shifted circle (phantom)
    radius = 6e-3;  % [m]
    x_pos = 25e-3 + x_shift_range(exp_num);  % [m]
    y_pos = 30e-3 + y_shift_range(exp_num);  % [m]
    scattering_circle = makeBall(Nx_tot, Ny_tot, Nz_tot, round(x_pos/dx), round(y_pos/dx), Nz_tot/2, round(radius/dx));
    sound_speed_map(scattering_circle == 1) = scattering_c0(scattering_circle == 1);
    density_map(scattering_circle == 1) = scattering_rho0(scattering_circle == 1);

    % =========================================================================
    % SAVE PHANTOM IMAGE AND DATA
    % =========================================================================

    % Create folder for saving data
    folder_name = sprintf('experiment_%d_results', exp_num);
    if ~exist(folder_name, 'dir')
        mkdir(folder_name);
    end

    % Plot the phantom (sound speed map)
    figure;
    imagesc((0:Nx_tot - 1) * dx * 1e3, (0:Ny_tot - 1) * dy * 1e3, sound_speed_map(:, :, Nz_tot / 2));
    axis image;
    colormap(gray);
    title('Phantom Sound Speed Map');
    xlabel('Horizontal Position [mm]');
    ylabel('Depth [mm]');
    saveas(gcf, fullfile(folder_name, sprintf('phantom_image_exp_%d.png', exp_num)));
    close(gcf);  % Close the figure to save memory

    % Save the phantom data (sound speed and density)
    save(fullfile(folder_name, sprintf('phantom_data_exp_%d.mat', exp_num)), 'sound_speed_map', 'density_map', 'scattering_circle');

    % =========================================================================
    % RUN THE SIMULATION
    % =========================================================================

    % Preallocate storage for scan lines
    scan_lines = zeros(number_scan_lines, kgrid.Nt);

    % Input arguments for simulation
    input_args = {
        'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], ...
        'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

    if RUN_SIMULATION
        % Medium position
        medium_position = 1;

        % Loop through scan lines
        for scan_line_index = 1:number_scan_lines
            % Display status
            disp(['Computing scan line ' num2str(scan_line_index) ' of ' num2str(number_scan_lines)]);

            % Load section of the medium
            medium.sound_speed = sound_speed_map(:, medium_position:medium_position + Ny - 1, :);
            medium.density = density_map(:, medium_position:medium_position + Ny - 1, :);

            % Run simulation
            sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});

            % Extract scan line
            scan_lines(scan_line_index, :) = transducer.scan_line(sensor_data);

            % Update medium position
            medium_position = medium_position + transducer.element_width;
        end

        % Create folder for saving data
        folder_name = sprintf('experiment_%d_results', exp_num);
        if ~exist(folder_name, 'dir')
            mkdir(folder_name);
        end

        % Save scan lines and phantom data
        save(fullfile(folder_name, sprintf('scan_lines_experiment_%d.mat', exp_num)), 'scan_lines', 'sound_speed_map', 'density_map', 'scattering_circle');

    else
        % Load scan lines from file
        load(sprintf('experiment_%d_results/scan_lines_experiment_%d.mat', exp_num));
    end

    % =========================================================================
    % PROCESS THE RESULTS
    % =========================================================================

    % -----------------------------
    % Remove Input Signal
    % -----------------------------
    scan_line_win = getWin(kgrid.Nt * 2, 'Tukey', 'Param', 0.05).';
    scan_line_win = [zeros(1, length(input_signal) * 2), scan_line_win(1:end/2 - length(input_signal) * 2)];
    scan_lines = bsxfun(@times, scan_line_win, scan_lines);
    scan_line_example(1, :) = scan_lines(end/2, :);

    % -----------------------------
    % Time Gain Compensation
    % -----------------------------
    t0 = length(input_signal) * kgrid.dt / 2;
    r = c0 * ( (1:length(kgrid.t_array)) * kgrid.dt - t0 ) / 2;    % [m]
    tgc_alpha_db_cm = medium.alpha_coeff * (tone_burst_freq * 1e-6)^medium.alpha_power;
    tgc_alpha_np_m = tgc_alpha_db_cm / 8.686 * 100;
    tgc = exp(tgc_alpha_np_m * 2 * r);
    scan_lines = bsxfun(@times, tgc, scan_lines);
    scan_line_example(2, :) = scan_lines(end/2, :);

    % -----------------------------
    % Frequency Filtering
    % -----------------------------
    scan_lines_fund = gaussianFilter(scan_lines, 1/kgrid.dt, tone_burst_freq, 100, true);
    scan_lines_harm = gaussianFilter(scan_lines, 1/kgrid.dt, 2 * tone_burst_freq, 30, true);
    scan_line_example(3, :) = scan_lines_fund(end/2, :);

    % -----------------------------
    % Envelope Detection
    % -----------------------------
    scan_lines_fund = envelopeDetection(scan_lines_fund);
    scan_lines_harm = envelopeDetection(scan_lines_harm);
    scan_line_example(4, :) = scan_lines_fund(end/2, :);

    % -----------------------------
    % Log Compression
    % -----------------------------
    compression_ratio = 30;
    scan_lines_fund = logCompression(scan_lines_fund, compression_ratio, true);
    scan_lines_harm = logCompression(scan_lines_harm, compression_ratio, true);
    scan_line_example(5, :) = scan_lines_fund(end/2, :);

    % -----------------------------
    % Scan Conversion
    % -----------------------------
    scale_factor = 2;
    scan_lines_fund = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_fund, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');
    scan_lines_harm = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_harm, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');

    % =========================================================================
    % VISUALISATION AND SAVING FIGURES
    % =========================================================================

    % B-mode Image
    horz_axis = (0:length(scan_lines_fund(:, 1)) - 1) * transducer.element_width * dy / scale_factor * 1e3;
    figure;
    imagesc(horz_axis, r * 1e3, scan_lines_fund.');
    axis image;
    colormap(gray);
    set(gca, 'YLim', [5, 40]);
    title('B-mode Image');
    xlabel('Horizontal Position [mm]');
    ylabel('Depth [mm]');
    saveas(gcf, fullfile(folder_name, sprintf('B-mode_image_exp_%d.png', exp_num)));

    % Harmonic Image
    figure;
    imagesc(horz_axis, r * 1e3, scan_lines_harm.');
    axis image;
    colormap(gray);
    set(gca, 'YLim', [5, 40]);
    title('Harmonic Image');
    xlabel('Horizontal Position [mm]');
    ylabel('Depth [mm]');
    saveas(gcf, fullfile(folder_name, sprintf('Harmonic_image_exp_%d.png', exp_num)));

    % Save all variables to a .mat file
    save(fullfile(folder_name, sprintf('experiment_%d_all_data.mat', exp_num)));

    % Close figures to save memory
    close all;
    clear transducer;
    clear kgrid;
    clear scan_line_example;
    clear scan_lines;

    % break;
end
