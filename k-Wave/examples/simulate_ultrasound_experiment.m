function [scan_lines, sound_speed_map, density_map, scattering_circle, B_mode_image, harmonic_image,horz_axis,transducer] = ...
    simulate_ultrasound_experiment(exp_num, x_shift, y_shift, radius, background_map_mean, background_map_std, scatter_map_mean)

    % =========================================================================
    % SIMULATION SETTINGS AND GRID DEFINITION
    % =========================================================================
    DATA_CAST = 'gpuArray-single'; % Set to 'single' or 'gpuArray-single' to speed up computations
    RUN_SIMULATION = true;         % Set to true to run simulation

    % Define the k-wave grid
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
    source_strength = 1e6;          % [Pa]
    tone_burst_freq = 1.5e6/sc;     % [Hz]
    tone_burst_cycles = 4;
    input_signal = toneBurst(1/kgrid.dt, tone_burst_freq, tone_burst_cycles);
    input_signal = (source_strength ./ (c0 * rho0)) .* input_signal;

    % =========================================================================
    % DEFINE THE ULTRASOUND TRANSDUCER
    % =========================================================================
    clear transducer;
    transducer.number_elements = 32/sc;  	% Total number of transducer elements
    transducer.element_width = 2;          % Width of each element [grid points]
    transducer.element_length = 24/sc;    	% Length of each element [grid points]
    transducer.element_spacing = 0;        % Spacing between the elements [grid points]
    transducer.radius = inf;               % Radius of curvature of the transducer [m]

    transducer_width = transducer.number_elements * transducer.element_width ...
        + (transducer.number_elements - 1) * transducer.element_spacing;

    % Position the transducer in the middle of the computational grid
    transducer.position = round([1, Ny/2 - transducer_width/2, Nz/2 - transducer.element_length/2]);
    transducer.sound_speed = c0;
    transducer.focus_distance = 20e-3;
    transducer.elevation_focus_distance = 19e-3;
    transducer.steering_angle = 0;
    transducer.transmit_apodization = 'Hanning';
    transducer.receive_apodization = 'Rectangular';
    transducer.active_elements = ones(transducer.number_elements, 1);
    transducer.input_signal = input_signal;
    transducer = kWaveTransducer(kgrid, transducer);

    % =========================================================================
    % DEFINE THE PHANTOM WITH SHIFTING CENTER AND CHANGING RADIUS
    % =========================================================================

    % Define image size to scan across
    number_scan_lines = 96/sc;
    Nx_tot = Nx;
    Ny_tot = Ny + number_scan_lines * transducer.element_width;
    Nz_tot = Nz;

    % Background scatterer with specified mean and standard deviation
    background_map = background_map_mean + background_map_std * randn([Nx_tot, Ny_tot, Nz_tot]);

    % Scattering map with the specified mean
    scattering_map = scatter_map_mean + randn([Nx_tot, Ny_tot, Nz_tot]);
    scattering_c0 = c0 + 25 + 75 * scattering_map;
    scattering_c0(scattering_c0 > 1600) = 1600;
    scattering_c0(scattering_c0 < 1400) = 1400;
    scattering_rho0 = scattering_c0 / 1.5;

    % Uniform medium
    sound_speed_map = c0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;
    density_map = rho0 * ones(Nx_tot, Ny_tot, Nz_tot) .* background_map;

    % Shifted circle (phantom)
    x_pos = 25e-3 + x_shift;
    y_pos = 30e-3 + y_shift;
    scattering_circle = makeBall(Nx_tot, Ny_tot, Nz_tot, round(x_pos/dx), round(y_pos/dx), Nz_tot/2, round(radius/dx));
    sound_speed_map(scattering_circle == 1) = scattering_c0(scattering_circle == 1);
    density_map(scattering_circle == 1) = scattering_rho0(scattering_circle == 1);

    % =========================================================================
    % RUN THE SIMULATION
    % =========================================================================
    scan_lines = zeros(number_scan_lines, kgrid.Nt);
    input_args = {'PMLInside', false, 'PMLSize', [pml_x_size, pml_y_size, pml_z_size], 'DataCast', DATA_CAST, 'DataRecast', true, 'PlotSim', false};

    if RUN_SIMULATION
        medium_position = 1;
        for scan_line_index = 1:number_scan_lines
            medium.sound_speed = sound_speed_map(:, medium_position:medium_position + Ny - 1, :);
            medium.density = density_map(:, medium_position:medium_position + Ny - 1, :);
            sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});
            scan_lines(scan_line_index, :) = transducer.scan_line(sensor_data);
            medium_position = medium_position + transducer.element_width;
        end
    end

    % =========================================================================
    % PROCESS THE RESULTS (Generate B-mode and Harmonic Images)
    % =========================================================================
    scan_line_win = getWin(kgrid.Nt * 2, 'Tukey', 'Param', 0.05).';
    scan_line_win = [zeros(1, length(input_signal) * 2), scan_line_win(1:end/2 - length(input_signal) * 2)];
    scan_lines = bsxfun(@times, scan_line_win, scan_lines);

    % Time gain compensation
    t0 = length(input_signal) * kgrid.dt / 2;
    r = c0 * ( (1:length(kgrid.t_array)) * kgrid.dt - t0 ) / 2;
    tgc_alpha_db_cm = medium.alpha_coeff * (tone_burst_freq * 1e-6)^medium.alpha_power;
    tgc_alpha_np_m = tgc_alpha_db_cm / 8.686 * 100;
    tgc = exp(tgc_alpha_np_m * 2 * r);
    scan_lines = bsxfun(@times, tgc, scan_lines);

    % Frequency filtering
    scan_lines_fund = gaussianFilter(scan_lines, 1/kgrid.dt, tone_burst_freq, 100, true);
    scan_lines_harm = gaussianFilter(scan_lines, 1/kgrid.dt, 2 * tone_burst_freq, 30, true);

    % Envelope detection and log compression
    scan_lines_fund = envelopeDetection(scan_lines_fund);
    scan_lines_harm = envelopeDetection(scan_lines_harm);
    compression_ratio = 30;
    scan_lines_fund = logCompression(scan_lines_fund, compression_ratio, true);
    scan_lines_harm = logCompression(scan_lines_harm, compression_ratio, true);

    % Scan conversion
    scale_factor = 2;
    scan_lines_fund = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_fund, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');
    scan_lines_harm = interp2(1:kgrid.Nt, (1:number_scan_lines).', scan_lines_harm, 1:kgrid.Nt, (1:1/scale_factor:number_scan_lines).');

    % B-mode and Harmonic images
    horz_axis = (0:length(scan_lines_fund(:, 1)) - 1) * transducer.element_width * dy / scale_factor * 1e3;
    B_mode_image = scan_lines_fund.';
    harmonic_image = scan_lines_harm.';
end