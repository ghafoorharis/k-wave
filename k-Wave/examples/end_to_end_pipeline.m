num_experiments = 3;  % Number of experiments
domain_size = 40e-3;  % Size of the domain (40 mm)
radii_range = [2e-3, 4e-3, 6e-3, 8e-3, 10e-3];  % Radii of the phantom

for exp_num = 1:num_experiments
    for radius_idx = 1:length(radii_range)
        % Adjust parameters as needed
        radius = radii_range(radius_idx);  % Get the radius for this experiment

        % Define sensible bounds for x_shift and y_shift
        % Shifts should ensure the phantom stays within the domain
        max_shift = domain_size - radius;  % Max possible shift ensuring the phantom stays in bounds
        min_shift = radius;                % Min shift to avoid phantom going out of bounds

        % Generate random x_shift and y_shift within these bounds
        x_shift = (max_shift - min_shift) * rand + min_shift;  % Random value between min_shift and max_shift
        y_shift = (max_shift - min_shift) * rand + min_shift;  % Random value between min_shift and max_shift
        
        background_map_mean = 1;        % Modify if you want to test different noise levels
        background_map_std = 0.008;     % Modify if you want to test different noise levels
        scatter_map_mean = 0;           % Modify to adjust scattering behavior
        
        % Call the simulation function and get results
        [scan_lines, sound_speed_map, density_map, scattering_circle, B_mode_image, harmonic_image,horz_axis,transducer] = ...
            simulate_ultrasound_experiment(exp_num, x_shift, y_shift, radius, background_map_mean, background_map_std, scatter_map_mean);
        
        % Create folder for saving results
        folder_name = sprintf('experiment_%d_radius_%.0fmm_xshift_%.0f_yshift_%.0f_results', exp_num, radius*1e3, x_shift*1e3, y_shift*1e3);
        if ~exist(folder_name, 'dir')
            mkdir(folder_name);
        end

        % Save phantom data
        save(fullfile(folder_name, sprintf('phantom_data_exp_%d.mat', exp_num)), 'sound_speed_map', 'density_map', 'scattering_circle');

        % Save scan lines data
        save(fullfile(folder_name, sprintf('scan_lines_experiment_%d.mat', exp_num)), 'scan_lines');

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
        
        % Save all variables in a .mat file
        save(fullfile(folder_name, sprintf('experiment_%d_all_data.mat', exp_num)), 'scan_lines', 'sound_speed_map', 'density_map', 'scattering_circle', 'B_mode_image', 'harmonic_image');
        % Close figures to save memory
        close all;
        clear transducer;
        clear kgrid;
        clear scan_line_example;
        clear scan_lines;
    end
end
