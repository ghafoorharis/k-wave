function scan_lines = run_simulation(kgrid, medium, transducer, ...
    sound_speed_map, density_map, input_args, folder_name, ...
    exp_num,Ny_tot,scattering_circle, ...
    number_scan_lines)

    scan_lines = zeros(number_scan_lines, kgrid.Nt);
    
    medium_position = 1;
    for scan_line_index = 1:number_scan_lines
        medium.sound_speed = sound_speed_map(:, medium_position:medium_position + Ny_tot - 1, :);
        medium.density = density_map(:, medium_position:medium_position + Ny_tot - 1, :);
        sensor_data = kspaceFirstOrder3D(kgrid, medium, transducer, transducer, input_args{:});
        scan_lines(scan_line_index, :) = transducer.scan_line(sensor_data);
        medium_position = medium_position + transducer.element_width;
    end
    
    save(fullfile(folder_name, sprintf('scan_lines_experiment_%d.mat', exp_num)), 'scan_lines', 'sound_speed_map', 'density_map', 'scattering_circle');
end
