function visualize_and_save_results(B_mode_image, harmonic_image, ...
    horz_axis, folder_name, exp_num, r, ...
    sound_speed_map, dx, dy, Ny, Nz, ...
    Nx_tot, transducer, number_scan_lines)
    % =========================================================================
    % PLOT THE SCATTERING PHANTOM (SOUND SPEED MAP)
    %=========================================================================
    figure;
    imagesc((0:number_scan_lines * transducer.element_width - 1) * dy * 1e3, ...
        (0:Nx_tot-1) * dx * 1e3, ...
        sound_speed_map(:, 1 + Ny/2:end - Ny/2, Nz/2));
    axis image;
    colormap(gray);
    set(gca, 'YLim', [5, 40]); % Adjust Y-axis limits if necessary
    title('Scattering Phantom');
    xlabel('Horizontal Position [mm]');
    ylabel('Depth [mm]');

    % Save the phantom image
    saveas(gcf, fullfile(folder_name, sprintf('phantom_image_exp_%d.png', exp_num)));
    close(gcf); % Close figure to save memory

    % =========================================================================
    % VISUALIZE AND SAVE B-MODE IMAGE
    % =========================================================================
    figure;
    imagesc(horz_axis, r * 1e3, B_mode_image); % Plot B-mode image
    colormap(gray);
    axis image;
    title('B-mode Image');
    xlabel('Horizontal Position [mm]');
    ylabel('Depth [mm]');
    
    % Save the B-mode image
    saveas(gcf, fullfile(folder_name, sprintf('B_mode_image_exp_%d.png', exp_num)));
    close(gcf); % Close figure to save memory

    % =========================================================================
    % VISUALIZE AND SAVE HARMONIC IMAGE
    % =========================================================================
    figure;
    imagesc(horz_axis, r * 1e3, harmonic_image); % Plot harmonic image
    colormap(gray);
    axis image;
    title('Harmonic Image');
    xlabel('Horizontal Position [mm]');
    ylabel('Depth [mm]');
    
    % Save the Harmonic image
    saveas(gcf, fullfile(folder_name, sprintf('harmonic_image_exp_%d.png', exp_num)));
    close(gcf); % Close figure to save memory
end
