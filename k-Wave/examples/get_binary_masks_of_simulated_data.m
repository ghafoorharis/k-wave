% Load the circle data from the specified file
circle_data = load("C:\Users\CMME3\Downloads\k-wave-toolbox-version-1.4\" + ...
    "k-Wave\examples\B_Mode_Simulation_Dataset_New\experiment_485_results\" + ...
    "scan_lines_experiment_485.mat");

% Extract the sound_speed_map from the loaded data
sound_speed_map = circle_data.scattering_circle;

% Select the slice you want to plot (for example, the 25th slice along the z-dimension)
slice_index = 25;

% Plot the slice using imagesc for better scaling
figure;
imagesc(sound_speed_map(:, :, slice_index));
axis image;  % Maintain aspect ratio
colormap("gray");  % Use jet colormap for better visualization
colorbar;  % Display a colorbar to show the scale of values
title(['Sound Speed Map Slice at Index ', num2str(slice_index)]);
xlabel('X (grid points)');
ylabel('Y (grid points)');
