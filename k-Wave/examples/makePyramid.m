function pyramid = makePyramid(Nx, Ny, Nz, x_pos, y_pos, z_pos, size)
    % Initialize the pyramid matrix as a 3D array of zeros
    pyramid = zeros(Nx, Ny, Nz);
    
    % Calculate half of the pyramid base size
    half_base = round(size / 2);
    
    % Loop through each layer of the pyramid
    for z = 0:half_base
        % Calculate the size of the square at the current layer
        current_size = half_base - z;
        
        % Define the boundaries of the square at the current layer
        x_min = max(x_pos - current_size, 1);
        x_max = min(x_pos + current_size, Nx);
        y_min = max(y_pos - current_size, 1);
        y_max = min(y_pos + current_size, Ny);
        z_layer = min(z_pos + z, Nz);  % Move upwards in z-direction
        
        % Fill in the square at the current layer
        pyramid(x_min:x_max, y_min:y_max, z_layer) = 1;
    end
end
