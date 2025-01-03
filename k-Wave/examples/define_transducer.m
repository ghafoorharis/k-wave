function transducer = define_transducer(kgrid, c0, input_signal,sc,Ny,Nz)


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
    % number_active_elements = 32;
    transducer.active_elements = ones(transducer.number_elements, 1);
    
    % append input signal used to drive the transducer
    transducer.input_signal = input_signal;
    
    % create the transducer using the defined settings
    transducer = kWaveTransducer(kgrid, transducer);
    
    % print out transducer properties
    transducer.properties;

end
