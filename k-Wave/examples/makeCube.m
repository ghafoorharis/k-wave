function cube = makeCube(Nx, Ny, Nz, x_pos, y_pos, z_pos, size)
    cube = zeros(Nx, Ny, Nz);
    cube(x_pos-size:x_pos+size, y_pos-size:y_pos+size, z_pos-size:z_pos+size) = 1;
end
