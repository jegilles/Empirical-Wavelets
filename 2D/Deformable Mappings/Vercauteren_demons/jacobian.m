function det_J = jacobian(sx,sy)

    % Gradients
    [gx_y,gx_x] = gradient(sx);
    [gy_y,gy_x] = gradient(sy);
    
    % Add identity
    gx_x = gx_x + 1;  % zero displacement should yield a transformation T = Identity (points keep their positions)
    gy_y = gy_y + 1;  % adding identity matrix here
    
    % Determinant
    det_J = gx_x.*gy_y - ...
            gy_x.*gx_y;
end