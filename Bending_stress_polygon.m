function Bending_stress_edge = Bending_stress_polygon(chordDist, r ,rho, Vr, phi, Cldes) 
    load("AerofoilCoord.mat",'Aerofoil_coord')

    polygon = struct('x', Aerofoil_coord(:, 1), 'y', Aerofoil_coord(:, 2));


    % Bending_stress finds the Bending stress dist of the wind turbine due
    % to aerodynamic forces
    %
    % Input arguments:
    %   chord_dist: the blades chord distribution
    %   r: Contains the radial distances of the blade sections corresponding to the chord dist in chord_dist.
    %   rho: Density
    %   CL_des: CL_des
    %   Vr: is the inlet velocity in the rotating (blade) frame of reference. [Lubbock Aero page 19 eq 18]
    %   phi: is the relative inflow angle with respect to the tangential direction.
    %
    % Output:
    %   Bending_stress: Distribution of the centrifugal stress along the
    %   blade.
    

    % Define Constants: 
    B = 2; % number of blades
    
    % Find aerodynamic force:
    p_y = B * 0.5 .* (Vr.^2) .* Cldes .* chordDist .* sin(phi);
    p_z = B * 0.5 .* (Vr.^2) .* Cldes .* chordDist .* cos(phi);
    
    % Find the shear force
    % T_z = integral( p_y(x)) between x to r using the trapezium rule
    % T_y = integral( p_y(x)) between x to r using the trapezium rule
    
    delta_x = r(2) - r(1); % x difference between each trapezium
    
    T_z = zeros(1, length(r)); % Initialize T_z with zeros
    T_z(1) = 8.91; % Initialize the first element of T_z since integral of a point = 0
    for n = 2:length(r)
        x = linspace(r(n-1), r(end), n); % Create equally spaced points from r(n-1) to r(end)
        y = p_z(1:n); % Extract values of p_z up to n
        
        T_z(n) = trapz(x, y); % Apply the trapezoidal rule
    end
    
    
    T_y = zeros(1, length(r)); % Initialize T_y with zeros
    T_y(1) = 4.12; % Initialize the first element of T_y (integral of a point = 0)
    for n = 2:length(r)
        x = linspace(r(n-1), r(end), n); % Create equally spaced points from r(n-1) to r(end)
        y = p_y(1:n); % Extract values of p_y up to n
        
        T_y(n) = trapz(x, y); % Apply the trapezoidal rule
    end
    
    % Find the Bending moments
    % M_y = - integral( T_z(x)) between x to R using the trapezium rule
    % M_z = - integral( T_y(x)) between x to R using the trapezium rule
    
    M_z = 0;
    for n = 1:max(size(r)) - 1
        M_z = (T_y(n)+T_y(n+1)) * 0.5 * delta_x + M_z;  % Sum of each trapezium
    end
    
    M_z = - M_z;

    M_y = 0;
    for n = 1:max(size(r)) - 1
        M_y = (T_z(n)+T_z(n+1)) * 0.5 * delta_x + M_y;  % Sum of each trapezium
    end

    M_y = -M_y;

    % Find Iz, Iy, Izy {see onenote for axis definition}
    % Idealize Aerofoil as a polygon
    % thickness(0.17* chord length)
    
    % Initialize Iz, Iy, Izy
    Iz = zeros(size(r));
    Iy = zeros(size(r));
    Izy = zeros(size(r));
    
    for i = 1:numel(r)
        % Scale polygon coordinates by r(i)
        scaled_polygon = struct('x', r(i) * polygon.x, 'y', r(i) * polygon.y);
        
        % Calculate centroid of the scaled polygon
        A = 0;
        C_x = 0;
        C_y = 0;
    
        for j = 1:numel(scaled_polygon.x)
            k = j + 1;
            if j == numel(scaled_polygon.x)
                k = 1;
            end
            A = A + (scaled_polygon.x(j) * scaled_polygon.y(k) - scaled_polygon.x(k) * scaled_polygon.y(j));
            C_x = C_x + (scaled_polygon.x(j) + scaled_polygon.x(k)) * (scaled_polygon.x(j) * scaled_polygon.y(k) - scaled_polygon.x(k) * scaled_polygon.y(j));
            C_y = C_y + (scaled_polygon.y(j) + scaled_polygon.y(k)) * (scaled_polygon.x(j) * scaled_polygon.y(k) - scaled_polygon.x(k) * scaled_polygon.y(j));
        end
    
        centroid_x = C_x / (3 * A);
        centroid_y = C_y / (3 * A);
    
        % Calculate moments of inertia
        Ix = 0;
        Iy_temp = 0;
        Iz_temp = 0;
        Ixy = 0;
    
        for j = 1:numel(scaled_polygon.x)
            k = j + 1;
            if j == numel(scaled_polygon.x)
                k = 1;
            end
    
            ai = (scaled_polygon.x(k) * scaled_polygon.y(j) - scaled_polygon.x(j) * scaled_polygon.y(k)) / 2;
            Ix = Ix + (scaled_polygon.y(j)^2 + scaled_polygon.y(j) * scaled_polygon.y(k) + scaled_polygon.y(k)^2) * ai;
            Iy_temp = Iy_temp + (scaled_polygon.x(j)^2 + scaled_polygon.x(j) * scaled_polygon.x(k) + scaled_polygon.x(k)^2) * ai;
            Iz_temp = Iz_temp + (scaled_polygon.x(j)^2 + scaled_polygon.x(j) * scaled_polygon.x(k) + scaled_polygon.x(k)^2 + scaled_polygon.y(j)^2 + scaled_polygon.y(j) * scaled_polygon.y(k) + scaled_polygon.y(k)^2) * ai;
            Ixy = Ixy + (scaled_polygon.x(j) * scaled_polygon.y(k) + 2 * scaled_polygon.x(j) * scaled_polygon.y(j) + 2 * scaled_polygon.x(k) * scaled_polygon.y(k) + scaled_polygon.x(k) * scaled_polygon.y(j)) * ai;
        end
    
        Iy(i) = Ix - (centroid_x^2 * A);
        Iz(i) = Iy_temp - (centroid_y^2 * A);
        Izy(i) = Ixy - (centroid_x * centroid_y * A);
    end

    % find angle theta for the principal axes
    two_theta = atan((-2 * Izy) ./ (Iz - Iy));

    % I1
    I1 = (Iy + Iz) / 2 - ((Iy - Iz) / 2) .* cos(two_theta) - Izy .* sin(two_theta)

    % I2
    I2 = (Iy + Iz) / 2 + ((Iy - Iz) / 2) .* cos(two_theta) - Izy .* sin(two_theta)


    % Bending stress distribution at the critical region
    M1 = M_y;
    M2 = M_z;
    y_val_edge = chordDist;
    z_val_edge = 0;

    % Scaling factor for bending moments
    scaling_factor = 1e-6;

    Bending_stress_edge = (M1 ./ I1) .* z_val_edge + (M2 ./ I2) .* y_val_edge;
    Bending_stress_edge = Bending_stress_edge * scaling_factor;
  
end
