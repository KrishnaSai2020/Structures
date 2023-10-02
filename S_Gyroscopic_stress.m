
function Gyroscopic_stress_edge = S_Gyroscopic_stress(chordDis, Beta, r ,rho, Vr, phi, Cldes) 

    % S_Gyroscopic_stress finds the stress dist of the wind turbine due
    % to Gyroscopic forces
    %
    % Input arguments:
    %   chord_dist: the blades chord distribution
    %   r: Contains the radial distances of the blade sections corresponding to the chord dist in chord_dist.
    %   rho: Density
    %   CL_des: CL_des
    %   W: is the inlet velocity in the rotating (blade) frame of reference. [Lubbock Aero page 19 eq 18]
    %   phi: is the relative inflow angle with respect to the tangential direction.
    %
    % Output:
    %   gyroscopic_stress: Distribution of the gyroscopic stress along the
    %   blade.
    
    % Define Constants: 
    B = 2; % number of blades
    
    % Find aerodynamic force:
    p_y = B * 0.5 .* (Vr.^2) .* Cldes .* chordDis .* sind(phi) * rho;
    p_z = B * 0.5 .* (Vr.^2) .* Cldes .* chordDis .* cosd(phi) * rho;
    
    % Find the shear force
    % T_z = integral( p_y(x)) between x to r using the trapezium rule
    % T_y = integral( p_y(x)) between x to r using the trapezium rule
    
    delta_x = r(2) - r(1); % x difference between each trapezium
    
    T_z = zeros(1, length(r)); % Initialize T_z with zeros    
    for n = 1:length(r)
        y = p_z(n:length(r));
        x = linspace(r(n),r(length(r)),length(y));
        if length(y) >1
            T_z(n) = trapz(x, y);
        else
            T_z(n) = 0;
        end
    end

    T_y = zeros(1, length(r)); % Initialize T_z with zeros
    for n = 1:length(r)
        y = p_y(n:length(r));
        x = linspace(r(n),r(length(r)),length(y));
        if length(y) >1
            T_y(n) = trapz(x, y);
        else
            T_y(n) = 0;
        end
    end

    % Find Iz, Iy, Izy {see onenote for axis definition}
    % Idealise Aerofoil as ellipse A = chord length / 2, B = max
    % thickness/2
    % 
    % thickness(0.17* chord length)

    % Iz
    Iz = [];
    for i = 1:max(size(r))
        B = 0.17*0.5*chordDis(i);
        A = chordDis(i)/2;
        Iz(i) = (B * A^3) * (pi/4) ;
    end

    % Iy
    Iy = [];
    for i = 1:max(size(r))
        B = chordDis(i)/2;
        A = 0.17 * 0.5 * chordDis(i);
        Iy(i) = (B * A^3) * (pi/4);
    end

    % Izy
    Izy = 0;

    omega_yaw = 3 - 0.01 * (pi * 0.25^2 - 2); 

    M_gyroscopic_y = zeros(1,length(r));
    for i = 1:length(r)
        M_gyroscopic_y(i) = 2 * Iy(i) * 628.3185 * omega_yaw * cos(0);
    end

    M_gyroscopic_z = zeros(1,length(r));
    for i = 1:length(r)
        M_gyroscopic_z(i) = 2 * Iz(i) * 628.3185 * omega_yaw * cos(0);
    end
    
    M_y = M_gyroscopic_y;
    M_z = M_gyroscopic_z;
    
    % find angle theta for the principle axes
    two_theta = atan(0);
    theta = two_theta/2;

    % include twist:
    theta = theta + deg2rad(Beta);
    two_theta = theta *2;
    
    I1 = (Iy+Iz)/2 - ((Iy-Iz)/2).* cos(two_theta); % moments of inertia in the principle axes
    I2 = (Iy+Iz)/2 + ((Iy-Iz)/2).* cos(two_theta); 

    % converting Moments into the principle axes
    M1 = M_y; % M_y.*cos(two_theta/2) - M_z.*sin(two_theta/2); % M1 = My because cos0 = 1 sin0 = 0
    M2 = M_z; % M_y.* sin(two_theta/3) + M_z.* cos(two_theta/2); % M2 = mz by the same logic as above

    % Bending stress distribution at the critical region (leftmost point of ellipse since the semi major axis is the largest point)
    y_val_edge = -chordDis/2;
    z_val_edge = 0;

    Gyroscopic_stress_edge = ((M1 ./ I1) .* z_val_edge) - ((M2 ./ (I2)) .* y_val_edge);
     
end
