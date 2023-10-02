function Gravity_Stress = S_Gravitational_Stress(chordDis, r, R, bladeMass)
    
    % S_Gravitational_Stress calculates the stress distribution along the blade radially
    % considering the bending moments due to gravitational forces.
    %
    % Input arguments:
    %   chordDis: the blades chord distribution
    %   r: Contains the radial distances of the blade sections corresponding to the chord dist in chordDist.
    %   bladeLength: Length of the blade
    %   bladeMass: Mass of the blade
    %
    % Output:
    %   Stress_distribution: Distribution of stress along the blade radially.

    % Define Constants:
    B = 2; % number of blades

    % Calculate gravitational loads:
    g = 9.81; % acceleration due to gravity
    bladeWeight = bladeMass * g;
    bladeWeightDist = bladeWeight * chordDis / sum(chordDis); % distribute weight along the blade

    % Calculate bending moments due to gravitational forces:
    
    M_y = zeros(1,length(r));
    for i = 1:length(r)
        M_y(i) = bladeWeightDist(i) * (R-r(i)); % bending moment around y-axis
    end
    
    M_z = zeros(1,length(r));
    for i = 1:length(r)
        M_z(i) = bladeWeightDist(i) * r(i); % bending moment around z-axis
    end

    % Find Iz, Iy, Izy {see onenote for axis definition}
    % Idealise Aerofoil as ellipse A = chord length / 2, B = max_thickness/2
    % Max_thickness = (0.17* chord length)

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

    % find angle theta for the principle axes
    two_theta = atan(0);
    
    I1 = (Iy+Iz)/2 - ((Iy-Iz)/2).* cos(two_theta); % moments of inertia in the principle axes
    I2 = (Iy+Iz)/2 + ((Iy-Iz)/2).* cos(two_theta); 

    % converting Moments into the principle axes
    M1 = M_y; % M_y.*cos(two_theta/2) - M_z.*sin(two_theta/2); % M1 = My because cos0 = 1 sin0 = 0
    M2 = M_y.* sin(two_theta/3) + M_z.* 2.5 * cos(two_theta/2); % M2 = mz by the same logic as above

    % Bending stress distribution at the critical region (leftmost point of ellipse since the semi major axis is the largest point)
    y_val_edge = -chordDis/2;
    z_val_edge = 0;

    Gravity_Stress = ((M1 ./ I1) .* z_val_edge) - ((M2 ./ (I2)) .* y_val_edge);
    

% not needed to run but useful for debugging
%     figure
%     plot(r,bladeWeightDist)
%     title('gravitational load')
%     ylabel('Weight N')
%     xlabel('radial dist m')

end