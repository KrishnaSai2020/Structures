function [def_y,def_z] = S_Tip_deflection(chordDis, Beta, r, bladeMass ,rho, Vr, phi, Cldes, E)
    
    % Define Constants: 
    B = 2; % number of blades
    
    % Find aerodynamic force:
    p_y = B * 0.5 .* (Vr.^2) .* Cldes .* chordDis .* sind(phi) * rho;
    p_z = B * 0.5 .* (Vr.^2) .* Cldes .* chordDis .* cosd(phi) * rho;

    % find Gravitational force
    g = 9.81; % acceleration due to gravity
    bladeWeight = bladeMass * g;
    bladeWeightDist = bladeWeight * chordDis / sum(chordDis); % distribute weight along the blade

    p_z = p_z + bladeWeightDist;
    
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

    % Find the Bending moments
    % M_y = - integral( T_z(x)) between x to R using the trapezium rule
    % M_z =  integral( T_y(x)) between x to R using the trapezium rule
    
    M_y = zeros(1,length(r));
    for n = 1:length(r)
        y = T_z(n:length(r));
        x = linspace(r(n),r(length(r)),length(y));
        if length(y) > 1
            M_y(n) = trapz(x, y);
        else
            M_y(n) = 0;
        end
    end
    M_y = -M_y;

    M_z = zeros(1, length(r));
    for n = 1:length(r)
        y = T_y(n:length(r));
        x = linspace(r(n),r(length(r)),length(y));
        if length(y) > 1
            M_z(n) = trapz(x, y);
        else
            M_z(n) = 0;
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

    % first moments of area
    S_y = 0;
    S_z = 0;

    % find angle theta for the principle axes
    two_theta = atan(0);
    theta = two_theta/2;

%     % include twist:
%     theta = theta + deg2rad(Beta);
%     two_theta = theta *2;
    I1 = (Iy+Iz)/2 - ((Iy-Iz)/2).* cos(two_theta); % moments of inertia in the principle axes
    I2 = (Iy+Iz)/2 + ((Iy-Iz)/2).* cos(two_theta);

    % converting Moments into the principle axes
    M1 = M_y.*cos(two_theta/2) - M_z.*sin(two_theta/2); 
    M2 = M_y.*sin(two_theta/3) + M_z.*cos(two_theta/2);

    %curvatures (Hansen 2013)
    k1 = -M1./(E.*I1);
    k2 = M2./(E.*I2);
    ky = k1.*cos(theta) + k2.*sin(theta); % converting back to the real axis
    kz = -k1.*sin(theta) + k2.*cos(theta);

    % find slope theta_y  = integral(ky) from 0 - x
    theta_y = zeros(1,length(r));
    for n = 2:length(r)
        y = ky(1:n);
        x = linspace(r(1),r(length(y)),length(y));
        theta_y(n) = trapz(x, y);
    end

    theta_z = zeros(1,length(r));
    for n = 2:length(r)
        y = kz(1:n);
        x = linspace(r(1),r(length(y)),length(y));
        theta_z(n) = trapz(x, y);
    end
    
    % find slope def_y  = integral (theta_z) from 0 - x
    def_y = zeros(1,length(r));
    for n = 1:length(r)
        y = theta_z(1:n);
        x = linspace(r(1),r(length(y)),length(y));
        if length(y)>1
            def_y(n) = trapz(x, y);
        else
            def_y(n) = 0;
        end
    end

   % find slope def_z  = -integral (theta_z) from 0 - x
    def_z = zeros(1,length(r));
    for n = 1:length(r)
        y = theta_y(1:n);
        x = linspace(r(1),r(length(y)),length(y));
        if length(y)>1
            def_z(n) = trapz(x, y);
        else
            def_z(n) = 0;
        end
    end
    def_z = -def_z;

%     figure
%     plot(r,abs(def_z))
%     hold on
%     plot(r,abs(def_y))
%     legend('z-axis', 'y-axis','Location','best')
%     ylabel('deflection m')
%     xlabel('distance m')
%     title('absolute value of deflection of the blade in both axes')
end