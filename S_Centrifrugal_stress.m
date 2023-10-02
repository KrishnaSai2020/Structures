function centrifrugal_stress = S_Centrifrugal_stress(chordDis, r ,rho_Material, Omega) 

    % centrifrugal_stress finds the centrifrugal stress dist of the wind turbine
    % stress = force/area => centrifrugal stress dist = centri F dist/area dist
    % centrifurgal Force = rho * omega^2 * integral of areas [lubbock powerpoint]
    
    % Input arguments:
    %   chord_dist: the blades chord distribution
    %   r_dist: Contains the radial distances of the blade sections corresponding to the chord dist in chord_dist.
    % Output:
    %   centrifrugal_stress: Distribution of the centrifrugal stress along the
    %   blade.    
    
    % Area distribution is found using coordinates from Aerofoil_plotter - data
    % under - AerofoilCoord.mat (for a 1m length aerofoil)
    
    load("AerofoilCoord.mat",'Aerofoil_coord'); % Loads the 1m data
    
    % CALCULATE AREA OF 1 METER CHORD LENGTH AIRFOIL 
    A_1m = 0; % set as zero and add each sectional area
    for n = 1: max( size(Aerofoil_coord) ) - 1
        A_1m = abs((Aerofoil_coord(n,1)+Aerofoil_coord(n+1,2)) * 0.5 * (Aerofoil_coord(n,2)-Aerofoil_coord(n+1,2))) + A_1m;
        % Area of a 1 meter chord length airfoil; sum of all individual sections of the blade
    end

    % Now find actual area dist from chord dist by scaling the 1m area down by
    % the correct scale factor (find the integral using the trapezium rule)
    Area_dist = A_1m * chordDis.^2;     % Area distribution is found by multiplying a respective scale factor (the chord at that point) at each section
    A_y = Area_dist.*r;              % y variable for the integral
    delta_x = r(2) - r(1);      % since delta_x is constant, it can be calculated first
    %mean_area = mean(Area_dist);
    
    A_Integral = zeros(1,length(r));

    for n = 1:length(r)
        y = Area_dist(n:length(r)) * r(n);
        x = linspace(r(n),r(length(r)),length(y));
        if length(y) > 1
            A_Integral(n) = trapz(x, y);
        else
            A_Integral(n) = 0;
        end
    end

    %A_Integral = -A_Integral;

    % Now sub into stress equation
    F_x = rho_Material * Omega^2 .* A_Integral; % Total centrifugal force
%     figure
%     plot(r,F_x)
    centrifrugal_stress =  F_x ./ Area_dist; % Stress distribution of blade due to centrifugal force
end
