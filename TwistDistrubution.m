function [p,Beta] = TwistDistrubution(r0,R,TSR,AOAdes,n,degreeB,graph,nonopt)
%This function determines the twist distrubution for the wing turbine blade
%using a known euqation from literature

%----------------------------INPUTS-------------------------------------
%r0 = distance of root chord from hub
%R = blade radius
%TSR = tip speed ratio
%AOAdes = design angle of attack based of off the aerofoil
%n = number of section to analyse
%---------------------------Outputs-------------------------------------
%Beta = twist distrubution
%p = twist distrubution polyfit coefficients



p=zeros(degreeB); %initializes an array to store the twist distrubution equations
r = linspace(r0, R, n); % Linear taper model
Beta=atand(2*R./(3*TSR.*r)) - AOAdes; %twist equation from literature


%if non optimal mode enabled it displays the twist distrubution for a
%linear distrubution
if nonopt==1  
    Beta=polyval(polyfit(r,Beta,1),r); 
end


%if graphs enables creates a figure
if graph==1
    figure 
    hold on
end


%generates the polyfit coefficients for all polyfits up until the degree
%chosen
for i = 1:degreeB
    poly = polyfit(r, Beta, i); % Coefficients of the polynomial
    y=polyval(poly, r);
    plot(r,y, 'LineWidth', 2)
    for k=1:length(poly)
        p(i,k)=poly(k);
    end
end

hold off %disables the hold on the graphs


%If graphs enabled plots the twist distrubtion for the turbine blade
if graph==1
    figure 
    hold on
    plot(r, Beta, 'b-', 'LineWidth', 2);
    xlabel('Normalized position along blade');
    ylabel('Twist angle (degrees)');
    title('Twist Distribution of a Wind Turbine Blade');
    grid on;
    hold off
end

end