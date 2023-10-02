function [LSolidity,RSolidity] = SolidityDistrubution(B,r0,R,n,chordDis,TSR,graph)
%Calculates the local and rotor solidity of the blade 

%-------------------------------INPUTS--------------------------------
%r0 = distance of root chord from hub
%R = blade radius
%TSR = tip speed ratio
%n = number of section to analyse
%B = number of blades
%chordDis= chord distrubution

%------------------------------OUTPUTS--------------------------------
%LSolidity = local solidity of the blade
%RSolidity = Rotor solidity of the blade


r=linspace(r0,R,n); %radius at set points 
LSolidity = (B./(2*pi*r)).*chordDis; %local solidity
RSolidity = (B/(pi*R^2))*(trapz(r,chordDis)); % Use the trapezoidal rule to integrate the chord disturbution curve, and then sub into the rotor solidity equationv


% If grpahs enabled, plots the local solidity against the radius
if graph==1
    figure
    hold on
    plot(r, LSolidity, 'b-', 'LineWidth', 2);
    xlabel('Position along blade (m)');
    ylabel('Solidity');
    title(sprintf('Solidity distrubution of a Wind Turbine Blade; TSR = %s', TSR));
    grid on;
    hold off
end

%Checks if the blades geometry is sutible for analyse via BEM
if LSolidity<1
    disp('Holds for BEM')
end

end