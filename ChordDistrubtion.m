function [p,j,chordDis] = ChordDistrubtion(TSR,Cldes,B,R,r0,n,degreeChord,graph,nonopt,linChange)
%This function finds the equaiton for the chord distrubtion of the wind
%turbine blade starting at distance r0 from the centre, the function will
%begin by finding the root chord cr using n the number of sections


% Set axis limits
xLimit = R;
yLimit = R;

%initialize an array to store all polynomial fit equations
p=zeros(degreeChord);

r=linspace(r0,R,n); %radius at set points 
chordDis = ((16*pi*(R^2))/(9*B*Cldes*TSR^2)./r); %chord distrubution
j=(16*pi*(R^2))/(9*B*Cldes*TSR^2);

for i = 1:degreeChord
    poly = polyfit(r, chordDis, i); % Coefficients of the polynomial
    for k=1:length(poly)
        p(i,k)=poly(k);
    end
end

%If wanting the non optimized chord distrubution returns a quadratic based
%of the nonopt input
if nonopt==1
    chordDis=linChange+polyval(polyfit(r,chordDis,degreeChord),r); %produces a quadratic chord
end

% Plots the chord distrubtion
if graph==1
    figure
    hold on
    plot(r, chordDis, 'b-', 'LineWidth', 2);
    xlabel('Position along blade (m)');
    ylabel('Chord length (m)');
    title(sprintf('Chord Distribution of a Wind Turbine Blade; TSR = %s', TSR));
    xlim([0 xLimit])
    ylim([0 yLimit])
    grid on;
    hold off
end


end