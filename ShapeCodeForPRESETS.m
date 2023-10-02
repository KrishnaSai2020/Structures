%--------------------------HOUSEKEEPING-----------------------------
clear 
clc
close all
%---------------------------FORMATING------------------------------
set(0,'defaultfigurecolor',[1 1 1])
set(groot,'defaultAxesFontSize',16)
set(groot,'defaulttextfontsize',20)
set(groot,'defaultLineMarkerSize',4)
set(groot,'defaultLineLineWidth',0.5)
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')
set(groot, 'DefaultAxesBox', 'on')
set(groot, 'defaultAxesFontName','Cambria Math')

%-------------------------CONTROL VARIABLES-------------------------
Cldes=0.899; % desgin lift coefficient
R=0.25; % blade radius
B=2;%no of blades
AOAdes=6.604;%the design angle of attack
CLpoly=[0.083066667	0.350333333]; %Lift curve equation
CDpoly=[3.20E-05	-0.000240488	-0.000572652	0.028337532]; %Drag curve equation
Omega=100*pi; %rotational velocity
rho=1.19;%Sea level air density
n=67; %number of sections to divide the blade
N=50; %amount of TSR designs to test for power curve
nn=10; % number of r0 to test for 
ch=0.02;%hub chord
SafetyFactor=2;
r0=0.048; 
rho_Material=163; % Change accordingly from CAD
ac=0.2; %critical induction factor

%===========================CHANGE-TO-TEST============================
TSRdes=8.08;

%============================MANIPULATION==============================
linChange=+0.005; %the amount by which you would like to linearly alter the chord
degreeChord=1; %the degree of the chord distrubution
degreeTwist=1; %the degree of the twist distrubution

%===========================CONTROLLERS===============================
corrections=0;%turns on aerodynamic corrections
graph=1; %if graph ==1 then graphs will be enabled if graph ==0 graphs will be disabled (BE CAREFUL TAKES A STRONG COMPUTER TO RUN)
structureGraph=1; %when in structure  mode if on 1 will plot graphs if on 0 will not
nonopt=1; %if nonopt ==1 then the code will produce manufacturable results for Beta and chord ie linear twist and quadratic chord
mode=2; %if mode == 0  then the code will run to try fun optimium combinations of r0 and TSR and return them in GoodComb, if it is 1 it will use the inputed r0 and TSRdes to return the shape and features of the blade, if it is 3 it will calculate the strucutral abaility of the blade
%=====================================================================

if r0<0.048
    disp("The blade will stall")
elseif r0>0.075
    disp("Blade will not produce enough power to spin")
end   

%=============FIND THE OPTIMUM COMBINATIONS OF TSR AND r0============

%only input variable to this controller is r0, this section will take your
%r0 and between the range of the maximum and minimum TSR it will find TSR's
% that produce an acceptable power output between the maximum and minimum TSR
% you have selected


if mode==0
    TSRtest=linspace(7.85,15.71,N); %The N number of TSR's that will tested across the range 
    r0test=linspace(0.048,0.08,nn); 
    GoodComb=[];
    r=linspace(r0,R,n); %radius at set points
    [P,GoodTSR,GoodComb] = TSRdesFinder(TSRtest,Cldes,B,R,r0,n,degreeTwist,AOAdes,Omega,N,CLpoly,CDpoly,r,rho,GoodComb,graph,nonopt,linChange,degreeChord,ac,corrections);
    

%==============TEST A SINGULAR TSR AND r0 COMBINATION================

elseif mode~=0
    %disabling graphs if structural mode is selected
   if mode==2
       graph=0; 
   end

    r=linspace(r0,R,n); %radius at set points 
    [p,j,chordDis] = ChordDistrubtion(TSRdes,Cldes,B,R,r0,n,degreeChord,graph,nonopt,linChange); %finds chord distrubution
    [q,Beta] = TwistDistrubution(r0,R,TSRdes,AOAdes,n,degreeTwist,graph,nonopt); %finds twist distrubution
    [LSolidity,RSolidity] = SolidityDistrubution(B,r0,R,n,chordDis,TSRdes,graph); %finds solidity 
        
    if mode==1

        Q=[];
        P=[]; %initialize the power array
        Cp=[];%initialize the power coefficient array
        TSR=linspace(7.85,15.71,10);
        
        tiledlayout(1,10);
        % Initialise Vectors to store values in
        for i = 1:length(TSR)
            Uinf=Omega*R./TSR(i);
            TSRsec=(TSR(i).*r)/R; %Sectional tip speed ratio
            [a,ad,Cx,Cy,phi] = InductionFactorFinder(TSR(i),TSRsec,Beta,CLpoly,CDpoly,LSolidity,n,r,N,graph,B,R,ac,corrections); 
            Vr=(((Uinf.*(1-a)).^2)+(Omega.*r.*(1+ad)).^2).^0.5; %Resultant velocities per section
            Fx=B*0.5.*((Vr).^2).*Cx.*chordDis*((R-r0)/n);
            Fy=B*0.5.*((Vr).^2).*Cy.*chordDis*((R-r0)/n);
            Q(i)=sum(r.*Fy); %torque
            P(i)=Omega*Q(i); %power
            Cp(i)=P(i)/(0.5*rho*(Uinf^3)*(pi*R^2));
        end

        if graph==1
            hold off
            figure
            hold on
            plot(TSR,Cp,'r','LineWidth',2)
            xlabel('TSR');
            ylabel('Coefficient of power');
            title(sprintf('Cp versus TSR using design; TSR = %s', TSRdes));
            grid on;
            hold off
        end
    
        if graph==1
            figure
            hold on
            plot(TSR,P,'r','LineWidth',2)
            xlabel('TSR');
            ylabel('Power(W)');
            title(sprintf('Power at different TSRs using design; TSR = %s', TSRdes));
            grid on;
            hold off
        end
%===============================STRUCTURAL================================


    elseif mode==2

            Omega=SafetyFactor*Omega;
            
            Uinf=Omega*R./TSRdes;
            TSRsec=(TSRdes.*r)/R; %Sectional tip speed ratio
            [a,ad,Cx,Cy,phi] = InductionFactorFinder(TSRdes,TSRsec,Beta,CLpoly,CDpoly,LSolidity,n,r,N,graph,B,R,ac,corrections); 
            Vr=(((Uinf.*(1-a)).^2)+(Omega.*r.*(1+ad)).^2).^0.5; %Resultant velocities per section
            
            Fx=B*0.5.*((Vr).^2).*Cx.*chordDis*((R)/n);
            Fy=B*0.5.*((Vr).^2).*Cy.*chordDis*((R)/n);

            Fxeq=polyfit(r,Fx,7);
            Fyeq=polyfit(r,Fy,7);

            if structureGraph==1
                figure
                hold on
                plot(r,polyval(polyfit(r,Fx,7),r))
                plot(linspace(0,r(1),11),zeros(11)+Fx(1))
                hold off
            end
       
            Fx=[(zeros(1,10)+Fx(1)),Fx]; %filling in the fx array to include forces at the hob to root assuming them to be the worst case scenario ie same as the root chord 
            Fy=[(zeros(1,10)+Fy(1)),Fy]; %filling in the fy array to include forces at the hob to root assuming them to be the worst case scenario ie same as the root chord 

%----------------ADDED BY KRISHNA-----------------------------------------%
            Bending_stress_edge = S_Bending_stress(chordDis, Beta, r ,rho, Vr, phi, Cldes);
            centrifrugal_stress = S_Centrifrugal_stress(chordDis, r ,rho_Material, Omega);
            
            bladeMass = 0.02; %example mass of 20g used please change with data from CAD

            % Extra corrections 
            Gravity_Stress = S_Gravitational_Stress(chordDis, r, R , bladeMass);
            Gyroscopic_Stress_edge = S_Gyroscopic_stress(chordDis, Beta, r ,rho, Vr, phi, Cldes);

            nexttile
            Total_stress = centrifrugal_stress + Bending_stress_edge + Gyroscopic_Stress_edge + Gravity_Stress;
            Total_stress_noextra = centrifrugal_stress + Bending_stress_edge;
            hold on
            plot(r,Total_stress,'r')
            yline(25000000,'-k','Yield stress Balsa')
            title('Total stress')
            ylabel('Stress pascals')
            xlabel('Radial distance metres')

%             nexttile
%             Total_stress_with_gyro = Total_stress - Gravity_Stress;
%             plot(r,Total_stress_with_gyro,'Color','b','LineStyle',':')
%             Total_stress_noextra = centrifrugal_stress + Bending_stress_edge;
%             plot(r,Total_stress_noextra,'k','LineStyle','--')
%             legend('with gravity stress','with just gyroscopic stress','without add ons','Location','best')
%             hold off
            
%---------------IF YOU WANT TO PLOT EACH COMPONENT SEPERATELY-------------%
%             nexttile
%             hold on
%             plot(r,Bending_stress_edge,'Color','r')
%             plot(r,centrifrugal_stress,'Color', 'M','LineStyle','--')
%             plot(r,Gravity_Stress,'color','k','LineStyle',':')
%             plot(r,Gyroscopic_Stress_edge, 'LineStyle','-.')
%             title('individual stresses')
%             xlabel('distance m')
%             ylabel('stress Pa')
%          
%             legend('Bending','centrifrugral','gravity','gyro','Location','best')
%-------------------------------------------------------------------------%
            
            % I HAVE USED THE VALUE FOR Balsa
            E = 0.5 *10^9 ; % Youngs modulus of blade

            [def_y,def_z] = S_Tip_deflection(chordDis, Beta, r, bladeMass, rho, Vr, phi, Cldes, E);

            nexttile
            hold on
            plot(r,abs(def_z),'r')
            plot(r,abs(def_y),'k','LineStyle','--')
            legend('z-axis', 'y-axis','Location','best')
            ylabel('deflection m')
            xlabel('distance m')
            title('Absolute value of deflection of the blade in both axes')
            hold off
            
            Hand_calc_centri_stress = [0.8197, 0.73763, 0.655218, 0.5738, 0.491421, 0.40986, 0.32788, 0.2459, 0.1639, 0.08197, 0];
            Hand_calc_centri_stress = Hand_calc_centri_stress * 10^6; % convert to pascals
            Hand_calc_x_val = linspace(0,0.228,11) + 0.022;

            percentage_change_stress = ( (Gravity_Stress )./Total_stress_noextra) * 100;
            percentage_change_stress_gyro = ( (Gyroscopic_Stress_edge )./Total_stress_noextra) * 100;
            
%             dist_new = 
%             percent_change_new = perce

%             figure
%             plot(r,percentage_change_stress,'Color','r')
%             title('Percentage change in total stress with Gravity add on')
%             xlabel('distance m')
%             ylabel('percentage change')
%             hold off

            
%             centrifrugal_stress = S_Centrifrugal_stress(chordDis, r ,rho_Material, Omega);
%             figure
%             hold on
%             plot(r, centrifrugal_stress, 'Color','r', 'LineStyle','-')
%             plot(Hand_calc_x_val, Hand_calc_centri_stress, 'Color','k', 'LineStyle','--')
%             
%             legend('Real blade', 'cantilever')
%             xlabel('distance m')
%             ylabel('stress Pa')
%             hold off
    end
end

