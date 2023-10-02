function [P,GoodTSR,GoodComb] = TSRdesFinder(TSRtest,Cldes,B,R,r0,n,degree,AOAdes,Omega,N,CLpoly,CDpoly,r,rho,GoodComb,graph,nonopt,linChange)

GoodTSR=zeros(1,N); %initialize array to store TSR's that past the test
Q=zeros(1,N); % initialize torque array
P=zeros(N); %initialize the power array
Cp=zeros(1,N);%initialize the power coefficient array

for ii=1:N
    TSRdes=TSRtest(ii); %The ii'th TSR design being tested within the loop
    [p,j,chordDis] = ChordDistrubtion(TSRdes,Cldes,B,R,r0,n,degree,graph,nonopt,linChange); %finds chord distrubution
    [q,Beta] = TwistDistrubution(r0,R,TSRdes,AOAdes,n,degree,graph,nonopt); %finds twist distrubution
    [LSolidity,RSolidity] = SolidityDistrubution(B,r0,R,n,chordDis,TSRdes,graph); %finds solidity 


    TSR=linspace(7.85,15.71,N); %Generate range of TSR's to plot CP curve
    
    tiledlayout(1,N); %Displaying graphs

    %Based off of the number of TSR's to test in order order to plot the
    %power curve runs BEM
    for i = 1:length(TSR)
        Uinf=Omega*R./TSR(i); %finds free stream velocity
        TSRsec=(TSR(i).*r)/R; %Sectional tip speed ratio
        [a,ad,Cx,Cy] = InductionFactorFinder(TSR(i),TSRsec,Beta,CLpoly,CDpoly,LSolidity,n,r,N,graph); %Calculates the induction factors and the normal and tangential forces
        Vr=(((Uinf.*(1-a)).^2)+(Omega.*r.*(1+ad)).^2).^0.5; %Resultant velocities per section
        Fx=B*0.5.*((Vr).^2).*Cx.*chordDis*((R-r0)/n); %Normal sectional force
        Fy=B*0.5.*((Vr).^2).*Cy.*chordDis*((R-r0)/n); %Tangential sectional force
        Q(i)=sum(r.*Fy); %torque array
        P(ii,i)=Omega*Q(i); %power array
        Cp(i)=P(ii,i)/(0.5*rho*(Uinf^3)*(pi*R^2)); %power coefficients
    end

    if isnan(P(ii,N))    
    else
        if P(ii,1)>30
            GoodComb=[GoodComb;{r0,TSRdes}];
        end
    end
end



end