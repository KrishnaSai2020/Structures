function [a,ad,Cx,Cy,phi] = InductionFactorFinder(TSR,TSRsec,Beta,CLpoly,CDpoly,LSolidity,n,r,N,graph,B,R,ac,corrections)

diffa=1; %initialize the difference between the old and new induction factors
diffad=1; %intialize the difference between the old and new tangential induction factors
a=zeros(1,n);
ad=zeros(1,n);

while (diffa >0.02)&&(diffad>0.02)

    phi=atand(((1-a)./(1-ad))./TSRsec);
    AOAsec=phi-Beta;

    f=(B*(R-r))./(2.*r.*sind(phi));
    F=(1/90)*acosd(exp(-f));
    
    if corrections==0
        F=1;
    end

    CLsec=polyval(CLpoly,AOAsec);
    CDsec=polyval(CDpoly,AOAsec);
    
    Cx= cosd(phi).*CLsec+sind(phi).*CDsec;
    Cy= sind(phi).*CLsec-cosd(phi).*CDsec;
    
    Consta=(LSolidity./(4.*F.*(sind(phi).^2))).*Cx;
    Constad=(LSolidity./(4.*F.*(sind(phi).*cosd(phi)))).*Cy;
    
    aold=a;
    adold=ad;

    a=Consta./(1+Consta);
    ad=Constad./(1-Constad);

    if corrections==1
        for i=length(n)
            if a(n)<=ac
                a(n) = 1/(((4*F(n)*(sin(phi(n))^2)/LSolidity(n)*Cx(n))) + 1);
            elseif a(n)>ac
                K=(4*F(n)*(sind(phi(n)))^2)/(LSolidity(n)*Cx(n));
                a(n)=0.5*(2+K*(1-2*ac)-(((K*(1-2*ac)+2)^2 +4*(K*ac^2 -1))*0.5));
            end
        end
    end    

    diffa = max(abs(aold - a)) / max(abs(aold)); % finds the change between the maximum values of the induction factor
    diffad = max(abs(adold - ad)) / max(abs(adold)); % finds the change between the maximum values of the tangential induction factor between iterations
end


if graph==1
    nexttile
    plot(r,a, 'b',r,ad,'r','LineWidth',2)
    title(TSR);
end
end