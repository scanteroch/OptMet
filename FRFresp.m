% Function: FRF response
function [pFRF]=FRFresp(num, total_massper, th, fmin, fmax)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
%   num  - Number of resonators
%   total_massper  - The percentage mass
%   th - Individual sample of the unknown model parameters vector \th
%   fmin - Minimum frequency in the design interval for vibration
%          attenuation
%   fmax - Maximum frequency in the design interval for vibration
%          attenuation
%
% Outputs:
%
%   pFRF - Sum of the FRF response in the interval [fmin, fmax]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ji=0;
fstep=1;
freq=fmin:fstep:fmax; %calculated frequency range

for f=freq
    ji=ji+1;
    w=2*pi*f;
    m=num; %total number of resonators added to the U beam
    %% physical parameters of the materials
    rou=th(2);
    poi=0.3;
    yita=0.02;
    E=th(1)*10^6*(1+i*yita);
    
    %% geometrical parameters of the beam
    L_t=0.3; %total length of the beam
    
    Lw=0.049; %width of the beam
    Lh=0.01;%height of the side beam
    tw1=0.002;%wall thickness of the beam 
    tf=0.005;%wall thickness of the face
    Au=Lw*tf+2*Lh*tw1;%cross section area U beam
    m_beam=rou*Au*L_t;%total mass of the beam
    
    Caxis=(Lw*tf*tf/2+tw1*Lh*(Lh+tf))/(Lw*tf+Lh*tw1*2); % centroid
    IyO=Lw*tf^3/3+tw1*((tf/2+Lh)^3-(tf/2)^3)*2/3; 
    IyC=IyO-Caxis^2*Au; %moment of inertial
    k=(rou*Au/E/IyC)^0.25*w^0.5; %wavenumber in each segment of U beam
    
    if m>0
    
       %% paramters of the connecting beam & the resonators
        df=L_t/m; %length of each part
        Ms1=total_massper*m_beam/2/m*ones(1,m);
        Ms2=total_massper*m_beam/2/m*ones(1,m);
        Js=zeros(1,m);
        p=0.5*ones(1,m); %position of the resonators  
        f0=320;
        hs1=0.0025*ones(1,m); %Side 1 height of connecting beam
        hs2=0.0025*ones(1,m); %side 2 height of connecting beam
        Ls=0.02*ones(1,m);%length of the connecting beam
%         bs1=0.0025*ones(1,m);%Side 1 width of connecting beam
%         bs2=0.0025*ones(1,m);%Side 2 width of connecting beam
       bs1=(2*pi*f0)^2*Ms1.*Ls.^3*12/3/E./hs1.^3;
       bs2=(2*pi*f0)^2*Ms2.*Ls.^3*12/3/E./hs2.^3;

        As1=bs1.*hs1;
        As2=bs2.*hs2;
        
        Is1=bs1.*hs1.^3/12;
        Is2=bs2.*hs2.^3/12;
        ks1=(rou*As1/E./Is1).^0.25*w^0.5;
        ks2=(rou*As2/E./Is2).^0.25*w^0.5;

  
        %% intial values  
        A1=zeros(1,m);
        A2=zeros(1,m);
        B1=zeros(1,m);
        B2=zeros(1,m);
        C1=zeros(1,m);
        C2=zeros(1,m);
        D1=zeros(1,m);
        D2=zeros(1,m);
        Ns1=zeros(1,m);
        Ns2=zeros(1,m);
        T=eye(4);
    
        R=[1 1 1 1;
           -i -1 i 1;
           -1 1 -1 1;
           i*E*IyC*k^3 -E*IyC*k^3 -i*E*IyC*k^3 E*IyC*k^3];   
        Rr=[1/4 i/4 -1/4 -i/4/(E*IyC*k^3);
            1/4 -1/4 1/4 -1/4/(E*IyC*k^3);
            1/4 -i/4 -1/4 i/4/(E*IyC*k^3);
            1/4 1/4 1/4 1/4/(E*IyC*k^3)];  %inv(R)
    
        Lc=[1 1 1 1;
            -i -1 i 1;
            -E*IyC*k^2 E*IyC*k^2 -E*IyC*k^2 E*IyC*k^2;
            i*E*IyC*k^3 -E*IyC*k^3 -i*E*IyC*k^3 E*IyC*k^3];
        Lcr=[1/4,  i/4, -1/(4*E*IyC*k^2), -i/(4*E*IyC*k^3);
            1/4, -1/4,  1/(4*E*IyC*k^2), -1/(4*E*IyC*k^3);
            1/4, -i/4, -1/(4*E*IyC*k^2),  i/(4*E*IyC*k^3);
            1/4,  1/4,  1/(4*E*IyC*k^2),  1/(4*E*IyC*k^3)];%inv(Lc)

        %%
        for n=1:m
            Ts11_1=[1 1;
                    -i*ks1(n) -ks1(n)];
            Ts11_2=[1 1;
                    -i*ks2(n) -ks2(n)];
            Ts12_1=[1 1;
                    i*ks1(n) ks1(n)];
            Ts12_2=[1 1;
                    i*ks2(n) ks2(n)];
            Ts21_1=[(i*E*Is1(n)*ks1(n)^3+Ms1(n)*w^2)*exp(-i*ks1(n)*Ls(n)) (-E*Is1(n)*ks1(n)^3+Ms1(n)*w^2)*exp(-ks1(n)*Ls(n));
                    (-E*Is1(n)*ks1(n)^2+i*ks1(n)*Js(n)*w^2)*exp(-i*ks1(n)*Ls(n)) (E*Is1(n)*ks1(n)^2+Js(n)*w^2*ks1(n))*exp(-ks1(n)*Ls(n))];
            Ts21_2=[(i*E*Is2(n)*ks2(n)^3+Ms2(n)*w^2)*exp(-i*ks2(n)*Ls(n)) (-E*Is2(n)*ks2(n)^3+Ms2(n)*w^2)*exp(-ks2(n)*Ls(n));
                    (-E*Is2(n)*ks2(n)^2+i*ks2(n)*Js(n)*w^2)*exp(-i*ks2(n)*Ls(n)) (E*Is2(n)*ks2(n)^2+Js(n)*w^2*ks2(n))*exp(-ks2(n)*Ls(n))];
            Ts22_1=[(-i*E*Is1(n)*ks1(n)^3+Ms1(n)*w^2)*exp(i*ks1(n)*Ls(n)) (E*Is1(n)*ks1(n)^3+Ms1(n)*w^2)*exp(ks1(n)*Ls(n));
                    (-E*Is1(n)*ks1(n)^2-i*ks1(n)*Js(n)*w^2)*exp(i*ks1(n)*Ls(n)) (E*Is1(n)*ks1(n)^2-ks1(n)*Js(n)*w^2)*exp(ks1(n)*Ls(n))];
            Ts22_2=[(-i*E*Is2(n)*ks2(n)^3+Ms2(n)*w^2)*exp(i*ks2(n)*Ls(n)) (E*Is2(n)*ks2(n)^3+Ms2(n)*w^2)*exp(ks2(n)*Ls(n));
                    (-E*Is2(n)*ks2(n)^2-i*ks2(n)*Js(n)*w^2)*exp(i*ks2(n)*Ls(n)) (E*Is2(n)*ks2(n)^2-ks2(n)*Js(n)*w^2)*exp(ks2(n)*Ls(n))];
            Ts11r_1=[0.5+i*0.5 (0.5+0.5*i)/ks1(n);
                    0.5-0.5*i (-0.5-0.5*i)/ks1(n)];%inv(Ts11_1)
            Ts11r_2=[0.5+i*0.5 (0.5+0.5*i)/ks2(n);
                    0.5-0.5*i (-0.5-0.5*i)/ks2(n)];%inv(Ts11_2)
            Tst_1=-Ts11r_1*Ts12_1;
            Tst_2=-Ts11r_2*Ts12_2;
            Tsx_1=-Ts21_1*Ts11r_1*Ts12_1+Ts22_1;
            Tsx_2=-Ts21_2*Ts11r_2*Ts12_2+Ts22_2;
            Tsr_1=inv(Tsx_1);
            Tsr_2=inv(Tsx_2);
            Tsm_1=Tst_1*Tsr_1;
            Tsm_2=Tst_2*Tsr_2;
            C1(n)=-Ms1(n)*w^2*Tsr_1(1,1);
            C2(n)=-Ms2(n)*w^2*Tsr_2(1,1);
            D1(n)=-Ms1(n)*w^2*Tsr_1(2,1);
            D2(n)=-Ms2(n)*w^2*Tsr_2(2,1);
            A1(n)=-Ms1(n)*w^2*Tsm_1(1,1);
            A2(n)=-Ms2(n)*w^2*Tsm_2(1,1);
            B1(n)=-Ms1(n)*w^2*Tsm_1(2,1);
            B2(n)=-Ms2(n)*w^2*Tsm_2(2,1);
            Ns1(n)=-E*Is1(n)*ks1(n)^3*(i*A1(n)-B1(n)-i*C1(n)+D1(n)); 
            Ns2(n)=-E*Is2(n)*ks2(n)^3*(i*A2(n)-B2(n)-i*C2(n)+D2(n)); 
            Ns(n)=Ns1(n)+Ns2(n);
            la_l=[exp(-i*k*p(n)*df) 0 0 0;
                  0 exp(-k*p(n)*df) 0 0;
                  0 0 exp(i*k*p(n)*df) 0;
                  0 0 0 exp(k*p(n)*df)];
            la_r=[exp(-i*k*(1-p(n))*df) 0 0 0;
                  0 exp(-k*(1-p(n))*df) 0 0;
                  0 0 exp(i*k*(1-p(n))*df) 0;
                  0 0 0 exp(k*(1-p(n))*df)];
            Rn=[1 1 1 1;
                -i -1 i 1;
                -1 1 -1 1;
                i*E*IyC*k^3+Ns(n) -E*IyC*k^3+Ns(n) -i*E*IyC*k^3+Ns(n) E*IyC*k^3+Ns(n)];
            Tn=la_r*Rr*Rn*la_l;
            T=Tn*T; 
        end

        H2=[i -1 -i 1;
            -1 1 -1 1]*T;
        H1=[i*k^3*E*IyC -k^3*E*IyC -i*k^3*E*IyC k^3*E*IyC;
            -1 1 -1 1];
        Th=[H1;H2];
        Tz=T/Th;
        alfam1=Tz(1,1);
        beltam1= Tz(2,1);
        xim1=Tz(3,1);
        detam1=Tz(4,1);
        disFRF(ji)=20*log10(abs(alfam1+beltam1+xim1+detam1));  
        
    elseif m==0
           la_T=[exp(-i*k*L_t) 0 0 0;
                  0 exp(-k*L_t) 0 0;
                  0 0 exp(i*k*L_t) 0;
                  0 0 0 exp(k*L_t)];
           H2=[i -1 -i 1;
            -1 1 -1 1]*la_T;
           H1=[i*k^3*E*IyC -k^3*E*IyC -i*k^3*E*IyC k^3*E*IyC;
               -1 1 -1 1];
           Th=[H1;H2];
           Tz=la_T*inv(Th);
           alfam1=Tz(1,1);
           beltam1= Tz(2,1);
           xim1=Tz(3,1);
           detam1=Tz(4,1);
           disFRF(ji)=20*log10(abs(alfam1+beltam1+xim1+detam1));
    end
end

% Performance index: \sum disFRF \in [freq_interest] ---> Minimize it
pFRF=sum(disFRF);
end