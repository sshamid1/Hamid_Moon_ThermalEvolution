function [TC, TS, TLi, Tdis, dTCdt, dridt, TC_nuc, ...
          Qsc, Qrc, Qgc, Qlc, Qic, Qad, TDM, c_oc, ...
          Psc, Prc, Pgc, Plc, Pic, Pk, P_total, Bs_core, B_ML, B_CIA, B_MAC,...
          Ro_ML, Ro_CIA, Ro_MAC, mag_ML, mag_CIA, mag_MAC] = ...
                    runEnergyCalc(t, QC, TC, ri, core_Kppm, c0, kc)
    
const_Moon;                       % Load constants from other script
TL0 = 1862 + dTLdc*(c0-0.06);   % Liquidus temperature at center (K).
                                %TL0 = 1920 + dTLdc*(c0-0.04);      
                                  % Was originally 1800 but I changed 
                                  % it based off of liquidus check we performed
TC_nuc = TL0 * (1-(rc/Lp)^2-Ap*(rc/Lp)^4)^gamm; % CMB T when IC starts (K)

%%% RUN THE MODEL
Qrc = Mc*h0C*core_Kppm*exp(-HlamC*t); % Radiogenic heating in the core (W)

% Create an inner core if there isn't one but there should be
if (ri < 1e4 && TC < TC_nuc)
    ri = 1e4;
end

% If there is an inner core, calculate TC self-consistently (ignore input TC)
if ri > 0    % YES inner core!
    Mic = Mc - 4/3*pi*rho0*Lp^3*(fc(rc/Lp,0,Ap) - fc(ri/Lp,0,Ap)); % Mass of the inner core (kg)
    TLi = TL0 - K0*dTLdP*(ri/Lp)^2 + dTLdc*c0*(ri/Lp)^3/fc(rc/Lp,0,Ap);  % Temperature at the inner core boundary
    TC = TLi * (1-(rc/Lp)^2 - Ap*(rc/Lp)^4)^gamm / (1-(ri/Lp)^2 - Ap*(ri/Lp)^4)^gamm; % Temperature at the CMB
    dTLdri = -2*(K0*dTLdP*ri/Lp^2) + (3*dTLdc*c0*ri^2 / (Lp^3*fc(rc/Lp,0,Ap))); % Growth rate of inner core
    rhoi = rho0 * (1 -(ri/Lp)^2 - Ap*(ri/Lp)^4);  % Density at the inner core boundary 
    gi = 4/3*pi*G*rho0*ri*(1-0.6*(ri/Lp)^2 - 3/7*Ap*(ri/Lp)^4); % Gravitational acceleration at inner core boundary

    dTadP = gamm*TLi/K0;     % Relationship between cooling rate and inner core growth rate
    dridTC = min(-(1/(dTLdP-dTadP)) * TLi/(TC*rhoi*gi),0);
        
    Psc = (-4/3*pi*rho0*Cc*Lp^3 * ...  
            (1-(ri/Lp)^2-Ap*(ri/Lp)^4)^(-gamm) * ...
            (dTLdri+2*gamm*TLi*ri/Lp^2 * ...
            (1+2*Ap*(ri/Lp)^2)/(1-(ri/Lp)^2-Ap*(ri/Lp)^4)) * ...  
            (fc(rc/Lp, gamm, Ap) - fc(ri/Lp, gamm, Ap)));
    Pgc = (8*pi^2*c0*G*rho0^2*bet*ri^2*Lp^2 / fc(rc/Lp,0,Ap) * ...
            (fchi(rc/Lp,ri,Lp) - fchi(ri/Lp,ri,Lp)));
    Plc = 4*pi*ri^2*rhoi*TLi*DSc;
    Pic = Cc*Mic*dTLdP*K0*(2*ri/Lp^2 + 3.2*ri/Lp^5);
        
    Qsc_til = Psc * dridTC;
    Qgc_til = Pgc * dridTC;
    Qlc_til = Plc * dridTC;
    Qic_til = Pic * dridTC;
    
% If there is no inner core, things are simpler...
else       
    Qsc_til = -4/3*(pi*rho0*Cc*Lp^3 * fc(rc/Lp, gamm, Ap) * ...
        (1-(rc/Lp)^2 - Ap*(rc/Lp)^4)^(-gamm)); 
    Qgc_til = 0;
    Qlc_til = 0;
    Qic_til = 0;
        
    TLi = TC*(1-(rc/Lp)^2 - Ap*(rc/Lp)^4)^(-gamm);
end
    
% Calculate cooling rate
dTCdt = -(QC - Qrc)/(Qsc_til + Qgc_til + Qlc_til + Qic_til); % Cooling rate of the core
if ri > 0
    dridt = -dridTC * dTCdt;
else
    dridt = 0;
end

% Absolute values of all energetic terms    
Qsc = -Qsc_til * dTCdt;
Qgc = -Qgc_til * dTCdt;
Qlc = -Qlc_til * dTCdt;
Qic = -Qic_til * dTCdt;
    
%%% ENERGY BUDGET FOR THE DYNAMO
EK = 16*pi*gamm^2*kc*Lp*(fk(rc/Lp,Ap)-fk(ri/Lp,Ap)); % Energy sink associated with thermal conductivity
    
TS = TLi*((1-(ri/Lp)^2-Ap*(ri/Lp)^4)^(-gamm) * ...         % Temperature associated with secular cooling. heat released as core cools
           (fc(rc/Lp,gamm,Ap) - fc(ri/Lp,gamm,Ap)) / ...
           (fc(rc/Lp,0,Ap) - fc(ri/Lp,0,Ap)));
Tdis = ((TLi / (1 - (ri/Lp)^2 - Ap*(ri/Lp)^4)^gamm) * ... % Temperature associated with average dissipation
        ((fc(rc/Lp,0,Ap) - fc(ri/Lp,0,Ap)) / (fc(rc/Lp,-gamm,Ap) - fc(ri/Lp,-gamm,Ap))));

% Dissipation budget (heats multiplied by efficiency terms)
Plc = (Tdis*(TLi-TC)/(TLi*TC))*Qlc;
Pic = (Tdis*(TLi-TC)/(TLi*TC))*Qic;
Pgc = (Tdis/TC)*Qgc;
Prc = ((Tdis-TC)/TC)*Qrc;
Psc = (Tdis*(TS-TC)/(TS*TC))*Qsc;
Pk = Tdis*EK;

Qad = (TS*TC)*EK/(TS-TC); % Adiabatic heat flow (W)

% Composition of the outer core
Mic = 4/3*pi*rho0*Lp^3*fc(ri/Lp, 0, Ap); % Mass of the inner core (kg)
c_oc = c0*Mc/(Mc-Mic);         % Abundance of sulfur in the outer core

% SCALING LAW: MAGNETIC FIELD STRENGTH (not to be excessively
% trusted...) Buoyancy-flux formulation
P_inner = Plc + Pic + Pgc; %disspation @ the inner= dissipation associated with grav. energy & latent heat from solidification
P_outer = Prc + Psc - Pk; %dissipation @ outer core= dissipation associated with rad., secular cooling, & sink associated w/ kc
Vc = 4/3*pi*(rc-ri)^3; % Volume of the core (m^3)
rho_av = Mc/Vc;        % Average density of the core (kg/m^3)
D = rc-ri;             % Thickness of the outer core (m)
phi_outer = rc*g/2;         % Gravitational potential at the CMB (m^2/s^2)
phi_inner = ri^2/rc * g/2;  % Gravitational potential at the inner core boundary (m^2/s^2)
phi_mean = 0.3*g/rc*((rc^5-ri^5)/(rc^3-ri^3)); % Average grav. potential in the outer core (m^2/s^2)

P_total = P_inner + P_outer;

if P_total > 0
    F_inner = P_inner/(phi_mean - phi_inner); %the potential energy at the inner core= grav. energy & latent heat from solidification
    F_outer = P_outer/(phi_outer - phi_mean);
    f_rat = F_inner/(F_outer+F_inner);

    powB = (P_inner + P_outer)/(Vc*rho_av*OmegaE^3*D^2);

    b_dip = 7.3*(1-ri/rc)*(1+f_rat);
    B_rms = powB^0.34*sqrt(rho_av*mu0)*OmegaE*D; 

    TDM = max(4*pi*rc^3/(sqrt(2)*mu0) * B_rms/b_dip, 0); %TDM=dipole moment
else
    TDM = 0;
end

Bs_core = mu0*TDM/(4*pi*rp^3); %magnetic field strength at surface of the moon. @ the equator


%% Scheinberg et al. 2018 Scaling Laws:

% Variables defined elsewhere in the code
% OmegaE = Angular rotation rate of the Moon: 2*pi/(27.3*24*3600)
% rho0 = 6477; density of the outer core/central density (kg/m3) 
% rc = 350e3; Radius of the core (m)

%Define variables
mu_0 = 4*pi*1e-7; %magnetic permeability
c = 0.63; % constant of proportionality
phi = P_total/Vc; %W/m3 Total dissipation/vol. of the core
rB = rc/rp; % the ratio of the top of the dynamo region to the planet radius.

% added an if statement since some dissipation terms are (-).
if phi <0 
    
    phi = 0;
end

%%Rossby # for buoyancy flux
%Ro_BF= U/(2*OmegaE*rc);
%rc=core radius=characteristic length scale of the flow

%Assuming a balance between inertial and buoyancy forces. Mixing length theory (ML)
B_ML = sqrt(2*mu_0*c*((rho0*rc^2*phi^2)^1/3))*rB^3/7;
U_ML =((phi*rc)/rho0)^(1/3); % fluid velocity ML (m/s)
Ro_ML = U_ML*((2*OmegaE*rc)^(-1)); %ML rossby #
mag_ML = mu0*1e5*U_ML*rc; %magnetic reynolds number

%Assuming a balance of Coriolis, inertial, and gravitational (Archimides)
%forces: CIA scaling law

B_CIA = sqrt(2*mu_0*c*((rho0^2*rc^4*OmegaE*phi^3)^1/5))*rB^3/7;
U_CIA =((phi/rho0)^(2/5))*((rc/OmegaE)^(1/5)); % fluid velocity CIA (m/s)
lengthscale=sqrt((U_CIA*rc)/OmegaE); %characteristic length scale of the flow for CIA
Ro_CIA = sqrt(U_CIA)*((4*OmegaE*rc)^(-1/2)); %CIA rossby #
mag_CIA = mu0*1e5*U_CIA*rc; %magnetic reynolds number

%CIA scaling was iffy b/c in some cases U_CIA was 0 & I was taking
%U_CIA/sqrt(U_CIA). This is why I was getting (-) numbers.

%Assuming a balance between Lorentz, gravitational and Coriolis forces: 
%(Magneto-Archimedes Coriolis) MAC scaling.

B_MAC = sqrt(2*mu_0*c*((rho0*rc^2*OmegaE*phi)^1/2))*rB^3/7;
U_MAC =(phi/(rho0*OmegaE))^(1/2); % fluid velocity MAC (m/s)
Ro_MAC = U_MAC*((2*OmegaE*rc)^(-1)); %MAC rossby #
mag_MAC = mu0*1e5*U_MAC*rc; %magnetic reynolds number
%The reason why Ro_MAC & Ro_ML looked flat was b/c Ro_CIA was super large &
%they were all plotted together on same axis

% Surface field
% B_surf = 1/7*B_CIA*rB^3; Commented out b/c I added it to each scaling law


end 
