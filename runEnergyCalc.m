function [TC, TS, TLi, Tdis, dTCdt, dridt, ...
          Qsc, Qrc, Qgc, Qlc, Qic, TDM, c_oc, ...
          Psc, Prc, Pgc, Plc, Pic, Pk, P_total, Bs_core] = ...
                    runEnergyCalc(t, QC, TC, ri, core_Kppm, c0)
    
const_Moon;                       % Load constants from other script
TL0 = 1800 + dTLdc*(c0-0.06);     % Liquidus temperature at center (K)
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
    dridTC = -(1/(dTLdP-dTadP)) * TLi/(TC*rhoi*gi);
        
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
    
TS = TLi*((1-(ri/Lp)^2-Ap*(ri/Lp)^4)^(-gamm) * ...         % Temperature associated with secular cooling
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

% Composition of the outer core
Mic = 4/3*pi*rho0*Lp^3*fc(ri/Lp, 0, Ap); % Mass of the inner core (kg)
c_oc = c0*Mc/(Mc-Mic);         % Abundance of sulfur in the outer core

% SCALING LAW: MAGNETIC FIELD STRENGTH (not to be excessively trusted...)
P_inner = Plc + Pic + Pgc;
P_outer = Prc + Psc + - Pk;
Vc = 4/3*pi*(rc-ri)^3; % Volume of the core (m^3)
rho_av = Mc/Vc;        % Average density of the core (kg/m^3)
D = rc-ri;             % Thickness of the outer core (m)
phi_outer = rc*g/2;         % Gravitational potential at the CMB (m^2/s^2)
phi_inner = ri^2/rc * g/2;  % Gravitational potential at the inner core boundary (m^2/s^2)
phi_mean = 0.3*g/rc*((rc^5-ri^5)/(rc^3-ri^3)); % Average grav. potential in the outer core (m^2/s^2)

P_total = P_inner + P_outer;

if P_total > 0
    F_inner = P_inner/(phi_mean - phi_inner);
    F_outer = P_outer/(phi_outer - phi_mean);
    f_rat = F_inner/(F_outer+F_inner);

    powB = (P_inner + P_outer)/(Vc*rho_av*OmegaE^3*D^2);

    b_dip = 7.3*(1-ri/rc)*(1+f_rat);
    B_rms = powB^0.34*sqrt(rho_av*mu0)*OmegaE*D;

    TDM = max(4*pi*rc^3/(sqrt(2)*mu0) * B_rms/b_dip, 0);
else
    TDM = 0;
end

Bs_core = mu0*TDM/(4*pi*rp^3);

end
