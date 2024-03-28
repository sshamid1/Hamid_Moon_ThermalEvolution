c0 = 0.06;         % Mass fraction of sulfur 
P0 = 5.15;         % Pressure (GPa), ~5.15 GPa at the center, ~4.43 GPa at the CMB. Pressure at top and bottom to find changes in liquidus

dP = .5;      % "Small" pressure change
dc = 0.02;    % "Small" composition change

[TL0_BW, x] = calcTL_BW(P0,c0);  % Liquidus temperature at the reference pressure and composition

dTLdP_BW = (calcTL_BW(P0+dP,c0)-TL0_BW)/dP;   % Numerical derivative w.r.t. pressure
dTLdc_BW = (calcTL_BW(P0,c0+dc)-TL0_BW)/dc;   % Numerical derivative w.r.t. composition

TL0_H = calcTL_H(P0,c0);  % Liquidus temperature at the reference pressure and composition

dTLdP_H = (calcTL_H(P0+dP,c0)-TL0_H)/dP;   % Numerical derivative w.r.t. pressure
dTLdc_H = (calcTL_H(P0,c0+dc)-TL0_H)/dc;   % Numerical derivative w.r.t. composition


function [TL, x] = calcTL_BW(P,c0)

    % atomic mass of each element (Fe and S, g/mol)
    am_Fe = 55.845;
    am_S = 32.065;

    % Number of moles in 1 g of core alloy
    N_S = c0/am_S;         % Moles of S
    N_Fe = (1-c0)/am_Fe;   % Moles of Fe

    % Fe-FeS system 
    N_FeS = N_S;
    Nx_FeS = N_FeS/(N_FeS + (N_Fe - N_S));

    x = Nx_FeS; % The variable x in Buono & Walker (2011, GCA)
    
    % Transcribe from Buono & Walker (2011) eq 29
    TL = (-2.4724*P^4 + 28.025*P^3 + 9.1404*P^2 + 581.71*P + 3394.8)*x^4 ...
          + (1.7978*P^4 - 6.7881*P^3 - 197.69*P^2 - 271.69*P - 8219.5)*x^3 ...
          + (-0.1702*P^4 - 9.3959*P^3 + 163.53*P^2 - 319.35*P + 5698.6)*x^2 ...
          + (-0.2308*P^4 + 7.1*P^3 - 64.118*P^2 + 105.98*P - 1621.9)*x ...
          + (0.2302*P^4 - 5.3688*P^3 + 38.124*P^2 - 46.681*P + 1813.8);
    
end


function TL = calcTL_H(P_GPa,c0)

    P = 1e9*P_GPa;

    Tm0 = 1809;
    Tm1 = 1.99e-11;
    Tm2 = -1.84e-22;

    Teu0 = 1261;
    Teu1 = -6.91e-12;

    chi_eu0 = 0.207;
    chi_eu1 = 0;

    alpha_c = (Tm0*(1+Tm1*P+Tm2*P^2)-Teu0*(1+Teu1*P)) ... 
               / (Tm0*chi_eu0*(1+chi_eu1*P));
       
    TL = Tm0 * (1+Tm1*P+Tm2*P^2)*(1-alpha_c*c0);

end

%-------------------------------------------------------------------------------------
%% Updates to central pressure, central density, effective modulus, & the
%% effective derivative of the effective modulus: 

% We were initially using too high of a central pressure for the Moon. 
% We had 5.8 GPa, but Laneuville (2014) cited 5.1 GPa from Garcia et al. (2011)... 
% which actually states 5.3 GPa... in an erratum published in 2012... As a compromise, 
% Joe recalculated the structural parameters for a central pressure of 5.15 GPa, which
% still can reproduce the core mass and radius used in Scheinberg et al. (2018). 

%% Liquidus Update:

% The solidus of the lower mantle is ~1700 K. So, if we use the Buono liquidus, then an 
% inner core tends to co-exist with a (partially) molten lower mantle. On the other hand, 
% the Hauck liquidus implies that the inner core does not nucleate until the lower mantle 
% is almost fully solidified. Of course, the sulfur content of the lunar core is also 
% unknown-the ~80 K difference between the liquidus formulations is equivalent to a ~2-5 wt% 
% difference in the sulfur content. In other words, models with the Hauck liquidus have the 
% same behavior as ones with the Buono liquidus if lower sulfur contents are used. I'm using 
%Buono & Walker parameterizations. The Hauck liquidus predicts that the
%final liquidus temperature is 300-400 K too low since we're flattening the
%derivative instead of taking into account the entire curve (liquidus
%temperature vs sulfur mass fraction). In this case, the model would be
%predicting that the core would grow faster than what it's supposed to. The
%core should be smaller with Buono % Walker liquidus.

% Hauck: TL0 = 1741 + dTLdc*(c0-0.06); dTLdc = -4.08e3; dTLdP = 3.2e-8; OR
% Buono: TL0 = 1862 + dTLdc*(c0-0.06); dTLdc = -1.63e3; dTLdP = 1.5e-8; 

