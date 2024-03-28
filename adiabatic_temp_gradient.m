%% Adiabatic thermal gradient in the core
rho0= 6477; %density at the center of the core (kg/m3) 
gamma = 1.65; %Grüneisen parameter for the core
Lp = 350e3; %length scale (m). Use radius of core instead of what's in appendix?
Ap = 1.59; %Constant in density profile
T0= 1745; %Center/ICB temperature (K). Used temp @ 4.5 Gyr from nominal model
r = 250e3; %inner core radius at reference T0 (m). 

rho_r= rho0*((1-((r/Lp)^2)) - (Ap*((r/Lp)^4))); %radial density (kg/m3)
Ta(r) = T0*((rho_r/rho0)^gamma); %Adiabatic thermal gradient in the core (K/GPa)
Ta(r)