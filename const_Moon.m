%% Stuff you could play with
dTLdc = -2.5e3;     % Compositional dependence of liquidus temperature (K)
dTLdP = 3e-8;   % Pressure dependence of liquidus temperature (K/Pa)

%% Stuff you should probably leave constant
mu0 = 4e-7 * pi;  % Vacuum pearmeability (SI)
G = 6.67e-11;     % Gravitational constant (SI)
y2s = 3.156e7;    % Seconds per Earth year 
R = 8.3145;       % Universal gas constant (J/K/mol)

rp = 1737e3;                     % Radius of the Moon (m)
OmegaE = 2*pi/(27.3*24*3600);    % Angular rotation rate (rad/s)

rc = 350e3;                         % Radius of the core (m)
K0 = 121.4e9;                       % Effective modulus (Pa) denotes the 
                                    % resistance to change of the material when a 
                                    % uniform pressure is applied. Given by the ratio 
                                    % of volume stress to volume strain
                                    
K1 = 5.7871; %6.2398                % Effective derivative of effective modulus 
rho0 = 6477;                        % Central density (kg/m^3)
Lp = sqrt(3*K0/(2*pi*G*rho0^2));    % Length scale (m)
Ap = 0.1*(5*K1-13);                 % Constant in density profile
P0 = 5.15e9;                         % Central pressure (Pa)
Mc = 4/3*pi*rho0*Lp^3*fc(rc/Lp, 0, Ap);                    % Mass of the core (kg)
g = 4/3*pi*G*rho0*rc*(1-0.6*(rc/Lp)^2 - 3/7*Ap*(rc/Lp)^4); % Gravitational acceleration near the CMB (m/s^2)
        
gamm = 1.65;      % Gruneisen parameter: Dimensionless combination of 
                  % bulk modulus, heat capacity at constant volume, thermal 
                  % expansion coefficient, and volume together. describes the 
                  % effect that changing the volume of a crystal lattice has 
                  % on its vibrational properties, and, as a consequence, 
                  % the effect that changing temperature has on the size or
                  % dynamics of the crystal lattice
    
Cc = 835;         % Specific heat of the core (J/kg/K)  
bet = 2.3;        % Coefficient of compositional expansion for inner core
DSc = 200;        % Entropy of melting (J/K/kg)
HlamC = 1.76e-17; % Average decay constant (1/s)
h0C = 4.1834e-14; % Heating per unit mass per ppm of K in the core (W/kg/ppm)
