%% Magnetic reynolds number
sigma = 6e-4; %electrical conductivity (S/m) Kuckes et al. 1971
mu_0 = 1.257e-6; %permeability of free space (H/m)
L =350e3; %characteristic length scale of the flow (m). Used radius of core (m)
Re_m = 50; %critical value for magnetic field generation

%Re_m = mu_0*sigma*U*L %magnetic reynolds number
U = Re_m/(mu_0*sigma*L) %characteristic flow velocity needed to produce Re_m=50 m/s
