% Specify initial conditions
rc = 325e3;     % Radius of the core (m)
c0 = 0.04;      % Mass fraction of sulfur in the core 
core_Kppm = 20; % Abundance of potassium in the core (ppm)

QC = 5e8;       % Core/mantle heat flow (W, note 1e9 W = 1000 MW)

% Note that TC is ignored if ri > 0! 
TC = 2000;      % Temperature at the core/mantle boundary (K)
ri = 0;         % Initial radius of the inner core (m)

t = 3e9*3.154e7; % Time (4 Gyr, semi-randomly)

% Calculate instantaneous energetics
[TC, TS, TLi, Tdis, dTCdt, dridt, ...
 Qsc, Qrc, Qgc, Qlc, Qic, TDM, c_oc, ...
 Psc, Prc, Pgc, Plc, Pic, Pk, P_total, Bs_core] = ...
     runEnergyCalc(t, QC, TC, ri, core_Kppm, c0);
   
% Output variables:
% TC, TS, TLi, Tdis = temperatures in the core
% dTCdt = cooling rate
% dridt = inner core growth rate
% Qsc, Qrc, Qgc, Qlc, Qic = heat budget
% Psc, Prc, Pgc, Plc, Pic, Pk, P_total = dissipation budget
% TDM = True Dipole Moment for the magnetic field
% Bs_core = Magnetic field intensity predicted at the surface