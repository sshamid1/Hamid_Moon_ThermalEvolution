%% Run the model

% Specify initial conditions
rc = 325e3;     % Radius of the core (m)
c0 = 0.035;      % Mass fraction of sulfur in the core
core_Kppm = 20; % Abundance of potassium in the core (ppm)

TC = 1800;      % Initial temperature of the core (K)
ri = 0;         % Initial radius of the inner core (m)

% Specify number of timesteps
NN = 4500;
y2s = 3.154e7;

% ...specifically, run from 0 to 4.5 Gyr in 1 Myr steps
t_all = linspace(0, 4.5e9*y2s, NN);  
dt = t_all(2)-t_all(1);

% Assume linear decay of heat flow over time <-- OBVIOUSLY TOO SIMPLE?
QC_all = 1e9 * linspace(5,2,NN);

% Create arrays to store model output
TC_all = zeros(1,NN); TS_all = zeros(1,NN); TLi_all = zeros(1,NN);
Tdis_all = zeros(1,NN); ri_all = zeros(1,NN); c_oc_all=zeros(1,NN);

Qsc_all = zeros(1,NN); Qrc_all = zeros(1,NN); Qgc_all = zeros(1,NN);
Qlc_all = zeros(1,NN); Qic_all = zeros(1,NN); TDM_all = zeros(1,NN);

Psc_all = zeros(1,NN); Prc_all = zeros(1,NN); Pgc_all = zeros(1,NN);
Plc_all = zeros(1,NN); Pic_all = zeros(1,NN); Pk_all = zeros(1,NN);
P_total_all = zeros(1,NN); Bs_core_all = zeros(1,NN);

% Run the model
for ii=1:NN
    
    % Calculate instantaneous energetics
    [TC, TS, TLi, Tdis, dTCdt, dridt, ...
          Qsc, Qrc, Qgc, Qlc, Qic, TDM, c_oc, ...
          Psc, Prc, Pgc, Plc, Pic, Pk, P_total, Bs_core] = ...
               runEnergyCalc(t_all(ii), QC_all(ii), TC, ri, core_Kppm, c0);
   
    % Store instantaneous energetics
    TC_all(1,ii) = TC; TS_all(1,ii) = TS; TLi_all(1,ii) = TLi;
    Tdis_all(1,ii) = Tdis; ri_all(1,ii) = ri; c_oc_all(1,ii) = c_oc;
    
    Qsc_all(1,ii) = Qsc; Qrc_all(1,ii) = Qrc; Qgc_all(1,ii) = Qgc;
    Qlc_all(1,ii) = Qlc; Qic_all(1,ii) = Qic; TDM_all(1,ii) = TDM;
    
    Psc_all(1,ii) = Psc; Prc_all(1,ii) = Prc; Pgc_all(1,ii) = Pgc;
    Plc_all(1,ii) = Plc; Pic_all(1,ii) = Pic; Pk_all(1,ii) = Pk;
    P_total_all(1,ii) = P_total; Bs_core_all(1,ii) = Bs_core;
    
    % Advance the temperature and inner core radius one timestep
    TC = TC - dTCdt*dt;
    ri = min(rc, ri + dridt*dt);   
    
end


%% Plots
t_Gyr = t_all/(1e9*y2s);   % Time in billion years for convenience

% Colors!
red = [228,26,28]./255; blue = [55,126,184]./255; green = [77,175,74]./255;
purple = [152,78,163]./255; orange = [255,127,0]./255; yellow = [255,255,51]./255;
brown = [166,86,40]./255; pink = [247,129,191]./255; grey= [33,33,33]./255;

% Fonts, etc. 
lw = 3;
fn = 'Arial';
fs = 18;

figure(1)
clf

% Temperatures
subplot(3,2,1)
hold on

plot(t_Gyr, TC_all, 'Color', grey, 'LineWidth', lw)
plot(t_Gyr, Tdis_all, 'Color', blue, 'LineWidth', lw)
plot(t_Gyr, TLi_all, 'Color', purple, 'LineWidth', lw)

hl = legend('At CMB','Average','Center/ICB');
legend boxoff

hx = xlabel('Time (Gyr)');
hy = ylabel('Temperature (K)');

set([gca hx hy hl], 'FontName', fn, 'FontSize', fs)
set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on', ...
    'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02])

% Heat budget
subplot(3,2,2)
hold on

plot(t_Gyr, QC_all./1e6, 'Color', grey, 'LineWidth', lw)
plot(t_Gyr, Qsc_all./1e6, 'Color', blue, 'LineWidth', lw)
plot(t_Gyr, Qrc_all./1e6, 'Color', green, 'LineWidth', lw)
plot(t_Gyr, Qgc_all./1e6, 'Color', purple, 'LineWidth', lw)
plot(t_Gyr, Qlc_all./1e6, 'Color', brown, 'LineWidth', lw)

hx = xlabel('Time (Gyr)');
hy = ylabel('Heat budget (MW)');

hl = legend({'Q_C','Q_S','Q_R','Q_G','Q_L'},'NumColumns',3);
legend boxoff

set([gca hx hy hl], 'FontName', fn, 'FontSize', fs)
set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on', ...
    'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02])

% Inner core
subplot(3,2,3)
hold on
plot([0 4.5], [325 325], 'k:', 'LineWidth', 3)
plot(t_Gyr, ri_all./1e3, 'Color', grey, 'LineWidth', lw)

hx = xlabel('Time (Gyr)');
hy = ylabel('Radius of the inner core (km)');
ht = text(0.1,350,'Radius of the core');

set([gca ht hx hy], 'FontName', fn, 'FontSize', fs)
set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on', ...
    'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02],'box','off')

% Sulfur in outer core
subplot(3,2,4)
hold on
plot([0 4.5], [23 23], 'k:', 'LineWidth', 3)
plot(t_Gyr, 1e2*c_oc_all, 'Color', grey, 'LineWidth', lw)

hx = xlabel('Time (Gyr)');
hy = ylabel('Sulfur abundance (wt%)');
ht = text(0.1,22,'Eutectic composition');

set([gca ht hx hy], 'FontName', fn, 'FontSize', fs)
set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on', ...
    'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02],'box','off')

% Dissipation budget
subplot(3,2,5)
hold on

plot(t_Gyr, P_total_all./1e6, 'Color', grey, 'LineWidth', lw)
plot(t_Gyr, Psc_all./1e6, 'Color', blue, 'LineWidth', lw)
plot(t_Gyr, Prc_all./1e6, 'Color', green, 'LineWidth', lw)
plot(t_Gyr, Pgc_all./1e6, 'Color', purple, 'LineWidth', lw)
plot(t_Gyr, Plc_all./1e6, 'Color', brown, 'LineWidth', lw)
plot(t_Gyr, -Pk_all./1e6, 'Color', pink, 'LineWidth', lw)

plot([0 4.5], [0 0], 'k:', 'LineWidth', 2)

hl = legend({'P_C','P_S','P_R','P_G','P_L','P_K'},'NumColumns',2);
legend boxoff

hx = xlabel('Time (Gyr)');
hy = ylabel('Dissipation budget (MW)');

set([gca hx hy], 'FontName', fn, 'FontSize', fs)
set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on', ...
    'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02])


% Magnetic field
subplot(3,2,6)

plot(t_Gyr, 1e6*Bs_core_all, 'Color', grey, 'LineWidth', lw)

hx = xlabel('Time (Gyr)');
hy = ylabel('Surface field intensity (\mu{T})');

set([gca hx hy], 'FontName', fn, 'FontSize', fs)
set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on', ...
    'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02],'box','off')
