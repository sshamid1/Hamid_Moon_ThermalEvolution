%% Specify model and initial conditions 
modName = 'V19K50p54'; % Scheinberg et al. (2018)
c0 = 0.075;              % Mass fraction of sulfur in the bulk core
core_Kppm = 25;         % Abundance of potassium in the bulk core (ppm) <50 ppm
QC_now = 0;           % CMB heat flow today (W) (vary from 0-2e9 W)
kc = 40;                % Thermal conductivity (W/m/K) 10-50, try in steps of 10

%% Run the model

load([modName,'.mat'])

% Specify timesteps
dt_Myr = 0.5;       % Timestep in Myr
t_Gyr_end = 4.2;    % Desired length of simulation in Gyr
t_Gyr_extra = ((t_Gyr(end)+1e-3*dt_Myr):(1e-3*dt_Myr):t_Gyr_end)'; % Gyr. So that it goes to present day

t_Gyr_all = 0.3+[t_Gyr; t_Gyr_extra];  % All the timesteps (should be 8400)

y2s = 3.156e7;                % seconds per year
t_all = 1e9*y2s*t_Gyr_all;    % time in seconds
dt = t_all(2)-t_all(1);       % constant timestep in seconds

% Append linear decline in heat flow after BMO solidification
QC_all = [Q_CMB; (linspace(Q_CMB(end),QC_now,length(t_Gyr_extra)))'];

% Initialize core radii and temperature
rc = 350e3;    % Radius of the core (m)
ri = 0;        % Initial radius of the inner core (m)
TC = TC0;      % Initial temperature at the CMB (K)

% Model temperature
t_Gyr_mod = 0.3+t_Gyr; %t_gyr timesteps from scheinberg paper.
TC_mod = (TC0-1700)+T_BMO;
i_mod = length(T_BMO);

% Determine the relationship b/t TC and ri
N_ri = 1e3;
ris = linspace(10e3,rc-1,N_ri);
TC_ri = zeros(1,N_ri); 
for jj = 1:N_ri
    [TC, ~, ~, ~, ~, ~, TC_nuc, ...
         ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~, ~] = ...
             runEnergyCalc(0, 0, 0, ris(jj), 0, c0, 0);
    TC_ri(1,jj) = TC;
end
ris = [0, ris]; TC_ri = [TC_nuc, TC_ri];

% Create arrays to store model output
NN = length(QC_all);
TC_all = zeros(1,NN); TS_all = zeros(1,NN); TLi_all = zeros(1,NN);
Tdis_all = zeros(1,NN); ri_all = zeros(1,NN); c_oc_all=zeros(1,NN);

Qsc_all = zeros(1,NN); Qrc_all = zeros(1,NN); Qgc_all = zeros(1,NN);
Qlc_all = zeros(1,NN); Qic_all = zeros(1,NN); Qad_all = zeros(1,NN);
TDM_all = zeros(1,NN);

Psc_all = zeros(1,NN); Prc_all = zeros(1,NN); Pgc_all = zeros(1,NN);
Plc_all = zeros(1,NN); Pic_all = zeros(1,NN); Pk_all = zeros(1,NN);
P_total_all = zeros(1,NN); Bs_core_all = zeros(1,NN); 
B_ML_all = zeros(1,NN); B_CIA_all = zeros(1,NN); B_MAC_all = zeros(1,NN); 
Ro_ML_all = zeros(1,NN); Ro_CIA_all = zeros(1,NN); Ro_MAC_all = zeros(1,NN);
mag_ML_all = zeros(1,NN); mag_CIA_all = zeros(1,NN); mag_MAC_all = zeros(1,NN);

% Run the model
for ii=1:NN
    
    %Thermal evolution BEFORE the BMO solidifies
    
    if ii < i_mod
        TC = TC_mod(ii);
        if TC < TC_nuc
            ri = interp1(TC_ri,ris,TC);
        else
            ri = 0;
        end
        [~, TS, TLi, Tdis, ~, ~, ~, ...
             Qsc, Qrc, Qgc, Qlc, Qic, Qad, TDM, c_oc, ...
             Psc, Prc, Pgc, Plc, Pic, Pk, P_total, Bs_core, B_ML, B_CIA, ...
             B_MAC,Ro_ML,Ro_CIA,Ro_MAC,mag_ML,mag_CIA,mag_MAC] = ...
          runEnergyCalc(t_all(ii), QC_all(ii), TC, ri, core_Kppm, c0, kc);    
    else
   
  %Thermal evolution AFTER the BMO solidifies
  
    % Calculate instantaneous energetics
    [TC, TS, TLi, Tdis, dTCdt, dridt, TC_nuc, ...
          Qsc, Qrc, Qgc, Qlc, Qic, Qad, TDM, c_oc, ...
          Psc, Prc, Pgc, Plc, Pic, Pk, P_total, Bs_core, B_ML, B_CIA, ...
          B_MAC,Ro_ML,Ro_CIA,Ro_MAC,mag_ML,mag_CIA,mag_MAC] = ...
      runEnergyCalc(t_all(ii), QC_all(ii), TC, ri, core_Kppm, c0, kc);
    
    % Advance the temperature and inner core radius one timestep
    TC = TC - dTCdt*dt;
    ri = min(rc, ri + dridt*dt);   
    end
    
    % Store instantaneous energetics
    TC_all(1,ii) = TC; TS_all(1,ii) = TS; TLi_all(1,ii) = TLi;
    Tdis_all(1,ii) = Tdis; ri_all(1,ii) = ri; c_oc_all(1,ii) = c_oc;
    
    Qsc_all(1,ii) = Qsc; Qrc_all(1,ii) = Qrc; Qgc_all(1,ii) = Qgc;
    Qlc_all(1,ii) = Qlc; Qic_all(1,ii) = Qic; Qad_all(1,ii) = Qad;
    TDM_all(1,ii) = TDM;
    
    Psc_all(1,ii) = Psc; Prc_all(1,ii) = Prc; Pgc_all(1,ii) = Pgc;
    Plc_all(1,ii) = Plc; Pic_all(1,ii) = Pic; Pk_all(1,ii) = Pk;
    P_total_all(1,ii) = P_total; Bs_core_all(1,ii) = Bs_core; 
    B_ML_all(1,ii)= B_ML; B_CIA_all(1,ii) = B_CIA; B_MAC_all(1,ii) = B_MAC;
    Ro_ML_all(1,ii) = Ro_ML; Ro_CIA_all(1,ii) = Ro_CIA; Ro_MAC_all(1,ii) = Ro_MAC;
    mag_ML_all(1,ii) = mag_ML; mag_CIA_all(1,ii) = mag_CIA; mag_MAC_all(1,ii) = mag_MAC;
    %I had an issue w/ Ro_MAC being 0 b/c I forgot to finish statement. I
    %put Ro_MAC(1,ii)=Ro_MAC instead of Ro_MAC_all(1,ii) = Ro_MAC

end


%% Plots
t_Gyr = t_all/(1e9*y2s);   % Time in billion years for convenience
i_ic = find(~ri_all,1,'last'); 


% Colors!
red = [228,26,28]./255; blue = [55,126,184]./255; green = [77,175,74]./255;
purple = [152,78,163]./255; orange = [255,127,0]./255; yellow = [255,255,51]./255;
brown = [166,86,40]./255; pink = [247,129,191]./255; grey= [33,33,33]./255; 
black= [0,0,0]; light_blue =[0.3010 0.7450 0.9330]; muted_green=[0.4660 0.6740 0.1880];
burgundy= [0.6350 0.0780 0.1840]; burnt_orange= [0.8500 0.3250 0.0980]; dark_purple= [0.4940 0.1840 0.5560];
purple_grey=[170,183,201]./255;rich_blue=[60,102,166]./255;pale_red=[228,122,119]./255;burnt_yellow=[230,173,87]./255;

% Fonts, etc. 
lw = 2;
fn = 'Courier';
fs = 20;

figure(1)
clf

%% Temperatures
% subplot(3,2,1) %rows, columns, position
% hold on
% 
% plot(t_Gyr, TC_all, 'Color', grey, 'LineWidth', lw)
% %plot(t_Gyr, Tdis_all, 'Color', light_blue, 'LineWidth', lw, 'LineStyle', '-.')
% plot(t_Gyr, TLi_all, 'Color', burgundy, 'LineWidth', lw, 'LineStyle', '--')
% plot(t_Gyr, TS_all, 'Color', burnt_orange, 'LineWidth', lw, 'LineStyle', '-.')
% 
% plot(t_Gyr_mod,TC_mod, 'Color', grey, 'LineWidth', lw)
% 
% hl = legend('At CMB','Average','Center/ICB');
% legend boxoff
% 
% hx = xlabel('Time (Gyr)');
% hy = ylabel('Temperature (K)');
% 
% set([gca hx hy hl], 'FontName', fn, 'FontSize', fs)
% set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on', ...
%     'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02])
% 
% % %% Heat budget
% subplot(3,2,2)
% hold on
% 
% plot(t_Gyr, smooth(QC_all./1e9,100), 'Color', grey, 'LineWidth', lw)
% plot(t_Gyr, smooth(Qsc_all./1e9,100), 'Color', rich_blue, 'LineWidth', lw,'LineStyle', '-')
% plot(t_Gyr, Qrc_all./1e9, 'Color', purple_grey, 'LineWidth', lw,'LineStyle', '-')
% plot(t_Gyr, smooth(Qgc_all./1e9,100), 'Color', orange, 'LineWidth', lw,'LineStyle', '-')
% plot(t_Gyr, smooth(Qlc_all./1e9,100), 'Color', dark_purple, 'LineWidth', lw,'LineStyle', '-')
% plot(t_Gyr(1:i_ic), Qad_all(1:i_ic)./1e9, 'Color', red, 'LineWidth', lw, 'LineStyle', '--')
% 
% hx = xlabel('Time (Gyr)');
% hy = ylabel({'Heat Budget', 'in the Core (GW)'});
% 
% hl = legend({'Q_{CMB}','Q_S','Q_R','Q_G','Q_L','Q_{AD}'},'NumColumns',2, 'FontSize', 13);
% legend boxoff
% 
% set([gca hx hy hl], 'FontName', fn, 'FontSize', fs)
% set(gca,'YLim',[0,ceil(max(Qad_all(1:i_ic))./1e9)],'XLim',[0 4.5],...
%     'TickDir','out','XMinorTick','on', ...
%     'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02])
% 
% %% Inner core
% subplot(3,2,3)
% hold on
% plot([0 4.5], [350 350], 'k:', 'LineWidth', 2, 'LineStyle', '--')
% plot(t_Gyr, ri_all./1e3, 'Color', grey, 'LineWidth', lw)
% 
% hx = xlabel('Time (Gyr)');
% hy = ylabel('Inner Core Radius (km)');
% ht = text(0.1,325,'Outer Core Radius'); %remember that the outer core radius includes inner core radius. Outer core is 350 km which includes a smaller inner core
% 
% set([gca ht hx hy], 'FontName', fn, 'FontSize', fs)
% set(gca,'XLim',[0 4.5],'YLim',[0,350],'TickDir','out','XMinorTick','on', ...
%     'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02],'box','off')
% 
% %% Sulfur in outer core
% subplot(3,2,4)
% hold on
% plot([0 4.5], [24 24], 'k:', 'LineWidth', 2,'LineStyle', '--') %eutectic [S] composition Fig. 4 Breuer et al. 2015
% plot(t_Gyr, 1e2*c_oc_all, 'Color', grey, 'LineWidth', lw)
% 
% hx = xlabel('Time (Gyr)');
% hy = ylabel({'Sulfur Abundance', 'in the Outer Core (wt%)'});
% ht = text(0.1,23,'Eutectic Composition','FontSize',fs, 'FontName',fn);
% 
% set([gca hx hy], 'FontName', fn, 'FontSize', fs)
% set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on', ...
%    'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02],'box','off')
% 
% %% Dissipation budget
% subplot(3,2,5)
% hold on
% 
% %plot(t_Gyr, smooth(P_total_all./1e6,100), 'Color', grey, 'LineWidth',lw) %commented out b/c the smooth function averages 
%                                                                           %dissipation over time which can make it look 
%                                                                           %like there isn?t positive dissipation, but in 
%                                                                           %actualuality there is. It?s just smoothed out. 
% plot(t_Gyr,P_total_all./1e6, 'Color', grey, 'LineWidth', lw) 
% plot(t_Gyr, smooth(Psc_all./1e6,100), 'Color', rich_blue, 'LineWidth', lw,'LineStyle','-')
% plot(t_Gyr, Prc_all./1e6, 'Color', purple_grey, 'LineWidth', lw,'LineStyle','-')
% plot(t_Gyr, smooth(Pgc_all./1e6,100), 'Color',orange , 'LineWidth', lw,'LineStyle','-')
% plot(t_Gyr, smooth(Plc_all./1e6,100), 'Color',dark_purple, 'LineWidth', lw, 'LineStyle','-')
% plot(t_Gyr, -Pk_all./1e6, 'Color', burgundy, 'LineWidth', lw)
% 
% plot([0 4.5], [0 0], 'k-', 'LineWidth', 2)
% 
% hl = legend({'\Phi_{CMB}','\Phi_S','\Phi_R','\Phi_G','\Phi_L','\Phi_K'},'NumColumns',3,'FontSize', 18);
% legend boxoff
% 
% hx = xlabel('Time (Gyr)');
% hy = ylabel({'Dissipation Budget', 'in the Core (MW)'});
% 
% set([gca hx hy], 'FontName', fn, 'FontSize', fs)
% set(gca,'YLim',[min(-Pk_all./1e6),...
%     ceil(1.1*max(max(smooth(P_total_all,100),max(smooth(Psc_all,100))))./1e6)],...
%     'XLim',[0 4.5],...
%     'TickDir','out','XMinorTick','on', ...
%     'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02])
% 
% %% Magnetic field. Note: the BMO exists in our model, we're just not
% %plotting or recalculating the field, just temperatures and energy/entropy budgets.
% subplot(3,2,6)
% 
% hold on
% plot(t_Gyr, smooth(1e6*B_MAC_all,100), 'Color', purple_grey, 'LineWidth', lw)
% plot(t_Gyr, smooth(1e6*Bs_core_all,100), 'Color', rich_blue, 'LineWidth', lw)
% plot(t_Gyr, smooth(1e6*B_ML_all,100), 'Color', black, 'LineWidth', lw)
% plot(t_Gyr, smooth(1e6*B_CIA_all,100), 'Color', burgundy, 'LineWidth', lw)
%  
% % plot([0.25 0.94], [35 35], 'color', black,'LineWidth', 2) %'LineStyle', '-.') %High field epoch: Evans 2018,Scheinberg 2018, Mighani 2020
% %  plot([3.5 4.5], [35 35], 'color', red,'LineWidth', lw, 'LineStyle', '-') %no dynamo Mighani 2020.
%   plot([0.3 2.9], [4.1 4.1], 'color', black, 'LineWidth', 2, 'LineStyle','--') %nominal V19K50p54 ML to get a sense of the "minimum" field possible
% % plot([0.3 2.3], [3.7 3.7], 'color', red,'LineWidth', lw, 'LineStyle', '-.') %V19K50p36 ML
% % plot([0.3 1.9], [3.4 3.4], 'color', red,'LineWidth', lw, 'LineStyle', '-.') %V19K50p27 ML
% % plot([0.3 3.2], [6.5 6.5], 'color', red,'LineWidth', lw, 'LineStyle', '-.') %V19K25p54 ML
% % plot([0.3 1.5], [4.2 4.2], 'color', red,'LineWidth', lw, 'LineStyle', '-.') %V18.5K50p54 ML 
% % plot([0.3 2.4], [12 12], 'color', red,'LineWidth', lw, 'LineStyle', '-.') %V18K00p100 ML
%  
%  
% hl = legend({'MAC', 'BF','ML','CIA'},'NumColumns',1,'FontSize', fs,'NumColumns',2);
% legend boxoff
% 
% hx = xlabel('Time (Gyr)');
% hy = ylabel({'Surface Magnetic', 'Field Intensity (\mu{T})'});
% 
% ht = text(0.3, 2.5,'Nominal BMO-Hosted Dynamo', 'FontSize', 13);
% %ht = text(2.9, 7,'BMO-Hosted Dynamo (CIA)', 'FontSize', 10);
% %ht = text(2.9, 16,'BMO-Hosted Dynamo (MAC)', 'FontSize', 10);
% %ht = text(0.25, 17,'High-Field Epoch', 'FontSize', 14);
% 
% 
% set([gca hx hy], 'FontName', fn, 'FontSize', fs)
% set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on', ...
%     'YScale', 'log', 'YLim', [1e-3 50], 'YMinorTick','on','LineWidth',1,'TickLength',[.02 .02],'box','off')

%% Rossby number 
%  Ro_MAC_all(Ro_MAC_all < 1e-7) = 1e-7; %anything < 1e-7 = 1e-7
%  Ro_CIA_all(Ro_CIA_all < 1e-7) = 1e-7;
%  Ro_ML_all(Ro_ML_all < 1e-7) = 1e-7;
% % 
% % subplot(2,2,4)
% 
% hold on
% plot(t_Gyr,Ro_CIA_all, 'Color', burgundy, 'LineWidth', lw)
% plot(t_Gyr, Ro_ML_all, 'Color', black, 'LineWidth', lw)
% plot(t_Gyr,Ro_MAC_all, 'Color', purple_grey, 'LineWidth', lw)
% 
% hl = legend({'CIA','ML','MAC'},'NumColumns',1,'FontSize', 40,'NumColumns',1);
% legend boxoff
% 
% hx = xlabel('Time (Gyr)');
% hy = ylabel({'Local Rossby Number'});
% set([gca hx hy], 'FontName', fn, 'FontSize', 25);
% set(gca,'XLim',[0 4.5],'YLim',[1e-5 1e-1],'Yscale','log','TickDir','out','XMinorTick','on','YMinorTick','on','LineWidth',1.5,'TickLength',[.02 .02],'box','off')
% hold on
% 
% %When finding the Rossby # I accidently smoothed out the plot & multiplied everything by 100 when
%copying the format from the magnetic field plots, so that was corrected. 

%% Magnetic reynolds number

hold on
plot(t_Gyr,mag_ML_all, 'Color', burgundy, 'LineWidth', lw)
plot(t_Gyr,mag_CIA_all, 'Color', black, 'LineWidth', lw)
plot(t_Gyr,mag_MAC_all, 'Color', purple_grey, 'LineWidth', lw)

hl = legend({'ML','CIA','MAC'},'NumColumns',1,'FontSize', 40,'NumColumns',1);
legend boxoff

title('\sigma = 10^{5}');
hx = xlabel('Time (Gyr)');
hy = ylabel({'Magnetic Reynolds Number'});
set([gca hx hy], 'FontName', fn, 'FontSize', 25);
set(gca,'XLim',[0 4.5],'TickDir','out','XMinorTick','on','YMinorTick','on','LineWidth',1.5,'TickLength',[.02 .02],'box','off')
hold on
