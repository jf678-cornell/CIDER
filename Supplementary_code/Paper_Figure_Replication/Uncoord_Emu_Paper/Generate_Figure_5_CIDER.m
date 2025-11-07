%% Intro
% fullclear
addAllPaths
addpath(genpath('Uncoord_Emu_Paper'))
p = [];
p.var = 'TREFHT';
p.year_for_comp = [2055 2069];
p.ensemble_numbers = [1 2 3];
p.units = '°C';
p.years_of_ss = [2050 2069];
p.wraparound = 0;
% Longitude and latitude
load get_lat_and_lon.mat
if p.wraparound ==1
lon = [lon;360];
end 
ww = cos(lat/180*pi);
p.ww = ww;
p.lat = lat;
p.lon = lon;
p.latbounds = [-inf inf];
p.lonbounds = [-inf inf];
p.latbounds = [-inf inf];
p.lonbounds = [-inf inf];
%% CO2
% load Variables_for_multilat_emulator_precipitation.mat
% load Variables_for_multilat_emulator_temperature.mat
load CESM_params.mat
param_T_all
% param_T_all(1,:) = param_T_all(2,:)
param_P_all
% param_P_all(1,:) = param_P_all(2,:)

% load Variables_for_multilat_emulator.mat
% param_AOD_all(:,3)=-param_AOD_all(:,3);
T_base = globalMean(pattern_T_base);
P_base = globalMean(pattern_P_base);
AOD_base = globalMean(pattern_AOD_base);

load CO2_concentrations.mat
CO2levels_2020_2100_ssp245 = CO2_SSP245(6:86);
CO2levels_2035_2070_ssp245 = CO2_SSP245(6+15:86-31);
CO2levels_2035_2100_ssp245 = CO2_SSP245(6+15:86);
CO2levels_2035_2075_ssp245 = CO2_SSP245(6+15:86-25);

CO2levels_2020_2100_ssp126 = CO2_SSP126(6:86);
CO2levels_2035_2070_ssp126 = CO2_SSP126(6+15:86-31);
CO2levels_2035_2100_ssp126 = CO2_SSP126(6+15:86);

CO2levels_2035_2100_constant = zeros(1000,1)+ CO2_SSP245(6);
CO2levels_2035_2100_held = [CO2levels_2035_2100_ssp245;zeros(1000,1)];
CO2levels_2035_2100_held(21:end) = CO2levels_2035_2100_held(20);


CO2_ref = CO2levels_2035_2070_ssp245(1);
CO2_ref = CO2_SSP245(6+14);
CO2_forcing_SSP245 = 5.35*log((CO2levels_2035_2075_ssp245)/CO2_ref);
CO2_forcing_SSP245_month = repeatElements(CO2_forcing_SSP245,12);
years_of_all_forcing = [2005	2010	2020	2030	2040	2050	2060	2070	2080	2090	2100];
all_forcing =          [1.871	2.137	2.622	3.017	3.470	3.922	4.395	4.897	5.421	5.983	6.561];
CO2_forcing_SSP245_month = interp1(years_of_all_forcing,all_forcing,annualToMonthly(2035:2075))-interp1(years_of_all_forcing,all_forcing,2035);

%%
% uncoord_AOD_data_1 = processRun('AODVISstdn',1,'GAUSS-UNCOORD','203501-203612.nc',[2035 2035],[1],p);
% uncoord_AOD_data_2 = processRun('AODVISstdn',1,'GAUSS-UNCOORD','203512-204512.nc',[2035 2069],[1],p);
% uncoord_AOD_data_3 = processRun('AODVISstdn',1,'GAUSS-UNCOORD','204412-205312.nc',[2035 2069],[1],p);
uncoord_AOD_data = processRun('AODVISstdn',1,'GAUSS-UNCOORD','203501-207512.nc',[2035 2075],[1 2 3],p);

uncoord_T_data = processRun('TREFHT',1,'GAUSS-UNCOORD','203501-207512.nc',[2035 2075],[1 2 3],p);

uncoord_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-UNCOORD','203501-207512.nc',[2035 2075],[1 2 3],p);

%%

% uncoord_P_data_1 = processRun('PRECT',1000*24*60*60,'GAUSS-UNCOORD','203501-203612.nc',[2035 2035],[1],p);
% uncoord_P_data_2 = processRun('PRECT',1000*24*60*60,'GAUSS-UNCOORD','203512-204512.nc',[2035 2069],[1],p);
% uncoord_P_data_3 = processRun('PRECT',1000*24*60*60,'GAUSS-UNCOORD','204412-205312.nc',[2035 2069],[1],p);
% spliced_P = [uncoord_P_data_1.ensemble_monthly_average(1:12);uncoord_P_data_2.ensemble_monthly_average(2:end);uncoord_P_data_3.ensemble_monthly_average(14:end)];


% default_T_data = processRun('TREFHT',1,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3],p);
% monthlies_T = getMonthlies(default_T_data,2050,2069);
% deseasonalize_T = @(monthly_data) reshape(reshape(monthly_data,[12 length(monthly_data)/12])-monthlies_T,[length(monthly_data) 1]);
% Uncoord_T_data_deseasonalized = deseasonalize_T(spliced_T-T_base);

% default_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3],p);
% monthlies_P = getMonthlies(default_P_data,2050,2069);
% deseasonalize_P = @(monthly_data) reshape(reshape(monthly_data,[12 length(monthly_data)/12])-monthlies_P,[length(monthly_data) 1]);
% Uncoord_P_data_deseasonalized = deseasonalize_P(spliced_P-P_base);

load yearly_injection_rate.mat
% feedback_inj_30 = inj_rate_arr(:,3);
% feedback_30_AOD_data = processRun('AODVISstdn',1,'GAUSS-30N_30S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
% feedback_30_T_data = processRun('TREFHT',1,'GAUSS-30N_30S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
% feedback_30_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-30N_30S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
% injection_feedback_30_monthly = [zeros(420,1),repeatElements(feedback_inj_30,12)/12/2,zeros(420,3),repeatElements(feedback_inj_30,12)/12/2,zeros(420,1)];
% T_emulated_3030 = pattern_scale_emulate(param_AOD_all,param_T_all,injection_feedback_30_monthly,CO2_forcing_SSP245_month,all_T_patterns_scaled,T_base_pattern);
% P_emulated_3030 = pattern_scale_emulate(param_AOD_all,param_P_all,injection_feedback_30_monthly,CO2_forcing_SSP245_month,all_P_patterns_scaled,P_base_pattern);


PIT = getPreIndustrialTModel();
T_base_above_PI = T_base-PIT;
%%
Actor_1_ann    = [1 2 3 4 5 6 7 8 7 6 5 4 3 2  1  0  0  0 0 4 4 4 4 4 4 4 4 3 2 1 0 0 0 0 0 0 0 0 0 0 0]';
Actor_2_ann    = [0 0 0 0 0 1 2 3 4 5 6 7 8 9 10 11 12 13 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
Actor_3_ann    = [0 0 0 0 0 0 0 1 2 3 4 5 6 7  8  8  8  8 8 8 8 8 8 8 8 8 8 8 8 8 8 8.12 8.24 8.36 8.48 8.60 8.72 8.84 8.96 9.08 9.2]'; 
Actor_1_ann_30 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0  0  0  0  0 0 0 0 0 0 1 2 3 4 5 6 7 8 8.12 8.24 8.36 8.48 8.60 8.72 8.84 8.96 9.08 9.2]'; 

% Actor_1_ann = [1 2 3 4 5 6 7 8 7 6 5 4 3 2 1 0 0 0 ]';
% Actor_2_ann = [0 0 0 0 0 1 2 3 4 5 6 7 8 9 10 11 13 0 ]';
% Actor_3_ann = [0 0 0 0 0 0 0 1 2 3 4 5 6 7 8 8 8 8 ]'; 

Actor_1 = [];
Actor_1_30 = [];
Actor_2 = [];
Actor_3 = [];
for i = 1:length(Actor_1_ann)
    Actor_1 =  [Actor_1; [0 0 Actor_1_ann(i)/3 Actor_1_ann(i)/3 Actor_1_ann(i)/3 0 0 0 0 0 0 0]'];
    Actor_1_30 = [Actor_1_30; ones(12,1)*Actor_1_ann_30(i)/12];
    Actor_2 = [Actor_2; ones(12,1)*Actor_2_ann(i)/12];
    Actor_3 = [Actor_3; ones(12,1)*Actor_3_ann(i)/12];
end
empty_set = 0*Actor_1;
injection_tseries = [Actor_1,Actor_1_30,empty_set,Actor_2,empty_set,Actor_3,empty_set];

%%
% T_emulated = pattern_scale_emulate(param_AOD_all,param_T_all,injection_tseries,CO2_forcing_SSP245_month,all_T_patterns_scaled,T_base_pattern);
% T_emulated_SSP245 = pattern_scale_emulate(param_AOD_all,param_T_all,0*injection_tseries,CO2_forcing_SSP245_month,all_T_patterns_scaled,T_base_pattern);
% P_emulated = pattern_scale_emulate(param_AOD_all,param_P_all,injection_tseries,CO2_forcing_SSP245_month,all_P_patterns_scaled,P_base_pattern);
% Uncoord_T_emu = globalMean(T_emulated.sum_change,p);
% Uncoord_P_emu = globalMean(P_emulated.sum_change,p);

Uncoord_AOD_emu = globalMean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,injection_tseries,pattern_AOD_all));
Uncoord_T_emu = CIDER_response_from_all_injections_and_CO2([injection_tseries CO2_forcing_SSP245_month],param_AOD_all,param_T_all);
Uncoord_T_emu_pattern = CIDER_pattern_from_all_injections_and_CO2([injection_tseries CO2_forcing_SSP245_month],param_AOD_all,param_T_all,pattern_T_all)+pattern_T_base;
Uncoord_P_emu = CIDER_response_from_all_injections_and_CO2([injection_tseries CO2_forcing_SSP245_month],param_AOD_all,param_P_all);
%%
colors2 = colormap(brewermap([],'Set1'));

% figure
% end_year = 2075;
% tiley = tiledlayout(1,3);
% set(gcf, 'Position', [200, 200, 1800,1200]) % Set figure size
% tiley.Padding = 'tight';
% my_tile = nexttile();
% hold on
% plot(2035:1/12:end_year+11/12, uncoord_AOD_data.ensemble_monthly_average,'-b','LineWidth',2)
% % plot(2035:1/12:2069+11/12, uncoord_AOD_data_2.ensemble_monthly_average,'-c','LineWidth',2)
% % plot(2035:1/12:end_year+11/12, uncoord_AOD_data_3.ensemble_monthly_average,'-y','LineWidth',2)
% % plot(2035:1/12:end_year+11/12,sum(emulate_all_injection_array_AOD(param_AOD_all,param_T_all,injection_tseries,CO2_forcing_SSP245_month),2)+AOD_base,':r','LineWidth',4)
% % plot(2035:1/12:2056+11/12, (uncoord_AOD_data.ensemble_monthly_average(1:22*12)+uncoord_AOD_data_2.ensemble_monthly_average(1:22*12)+uncoord_AOD_data_3.ensemble_monthly_average)/3,'-k','LineWidth',4)
% plot(2035:1/12:end_year+11/12, Uncoord_AOD_emu+AOD_base,'-r','LineWidth',4)
% plot(2035:1/12:end_year+11/12, uncoord_AOD_data.individual_monthly_average,'-c','LineWidth',1)
% plot(2035:1/12:end_year+11/12, Uncoord_AOD_emu+AOD_base,'-r','LineWidth',4)
% plot(2035:1/12:end_year+11/12, uncoord_AOD_data.ensemble_monthly_average,'-b','LineWidth',2)
% 
% % plot(2035:1/12:2069+11/12,feedback_30_AOD_data.individual_monthly_average(:,1),'-k','LineWidth',4)
% % plot(2035:1/12:2069+11/12,sum(emulate_all_injection_array_AOD(param_AOD_all,param_T_all,injection_feedback_30_monthly,CO2_forcing_SSP245_month),2)+AOD_base,':k','LineWidth',4)
% title('UNCOORD AOD')
% box on
% grid on
% xlabel("Year")
% ylabel("AOD above baseline")
% legend("Simulation Ensemble Mean","CIDER","Ensemble Members")
% xlim([2035 2075])
% 
% set(gca,"FontSize",20)
% my_tile = nexttile();
% hold on
% plot(2035:end_year, uncoord_T_data.ensemble_annual_average-PIT,'-b','LineWidth',2)
% % plot(2035:2069, uncoord_T_data_2.individual_annual_average-PIT,'-c','LineWidth',2)
% % plot(2035:end_year, uncoord_T_data_3.individual_annual_average-PIT,'-y','LineWidth',2)
% % plot(2035:2069, (uncoord_T_data.individual_annual_average(1:22)+uncoord_T_data_2.individual_annual_average(1:22)+uncoord_T_data_3.individual_annual_average)/3-PIT,'-k','LineWidth',4)
% plot(2035:end_year, averageEvery2d(12,1,Uncoord_T_emu)+T_base-PIT,'-r','LineWidth',4)
% plot(2035:end_year, uncoord_T_data.individual_annual_average-PIT,'-c','LineWidth',1)
% plot(2035:end_year, averageEvery2d(12,1,Uncoord_T_emu)+T_base-PIT,'-r','LineWidth',4)
% plot(2035:end_year, uncoord_T_data.ensemble_annual_average-PIT,'-b','LineWidth',2)
% 
% % plot(2035:2069, feedback_30_T_data.individual_annual_average(:,1)-PIT,'-k','LineWidth',4)
% % plot(2035:1/12:2069+11/12, globalMean(T_emulated_3030.sum_change,p)-PIT,':k','LineWidth',4)
% title('UNCOORD Temperature')
% box on
% grid on
% xlabel("Year")
% ylabel("T above PI (°C)")
% xlim([2035 2075])
% legend("Simulation Ensemble Mean","CIDER","Ensemble Members")
% 
% set(gca,"FontSize",20)
% my_tile = nexttile();
% hold on
% plot(2035:end_year, uncoord_P_data.ensemble_annual_average,'-b','LineWidth',2)
% % plot(2035:2069, uncoord_P_data_2.individual_annual_average,'-c','LineWidth',2)
% % plot(2035:end_year, uncoord_P_data_3.individual_annual_average,'-y','LineWidth',2)
% % plot(2035:2056, (uncoord_P_data.individual_annual_average(1:22)+uncoord_P_data_2.individual_annual_average(1:22)+uncoord_P_data_3.individual_annual_average)/3,'-k','LineWidth',4)
% plot(2035:end_year, averageEvery2d(12,1,Uncoord_P_emu)+P_base,'-r','LineWidth',4)
% plot(2035:end_year, uncoord_P_data.individual_annual_average,'-c','LineWidth',1)
% plot(2035:end_year, averageEvery2d(12,1,Uncoord_P_emu)+P_base,'-r','LineWidth',4)
% plot(2035:end_year, uncoord_P_data.ensemble_annual_average,'-b','LineWidth',2)
% 
% % plot(2035:2069, feedback_30_P_data.individual_annual_average(:,1),'-k','LineWidth',4)
% % plot(2035:1/12:2069+11/12, globalMean(P_emulated_3030.sum_change,p),':k','LineWidth',4)
% title('UNCOORD Precipitation')
% box on
% grid on
% xlabel("Year")
% ylabel("Precipitation")
% legend("Simulation Ensemble Mean","CIDER","Ensemble Members")
% xlim([2035 2075])
% 
% set(gca,"FontSize",20)
% set(gcf,'renderer','painters')
% print(gcf,'-dpng',["Uncoord_Emu_Paper/Uncoord_Plots/Figure_UNCOORD_" + getNow() + ".png"],'-r300')
% close all
%%
colors = colormap(brewermap([],'Paired'));
colors2 = colormap(brewermap([],'Set1'));
colors3 = colormap(brewermap([],'Set2'));
figure
end_year = 2075;
tiley = tiledlayout(2,3);
set(gcf, 'Position', [200, 200, 1800,1200]) % Set figure size
tiley.Padding = 'tight';
my_tile = nexttile();
hold on
plot(2035:end_year,averageEvery2d(12,1,injection_tseries(:,1)*12),'LineWidth',4,"Color",colors(2,:))
plot(2035:end_year,averageEvery2d(12,1,injection_tseries(:,2)*12),'LineWidth',4,"Color",colors(1,:))
plot(2035:end_year,averageEvery2d(12,1,injection_tseries(:,4)*12),'LineWidth',4,"Color",colors(4,:))
plot(2035:end_year,averageEvery2d(12,1,injection_tseries(:,6)*12),'LineWidth',4,"Color",colors(6,:))
legend("Actor 1 @ 60°N","Actor 1 @ 30°N","Actor 2 @ 0°N","Actor 3 @ 30°S")
grid on
box on
set(gca,"FontSize",20)
title("(a)")
xlabel("Year")
ylabel("Injection Rates (Tg SO2/yr)")
xlim([2035 2075])
my_tile.TitleHorizontalAlignment = 'left';

my_tile = nexttile();
hold on
plot(2035:1/12:end_year+11/12, uncoord_AOD_data.ensemble_monthly_average,'-b','LineWidth',2)
% plot(2035:1/12:2069+11/12, uncoord_AOD_data_2.ensemble_monthly_average,'-c','LineWidth',2)
% plot(2035:1/12:end_year+11/12, uncoord_AOD_data_3.ensemble_monthly_average,'-y','LineWidth',2)
% plot(2035:1/12:end_year+11/12,sum(emulate_all_injection_array_AOD(param_AOD_all,param_T_all,injection_tseries,CO2_forcing_SSP245_month),2)+AOD_base,':r','LineWidth',4)
% plot(2035:1/12:2056+11/12, (uncoord_AOD_data.ensemble_monthly_average(1:22*12)+uncoord_AOD_data_2.ensemble_monthly_average(1:22*12)+uncoord_AOD_data_3.ensemble_monthly_average)/3,'-k','LineWidth',4)
plot(2035:1/12:end_year+11/12, Uncoord_AOD_emu+AOD_base,'-r','LineWidth',4)
plot(2035:1/12:end_year+11/12, uncoord_AOD_data.individual_monthly_average,'-c','LineWidth',1)
plot(2035:1/12:end_year+11/12, Uncoord_AOD_emu+AOD_base,'-r','LineWidth',4)
plot(2035:1/12:end_year+11/12, uncoord_AOD_data.ensemble_monthly_average,'-b','LineWidth',2)

% plot(2035:1/12:2069+11/12,feedback_30_AOD_data.individual_monthly_average(:,1),'-k','LineWidth',4)
% plot(2035:1/12:2069+11/12,sum(emulate_all_injection_array_AOD(param_AOD_all,param_T_all,injection_feedback_30_monthly,CO2_forcing_SSP245_month),2)+AOD_base,':k','LineWidth',4)
title('(b)')
box on
grid on
xlabel("Year")
ylabel("Global Mean AOD")
legend("CESM2 Ensemble Mean","CIDER","CESM2 Ensemble Members","Location","se")
xlim([2035 2075])
my_tile.TitleHorizontalAlignment = 'left';

set(gca,"FontSize",20)
my_tile = nexttile();
hold on
plot(2035:end_year, uncoord_T_data.ensemble_annual_average-PIT,'-b','LineWidth',2)
% plot(2035:2069, uncoord_T_data_2.individual_annual_average-PIT,'-c','LineWidth',2)
% plot(2035:end_year, uncoord_T_data_3.individual_annual_average-PIT,'-y','LineWidth',2)
% plot(2035:2069, (uncoord_T_data.individual_annual_average(1:22)+uncoord_T_data_2.individual_annual_average(1:22)+uncoord_T_data_3.individual_annual_average)/3-PIT,'-k','LineWidth',4)
plot(2035:end_year, averageEvery2d(12,1,Uncoord_T_emu)+T_base-PIT,'-r','LineWidth',4)
plot(2035:end_year, uncoord_T_data.individual_annual_average-PIT,'-c','LineWidth',1)
plot(2035:end_year, averageEvery2d(12,1,Uncoord_T_emu)+T_base-PIT,'-r','LineWidth',4)
plot(2035:end_year, uncoord_T_data.ensemble_annual_average-PIT,'-b','LineWidth',2)
line([2040 2058],[.5 .5],"LineStyle","--","Color",colors(4,:),"LineWidth",4)
line([2058 2075],[1 1],"LineStyle","--","Color","m","LineWidth",4)
% plot(2035:2069, feedback_30_T_data.individual_annual_average(:,1)-PIT,'-k','LineWidth',4)
% plot(2035:1/12:2069+11/12, globalMean(T_emulated_3030.sum_change,p)-PIT,':k','LineWidth',4)
title('(c)')
box on
grid on
xlabel("Year")
ylabel("Global Mean Temperature above PI (°C)")
xlim([2035 2075])
ylim([0.4 1.9])
% legend("CESM2 Ensemble Mean","CIDER","CESM2 Ensemble Members","Location","se")
my_tile.TitleHorizontalAlignment = 'left';

set(gca,"FontSize",20)
my_tile = nexttile();
hold on
plot(2035:end_year, regionalMean(uncoord_T_data.ensemble_annual_matrix,[60 90],[0 360])-273.15,'-b','LineWidth',2)
% plot(2035:2069, uncoord_T_data_2.individual_annual_average-PIT,'-c','LineWidth',2)
% plot(2035:end_year, uncoord_T_data_3.individual_annual_average-PIT,'-y','LineWidth',2)
% plot(2035:2069, (uncoord_T_data.individual_annual_average(1:22)+uncoord_T_data_2.individual_annual_average(1:22)+uncoord_T_data_3.individual_annual_average)/3-PIT,'-k','LineWidth',4)
plot(2035:end_year, regionalMean(averageEvery(12,1,Uncoord_T_emu_pattern),[60 90],[0 360])-273.15,'-r','LineWidth',4)
plot(2035:end_year,  regionalMean(uncoord_T_data.individual_annual_matrix,[60 90],[0 360])-273.15,'-c','LineWidth',1)
plot(2035:end_year, regionalMean(averageEvery(12,1,Uncoord_T_emu_pattern),[60 90],[0 360])-273.15,'-r','LineWidth',4)
plot(2035:end_year, regionalMean(uncoord_T_data.ensemble_annual_matrix,[60 90],[0 360])-273.15,'-b','LineWidth',2)
line([2035 2058],[-7.9 -7.9],"LineStyle","--","Color",colors(2,:),"LineWidth",4)

% plot(2035:2069, feedback_30_T_data.individual_annual_average(:,1)-PIT,'-k','LineWidth',4)
% plot(2035:1/12:2069+11/12, globalMean(T_emulated_3030.sum_change,p)-PIT,':k','LineWidth',4)
title('(d)')
box on
grid on
xlabel("Year")
ylabel("Arctic Mean Temperature (°C)")
xlim([2035 2075])
% legend("CESM2 Ensemble Mean","CIDER","CESM2 Ensemble Members","Location","se")
my_tile.TitleHorizontalAlignment = 'left';

set(gca,"FontSize",20)
my_tile =  nexttile();
hold on
plot(2035:end_year, calculateT1(uncoord_T_data.ensemble_annual_matrix),'-b','LineWidth',2)
% plot(2035:2069, uncoord_T_data_2.individual_annual_average-PIT,'-c','LineWidth',2)
% plot(2035:end_year, uncoord_T_data_3.individual_annual_average-PIT,'-y','LineWidth',2)
% plot(2035:2069, (uncoord_T_data.individual_annual_average(1:22)+uncoord_T_data_2.individual_annual_average(1:22)+uncoord_T_data_3.individual_annual_average)/3-PIT,'-k','LineWidth',4)
plot(2035:end_year, calculateT1(averageEvery(12,1,Uncoord_T_emu_pattern)),'-r','LineWidth',4)
for i = 1:3
    plot(2035:end_year,  calculateT1(uncoord_T_data.individual_annual_matrix(:,:,:,i)),'-c','LineWidth',1)
end
plot(2035:end_year, calculateT1(averageEvery(12,1,Uncoord_T_emu_pattern)),'-r','LineWidth',4)
plot(2035:end_year, calculateT1(uncoord_T_data.ensemble_annual_matrix),'-b','LineWidth',2)
line([2042 2058],[.9 .9],"LineStyle","--","Color",colors(6,:),"LineWidth",4)

% plot(2035:2069, feedback_30_T_data.individual_annual_average(:,1)-PIT,'-k','LineWidth',4)
% plot(2035:1/12:2069+11/12, globalMean(T_emulated_3030.sum_change,p)-PIT,':k','LineWidth',4)
title('(e)')
box on
grid on
xlabel("Year")
ylabel("N-S Temperature Gradient (°C)")
xlim([2035 2075])
% legend("CESM2 Ensemble Mean","CIDER","CESM2 Ensemble Members","Location","se")

set(gca,"FontSize",20)
my_tile.TitleHorizontalAlignment = 'left';


my_tile =  nexttile();
hold on
plot(2035:end_year, uncoord_P_data.ensemble_annual_average,'-b','LineWidth',2)
% plot(2035:2069, uncoord_P_data_2.individual_annual_average,'-c','LineWidth',2)
% plot(2035:end_year, uncoord_P_data_3.individual_annual_average,'-y','LineWidth',2)
% plot(2035:2056, (uncoord_P_data.individual_annual_average(1:22)+uncoord_P_data_2.individual_annual_average(1:22)+uncoord_P_data_3.individual_annual_average)/3,'-k','LineWidth',4)
plot(2035:end_year, averageEvery2d(12,1,Uncoord_P_emu)+P_base,'-r','LineWidth',4)
plot(2035:end_year, uncoord_P_data.individual_annual_average,'-c','LineWidth',1)
plot(2035:end_year, averageEvery2d(12,1,Uncoord_P_emu)+P_base,'-r','LineWidth',4)
plot(2035:end_year, uncoord_P_data.ensemble_annual_average,'-b','LineWidth',2)

% plot(2035:2069, feedback_30_P_data.individual_annual_average(:,1),'-k','LineWidth',4)
% plot(2035:1/12:2069+11/12, globalMean(P_emulated_3030.sum_change,p),':k','LineWidth',4)
title('(f)')
box on
grid on
xlabel("Year")
ylabel("Precipitation (mm day^{-1})")
% legend("CESM2 Ensemble Mean","CIDER","CESM2 Ensemble Members","Location","se")
xlim([2035 2075])
my_tile.TitleHorizontalAlignment = 'left';

set(gca,"FontSize",20)
set(gcf,'renderer','painters')
print(gcf,'-dpng',["Uncoord_Emu_Paper/Uncoord_Plots/Figure_5_" + getNow() + ".png"],'-r300')

