addAllPaths
close all
load UKESM_params.mat
%% Load/create 
injection_12 = ones(420,1);
injection_12_60N = [];
injection_12_60S = [];
for i = 1:35
injection_12_60N = [injection_12_60N; [0 0 4 4 4 0 0 0 0 0 0 0]'];
injection_12_60S = [injection_12_60S; [0 0 0 0 0 0 0 0 4 4 4 0]'];
end

load yearly_injection_rate.mat
% 0N_1.0, 15N_15S, 30N_30S, 60N_60S_1.0, Global+1C
% feedback_inj_0N = inj_rate_arr(:,1);
feedback_inj_0N = (readtable("ukesm_EQ_inj_1.log").x30S_Tg_(2:end)+readtable("ukesm_EQ_inj_2.log").x30S_Tg_(2:end)++readtable("ukesm_EQ_inj_3.log").x30S_Tg_(2:end))/3;
feedback_inj_15 = 2*(readtable("ukesm_15NS_inj_1.log").x30S_Tg_(2:end-1)+readtable("ukesm_15NS_inj_2.log").x30S_Tg_(2:end)+readtable("ukesm_15NS_inj_3.log").x30S_Tg_(2:end))/3;
feedback_inj_30 = 2*(readtable("ukesm_30NS_inj_1.log").x30S_Tg_(2:end)+readtable("ukesm_30NS_inj_2.log").x30S_Tg_(2:end)+readtable("ukesm_30NS_inj_3.log").x30S_Tg_(2:end))/3;
feedback_inj_60 = 2*(readtable("ukesm_60NS_inj_1.log").x30S_Tg_(2:end)+readtable("ukesm_60NS_inj_2.log").x30S_Tg_(2:end)+readtable("ukesm_60NS_inj_3.log").x30S_Tg_(2:end))/3;
ARISE_inj_30N = (readtable("ukesm_ARISE_inj_1.log").x30N_Tg_(2:end-1)+readtable("ukesm_ARISE_inj_2.log").x30N_Tg_(2:end-1)+readtable("ukesm_ARISE_inj_3.log").x30N_Tg_(2:end-1))/3;
ARISE_inj_15N = (readtable("ukesm_ARISE_inj_1.log").x15N_Tg_(2:end-1)+readtable("ukesm_ARISE_inj_2.log").x15N_Tg_(2:end-1)+readtable("ukesm_ARISE_inj_3.log").x15N_Tg_(2:end-1))/3;
ARISE_inj_15S = (readtable("ukesm_ARISE_inj_1.log").x15S_Tg_(2:end-1)+readtable("ukesm_ARISE_inj_2.log").x15S_Tg_(2:end-1)+readtable("ukesm_ARISE_inj_3.log").x15S_Tg_(2:end-1))/3;
ARISE_inj_30S = (readtable("ukesm_ARISE_inj_1.log").x30S_Tg_(2:end-1)+readtable("ukesm_ARISE_inj_2.log").x30S_Tg_(2:end-1)+readtable("ukesm_ARISE_inj_3.log").x30S_Tg_(2:end-1))/3;
feedback_inj_0N_monthly = repeatElements(feedback_inj_0N,12)/12;
feedback_inj_15_monthly = repeatElements(feedback_inj_15,12)/12/2;
feedback_inj_30_monthly = repeatElements(feedback_inj_30,12)/12/2;
feedback_inj_60_monthly_base = repeatElements(feedback_inj_60,12)/3/2;
feedback_inj_60N_monthly = injection_12_60N/4.*feedback_inj_60_monthly_base;
feedback_inj_60S_monthly = injection_12_60S/4.*feedback_inj_60_monthly_base;
feedback_inj_60NS_monthly = feedback_inj_60N_monthly+feedback_inj_60S_monthly;
% feedback_inj_60N_monthly = feedback_inj_60NS_monthly/2;
% feedback_inj_60S_monthly = feedback_inj_60NS_monthly/2;

all_step_injections = [injection_12_60N ones(420,5)];
all_feedback_injections = [feedback_inj_60N_monthly feedback_inj_30_monthly feedback_inj_15_monthly feedback_inj_0N_monthly feedback_inj_15_monthly feedback_inj_30_monthly];

injection_matrix_fdbk_0N = [zeros(420,3) feedback_inj_0N_monthly zeros(420,3)];
injection_matrix_fdbk_15NS = [zeros(420,2) feedback_inj_15_monthly zeros(420,1) feedback_inj_15_monthly zeros(420,2)];
injection_matrix_fdbk_30NS = [zeros(420,1) feedback_inj_30_monthly zeros(420,3) feedback_inj_30_monthly zeros(420,1)];
injection_matrix_fdbk_60NS = [zeros(420,0) feedback_inj_60N_monthly zeros(420,5) feedback_inj_60S_monthly zeros(420,0)];
injection_matrix_ARISE = [zeros(420,1) repeatElements(ARISE_inj_30N,12)/12 repeatElements(ARISE_inj_15N,12)/12 zeros(420,1) repeatElements(ARISE_inj_15S,12)/12 repeatElements(ARISE_inj_30S,12)/12 zeros(420,1)];

%% CO2
load CO2_concentrations.mat
CO2levels_2035_2070_ssp245 = CO2_SSP245(6+15:86-31);

CO2_ref = CO2_SSP245(6+14);
CO2_forcing_SSP245 = 5.35*log((CO2levels_2035_2070_ssp245)/CO2_ref);
CO2_forcing_SSP245_month = repeatElements(CO2_forcing_SSP245,12);
years_of_all_forcing = [2005	2010	2020	2030	2040	2050	2060	2070	2080	2090	2100];
all_forcing =          [1.871	2.137	2.622	3.017	3.470	3.922	4.395	4.897	5.421	5.983	6.561];
CO2_forcing_SSP245_month = interp1(years_of_all_forcing,all_forcing,annualToMonthly(2035:2069))-interp1(years_of_all_forcing,all_forcing,2035);

%% Load SP Data
AOD_0N_matrix  = ncread("u-ch622_AOD.nc","AOD550_delta");
AOD_15N_matrix = ncread("u-cd352_AOD.nc","AOD550_delta");
AOD_15S_matrix = ncread("u-cd353_AOD.nc","AOD550_delta");
AOD_30N_matrix = ncread("u-cd297_AOD.nc","AOD550_delta");
AOD_30S_matrix = ncread("u-cd354_AOD.nc","AOD550_delta");
AOD_60N_matrix = ncread("u-cm174_AOD.nc","AOD550_delta");
AOD_60S_matrix = ncread("u-cm175_AOD.nc","AOD550_delta");


T_0N_matrix = ncread("u-ch622_tas.nc","tas_delta");
T_15N_matrix = ncread("u-cd352_tas.nc","tas_delta");
T_15S_matrix = ncread("u-cd353_tas.nc","tas_delta");
T_30N_matrix = ncread("u-cd297_tas.nc","tas_delta");
T_30S_matrix = ncread("u-cd354_tas.nc","tas_delta");
T_60N_matrix = ncread("u-cm174_tas.nc","tas_delta");
T_60S_matrix = ncread("u-cm175_tas.nc","tas_delta");

T_SSP245_matrix = regridUKESMtoCESM((ncread("tas_ssp_ukesm_1.nc","tas")+ncread("tas_ssp_ukesm_2.nc","tas")+ncread("tas_ssp_ukesm_3.nc","tas"))/3);
T_SSP245_matrix_base = mean(T_SSP245_matrix(:,:,12*10+1:12*30),3);
T_SSP245_matrix = T_SSP245_matrix(:,:,12*20+1:12*55) - T_SSP245_matrix_base;

pdimchange = 24*60*60;
P_60N_matrix = pdimchange*ncread("u-cm174_pr.nc","pr_delta");
P_30N_matrix = pdimchange*ncread("u-cd297_pr.nc","pr_delta");
P_15N_matrix = pdimchange*ncread("u-cd352_pr.nc","pr_delta");
P_0N_matrix  = pdimchange*ncread("u-ch622_pr.nc","pr_delta");
P_15S_matrix = pdimchange*ncread("u-cd353_pr.nc","pr_delta");
P_30S_matrix = pdimchange*ncread("u-cd354_pr.nc","pr_delta");
P_60S_matrix = pdimchange*ncread("u-cm175_pr.nc","pr_delta");
P_SSP245_matrix = pdimchange*regridUKESMtoCESM((ncread("pr_ssp_ukesm_1.nc","pr")+ncread("pr_ssp_ukesm_2.nc","pr")+ncread("pr_ssp_ukesm_3.nc","pr")+ncread("pr_ssp_ukesm_4.nc","pr")+ncread("pr_ssp_ukesm_5.nc","pr"))/5);
P_SSP245_matrix_base = mean(P_SSP245_matrix(:,:,12*10+1:12*30),3);
P_SSP245_matrix = P_SSP245_matrix(:,:,12*20+1:12*55) - P_SSP245_matrix_base;

%% Parse SP data
AOD_0N = globalMean(AOD_0N_matrix);
AOD_15N = globalMean(AOD_15N_matrix);
AOD_15S = globalMean(AOD_15S_matrix);
AOD_30N = globalMean(AOD_30N_matrix);
AOD_30S = globalMean(AOD_30S_matrix);
AOD_60N = globalMean(AOD_60N_matrix);
AOD_60S = globalMean(AOD_60S_matrix);
all_AOD_step_responses = [AOD_60N AOD_30N AOD_15N AOD_0N AOD_15S AOD_30S];

T_0N = globalMean(T_0N_matrix);
T_15N = globalMean(T_15N_matrix);
T_15S = globalMean(T_15S_matrix);
T_30N = globalMean(T_30N_matrix);
T_30S = globalMean(T_30S_matrix);
T_60N = globalMean(T_60N_matrix);
T_60S = globalMean(T_60S_matrix);
all_T_step_responses = [T_60N T_30N T_15N T_0N T_15S T_30S];
T_base = globalMean(T_SSP245_matrix_base);
T_SSP245 = globalMean(T_SSP245_matrix);

P_0N = globalMean(P_0N_matrix);
P_15N = globalMean(P_15N_matrix);
P_15S = globalMean(P_15S_matrix);
P_30N = globalMean(P_30N_matrix);
P_30S = globalMean(P_30S_matrix);
P_60N = globalMean(P_60N_matrix);
P_60S = globalMean(P_60S_matrix);
all_P_step_responses = [P_60N P_30N P_15N P_0N P_15S P_30S];

P_SSP245 = globalMean(P_SSP245_matrix);
P_base = globalMean(P_SSP245_matrix_base);
%% Load feedback data
ukesm_SSP245_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_SSP245_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_SSP245_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_SSP245_strataod_3.nc","__xarray_dataarray_variable__"))/3);
ukesm_SSP245_strataod_matrix = ukesm_SSP245_strataod_matrix(:,:,12*20+1:12*55);
ukesm_60NS_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_60NS_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_60NS_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_60NS_strataod_3.nc","__xarray_dataarray_variable__"))/3);
ukesm_60NS_strataod_matrix = ukesm_60NS_strataod_matrix(:,:,25:end);
ukesm_30NS_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_30NS_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_30NS_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_30NS_strataod_3.nc","__xarray_dataarray_variable__"))/3);
ukesm_30NS_strataod_matrix = ukesm_30NS_strataod_matrix(:,:,25:end);
ukesm_15NS_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_15NS_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_15NS_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_15NS_strataod_3.nc","__xarray_dataarray_variable__"))/3);
ukesm_15NS_strataod_matrix = ukesm_15NS_strataod_matrix(:,:,25:end);
ukesm_EQ_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_EQ_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_EQ_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_EQ_strataod_3.nc","__xarray_dataarray_variable__"))/3);
ukesm_EQ_strataod_matrix = ukesm_EQ_strataod_matrix(:,:,25:end);
AOD_base = mean(globalMean(ukesm_SSP245_strataod_matrix));
AOD_SSP245 = globalMean(ukesm_SSP245_strataod_matrix);
AOD_60NS_fdbk = globalMean(ukesm_60NS_strataod_matrix-ukesm_SSP245_strataod_matrix);
AOD_30NS_fdbk = globalMean(ukesm_30NS_strataod_matrix-ukesm_SSP245_strataod_matrix);
AOD_15NS_fdbk = globalMean(ukesm_15NS_strataod_matrix-ukesm_SSP245_strataod_matrix);
AOD_EQ_fdbk = globalMean(ukesm_EQ_strataod_matrix-ukesm_SSP245_strataod_matrix);
all_AOD_fdbk_responses = [AOD_60NS_fdbk AOD_30NS_fdbk AOD_15NS_fdbk AOD_EQ_fdbk AOD_15NS_fdbk AOD_30NS_fdbk];

T_30NS_fdbk = globalMean(regridUKESMtoCESM((ncread("ukesm_30NS_tas_1.nc","air_temperature")+ncread("ukesm_30NS_tas_2.nc","air_temperature")+ncread("ukesm_30NS_tas_3.nc","air_temperature"))/3))-T_base;
T_15NS_fdbk =  globalMean(regridUKESMtoCESM((ncread("ukesm_15NS_tas_1.nc","air_temperature")+ncread("ukesm_15NS_tas_2.nc","air_temperature")+ncread("ukesm_15NS_tas_3.nc","air_temperature"))/3))-T_base;
T_60NS_fdbk =  globalMean(regridUKESMtoCESM((ncread("ukesm_60NS_tas_1.nc","air_temperature")+ncread("ukesm_60NS_tas_2.nc","air_temperature")+ncread("ukesm_60NS_tas_3.nc","air_temperature"))/3))-T_base;
T_EQ_fdbk =  globalMean(regridUKESMtoCESM((ncread("ukesm_EQ_tas_1.nc","air_temperature")+ncread("ukesm_EQ_tas_2.nc","air_temperature")+ncread("ukesm_EQ_tas_3.nc","air_temperature"))/3))-T_base;

P_30NS_fdbk = pdimchange*globalMean(regridUKESMtoCESM((ncread("ukesm_30NS_pr_1.nc","precipitation_flux")+ncread("ukesm_30NS_pr_2.nc","precipitation_flux")+ncread("ukesm_30NS_pr_3.nc","precipitation_flux"))/3))-P_base;
P_15NS_fdbk =  pdimchange*globalMean(regridUKESMtoCESM((ncread("ukesm_15NS_pr_1.nc","precipitation_flux")+ncread("ukesm_15NS_pr_2.nc","precipitation_flux")+ncread("ukesm_15NS_pr_3.nc","precipitation_flux"))/3))-P_base;
P_60NS_fdbk =  pdimchange*globalMean(regridUKESMtoCESM((ncread("ukesm_60NS_pr_1.nc","precipitation_flux")+ncread("ukesm_60NS_pr_2.nc","precipitation_flux")+ncread("ukesm_60NS_pr_3.nc","precipitation_flux"))/3))-P_base;
P_EQ_fdbk =  pdimchange*globalMean(regridUKESMtoCESM((ncread("ukesm_EQ_pr_1.nc","precipitation_flux")+ncread("ukesm_EQ_pr_2.nc","precipitation_flux")+ncread("ukesm_EQ_pr_3.nc","precipitation_flux"))/3))-P_base;

%% Load ARISE data
ukesm_ARISE_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_arise_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_arise_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_arise_strataod_3.nc","__xarray_dataarray_variable__"))/3);
ukesm_ARISE_strataod_matrix_full = cat(4,regridUKESMtoCESM((ncread("ukesm_arise_strataod_1.nc","__xarray_dataarray_variable__"))),regridUKESMtoCESM((ncread("ukesm_arise_strataod_2.nc","__xarray_dataarray_variable__"))),regridUKESMtoCESM((ncread("ukesm_arise_strataod_3.nc","__xarray_dataarray_variable__"))),regridUKESMtoCESM((ncread("ukesm_arise_strataod_4.nc","__xarray_dataarray_variable__"))),regridUKESMtoCESM((ncread("ukesm_arise_strataod_5.nc","__xarray_dataarray_variable__"))));
ukesm_ARISE_strataod_matrix_lat = squeeze(mean(ukesm_ARISE_strataod_matrix_full,1));
ukesm_ARISE_T_matrix = regridUKESMtoCESM((ncread("ukesm_arise_tas_1.nc","air_temperature")+ncread("ukesm_arise_tas_2.nc","air_temperature")+ncread("ukesm_arise_tas_3.nc","air_temperature"))/3);
ukesm_ARISE_P_matrix = pdimchange*regridUKESMtoCESM((ncread("ukesm_arise_pr_1.nc","precipitation_flux")+ncread("ukesm_arise_pr_2.nc","precipitation_flux")+ncread("ukesm_arise_pr_3.nc","precipitation_flux"))/3);

AOD_ARISE_ann_mean = averageEvery2d(12,1,globalMean(ukesm_ARISE_strataod_matrix(:,:,25:end)));
T_ARISE_ann_mean = averageEvery2d(12,1,globalMean(ukesm_ARISE_T_matrix(:,:,25:end)));
P_ARISE_ann_mean = averageEvery2d(12,1,globalMean(ukesm_ARISE_P_matrix(:,:,25:end)));

%% Initiate plot
colors = colormap(brewermap([],'Dark2'));
colors2 = colormap(brewermap([],'Set1'));
colors3 = colormap(brewermap([],'Set2'));
ann_time = (2035:2069)+.5;
% figure
tiley = tiledlayout(3,3);
set(gcf, 'Position', [200, 200, 1200,1200]) % Set figure size
tiley.Padding = 'tight';

%% Plot single point AOD
sp_legend = ["SSP2-4.5"; "60°N"; "30°N"; "15°N"; "0°N";"15°S";"30°S";"60°S"];
my_tile = nexttile(1);
hold on
plot(ann_time, averageEvery2d(12,1,AOD_SSP245),"LineWidth",2,"Color",'k')
plot(ann_time, averageEvery2d(12,1,AOD_60N+AOD_base),"LineWidth",2,"Color",colors(1,:))
plot(ann_time, averageEvery2d(12,1,AOD_30N+AOD_base),"LineWidth",2,"Color",colors(2,:))
plot(ann_time, averageEvery2d(12,1,AOD_15N+AOD_base),"LineWidth",2,"Color",colors(3,:))
plot(ann_time, averageEvery2d(12,1,AOD_0N+AOD_base),"LineWidth",2,"Color",colors(4,:))
plot(ann_time, averageEvery2d(12,1,AOD_15S+AOD_base),"LineWidth",2,"Color",colors(5,:))
plot(ann_time, averageEvery2d(12,1,AOD_30S+AOD_base),"LineWidth",2,"Color",colors(6,:))
plot(ann_time(1:14), averageEvery2d(12,1,AOD_60S+AOD_base),"LineWidth",2,"Color",colors(7,:))
for i = 1:7
    if i == 1
        injection_12_60N = [];
        for j = 1:35
            injection_12_60N = [injection_12_60N; [0 0 4 4 4 0 0 0 0 0 0 0]'];
        end
        injections = [injection_12_60N zeros(420,6)];
    elseif i == 7
        injection_12_60S = [];
        for j = 1:35
            injection_12_60S = [injection_12_60S; [0 0 0 0 0 0 0 0 4 4 4 0]'];
        end
        injections = [zeros(420,6) injection_12_60S];
    else
        injections = [zeros(420,i-1) ones(420,1) zeros(420,7-i)];
    end
    AOD_emu_matrix = CIDER_AOD_pattern_from_all_injections(param_AOD_all,injections,pattern_AOD_all);
    AOD_emu = globalMean(AOD_emu_matrix);
    plot(ann_time, averageEvery2d(12,1,AOD_emu+AOD_base),"LineWidth",2,"Color",colors(i,:),"LineStyle","--")
end
hold off
ylabel("AOD")
my_tile.TitleHorizontalAlignment = 'left';
my_title = "(a) AOD";
title(my_title,"Fontsize",20)
set(gca,'FontSize', 20)
legend(sp_legend,"Location","se","FontSize",12)
xlim([2035 2070])
box on 
grid on
%% Plot single point T
my_tile = nexttile(2);
hold on
plot(ann_time, averageEvery2d(12,1,T_SSP245+T_base),"LineWidth",2,"Color",'k')
plot(ann_time, averageEvery2d(12,1,T_60N+T_SSP245+T_base),"LineWidth",2,"Color",colors(1,:))
plot(ann_time, averageEvery2d(12,1,T_30N+T_SSP245+T_base),"LineWidth",2,"Color",colors(2,:))
plot(ann_time, averageEvery2d(12,1,T_15N+T_SSP245+T_base),"LineWidth",2,"Color",colors(3,:))
plot(ann_time, averageEvery2d(12,1,T_0N+T_SSP245+T_base),"LineWidth",2,"Color",colors(4,:))
plot(ann_time, averageEvery2d(12,1,T_15S+T_SSP245+T_base),"LineWidth",2,"Color",colors(5,:))
plot(ann_time, averageEvery2d(12,1,T_30S+T_SSP245+T_base),"LineWidth",2,"Color",colors(6,:))
plot(ann_time, averageEvery2d(12,1,T_60S+T_SSP245+T_base),"LineWidth",2,"Color",colors(7,:))
for i = 1:8
    if i == 1
        injection_12_60N = [];
        for j = 1:35
            injection_12_60N = [injection_12_60N; [0 0 4 4 4 0 0 0 0 0 0 0]'];
        end
        injections = [injection_12_60N zeros(420,6)];
        color = colors(i,:);
    elseif i == 7
        injection_12_60S = [];
        for j = 1:35
            injection_12_60S = [injection_12_60S; [0 0 0 0 0 0 0 0 4 4 4 0]'];
        end
        injections = [zeros(420,6) injection_12_60S];
        color = colors(i,:);
    elseif i ==8 
        injections = zeros(420,7);
        color = 'k';
    else
        injections = [zeros(420,i-1) ones(420,1) zeros(420,7-i)];
        color = colors(i,:);
    end
    all_injection_and_CO2 = [injections CO2_forcing_SSP245_month];
    T_emu_matrix = CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2,param_AOD_all,param_T_all,pattern_T_all);
    T_emu = globalMean(T_emu_matrix);
    plot(ann_time, averageEvery2d(12,1,T_emu+T_base),"LineWidth",2,"Color",color,"LineStyle","--")
end
hold off
% legend(sp_legend)
ylabel("Temperature (K)")
% title("Figure 1 but British (in progress)")
my_tile.TitleHorizontalAlignment = 'left';
my_title = "(b) Temperature";
title(my_title,"Fontsize",20)
set(gca,'FontSize', 20)
xlim([2035 2070])

box on 
grid on
%% Plot single point P
my_tile = nexttile(3);
hold on
plot(ann_time, averageEvery2d(12,1,P_SSP245)+P_base,"LineWidth",2,"Color",'k')
plot(ann_time, averageEvery2d(12,1,P_60N+P_SSP245)+P_base,"LineWidth",2,"Color",colors(1,:))
plot(ann_time, averageEvery2d(12,1,P_30N+P_SSP245)+P_base,"LineWidth",2,"Color",colors(2,:))
plot(ann_time, averageEvery2d(12,1,P_15N+P_SSP245)+P_base,"LineWidth",2,"Color",colors(3,:))
plot(ann_time, averageEvery2d(12,1,P_0N+P_SSP245)+P_base,"LineWidth",2,"Color",colors(4,:))
plot(ann_time, averageEvery2d(12,1,P_15S+P_SSP245)+P_base,"LineWidth",2,"Color",colors(5,:))
plot(ann_time, averageEvery2d(12,1,P_30S+P_SSP245)+P_base,"LineWidth",2,"Color",colors(6,:))
plot(ann_time, averageEvery2d(12,1,P_60S+P_SSP245)+P_base,"LineWidth",2,"Color",colors(7,:))
for i = 1:8
    if i == 1
        injection_12_60N = [];
        for j = 1:35
            injection_12_60N = [injection_12_60N; [0 0 4 4 4 0 0 0 0 0 0 0]'];
        end
        injections = [injection_12_60N zeros(420,6)];
        color = colors(i,:);
    elseif i == 7
        injection_12_60S = [];
        for j = 1:35
            injection_12_60S = [injection_12_60S; [0 0 0 0 0 0 0 0 4 4 4 0]'];
        end
        injections = [zeros(420,6) injection_12_60S];
        color = colors(i,:);
    elseif i ==8 
        injections = zeros(420,7);
        color = 'k';
    else
        injections = [zeros(420,i-1) ones(420,1) zeros(420,7-i)];
        color = colors(i,:);
    end
    all_injection_and_CO2 = [injections CO2_forcing_SSP245_month];
    P_emu_matrix = CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2,param_AOD_all,param_P_all,pattern_P_all);
    P_emu = globalMean(P_emu_matrix);
    plot(ann_time, averageEvery2d(12,1,P_emu)+P_base,"LineWidth",2,"Color",color,"LineStyle","--")
end
hold off
% legend(sp_legend)
ylabel("Precipitation (mm day^{-1})")
my_tile.TitleHorizontalAlignment = 'left';
my_title = "(c) Precipitation";
title(my_title,"Fontsize",20)
set(gca,'FontSize', 20)
xlim([2035 2070])

box on 
grid on
%% Plot feedback AOD
my_tile = nexttile(4);
hold on
plot(ann_time, averageEvery2d(12,1,AOD_SSP245),"LineWidth",2,"Color",'k')
plot(ann_time, averageEvery2d(12,1,AOD_EQ_fdbk+AOD_base),"LineWidth",2,"Color",colors2(1,:))
plot(ann_time, averageEvery2d(12,1,AOD_15NS_fdbk+AOD_base),"LineWidth",2,"Color",colors2(2,:))
plot(ann_time, averageEvery2d(12,1,AOD_30NS_fdbk+AOD_base),"LineWidth",2,"Color",colors2(3,:))
plot(ann_time, averageEvery2d(12,1,AOD_60NS_fdbk+AOD_base),"LineWidth",2,"Color",colors2(4,:))

AOD_fdbk_EQ = globalMean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,injection_matrix_fdbk_0N,pattern_AOD_all));
AOD_fdbk_15NS = globalMean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,injection_matrix_fdbk_15NS,pattern_AOD_all));
AOD_fdbk_30NS = globalMean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,injection_matrix_fdbk_30NS,pattern_AOD_all));
AOD_fdbk_60NS = globalMean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,injection_matrix_fdbk_60NS,pattern_AOD_all));

plot(ann_time, averageEvery2d(12,1,AOD_fdbk_EQ+AOD_base),"LineWidth",2,"Color",colors2(1,:),"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,AOD_fdbk_15NS+AOD_base),"LineWidth",2,"Color",colors2(2,:),"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,AOD_fdbk_30NS+AOD_base),"LineWidth",2,"Color",colors2(3,:),"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,AOD_fdbk_60NS+AOD_base),"LineWidth",2,"Color",colors2(4,:),"LineStyle","--")

hold off
legend("SSP2-4.5","0°N Feedback for 1.0°C","15°N+15°S Feedback for 1.0°C","30°N+30°S Feedback for 1.0°C","60°N+60°S Feedback for 1.0°C","Location",'se',"Fontsize",12)
% xlabel("Year")
ylabel("AOD")
box on 
grid on
% set(gca,'Linewidth',2)
set(gca,'FontSize', 20)
% xlim([2035 2070])
% my_tile.TitleHorizontalAlignment = 'left';
% my_title = "(d)";
% title(my_title,"Fontsize",20)
my_tile.TitleHorizontalAlignment = 'left';
my_title = "(d)";
title(my_title,"Fontsize",20)
box on 
grid on
xlim([2035 2070])

%% Plot Feedback T
my_tile = nexttile(5);
hold on
plot(ann_time, averageEvery2d(12,1,T_SSP245+T_base),"LineWidth",2,"Color",'k')
plot(ann_time, averageEvery2d(12,1,T_EQ_fdbk(25:end)+T_base),"LineWidth",2,"Color",colors2(1,:))
plot(ann_time, averageEvery2d(12,1,T_15NS_fdbk(25:end)+T_base),"LineWidth",2,"Color",colors2(2,:))
plot(ann_time, averageEvery2d(12,1,T_30NS_fdbk(25:end)+T_base),"LineWidth",2,"Color",colors2(3,:))
plot(ann_time, averageEvery2d(12,1,T_60NS_fdbk(25:end)+T_base),"LineWidth",2,"Color",colors2(4,:))

T_emu_SSP245 = globalMean(CIDER_pattern_from_all_injections_and_CO2([0*injection_matrix_fdbk_0N CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
T_emu_fdbk_EQ = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_0N CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
T_emu_fdbk_15NS = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_15NS CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
T_emu_fdbk_30NS = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_30NS CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
T_emu_fdbk_60NS = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_60NS CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));

plot(ann_time, averageEvery2d(12,1,T_emu_SSP245+T_base),"LineWidth",2,"Color",'k',"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,T_emu_fdbk_EQ+T_base),"LineWidth",2,"Color",colors2(1,:),"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,T_emu_fdbk_15NS+T_base),"LineWidth",2,"Color",colors2(2,:),"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,T_emu_fdbk_30NS+T_base),"LineWidth",2,"Color",colors2(3,:),"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,T_emu_fdbk_60NS+T_base),"LineWidth",2,"Color",colors2(4,:),"LineStyle","--")

hold off
% legend("SSP2-4.5","0°N Feedback for 1.0°C","15°N+15°S Feedback for 1.0°C","30°N+30°S Feedback for 1.0°C","60°N+60°S Feedback for 1.0°C","Location",'se')
% xlabel("Year")
set(gca,'FontSize', 20)
xlim([2035 2070])

ylabel("Temperature (K)")
box on 
grid on
my_tile.TitleHorizontalAlignment = 'left';
my_title = "(e)";
title(my_title,"Fontsize",20)
%% Plot Feedback P
my_tile = nexttile(6);
hold on
plot(ann_time, averageEvery2d(12,1,P_SSP245)+P_base,"LineWidth",2,"Color",'k')
plot(ann_time, averageEvery2d(12,1,P_EQ_fdbk(25:end))+P_base,"LineWidth",2,"Color",colors2(1,:))
plot(ann_time, averageEvery2d(12,1,P_15NS_fdbk(25:end))+P_base,"LineWidth",2,"Color",colors2(2,:))
plot(ann_time, averageEvery2d(12,1,P_30NS_fdbk(25:end))+P_base,"LineWidth",2,"Color",colors2(3,:))
plot(ann_time, averageEvery2d(12,1,P_60NS_fdbk(25:end))+P_base,"LineWidth",2,"Color",colors2(4,:))

P_emu_SSP245 = globalMean(CIDER_pattern_from_all_injections_and_CO2([0*injection_matrix_fdbk_0N CO2_forcing_SSP245_month], param_AOD_all,param_P_all,pattern_P_all));
P_emu_fdbk_EQ = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_0N CO2_forcing_SSP245_month], param_AOD_all,param_P_all,pattern_P_all));
P_emu_fdbk_15NS = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_15NS CO2_forcing_SSP245_month], param_AOD_all,param_P_all,pattern_P_all));
P_emu_fdbk_30NS = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_30NS CO2_forcing_SSP245_month], param_AOD_all,param_P_all,pattern_P_all));
P_emu_fdbk_60NS = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_60NS CO2_forcing_SSP245_month], param_AOD_all,param_P_all,pattern_P_all));

plot(ann_time, averageEvery2d(12,1,P_emu_SSP245+P_base),"LineWidth",2,"Color",'k',"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,P_emu_fdbk_EQ+P_base),"LineWidth",2,"Color",colors2(1,:),"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,P_emu_fdbk_15NS+P_base),"LineWidth",2,"Color",colors2(2,:),"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,P_emu_fdbk_30NS+P_base),"LineWidth",2,"Color",colors2(3,:),"LineStyle","--")
plot(ann_time, averageEvery2d(12,1,P_emu_fdbk_60NS+P_base),"LineWidth",2,"Color",colors2(4,:),"LineStyle","--")

hold off
% legend("SSP2-4.5","0°N Feedback for 1.0°C","15°N+15°S Feedback for 1.0°C","30°N+30°S Feedback for 1.0°C","60°N+60°S Feedback for 1.0°C","Location",'se')
% xlabel("Year")
ylabel("Precipitation (mm day^{-1})")
box on 
grid on
my_tile.TitleHorizontalAlignment = 'left';
my_title = "(f)";
title(my_title,"Fontsize",20)
set(gca,'FontSize', 20)
xlim([2035 2070])

%% Plot ARISE AOD
my_tile = nexttile(7);
hold on
plot(ann_time, averageEvery2d(12,1,AOD_SSP245),"LineWidth",2,"Color",'k')
plot(ann_time, AOD_ARISE_ann_mean,"LineWidth",2,"Color",colors3(1,:))
AOD_fdbk_arise = globalMean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,injection_matrix_ARISE,pattern_AOD_all));
plot(ann_time, averageEvery2d(12,1,AOD_fdbk_arise+AOD_base),"LineWidth",2,"Color",colors3(1,:),"LineStyle","--")
ylabel("AOD")
hold off
my_tile.TitleHorizontalAlignment = 'left';
my_title = "(g)";
title(my_title,"Fontsize",20)
set(gca,'FontSize', 20)

legend("SSP2-4.5","Multi-Objective SAI","Fontsize",12,"Location","se")
box on 
grid on
xlim([2035 2070])
CIDER_UKESM_ARISE_AOD_lat = squeeze(mean(CIDER_AOD_pattern_from_all_injections(param_AOD_all,injection_matrix_ARISE,pattern_AOD_all),1));
model_UKESM_ARISE_AOD_lat = squeeze(mean(ukesm_ARISE_strataod_matrix(:,:,25:end),1));
save("Uncoord_Emu_Paper/UKESM_ARISE_comparison.mat","model_UKESM_ARISE_AOD_lat","CIDER_UKESM_ARISE_AOD_lat")
%% Plot ARISE T
my_tile = nexttile(8);
hold on
plot(ann_time, averageEvery2d(12,1,T_SSP245)+T_base,"LineWidth",2,"Color",'k')
plot(ann_time, T_ARISE_ann_mean,"LineWidth",2,"Color",colors3(1,:))
plot(ann_time, averageEvery2d(12,1,T_emu_SSP245+T_base),"LineWidth",2,"Color",'k',"LineStyle","--")
T_emu_arise = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_ARISE CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
plot(ann_time, averageEvery2d(12,1,T_emu_arise+T_base),"LineWidth",2,"Color",colors3(1,:),"LineStyle","--")
ylabel("Temperature (K)")
hold off
my_tile.TitleHorizontalAlignment = 'left';
my_title = "(h)";
title(my_title,"Fontsize",20)
set(gca,'FontSize', 20)
xlim([2035 2070])

box on 
grid on
%% Plot ARISE P
my_tile = nexttile(9);
hold on
plot(ann_time, averageEvery2d(12,1,P_SSP245)+P_base,"LineWidth",2,"Color",'k')
plot(ann_time, P_ARISE_ann_mean,"LineWidth",2,"Color",colors3(1,:))
plot(ann_time, averageEvery2d(12,1,P_emu_SSP245+P_base),"LineWidth",2,"Color",'k',"LineStyle","--")
P_emu_arise = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_ARISE CO2_forcing_SSP245_month], param_AOD_all,param_P_all,pattern_P_all));
plot(ann_time, averageEvery2d(12,1,P_emu_arise+P_base),"LineWidth",2,"Color",colors3(1,:),"LineStyle","--")
ylabel("Precipitation (mm day^{-1})")
hold off
my_tile.TitleHorizontalAlignment = 'left';
my_title = "(i)";
title(my_title,"Fontsize",20)
box on 
set(gca,'FontSize', 20)
xlim([2035 2070])

grid on
%% Save Figure
set(gcf,'renderer','painters')
print(gcf,'-dpng',["Uncoord_Emu_Paper/Uncoord_Plots/Figure_4_" + getNow() + ".png"],'-r300')

