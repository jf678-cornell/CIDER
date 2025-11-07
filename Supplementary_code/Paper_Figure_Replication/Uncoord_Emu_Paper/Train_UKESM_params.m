clc
addAllPaths
AOD_weightings = [1 1];

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

feedback_inj_0N_monthly = repeatElements(feedback_inj_0N,12)/12;
feedback_inj_15_monthly = repeatElements(feedback_inj_15,12)/12/2;
feedback_inj_30_monthly = repeatElements(feedback_inj_30,12)/12/2;
feedback_inj_60_monthly_base = repeatElements(feedback_inj_60,12)/3/2;
feedback_inj_60N_monthly = injection_12_60N/4.*feedback_inj_60_monthly_base;
feedback_inj_60S_monthly = injection_12_60S/4.*feedback_inj_60_monthly_base;
feedback_inj_60NS_monthly = feedback_inj_60N_monthly+feedback_inj_60S_monthly;
feedback_inj_60N_monthly = feedback_inj_60NS_monthly/2;
feedback_inj_60S_monthly = feedback_inj_60NS_monthly/2;


all_step_injections = [injection_12_60N ones(420,5)];
all_feedback_injections = [feedback_inj_60N_monthly feedback_inj_30_monthly feedback_inj_15_monthly feedback_inj_0N_monthly feedback_inj_15_monthly feedback_inj_30_monthly];
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
pattern_T_base = T_SSP245_matrix_base;

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
pattern_P_base = P_SSP245_matrix_base;

%% Parse SP data
AOD_0N = globalMean(AOD_0N_matrix);
AOD_15N = globalMean(AOD_15N_matrix);
AOD_15S = globalMean(AOD_15S_matrix);
AOD_30N = globalMean(AOD_30N_matrix);
AOD_30S = globalMean(AOD_30S_matrix);
AOD_60N = globalMean(AOD_60N_matrix);
AOD_60S = globalMean(AOD_60S_matrix);
all_AOD_step_responses = [AOD_60N AOD_30N AOD_15N AOD_0N AOD_15S AOD_30S];

ssl = 240;
pattern_AOD_all = cat(3,CIDER_get_pattern(AOD_60N_matrix,ssl),CIDER_get_pattern(AOD_30N_matrix,ssl),CIDER_get_pattern(AOD_15N_matrix,ssl),CIDER_get_pattern(AOD_0N_matrix,ssl),CIDER_get_pattern(AOD_15S_matrix,ssl),CIDER_get_pattern(AOD_30S_matrix,ssl),CIDER_get_pattern(AOD_60S_matrix,ssl/2));


T_0N = globalMean(T_0N_matrix);
T_15N = globalMean(T_15N_matrix);
T_15S = globalMean(T_15S_matrix);
T_30N = globalMean(T_30N_matrix);
T_30S = globalMean(T_30S_matrix);
T_60N = globalMean(T_60N_matrix);
T_60S = globalMean(T_60S_matrix);
all_T_step_responses = [T_60N T_30N T_15N T_0N T_15S T_30S];

T_SSP245 = globalMean(T_SSP245_matrix);
ssl = 240;
pattern_T_all = cat(3,CIDER_get_pattern(T_60N_matrix,ssl),CIDER_get_pattern(T_30N_matrix,ssl),CIDER_get_pattern(T_15N_matrix,ssl),CIDER_get_pattern(T_0N_matrix,ssl),CIDER_get_pattern(T_15S_matrix,ssl),CIDER_get_pattern(T_30S_matrix,ssl),CIDER_get_pattern(T_60S_matrix,ssl),CIDER_get_pattern(T_SSP245_matrix,ssl));

P_0N = globalMean(P_0N_matrix);
P_15N = globalMean(P_15N_matrix);
P_15S = globalMean(P_15S_matrix);
P_30N = globalMean(P_30N_matrix);
P_30S = globalMean(P_30S_matrix);
P_60N = globalMean(P_60N_matrix);
P_60S = globalMean(P_60S_matrix);
all_P_step_responses = [P_60N P_30N P_15N P_0N P_15S P_30S];

P_SSP245 = globalMean(P_SSP245_matrix);
ssl = 240;
pattern_P_all = cat(3,CIDER_get_pattern(P_60N_matrix,ssl),CIDER_get_pattern(P_30N_matrix,ssl),CIDER_get_pattern(P_15N_matrix,ssl),CIDER_get_pattern(P_0N_matrix,ssl),CIDER_get_pattern(P_15S_matrix,ssl),CIDER_get_pattern(P_30S_matrix,ssl),CIDER_get_pattern(P_60S_matrix,ssl),CIDER_get_pattern(P_SSP245_matrix,ssl));

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

AOD_60NS_fdbk = globalMean(ukesm_60NS_strataod_matrix-ukesm_SSP245_strataod_matrix);
AOD_30NS_fdbk = globalMean(ukesm_30NS_strataod_matrix-ukesm_SSP245_strataod_matrix);
AOD_15NS_fdbk = globalMean(ukesm_15NS_strataod_matrix-ukesm_SSP245_strataod_matrix);
AOD_EQ_fdbk = globalMean(ukesm_EQ_strataod_matrix-ukesm_SSP245_strataod_matrix);

ratio_60N = mean(AOD_60N(181:end))/(mean(AOD_60N(181:end))+mean(AOD_60S(48:end)));
ratio_30N = mean(AOD_30N(181:end))/(mean(AOD_30N(181:end))+mean(AOD_30S(181:end)));
ratio_15N = mean(AOD_15N(181:end))/(mean(AOD_15N(181:end))+mean(AOD_15S(181:end)));

all_AOD_fdbk_responses = 1/2*[AOD_60NS_fdbk AOD_30NS_fdbk AOD_15NS_fdbk 2*AOD_EQ_fdbk AOD_15NS_fdbk AOD_30NS_fdbk];
all_AOD_fdbk_responses = [ratio_60N*AOD_60NS_fdbk ratio_30N*AOD_30NS_fdbk ratio_15N*AOD_15NS_fdbk AOD_EQ_fdbk (1-ratio_15N)*AOD_15NS_fdbk (1-ratio_30N)*AOD_30NS_fdbk];

%% Train AOD
weightings = AOD_weightings;
all_params_AOD_except_60S = CIDER_train_AOD_params(all_step_injections,all_AOD_step_responses,all_feedback_injections,all_AOD_fdbk_responses,weightings);
params_AOD_60S = CIDER_train_AOD_params(injection_12_60S(1:172),AOD_60S,feedback_inj_60S_monthly,(1-ratio_60N)*AOD_60NS_fdbk,weightings);
param_AOD_all = [all_params_AOD_except_60S;params_AOD_60S];
% plot(CIDER_AOD_from_injection(param_AOD_all(4,:),injection_12))

%% Train T
param_T_all_inj_except_60S = CIDER_train_climate_params(all_AOD_step_responses,all_T_step_responses);
AOD_60S_emu = CIDER_AOD_from_injection(param_AOD_all(7,:),injection_12_60S);
param_T_60S = CIDER_train_climate_params(AOD_60S ,T_60S(1:172));
% param_T_60S_emu = CIDER_train_climate_params(AOD_60S_emu ,T_60S);
% plot(CIDER_response_from_1_forcing(param_T_60S,AOD_60S_emu))
% hold on
% plot(T_60S)
param_T_CO2 = CIDER_train_climate_params( CO2_forcing_SSP245_month,T_SSP245);
param_T_all = [param_T_all_inj_except_60S;param_T_60S;param_T_CO2];

%% Train P
lb = [0;-1];
ub = [5000;1];

param_P_all_inj_except_60S = CIDER_train_climate_params(all_AOD_step_responses,all_P_step_responses,param_T_all_inj_except_60S,lb,ub);
param_P_60S = CIDER_train_climate_params(AOD_60S ,P_60S(1:172),param_T_60S,lb,ub);
param_P_CO2 = CIDER_train_climate_params(CO2_forcing_SSP245_month,P_SSP245,param_T_CO2,lb,ub);
param_P_all = [param_P_all_inj_except_60S;param_P_60S;param_P_CO2];

save("Uncoord_Emu_Paper/UKESM_params.mat","param_AOD_all","param_P_all","param_T_all","pattern_P_all","pattern_T_all","pattern_AOD_all","pattern_T_base","pattern_P_base");
clc
disp("Done training!")
%%
% all_injection_and_CO2 = [zeros(420,6) ones(420,1) zeros(420,0) CO2_forcing_SSP245_month];
% sample_response = CIDER_response_from_all_injections_and_CO2(all_injection_and_CO2,param_AOD_all,param_T_all);
% figure
% hold on
% plot(1:length(sample_response),sample_response)
% plot(1:length(sample_response),T_60S+T_SSP245)
%%
figure