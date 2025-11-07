addAllPaths
close all
load CESM_params.mat
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

% feedback_inj_0N = readtable("ukesm_EQ_inj_1.log").x30S_Tg_(2:end);
% feedback_inj_15 = 2*readtable("ukesm_15NS_inj_1.log").x30S_Tg_(2:end-1);
% feedback_inj_30 = 2*readtable("ukesm_30NS_inj_1.log").x30S_Tg_(2:end);
% feedback_inj_60 = 2*readtable("ukesm_60NS_inj_1.log").x30S_Tg_(2:end);

feedback_inj_0N = inj_rate_arr(:,1);
feedback_inj_15 = inj_rate_arr(:,2);
feedback_inj_30 = inj_rate_arr(:,3);
feedback_inj_60 = inj_rate_arr(:,4);

feedback_inj_0N_monthly = repeatElements(feedback_inj_0N,12)/12;
feedback_inj_15_monthly = repeatElements(feedback_inj_15,12)/12/2;
feedback_inj_30_monthly = repeatElements(feedback_inj_30,12)/12/2;
feedback_inj_60_monthly_base = repeatElements(feedback_inj_60,12)/3/2;
feedback_inj_60N_monthly = injection_12_60N/4.*feedback_inj_60_monthly_base;
feedback_inj_60S_monthly = injection_12_60S/4.*feedback_inj_60_monthly_base;
feedback_inj_60NS_monthly = feedback_inj_60N_monthly+feedback_inj_60S_monthly;
feedback_inj_60N_monthly = feedback_inj_60NS_monthly/2;
feedback_inj_60S_monthly = feedback_inj_60NS_monthly/2;

injection_matrix_fdbk_0N = [zeros(420,3) feedback_inj_0N_monthly zeros(420,3)];
injection_matrix_fdbk_15NS = [zeros(420,2) feedback_inj_15_monthly zeros(420,1) feedback_inj_15_monthly zeros(420,2)];
injection_matrix_fdbk_30NS = [zeros(420,1) feedback_inj_30_monthly zeros(420,3) feedback_inj_30_monthly zeros(420,1)];
injection_matrix_fdbk_60NS = [zeros(420,0) feedback_inj_60N_monthly zeros(420,5) feedback_inj_60S_monthly zeros(420,0)];


all_step_injections = [injection_12_60N ones(420,5)];
all_feedback_injections = [feedback_inj_60N_monthly feedback_inj_30_monthly feedback_inj_15_monthly feedback_inj_0N_monthly feedback_inj_15_monthly feedback_inj_30_monthly];
%%
L_DEFAULT = readLog('GAUSS-DEFAULT',[1 2 3]);
injection_default = [mean(L_DEFAULT.S30,2) mean(L_DEFAULT.S15,2) mean(L_DEFAULT.N15,2) mean(L_DEFAULT.N30,2)];
injection_default_monthly = zeros(length(injection_default(1:end-1,1))*12,7);
injection_default_monthly(:,2) = repeatElements(mean(L_DEFAULT.N30(1:end-1,:),2),12)/12;
injection_default_monthly(:,3) = repeatElements(mean(L_DEFAULT.N15(1:end-1,:),2),12)/12;
injection_default_monthly(:,5) = repeatElements(mean(L_DEFAULT.S15(1:end-1,:),2),12)/12;
injection_default_monthly(:,6) = repeatElements(mean(L_DEFAULT.S30(1:end-1,:),2),12)/12;

L_LOWER_5 = readLog('GAUSS-LOWER-0.5',[1 2 3]);

injection_lower_5 = [mean(L_LOWER_5.S30,2) mean(L_LOWER_5.S15,2) mean(L_LOWER_5.N15,2) mean(L_LOWER_5.N30,2)];
injection_lower_5_monthly = zeros(length(injection_lower_5(1:end-1,1))*12,7);
injection_lower_5_monthly(:,2) = repeatElements(mean(L_LOWER_5.N30(1:end-1,:),2),12)/12;
injection_lower_5_monthly(:,3) = repeatElements(mean(L_LOWER_5.N15(1:end-1,:),2),12)/12;
injection_lower_5_monthly(:,5) = repeatElements(mean(L_LOWER_5.S15(1:end-1,:),2),12)/12;
injection_lower_5_monthly(:,6) = repeatElements(mean(L_LOWER_5.S30(1:end-1,:),2),12)/12;

L_LOWER_1 = readLog('GAUSS-LOWER-1.0',[1 2 3]);

injection_lower_1 = [mean(L_LOWER_1.S30,2) mean(L_LOWER_1.S15,2) mean(L_LOWER_1.N15,2) mean(L_LOWER_1.N30,2)];
injection_lower_1_monthly = zeros(length(injection_lower_1(1:end-1,1))*12,7);
injection_lower_1_monthly(:,2) = repeatElements(mean(L_LOWER_1.N30(1:end-1,:),2),12)/12;
injection_lower_1_monthly(:,3) = repeatElements(mean(L_LOWER_1.N15(1:end-1,:),2),12)/12;
injection_lower_1_monthly(:,5) = repeatElements(mean(L_LOWER_1.S15(1:end-1,:),2),12)/12;
injection_lower_1_monthly(:,6) = repeatElements(mean(L_LOWER_1.S30(1:end-1,:),2),12)/12;

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
% AOD_0N_matrix  = ncread("u-ch622_AOD.nc","AOD550_delta");
% AOD_15N_matrix = ncread("u-cd352_AOD.nc","AOD550_delta");
% AOD_15S_matrix = ncread("u-cd353_AOD.nc","AOD550_delta");
% AOD_30N_matrix = ncread("u-cd297_AOD.nc","AOD550_delta");
% AOD_30S_matrix = ncread("u-cd354_AOD.nc","AOD550_delta");
% AOD_60N_matrix = ncread("u-cm174_AOD.nc","AOD550_delta");
% AOD_60S_matrix = ncread("u-cm175_AOD.nc","AOD550_delta");

% ssp245_AOD_data = processRun('AODVISstdn',1,'SSP245','203501-206912.nc',[2035 2069],[1 2 3]);
% ssp245_AOD_base = mean(ssp245_AOD_data.ensemble_annual_average);
% 
% AOD_60N_matrix = readSinglePoint(60,'N',1,'AODVISstdn');
% AOD_30N_matrix = mean(readSinglePoint(30,'N',[1 2 3],'AODVISstdn'),4);
% AOD_15N_matrix = mean(readSinglePoint(15,'N',[1 2 3],'AODVISstdn'),4);
% AOD_0N_matrix  = mean(readSinglePoint( 0,'N',[1 2 3],'AODVISstdn'),4);
% AOD_15S_matrix = mean(readSinglePoint(15,'S',[1 2 3],'AODVISstdn'),4);
% AOD_30S_matrix = mean(readSinglePoint(30,'S',[1 2 3],'AODVISstdn'),4);
% AOD_60S_matrix = readSinglePoint(60,'S',1,'AODVISstdn');


% T_0N_matrix = ncread("u-ch622_tas.nc","tas_delta");
% T_15N_matrix = ncread("u-cd352_tas.nc","tas_delta");
% T_15S_matrix = ncread("u-cd353_tas.nc","tas_delta");
% T_30N_matrix = ncread("u-cd297_tas.nc","tas_delta");
% T_30S_matrix = ncread("u-cd354_tas.nc","tas_delta");
% T_60N_matrix = ncread("u-cm174_tas.nc","tas_delta");
% T_60S_matrix = ncread("u-cm175_tas.nc","tas_delta");



ssp245_T_data = processRun('TREFHT',1,'SSP245','203501-206912.nc',[2015 2069],[1 2 3]);
% T_SSP245_matrix = regridUKESMtoCESM((ncread("tas_ssp_ukesm_1.nc","tas")+ncread("tas_ssp_ukesm_2.nc","tas")+ncread("tas_ssp_ukesm_3.nc","tas"))/3);
T_SSP245_matrix = ssp245_T_data.ensemble_monthly_matrix;
T_SSP245_matrix_base = mean(T_SSP245_matrix(:,:,12*10+1:12*30),3);
T_SSP245_matrix = T_SSP245_matrix(:,:,12*20+1:12*55) - T_SSP245_matrix_base;
T_base = globalMean(T_SSP245_matrix_base);
pattern_T_base = T_SSP245_matrix_base;
T_60N_matrix = readSinglePoint(60,'N',1,'TREFHT')-T_SSP245_matrix- T_SSP245_matrix_base;
T_30N_matrix = mean(readSinglePoint(30,'N',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_15N_matrix = mean(readSinglePoint(15,'N',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_0N_matrix  = mean(readSinglePoint( 0,'N',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_15S_matrix = mean(readSinglePoint(15,'S',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_30S_matrix = mean(readSinglePoint(30,'S',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_60S_matrix = readSinglePoint(60,'S',1,'TREFHT')-T_SSP245_matrix- T_SSP245_matrix_base;
pdimchange = 1000*24*60*60;
% P_60N_matrix = pdimchange*ncread("u-cm174_pr.nc","pr_delta");
% P_30N_matrix = pdimchange*ncread("u-cd297_pr.nc","pr_delta");
% P_15N_matrix = pdimchange*ncread("u-cd352_pr.nc","pr_delta");
% P_0N_matrix  = pdimchange*ncread("u-ch622_pr.nc","pr_delta");
% P_15S_matrix = pdimchange*ncread("u-cd353_pr.nc","pr_delta");
% P_30S_matrix = pdimchange*ncread("u-cd354_pr.nc","pr_delta");
% P_60S_matrix = pdimchange*ncread("u-cm175_pr.nc","pr_delta");

% ssp245_P_data = processRun('PRECT',pdimchange,'SSP245','203501-206912.nc',[2015 2069],[1 2 3]);
% % P_SSP245_matrix = pdimchange*regridUKESMtoCESM((ncread("pr_ssp_ukesm_1.nc","pr")+ncread("pr_ssp_ukesm_2.nc","pr")+ncread("pr_ssp_ukesm_3.nc","pr")+ncread("pr_ssp_ukesm_4.nc","pr")+ncread("pr_ssp_ukesm_5.nc","pr"))/5);
% P_SSP245_matrix = ssp245_P_data.ensemble_monthly_matrix;
% P_SSP245_matrix_base = mean(P_SSP245_matrix(:,:,12*10+1:12*30),3);
% P_SSP245_matrix = P_SSP245_matrix(:,:,12*20+1:12*55) - P_SSP245_matrix_base;
% pattern_P_base = P_SSP245_matrix_base;
% P_base = globalMean(P_SSP245_matrix_base);
% P_60N_matrix = pdimchange*readSinglePoint(60,'N',1,'PRECT')-P_SSP245_matrix- P_SSP245_matrix_base;
% P_30N_matrix = mean(pdimchange*readSinglePoint(30,'N',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
% P_15N_matrix = mean(pdimchange*readSinglePoint(15,'N',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
% P_0N_matrix  = mean(pdimchange*readSinglePoint( 0,'N',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
% P_15S_matrix = mean(pdimchange*readSinglePoint(15,'S',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
% P_30S_matrix = mean(pdimchange*readSinglePoint(30,'S',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
% P_60S_matrix = pdimchange*readSinglePoint(60,'S',1,'PRECT')-P_SSP245_matrix- P_SSP245_matrix_base;

%% Parse SP data
% AOD_0N = globalMean(AOD_0N_matrix)-ssp245_AOD_base;
% AOD_15N = globalMean(AOD_15N_matrix)-ssp245_AOD_base;
% AOD_15S = globalMean(AOD_15S_matrix)-ssp245_AOD_base;
% AOD_30N = globalMean(AOD_30N_matrix)-ssp245_AOD_base;
% AOD_30S = globalMean(AOD_30S_matrix)-ssp245_AOD_base;
% AOD_60N = globalMean(AOD_60N_matrix)-ssp245_AOD_base;
% AOD_60S = globalMean(AOD_60S_matrix)-ssp245_AOD_base;
% 
% 
% all_AOD_step_responses = [AOD_60N AOD_30N AOD_15N AOD_0N AOD_15S AOD_30S];
% 
% ssl = 240;
% pattern_AOD_all = cat(3,CIDER_get_pattern(AOD_60N_matrix,ssl),CIDER_get_pattern(AOD_30N_matrix,ssl),CIDER_get_pattern(AOD_15N_matrix,ssl),CIDER_get_pattern(AOD_0N_matrix,ssl),CIDER_get_pattern(AOD_15S_matrix,ssl),CIDER_get_pattern(AOD_30S_matrix,ssl),CIDER_get_pattern(AOD_60S_matrix,ssl/2));


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

% P_0N = globalMean(P_0N_matrix);
% P_15N = globalMean(P_15N_matrix);
% P_15S = globalMean(P_15S_matrix);
% P_30N = globalMean(P_30N_matrix);
% P_30S = globalMean(P_30S_matrix);
% P_60N = globalMean(P_60N_matrix);
% P_60S = globalMean(P_60S_matrix);
% all_P_step_responses = [P_60N P_30N P_15N P_0N P_15S P_30S];
% 
% P_SSP245 = globalMean(P_SSP245_matrix);
% ssl = 240;
% pattern_P_all = cat(3,CIDER_get_pattern(P_60N_matrix,ssl),CIDER_get_pattern(P_30N_matrix,ssl),CIDER_get_pattern(P_15N_matrix,ssl),CIDER_get_pattern(P_0N_matrix,ssl),CIDER_get_pattern(P_15S_matrix,ssl),CIDER_get_pattern(P_30S_matrix,ssl),CIDER_get_pattern(P_60S_matrix,ssl),CIDER_get_pattern(P_SSP245_matrix,ssl));

%% Load feedback data
% ukesm_SSP245_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_SSP245_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_SSP245_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_SSP245_strataod_3.nc","__xarray_dataarray_variable__"))/3);
% ukesm_SSP245_strataod_matrix = ukesm_SSP245_strataod_matrix(:,:,12*20+1:12*55);
% ukesm_60NS_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_60NS_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_60NS_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_60NS_strataod_3.nc","__xarray_dataarray_variable__"))/3);
% ukesm_60NS_strataod_matrix = ukesm_60NS_strataod_matrix(:,:,25:end);
% ukesm_30NS_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_30NS_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_30NS_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_30NS_strataod_3.nc","__xarray_dataarray_variable__"))/3);
% ukesm_30NS_strataod_matrix = ukesm_30NS_strataod_matrix(:,:,25:end);
% ukesm_15NS_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_15NS_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_15NS_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_15NS_strataod_3.nc","__xarray_dataarray_variable__"))/3);
% ukesm_15NS_strataod_matrix = ukesm_15NS_strataod_matrix(:,:,25:end);
% ukesm_EQ_strataod_matrix = regridUKESMtoCESM((ncread("ukesm_EQ_strataod_1.nc","__xarray_dataarray_variable__")+ncread("ukesm_EQ_strataod_2.nc","__xarray_dataarray_variable__")+ncread("ukesm_EQ_strataod_3.nc","__xarray_dataarray_variable__"))/3);
% ukesm_EQ_strataod_matrix = ukesm_EQ_strataod_matrix(:,:,25:end);
% 
% 
% AOD_60NS_fdbk = globalMean(ukesm_60NS_strataod_matrix-ukesm_SSP245_strataod_matrix);
% AOD_30NS_fdbk = globalMean(ukesm_30NS_strataod_matrix-ukesm_SSP245_strataod_matrix);
% AOD_15NS_fdbk = globalMean(ukesm_15NS_strataod_matrix-ukesm_SSP245_strataod_matrix);
% AOD_EQ_fdbk = globalMean(ukesm_EQ_strataod_matrix-ukesm_SSP245_strataod_matrix);
% 
% ssp245_AOD_data = processRun('AODVISstdn',1,'SSP245','203501-206912.nc',[2035 2098],[1 2 3]);
% default_AOD_data = processRun('AODVISstdn',1,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3]);
% lower_5_AOD_data = processRun('AODVISstdn',1,'GAUSS-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% % feedback_60_AOD_data = processRun('AODVISstdn',1,'GAUSS-60N_60S-LOWER-0.5','203501-206912.nc',[2035 2069],[1]);
% feedback_30_AOD_data = processRun('AODVISstdn',1,'GAUSS-30N_30S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% feedback_15_AOD_data = processRun('AODVISstdn',1,'GAUSS-15N_15S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% feedback_0N_AOD_data = processRun('AODVISstdn',1,'GAUSS-0N-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% AOD_base = mean(ssp245_AOD_data.ensemble_monthly_average);
% ssp245_AOD_data_above_base = ssp245_AOD_data.ensemble_monthly_average - AOD_base;
% default_AOD_data_above_base = default_AOD_data.ensemble_monthly_average - AOD_base;
% lower_5_AOD_data_above_base = lower_5_AOD_data.ensemble_monthly_average - AOD_base;
% AOD_EQ_fdbk = feedback_0N_AOD_data.ensemble_monthly_average - AOD_base;
% AOD_15NS_fdbk = feedback_15_AOD_data.ensemble_monthly_average - AOD_base;
% AOD_30NS_fdbk = feedback_30_AOD_data.ensemble_monthly_average - AOD_base;
% 
% all_AOD_fdbk_responses = 1/2*[AOD_60NS_fdbk AOD_30NS_fdbk AOD_15NS_fdbk 2*AOD_EQ_fdbk AOD_15NS_fdbk AOD_30NS_fdbk];

feedback_30_T_data = processRun('TREFHT',1,'GAUSS-30N_30S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
feedback_15_T_data = processRun('TREFHT',1,'GAUSS-15N_15S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
feedback_0N_T_data = processRun('TREFHT',1,'GAUSS-0N-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
default_T_data = processRun('TREFHT',1,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3]);
lower_5_T_data = processRun('TREFHT',1,'GAUSS-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
lower_1_T_data = processRun('TREFHT',1,'GAUSS-LOWER-1.0','203501-206912.nc',[2035 2069],[1 2 3]);
% T_base = mean(ssp245_T_data.ensemble_monthly_average);
ssp245_T_data_above_base = ssp245_T_data.ensemble_monthly_average - T_base;
default_T_data_above_base = default_T_data.ensemble_monthly_average - T_base;
lower_5_T_data_above_base = lower_5_T_data.ensemble_monthly_average - T_base;
T_EQ_fdbk = feedback_0N_T_data.ensemble_monthly_average - T_base;
T_15NS_fdbk = feedback_15_T_data.ensemble_monthly_average - T_base;
T_30NS_fdbk = feedback_30_T_data.ensemble_monthly_average - T_base;
T_EQ_fdbk = feedback_0N_T_data.individual_monthly_average - T_base;
T_15NS_fdbk = feedback_15_T_data.individual_monthly_average - T_base;
T_30NS_fdbk = feedback_30_T_data.individual_monthly_average - T_base;
% default_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3]);
% lower_5_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% lower_1_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-LOWER-1.0','203501-206912.nc',[2035 2069],[1 2 3]);
% feedback_30_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-30N_30S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% feedback_15_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-15N_15S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% feedback_0N_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-0N-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% % P_base = mean(ssp245_P_data.ensemble_monthly_average);
% ssp245_P_data_above_base = ssp245_P_data.ensemble_monthly_average - P_base;
% default_P_data_above_base = default_P_data.ensemble_monthly_average - P_base;
% lower_5_P_data_above_base = lower_5_P_data.ensemble_monthly_average - P_base;
% P_EQ_fdbk = feedback_0N_P_data.ensemble_monthly_average - P_base;
% P_15NS_fdbk = feedback_15_P_data.ensemble_monthly_average - P_base;
% P_30NS_fdbk = feedback_30_P_data.ensemble_monthly_average - P_base;

%% Load GAUSS data
% ssp245_AOD_data = processRun('AODVISstdn',1,'SSP245','203501-206912.nc',[2035 2069],[1 2 3]);
% default_AOD_data = processRun('AODVISstdn',1,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3]);
% lower_5_AOD_data = processRun('AODVISstdn',1,'GAUSS-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% lower_1_AOD_data = processRun('AODVISstdn',1,'GAUSS-LOWER-1.0','203501-206912.nc',[2035 2069],[1 2 3]);
ssp245_T_data = processRun('TREFHT',1,'SSP245','203501-206912.nc',[2035 2069],[1 2 3]);
default_T_data = processRun('TREFHT',1,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3]);
lower_5_T_data = processRun('TREFHT',1,'GAUSS-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
lower_1_T_data = processRun('TREFHT',1,'GAUSS-LOWER-1.0','203501-206912.nc',[2035 2069],[1 2 3]);
% ssp245_P_data = processRun('PRECT',1000*24*60*60,'SSP245','203501-206912.nc',[2035 2069],[1 2 3]);
% default_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3]);
% lower_5_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
% lower_1_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-LOWER-1.0','203501-206912.nc',[2035 2069],[1 2 3]);
% 


%% Initiate plot
colors = colormap(brewermap([],'Dark2'));
ann_time = (2035:2069)+.5;
% figure
tiley = tiledlayout(1,1);
set(gcf, 'Position', [200, 200, 1200,1200]) % Set figure size
tiley.Padding = 'tight';

% %% Plot single point T
% my_tile = nexttile(2);
% hold on
% plot(ann_time, ssp245_T_data.ensemble_annual_average,"LineWidth",2,"Color",'k')
% plot(ann_time, averageEvery2d(12,1,T_60N+T_SSP245+T_base),"LineWidth",2,"Color",colors(1,:))
% plot(ann_time, averageEvery2d(12,1,T_30N+T_SSP245+T_base),"LineWidth",2,"Color",colors(2,:))
% plot(ann_time, averageEvery2d(12,1,T_15N+T_SSP245+T_base),"LineWidth",2,"Color",colors(3,:))
% plot(ann_time, averageEvery2d(12,1,T_0N+T_SSP245+T_base),"LineWidth",2,"Color",colors(4,:))
% plot(ann_time, averageEvery2d(12,1,T_15S+T_SSP245+T_base),"LineWidth",2,"Color",colors(5,:))
% plot(ann_time, averageEvery2d(12,1,T_30S+T_SSP245+T_base),"LineWidth",2,"Color",colors(6,:))
% plot(ann_time, averageEvery2d(12,1,T_60S+T_SSP245+T_base),"LineWidth",2,"Color",colors(7,:))
% for i = 1:8
%     if i == 1
%         injection_12_60N = [];
%         for j = 1:35
%             injection_12_60N = [injection_12_60N; [0 0 4 4 4 0 0 0 0 0 0 0]'];
%         end
%         injections = [injection_12_60N zeros(420,6)];
%         color = colors(i,:);
%     elseif i == 7
%         injection_12_60S = [];
%         for j = 1:35
%             injection_12_60S = [injection_12_60S; [0 0 0 0 0 0 0 0 4 4 4 0]'];
%         end
%         injections = [zeros(420,6) injection_12_60S];
%         color = colors(i,:);
%     elseif i ==8 
%         injections = zeros(420,7);
%         color = 'k';
%     else
%         injections = [zeros(420,i-1) ones(420,1) zeros(420,7-i)];
%         color = colors(i,:);
%     end
%     all_injection_and_CO2 = [injections CO2_forcing_SSP245_month];
%     T_emu_matrix = CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2,param_AOD_all,param_T_all,pattern_T_all);
%     T_emu = globalMean(T_emu_matrix);
%     plot(ann_time, averageEvery2d(12,1,T_emu+T_base),"LineWidth",2,"Color",color,"LineStyle","--")
% end
% hold off
% % legend(sp_legend)
% ylabel("Temperature (K)")
% title("Figure 1 but British (in progress)")
all_variability = get_global_mean_ann_variability("TREFHT",35,12);

%% Plot Feedback T
colors2 = colormap(brewermap([],'Set1'));
ens_variability_1 = mean(all_variability(:,1:3),2);
ens_variability_2 = mean(all_variability(:,4:6),2);
ens_variability_3 = mean(all_variability(:,7:9),2);
ens_variability_4 = mean(all_variability(:,10:12),2);
ens_variability_1 = (all_variability(:,1:3));
ens_variability_2 = (all_variability(:,4:6));
ens_variability_3 = (all_variability(:,7:9));
ens_variability_4 = (all_variability(:,10:12));
PIT = getPreIndustrialTModel();
my_tile = nexttile(1);
hold on
% plot(ann_time, ssp245_T_data.ensemble_annual_average,"LineWidth",4,"Color",'k')
% plot(ann_time, mean(averageEvery2d(12,1,T_30NS_fdbk+T_base),2),"LineWidth",4,"Color",colors(9,:))



T_emu_SSP245 = globalMean(CIDER_pattern_from_all_injections_and_CO2([0*injection_matrix_fdbk_0N CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
T_emu_fdbk_EQ = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_0N CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
T_emu_fdbk_15NS = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_15NS CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
T_emu_fdbk_30NS = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_30NS CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
% T_emu_fdbk_60NS = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_matrix_fdbk_60NS CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
plot(ann_time, ssp245_T_data.ensemble_annual_average-PIT,"LineWidth",4,"Color",'k')
plot(ann_time, mean(averageEvery2d(12,1,T_emu_SSP245+T_base)+ens_variability_1,2)-PIT,"LineWidth",4,"Color",colors(3,:),"LineStyle","-")
plot(ann_time, mean(averageEvery2d(12,1,T_30NS_fdbk+T_base),2)-PIT,"LineWidth",4,"Color",colors(9,:))
plot(ann_time, mean(averageEvery2d(12,1,T_emu_fdbk_30NS+T_base)+ens_variability_4,2)-PIT,"LineWidth",4,"Color",colors(2,:),"LineStyle","-")

plot(ann_time, ssp245_T_data.individual_annual_average-PIT,"LineWidth",2,"Color",'k',"LineStyle",":")
plot(ann_time, averageEvery2d(12,1,T_30NS_fdbk+T_base)-PIT,"LineWidth",2,"Color",colors(9,:),"LineStyle",":")
plot(ann_time, averageEvery2d(12,1,T_emu_SSP245+T_base)+ens_variability_1-PIT,"LineWidth",2,"Color",colors(3,:),"LineStyle",":")
plot(ann_time, averageEvery2d(12,1,T_emu_fdbk_30NS+T_base)+ens_variability_4-PIT,"LineWidth",2,"Color",colors(2,:),"LineStyle",":")
plot(ann_time, ssp245_T_data.ensemble_annual_average-PIT,"LineWidth",4,"Color",'k')
plot(ann_time, mean(averageEvery2d(12,1,T_emu_SSP245+T_base)+ens_variability_1,2)-PIT,"LineWidth",4,"Color",colors(3,:),"LineStyle","-")
plot(ann_time, mean(averageEvery2d(12,1,T_30NS_fdbk+T_base),2)-PIT,"LineWidth",4,"Color",colors(9,:))
plot(ann_time, mean(averageEvery2d(12,1,T_emu_fdbk_30NS+T_base)+ens_variability_4,2)-PIT,"LineWidth",4,"Color",colors(2,:),"LineStyle","-")
hold off
legend("SSP2-4.5, CESM2-WACCM6","SSP2-4.5, CIDER","30°N+30°S Feedback for 1.0°C, CESM2-WACCM6","30°N+30°S Feedback for 1.0°C, CIDER","Location",'nw')
xlabel("Year")
ylabel("Temperature (°C above PI)")
box on 
grid on
set(gca,'FontSize', 20)

print(gcf,'-dpng',["Uncoord_Emu_Paper/Uncoord_Plots/Figure_6_" + getNow() + ".png"],'-r300')

% %% Plot GAUSS T
% my_tile = nexttile(8);
% hold on
% plot(ann_time, ssp245_T_data.ensemble_annual_average,"LineWidth",2,"Color",'k')
% plot(ann_time, default_T_data.ensemble_annual_average,"LineWidth",2,"Color",'b')
% plot(ann_time, lower_5_T_data.ensemble_annual_average,"LineWidth",2,"Color",'b')
% plot(ann_time, lower_1_T_data.ensemble_annual_average,"LineWidth",2,"Color",'b')
% plot(ann_time, averageEvery2d(12,1,T_emu_SSP245+T_base),"LineWidth",2,"Color",'k',"LineStyle","--")
% T_emu_arise = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_default_monthly CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
% T_emu_lower_5 = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_lower_5_monthly CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
% T_emu_lower_1 = globalMean(CIDER_pattern_from_all_injections_and_CO2([injection_lower_1_monthly CO2_forcing_SSP245_month], param_AOD_all,param_T_all,pattern_T_all));
% plot(ann_time, averageEvery2d(12,1,T_emu_arise+T_base),"LineWidth",2,"Color",'b',"LineStyle","--")
% plot(ann_time, averageEvery2d(12,1,T_emu_lower_5+T_base),"LineWidth",2,"Color",'b',"LineStyle","--")
% plot(ann_time, averageEvery2d(12,1,T_emu_lower_1+T_base),"LineWidth",2,"Color",'b',"LineStyle","--")
% 
% hold off

