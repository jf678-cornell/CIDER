clc
clear all
addAllPaths

AOD_weightings = [1 25];
AOD_weightings = [1 15];
reweight = 420;

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
% feedback_inj_60N_monthly = feedback_inj_60NS_monthly/2;
% feedback_inj_60S_monthly = feedback_inj_60NS_monthly/2;


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
% AOD_0N_matrix  = ncread("u-ch622_AOD.nc","AOD550_delta");
% AOD_15N_matrix = ncread("u-cd352_AOD.nc","AOD550_delta");
% AOD_15S_matrix = ncread("u-cd353_AOD.nc","AOD550_delta");
% AOD_30N_matrix = ncread("u-cd297_AOD.nc","AOD550_delta");
% AOD_30S_matrix = ncread("u-cd354_AOD.nc","AOD550_delta");
% AOD_60N_matrix = ncread("u-cm174_AOD.nc","AOD550_delta");
% AOD_60S_matrix = ncread("u-cm175_AOD.nc","AOD550_delta");

ssp245_AOD_data = processRun('AODVISstdn',1,'SSP245','203501-206912.nc',[2035 2069],[1 2 3]);
ssp245_AOD_base = mean(ssp245_AOD_data.ensemble_annual_average);

AOD_60N_matrix = readSinglePoint(60,'N',1,'AODVISstdn');
AOD_30N_matrix = mean(readSinglePoint(30,'N',[1 2 3],'AODVISstdn'),4);
AOD_15N_matrix = mean(readSinglePoint(15,'N',[1 2 3],'AODVISstdn'),4);
AOD_0N_matrix  = mean(readSinglePoint( 0,'N',[1 2 3],'AODVISstdn'),4);
AOD_15S_matrix = mean(readSinglePoint(15,'S',[1 2 3],'AODVISstdn'),4);
AOD_30S_matrix = mean(readSinglePoint(30,'S',[1 2 3],'AODVISstdn'),4);
AOD_60S_matrix = readSinglePoint(60,'S',1,'AODVISstdn');


% T_0N_matrix = ncread("u-ch622_tas.nc","tas_delta");
% T_15N_matrix = ncread("u-cd352_tas.nc","tas_delta");
% T_15S_matrix = ncread("u-cd353_tas.nc","tas_delta");
% T_30N_matrix = ncread("u-cd297_tas.nc","tas_delta");
% T_30S_matrix = ncread("u-cd354_tas.nc","tas_delta");
% T_60N_matrix = ncread("u-cm174_tas.nc","tas_delta");
% T_60S_matrix = ncread("u-cm175_tas.nc","tas_delta");



ssp245_T_data = processRun('TREFHT',1,'SSP245','203501-206912.nc',[2015 2069],[1 2 3]);
T_SSP245_matrix = ssp245_T_data.ensemble_monthly_matrix;
T_SSP245_matrix_base = mean(T_SSP245_matrix(:,:,12*10+1:12*30),3);
% T_SSP245_matrix_base = mean(T_SSP245_matrix(:,:,12*15+1:12*24),3);
T_SSP245_matrix = T_SSP245_matrix(:,:,12*20+1:12*55) - T_SSP245_matrix_base;
pattern_T_base = T_SSP245_matrix_base;
T_60N_matrix = readSinglePoint(60,'N',1,'TREFHT')-T_SSP245_matrix- T_SSP245_matrix_base;
T_30N_matrix = mean(readSinglePoint(30,'N',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_15N_matrix = mean(readSinglePoint(15,'N',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_0N_matrix  = mean(readSinglePoint( 0,'N',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_15S_matrix = mean(readSinglePoint(15,'S',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_30S_matrix = mean(readSinglePoint(30,'S',[1 2 3],'TREFHT'),4)-T_SSP245_matrix- T_SSP245_matrix_base;
T_60S_matrix = readSinglePoint(60,'S',1,'TREFHT')-T_SSP245_matrix- T_SSP245_matrix_base;



ssp245_RH_data = processRun('RHREFHT',1,'SSP245','203501-206912.nc',[2015 2069],[1]);
RH_SSP245_matrix = ssp245_RH_data.ensemble_monthly_matrix;
RH_SSP245_matrix_base = mean(RH_SSP245_matrix(:,:,12*10+1:12*30),3);
RH_SSP245_matrix_base = mean(RH_SSP245_matrix(:,:,  12*15+1:12*24),3);
RH_SSP245_matrix = RH_SSP245_matrix(:,:,12*20+1:12*55) - RH_SSP245_matrix_base;
pattern_RH_base = RH_SSP245_matrix_base;
RH_60N_matrix = readSinglePoint(60,'N',1,'RHREFHT')-RH_SSP245_matrix- RH_SSP245_matrix_base;
RH_30N_matrix = mean(readSinglePoint(30,'N',[1 2 3],'RHREFHT'),4)-RH_SSP245_matrix- RH_SSP245_matrix_base;
RH_15N_matrix = mean(readSinglePoint(15,'N',[1 2 3],'RHREFHT'),4)-RH_SSP245_matrix- RH_SSP245_matrix_base;
RH_0N_matrix  = mean(readSinglePoint( 0,'N',[1 2 3],'RHREFHT'),4)-RH_SSP245_matrix- RH_SSP245_matrix_base;
RH_15S_matrix = mean(readSinglePoint(15,'S',[1 2 3],'RHREFHT'),4)-RH_SSP245_matrix- RH_SSP245_matrix_base;
RH_30S_matrix = mean(readSinglePoint(30,'S',[1 2 3],'RHREFHT'),4)-RH_SSP245_matrix- RH_SSP245_matrix_base;
RH_60S_matrix = readSinglePoint(60,'S',1,'RHREFHT')-RH_SSP245_matrix- RH_SSP245_matrix_base;

ssp245_U10_data = processRun('U10',1,'SSP245','203501-206912.nc',[2015 2069],[1]);
U10_SSP245_matrix = ssp245_U10_data.ensemble_monthly_matrix;
U10_SSP245_matrix_base = mean(U10_SSP245_matrix(:,:,12*10+1:12*30),3);
U10_SSP245_matrix_base = mean(U10_SSP245_matrix(:,:,  12*15+1:12*24),3);
U10_SSP245_matrix = U10_SSP245_matrix(:,:,12*20+1:12*55) - U10_SSP245_matrix_base;
pattern_U10_base = U10_SSP245_matrix_base;
U10_60N_matrix = readSinglePoint(60,'N',1,'U10')-U10_SSP245_matrix- U10_SSP245_matrix_base;
U10_30N_matrix = mean(readSinglePoint(30,'N',[1 2 3],'U10'),4)-U10_SSP245_matrix- U10_SSP245_matrix_base;
U10_15N_matrix = mean(readSinglePoint(15,'N',[1 2 3],'U10'),4)-U10_SSP245_matrix- U10_SSP245_matrix_base;
U10_0N_matrix  = mean(readSinglePoint( 0,'N',[1 2 3],'U10'),4)-U10_SSP245_matrix- U10_SSP245_matrix_base;
U10_15S_matrix = mean(readSinglePoint(15,'S',[1 2 3],'U10'),4)-U10_SSP245_matrix- U10_SSP245_matrix_base;
U10_30S_matrix = mean(readSinglePoint(30,'S',[1 2 3],'U10'),4)-U10_SSP245_matrix- U10_SSP245_matrix_base;
U10_60S_matrix = readSinglePoint(60,'S',1,'U10')-U10_SSP245_matrix- U10_SSP245_matrix_base;

pdimchange = 1000*24*60*60;
ssp245_P_data = processRun('PRECT',pdimchange,'SSP245','203501-206912.nc',[2015 2069],[1 2 3]);
P_SSP245_matrix = ssp245_P_data.ensemble_monthly_matrix;
P_SSP245_matrix_base = mean(P_SSP245_matrix(:,:,12*10+1:12*30),3);
P_SSP245_matrix_base = mean(P_SSP245_matrix(:,:,12*15:12*25),3);
P_SSP245_matrix = P_SSP245_matrix(:,:,12*20+1:12*55) - P_SSP245_matrix_base;
pattern_P_base = P_SSP245_matrix_base;

P_60N_matrix = pdimchange*readSinglePoint(60,'N',1,'PRECT')-P_SSP245_matrix- P_SSP245_matrix_base;
P_30N_matrix = mean(pdimchange*readSinglePoint(30,'N',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
P_15N_matrix = mean(pdimchange*readSinglePoint(15,'N',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
P_0N_matrix  = mean(pdimchange*readSinglePoint( 0,'N',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
P_15S_matrix = mean(pdimchange*readSinglePoint(15,'S',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
P_30S_matrix = mean(pdimchange*readSinglePoint(30,'S',[1 2 3],'PRECT'),4)-P_SSP245_matrix- P_SSP245_matrix_base;
P_60S_matrix = pdimchange*readSinglePoint(60,'S',1,'PRECT')-P_SSP245_matrix- P_SSP245_matrix_base;


qdimchange = 24*60*60;

ssp245_Q_data = processRun('QFLX',qdimchange,'SSP245','203501-206912.nc',[2015 2069],[1 2 3]);
Q_SSP245_matrix = ssp245_Q_data.ensemble_monthly_matrix;
Q_SSP245_matrix_base = mean(Q_SSP245_matrix(:,:,12*10+1:12*30),3);
Q_SSP245_matrix_base = mean(Q_SSP245_matrix(:,:,12*15:12*25),3);
Q_SSP245_matrix = Q_SSP245_matrix(:,:,12*20+1:12*55) - Q_SSP245_matrix_base;
pattern_Q_base = Q_SSP245_matrix_base;

Q_60N_matrix = qdimchange*readSinglePoint(60,'N',1,'QFLX')-Q_SSP245_matrix- Q_SSP245_matrix_base;
Q_30N_matrix = mean(qdimchange*readSinglePoint(30,'N',[1 2 3],'QFLX'),4)-Q_SSP245_matrix- Q_SSP245_matrix_base;
Q_15N_matrix = mean(qdimchange*readSinglePoint(15,'N',[1 2 3],'QFLX'),4)-Q_SSP245_matrix- Q_SSP245_matrix_base;
Q_0N_matrix  = mean(qdimchange*readSinglePoint( 0,'N',[1 2 3],'QFLX'),4)-Q_SSP245_matrix- Q_SSP245_matrix_base;
Q_15S_matrix = mean(qdimchange*readSinglePoint(15,'S',[1 2 3],'QFLX'),4)-Q_SSP245_matrix- Q_SSP245_matrix_base;
Q_30S_matrix = mean(qdimchange*readSinglePoint(30,'S',[1 2 3],'QFLX'),4)-Q_SSP245_matrix- Q_SSP245_matrix_base;
Q_60S_matrix = qdimchange*readSinglePoint(60,'S',1,'QFLX')-Q_SSP245_matrix- Q_SSP245_matrix_base;


%% Parse SP data
AOD_0N = globalMean(AOD_0N_matrix)-ssp245_AOD_base;
AOD_15N = globalMean(AOD_15N_matrix)-ssp245_AOD_base;
AOD_15S = globalMean(AOD_15S_matrix)-ssp245_AOD_base;
AOD_30N = globalMean(AOD_30N_matrix)-ssp245_AOD_base;
AOD_30S = globalMean(AOD_30S_matrix)-ssp245_AOD_base;
AOD_60N = globalMean(AOD_60N_matrix)-ssp245_AOD_base;
AOD_60S = globalMean(AOD_60S_matrix)-ssp245_AOD_base;


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
ssl2 = 420;

pattern_T_all = cat(3,CIDER_get_pattern(T_60N_matrix,ssl),CIDER_get_pattern(T_30N_matrix,ssl),CIDER_get_pattern(T_15N_matrix,ssl),CIDER_get_pattern(T_0N_matrix,ssl),CIDER_get_pattern(T_15S_matrix,ssl),CIDER_get_pattern(T_30S_matrix,ssl),CIDER_get_pattern(T_60S_matrix,ssl),CIDER_get_pattern(T_SSP245_matrix,ssl2));

RH_0N = globalMean(RH_0N_matrix);
RH_15N = globalMean(RH_15N_matrix);
RH_15S = globalMean(RH_15S_matrix);
RH_30N = globalMean(RH_30N_matrix);
RH_30S = globalMean(RH_30S_matrix);
RH_60N = globalMean(RH_60N_matrix);
RH_60S = globalMean(RH_60S_matrix);
all_RH_step_responses = [RH_60N RH_30N RH_15N RH_0N RH_15S RH_30S];

RH_SSP245 = globalMean(RH_SSP245_matrix);
ssl = 240;
pattern_RH_all = cat(3,CIDER_get_pattern(RH_60N_matrix,ssl),CIDER_get_pattern(RH_30N_matrix,ssl),CIDER_get_pattern(RH_15N_matrix,ssl),CIDER_get_pattern(RH_0N_matrix,ssl),CIDER_get_pattern(RH_15S_matrix,ssl),CIDER_get_pattern(RH_30S_matrix,ssl),CIDER_get_pattern(RH_60S_matrix,ssl),CIDER_get_pattern(RH_SSP245_matrix,ssl));

U10_0N = globalMean(U10_0N_matrix);
U10_15N = globalMean(U10_15N_matrix);
U10_15S = globalMean(U10_15S_matrix);
U10_30N = globalMean(U10_30N_matrix);
U10_30S = globalMean(U10_30S_matrix);
U10_60N = globalMean(U10_60N_matrix);
U10_60S = globalMean(U10_60S_matrix);
all_U10_step_responses = [U10_60N U10_30N U10_15N U10_0N U10_15S U10_30S];

U10_SSP245 = globalMean(U10_SSP245_matrix);
ssl = 240;
pattern_U10_all = cat(3,CIDER_get_pattern(U10_60N_matrix,ssl),CIDER_get_pattern(U10_30N_matrix,ssl),CIDER_get_pattern(U10_15N_matrix,ssl),CIDER_get_pattern(U10_0N_matrix,ssl),CIDER_get_pattern(U10_15S_matrix,ssl),CIDER_get_pattern(U10_30S_matrix,ssl),CIDER_get_pattern(U10_60S_matrix,ssl),CIDER_get_pattern(U10_SSP245_matrix,ssl));


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


Q_0N = globalMean(Q_0N_matrix);
Q_15N = globalMean(Q_15N_matrix);
Q_15S = globalMean(Q_15S_matrix);
Q_30N = globalMean(Q_30N_matrix);
Q_30S = globalMean(Q_30S_matrix);
Q_60N = globalMean(Q_60N_matrix);
Q_60S = globalMean(Q_60S_matrix);
all_Q_step_responses = [Q_60N Q_30N Q_15N Q_0N Q_15S Q_30S];

Q_SSP245 = globalMean(Q_SSP245_matrix);
ssl = 240;
pattern_Q_all = cat(3,CIDER_get_pattern(Q_60N_matrix,ssl),CIDER_get_pattern(Q_30N_matrix,ssl),CIDER_get_pattern(Q_15N_matrix,ssl),CIDER_get_pattern(Q_0N_matrix,ssl),CIDER_get_pattern(Q_15S_matrix,ssl),CIDER_get_pattern(Q_30S_matrix,ssl),CIDER_get_pattern(Q_60S_matrix,ssl),CIDER_get_pattern(Q_SSP245_matrix,ssl));

%% Load feedback data
ssp245_AOD_data = processRun('AODVISstdn',1,'SSP245','203501-206912.nc',[2035 2098],[1 2 3]);
default_AOD_data = processRun('AODVISstdn',1,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3]);
lower_5_AOD_data = processRun('AODVISstdn',1,'GAUSS-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
feedback_60_AOD_data = processRun('AODVISstdn',1,'GAUSS-60N_60S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
feedback_30_AOD_data = processRun('AODVISstdn',1,'GAUSS-30N_30S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
feedback_15_AOD_data = processRun('AODVISstdn',1,'GAUSS-15N_15S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
feedback_0N_AOD_data = processRun('AODVISstdn',1,'GAUSS-0N-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3]);
AOD_base = mean(ssp245_AOD_data.ensemble_monthly_average);
pattern_AOD_base = mean(ssp245_AOD_data.ensemble_monthly_matrix,3);
ssp245_AOD_data_above_base = ssp245_AOD_data.ensemble_monthly_average - AOD_base;
default_AOD_data_above_base = default_AOD_data.ensemble_monthly_average - AOD_base;
lower_5_AOD_data_above_base = lower_5_AOD_data.ensemble_monthly_average - AOD_base;
AOD_EQ_fdbk = feedback_0N_AOD_data.ensemble_monthly_average - AOD_base;
AOD_15NS_fdbk = feedback_15_AOD_data.ensemble_monthly_average - AOD_base;
AOD_30NS_fdbk = feedback_30_AOD_data.ensemble_monthly_average - AOD_base;
AOD_60NS_fdbk = feedback_60_AOD_data.ensemble_monthly_average - AOD_base;


ratio_60N = mean(AOD_60N(181:end))/(mean(AOD_60N(181:end))+mean(AOD_60S(181:end)));
ratio_30N = mean(AOD_30N(181:end))/(mean(AOD_30N(181:end))+mean(AOD_30S(181:end)));
ratio_15N = mean(AOD_15N(181:end))/(mean(AOD_15N(181:end))+mean(AOD_15S(181:end)));
all_AOD_fdbk_responses = 1/2*[AOD_60NS_fdbk AOD_30NS_fdbk AOD_15NS_fdbk 2*AOD_EQ_fdbk AOD_15NS_fdbk AOD_30NS_fdbk];
all_AOD_fdbk_responses = [ratio_60N*AOD_60NS_fdbk ratio_30N*AOD_30NS_fdbk ratio_15N*AOD_15NS_fdbk AOD_EQ_fdbk (1-ratio_15N)*AOD_15NS_fdbk (1-ratio_30N)*AOD_30NS_fdbk];

%% Train AOD
weightings = AOD_weightings;
all_step_injections_short = all_step_injections(1:60,:);
all_AOD_step_responses_short = all_AOD_step_responses(1:60,:);
% all_params_AOD_except_60S = CIDER_train_AOD_params(all_step_injections,all_AOD_step_responses,all_feedback_injections,all_AOD_fdbk_responses,weightings);
all_params_AOD_except_60S = CIDER_train_AOD_params(all_step_injections,all_AOD_step_responses,all_feedback_injections,all_AOD_fdbk_responses,weightings,all_step_injections_short,all_AOD_step_responses_short);
params_AOD_60S = CIDER_train_AOD_params(injection_12_60S,AOD_60S,feedback_inj_60S_monthly,AOD_60NS_fdbk*(1-ratio_60N),weightings);
param_AOD_all = [all_params_AOD_except_60S;params_AOD_60S];
% plot(CIDER_AOD_from_injection(param_AOD_all(4,:),injection_12))

%% Train T
param_T_all_inj_except_60S = CIDER_train_climate_params(all_AOD_step_responses,all_T_step_responses);
param_T_all_inj_except_60S = CIDER_train_climate_params(all_AOD_step_responses(1:reweight,:),all_T_step_responses(1:reweight,:));
AOD_60S_emu = CIDER_AOD_from_injection(param_AOD_all(7,:),injection_12_60S);
param_T_60S = CIDER_train_climate_params(AOD_60S ,T_60S);
param_T_CO2 = CIDER_train_climate_params( CO2_forcing_SSP245_month,T_SSP245)
param_T_CO2_hist = CIDER_train_climate_params_with_history( CO2_forcing_SSP245_month,T_SSP245)
param_T_all = [param_T_all_inj_except_60S;param_T_60S;param_T_CO2];
%% Train RH
param_RH_all_inj_except_60S = CIDER_train_climate_params(all_AOD_step_responses,all_RH_step_responses);
AOD_60S_emu = CIDER_AOD_from_injection(param_AOD_all(7,:),injection_12_60S);
param_RH_60S = CIDER_train_climate_params(AOD_60S ,RH_60S);
param_RH_CO2 = CIDER_train_climate_params( CO2_forcing_SSP245_month,RH_SSP245);
param_RH_all = [param_RH_all_inj_except_60S;param_RH_60S;param_RH_CO2];
%% Train U10
param_U10_all_inj_except_60S = CIDER_train_climate_params(all_AOD_step_responses,all_U10_step_responses);
AOD_60S_emu = CIDER_AOD_from_injection(param_AOD_all(7,:),injection_12_60S);
param_U10_60S = CIDER_train_climate_params(AOD_60S ,U10_60S);
param_U10_CO2 = CIDER_train_climate_params( CO2_forcing_SSP245_month,U10_SSP245);
param_U10_all = [param_U10_all_inj_except_60S;param_U10_60S;param_U10_CO2];
%% Train P
param_P_all_inj_except_60S = CIDER_train_climate_params(all_AOD_step_responses,all_P_step_responses,[240 2e-4],[0 -1],[3000 1]);
param_P_60S = CIDER_train_climate_params(AOD_60S ,P_60S,[240 2e-4],[0 -1],[3000 1]);
param_P_CO2 = CIDER_train_climate_params(CO2_forcing_SSP245_month,P_SSP245,[240 2e-4],[0 -1],[3000 1]);
param_P_all = [param_P_all_inj_except_60S;param_P_60S;param_P_CO2];

%% Train Q
param_Q_all_inj_except_60S = CIDER_train_climate_params(all_AOD_step_responses,all_Q_step_responses,[240 2e-4],[-3000 -1],[3000 1]);
param_Q_60S = CIDER_train_climate_params(AOD_60S ,Q_60S,[240 .002],[0 -1],[3000 1]);
param_Q_CO2 = CIDER_train_climate_params(CO2_forcing_SSP245_month,Q_SSP245,[120 0.002],[0 -1],[3000 1]);
param_Q_all = [param_Q_all_inj_except_60S;param_Q_60S;param_Q_CO2];
%%
T_base = globalMean(pattern_T_base);
P_base = globalMean(pattern_P_base);
Q_base = globalMean(pattern_Q_base);
RH_base = globalMean(pattern_RH_base);
U10_base = globalMean(pattern_U10_base);
AOD_base = globalMean(pattern_AOD_base);
save("Uncoord_Emu_Paper/CESM_params.mat","param_AOD_all","param_P_all","param_T_all","param_Q_all","param_RH_all","param_U10_all","pattern_P_all","pattern_T_all","pattern_Q_all","pattern_RH_all","pattern_U10_all","pattern_AOD_all","pattern_AOD_base","pattern_T_base","pattern_P_base","pattern_Q_base","pattern_U10_base","pattern_RH_base");
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