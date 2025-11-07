%% Intro
fullclear
addAllPaths
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
load CO2_concentrations.mat
load CESM_params.mat
T_base = globalMean(pattern_T_base);
P_base = globalMean(pattern_P_base);
CO2levels_2020_2100_ssp245 = CO2_SSP245(6:86);
CO2levels_2035_2070_ssp245 = CO2_SSP245(6+15:86-31);
CO2levels_2035_2100_ssp245 = CO2_SSP245(6+15:86);

CO2levels_2020_2100_ssp126 = CO2_SSP126(6:86);
CO2levels_2035_2070_ssp126 = CO2_SSP126(6+15:86-31);
CO2levels_2035_2100_ssp126 = CO2_SSP126(6+15:86);

CO2levels_2035_2100_constant = zeros(1000,1)+ CO2_SSP245(6);
CO2levels_2035_2100_held = [CO2levels_2035_2100_ssp245;zeros(1000,1)];
CO2levels_2035_2100_held(21:end) = CO2levels_2035_2100_held(20);


CO2_ref = CO2levels_2035_2070_ssp245(1);
CO2_ref = CO2_SSP245(6+14);
CO2_forcing_SSP245 = 5.35*log((CO2levels_2035_2070_ssp245)/CO2_ref);
CO2_forcing_SSP245_month = repeatElements(CO2_forcing_SSP245,12);
years_of_all_forcing = [2005	2010	2020	2030	2040	2050	2060	2070	2080	2090	2100];
all_forcing =          [1.871	2.137	2.622	3.017	3.470	3.922	4.395	4.897	5.421	5.983	6.561];
CO2_forcing_SSP245_month = interp1(years_of_all_forcing,all_forcing,annualToMonthly(2035:2069))-interp1(years_of_all_forcing,all_forcing,2035);

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
injection_lower_5_1 = [L_LOWER_5.S30(1:end-1,1) L_LOWER_5.S15(1:end-1,1) L_LOWER_5.N15(1:end-1,1) L_LOWER_5.N30(1:end-1,1)];
injection_lower_5_2 = [L_LOWER_5.S30(1:end-1,2) L_LOWER_5.S15(1:end-1,2) L_LOWER_5.N15(1:end-1,2) L_LOWER_5.N30(1:end-1,2)];
injection_lower_5_3 = [L_LOWER_5.S30(1:end-1,3) L_LOWER_5.S15(1:end-1,3) L_LOWER_5.N15(1:end-1,3) L_LOWER_5.N30(1:end-1,3)];

injection_lower_5_monthly = zeros(length(injection_lower_5(1:end-1,1))*12,7);
injection_lower_5_monthly_1 = injection_lower_5_monthly;
injection_lower_5_monthly_2 = injection_lower_5_monthly;
injection_lower_5_monthly_3 = injection_lower_5_monthly;
for i = 1:2
    injection_lower_5_monthly_1(:,i+1) = repeatElements(injection_lower_5_1(:,5-i),12)/12;
    injection_lower_5_monthly_2(:,i+1) = repeatElements(injection_lower_5_2(:,5-i),12)/12;
    injection_lower_5_monthly_3(:,i+1) = repeatElements(injection_lower_5_3(:,5-i),12)/12;

end
for i = 3:4
    injection_lower_5_monthly_1(:,i+2) = repeatElements(injection_lower_5_1(:,5-i),12)/12;
    injection_lower_5_monthly_2(:,i+2) = repeatElements(injection_lower_5_2(:,5-i),12)/12;
    injection_lower_5_monthly_3(:,i+2) = repeatElements(injection_lower_5_3(:,5-i),12)/12;
end


injection_lower_5_monthly(:,2) = repeatElements(mean(L_LOWER_5.N30(1:end-1,:),2),12)/12;
injection_lower_5_monthly(:,3) = repeatElements(mean(L_LOWER_5.N15(1:end-1,:),2),12)/12;
injection_lower_5_monthly(:,5) = repeatElements(mean(L_LOWER_5.S15(1:end-1,:),2),12)/12;
injection_lower_5_monthly(:,6) = repeatElements(mean(L_LOWER_5.S30(1:end-1,:),2),12)/12;


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
%% injections symmetric 
load yearly_injection_rate.mat
% 0N_1.0, 15N_15S, 30N_30S, 60N_60S_1.0, Global+1C
feedback_inj_0N = inj_rate_arr(:,1);
feedback_inj_15 = inj_rate_arr(:,2);
feedback_inj_30 = inj_rate_arr(:,3);

%% Fit INJ to AOD

injection_12 = ones(420,1);
injection_12_60N = [];
injection_12_60S = [];
for i = 1:35
injection_12_60N = [injection_12_60N; [0 0 4 4 4 0 0 0 0 0 0 0]'];
injection_12_60S = [injection_12_60S; [0 0 0 0 0 0 0 0 4 4 4 0]'];
end


injection_feedback_0N_monthly_single = repeatElements(feedback_inj_0N,12)/12;
injection_feedback_0N_monthly = [zeros(420,3),repeatElements(feedback_inj_0N,12)/12,zeros(420,2)];
injection_feedback_15_monthly = [zeros(420,2),repeatElements(feedback_inj_15,12)/12/2,zeros(420,1),repeatElements(feedback_inj_15,12)/12/2,zeros(420,1)];
injection_feedback_30_monthly = [zeros(420,1),repeatElements(feedback_inj_30,12)/12/2,zeros(420,3),repeatElements(feedback_inj_30,12)/12/2];


%% Fit AOD -> Temp timescale

T_base_above_PI = T_base-getPreIndustrialTModel();

ssp245_T_data = processRun('TREFHT',1,'SSP245','203501-206912.nc',[2035 2069],[1 2 3],p);
default_T_data = processRun('TREFHT',1,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3],p);
lower_5_T_data = processRun('TREFHT',1,'GAUSS-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
lower_5_T_data_individual_annual_matrix = lower_5_T_data.individual_annual_matrix;
lower_5_T_data_individual_annual_matrix_dif = lower_5_T_data.individual_annual_matrix-ssp245_T_data.individual_annual_matrix;

% feedback_30_T_data = processRun('TREFHT',1,'GAUSS-30N_30S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
% feedback_15_T_data = processRun('TREFHT',1,'GAUSS-15N_15S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
% feedback_0N_T_data = processRun('TREFHT',1,'GAUSS-0N-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
monthlies = getMonthlies(default_T_data,2050,2069);
deseasonalize = @(monthly_data) reshape(reshape(monthly_data,[12 length(monthly_data)/12])-monthlies,[length(monthly_data) 1]);

ssp245_T_data_deseasonalized = deseasonalize(ssp245_T_data.ensemble_monthly_average-T_base);
default_T_data_deseasonalized = deseasonalize(default_T_data.ensemble_monthly_average-T_base);
lower_5_T_data_deseasonalized = deseasonalize(lower_5_T_data.ensemble_monthly_average-T_base);
% feedback_30_T_data_deseasonalized = deseasonalize(feedback_30_T_data.ensemble_monthly_average-T_base);
% feedback_15_T_data_deseasonalized = deseasonalize(feedback_15_T_data.ensemble_monthly_average-T_base);
% feedback_0N_T_data_deseasonalized = deseasonalize(feedback_0N_T_data.ensemble_monthly_average-T_base);

ensemble_mean_sp = @(lat,NS,var) (readSinglePoint(lat,NS,1,var)+readSinglePoint(lat,NS,2,var))/2;


default_T_pattern = ssp245_T_data.individual_annual_matrix-default_T_data.individual_annual_matrix;
default_T_pattern = mean(mean(default_T_pattern(:,:,16:end,:),4),3)/globalMean(mean(mean(default_T_pattern(:,:,16:end,:),4),3));

% T_0N_2 = globalMean(readSinglePoint(0,'N',2,'TREFHT'),p);
% T_0N_1 = globalMean(readSinglePoint(0,'N',1,'TREFHT'),p);
% T_0N = deseasonalize(mean([T_0N_1 T_0N_2],2)-T_base);
% 
% T_15N_2 = globalMean(readSinglePoint(15,'N',2,'TREFHT'),p);
% T_15N_1 = globalMean(readSinglePoint(15,'N',1,'TREFHT'),p);
% T_15N = deseasonalize(mean([T_15N_1 T_15N_2],2)-T_base);
% 
% T_15S_2 = globalMean(readSinglePoint(15,'S',2,'TREFHT'),p);
% T_15S_1 = globalMean(readSinglePoint(15,'S',1,'TREFHT'),p);
% T_15S = deseasonalize(mean([T_15S_1 T_15S_2],2)-T_base);
% 
% T_30N_2 = globalMean(readSinglePoint(30,'N',2,'TREFHT'),p);
% T_30N_1 = globalMean(readSinglePoint(30,'N',1,'TREFHT'),p);
% T_30N = deseasonalize(mean([T_30N_1 T_30N_2],2)-T_base);
% 
% T_30S_2 = globalMean(readSinglePoint(30,'S',2,'TREFHT'),p);
% T_30S_1 = globalMean(readSinglePoint(30,'S',1,'TREFHT'),p);
% T_30S = deseasonalize(mean([T_30S_1 T_30S_2],2)-T_base);
% 
% T_60N = deseasonalize(globalMean(readSinglePoint(60,'N',1,'TREFHT'),p)-T_base);
% T_60S = deseasonalize(globalMean(readSinglePoint(60,'S',1,'TREFHT'),p)-T_base);

weighting = [1 1 1];
funSSP = @(params) weighting(1)*optimize_params_forcing_to_T(params,ssp245_T_data_deseasonalized(ssp245_T_data.months>2035&ssp245_T_data.months<2070),CO2_forcing_SSP245_month,ssp245_AOD_data_above_base(1:end));
funSSPsinglepoint = @(params,AOD,T) optimize_params_forcing_to_T(params,T,CO2_forcing_SSP245_month,AOD)*weighting(2)+funSSP(params);

clear ssp245_T_data default_T_data lower_5_T_data   
%% Fit AOD -> Temp timescale

P_base_above_PI = 0;

ssp245_P_data = processRun('PRECT',1000*24*60*60,'SSP245','203501-206912.nc',[2035 2069],[1 2 3],p);
default_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-DEFAULT','203501-206912.nc',[2035 2069],[1 2 3],p);
lower_5_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
lower_5_P_data_individual_annual_matrix = lower_5_P_data.individual_annual_matrix;
lower_5_P_data_individual_annual_matrix_dif = lower_5_P_data.individual_annual_matrix-ssp245_P_data.individual_annual_matrix;
% feedback_30_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-30N_30S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
% feedback_15_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-15N_15S-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
% feedback_0N_P_data = processRun('PRECT',1000*24*60*60,'GAUSS-0N-LOWER-0.5','203501-206912.nc',[2035 2069],[1 2 3],p);
monthlies = getMonthlies(default_P_data,2050,2069);
deseasonalize = @(monthly_data) reshape(reshape(monthly_data,[12 length(monthly_data)/12])-monthlies,[length(monthly_data) 1]);

% ssp245_P_data_deseasonalized = deseasonalize(ssp245_P_data.ensemble_monthly_average-P_base);
% default_P_data_deseasonalized = deseasonalize(default_P_data.ensemble_monthly_average-P_base);
lower_5_P_data_deseasonalized = deseasonalize(lower_5_P_data.ensemble_monthly_average-P_base);
% feedback_30_P_data_deseasonalized = deseasonalize(feedback_30_P_data.ensemble_monthly_average-P_base);
% feedback_15_P_data_deseasonalized = deseasonalize(feedback_15_P_data.ensemble_monthly_average-P_base);
% feedback_0N_P_data_deseasonalized = deseasonalize(feedback_0N_P_data.ensemble_monthly_average-P_base);


default_P_pattern = ssp245_P_data.individual_annual_matrix-default_P_data.individual_annual_matrix;
default_P_pattern = mean(mean(default_P_pattern(:,:,16:end,:),4),3)/globalMean(mean(mean(default_P_pattern(:,:,16:end,:),4),3));

ensemble_mean_sp = @(lat,NS,var) (readSinglePoint(lat,NS,1,var)+readSinglePoint(lat,NS,2,var))/2;
% 
% P_0N_2 = globalMean(readSinglePoint(0,'N',2,'PRECT'),p)*1000*24*60*60;
% P_0N_1 = globalMean(readSinglePoint(0,'N',1,'PRECT'),p)*1000*24*60*60;
% P_0N = deseasonalize(mean([P_0N_1 P_0N_2],2)-P_base);
% 
% P_15N_2 = globalMean(readSinglePoint(15,'N',2,'PRECT'),p)*1000*24*60*60;
% P_15N_1 = globalMean(readSinglePoint(15,'N',1,'PRECT'),p)*1000*24*60*60;
% P_15N = deseasonalize(mean([P_15N_1 P_15N_2],2)-P_base);
% 
% P_15S_2 = globalMean(readSinglePoint(15,'S',2,'PRECT'),p)*1000*24*60*60;
% P_15S_1 = globalMean(readSinglePoint(15,'S',1,'PRECT'),p)*1000*24*60*60;
% P_15S = deseasonalize(mean([P_15S_1 P_15S_2],2)-P_base);
% 
% P_30N_2 = globalMean(readSinglePoint(30,'N',2,'PRECT'),p)*1000*24*60*60;
% P_30N_1 = globalMean(readSinglePoint(30,'N',1,'PRECT'),p)*1000*24*60*60;
% P_30N = deseasonalize(mean([P_30N_1 P_30N_2],2)-P_base);
% 
% P_30S_2 = globalMean(readSinglePoint(30,'S',2,'PRECT'),p)*1000*24*60*60;
% P_30S_1 = globalMean(readSinglePoint(30,'S',1,'PRECT'),p)*1000*24*60*60;
% P_30S = deseasonalize(mean([P_30S_1 P_30S_2],2)-P_base);
% 
% P_60N = deseasonalize(globalMean(readSinglePoint(60,'N',1,'PRECT'),p)*1000*24*60*60-P_base);
% P_60S = deseasonalize(globalMean(readSinglePoint(60,'S',1,'PRECT'),p)*1000*24*60*60-P_base);

% weighting = [1 1 1];
% funSSP = @(params) weighting(1)*optimize_params_forcing_to_T(params,ssp245_T_data_deseasonalized(ssp245_T_data.months>2035&ssp245_T_data.months<2070),CO2_forcing_SSP245_month,ssp245_AOD_data_above_base(1:end));
% funSSPsinglepoint = @(params,AOD,T) optimize_params_forcing_to_T(params,T,CO2_forcing_SSP245_month,AOD)*weighting(2)+funSSP(params);

clear  default_P_data lower_5_P_data   
%%
all_injection_and_CO2 = [injection_lower_5_monthly CO2_forcing_SSP245_month];
all_injection_and_CO2_1 = [injection_lower_5_monthly_1 0*CO2_forcing_SSP245_month];
all_injection_and_CO2_2 = [injection_lower_5_monthly_2 0*CO2_forcing_SSP245_month];
all_injection_and_CO2_3 = [injection_lower_5_monthly_3 0*CO2_forcing_SSP245_month];
% patternscaled_multiobj_T_all = pattern_scale_emulate(param_AOD_all,param_T_all,injection_lower_5_monthly,CO2_forcing_SSP245_month,all_T_patterns_scaled,T_base_pattern);
patternscaled_multiobj_T = averageEvery(12,1,CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2,param_AOD_all,param_T_all,pattern_T_all)+0*pattern_T_base);
patternscaled_multiobj_T_1 = averageEvery(12,1,CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2_1,param_AOD_all,param_T_all,pattern_T_all)+0*pattern_T_base);
patternscaled_multiobj_T_2 = averageEvery(12,1,CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2_2,param_AOD_all,param_T_all,pattern_T_all)+0*pattern_T_base);
patternscaled_multiobj_T_3 = averageEvery(12,1,CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2_3,param_AOD_all,param_T_all,pattern_T_all)+0*pattern_T_base);
patternscaled_multiobj_T = cat(4,patternscaled_multiobj_T_1,patternscaled_multiobj_T_2,patternscaled_multiobj_T_3);
% patternscaled_multiobj_T = patternscaled_multiobj_T_all.sum_change;
patternscaled_multiobj_T_def_scale = patternscaled_multiobj_T;
for i = 1:3
    patternscaled_multiobj_T_def_scale(:,:,:,i) = CIDER_pattern_scale(globalMean(patternscaled_multiobj_T(:,:,:,i)),default_T_pattern);
end


lower_steady_state_T = cat(3,lower_5_T_data_individual_annual_matrix_dif(:,:,16:end,1),lower_5_T_data_individual_annual_matrix_dif(:,:,16:end,2),lower_5_T_data_individual_annual_matrix_dif(:,:,16:end,3));
grid_mean_T = mean(lower_steady_state_T,3);
for i = 1:288
    for j = 1:192
        grid_ss_T = squeeze(lower_steady_state_T(i,j,:));
        grid_deviations_T = grid_ss_T - grid_mean_T(i,j);
        % grid_deviations_reshaped_T = reshape(grid_deviations_T,12,[]);
        % grid_monthlies_T = mean(grid_deviations_reshaped_T,2);
        % grid_deviations_reshaped_deseasonalized_T = grid_deviations_reshaped_T-grid_monthlies_T;
        % grid_deviations_deseasonalized_T = reshape(grid_deviations_reshaped_deseasonalized_T,1,[]);
        % grid_std_T(i,j) = sqrt(sum((grid_deviations_deseasonalized_T).^2)/(length(grid_ss_T)-1));
        grid_std_T(i,j) = sqrt(sum((grid_deviations_T).^2)/(length(grid_ss_T)-1));
        grid_ste_T(i,j) = sqrt(sum((grid_deviations_T).^2))/(length(grid_ss_T)-1);
    end
end

grid_std_T = getPIstd("TREFHT",1)*sqrt(2);
grid_ste_T = grid_std_T/sqrt(60);
difference_from_MO_T =lower_5_T_data_individual_annual_matrix_dif -patternscaled_multiobj_T;
difference_from_MO_deseasonalized_T = -0*(patternscaled_multiobj_T-lower_5_T_data_individual_annual_matrix_dif);
difference_from_MO_T_def_scale =lower_5_T_data_individual_annual_matrix_dif -patternscaled_multiobj_T_def_scale;

% aa_T = sqrt(mean(mean(globalMean(difference_from_MO_T(:,:,16:end,:).^2,p),2),1));
% bb_T = sqrt(globalMean(mean(mean(difference_from_MO_T(:,:,16:end,:),3),4).^2,p));

% for i = 1:288
%     for j = 1:192
%         for k = 1:3
%             grid_tseries_T = lower_5_T_data_individual_annual_matrix(i,j,:,k);
%             grid_tseries_reshaped_T = reshape(grid_tseries_T,12,[]);
%             grid_tseries_reshaped_ss_T = grid_tseries_reshaped_T(:,16:end);
%             monthlies_T = mean(grid_tseries_reshaped_T,2)-mean(mean(grid_tseries_reshaped_T,2),1);
%             grid_tseries_reshaped_deseasonalized_T = grid_tseries_reshaped_T-monthlies_T;
%             grid_tseries_deseasonalized_T = reshape(grid_tseries_reshaped_deseasonalized_T,1,[]);
%             difference_from_MO_deseasonalized_T(i,j,:,k) = squeeze(patternscaled_multiobj_T(i,j,:,1))-grid_tseries_deseasonalized_T';
%         end
%     end
% end % Victorian Albert museum
% difference_from_MO_rmse_T = mean(mean(difference_from_MO_deseasonalized_T.^2,4),3).^0.5;
difference_from_MO_std_T = mean(mean(difference_from_MO_T./grid_std_T,4),3);
difference_from_MO_std_T = mean(mean(difference_from_MO_T(:,:,16:end,:)./grid_std_T,4),3);
difference_from_MO_regular_T = mean(mean(difference_from_MO_T(:,:,16:end,:),4),3);
difference_from_MO_ste_T_def = mean(mean(difference_from_MO_T_def_scale(:,:,16:end,:)./grid_ste_T,4),3);

% difference_from_MO_ste_T = mean(mean(difference_from_MO_T(:,:,16:end,:)./grid_std_T,4),3).*sqrt(60);
difference_from_MO_ste_T = mean(mean(difference_from_MO_T(:,:,16:end,:)./grid_ste_T,4),3);
% difference_from_MO_ste_T = mean(mean(difference_from_MO_T(:,:,16:end,:)./grid_std_T,4),3).*sqrt(60-1);
difference_from_MO_pct_T = mean(mean(difference_from_MO_T(:,:,16:end,:)./mean(mean(lower_5_T_data_individual_annual_matrix_dif(:,:,16:end,:),4),3),4),3)*100;
% difference_from_MO_ste = mean(mean(difference_from_MO_rmse(:,:,:,:)./grid_std,4),3).*sqrt(20);
difference_from_MO_to_plot = mean(mean(difference_from_MO_T,4),3);
lat_to_stipple_T = [];
lon_to_stipple_T = [];
please_stipple_T = zeros(288,192);
for j = 1:192
    for i = 1:288
        please_stipple_T(i,j) = abs(difference_from_MO_ste_T(i,j))>2;
        if please_stipple_T(i,j)
            lat_to_stipple_T = [lat_to_stipple_T;lat(j)];
            lon_to_stipple_T = [lon_to_stipple_T;lon(i)];
        end
    end
end

% save("Uncoord_Emu_Paper/Testing_Vars","patternscaled_multiobj_T_all","patternscaled_multiobj_T","param_AOD_all","param_T_all","injection_lower_5_monthly","CO2_forcing_SSP245_month","all_T_patterns_scaled","T_base_pattern")

difference_from_MO_regular_T = [difference_from_MO_regular_T;difference_from_MO_regular_T(1,:)];
difference_from_MO_ste_T = [difference_from_MO_ste_T;difference_from_MO_ste_T(1,:)];
difference_from_MO_pct_T = [difference_from_MO_pct_T;difference_from_MO_pct_T(1,:)];
difference_from_MO_ste_T_def = [difference_from_MO_ste_T_def;difference_from_MO_ste_T_def(1,:)];

% difference_from_MO_rmse_T = [difference_from_MO_rmse_T;difference_from_MO_rmse_T(1,:)];
please_stipple_T = [please_stipple_T;please_stipple_T(1,:)];



%%
% patternscaled_multiobj_P = pattern_scale_emulate(param_AOD_all,param_P_all,injection_lower_5_monthly,CO2_forcing_SSP245_month,all_P_patterns_scaled,P_base_pattern);
% patternscaled_multiobj_P = patternscaled_multiobj_P.sum_change;
patternscaled_multiobj_P = averageEvery(12,1,CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2,param_AOD_all,param_P_all,pattern_P_all)+0*pattern_P_base);
patternscaled_multiobj_P_1 = averageEvery(12,1,CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2_1,param_AOD_all,param_P_all,pattern_P_all)+0*pattern_P_base);
patternscaled_multiobj_P_2 = averageEvery(12,1,CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2_2,param_AOD_all,param_P_all,pattern_P_all)+0*pattern_P_base);
patternscaled_multiobj_P_3 = averageEvery(12,1,CIDER_pattern_from_all_injections_and_CO2(all_injection_and_CO2_3,param_AOD_all,param_P_all,pattern_P_all)+0*pattern_P_base);
patternscaled_multiobj_P = cat(4,patternscaled_multiobj_P_1,patternscaled_multiobj_P_2,patternscaled_multiobj_P_3);
%
patternscaled_multiobj_P_def_scale = patternscaled_multiobj_P;
for i = 1:3
    patternscaled_multiobj_P_def_scale(:,:,:,i) = CIDER_pattern_scale(globalMean(patternscaled_multiobj_P(:,:,:,i)),default_P_pattern);
end

lower_steady_state_P = cat(3,lower_5_P_data_individual_annual_matrix_dif(:,:,16:end,1),lower_5_P_data_individual_annual_matrix_dif(:,:,16:end,2),lower_5_P_data_individual_annual_matrix_dif(:,:,16:end,3));
grid_mean_P = mean(lower_steady_state_P,3);
for i = 1:288
    for j = 1:192
        grid_ss_P = squeeze(lower_steady_state_P(i,j,:));
        grid_deviations_P = grid_ss_P - grid_mean_P(i,j);
        % grid_deviations_reshaped_P = reshape(grid_deviations_P,12,[]);
        % grid_monthlies_P = mean(grid_deviations_reshaped_P,2);
        % grid_deviations_reshaped_deseasonalized_P = grid_deviations_reshaped_P-grid_monthlies_P;
        % grid_deviations_deseasonalized_P = reshape(grid_deviations_reshaped_deseasonalized_P,1,[]);
        % grid_std_P(i,j) = sqrt(sum((grid_deviations_deseasonalized_P).^2)/(length(grid_ss_P)-1));
        grid_std_P(i,j) = sqrt(sum((grid_deviations_P).^2)/(length(grid_ss_P)-1));
        grid_ste_P(i,j) = sqrt(sum((grid_deviations_P).^2))/(length(grid_ss_P)-1);
    end
end

grid_std_P = getPIstd("PRECT",1000*24*60*60)*sqrt(2);
grid_ste_P = grid_std_P/sqrt(60);

difference_from_MO_P = -(patternscaled_multiobj_P-lower_5_P_data_individual_annual_matrix_dif);
difference_from_MO_P_def_scale =lower_5_P_data_individual_annual_matrix_dif - patternscaled_multiobj_P_def_scale;

difference_from_MO_deseasonalized_P = 0*(patternscaled_multiobj_P-lower_5_P_data_individual_annual_matrix_dif);
difference_from_MO_ste_P_def = mean(mean(difference_from_MO_P_def_scale(:,:,16:end,:)./grid_ste_P,4),3);

aa_P = sqrt(mean(mean(globalMean(difference_from_MO_P(:,:,16:end,:).^2,p),2),1));
bb_P = sqrt(globalMean(mean(mean(difference_from_MO_P(:,:,16:end,:),3),4).^2,p));

% for i = 1:288
%     for j = 1:192
%         for k = 1:3
%             grid_tseries_P = lower_5_P_data_individual_annual_matrix(i,j,:,k);
%             grid_tseries_reshaped_P = reshape(grid_tseries_P,12,[]);
%             grid_tseries_reshaped_ss_P = grid_tseries_reshaped_P(:,16:end);
%             monthlies_P = mean(grid_tseries_reshaped_P,2)-mean(mean(grid_tseries_reshaped_P,2),1);
%             grid_tseries_reshaped_deseasonalized_P = grid_tseries_reshaped_P-monthlies_P;
%             grid_tseries_deseasonalized_P = reshape(grid_tseries_reshaped_deseasonalized_P,1,[]);
%             difference_from_MO_deseasonalized_P(i,j,:,k) = squeeze(patternscaled_multiobj_P(i,j,:,1))-grid_tseries_deseasonalized_P';
%         end
%     end
% end 
% difference_from_MO_rmse_P = mean(mean(difference_from_MO_deseasonalized_P.^2,4),3).^0.5;
difference_from_MO_std_P = mean(mean(difference_from_MO_P./grid_std_P,4),3);
difference_from_MO_std_P = mean(mean(difference_from_MO_P(:,:,16:end,:)./grid_std_P,4),3);
difference_from_MO_regular_P = mean(mean(difference_from_MO_P(:,:,16:end,:),4),3);

% difference_from_MO_ste_P = mean(mean(difference_from_MO_P(:,:,16:end,:)./grid_std_P,4),3).*sqrt(60);
difference_from_MO_ste_P = mean(mean(difference_from_MO_P(:,:,16:end,:)./grid_ste_P,4),3);
% difference_from_MO_ste_P = mean(mean(difference_from_MO_P(:,:,16:end,:)./grid_std_P,4),3).*sqrt(60-1);
difference_from_MO_pct_P = mean(mean(difference_from_MO_P(:,:,16:end,:)./mean(mean(lower_5_P_data_individual_annual_matrix_dif(:,:,16:end,:),4),3),4),3)*100;
% difference_from_MO_ste = mean(mean(difference_from_MO_rmse(:,:,:,:)./grid_std,4),3).*sqrt(20);
difference_from_MO_to_plot = mean(mean(difference_from_MO_P,4),3);
lat_to_stipple_P = [];
lon_to_stipple_P = [];
please_stipple_P = zeros(288,192);
for j = 1:192
    for i = 1:288
        please_stipple_P(i,j) = abs(difference_from_MO_ste_P(i,j))>2;
        if please_stipple_P(i,j)
            lat_to_stipple_P = [lat_to_stipple_P;lat(j)];
            lon_to_stipple_P = [lon_to_stipple_P;lon(i)];
        end
    end
end


difference_from_MO_regular_P = [difference_from_MO_regular_P;difference_from_MO_regular_P(1,:)];
difference_from_MO_ste_P = [difference_from_MO_ste_P;difference_from_MO_ste_P(1,:)];
difference_from_MO_pct_P = [difference_from_MO_pct_P;difference_from_MO_pct_P(1,:)];
difference_from_MO_ste_P_def= [difference_from_MO_ste_P_def;difference_from_MO_ste_P_def(1,:)];
% difference_from_MO_rmse_P = [difference_from_MO_rmse_P;difference_from_MO_rmse_P(1,:)];
please_stipple_P = [please_stipple_P;please_stipple_P(1,:)];


%%
[LAT,LON] = meshgrid(lat,[lon;360]);
injection_all_60N = [injection_12_60N,zeros(420,6)];
injection_all = injection_default_monthly;
time = annualToMonthly(2035:2069);
create_SP = @(i) [zeros(420,i-1),ones(420,1),zeros(420,7-i)];
create_SP_pulse = @(i) [zeros(420,i-1),injection_12_60N,zeros(420,7-i)];
injection_feedback_0N_monthly = [zeros(420,3),repeatElements(feedback_inj_0N,12)/12,zeros(420,3)];
injection_feedback_15_monthly = [zeros(420,2),repeatElements(feedback_inj_15,12)/12/2,zeros(420,1),repeatElements(feedback_inj_15,12)/12/2,zeros(420,2)];
injection_feedback_30_monthly = [zeros(420,1),repeatElements(feedback_inj_30,12)/12/2,zeros(420,3),repeatElements(feedback_inj_30,12)/12/2,zeros(420,1)];
T_base_above_PI = T_base-getPreIndustrialTModel();

CO2_forcing_SSP245_month = repeatElements(CO2_forcing_SSP245,12);
colors = colormap(brewermap([],'Dark2'));
% figure
% tiley = tiledlayout(1,2);
tiley = tiledlayout(3,2);
% set(gcf, 'Position', [200, 200, 1800,1200]) % Set figure size
set(gcf, 'Position', [200, 200, 1200,675]) % Set figure size
tiley.TileSpacing = 'compact';
tiley.Padding = 'compact';

ax3 = nexttile;

load coastlines
worldmap('World');
box on
hold on
mmin = -1;
mmax = 1;
l_colb = 60;
key_for_color = linspace(mmin,mmax,l_colb);
fc = brewermap(l_colb,'*RdBu');
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
% fc(key_for_color<2&key_for_color>-2,:) = brighten(fc(key_for_color<2&key_for_color>-2,:),.5);
number_of_levels=6;
v2 = mmin:(mmax-mmin)/(l_colb/5):mmax;
hold on
% contourfm(lat,[lon;360],double(difference_from_MO_ste_T'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
% contourfm(lat,[lon;360],double(difference_from_MO_pct_T'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
% contourfm(lat,[lon;360],double((difference_from_MO_ste.*please_stipple)'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
contourfm(lat,[lon;360],double(difference_from_MO_regular_T'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
colormap(ax3,fc)

geoshow(coastlat,coastlon,'Color','k')
% scatterm(lat_to_stipple_T,lon_to_stipple_T,4,'k','filled')
mlabel off; plabel off; gridm off
hold off
caxis([mmin mmax])
cbmap1=flipud(cbrewer('div','*RdBu',number_of_levels));
hl = colorbar('YTick',[-1  0  1]);
hl.FontSize = 16;
title("(a) CIDER, T Error, °C")
set(gca,"Fontsize",16)

ax4 = nexttile;
load coastlines
worldmap('World');
box on
hold on
mmin = -2;
mmax = 2;
l_colb = 60;
key_for_color = linspace(mmin,mmax,l_colb);
fc2 = brewermap(l_colb,'BrBG');
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
% fc2(key_for_color<2&key_for_color>-2,:) = brighten(fc2(key_for_color<2&key_for_color>-2,:),.5);
number_of_levels=6;
v2 = mmin:(mmax-mmin)/(l_colb/5):mmax;
hold on
% contourfm(lat,[lon;360],double(difference_from_MO_ste_P'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
% contourfm(lat,[lon;360],double(difference_from_MO_pct_P'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
% contourfm(lat,[lon;360],double((difference_from_MO_ste.*please_stipple)'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
contourfm(lat,[lon;360],double(difference_from_MO_regular_P'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
colormap(ax4, fc2)

geoshow(coastlat,coastlon,'Color','k')
% scatterm(lat_to_stipple_P,lon_to_stipple_P,4,'k','filled')
mlabel off; plabel off; gridm off
hold off
caxis([mmin mmax])
cbmap=flipud(cbrewer('div','BrBG',number_of_levels));
h2 = colorbar('YTick',[-2  0  2]);
h2.FontSize = 16;
title("(b) CIDER, P Error, mm/day")
set(gca,"Fontsize",16)


ax1 = nexttile;

load coastlines
worldmap('World');
box on
hold on
mmin = -3;
mmax = 3;
l_colb = 60;
key_for_color = linspace(mmin,mmax,l_colb);
fc = brewermap(l_colb,'*RdBu');
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
fc(key_for_color<2&key_for_color>-2,:) = brighten(fc(key_for_color<2&key_for_color>-2,:),.8);
% fc(key_for_color<-0.95&key_for_color>-1.05,:) = [1 0 1];
% fc(key_for_color>0.95&key_for_color<1.05,:) = [1 0 1];
% fc(key_for_color<-1.95&key_for_color>-2.05,:) = [.5 0 1];
% fc(key_for_color>1.95&key_for_color<2.05,:) = [.5 0 1];
number_of_levels=6;
v2 = mmin:(mmax-mmin)/(l_colb/5):mmax;
hold on
contourfm(lat,[lon;360],double(difference_from_MO_ste_T'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
% contourfm(lat,[lon;360],double(difference_from_MO_pct_T'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
% contourfm(lat,[lon;360],double((difference_from_MO_ste.*please_stipple)'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
%contourfm(lat,[lon;360],double(difference_from_MO_regular'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
% contourfm(lat,[lon;360],abs(double(difference_from_MO_ste_T')),[1.9 2.1],'k', 'LineWidth', 5)
colormap(ax1,fc)

geoshow(coastlat,coastlon,'Color','k')
% scatterm(lat_to_stipple_T,lon_to_stipple_T,4,'k','filled')
mlabel off; plabel off; gridm off
hold off
caxis([mmin mmax])
cbmap1=flipud(cbrewer('div','*RdBu',number_of_levels));
hl = colorbar('YTick',[-3 -2 -1 0 1 2 3]);
hl.FontSize = 16;

hold on 
% contourm(lat,[lon;360],abs(double(difference_from_MO_ste_T')),[0.99 1.01],'Color',[1 0 1], 'LineWidth', 1);
% contourm(lat,[lon;360],abs(double(difference_from_MO_ste_T')),[1.99 2.01],'Color',[.5 0 1], 'LineWidth', 1);
hold off
title("(c) CIDER, T Error, SE-Normalized")
set(gca,"Fontsize",16)
ax2 = nexttile;

load coastlines
worldmap('World');
box on
hold on
mmin = -3;
mmax = 3;
l_colb = 60;
key_for_color = linspace(mmin,mmax,l_colb);
fc2 = brewermap(l_colb,'BrBG');
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
fc2(key_for_color<2&key_for_color>-2,:) = brighten(fc2(key_for_color<2&key_for_color>-2,:),.8);

number_of_levels=6;
v2 = mmin:(mmax-mmin)/(l_colb/5):mmax;
hold on
contourfm(lat,[lon;360],double(difference_from_MO_ste_P'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
colormap(ax2, fc2)

geoshow(coastlat,coastlon,'Color','k')
mlabel off; plabel off; gridm off
hold off
% caxis([mmin mmax])
clim([mmin mmax])
cbmap=flipud(cbrewer('div','BrBG',number_of_levels));
h2 = colorbar('YTick',[-3 -2 -1 0 1 2 3]);
h2.FontSize = 16;
title("(d) CIDER, P Error, SE-Normalized")
set(gca,"Fontsize",16)
hold on
% contourm(lat,[lon;360],abs(double(difference_from_MO_ste_P')),[0.99 1.01],'Color',[1 0 1], 'LineWidth', 1);
% contourm(lat,[lon;360],abs(double(difference_from_MO_ste_P')),[1.99 2.01],'Color',[.5 0 1], 'LineWidth', 1);
hold off

ax1 = nexttile;

load coastlines
worldmap('World');
box on
hold on
mmin = -3;
mmax = 3;
l_colb = 60;
key_for_color = linspace(mmin,mmax,l_colb);
fc = brewermap(l_colb,'*RdBu');
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
fc(key_for_color<2&key_for_color>-2,:) = brighten(fc(key_for_color<2&key_for_color>-2,:),.8);
% fc(key_for_color<-0.95&key_for_color>-1.05,:) = [1 0 1];
% fc(key_for_color>0.95&key_for_color<1.05,:) = [1 0 1];
% fc(key_for_color<-1.95&key_for_color>-2.05,:) = [.5 0 1];
% fc(key_for_color>1.95&key_for_color<2.05,:) = [.5 0 1];
number_of_levels=6;
v2 = mmin:(mmax-mmin)/(l_colb/5):mmax;
hold on
contourfm(lat,[lon;360],double(difference_from_MO_ste_T_def'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
colormap(ax1,fc)

geoshow(coastlat,coastlon,'Color','k')
% scatterm(lat_to_stipple_T,lon_to_stipple_T,4,'k','filled')
mlabel off; plabel off; gridm off
hold off
caxis([mmin mmax])
cbmap1=flipud(cbrewer('div','*RdBu',number_of_levels));
hl = colorbar('YTick',[-3 -2 -1 0 1 2 3]);
hl.FontSize = 16;

hold on 
% contourm(lat,[lon;360],abs(double(difference_from_MO_ste_T')),[0.99 1.01],'Color',[1 0 1], 'LineWidth', 1);
% contourm(lat,[lon;360],abs(double(difference_from_MO_ste_T')),[1.99 2.01],'Color',[.5 0 1], 'LineWidth', 1);
hold off
title("(e) Scaling Multi-Obj-1.5, T Error, SE-Normalized")
set(gca,"Fontsize",16)
ax2 = nexttile;

load coastlines
worldmap('World');
box on
hold on
mmin = -3;
mmax = 3;
l_colb = 60;
key_for_color = linspace(mmin,mmax,l_colb);
fc2 = brewermap(l_colb,'BrBG');
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
% fc2(key_for_color<2&key_for_color>-2,:) = 1;
fc2(key_for_color<2&key_for_color>-2,:) = brighten(fc2(key_for_color<2&key_for_color>-2,:),.8);

number_of_levels=6;
v2 = mmin:(mmax-mmin)/(l_colb/5):mmax;
hold on
contourfm(lat,[lon;360],double(difference_from_MO_ste_P_def'),mmin:(mmax-mmin)/l_colb:mmax,'LineColor','none');
colormap(ax2, fc2)

geoshow(coastlat,coastlon,'Color','k')
mlabel off; plabel off; gridm off
hold off
% caxis([mmin mmax])
clim([mmin mmax])
cbmap=flipud(cbrewer('div','BrBG',number_of_levels));
h2 = colorbar('YTick',[-3 -2 -1 0 1 2 3]);
h2.FontSize = 16;
title("(f) Scaling Multi-Obj-1.5, P Error, SE-Normalized")
set(gca,"Fontsize",16)
hold on
% contourm(lat,[lon;360],abs(double(difference_from_MO_ste_P')),[0.99 1.01],'Color',[1 0 1], 'LineWidth', 1);
% contourm(lat,[lon;360],abs(double(difference_from_MO_ste_P')),[1.99 2.01],'Color',[.5 0 1], 'LineWidth', 1);
hold off
set(gcf,'renderer','painters')
print(gcf,'-dpng',["Uncoord_Emu_Paper/Uncoord_Plots/Figure_3_" + getNow() + ".png"],'-r300')


function pattern = pattern_scale(mean_temperatures,all_T_patterns_scaled,T_base_pattern)
pattern = zeros(288,192,length(mean_temperatures(:,1)));
for i = 1:8
    pattern = pattern+all_T_patterns_scaled(:,:,i).*reshape(mean_temperatures(:,i),1,1,[]);      
end
pattern = pattern+T_base_pattern;

end
function T_emulated_all = emulate_all_injection_array(param_AOD_all,param_T_all,injection_all,CO2_forcing)
% [60N; 30N; 15N; 0N; 15S; 30S]
inj_points = 7;
size_injection = size(injection_all);
if size_injection(1) == inj_points
    injection_all = injection_all';
end
T_emulated_all = zeros(length(injection_all(:,1)),inj_points+1);
for i = 1:inj_points+1
    if i<inj_points+1
    T_emulated_all(:,i) = emulate_sp_injection(param_AOD_all(i,:),param_T_all(i,:),injection_all(:,i),0*CO2_forcing);
    else
    T_emulated_all(:,i) = emulate_sp_injection(param_AOD_all(4,:),param_T_all(4,:),0*injection_all(:,4),CO2_forcing);

    end
end
end

function AOD_emulated_all = emulate_all_injection_array_AOD(param_AOD_all,param_T_all,injection_all,CO2_forcing)
% [60N; 30N; 15N; 0N; 15S; 30S]
inj_points = 7;
size_injection = size(injection_all);
if size_injection(1) == inj_points
    injection_all = injection_all';
end
AOD_emulated_all = zeros(length(injection_all(:,1)),inj_points+1);
for i = 1:inj_points+1
    if i<inj_points+1
    AOD_emulated_all(:,i) = emulate_sp_injection_AOD(param_AOD_all(i,:),param_T_all(i,:),injection_all(:,i),0*CO2_forcing);
    else
    AOD_emulated_all(:,i) = emulate_sp_injection_AOD(param_AOD_all(4,:),param_T_all(4,:),0*injection_all(:,4),CO2_forcing);

    end
end
end

function T_emulated = emulate_sp_injection(param_AOD,param_T,injection,CO2_forcing)
beta = param_AOD(1);
alpha = param_AOD(2);
gamma = param_AOD(3);
AOD_emulated = zeros(length(injection),1);
for k = 1:length(AOD_emulated)
    AOD_emulated(k) = 0;
    for j = 1:k
        AOD_emulated(k) = AOD_emulated(k) + impulse_firstOdiff_wGammaInj(beta, alpha,gamma,injection(k-j+1), j)*injection(k-j+1);
    end
end

% for k = 1:length(AOD_emulated)
%     AOD_emulated(k) = 0;
%     for j = 0:k-1
%         AOD_emulated(k) = AOD_emulated(k) + impulse_firstOdiff_wGammaInj(beta, alpha,gamma,injection(k-j), j)*injection(k-j);
%     end
% end

% Load in parameters for GHG
impulse_p_CO2 = [];
impulse_p_CO2.tau = param_T(1);
impulse_p_CO2.mu = param_T(2);
% Load in parameters for SAI
impulse_p_SAI = [];
impulse_p_SAI.tau = param_T(1);
impulse_p_SAI.mu = param_T(3);
T_emulated = zeros(length(injection),1);
% Convolve the impulse response
for k = 1:length(T_emulated)
    for j = 1:k
        T_emulated(k) = T_emulated(k) + impulse_semiInfDiff(j,impulse_p_CO2)*CO2_forcing(k-j+1)+ impulse_semiInfDiff(j,impulse_p_SAI)*AOD_emulated(k-j+1);
    end
end
% for k = 1:length(T_emulated)
%     for j = 0:k-1
%         T_emulated(k) = T_emulated(k) + impulse_semiInfDiff(j,impulse_p_CO2)*CO2_forcing(k-j)+ impulse_semiInfDiff(j,impulse_p_SAI)*AOD_emulated(k-j);
%     end
% end


end


function AOD_emulated = emulate_sp_injection_AOD(param_AOD,param_T,injection,CO2_forcing)
beta = param_AOD(1);
alpha = param_AOD(2);
gamma = param_AOD(3);
AOD_emulated = zeros(length(injection),1);
for k = 1:length(AOD_emulated)
    AOD_emulated(k) = 0;
    for j = 1:k
        AOD_emulated(k) = AOD_emulated(k) + impulse_firstOdiff_wGammaInj(beta, alpha,gamma,injection(k-j+1), j)*injection(k-j+1);
    end
end

% for k = 1:length(AOD_emulated)
%     AOD_emulated(k) = 0;
%     for j = 0:k-1
%         AOD_emulated(k) = AOD_emulated(k) + impulse_firstOdiff_wGammaInj(beta, alpha,gamma,injection(k-j), j)*injection(k-j);
%     end
% end



end

function T_emulated = emulate_sp_injection_old(param_AOD,param_T,injection,CO2_forcing)
beta = param_AOD(1);
alpha = param_AOD(2);
AOD_emulated = zeros(length(injection),1);
for k = 1:length(AOD_emulated)
    AOD_emulated(k) = 0;
    for j = 1:k
        AOD_emulated(k) = AOD_emulated(k) + impulse_firstOdiff(beta, alpha, j)*injection(k-j+1);
    end
end

% Load in parameters for GHG
impulse_p_CO2 = [];
impulse_p_CO2.tau = param_T(1);
impulse_p_CO2.mu = param_T(2);
% Load in parameters for SAI
impulse_p_SAI = [];
impulse_p_SAI.tau = param_T(1);
impulse_p_SAI.mu = param_T(3);
T_emulated = zeros(length(injection),1);
% Convolve the impulse response
for k = 1:length(T_emulated)
    for j = 1:k
        T_emulated(k) = T_emulated(k) + impulse_semiInfDiff(j,impulse_p_CO2)*CO2_forcing(k-j+1)+ impulse_semiInfDiff(j,impulse_p_SAI)*AOD_emulated(k-j+1);
    end
end

% for k = 1:length(T_emulated)
%     for j = 0:k-1
%         T_emulated(k) = T_emulated(k) + impulse_semiInfDiff(j,impulse_p_CO2)*CO2_forcing(k-j)+ impulse_semiInfDiff(j,impulse_p_SAI)*AOD_emulated(k-j);
%     end
% end

end

function AOD_emulated = emulate_sp_injection_AOD_old(param_AOD,param_T,injection,CO2_forcing)
beta = param_AOD(1);
alpha = param_AOD(2);
AOD_emulated = zeros(length(injection),1);
for k = 1:length(AOD_emulated)
    AOD_emulated(k) = 0;
    for j = 1:k
        AOD_emulated(k) = AOD_emulated(k) + impulse_firstOdiff(beta, alpha, j)*injection(k-j+1);
    end
end

% Load in parameters for GHG
impulse_p_CO2 = [];
impulse_p_CO2.tau = param_T(1);
impulse_p_CO2.mu = param_T(2);
% Load in parameters for SAI
impulse_p_SAI = [];
impulse_p_SAI.tau = param_T(1);
impulse_p_SAI.mu = param_T(3);
T_emulated = zeros(length(injection),1);
% Convolve the impulse response
for k = 1:length(T_emulated)
    for j = 1:k
        T_emulated(k) = T_emulated(k) + impulse_semiInfDiff(j,impulse_p_CO2)*CO2_forcing(k-j+1)+ impulse_semiInfDiff(j,impulse_p_SAI)*AOD_emulated(k-j+1);
    end
end

end