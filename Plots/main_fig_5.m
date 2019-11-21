%% Plot change in dynamics and age profile of infection post-vaccination
clear
D_base = load('Baseline_incidence_R1.5.mat');
D_ctrl = load('fifty_reduction_R1.5.mat');
load('predictions_for_plotting_model_schools_R_1.5.mat');
Eff_struct = load('EffectivenessData_R1.5.mat');
X15 = Eff_struct.data;
A_base = load('Baseline_age_distribution_R1.5.mat');
base_age_distrib = A_base.age_distribution;
A_ctrl = load('fifty_reduction_age_distribution_R1.5.mat');
ctrl_age_distrib = A_ctrl.age_distribution;
R = 3.336;
load('vaccine_efficiency_1.5.mat');
load('baseline_inf_numbers_1.5.mat');
%%
I_base = D_base.incidence;
t_base = D_base.times;

I_ctrl = D_ctrl.incidence;
t_ctrl = D_ctrl.times;
t = t_base;
%% Calculate median and 2.5-97.5% prediction intervals
I_perc_base = prctile(I_base,[2.5,50,97.5],1);
I_av_base = mean(I_base);
I_sig_base = std(I_base);

I_perc_ctrl = prctile(I_ctrl,[2.5,50,97.5],1);
I_av_ctrl = mean(I_ctrl);
I_sig_ctrl = std(I_ctrl);
%%
t_ten_years = 0:7:(365*10);
ind = 1:length(t_ten_years);

I_b_05 = prctile(R*I_base,5,1);
I_b_50 = prctile(R*I_base,50,1);
I_b_95 = prctile(R*I_base,95,1);
%% Data processing for retrospective analysis
dates = datestr(pred_times + datenum(2000,1,1));
TotalTrueIncidence = sum(true_incidence,2);
TotalModelIncidence = sum(model_pred,2);

Upper95Pred = poissinv(0.995*ones(size(TotalModelIncidence)),TotalModelIncidence) - TotalModelIncidence;
Lower95Pred = TotalModelIncidence - poissinv(0.005*ones(size(TotalModelIncidence)),TotalModelIncidence);

StandErr = (TotalTrueIncidence - TotalModelIncidence)./sqrt(TotalModelIncidence);
[mu_hat,std_hat] = normfit(StandErr);

%Age distribution of cases
TrueAgeDistribution = sum(true_incidence,1);


ModelAgeDistribution = sum(model_pred,1);




%%
figure(1)
clf
hold on
%%
%Bin the projected age distributions -- mean per year
B2 = R*[sum(base_age_distrib(1:3,:));sum(base_age_distrib(4:6,:));sum(base_age_distrib(7:9,:));sum(base_age_distrib(10:12,:));base_age_distrib(13:end,:)];
Base2 = mean(B2,2);
C2 =  R*[sum(ctrl_age_distrib(1:3,:));sum(ctrl_age_distrib(4:6,:));sum(ctrl_age_distrib(7:9,:));sum(ctrl_age_distrib(10:12,:));ctrl_age_distrib(13:end,:)];
Ctrl2 = mean(C2,2);
%%
%Prediction intervals for future age distributions
pred_intervals_B2 = zeros(length(Base2),2);
pred_intervals_C2 = zeros(length(Ctrl2),2);
for i = 1:length(Ctrl2)
    pred_intervals_B2(i,:) = prediction_interval_for_mixed_poisson(B2(i,:),0.025,0.975);
    pred_intervals_C2(i,:) = prediction_interval_for_mixed_poisson(C2(i,:),0.025,0.975);
end

%%
figure(1)
clf
subplot(1,2,2)
hold on
xticklocs = 1:5:40;

bar(xticklocs-1,Base2,0.375,'FaceColor',[0 0 1],'LineStyle','none')
bar(xticklocs+1,Ctrl2,0.375,'FaceColor',[1 0 0],'LineStyle','none')
errorbar(xticklocs-1,Base2,Base2 -pred_intervals_B2(:,1)  ,pred_intervals_B2(:,2) - Base2,'LineStyle','none','color','black','LineWidth',2)
errorbar(xticklocs+1,Ctrl2,Ctrl2 -pred_intervals_C2(:,1)  ,pred_intervals_C2(:,2) - Ctrl2,'LineStyle','none','color','black','LineWidth',2)


xticks(xticklocs)
xticklabels({'0-3m','3-6m','6-9m','9m-1y','1-2y','2-3y','3-4y','4-5y'});
XtickAgeDistrib2 = {'0-3m','3-6m','6-9m','9m-1y','1-2y','2-3y','3-4y','4-5y'};
l = legend('Baseline','Post-intervention');


l.Box = 'off';
set(gca,'FontSize',16)
ylabel('Hospitalisations','FontSize',22)
title('10 Year forecast: age distribution','FontSize',19)
% ylim([0,75])

subplot(1,2,1)


hold on
shadedErrorBar( t_ten_years ,I_b_50(ind),[I_b_95(ind)-I_b_50(ind) ; I_b_50(ind)-I_b_05(ind)],...
            {'Color',[0 0 1],'LineWidth',3},0.5);
shadedErrorBar( t_ten_years ,0.5*I_b_50(ind),0.5*[I_b_95(ind)-I_b_50(ind) ; I_b_50(ind)-I_b_05(ind)],...
            {'Color',[1 0 0],'LineWidth',3},0.5);

ylim([0,35])
xlim([-20,10*365+20])
set(gca,'FontSize',18,'XTick',365:365:365*10,'XTickLabel',1:10)

xlabel('Years after vaccination implementation','FontSize',18);
ylabel('Weekly hospitalisations','FontSize',22);

title('10 Year forecast: time series','FontSize',19)

% % Legend for projection plot
line([1*365 1.5*365],[33 33],'Color',[0 0 1],'LineWidth',3)
line([1*365 1.5*365],[30 30],'Color',[1 0.3 0.3],'LineWidth',3)
text(1.6*365,33.2,'Baseline','FontSize',18);
text(1.6*365,30.2,'Post-intervention','FontSize',18);


%% 
annotation('textbox',[0.03,0.86,0.1,0.1],'String','A','FitBoxToText','on','FontSize',40,'LineStyle','none');
annotation('textbox',[0.48,0.86,0.1,0.1],'String','B','FitBoxToText','on','FontSize',40,'LineStyle','none');
