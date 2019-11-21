%% Figure 2
% load data
D_base = load('Baseline_incidence_R1.5.mat');
D_ctrl = load('fifty_reduction_R1.5.mat');
load('predictions_for_plotting_model_schools_R_1.5.mat');
Eff_struct = load('EffectivenessData_R1.5.mat');
X15 = Eff_struct.data;
A_base = load('Baseline_age_distribution_R1.5.mat');
base_age_distrib = A_base.age_distribution;
A_ctrl = load('fifty_reduction_age_distribution_R1.5.mat');
ctrl_age_distrib = A_ctrl.age_distribution;
load('underlying_latent_vars_1.5.mat');
load('N1_sim.mat');
%%
I_base = D_base.incidence;
t_base = D_base.times;

I_ctrl = D_ctrl.incidence;
t_ctrl = D_ctrl.times;
t = t_base;
%% Calculate median and 1-99% prediction intervals
I_perc_base = prctile(I_base,[1,50,99],1);
I_av_base = mean(I_base);
I_sig_base = std(I_base);

I_perc_ctrl = prctile(I_ctrl,[1,50,99],1);
I_av_ctrl = mean(I_ctrl);
I_sig_ctrl = std(I_ctrl);

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



%% Find dates for ticks
TimesForTicks = [];
for y = 2002:2016
    TimesForTicks = [TimesForTicks;datenum(y,1,1) - datenum(2000,1,1)+1];
end

DatesForTicks = datevec(TimesForTicks);
DatesForTicks = DatesForTicks(:,1) + 2000;   


%% Main plot

figure(1)
clf

% top plot
subplot(2,1,1)

hold on


scatter(pred_times(2:end),sum(true_incidence,2),70,'filled','MarkerFaceColor','red')
shadedErrorBar(pred_times(2:end),sum(model_pred,2),[Upper95Pred';Lower95Pred'],{'LineWidth',3,'color',[0 0 0]},1)
xlim([pred_times(2)-150,pred_times(end)+150])
ylim([0,35])
set(gca,'XTick',TimesForTicks(1:2:end),'XTickLabel',DatesForTicks(1:2:end));
set(gca,'FontSize',20)
ylabel('Weekly hospitalisations','FontSize',20)
title('Regression on past hospitalisations')
% Legend for top left
line([4384+380 4384+(365+ 380)],[32 32],'LineWidth',3,'color',[0 0 0])
scatter(4384+(365/2)+380,29,70,'filled','MarkerFaceColor','red')
text(4384+(465 + 380),32.2,'Model','FontSize',15);
text(4384+(465 + 380),29.2,'Data','FontSize',15);

yyaxis right
plot(times_N1,N1_sim,'LineWidth',3,'LineStyle','--')
ylabel('Number of under ones','Rotation',270)

 

%% bottom left plot

subplot(2,1,2)
hold on
M2 = [sum(ModelAgeDistribution(1:3)),sum(ModelAgeDistribution(4:6)),sum(ModelAgeDistribution(7:9)),sum(ModelAgeDistribution(10:12)),ModelAgeDistribution(13:end)];
T2 = [sum(TrueAgeDistribution(1:3)),sum(TrueAgeDistribution(4:6)),sum(TrueAgeDistribution(7:9)),sum(TrueAgeDistribution(10:12)),TrueAgeDistribution(13:end)];

%Bin the projected age distributions
B2 = [sum(base_age_distrib(1:3,:));sum(base_age_distrib(4:6,:));sum(base_age_distrib(7:9,:));sum(base_age_distrib(10:12,:));base_age_distrib(13:end,:)];
Base2 = mean(B2,2);
C2 =  [sum(ctrl_age_distrib(1:3,:));sum(ctrl_age_distrib(4:6,:));sum(ctrl_age_distrib(7:9,:));sum(ctrl_age_distrib(10:12,:));ctrl_age_distrib(13:end,:)];
Ctrl2 = mean(C2,2);

%Prediction intervals for future age distributions
pred_intervals_C2 = zeros(length(Ctrl2),2);
for i = 1:length(Ctrl2)
    pred_intervals_C2(i,:) = prediction_interval_for_mixed_poisson(C2(i,:),0.025,0.975);
end


BarUpper95Pred = poissinv(0.9995*ones(size(M2)),M2) - M2;
BarLower95Pred = M2 - poissinv(0.0005*ones(size(M2)),M2);

BarUpper95Pred_base = poissinv(0.975*ones(size(base_age_distrib)),base_age_distrib) - base_age_distrib;
BarLower95Pred_base = base_age_distrib - poissinv(0.025*ones(size(base_age_distrib)),base_age_distrib);

%     

hold on
xticklocs = 1:5:40;



bar(xticklocs-1,T2,0.15,'FaceColor',[1 0 0],'LineStyle','none')
bar(xticklocs,M2,0.15,'FaceColor',[0.7 0.7 1],'LineStyle','none')
errorbar(xticklocs,M2,BarLower95Pred,BarUpper95Pred,'LineStyle','none','color','black','LineWidth',2)



xticks(xticklocs)
xticklabels({'0-3m','3-6m','6-9m','9m-1y','1-2y','2-3y','3-4y','4-5y'});
XtickAgeDistrib2 = {'0-3m','3-6m','6-9m','9m-1y','1-2y','2-3y','3-4y','4-5y'};
l = legend('Data','Model');

l.Box = 'off';
set(gca,'FontSize',18)
ylabel('Hospitalisations','FontSize',22)
title('Hospitalisation age profile','FontSize',19)

%%
annotation('textbox',[0.03,0.86,0.1,0.1],'String','A','FitBoxToText','on','FontSize',40,'LineStyle','none');
annotation('textbox',[0.03,0.4,0.1,0.1],'String','B','FitBoxToText','on','FontSize',40,'LineStyle','none');
