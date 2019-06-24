%% alternate plot 3
D_base = load('Baseline_incidence_R1.5_finer_DT.mat');
D_ctrl = load('fifty_reduction_R1.5_finer_DT.mat');
load('predictions_for_plotting_model_schools_R_1.5.mat');
Eff_struct = load('EffectivenessData_R1.5.mat');
X15 = Eff_struct.data;
A_base = load('Baseline_age_distribution_R1.5.mat');
base_age_distrib = A_base.age_distribution;
A_ctrl = load('fifty_reduction_age_distribution_R1.5.mat');
ctrl_age_distrib = A_ctrl.age_distribution;
load('underlying_latent_vars_1.5.mat');
% load('Number_U1s_data.mat');
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
% TrueAgeDistribution = sum(true_incidence(end-5*52:end,:),1);
% TrueAgeDistribution = TrueAgeDistribution/sum(TrueAgeDistribution);

ModelAgeDistribution = sum(model_pred,1);
% ModelAgeDistribution = sum(model_pred(end-5*52:end,:),1);
% ModelAgeDistribution = ModelAgeDistribution/sum(ModelAgeDistribution);



%% Find dates for ticks
TimesForTicks = [];
for y = 2002:2016
    TimesForTicks = [TimesForTicks;datenum(y,1,1) - datenum(2000,1,1)+1];
end

% DatesForTicks = datestr(TimesForTicks);
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
% xlabel('Year','FontSize',16);
title('Regression on past hospitalisations')
% Legend for top left
line([4384+380 4384+(365+ 380)],[32 32],'LineWidth',3,'color',[0 0 0])
scatter(4384+(365/2)+380,29,70,'filled','MarkerFaceColor','red')
text(4384+(465 + 380),32.2,'Model','FontSize',15);
text(4384+(465 + 380),29.2,'Data','FontSize',15);

yyaxis right
% plot(days,beta_t,times_N1,N1_sim)
plot(times_N1,N1_sim,'LineWidth',3,'LineStyle','--')
ylabel('Number of under ones','Rotation',270)
% Legend for projection plot
% line([1*365 1.5*365],[32 32],'Color',[0 0 1 0.2],'LineWidth',2,'Parent',ax2)
% line([1*365 1.5*365],[29 29],'Color',[1 0.3 0.3 0.2],'LineWidth',2,'Parent',ax2)
% text(1.6*365,32.2,'Projection: No intervention','FontSize',15);
% text(1.6*365,29.2,'Projection: Post-intervention','FontSize',15);

% seasonality bottom left plot
% subplot(2,1,2)
% plot(days,beta_t)
% set(gca,'XTick',TimesForTicks(1:2:end),'XTickLabel',DatesForTicks(1:2:end));
% set(gca,'FontSize',20)
% 
% ax2 = gca; % current axes
% % ax1.Position = ax1.Position + [-0.05,0,-0.3,0]; % position of first axes
% ax2.Position = ax2.Position + [0,0,-0.3,0]; % position of first axes
% xlim([pred_times(2)-150,pred_times(end)+150])
% 

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

% bar(xticklocs-1,T2,0.15,'FaceColor',[1 0 0],'LineStyle','none')
% bar(xticklocs,M2,0.15,'FaceColor',[0.7 0.7 1],'LineStyle','none')
% bar(xticklocs+1,Ctrl2,0.15,'FaceColor',[0 1 0],'LineStyle','none')
% errorbar(xticklocs,M2,BarLower95Pred,BarUpper95Pred,'LineStyle','none','color','black','LineWidth',2)

bar(xticklocs-1,100*T2/sum(T2),0.15,'FaceColor',[1 0 0],'LineStyle','none')
bar(xticklocs,100*M2/sum(M2),0.15,'FaceColor',[0.7 0.7 1],'LineStyle','none')
% bar(xticklocs+1,100*Ctrl2/sum(Ctrl2),0.15,'FaceColor',[0 1 0],'LineStyle','none')
errorbar(xticklocs,100*M2/sum(M2),100*BarLower95Pred/sum(M2),100*BarUpper95Pred/sum(M2),'LineStyle','none','color','black','LineWidth',2)
% errorbar(xticklocs+1,100*Ctrl2/sum(Ctrl2),100*(Ctrl2 - pred_intervals_C2(:,1))/sum(Ctrl2),100*(pred_intervals_C2(:,2)-Ctrl2)/sum(Ctrl2),'LineStyle','none','color','black','LineWidth',2)

% bar(xticklocs-1,T2,0.15,'FaceColor',[1 0 0],'LineStyle','none')
% bar(xticklocs,M2,0.15,'FaceColor',[0.7 0.7 1],'LineStyle','none')
% errorbar(xticklocs,M2,BarLower95Pred,BarUpper95Pred,'LineStyle','none','color','black','LineWidth',2)


xticks(xticklocs)
xticklabels({'0-3m','3-6m','6-9m','9m-1y','1-2y','2-3y','3-4y','4-5y'});
XtickAgeDistrib2 = {'0-3m','3-6m','6-9m','9m-1y','1-2y','2-3y','3-4y','4-5y'};
% l = legend('Data: 2012-2016','Model: 5 yrs Pre-intervention','Model: 5 yrs Post-intervention');
l = legend('Data','Model');

l.Box = 'off';
set(gca,'FontSize',18)
ylabel('% hospitalisations','FontSize',22)
title('Hospitalisation age profile','FontSize',19)
% ylim([0,75])
% bottom right

%%
annotation('textbox',[0.03,0.86,0.1,0.1],'String','A','FitBoxToText','on','FontSize',40,'LineStyle','none');
annotation('textbox',[0.03,0.4,0.1,0.1],'String','B','FitBoxToText','on','FontSize',40,'LineStyle','none');
