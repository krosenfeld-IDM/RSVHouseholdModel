%% Force of infection plot

load('force_of_infection_1.5.mat');
load('mean_inf_rate.mat');
load('In_household_infection_rates_by_HHSize.mat');

%%
figure(1)
clf


subplot(2,2,1)
U1_foi = sum(mean_foi_A(1:12));
School_age_foi = sum(mean_foi_A(16:28));
other_foi = sum(mean_foi_A([13:15 29:30]));
xticklocs = [1,4,7];
labels = {'Under-ones', 'School age', 'Other over-ones' };
bar(xticklocs,[U1_foi,School_age_foi,other_foi]);
set(gca,'FontSize',16,'XTick',xticklocs,'XTickLabel',labels,'XTickLabelRotation',-45);
title('Force of inf. on age groups')%,'FontSize',14)


subplot(2,2,2)
stacked_HH_bars = [mean_foi_H_U1 mean_foi_H_O1_1+mean_foi_H_O1_0];
bar(stacked_HH_bars,'stacked');
set(gca,'FontSize',16)
title('Force of inf. on households')
xlabel('HH size')
legend('U1s','O1s')



subplot(2,2,3)
Meshendpts = [linspace(30.4,365.25,12), linspace(2*365.25,365.25*17,16),18*365.25,150*365.25];
Meshstartpts = [0,Meshendpts(1:(end-1))];

Mesh = mean( [Meshstartpts ; Meshendpts], 1)/365.25;

hold on

% plot(Mesh(1:(end-1)),mean_transmission_rate_in_HH(1:(end-1)) ) ;
 plot([1:30,40],[mean_transmission_rate_in_HH;mean_transmission_rate_in_HH(end)],'color',[ 0 0 0],'LineWidth',7 ) ;

% bar(33,mean_transmission_rate_in_HH(end),5)
set(gca,'XTick',[1:4:29,33],'Box','on')
set(gca,'XTickLabel',{'0-1m','4-5m','8-9m','1yr','5yr','9yr','13yr','17yr','18+yrs'})
set(gca,'FontSize',16,'XTickLabelRotation',-45);
title('Per-capita within HH infection rate')
xlim([0,35])


subplot(2,2,4)
bar([mean(InHH_inf_rate_u1,2) mean(InHH_inf_rate_o1,2)],'stacked')
set(gca,'FontSize',16)
title('Total infection rate within HHs')
xlabel('HH size')
legend('U1s','O1s')
%%

annotation('textbox',[0.04,0.87,0.1,0.1],'String','A','FitBoxToText','on','FontSize',45,'LineStyle','none');
annotation('textbox',[0.48,0.87,0.1,0.1],'String','B','FitBoxToText','on','FontSize',45,'LineStyle','none');
annotation('textbox',[0.04,0.4,0.1,0.1],'String','C','FitBoxToText','on','FontSize',45,'LineStyle','none');
annotation('textbox',[0.48,0.4,0.1,0.1],'String','D','FitBoxToText','on','FontSize',45,'LineStyle','none');
