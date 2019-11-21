%% Script for plotting heat maps of the mean effectiveness of interventions

load('Eff_inf_Mat_R1.5_MABcov_0.5.mat');
load('Eff_inf_Mat_R1.5_MABcov_1.0.mat');
load('EffMat_R1.5_MABcov_0.5.mat');
load('EffMat_R1.5_MABcov_1.0.mat');
%% Plot
figure(1)
clf


subplot(2,2,1)

HM = generate_heatmap(X15_1,'% reduction in hospitalisation',60);
HM.FontSize = 18;

title('% reduction in hospitalisation');



subplot(2,2,3)

HM = generate_heatmap(IX15_1,'% reduction in total infections',10);
HM.FontSize = 18;

annotation('textbox',[0.03,0.89,0.1,0.1],'String','A','FitBoxToText','on','FontSize',40,'LineStyle','none');
annotation('textbox',[0.48,0.89,0.1,0.1],'String','B','FitBoxToText','on','FontSize',40,'LineStyle','none');

subplot(2,2,2)

HM = generate_heatmap(X15_50,'% reduction in hospitalisation',60);
HM.FontSize = 18;

title('% reduction in hospitalisation');



subplot(2,2,4)

HM = generate_heatmap(IX15_50,'% reduction in total infections',10);
HM.FontSize = 18;

annotation('textbox',[0.03,0.42,0.1,0.1],'String','C','FitBoxToText','on','FontSize',40,'LineStyle','none');
annotation('textbox',[0.48,0.42,0.1,0.1],'String','D','FitBoxToText','on','FontSize',40,'LineStyle','none');


annotation('textbox',[0.15,0.9,0.1,0.1],'String','100% maternal vaccine coverage','FitBoxToText','on','FontSize',30,'LineStyle','none');
annotation('textbox',[0.65,0.9,0.1,0.1],'String','50% maternal vaccine coverage','FitBoxToText','on','FontSize',30,'LineStyle','none');


