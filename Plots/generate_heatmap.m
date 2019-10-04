function [ HM,heatmap_data0,house_prot,mat_prot ] = generate_heatmap( X0,title_str ,col_lim)
%% Collate the pairs of interventions
mat_prot = unique(X0(1,:));
house_prot = unique(X0(2,:));

%% Find median protection for the intervention combination

heatmap_data0 = zeros(length(mat_prot),length(house_prot));
i = 1;
j=1;
for mp = mat_prot
    for hp = house_prot
        if ~(mp == 0 && hp == 0)
            F = find( (X0(1,:)==mp)&(X0(2,:)==hp));
            d = X0(3:end,F);
            heatmap_data0(i,j) = median(d);
        end
        j = j+1;
    end
    j = 1;
    i = i+1;
end
xlabels = {'0%','25%','50%','75%','100%'};
HM = heatmap(xlabels,mat_prot,heatmap_data0*100);
HM.Colormap =  flipud(jet(1200));
HM.Title = title_str;
HM.FontSize = 24;
HM.ColorLimits = [0,col_lim];
HM.XLabel = 'IRP household coverage';
HM.YLabel = 'Duration of MAB protection (days)';
HM.ColorbarVisible = 'off';

end

