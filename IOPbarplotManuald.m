
clear all
close all
cd ('C:\Users\mbahranifard3\Documents\Master_2019\Matlab\Stats Analysis\IOPmeasurementTransgenic\')
        raw = readtable ('alltable.csv');
%
    Gsum = groupsummary(raw,["Group","Tg"])
    %change the reference to 100 for visualizaiton
%     raw.Var2 = raw.Var2-100;
    
    figure('Units','pixels','WindowStyle','normal','Position',[200,200,325,250]);
    
    raw = rmmissing(raw, 1);
    twosum = groupsummary(raw,["Group","Tg"])
    [secondarygrp,gpnum] = findgroups (raw.Group, raw.Tg);
        splitknot =splitapply(@mean,raw.IOP,secondarygrp);
        amatknot =splitapply(@std,raw.IOP,secondarygrp);

     grps=unique(secondarygrp);
     
     for i = 1:numel(grps)
         j = twosum.GroupCount(i)
         normrand = randn(1,j);
            jitterrand = 0.3*normrand/max(abs(normrand));
         secondarygrp (secondarygrp == i) = repelem(i,j) + jitterrand;
     end
     
         
    secgrptxt = string(raw.ID);
%     for k = 1:size(grps)

    bp=bar(grps, splitknot,'FaceColor','flat');
    hold all
    errorbar(grps,splitknot,amatknot,'linestyle','none','color','k','linewidth',1.2)
%     end
%    scatter(secondarygrp,raw.AllMat_4,20,'jitter','on','jitterAmount',0.25,'markerfacecolor','k','markeredgecolor','none')
  textscatter(secondarygrp,raw.IOP, secgrptxt)


         
 colormap = [0 0.4470 0.7410;0.8500 0.3250 0.0980	;0.4660 0.6740 0.1880	;0.9290 0.6940 0.1250	]
 for k = 1:size(grps)
    bp.CData(k,:) = colormap(k,:);
 end
 ylim(gca,[0 22])
ylabel('IOP (mmHg)')
set(gca, 'xtick',[1 2 3  4],'xticklabel',[])
        
%     scatter(raw.grouping, raw.Var3,'jitter', 'on','jitterAmount',0.25)
% 