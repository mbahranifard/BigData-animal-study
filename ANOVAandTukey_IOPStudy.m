close all
clear all
cd('C:\Users\mbahranifard3\Documents\Master_2019\Matlab\Stats Analysis\IOPmeasurementTransgenic\');
RawData = readtable ('AvgIOP2.csv');
% figure(1)
% outlier detection and removal
% OutElem = isoutlier(RawData(:,1:3),1);
%for mike it was logarithmic
% OutData = RawData(:,1:3);
% OutData(OutElem) = NaN;
OutData = RawData;

%% SW Test Normality

% for i = 1:3
% [hk,pk] = swtest(OutData(:,i))
% end


%%
% yMean = mean(OutData,1,'omitnan');
% err_down = prctile(OutData,2.5,1);
% err_up = prctile(OutData,97.5,1);
% figure
% hold on
% bar([1,2,3],yMean);
% errorbar([1,2,3],yMean,abs(yMean-err_down),abs(yMean-err_up),'color','k')

%%
%one way anova
% [p,tbl,stats] = anova1(OutData);
% months = num2str(OutData(:,1));
% y = round(OutData(:,2),2)
% y = ceil(OutData(:,2))
y = table2array(OutData(:,2));
months = table2array(OutData(:,1));
Het = (OutData{:,3});
[p,tbl,stats] = anovan(y,{months,Het},'model',2,'varnames',["age","Het"])
figure
[results,means,~,gnames] = multcompare(stats,'Dimension',[1 2],"CType","hsd"); %"hsd"

%% confidence interval on means - from ANOVA
Ameans = stats.means;
n = stats.n;
s = stats.s;
df = numel(OutData(~isnan(OutData)))-numel(Ameans);
CV = tinv (.975,df)*s*sqrt(1./n);
