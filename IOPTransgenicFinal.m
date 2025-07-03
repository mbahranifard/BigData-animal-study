close all hidden
clear all
clc

cd('C:\Users\mbahranifard3\Documents\Master_2019\Matlab\Stats Analysis\IOPmeasurementTransgenic\');

% Mouse.xlsx is original figures, Mouse1 was updated on 6.30.2023 with
% iPSC-TMs
AnimalRaw = readtable ('Mouse4.xlsx','Sheet',1,'Format','auto');
IOPRaw = readtable ('Mouse4.xlsx','Sheet',3,'Format','auto','ReadVariableNames',false);
Injected = readtable ('Mouse4.xlsx','Sheet',2,'Format','auto');
%added for iPSCTM
ipsc = readtable ('ipsctm.xlsx','Sheet',1,'Format','auto','ReadVariableNames',false);




% Gadjet = NaN(500,4);
% GroupVec = [2.5,3.5,4.5,5.5,6];
GroupVec = [4.5,6];
IOPnum = transpose(IOPRaw{1,:});
IOPnum = rmmissing(IOPnum);
AnimalRmInd = find(contains(IOPnum,'Ani'));
IOPnum(AnimalRmInd) = [];
IOPmat = str2num(cell2mat(IOPnum));

HetWT = find(contains(AnimalRaw.Het_WT,'W'));
HetWTnum = ones(length(AnimalRaw.AnimalID),1);
HetWTnum(HetWT) = 0;

for i= 1:length(ShamS)
    x = find(AnimalRaw.AnimalID==ShamS(i));
    HetWTnum(x)= 2;
end


for i= 1:length(ShamL)
    x = find(AnimalRaw.AnimalID==ShamL(i));
    HetWTnum(x)= 3;
end



rmtemp = [];






%     for j=1:length(IOPnum)
%     IOPmat(j) = str2num(IOPnum{j});
%     end
% 
% for col = 12:size(AnimalRaw,2)

    
for col = 12:28
    if numel(rmmissing(AnimalRaw{:,col}))==0
        continue
    end
    Gadjet(1:numel(AnimalRaw.AnimalID),1) = AnimalRaw.AnimalID;
    Gadjet(1:numel(AnimalRaw.AnimalID),2) = hours((AnimalRaw{:,col}-AnimalRaw.DOB))/(24*30);
    Gadjet(1:numel(AnimalRaw.AnimalID),5) = HetWTnum;
    
    GadjetDates(1:numel(AnimalRaw.AnimalID)) = AnimalRaw{:,col};

    GadjetDates(isnan(Gadjet(:,2))) = []; 
    Gadjet(isnan(Gadjet(:,2)),:) = []; 
    
    
        for i = 1:numel(Gadjet(:,1))
            if Gadjet(i,2) > 6 
                Gadjet(i,3) = GroupVec(end);
            elseif (Gadjet(i,2)-GroupVec(1)) < 0.6 & (Gadjet(i,2)-GroupVec(1)) >= 0
                Gadjet(i,3) = GroupVec(1);
            else 
                rmtemp = [rmtemp i]; 
                continue

%                 [~,I] = min(abs(GroupVec-Gadjet(i,2)));
%                 Gadjet(i,3)= GroupVec(I);
                
            end
            
            % removing the ones injected
            injrmidx = find(Injected.AnimalID==Gadjet(i,1));
            if ~isempty(injrmidx)
               injrm = hours(GadjetDates(i)-Injected.InjectionDate(injrmidx));
               if injrm >0
               rmtemp = [rmtemp i]; 
                continue
               end
               
            end
            
            %getting the IOP
            K = find(IOPmat == Gadjet(i,1));
            
            if ~isempty(K)
% %             Gadjet(i,4) = mean(IOPRaw{(col-11)*14,[2*K-1,2*K]});
%                 DateFind = find(IOPRaw{:,2*K(1)-1} == char(GadjetDates(i)));
%                     GadDateChar = char(datestr(GadjetDates(i),'mm/dd/yyyy'));
                    GadDateStr = string(datestr(GadjetDates(i),'mm/dd/yyyy'));

%                     if GadDateChar(1)== '0'
%                         GadDateStr = string(GadDateChar(2:end));
%                     end
                    
                    IOPDATE = datetime(IOPRaw{:,2*K(1)-1},'InputFormat','MM/dd/yyyy');
                    IOPDATE = string(datestr(IOPDATE,'mm/dd/yyyy'));
                    DateFind = find(contains(IOPDATE,GadDateStr));


                if ~isempty(DateFind)
                                        
                    pone = str2num(cell2mat(IOPRaw{DateFind+12,2*K(1)-1}));
                    ptwo = str2num(cell2mat(IOPRaw{DateFind+12,2*K(1)}));
                    stdone = str2num(cell2mat(IOPRaw{DateFind+13,2*K(1)-1}));
                    stdtwo = str2num(cell2mat(IOPRaw{DateFind+13,2*K(1)}));
                    Gadjet(i,4) = mean([pone, ptwo]);
                    Gadjet(i,6) = mean([stdone, stdtwo]);

                else 
                 formatSpec = 'Date %s for animal %d missing\n';
                 %remove entire row
                 rmtemp = [rmtemp i]; 

                 fprintf(formatSpec,datestr(GadjetDates(i),'mm/dd/yyyy'),Gadjet(i,1))
                end
            
            end
        
        end
        Gadjet(rmtemp,:)=[];
     AvgIOP{col-11} = Gadjet;
     rmtemp = [];
     clear Gadjet K GadjetDates 
end

AllMat = [AvgIOP{1}];
for i = 2:size(AvgIOP,2)
AllMat = [AllMat;AvgIOP{i}];
end

AllMat(isnan(AllMat(:,1)),:)=[];

%% find non-unique items and remove the samller ages
idxtot = [];
idxtotlong = [];
for t = 1:numel(GroupVec)
    idxhalftot = [];
    maxlongtot = [];
%     maxloc2 =[];
    agematidx = find(AllMat(:,3)==GroupVec(t));
    agemat = AllMat(agematidx,:);
    u = unique(agemat(:,1));
    [n,bin]=histc(agemat(:,1),u);
    v = find(n>1);
    for s = 1: numel(v)
        maxloc2=[];
        idx=find(bin==v(s));
        [~,maxidx]= max(agemat(idx,2));
        % finding short term measurements from long term measurements
        maxloc = idx(maxidx);
        idx(maxidx) = [];
        longidx = find(ShamL == agemat(maxidx,1));
%         agemat(maxidx,5)
        if agemat(maxloc,5)==3
                    [~,maxidx2]= max(agemat(idx,2));
                    maxloc2 = idx(maxidx2);
                    idx(maxidx2) = [];        
        end
        maxlongtot = [maxlongtot; maxloc2]; %says which ones are short terms from long terms
        idxhalftot = [idxhalftot; idx];

    end 
            idxtot = [idxtot ;agematidx(idxhalftot)];
            idxtotlong = [idxtotlong;agematidx(maxlongtot)];
            
end
            AllMat(idxtotlong,5)=2; %switches shorter term reps to short term
            AllMat(idxtot,:)=[];
            
            AllMatSupp=[];
            
            %% iPSC-TM resutls
           
            ipscgp = unique(ipsc.Var6);
            for j=1:numel(ipscgp)
                 ioptempmat =[NaN NaN NaN];
                junktab = ipsc(string(ipsc.Var6)==string(ipscgp(j)),:);
 for i = 1:numel(junktab.Var1)
        K = find(IOPmat == junktab.Var1(i));

if ~isempty(K)
% %             Gadjet(i,4) = mean(IOPRaw{(col-11)*14,[2*K-1,2*K]});
%                 DateFind = find(IOPRaw{:,2*K(1)-1} == char(GadjetDates(i)));
%                     GadDateChar = char(datestr(GadjetDates(i),'mm/dd/yyyy'));
    if ~isnat(junktab.Var4(i))
    GadDateStr = string(datestr(junktab.Var4(i),'mm/dd/yyyy'));
    else
         formatSpec = 'Date %s for animal %d is not time\n';
        fprintf(formatSpec,datestr(GadDateStr,'mm/dd/yyyy'),junktab.Var1(i))
        continue
    end
%                     if GadDateChar(1)== '0'
%                         GadDateStr = string(GadDateChar(2:end));
%                     end

    IOPDATE = datetime(IOPRaw{:,2*K(1)-1},'InputFormat','MM/dd/yyyy');
    IOPDATE = string(datestr(IOPDATE,'mm/dd/yyyy'));
    DateFind = find(contains(IOPDATE,GadDateStr));


if ~isempty(DateFind)

    pone = str2num(cell2mat(IOPRaw{DateFind+12,2*K(1)-1}));
    if isempty(pone)
        pone=NaN;
    end
    ptwo = str2num(cell2mat(IOPRaw{DateFind+12,2*K(1)}));
    stdone = str2num(cell2mat(IOPRaw{DateFind+13,2*K(1)-1}));
    stdtwo = str2num(cell2mat(IOPRaw{DateFind+13,2*K(1)}));
%     Gadjet(i,4) = mean([pone, ptwo]);
%     Gadjet(i,6) = mean([stdone, stdtwo]);
    if junktab.Var2{i} == "OS"
        ipscP = pone;
        ipscSTD = stdone;
    else
        ipscP = ptwo;
        ipscSTD = stdtwo; 
    end

else 
 formatSpec = 'Date %s for animal %d missing\n';
 %remove entire row
%  rmtemp = [rmtemp i]; 
 fprintf(formatSpec,datestr(GadDateStr,'mm/dd/yyyy'),junktab.Var1(i))
 continue
end
idvec = [6 7 8 9];
%outlier removal
if j==2 & junktab.Var1(i)==269
    continue
end

AllMat=[AllMat;[junktab.Var1(i) 0 6 ipscP idvec(j) ipscSTD]];
AllMatSupp = [AllMatSupp;[num2str(junktab.Var1(i)),junktab.Var2{i}]];

ioptempmat = [ioptempmat;[junktab.Var1(i) pone ptwo]];

end

end
            end
            
           
            %%

            
            %% removing certain measurements
%             AllMat(AllMat(:,1)<149,:)=[];
            %%
            
            SaveName = ('C:\Users\mbahranifard3\Documents\Master_2019\Matlab\Stats Analysis\IOPmeasurementTransgenic\alltable.csv');
            varnames = {'ID','Age','Group','IOP','Tg','std'};
            alltable = array2table(AllMat(:,1:6),'VariableNames',varnames);
            writetable(alltable,SaveName);



% for injected animals loaded from totomat
cellomat=readmatrix('C:\Users\mbahranifard3\Desktop\InjectRaw\totomat.csv');
cellomatidx = find(isnan(cellomat(:)));
cellomatidx(cellomatidx>size(cellomat,1))=[];
cellomat(cellomatidx,:)=[];
cellocell=readtable('C:\Users\mbahranifard3\Desktop\InjectRaw\cellocell.csv');
for i = 1:size(cellomat,2)/3
%     grpidx = find(cellomat(:,4)==groups(i));
%     OutElem = isoutlier(DataRaw.Var3(grpidx),1);
%     grpidx(grpidx == OutElem) = [];
%     grpidxcell = {grpidxcell grpidx};
    Yi{i}  = cellomat(:,3*i-1);
    NameTag{i}=cellocell.cellocell;
%     [H, pValue, SWstatistic]=swtest(Yi);
%     pValue
    sYi{i} = cellomat(:,3*i); %std of measurement
    rmindx = find(isnan(Yi{i}));
    Yi{i}(rmindx)=[];
    sYi{i}(rmindx)=[];
    NameTag{i}(rmindx)=[];
    
    


end

AllMat(AllMat(:,3)==4.5,:)=[];
groups = unique(AllMat(:,5));
grpidxcell = [];
%for non-injected
for i = size(cellomat,2)/3+1:size(cellomat,2)/3+numel(groups)
    grpidx = find(AllMat(:,5)==groups(i-size(cellomat,2)/3));
%     OutElem = isoutlier(DataRaw.Var3(grpidx),1);
%     grpidx(grpidx == OutElem) = [];
    grpidxcell = {grpidxcell grpidx};
    Yi{i}  = AllMat(grpidx,4);
    NameTag{i}=string(AllMat(grpidx,1));
%     [H, pValue, SWstatistic]=swtest(Yi);
%     pValue
    sYi{i} = AllMat(grpidx,6); %std of measurement
    rmindx = find(isnan(Yi{i}));
    Yi{i}(rmindx)=[];
    sYi{i}(rmindx)=[];
    NameTag{i}(rmindx)=[];
    
    


end

f3=figure(100);set(gca,'fontsize',12)
    ha=tight_subplot(1,1,[.02 .125],[.14 .06],[.15 0.05]);
    set(f3,'position',[100 100 500 500],'color','w')
%     title(Names)
    set(gca,'fontsize',11) 

    %     Acol={[256 0 0]/256 [0 0 256]/256 [0 192 0]/256 [255 160 64]/256 };                         % Light versions of red and bluw 
% colcell = lines(10);
% Acol = colcell([1:6,9,10],:);
% Mcol={[256 0 0]/256*.6 [0 0 256]/256*.6 [0 192 0]/256*.6 [255 160 64]/256*.6};  
Mcol = Acol*.3;

Acol = num2cell(Acol,2);
Mcol = num2cell(Mcol,2);
    
%% re-sorting the groups
tempcell={Ybar,MEY,sPop2,Yi,sYi,NameTag};
k = [3,4,5,1,6,2,10,9,7,8];




Captions=["WT","Tg","Sham Short","MSC Short","Sham Mid","MSC Mid","Sham Long","MSC Long","iPSC Short","iPSC Mid"];



   %% stats
%        [p,tbl,stats] = anovan(Yi,{alltable.Group,alltable.Tg},'model',2,'varnames',{'Group','Tg'})
y=[];
vec=[];
for i = 1:numel(Yi)
    y=[y;Yi{i}];
    vec = [vec;transpose(repelem(i,numel(Yi{i})))];
end

       [p,tbl,stats] = anovan(y,vec,'model',1,'varnames',{'Group'})

        statfig=figure(590);
        set(statfig,'Position',[100 100 800 800])
        
       [results,means,~,gnames] = multcompare(stats,"CType","hsd"); %"hsd"
       
       %% multiple comparison without correction. Total of 15 comparisons are being made so critical alpha = 0.05/15
        unique_numbers =unique(vec); 
       counts = histc(vec, unique_numbers);
        noposthocMat=[];

        for i = 1:numel(Ybar)
           for j = i+1:numel(Ybar)
               % Calculate pooled standard deviation
                pooled_std = sqrt(((counts(i) - 1) * sPop2(i) + (counts(j) - 1) * sPop2(j)) / (counts(i) + counts(j) - 2));
                % Calculate t-statistic
                t_statistic = (Ybar(i) - Ybar(j)) / (pooled_std * sqrt(1/counts(i) + 1/counts(j)));
                % Degrees of freedom
                df = counts(i) + counts(j) - 2;
                % Two-tailed t-test
                alpha = 0.05; % Significance level
                p_value = 1-tcdf(abs(t_statistic), df); % Calculate p-value for two-tailed test
                %confidence bounds
%                 t_score = tinv(1 - alpha / 2, df);
                % Calculate the margin of error
%                 margin_of_error = t_score * (pooled_std * sqrt(1/numdif(i) + 1/numdif(j)));
%                 meanci = [abs(avgdif(i) - avgdif(j))];
                noposthocMat = [noposthocMat;[i j p_value]];
           end
       end
       %%
       
%        [results,means,~,gnames] = multcompare(stats,"CType","bonferroni",'alpha',alpha); 
       
       
             % stats for the difference between treatment vs control over time
       %for MSCs
       alpha = 0.05;
       injct = [3,5,7]; injmsc = [4,6,8]; injipsc = [9,10];
       avgdif = [means(injmsc),means(injipsc)]-[means(injct),means(injct(1:2))];
       stddif = sqrt([sPop2(injct),sPop2(injct(1:2))]+[sPop2(injmsc),sPop2(injipsc)]);
       unique_numbers = unique(vec);
       counts = histc(vec, unique_numbers);
       numdif = round(([counts(injct);counts(injct(1:2))]+[counts(injmsc);counts(injipsc)])/2);
       MEYmscratio = sqrt((stddif'.^2)./numdif).*tinv(0.975,numdif-1); 
       ratiostats = [[injmsc,injipsc]; [injct,injct(1:2)]; ([avgdif;avgdif-MEYmscratio';avgdif+MEYmscratio'])];
       diffcompMat=[];
       for i = 1:numel(avgdif)
           for j = i+1:numel(avgdif)
               % Calculate pooled standard deviation
                pooled_std = sqrt(((numdif(i) - 1) * stddif(i)^2 + (numdif(j) - 1) * stddif(j)^2) / (numdif(i) + numdif(j) - 2));
                % Calculate t-statistic
                t_statistic = (avgdif(i) - avgdif(j)) / (pooled_std * sqrt(1/numdif(i) + 1/numdif(j)));
                % Degrees of freedom
                df = numdif(i) + numdif(j) - 2;
                % Two-tailed t-test
                alpha = 0.05; % Significance level
                p_value = 1-tcdf(abs(t_statistic), df); % Calculate p-value for two-tailed test
                %confidence bounds
                t_score = tinv(1 - alpha / 2, df);
                % Calculate the margin of error
                margin_of_error = t_score * (pooled_std * sqrt(1/numdif(i) + 1/numdif(j)));
                meanci = [abs(avgdif(i) - avgdif(j)) margin_of_error];
                diffcompMat = [diffcompMat;[i j meanci p_value]];
           end
       end
      %for ipsc
%       avgdifipsc = means(injipsc)-means(injct(1:2));
%       stddif = sqrt(sPop2(injct(1:2))+sPop2(injipsc));
%       numdif = round((counts(injct(1:2))+counts(injipsc))/2);
%       i=1;j=2;
%     pooled_std = sqrt(((numdif(i) - 1) * stddif(i)^2 + (numdif(j) - 1) * stddif(j)^2) / (numdif(i) + numdif(j) - 2));
%     % Calculate t-statistic
%     t_statistic = (avgdifipsc(i) - avgdifipsc(j)) / (pooled_std * sqrt(1/numdif(i) + 1/numdif(j)));
%     % Degrees of freedom
%     df = numdif(i) + numdif(j) - 2;
%     % Two-tailed t-test
%     alpha = 0.05; % Significance level
%     p_value = 1-tcdf(abs(t_statistic), df); % Calculate p-value for two-tailed test
%     %confidence bounds
%     t_score = tinv(1 - alpha / 2, df);
%     % Calculate the margin of error
%     margin_of_error = t_score * (pooled_std * sqrt(1/numdif(i) + 1/numdif(j)));
%     meanci = [abs(avgdifipsc(i) - avgdifipsc(j)) margin_of_error];
%     MEYipscratio = sqrt((stddif'.^2)./numdif).*tinv(0.975,numdif-1); 
%     ratiostats = [ratiostats,[injipsc; injct(1:2); ([avgdifipsc;avgdifipsc-MEYipscratio';avgdifipsc+MEYipscratio'])]];
%     diffcompMat = [diffcompMat;[9 10 meanci p_value]];
 %%          
 for j = 1:numel(Ybar)
YCI(j,1:3) = ([Ybar(j) Ybar(j)-MEY(j) Ybar(j)+MEY(j)]);
end
round(YCI,1)           
            