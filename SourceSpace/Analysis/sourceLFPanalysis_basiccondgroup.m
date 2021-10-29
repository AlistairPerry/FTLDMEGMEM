function [MMNpall] = sourceLFPanalysis_basiccondgroup

%% Setup
%Open up SPM dir (and field trip)

%restoredefaultpath


addpath('/imaging/rowe/users/ap09/Toolbox')

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest')

addpath('/imaging/rowe/users/ap09/Toolbox/fieldtrip')

%addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/fieldtrip/'); %EK

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/boundedline-pkg'))

addpath('/imaging/rowe/users/ap09/Toolbox/ekMEG/ntadscripts')

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/VBA-toolbox'))


ft_defaults()


%% Subjects

preprocdir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/analysis/SourceSpace/LFPs_icafixes/LFPs_COH_wrad_allMEGch_wfids/B_Data';

preprocsteps = 'fmraeMaffffdS';

LFPprefix = 'LFP6inv1';

LFPpostfix = 'rov123_sss.mat'; %for now


subjs = {'c2','c3','c7','c8','c9','c11','c13','c14','c18','c20','c23','c24','c1','c4','c5','c6','c10','c12','c15','c16','c17','c19','c21','c22', 'p1','p3','p6','p7','p9','p12','p13','p14','p17','p18','p21','p24','p2','p4','p5','p8','p10','p11','p15','p16','p19','p20','p22','p23','p25'};

sessions = {'s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2'};


%Final subjects - those w/both sessions

subjs_new = {'c2';'c3';'c7';'c8';'c9';'c11';'c13';'c14';'c23';'c24';'c1';'c4';'c5';'c6';'c10';'c12';'c15';'c16';'c21';'c22'; 'p1';'p6';'p7';'p9';'p12';'p13';'p14';'p17';'p18';'p21';'p24';'p2';'p4';'p5';'p8';'p11';'p16';'p22';'p23';'p25'};


%Now find their sessions

for i = 1:length(subjs_new)
    
    match = strcmp(subjs, subjs_new{i});
    
    
    if nnz(match)
        
        sessions_new(i,1) = sessions(match);
        
    end
    
end


sessions = sessions_new;

subjs = subjs_new;



%Pull out cons and pats

for i = 1:length(subjs)
    
   
    cmatch = contains(subjs{i},'c');
    
    
    if cmatch == 1
    
        grpst(i,1) = 1;
       
        
    else
        
        grpst(i,1) = -1;
        
        
    end
    
    
end


subjs_con = subjs(grpst==1);

sessions_con = sessions(grpst==1);


subjs_pat = subjs(grpst==-1);

sessions_pat = sessions(grpst==-1);



%% Other setup


devtrl = 1;
std1trl = 2;
std3trl = 4;
stdtrl = 7;


ds = 1;


if ds == 1
    
    %Determine MMN start/finish
    
    MMNs = 226./2;
    MMNf = 276./2;
    
    xticks=[1 50 100 150 200 251];
    xlabels={'-100' '0' '100' '200' '300' '400'};
    
else
    
    MMNs = 226;
    MMNf = 276;
    
    xticks=[1 100 200 300 400];
    xlabels={'-100' '0' '100' '200' '300' '400'};
    
end


nsources = 6;


%Figure output dir

FigOutDir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/analysis/SourceSpace/LFPs_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots';

mkdir(FigOutDir);


figoutprefix = 'commconpatplac_remo_wbothsess';


%Plus disease group information

% Can call it directly as it will never change

DiagData = readcell('/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/PatDiagInfo.txt');



%% Read in placebo and drug data - controls


%Load data and setup

% fid = fopen(LFPplaclistcon);
% Data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');

nfilesc=length(subjs_con);


%Create subject cells

avgdev_all=cell(nfilesc,1);

avgstd_all=cell(nfilesc,1);


% Read in patient data


% fid = fopen(LFPplaclistpat);
% Data_Pat = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');

nfilesp=length(subjs_pat);



%Create subject cells

avgdev_all_pat=cell(nfilesp,1);

avgstd_all_pat=cell(nfilesp,1);



%% Setup design matrix(s)

for s = 1:nfilesc
    
S(s,1) = s;
    
end

F1 = repmat([1;1], [nfilesc 1]);

F2 = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesc 1]));


%Now patients

for s = 1:nfilesp
    
        Sp(s,1) = s+nfilesc;
    
end

F1p = repmat([2;2], [nfilesp 1]);
F2p = cat(1, repmat(1, [nfilesp 1]), repmat(2, [nfilesp 1]));


%Combine design matrix(s)

Sall = cat(1, S, S, Sp, Sp);

F1all = cat(1, F1, F1p);
F2all = cat(1, F2, F2p);



%% Loop through subjects/files for controls
%Placebo


for subj = 1:nfilesc
    
    
    D = spm_eeg_load([preprocdir '/' subjs_con{subj} '/' LFPprefix '_' preprocsteps '_' subjs_con{subj} '_' sessions_con{subj} '_' LFPpostfix]);
    
    
    %Other processing, including baseline correction and smoothing
    
    timelock=D.fttimelock;
    
    %Baseline correct
    cfg = [];
    cfg.baseline = [-0.100 -0.001];
    [baselinecorr] = ft_timelockbaseline(cfg, timelock);

    
    %Go back into sources
    
    %Std
    cfg = [];
    cfg.trials = timelock.trialinfo == 7;
    data_std = ft_redefinetrial(cfg, baselinecorr);
    
    
    %Std3
    cfg = [];
    cfg.trials = timelock.trialinfo == 4;
    data_std3 = ft_redefinetrial(cfg, baselinecorr);
    
    %Std1
    cfg = [];
    cfg.trials = timelock.trialinfo == 2;
    data_std1 = ft_redefinetrial(cfg, baselinecorr);
    
    %Dev
    cfg = [];
    cfg.trials = timelock.trialinfo == 1;
    data_dev = ft_redefinetrial(cfg, baselinecorr);
    
 
    
    
    for source = 1:nsources
        
        %Std
        
        data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
        
        
        %Std3
        
        data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
        
        
        %Std1
        
        data_std1.avg(source,:) = smooth(data_std1.avg(source,:), 20, 'moving');
        
        
        %Dev
        
        data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
        
    end
    
    
    %Combine into a group cell
    
    avgdev_all{subj,1}=data_dev;
    
    avgstd_all{subj,1}=data_std;
    
    avgstd1_all{subj,1}=data_std1;
    
    avgstd3_all{subj,1}=data_std3;
    
    
end


% Group calculations

cfg=[];

grandavgdev  = ft_timelockgrandaverage(cfg, avgdev_all{:});

grandavgstd   = ft_timelockgrandaverage(cfg, avgstd_all{:});


%% Repeat for patients


submatcount_pat = 0;


for pat = 1:nfilesp
    
    
    filen = [preprocdir '/' subjs_pat{pat} '/' LFPprefix '_' preprocsteps '_' subjs_pat{pat} '_' sessions_pat{pat} '_' LFPpostfix];
    
    
    if isfile(filen)
        
        submatcount_pat = submatcount_pat + 1;
        
        D = spm_eeg_load(filen);
        
        
        %Pull out subject IDs first
        
        [filepath,~,~] = fileparts(filen);
        
        
        [~,psubj,~] = fileparts(filepath);
        
        
        %And also find diag info
        
        PatMatch=strcmp(DiagData(:,1), psubj);
        
        extdiaginfo(submatcount_pat,1) = DiagData(PatMatch==1,2);
        
        
        
        %Other processing, including baseline correction and smoothing
        
        timelock=D.fttimelock;
        
        %Baseline correct
        cfg = [];
        cfg.baseline = [-0.100 -0.001];
        [baselinecorr] = ft_timelockbaseline(cfg, timelock);
        
        
        %Go back into sources
        
        %Std
        cfg = [];
        cfg.trials = timelock.trialinfo == 7;
        data_std = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Std3
        cfg = [];
        cfg.trials = timelock.trialinfo == 4;
        data_std3 = ft_redefinetrial(cfg, baselinecorr);
        
        %Std1
        cfg = [];
        cfg.trials = timelock.trialinfo == 2;
        data_std1 = ft_redefinetrial(cfg, baselinecorr);


        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
            
            
            %Std1
            
            data_std1.avg(source,:) = smooth(data_std1.avg(source,:), 20, 'moving');
            
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into a group cell
        
        avgdev_all_pat{pat,1}=data_dev;
        
        avgstd_all_pat{pat,1}=data_std;
        
        avgstd1_all_pat{pat,1}=data_std1;
        
        avgstd3_all_pat{pat,1}=data_std3;
        
        
    end
    
end


%Calculate grand average over all pats


cfg=[];

grandavgdev_pat  = ft_timelockgrandaverage(cfg, avgdev_all_pat{:});

grandavgstd_pat  = ft_timelockgrandaverage(cfg, avgstd_all_pat{:});



%% Now we can loop through group and condition plots

nsources = 6;

meanMMNstats = struct;

meanMMNcol_all_allnodes = [];


for source=1:nsources
    


% First controls
    
    tempdevcol = zeros(nfilesc, numel(grandavgdev.time));
    
    tempstdcol = zeros(nfilesc, numel(grandavgdev.time));
    
    tempdiffcol = zeros(nfilesc, numel(grandavgdev.time));
    
    
    for isub = 1:nfilesc
        
        % Collate all individual data for variance calcs
        
        tempdevcol(isub,:) = avgdev_all{isub,1}.avg(source,:);
        
        tempstdcol(isub,:) = avgstd_all{isub,1}.avg(source,:);
        
        tempdiffcol(isub,:) = tempstdcol(isub,:)-tempdevcol(isub,:);
        
    end
    
    
    %Concat for later mean calcs
    
    tempdatacol = cat(1, tempstdcol, tempdevcol);
    
    %tempdatacol_wide = cat(2, tempstdcol, tempdevcol);

    
    
    devvar(source,:)=std(tempdevcol, [], 1)./sqrt(nfilesc);
    
    stdvar(source,:)=std(tempstdcol, [], 1)./sqrt(nfilesc);
    
    diffvar(source,:)=std(tempdiffcol, [], 1)./sqrt(nfilesc);
    
    
   
    
    
    
    %% Now Patients
    
    tempdevcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    tempstdcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    tempdiffcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    
    
    for isub = 1:nfilesp
        
        % Collate all individual data for variance calcs
        
        tempdevcol_pat(isub,:) = avgdev_all_pat{isub,1}.avg(source,:);
        
        tempstdcol_pat(isub,:) = avgstd_all_pat{isub,1}.avg(source,:);
        
        tempdiffcol_pat(isub,:) = tempstdcol_pat(isub,:)-tempdevcol_pat(isub,:);
        
    end
    
    
    %Concat for later mean calcs
    
    tempdatacol_pat = cat(1, tempstdcol_pat, tempdevcol_pat);
    
    %tempdatacol_pat_wide = cat(2, tempstdcol_pat, tempdevcol_pat);
    
    
    devvar_pat(source,:)=std(tempdevcol_pat, [], 1)./sqrt(nfilesp);
    
    stdvar_pat(source,:)=std(tempstdcol_pat, [], 1)./sqrt(nfilesp);
    
    diffvar_pat(source,:)=std(tempdiffcol_pat, [], 1)./sqrt(nfilesp);  
    
    
    % Start plots..
    
    
    % Standard in placebo, both groups
    
    figure1 = figure(1), subplot(2,3,source)
    
    hold on
    
    
    figure(1), boundedline(grandavgdev.time, mean(tempstdcol), stdvar(source,:), 'alpha', 'b'); % alpha makes bounds transparent
    
    figure(1), boundedline(grandavgdev_pat.time, mean(tempstdcol_pat), stdvar_pat(source,:), 'alpha', 'r'); % alpha makes bounds transparent
    
    
    xlim([-0.100 0.400])
    
    ax = gca;
    
    ax.YAxisLocation = 'origin';
    
    ax.YTick = [];
    
    ax.TickDir = 'out';
    
    title(grandavgdev.label(source))
    
    box off;
    
    
    %saveas(figure(2), [FigOutDir '/' 'LFPdiffcalcs_fullsample_plac_wsem_ConPat.tif']);
    
    
    
    
    % Now deviant in placebo, both groups
    
    figure2 = figure(2), subplot(2,3,source)
    
    hold on
    
    figure(2), boundedline(grandavgdev.time, mean(tempdevcol), devvar(source,:), 'alpha', 'b'); % alpha makes bounds transparent
    
    figure(2), boundedline(grandavgdev.time, mean(tempdevcol_pat), devvar_pat(source,:), 'alpha', 'r'); % alpha makes bounds transparent
    
    
    xlim([-0.100 0.400])
    
    ax = gca;
    
    ax.YAxisLocation = 'origin';
    
    ax.YTick = [];
    
    ax.TickDir = 'out';
    
    title(grandavgdev.label(source))
    
    box off;
    
    
    %saveas(figure(3), [FigOutDir '/' 'LFPdiffcalcs_fullsample_drug_wsem_ConPat.tif']);
    
    
    
    % Repeat with difference (std - dev) across conditions - important for adding significance
    %
    figure3 = figure(3), subplot(2,3,source)

%     hold on
% 
%     figure(3), boundedline(grandavgdev.time, mean(tempdiffcol), diffvar(source,:), 'alpha', 'b'); % alpha makes bounds transparent
%  
%     figure(3), boundedline(grandavgdev.time, mean(tempdiffcol_pat), diffvar_pat(source,:), 'alpha', 'r'); % alpha makes bounds transparent
%     
%     
%     xlim([-0.100 0.400])
%     
%     ax = gca;
%     
%     ax.YAxisLocation = 'origin';
%     
%     ax.YTick = 0;
%     
%     ax.TickDir = 'out';


    %Note, independent ttest now
    %First means and stds of difference
    m1diff_con=squeeze(nanmean(tempdiffcol,1));
    s1diff_con=squeeze(nanstd(tempdiffcol)./sqrt(size(tempdiffcol,1)));
    
    m1diff_pat=squeeze(nanmean(tempdiffcol_pat,1));
    s1diff_pat=squeeze(nanstd(tempdiffcol_pat)./sqrt(size(tempdiffcol_pat,1)));
    
    [h,p,ci,stats] = ttest2(tempdiffcol,tempdiffcol_pat);
    
    
    %FDR-corrected
    
    [~, ~, ~, adj_p]=fdr_bh(p(52:end));
    adj_p_wbl = cat(2, NaN(1, 51), adj_p);
    
    
    %Uncorrected 
    
    k=max(max([m1diff_con; m1diff_pat]))+max(max([s1diff_con;s1diff_pat]));k=k*1.1;
    %p=0.9.*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;

 
    
%   Don't need to do the below 
%     temp_var = strcat( 'figure',num2str(source) );
%     
%     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(3);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_con,'LineWidth',4,'Color',[0/255 0/255 255/255]); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[255 0 0]./255); legendflex({'CON','FTLD'},'fontsize',14,'fontname','Arial','fontweight','bold');
    boundedline([1:size(m1diff_con,2)],m1diff_con(1,:),s1diff_con(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[255 0 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    plot(p(adj_p_wbl<0.05),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev.label(source)]); 

    
    if ds == 1
    
    xlim([0 251]); %ylim([0 1]);%axis square
    
    else
        
    xlim([0 500])
        
    end

    %title(grandavgdev.label(source))
    
    
    % Remove the box around the plot, while we're at it:
    
    box off;
     
    
    
    %% Calculate means for MMN ANOVA
    
    %Setup mean MMN
    
    %Need to concatenate controls and patients
    
    meanMMNcol = mean(tempdatacol(:,MMNs:MMNf),2);
    
    meanMMNcol_pat = mean(tempdatacol_pat(:,MMNs:MMNf),2);
    
    
    
    %ttest
    
%     [~,p,~,stats] = ttest2(meanMMNcol, meanMMNcol_pat, 'tail', 'both');
%     
%     
%     meanMMNstats.p(source,1) = p;
%     
%     meanMMNstats.t(source,1) = stats.tstat;
%     
%     meanMMNstats.df = stats.df;
    
    
    
    %Start ANOVA - redundant now
    
%     meanMMNcol_all = cat(1, meanMMNcol, meanMMNcol_pat);
%     
%     meanMMNcol_all_allnodes = cat(2, meanMMNcol_all_allnodes, meanMMNcol_all);
%     
    
    %meanMMNcolall(1:nsubjs*4, source) = meanMMNcol;
    
%     X = cat(2, meanMMNcol_all, F1all, F2all, Sall);
%     
%     
%     [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
%     
%     
%     pgrp = Ps(1);
%     pcond = Ps(3);
%     pint = Ps(4);
%     
%     MMNpall(source,1)=pgrp;
%     MMNpall(source,2)=pcond;
%     MMNpall(source,3)=pint;
    
    
    % Mean plots
    
    
    figure4 = figure(4), subplot(2,3,source)
    
    hold on
    
    figure(4), plot([1 2], [mean(meanMMNcol(1:nfilesc)) mean(meanMMNcol(nfilesc+1:nfilesc+nfilesc))], '-b','LineWidth',1,'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',6)
    
    figure(4), plot([1 2], [mean(meanMMNcol_pat(1:nfilesp)) mean(meanMMNcol_pat(nfilesp+1:nfilesp+nfilesp))], '-r','LineWidth',1,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',6)

    
    title([data_std.label{source, 1}])

    
    %xlabel('Condition', 'FontSize', 10, 'FontWeight', 'Bold', 'Color', 'black')
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    ax.XLim = [0.5 2.5];
    ax.XTick = [1 2];
    ax.XTickLabel = {'Std' 'Dev'};
    
    
    %Write out mean data for calculating in JASP RM ANOVA
    
    
%     m1_std = nanmean(tempstdcol(:,MMNs:MMNf),2);
%     
%     m1_dev = nanmean(tempdevcol(:,MMNs:MMNf),2);
%     
%     m1_std_pat = nanmean(tempstdcol_pat(:,MMNs:MMNf),2);
%     
%     m1_dev_pat = nanmean(tempdevcol_pat(:,MMNs:MMNf),2);
%     
%     
%     
%     meanMMNcol_lw = cat(2, m1_std, m1_dev);
%     
%     meanMMNcol_lw_pat = cat(2, m1_std_pat, m1_dev_pat);
%     
%     meanMMNcol_lw_all = cat(1, meanMMNcol_lw, meanMMNcol_lw_pat);
% 
%     meanMMNcol_lw_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));
%     
%     MMNtable_lw = table(meanMMNcol_lw_all(:,1), meanMMNcol_lw_all(:,2), meanMMNcol_lw_grp);
%     
%     
%     MMNtable_lw.Properties.VariableNames = {'mean_STD', 'mean_DEV', 'Group'};
% 
%     writetable(MMNtable_lw, [FigOutDir '/' figoutprefix '_' 'LFP_' grandavgdev.label{source} '_STDDEVmean' '_' 'ConsandPats' '_GrpCondInt_' 'fullsample_forJASP.txt'], 'Delimiter','tab')
%     
    

    %Write out mean data for ttest and calculating in JASP
    
    
    conmmndiff = nanmean(tempdiffcol(:,MMNs:MMNf),2);
    
    patmmndiff = nanmean(tempdiffcol_pat(:,MMNs:MMNf),2);
    
    
    [~,p,~,stats] = ttest2(conmmndiff,patmmndiff,'tail','both');
    
    MMNpall(source,1) = p;
    MMNtall(source,1) = stats.tstat;
    
    
    %Collate together - can do all source regions at once - no problem.
    
    meanMMNdiffcol_all(:, source) = cat(1, conmmndiff, patmmndiff);


end


%Output

%Save ttest

MMNstats = struct;

MMNstats.t = MMNtall;

MMNstats.p = MMNpall;


save([FigOutDir '/' figoutprefix '_' 'LFPs_'  'MMNdiffmean_' 'ConsandPats_ttest.mat'], 'MMNstats');



%Mean data for JASP

meanMMNdiffcol_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));

meanMMNdiffcol_subgrp = cat(1, repmat(3, [nfilesc 1]), cell2mat(extdiaginfo));


MMNdifftable = table(meanMMNdiffcol_all(:,1), meanMMNdiffcol_all(:,2), meanMMNdiffcol_all(:,3), meanMMNdiffcol_all(:,4), meanMMNdiffcol_all(:,5), meanMMNdiffcol_all(:,6), meanMMNdiffcol_grp, meanMMNdiffcol_subgrp);

MMNdifftable.Properties.VariableNames = {['meanMMNcol_1'], ['meanMMNcol_2'], ['meanMMNcol_3'], ['meanMMNcol_4'], ['meanMMNcol_5'], ['meanMMNcol_6'], 'Group', 'PatSubgrp'};
%MMNdifftable.Properties.VariableNames = {join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), 'Group'};

writetable(MMNdifftable, [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MMNdiffmean_' 'ConsandPats_' 'forJASP_' 'fullsample.txt'], 'Delimiter','tab')



%Mean data for JASP

% meanMMNdiffcol_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1])); 
% 
% 
% MMNdifftable = table(meanMMNdiffcol_all(:,1), meanMMNdiffcol_all(:,2), meanMMNdiffcol_all(:,3), meanMMNdiffcol_all(:,4), meanMMNdiffcol_all(:,5), meanMMNdiffcol_all(:,6), meanMMNdiffcol_grp);
% 
% MMNdifftable.Properties.VariableNames = {['meanMMNcol_1'], ['meanMMNcol_2'], ['meanMMNcol_3'], ['meanMMNcol_4'], ['meanMMNcol_5'], ['meanMMNcol_6'], 'Group'};
% MMNdifftable.Properties.VariableNames = {join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), 'Group'};
% 
% writetable(MMNdifftable, [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MMNdiffmean_' 'ConsandPats_' 'forJASP_' 'fullsample.txt'], 'Delimiter','tab')


%Figures

saveas(figure1, [FigOutDir '/' figoutprefix '_LFPstd_fullsample_plac_wsem_ConPat.tif']);

saveas(figure2, [FigOutDir '/' figoutprefix '_LFPdev_fullsample_plac_wsem_ConPat.tif']);

saveas(figure3, [FigOutDir '/' figoutprefix '_LFPdiffcalcs_fullsample_plac_wsem_ConPat.tif']);

saveas(figure4, [FigOutDir '/' figoutprefix '_LFPmmn_fullsample_condgroup_mean_int.tif']);


%Save output table

% summrestable = table(MMNpall(:,1), MMNpall(:,2), MMNpall(:,3));
% summrestable.Properties.VariableNames = {'Group' 'Cond' 'Int'};
% 
% writetable(summrestable, [FigOutDir '/' figoutprefix '_restable_2wayANOVA_CondGrps_wInt.txt']);


%And full data table

%Can push out into rainclouds for plotting in R :)

% MMNtable = table(meanMMNcol_all_allnodes(:,1), meanMMNcol_all_allnodes(:,2), meanMMNcol_all_allnodes(:,3), meanMMNcol_all_allnodes(:,4), meanMMNcol_all_allnodes(:,5), meanMMNcol_all_allnodes(:,6), F1all, F2all, Sall);
% 
% MMNtable.Properties.VariableNames = {'meanMMNcol_LIFG', 'meanMMNcol_LSTG', 'meanMMNcol_LAUD', 'meanMMNcol_RIFG', 'meanMMNcol_RSTG', 'meanMMNcol_RAUD', 'Group', 'Trial', 'ID'};
% 
% writetable(MMNtable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'MMNmean' '_' 'ConsandPats' '_' 'fullsample.txt'], 'Delimiter','tab')
% 


%% Do again for rep3

figcount = 4;


MMNstats = struct;

meanMMNdiffcol_all = [];

MMNpall = [];

MMNtall = [];


for source=1:nsources
    
    figcount = figcount + 1;
    

% First controls
    
    tempdevcol = zeros(nfilesc, numel(grandavgdev.time));
    
    tempstdcol = zeros(nfilesc, numel(grandavgdev.time));
    
    tempdiffcol = zeros(nfilesc, numel(grandavgdev.time));
    
    
    for isub = 1:nfilesc
        
        % Collate all individual data for variance calcs
        
        tempdevcol(isub,:) = avgdev_all{isub,1}.avg(source,:);
        
        tempstdcol(isub,:) = avgstd3_all{isub,1}.avg(source,:);
        
        tempdiffcol(isub,:) = tempstdcol(isub,:)-tempdevcol(isub,:);
        
    end
    
    
    %Concat for later mean calcs
    
    tempdatacol = cat(1, tempstdcol, tempdevcol);
    
    tempdatacol_wide = cat(2, tempstdcol, tempdevcol);

    
    
    devvar(source,:)=std(tempdevcol, [], 1)./sqrt(nfilesc);
    
    stdvar(source,:)=std(tempstdcol, [], 1)./sqrt(nfilesc);
    
    diffvar(source,:)=std(tempdiffcol, [], 1)./sqrt(nfilesc);
    

    
    %% Now Patients
    
    tempdevcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    tempstdcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    tempdiffcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    
    
    for isub = 1:nfilesp
        
        % Collate all individual data for variance calcs
        
        tempdevcol_pat(isub,:) = avgdev_all_pat{isub,1}.avg(source,:);
        
        tempstdcol_pat(isub,:) = avgstd3_all_pat{isub,1}.avg(source,:);
        
        tempdiffcol_pat(isub,:) = tempstdcol_pat(isub,:)-tempdevcol_pat(isub,:);
        
    end
    
    
    %Concat for later mean calcs
    
    tempdatacol_pat = cat(1, tempstdcol_pat, tempdevcol_pat);
    
    tempdatacol_pat_wide = cat(2, tempstdcol_pat, tempdevcol_pat);
    
    
    devvar_pat(source,:)=std(tempdevcol_pat, [], 1)./sqrt(nfilesp);
    
    stdvar_pat(source,:)=std(tempstdcol_pat, [], 1)./sqrt(nfilesp);
    
    diffvar_pat(source,:)=std(tempdiffcol_pat, [], 1)./sqrt(nfilesp);  
    
    

    % Repeat with difference (std - dev) across conditions - important for adding significance
    %figure5 = figure(5), subplot(2,3,source)

%     hold on
% 
%     figure(3), boundedline(grandavgdev.time, mean(tempdiffcol), diffvar(source,:), 'alpha', 'b'); % alpha makes bounds transparent
%  
%     figure(3), boundedline(grandavgdev.time, mean(tempdiffcol_pat), diffvar_pat(source,:), 'alpha', 'r'); % alpha makes bounds transparent
%     
%     
%     xlim([-0.100 0.400])
%     
%     ax = gca;
%     
%     ax.YAxisLocation = 'origin';
%     
%     ax.YTick = 0;
%     
%     ax.TickDir = 'out';


    %Note, independent ttest now
    %First means and stds of difference
    m1diff_con=squeeze(nanmean(tempdiffcol,1));
    s1diff_con=squeeze(nanstd(tempdiffcol)./sqrt(size(tempdiffcol,1)));
    
    m1diff_pat=squeeze(nanmean(tempdiffcol_pat,1));
    s1diff_pat=squeeze(nanstd(tempdiffcol_pat)./sqrt(size(tempdiffcol_pat,1)));
    
    [h,p,ci,stats] = ttest2(tempdiffcol,tempdiffcol_pat);
    
    
    %FDR-corrected
    
    [~, ~, ~, adj_p]=fdr_bh(p(52:end));
    adj_p_wbl = cat(2, NaN(1, 51), adj_p);
    
    
    %Uncorrected 
    
    k=max(max([m1diff_con; m1diff_pat]))+max(max([s1diff_con;s1diff_pat]));k=k*1.1;
    %p=0.9.*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;

 
    
%   Don't need to do the below 
%     temp_var = strcat( 'figure',num2str(source) );
%     
%     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 0.33 0.25]);
    plot(m1diff_con,'LineWidth',4,'Color',[0/255 0/255 255/255]); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[255 0 0]./255); %legendflex({'CON','FTLD'},'fontsize',14, 'FontWeight', 'bold');
    boundedline([1:size(m1diff_con,2)],m1diff_con(1,:),s1diff_con(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[255 0 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    %plot(p(adj_p_wbl<0.05),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    
    %According to RH nodes
    ylim([-0.2 0.1])
    
    
    %set(gca, 'YTick', 'Fontsize',14, 'FontWeight', 'bold')
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',13, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial'); ylabel({'Mismatch response'; '(rep3-dev)'}, 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    %title(['rep3' '' grandavgdev.label(source)]);
    
    
    if ds == 1
        
        xlim([0 251]); %ylim([0 1]);%axis square
        
    else
        
        xlim([0 500])
        
    end
    
    %title(grandavgdev.label(source))
    
    
    % Remove the box around the plot, while we're at it:
    
    box off

    set(gca, 'TickDir', 'Out')
    
    
    saveas(figure(figcount), join([FigOutDir '/' figoutprefix '_' 'LFPs_'  'MMN3diffmean_' grandavgdev.label{source} '_ConsandPats_fullsample.tif']));

    
    
    %% Calculate means for MMN ANOVA

    
    %Now ANOVA
    
    %Need to concatenate controls and patients
    
    meanMMNcol = mean(tempdatacol(:,MMNs:MMNf),2);
    
    meanMMNcol_pat = mean(tempdatacol_pat(:,MMNs:MMNf),2);
    
    meanMMNcol_all = cat(1, meanMMNcol, meanMMNcol_pat);
    
    meanMMNcol_all_allnodes = cat(2, meanMMNcol_all_allnodes, meanMMNcol_all);
    
    
    %meanMMNcolall(1:nsubjs*4, source) = meanMMNcol;
    
%     X = cat(2, meanMMNcol_all, F1all, F2all, Sall);
%     
%     
%     [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
%     
%     
%     pgrp = Ps(1);
%     pcond = Ps(3);
%     pint = Ps(4);
%     
%     MMNpall(source,1)=pgrp;
%     MMNpall(source,2)=pcond;
%     MMNpall(source,3)=pint;
    
    
    % Mean plots
    
    
%     figure6 = figure(6), subplot(2,3,source)
%     
%     hold on
%     
%     figure(6), plot([1 2], [mean(meanMMNcol(1:nfilesc)) mean(meanMMNcol(nfilesc+1:nfilesc+nfilesc))], '-b','LineWidth',1,'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',6)
%     
%     figure(6), plot([1 2], [mean(meanMMNcol_pat(1:nfilesp)) mean(meanMMNcol_pat(nfilesp+1:nfilesp+nfilesp))], '-r','LineWidth',1,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',6)
%     
%     
%     title(['rep3' '' data_std.label{source, 1}])
%     
%     
%     %xlabel('Condition', 'FontSize', 10, 'FontWeight', 'Bold', 'Color', 'black')
%     
%     ax=gca;
%     
%     ax.YColor = 'black';
%     ax.XColor = 'black';
%     ax.FontWeight = 'bold';
%     ax.FontSize = 10;
%     
%     ax.XLim = [0.5 2.5];
%     ax.XTick = [1 2];
%     ax.XTickLabel = {'Std3' 'Dev'};
    
    
    %Write out mean data for ttest and calculating in JASP
    
    
    conmmndiff = nanmean(tempdiffcol(:,MMNs:MMNf),2);
    
    patmmndiff = nanmean(tempdiffcol_pat(:,MMNs:MMNf),2);
    
    
    [~,p,~,stats] = ttest2(conmmndiff,patmmndiff,'tail','both');
    
    MMNpall(source,1) = p;
    MMNtall(source,1) = stats.tstat;
    
    
    %Collate together - can do all source regions at once - no problem.
    
    meanMMNdiffcol_all(:, source) = cat(1, conmmndiff, patmmndiff);


end


%Output

%Save ttest

MMNstats = struct;

MMNstats.t = MMNtall;

MMNstats.p = MMNpall;


save([FigOutDir '/' figoutprefix '_' 'LFPs_'  'MMN3diffmean_' 'ConsandPats_ttest.mat'], 'MMNstats');



%Mean data for JASP

meanMMNdiffcol_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));

meanMMNdiffcol_subgrp = cat(1, repmat(3, [nfilesc 1]), cell2mat(extdiaginfo));


MMNdifftable = table(meanMMNdiffcol_all(:,1), meanMMNdiffcol_all(:,2), meanMMNdiffcol_all(:,3), meanMMNdiffcol_all(:,4), meanMMNdiffcol_all(:,5), meanMMNdiffcol_all(:,6), meanMMNdiffcol_grp, meanMMNdiffcol_subgrp);

MMNdifftable.Properties.VariableNames = {['meanMMNcol_1'], ['meanMMNcol_2'], ['meanMMNcol_3'], ['meanMMNcol_4'], ['meanMMNcol_5'], ['meanMMNcol_6'], 'Group', 'PatSubgrp'};
%MMNdifftable.Properties.VariableNames = {join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), join(['meanMMNcol_' grandavgdev.label(1)]), 'Group'};

writetable(MMNdifftable, [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MMN3diffmean_' 'ConsandPats_' 'forJASP_' 'fullsample.txt'], 'Delimiter','tab')



%Figures

%saveas(figure5, [FigOutDir '/' figoutprefix '_LFPdiffcalcs_MMN3_fullsample_plac_wsem_ConPat.tif']);

%saveas(figure6, [FigOutDir '/' figoutprefix '_LFPmmn3_fullsample_condgroup_mean_int.tif']);


%Save output table

% summrestable = table(MMNpall(:,1), MMNpall(:,2), MMNpall(:,3));
% summrestable.Properties.VariableNames = {'Group' 'Cond' 'Int'};
% 
% writetable(summrestable, [FigOutDir '/' figoutprefix 'MMN3_restable_2wayANOVA_CondGrps_wInt.txt']);


%And full data table

%Can push out into rainclouds for plotting in R :)

% MMNtable = table(meanMMNcol_all_allnodes(:,1), meanMMNcol_all_allnodes(:,2), meanMMNcol_all_allnodes(:,3), meanMMNcol_all_allnodes(:,4), meanMMNcol_all_allnodes(:,5), meanMMNcol_all_allnodes(:,6), F1all, F2all, Sall);
% 
% MMNtable.Properties.VariableNames = {'meanMMNcol_LIFG', 'meanMMNcol_LSTG', 'meanMMNcol_LAUD', 'meanMMNcol_RIFG', 'meanMMNcol_RSTG', 'meanMMNcol_RAUD', 'Group', 'Trial', 'ID'};
% 
% writetable(MMNtable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'MMN3mean' '_' 'ConsandPats' '_' 'fullsample.txt'], 'Delimiter','tab')





%% RFT
%%=========================================================================

% figcount = 10;
% 
% %Setup DM
% 
% %First, set up design matrix for RFT
% 
% DM = zeros(nfilesc+nfilesp,2);
% 
% DM(1:nfilesc,1) = repmat(1, [nfilesc 1]);
% 
% DM(nfilesc+1:end,2) = repmat(1, [nfilesp 1]);
% 
% 
% c = [1; -1];
% 
% 
% %Just do rep3
% 
% 
% for source = 4:nsources
%    
%     
%     %Start with controls
%     
%     
%     for isub = 1:nfilesc
%         
%         % Collate all individual data for variance calcs
%         
%         tempdevcol(isub,:) = avgdev_all{isub,1}.avg(source,:);
%         
%         tempstdcol(isub,:) = avgstd3_all{isub,1}.avg(source,:);
%         
%         tempdiffcol(isub,:) = tempstdcol(isub,:)-tempdevcol(isub,:);
%         
%     end
%     
%     
%     %Patients
%     
%     for isub = 1:nfilesp
%         
%         % Collate all individual data for variance calcs
%         
%         tempdevcol_pat(isub,:) = avgdev_all_pat{isub,1}.avg(source,:);
%         
%         tempstdcol_pat(isub,:) = avgstd3_all_pat{isub,1}.avg(source,:);
%         
%         tempdiffcol_pat(isub,:) = tempstdcol_pat(isub,:)-tempdevcol_pat(isub,:);
%         
%     end
%         
%     
%     
%     figcount = figcount + 1;
%         
%         
%     mmndiffall = cat(1, tempdiffcol, tempdiffcol_pat);
%     
%     
%     y = mmndiffall;
%     
%     
%     [stat,out] = RFT_GLM_contrast(DM,y,c,'t',1,1);
%     
%     
%     
%     %Save figure
%     
%     saveas(figure(figcount), [FigOutDir '/' figoutprefix '_' 'LFP_' 'rep3_' 'ControlsandPats_' grandavgdev.label{source} '_Diff_fullsample_RFT.tif']);
%     
%     
%     %Saveoutput
%     
%     save([FigOutDir '/' figoutprefix '_' 'LFP_' 'rep3_' 'ControlsandPats_' grandavgdev.label{source} '_Diff_fullsample_RFT.mat'], 'stat');
%     
%     
% end
    



%%=========================================================================
%% Disease group info
%%=========================================================================

%First determine number of bv/psp

extdiaginfo = cell2mat(extdiaginfo);

nfilesp_bv = length(find(extdiaginfo==1));
nfilesp_psp = length(find(extdiaginfo==2));


%Setup design matrix for ANOVA

%Design matrix for ANOVA..

S = [];

for s = 1:nfilesc
    
    S(s,1) = s;
    
end


F1 = repmat([1;1], [nfilesc 1]);

F2 = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesc 1]));


%Now patients

for s = 1:nfilesp_bv
    
    Sp_bv(s,1) = s+nfilesc;
    
end


for s = 1:nfilesp_psp
    
    Sp_psp(s,1) = s+nfilesp_bv+nfilesc;
    
end


F1p_bv = repmat([2;2], [nfilesp_bv 1]);
F2p_bv = cat(1, repmat(1, [nfilesp_bv 1]), repmat(2, [nfilesp_bv 1]));


F1p_psp = repmat([3;3], [nfilesp_psp 1]);
F2p_psp = cat(1, repmat(1, [nfilesp_psp 1]), repmat(2, [nfilesp_psp 1]));


%Combine design matrix(s)

Sall = cat(1, S, S, Sp_bv, Sp_bv, Sp_psp, Sp_psp);

F1all = cat(1, F1, F1p_bv, F1p_psp);
F2all = cat(1, F2, F2p_bv, F2p_psp);


MMNpall = [];


%First rep6


for source = 1:nsources
    
    % First controls
    
    tempdevcol = zeros(nfilesc, numel(grandavgdev.time));
    
    tempstdcol = zeros(nfilesc, numel(grandavgdev.time));
    
    tempdiffcol = zeros(nfilesc, numel(grandavgdev.time));
    
    
    for isub = 1:nfilesc
        
        % Collate all individual data for variance calcs
        
        tempdevcol(isub,:) = avgdev_all{isub,1}.avg(source,:);
        
        tempstdcol(isub,:) = avgstd_all{isub,1}.avg(source,:);
        
        tempdiffcol(isub,:) = tempstdcol(isub,:)-tempdevcol(isub,:);
        
    end
    
    
    %Concat for later mean calcs
    
    tempdatacol = cat(1, tempstdcol, tempdevcol);
    
    
    devvar(source,:)=std(tempdevcol, [], 1)./sqrt(nfilesc);
    
    stdvar(source,:)=std(tempstdcol, [], 1)./sqrt(nfilesc);
    
    diffvar(source,:)=std(tempdiffcol, [], 1)./sqrt(nfilesc);
    
    
    
    %% Now Patients
    
    tempdevcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    tempstdcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    tempdiffcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    
    
    for isub = 1:nfilesp
        
        % Collate all individual data for variance calcs
        
        tempdevcol_pat(isub,:) = avgdev_all_pat{isub,1}.avg(source,:);
        
        tempstdcol_pat(isub,:) = avgstd_all_pat{isub,1}.avg(source,:);
        
        tempdiffcol_pat(isub,:) = tempstdcol_pat(isub,:)-tempdevcol_pat(isub,:);
        
    end
    
    
    %Concat for later mean calcs
    
    %tempdatacol_pat = cat(1, tempstdcol_pat, tempdevcol_pat);
    
    
    %     devvar_pat(source,:)=std(tempdevcol_pat, [], 1)./sqrt(nfilesp);
    %
    %     stdvar_pat(source,:)=std(tempstdcol_pat, [], 1)./sqrt(nfilesp);
    %
    %     diffvar_pat(source,:)=std(tempdiffcol_pat, [], 1)./sqrt(nfilesp);
    
    
    
    % Now separate the groups (fingers crossed)
    
    
    tempdevcol_pat_bv = tempdevcol_pat(extdiaginfo==1,:);
    
    tempstdcol_pat_bv = tempstdcol_pat(extdiaginfo==1,:);
    
    tempdiffcol_pat_bv = tempdiffcol_pat(extdiaginfo==1,:);
    
    tempdatacol_pat_bv = cat(1, tempstdcol_pat_bv, tempdevcol_pat_bv);
    
    
    tempdevcol_pat_psp = tempdevcol_pat(extdiaginfo==2,:);
    
    tempstdcol_pat_psp = tempstdcol_pat(extdiaginfo==2,:);
    
    tempdiffcol_pat_psp = tempdiffcol_pat(extdiaginfo==2,:);
    
    tempdatacol_pat_psp = cat(1, tempstdcol_pat_psp, tempdevcol_pat_psp);
    
    
    %Variance calcs
    
    devvar_pat_bv(source,:)=std(tempdevcol_pat_bv, [], 1)./sqrt(nfilesp_bv);
    
    stdvar_pat_bv(source,:)=std(tempstdcol_pat_bv, [], 1)./sqrt(nfilesp_bv);
    
    diffvar_pat_bv(source,:)=std(tempdiffcol_pat_bv, [], 1)./sqrt(nfilesp_bv);
    
    
    devvar_pat_psp(source,:)=std(tempdevcol_pat_psp, [], 1)./sqrt(nfilesp_psp);
    
    stdvar_pat_psp(source,:)=std(tempstdcol_pat_psp, [], 1)./sqrt(nfilesp_psp);
    
    diffvar_pat_psp(source,:)=std(tempdiffcol_pat_psp, [], 1)./sqrt(nfilesp_psp);
    
    
    %Means across three groups for difference
    
    m1diff_con=squeeze(nanmean(tempdiffcol,1));
    s1diff_con=squeeze(nanstd(tempdiffcol)./sqrt(size(tempdiffcol,1)));
    
    m1diff_pat_bv=squeeze(nanmean(tempdiffcol_pat_bv,1));
    s1diff_pat_bv=squeeze(nanstd(tempdiffcol_pat_bv)./sqrt(size(tempdiffcol_pat_bv,1)));
    
    m1diff_pat_psp=squeeze(nanmean(tempdiffcol_pat_psp,1));
    s1diff_pat_psp=squeeze(nanstd(tempdiffcol_pat_psp)./sqrt(size(tempdiffcol_pat_psp,1)));
    
    
    %Get to plots
    
    %Start with subgroup plotting
    
    [h,p,ci,stats]=ttest2(squeeze(tempdiffcol_pat_bv),squeeze(tempdiffcol_pat_psp));
    k=max(max([m1diff_pat_bv(1,:); m1diff_pat_psp(1,:)]))+max(max([s1diff_pat_bv(1,:);s1diff_pat_psp(1,:)])); k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    figure(11), subplot(2,3,source);
    
    figure(11); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1diff_pat_bv(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1diff_pat_psp(1,:),'LineWidth',4,'Color',[0 255 0]./255); legendflex({'bvFTD','PSP'},'fontsize',10)
    boundedline([1:size(m1diff_pat_bv,2)],m1diff_pat_bv(1,:),s1diff_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_psp,2)],m1diff_pat_psp(1,:),s1diff_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN')
    title([grandavgdev.label(source) ' MMN in bv/PSP'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    %Now including three groups (i.e. adding controls)
    %Can go straight into analysis
    
    diff_allgrps = cat(1, tempdiffcol, tempdiffcol_pat_bv, tempdiffcol_pat_psp);
    
    grpst_allgrps = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp_bv 1]), repmat(3, [nfilesp_psp 1]));
    
    
    %Note, one-way anova now
    for t_win = 1:length(diff_allgrps(1,:))
        p(1,t_win) = anova1(diff_allgrps(:,t_win),grpst_allgrps, 'off');
    end
    
    k=max(max([m1diff_con; m1diff_pat_bv; m1diff_pat_psp]))+max(max([s1diff_con; s1diff_pat_bv; s1diff_pat_psp])); k=k*1.1;
    
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %Start figure
    
    figure(12), subplot(2,3,source);
    
    
    figure(12);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_con,'LineWidth',4,'Color',[0 0 255]./255); hold on
    plot(m1diff_pat_bv,'LineWidth',4,'Color',[255 0 0]./255);
    plot(m1diff_pat_psp,'LineWidth',4,'Color',[0/255 255/255 0/255]); legendflex({'CON', 'bvFTD','PSP'},'fontsize',10)
    
    boundedline([1:size(m1diff_con,2)],m1diff_con(1,:),s1diff_con(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_bv,2)],m1diff_pat_bv(1,:),s1diff_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_psp,2)],m1diff_pat_psp(1,:),s1diff_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN')
    title([grandavgdev.label(source) ' MMN'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
        

    
    
    
end

%Save output

%Figures

saveas(figure(11), [FigOutDir '/' figoutprefix '_LFPdiffcalcs_MMN_fullsample_plac_wsem_PatSubgrps.tif']);

saveas(figure(12), [FigOutDir '/' figoutprefix '_LFPdiffcalcs_MMN_fullsample_drug_wsem_ConplusSubgrps.tif']);




%Now rep3

MMNpall = [];

meanMMNstats = struct;


for source = 1:nsources
    
    % First controls
    
    tempdevcol = zeros(nfilesc, numel(grandavgdev.time));
    
    tempstdcol = zeros(nfilesc, numel(grandavgdev.time));
    
    tempdiffcol = zeros(nfilesc, numel(grandavgdev.time));
    
    
    for isub = 1:nfilesc
        
        % Collate all individual data for variance calcs
        
        tempdevcol(isub,:) = avgdev_all{isub,1}.avg(source,:);
        
        tempstdcol(isub,:) = avgstd3_all{isub,1}.avg(source,:);
        
        tempdiffcol(isub,:) = tempstdcol(isub,:)-tempdevcol(isub,:);
        
    end
    
    
    %Concat for later mean calcs
    
    tempdatacol = cat(1, tempstdcol, tempdevcol);
    
    
    devvar(source,:)=std(tempdevcol, [], 1)./sqrt(nfilesc);
    
    stdvar(source,:)=std(tempstdcol, [], 1)./sqrt(nfilesc);
    
    diffvar(source,:)=std(tempdiffcol, [], 1)./sqrt(nfilesc);
    
    
    
    %% Now Patients
    
    tempdevcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    tempstdcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    tempdiffcol_pat = zeros(nfilesp, numel(grandavgdev.time));
    
    
    
    for isub = 1:nfilesp
        
        % Collate all individual data for variance calcs
        
        tempdevcol_pat(isub,:) = avgdev_all_pat{isub,1}.avg(source,:);
        
        tempstdcol_pat(isub,:) = avgstd3_all_pat{isub,1}.avg(source,:);
        
        tempdiffcol_pat(isub,:) = tempstdcol_pat(isub,:)-tempdevcol_pat(isub,:);
        
    end
    
    
    %Concat for later mean calcs
    
    %tempdatacol_pat = cat(1, tempstdcol_pat, tempdevcol_pat);
    
    
    %     devvar_pat(source,:)=std(tempdevcol_pat, [], 1)./sqrt(nfilesp);
    %
    %     stdvar_pat(source,:)=std(tempstdcol_pat, [], 1)./sqrt(nfilesp);
    %
    %     diffvar_pat(source,:)=std(tempdiffcol_pat, [], 1)./sqrt(nfilesp);
    
    
    
    % Now separate the groups (fingers crossed)
    
    
    tempdevcol_pat_bv = tempdevcol_pat(extdiaginfo==1,:);
    
    tempstdcol_pat_bv = tempstdcol_pat(extdiaginfo==1,:);
    
    tempdiffcol_pat_bv = tempdiffcol_pat(extdiaginfo==1,:);
    
    tempdatacol_pat_bv = cat(1, tempstdcol_pat_bv, tempdevcol_pat_bv);
    
    
    tempdevcol_pat_psp = tempdevcol_pat(extdiaginfo==2,:);
    
    tempstdcol_pat_psp = tempstdcol_pat(extdiaginfo==2,:);
    
    tempdiffcol_pat_psp = tempdiffcol_pat(extdiaginfo==2,:);
    
    tempdatacol_pat_psp = cat(1, tempstdcol_pat_psp, tempdevcol_pat_psp);
    
    
    %Variance calcs
    
    devvar_pat_bv(source,:)=std(tempdevcol_pat_bv, [], 1)./sqrt(nfilesp_bv);
    
    stdvar_pat_bv(source,:)=std(tempstdcol_pat_bv, [], 1)./sqrt(nfilesp_bv);
    
    diffvar_pat_bv(source,:)=std(tempdiffcol_pat_bv, [], 1)./sqrt(nfilesp_bv);
    
    
    devvar_pat_psp(source,:)=std(tempdevcol_pat_psp, [], 1)./sqrt(nfilesp_psp);
    
    stdvar_pat_psp(source,:)=std(tempstdcol_pat_psp, [], 1)./sqrt(nfilesp_psp);
    
    diffvar_pat_psp(source,:)=std(tempdiffcol_pat_psp, [], 1)./sqrt(nfilesp_psp);
    
    
    %Means across three groups for difference
    
    m1diff_con=squeeze(nanmean(tempdiffcol,1));
    s1diff_con=squeeze(nanstd(tempdiffcol)./sqrt(size(tempdiffcol,1)));
    
    m1diff_pat_bv=squeeze(nanmean(tempdiffcol_pat_bv,1));
    s1diff_pat_bv=squeeze(nanstd(tempdiffcol_pat_bv)./sqrt(size(tempdiffcol_pat_bv,1)));
    
    m1diff_pat_psp=squeeze(nanmean(tempdiffcol_pat_psp,1));
    s1diff_pat_psp=squeeze(nanstd(tempdiffcol_pat_psp)./sqrt(size(tempdiffcol_pat_psp,1)));
    
    
    %Get to plots
    
    %Start with subgroup plotting
    
    [h,p,ci,stats]=ttest2(squeeze(tempdiffcol_pat_bv),squeeze(tempdiffcol_pat_psp));
    k=max(max([m1diff_pat_bv(1,:); m1diff_pat_psp(1,:)]))+max(max([s1diff_pat_bv(1,:);s1diff_pat_psp(1,:)])); k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    figure(13), subplot(2,3,source);
    
    figure(13); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1diff_pat_bv(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1diff_pat_psp(1,:),'LineWidth',4,'Color',[0 255 0]./255); legendflex({'bvFTD','PSP'},'fontsize',10)
    boundedline([1:size(m1diff_pat_bv,2)],m1diff_pat_bv(1,:),s1diff_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_psp,2)],m1diff_pat_psp(1,:),s1diff_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN3')
    title([grandavgdev.label(source) ' MMN3 in bv/PSP'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
   
    
    %Now including three groups (i.e. adding controls)
    %Can go straight into analysis
    
    diff_allgrps = cat(1, tempdiffcol, tempdiffcol_pat_bv, tempdiffcol_pat_psp);
    
    grpst_allgrps = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp_bv 1]), repmat(3, [nfilesp_psp 1]));
    
    
    %Note, one-way anova now
    for t_win = 1:length(diff_allgrps(1,:))
        p(1,t_win) = anova1(diff_allgrps(:,t_win),grpst_allgrps, 'off');
    end
    
    k=max(max([m1diff_con; m1diff_pat_bv; m1diff_pat_psp]))+max(max([s1diff_con; s1diff_pat_bv; s1diff_pat_psp])); k=k*1.1;
    
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %Start figure
    
    figure(14), subplot(2,3,source);
    
    
    figure(14);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_con,'LineWidth',4,'Color',[0 0 255]./255); hold on
    plot(m1diff_pat_bv,'LineWidth',4,'Color',[255 0 0]./255);
    plot(m1diff_pat_psp,'LineWidth',4,'Color',[0/255 255/255 0/255]); legendflex({'CON', 'bvFTD','PSP'},'fontsize',10)
    
    boundedline([1:size(m1diff_con,2)],m1diff_con(1,:),s1diff_con(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_bv,2)],m1diff_pat_bv(1,:),s1diff_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_psp,2)],m1diff_pat_psp(1,:),s1diff_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN3')
    title([grandavgdev.label(source) ' MMN3'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Need to concatenate controls and patients
    
    meanMMNcol = mean(tempdiffcol(:,MMNs:MMNf),2);
    
    meanMMNcol_pat_bv = mean(tempdiffcol_pat_bv(:,MMNs:MMNf),2);
    
    meanMMNcol_pat_psp = mean(tempdiffcol_pat_psp(:,MMNs:MMNf),2);
    
    
    meanMMNcol_all = cat(1, meanMMNcol, meanMMNcol_pat_bv, meanMMNcol_pat_psp);

    
    %ANOVA1
    
    [p,anovatab,stats] = anova1(meanMMNcol_all,grpst_allgrps, 'off');
    
    
    meanMMNstats.p(source) = p;
    
    meanMMNstats.F(source) = anovatab{2,6};
    
    meanMMNstats.d(1,1) = anovatab{2,3};

    meanMMNstats.d(1,2) = stats.df;
    
    
    %For ANOVA factorial
    
    
    %meanMMNcol_all_allnodes = cat(2, meanMMNcol_all_allnodes, meanMMNcol_all);
    
    
    %meanMMNcolall(1:nsubjs*4, source) = meanMMNcol;
    
%     X = cat(2, meanMMNcol_all, F1all, F2all, Sall);
%     
%     
%     [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
%     
%     
%     pgrp = Ps(1);
%     pcond = Ps(3);
%     pint = Ps(4);
%     
%     MMNpall(source,1)=cell2mat(pgrp);
%     MMNpall(source,2)=cell2mat(pcond);
%     MMNpall(source,3)=cell2mat(pint);
    
    
end


%Save output

save([FigOutDir '/' figoutprefix 'MMN3_restable_ANOVA1_ConplusSubgrps.mat'], 'meanMMNstats');


%Tables

%Save output table

% summrestable = table(MMNpall(:,1), MMNpall(:,2), MMNpall(:,3));
% summrestable.Properties.VariableNames = {'Group' 'Cond' 'Int'};
% 
% writetable(summrestable, [FigOutDir '/' figoutprefix 'MMN3_restable_2wayANOVA_CondGrps_wInt_ConplusSubgrps.txt']);
% 

%Figures

saveas(figure(13), [FigOutDir '/' figoutprefix '_LFPdiffcalcs_MMN3_fullsample_plac_wsem_PatSubgrps.tif']);

saveas(figure(14), [FigOutDir '/' figoutprefix '_LFPdiffcalcs_MMN3_fullsample_drug_wsem_ConplusSubgrps.tif']);



%% Wow, Finished

sprintf('finito, output saved to %s \n', FigOutDir)


end