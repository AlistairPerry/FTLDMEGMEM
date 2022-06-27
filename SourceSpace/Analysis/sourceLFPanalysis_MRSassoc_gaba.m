function [MMNstats] = sourceLFPanalysis_MRSassoc_gaba

%% Setup

% mutualsubjscon='Con_BothSess_subjs.txt';
% mutualsubjspat='Pat_BothSess_subjs.txt';
% LFPplaclistcon, LFPdruglistcon, LFPplaclistpat, LFPdruglistpat)
% [MMNpall] = sourceLFPanalysis_basiccond_druggroup_rmanov('Con_BothSess_subjs.txt', 'Pat_BothSess_subjs.txt', 'Con_AllP_Sess.txt', 'Con_AllD_Sess.txt', 'Pat_P_fnames.txt', 'Pat_D_fnames.txt')

%Open up SPM dirs (and field trip), and others

%Add paths


addpath('/imaging/rowe/users/ap09/Toolbox')

addpath('/imaging/rowe/users/ap09/Toolbox/ekMEG/ntadscripts')

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest')

addpath('/imaging/rowe/users/ap09/Toolbox/fieldtrip')

%addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/fieldtrip/'); %EK


addpath(genpath('/imaging/rowe/users/ap09/Toolbox/boundedline-pkg'))

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/VBA-toolbox'))




ft_defaults()

%spm('defaults', 'eeg')


%% Directory structure

FigOutDir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/analysis/SourceSpace/LFPs_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis';

mkdir(FigOutDir)


%figoutprefix = 'adv_ssst_newmaxf_fixICA_wfids_nodipcor';


%When removing that outlier
figoutprefix = 'adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19_FIX';



%% MRS, Demographics, and diagnostic info

% Load MRS table and diagnosis information

MRStable = readtable('/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/Clinical/FTD-MEM_testscores.xls', 'FileType', 'spreadsheet', 'Sheet', 'MRS', 'Range', 'A1:T46');


Demtable = readtable('/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/Clinical/FTD-MEM_testscores.xls', 'FileType', 'spreadsheet', 'Sheet', 'Demographics', 'Range', 'A2:T47');


% Can call it directly as it wont change
DiagData = readcell('/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/PatDiagInfo.txt');


%Other parameter setup

nsources = 6;
%nsources = 8; %w IPC

ds = 1;


devtrl = 1;

std3trl = 4;

stdtrl = 7;


xticks=[1 50 100 150 200 251];

xlabels={'-100' '0' '100' '200' '300' '400'};



%% Subjects
%Going to be hard to compile them altogether, but oh well

preprocdir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/analysis/SourceSpace/LFPs_icafixes/LFPs_COH_wrad_allMEGch_wfids/B_Data';

preprocsteps = 'fmraeMaffffdS';

LFPprefix = 'LFP6inv1';

LFPpostfix = 'rov123_sss.mat'; %note i've put it to the one for comb sess



Pat_P = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Pat_P.txt';

%Pat_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Pat_D.txt';


%And when removing outlier

%Pat_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Pat_D_nop19.txt';


%Both P10 (small # trials) and P19 (noisy)

Pat_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Pat_D_nop10p19.txt';


%Just p10 to remove

%Pat_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Pat_D_nop10.txt';


fid = fopen(Pat_P);
PlacData_pat = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
PlacData_pat = PlacData_pat{1,1};

fid = fopen(Pat_D);
DrugData_pat = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
DrugData_pat = DrugData_pat{1,1};



%% Patients


submatcount_pat = 0;



for pplasubj = 1:length(PlacData_pat(:,1))
    
    
    filen = [preprocdir '/' PlacData_pat{pplasubj} '/' LFPprefix '_' preprocsteps '_' PlacData_pat{pplasubj} '_' 's1' '_' LFPpostfix];
    
    
    filen_d = [preprocdir '/' PlacData_pat{pplasubj} '/' LFPprefix '_' preprocsteps '_' PlacData_pat{pplasubj} '_' 's2' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount_pat = submatcount_pat + 1;
        
        
        allpats_bothsess{submatcount_pat,1} = PlacData_pat{pplasubj};
        
        
        %Also find diag info
        
        PatMatch=strcmpi(DiagData(:,1), PlacData_pat{pplasubj});

        extdiaginfo(submatcount_pat,1) = DiagData(PatMatch==1,2);
        
        
        D=spm_eeg_load(filen);
        
        
        D2=spm_eeg_load(filen_d);
        
        
        
        plasubj = PlacData_pat{pplasubj};
        
        
        
        %Standards
        
%         dipmatch = strcmp(plasubj, PlacToFlip(:,1));
        
        
        
%         if nnz(dipmatch)==1
%             
%             nodetoflip = cell2mat(PlacToFlip(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 flipndtrl = -(D(nodetoflip, :, stdtrl));
%                 
%                 D(nodetoflip, :, stdtrl) = flipndtrl;
%                 
%                 %D(nodetoflip, :, devtrl) = -(D(nodetoflip, :, devtrl));
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     flipndtrl = -(D(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     D(nodetoflip(nodeflip), :, stdtrl) = flipndtrl;
%                     
%                     %D(nodetoflip(nodeflip), :, devtrl) = -(D(nodetoflip(nodeflip), :, devtrl));
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
%         
%         
%         
%         %Deviants
%         
%         dipmatch = strcmp(plasubj, PlacToFlip_Dev(:,1));
%         
%         
%         
%         if nnz(dipmatch)==1
%             
%             
%             nodetoflip = cell2mat(PlacToFlip_Dev(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 %D(nodetoflip, :, stdtrl) = -(D(nodetoflip, :, stdtrl));
%                 
%                 flipndtrl = -(D(nodetoflip, :, devtrl));
%                 
%                 D(nodetoflip, :, devtrl) = flipndtrl;
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     %D(nodetoflip(nodeflip), :, stdtrl) = -(D(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     flipndtrl = -(D(nodetoflip(nodeflip), :, devtrl));                    
%                     
%                     D(nodetoflip(nodeflip), :, devtrl) = flipndtrl;
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
        
        
        %Other processing, including baseline correction and smoothing
        
        timelock=D.fttimelock;
        
        
        %Baseline correct
        cfg = [];
        cfg.baseline = [-0.100 -0.001];
        [baselinecorr] = ft_timelockbaseline(cfg, timelock);
        
        
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
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20,'moving');
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20,'moving');
            
            
            %Std3
            
            data_std1.avg(source,:) = smooth(data_std1.avg(source,:), 20,'moving');
            
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');           
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgstd1_pat_col(submatcount_pat,:, :) = permute(data_std1.avg, [3 2 1]);
        
        avgdev_pat_col(submatcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev_pat{submatcount_pat,1}=data_dev;
        
        avgstd_pat{submatcount_pat,1}=data_std;
        
        avgstd3_pat{submatcount_pat,1}=data_std3;
        
        
        
        % Now with automatic correction
        
%         timelockac = timelock;
%         
%         for source = 1:nsources
%             
%             if timelockac.trial(devtrl, source, 85) > 0
%                 
%                 timelockac.trial(devtrl, source, :) = -(timelockac.trial(devtrl, source, :));
%                 
%             end
%             
%             
%             if timelockac.trial(stdtrl, source, 85) > 0
%                 
%                 timelockac.trial(stdtrl, source, :) = -(timelockac.trial(stdtrl, source, :));
%                 
%             end
%             
%         end
%         
%         
%         %Baseline correct
%         cfg = [];
%         cfg.baseline = [-0.100 -0.001];
%         [baselinecorr_ac] = ft_timelockbaseline(cfg, timelockac);
%         
%         
%         %Std
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 7;
%         data_std_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         %Dev
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 1;
%         data_dev_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         
%         %Smoothing
%         
%         for source = 1:nsources
%             
%             %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             %Dev
%             
%             data_dev_ac.avg(source,:) = smooth(data_dev_ac.avg(source,:),20,'moving');
%             
%         end
%         
%         
%         %Combine into 3D matrix
%         
%         avgstd_pat_col_ac(submatcount_pat,:, :) = permute(data_std_ac.avg, [3 2 1]);
%         
%         avgdev_pat_col_ac(submatcount_pat,:, :) = permute(data_dev_ac.avg, [3 2 1]);
        

        
        
        %And repeat for drug session
        

        %Other processing, including baseline correction and smoothing
        
        timelock=D2.fttimelock;
        
        
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
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20,'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20,'moving');
            
            
            %Std1
            
            data_std1.avg(source,:) = smooth(data_std1.avg(source,:), 20,'moving');
            
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_d_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_d_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgstd1_d_pat_col(submatcount_pat,:, :) = permute(data_std1.avg, [3 2 1]);
        
        avgdev_d_pat_col(submatcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev_d_pat{submatcount_pat,1}=data_dev;
        
        avgstd_d_pat{submatcount_pat,1}=data_std;
        
        avgstd3_d_pat{submatcount_pat,1}=data_std3;
        
        
        
        % Auto-correction
%         
%         timelockac = D2.fttimelock;
%         
%         for source = 1:nsources
%             
%             if timelockac.trial(devtrl, source, 85) > 0
%                 
%                 timelockac.trial(devtrl, source, :) = -(timelockac.trial(devtrl, source, :));
%                 
%             end
%             
%             
%             if timelockac.trial(stdtrl, source, 85) > 0
%                 
%                 timelockac.trial(stdtrl, source, :) = -(timelockac.trial(stdtrl, source, :));
%                 
%             end
%             
%         end
%         
%         
%         %Baseline correct
%         cfg = [];
%         cfg.baseline = [-0.100 -0.001];
%         [baselinecorr_ac] = ft_timelockbaseline(cfg, timelockac);
%         
%         
%         %Std
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 7;
%         data_std_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         %Dev
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 1;
%         data_dev_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         
%         %Smoothing
%         
%         for source = 1:nsources
%             
%             %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             %Dev
%             
%             data_dev_ac.avg(source,:) = smooth(data_dev_ac.avg(source,:),20,'moving');
%             
%         end
%         
%         
%         %Combine into 3D matrix
%         
%         avgstd_d_pat_col_ac(submatcount_pat,:, :) = permute(data_std_ac.avg, [3 2 1]);
%         
%         avgdev_d_pat_col_ac(submatcount_pat,:, :) = permute(data_dev_ac.avg, [3 2 1]);
        
        
    end
    
    
end



for pdrugsubj = 1:length(DrugData_pat(:,1))
    
    
    filen = [preprocdir '/' DrugData_pat{pdrugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData_pat{pdrugsubj} '_' 's2' '_' LFPpostfix];
    
    filen_d = [preprocdir '/' DrugData_pat{pdrugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData_pat{pdrugsubj} '_' 's1' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount_pat = submatcount_pat + 1;
        
        
        allpats_bothsess{submatcount_pat,1} = DrugData_pat{pdrugsubj};
        
        
        
        %Ext diag info again
        
        PatMatch=strcmpi(DiagData(:,1), DrugData_pat{pdrugsubj});

        extdiaginfo(submatcount_pat,1) = DiagData(PatMatch==1,2);
        
        
        %Start again
        
        D=spm_eeg_load(filen);
        
        D2=spm_eeg_load(filen_d);
        
        
        
        drugsubj = DrugData_pat{pdrugsubj};
        
        
%         dipmatch = strcmp(drugsubj, PlacToFlip(:,1));
%         
%         
%         
%         %Standards
%         
%         if nnz(dipmatch)==1
%             
%             
%             nodetoflip = cell2mat(PlacToFlip(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 flipndtrl = -(D(nodetoflip, :, stdtrl));
%                 
%                 D(nodetoflip, :, stdtrl) = flipndtrl;
%                 
%                 %D(nodetoflip, :, devtrl) = -(D(nodetoflip, :, devtrl));
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     flipndtrl = -(D(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     D(nodetoflip(nodeflip), :, stdtrl) = flipndtrl;
%                     
%                     %D(nodetoflip(nodeflip), :, devtrl) = -(D(nodetoflip(nodeflip), :, devtrl));
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
%         
%         
%         
%         %Deviants
%         
%         
%         dipmatch = strcmp(drugsubj, PlacToFlip_Dev(:,1));
%         
%         
%         if nnz(dipmatch)==1
%             
%             
%             nodetoflip = cell2mat(PlacToFlip_Dev(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 %D(nodetoflip, :, stdtrl) = -(D(nodetoflip, :, stdtrl));
%                 
%                 flipndtrl = -(D(nodetoflip, :, devtrl));
%                 
%                 D(nodetoflip, :, devtrl) = flipndtrl;
%                 
%             
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     %D(nodetoflip(nodeflip), :, stdtrl) = -(D(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     flipndtrl = -(D(nodetoflip(nodeflip), :, devtrl));
%                     
%                     D(nodetoflip(nodeflip), :, devtrl) = flipndtrl;
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
        
        
        %Other processing, including baseline correction and smoothing
        
        timelock=D.fttimelock;
        
        %Baseline correct
        cfg = [];
        cfg.baseline = [-0.100 -0.001];
        [baselinecorr] = ft_timelockbaseline(cfg, timelock);
        
        
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
        
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20,'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20,'moving');
            
            
            %Std1
            
            data_std1.avg(source,:) = smooth(data_std1.avg(source,:), 20,'moving');
            
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgstd1_pat_col(submatcount_pat,:, :) = permute(data_std1.avg, [3 2 1]);
        
        avgdev_pat_col(submatcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        %Combine into a group cell
        
        avgdev_pat{submatcount_pat,1}=data_dev;
        
        avgstd_pat{submatcount_pat,1}=data_std;
        
        avgstd3_pat{submatcount_pat,1}=data_std3;
        
        
        %Autocorrect

        
%         timelockac = timelock;
%         
%         for source = 1:nsources
%             
%             if timelockac.trial(devtrl, source, 85) > 0
%                 
%                 timelockac.trial(devtrl, source, :) = -(timelockac.trial(devtrl, source, :));
%                 
%             end
%             
%             
%             if timelockac.trial(stdtrl, source, 85) > 0
%                 
%                 timelockac.trial(stdtrl, source, :) = -(timelockac.trial(stdtrl, source, :));
%                 
%             end
%             
%         end
%         
%         
%         %Baseline correct
%         cfg = [];
%         cfg.baseline = [-0.100 -0.001];
%         [baselinecorr_ac] = ft_timelockbaseline(cfg, timelockac);
%         
%         
%         %Std
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 7;
%         data_std_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         %Dev
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 1;
%         data_dev_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         
%         %Smoothing
%         
%         for source = 1:nsources
%             
%             %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             
%                         %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             %Dev
%             
%             data_dev_ac.avg(source,:) = smooth(data_dev_ac.avg(source,:),20,'moving');
%             
%         end
%         
%         
%         %Combine into 3D matrix
%         
%         avgstd_pat_col_ac(submatcount_pat,:, :) = permute(data_std_ac.avg, [3 2 1]);
%         
%         avgdev_pat_col_ac(submatcount_pat,:, :) = permute(data_dev_ac.avg, [3 2 1]);
        

        
        
        %And repeat for drug session
        
        
%         dipmatch = strcmp(drugsubj, DrugToFlip(:,1));
%         
%         
%         if nnz(dipmatch)==1
%             
%             nodetoflip = cell2mat(DrugToFlip(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 flipndtrl = -(D2(nodetoflip, :, stdtrl));
%                 
%                 D2(nodetoflip, :, stdtrl) = flipndtrl;
%                 
%                 %D2(nodetoflip, :, devtrl) = -(D2(nodetoflip, :, devtrl));
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     flipndtrl = -(D2(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     D2(nodetoflip(nodeflip), :, stdtrl) = flipndtrl;
%                     
%                     %D2(nodetoflip(nodeflip), :, devtrl) = -(D2(nodetoflip(nodeflip), :, devtrl));
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
%         
%         
%         %Deviants
%         
%         dipmatch = strcmp(drugsubj, DrugToFlip_Dev(:,1));
%         
%         
%         if nnz(dipmatch)==1
%             
%             nodetoflip = cell2mat(DrugToFlip_Dev(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 %D2(nodetoflip, :, stdtrl) = -(D2(nodetoflip, :, stdtrl));
%                 
%                 flipndtrl = -(D2(nodetoflip, :, devtrl));
%                 
%                 D2(nodetoflip, :, devtrl) = flipndtrl;
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     %D2(nodetoflip(nodeflip), :, stdtrl) = -(D2(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     flipndtrl = -(D2(nodetoflip(nodeflip), :, devtrl));
%                     
%                     D2(nodetoflip(nodeflip), :, devtrl) = flipndtrl;
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
        
        
        
        %Other processing, including baseline correction and smoothing
        
        timelock=D2.fttimelock;
        
        
        %Baseline correct
        cfg = [];
        cfg.baseline = [-0.100 -0.001];
        [baselinecorr] = ft_timelockbaseline(cfg, timelock);
        
        
        
        %Go back into sources
        
        %Std
        cfg = [];
        cfg.trials = timelock.trialinfo == 7;
        data_std = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Std3 - oops
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 4;
%         data_std = ft_redefinetrial(cfg, baselinecorr);
        
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
        
        
        %Combine into 3D matrix
        
        avgstd_d_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_d_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgstd1_d_pat_col(submatcount_pat,:, :) = permute(data_std1.avg, [3 2 1]);
        
        avgdev_d_pat_col(submatcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev_d_pat{submatcount_pat,1}=data_dev;
        
        avgstd_d_pat{submatcount_pat,1}=data_std;
        
        avgstd3_d_pat{submatcount_pat,1}=data_std3;
        
        
        
        % Auto-correction
        
%         timelockac = D2.fttimelock;
%         
%         for source = 1:nsources
%         
%         if timelockac.trial(devtrl, source, 85) > 0
%             
%             timelockac.trial(devtrl, source, :) = -(timelockac.trial(devtrl, source, :));
%             
%         end
%         
%         
%         if timelockac.trial(stdtrl, source, 85) > 0
%             
%             timelockac.trial(stdtrl, source, :) = -(timelockac.trial(stdtrl, source, :));
%             
%         end
%         
%         end
%         
%         %Baseline correct
%         cfg = [];
%         cfg.baseline = [-0.100 -0.001];
%         [baselinecorr_ac] = ft_timelockbaseline(cfg, timelockac);
%         
%         
%         %Std
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 7;
%         data_std_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         %Dev
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 1;
%         data_dev_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         
%         %Smoothing
%         
%         for source = 1:nsources
%             
%             %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             %Dev
%             
%             data_dev_ac.avg(source,:) = smooth(data_dev_ac.avg(source,:),20,'moving');
%             
%         end
%         
%         
%         %Combine into 3D matrix
%         
%         avgstd_d_pat_col_ac(submatcount_pat,:, :) = permute(data_std_ac.avg, [3 2 1]);
%         
%         avgdev_d_pat_col_ac(submatcount_pat,:, :) = permute(data_dev_ac.avg, [3 2 1]);
%         
        
        
    end
    
    
end



%% Find matching data in MRS table
% Remember to include only those with both MEG and MRS - likely some will
% not have MRS


MRStable_pats = MRStable(21:end,:);

MRStable_patids_lc = lower(MRStable_pats.Participant);


DEMtable_pats = Demtable(21:end,:);


%Now those w/ both sess

match = 0;

npats_bothsess = size(avgstd_pat,1);

pats_bothsess_keep = zeros(1, npats_bothsess);



for isub = 1:npats_bothsess
    
    
    % subj
    
    subj = allpats_bothsess{isub};
    
    SubjMatch=strcmp(MRStable_patids_lc, [subj]);
    

    if ~isempty(SubjMatch) && ~isnan(MRStable_pats.gaba_ifg_corr(SubjMatch))
        
        
    pats_bothsess_keep(isub) = 1;
    
    match = match + 1;
    
    
    ExtTable_gabaifg_bothsess(match,:) = MRStable_pats.gaba_ifg(SubjMatch==1);
    
    ExtTable_gabaifgcor_bothsess(match,:) = MRStable_pats.gaba_ifg_corr(SubjMatch==1);
    
    ExtTable_gabaocccor_bothsess(match,:) = MRStable_pats.gaba_occ_corr(SubjMatch==1);
    
    
    ExtTable_Dems_Age(match,:) = DEMtable_pats.AgeVisit1(SubjMatch==1);
    
    ExtTable_Dems_Sex{match,:} = DEMtable_pats.Gender(SubjMatch==1);

    
    end
    
    
end



pats_bothsess_keep = pats_bothsess_keep';


%Remove string from sex info

ExtTable_Dems_Sex_mat = zeros(length(ExtTable_Dems_Age),1);

ExtTable_Dems_Sex_mat(string(ExtTable_Dems_Sex)=='F', :) = 1; %Code F as 1

ExtTable_Dems_Sex_mat(string(ExtTable_Dems_Sex)=='M', :) = -1;





%% Remove subjects without MRS

%Grand averages

cfg=[];


%Plac

grandavgdev_pat  = ft_timelockgrandaverage(cfg, avgdev_pat{pats_bothsess_keep==1});

grandavgstd_pat  = ft_timelockgrandaverage(cfg, avgstd_pat{pats_bothsess_keep==1});


%Drug

grandavgdev_d_pat  = ft_timelockgrandaverage(cfg, avgdev_d_pat{pats_bothsess_keep==1});

grandavgstd_d_pat  = ft_timelockgrandaverage(cfg, avgstd_d_pat{pats_bothsess_keep==1});




% And collapsed data

%Non-corrected

avgstd_pat_col_keep = avgstd_pat_col(pats_bothsess_keep==1,:,:);

avgstd3_pat_col_keep = avgstd3_pat_col(pats_bothsess_keep==1,:,:);

avgstd1_pat_col_keep = avgstd1_pat_col(pats_bothsess_keep==1,:,:);

avgdev_pat_col_keep = avgdev_pat_col(pats_bothsess_keep==1,:,:);


avgstd_d_pat_col_keep = avgstd_d_pat_col(pats_bothsess_keep==1,:,:);

avgstd3_d_pat_col_keep = avgstd3_d_pat_col(pats_bothsess_keep==1,:,:);

avgstd1_d_pat_col_keep = avgstd1_d_pat_col(pats_bothsess_keep==1,:,:);

avgdev_d_pat_col_keep = avgdev_d_pat_col(pats_bothsess_keep==1,:,:);



%Auto-corrected

% avgstd_pat_col_ac_keep = avgstd_pat_col_ac(pats_bothsess_keep==1,:,:);
% 
% avgdev_pat_col_ac_keep = avgdev_pat_col_ac(pats_bothsess_keep==1,:,:);
% 
% 
% avgstd_d_pat_col_ac_keep = avgstd_d_pat_col_ac(pats_bothsess_keep==1,:,:);
% 
% avgdev_d_pat_col_ac_keep = avgdev_d_pat_col_ac(pats_bothsess_keep==1,:,:);



% And lastly diagnostic info

extdiaginfoold = extdiaginfo;
extdiaginfo = extdiaginfo(pats_bothsess_keep==1,:);




%% Basic plotting and analyses - association with MRS scores


%First setup

ds = 1;

if ds == 1
    
    %Determine MMN start/finish
    
    MMNs = 226./2;
    MMNf = 276./2;
    
    %And also p3
    
    p3s = 375./2;
    p3f = 425./2;
    
else
    
    MMNs = 226;
    MMNf = 276;
    
    p3s = 275;
    p3f = 325;
    
end




%%=========================================================================
%% Now MRS assoc with Placebo
%%=========================================================================


MMNrall = [];

MMNpall = [];


MMNrall_cor = [];

MMNpall_cor = [];


%%=========================================================================
%% First rep6
%%=========================================================================


for source = 1:nsources
   
    
    %Placebo - Pats
    
    diff_pat = avgstd_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
    
    m1diff_pat=squeeze(nanmean(diff_pat,1));
    s1diff_pat=squeeze(nanstd(diff_pat)./sqrt(size(diff_pat,1)));
    
        
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifg_bothsess);
%         
%         pall(1,t_win) = p(1,2);


        [b, stats] = robustfit(ExtTable_gabaifg_bothsess, diff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(1), subplot(2,3,source);
    
    hold on
    
    figure(1);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    %plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); %legendflex({'MEM','PLA'},'fontsize',14)
    %boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) 'PLA MRSIFGgaba Assoc in MMN Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
   
    
    
    % Corrected values
    
    ball = [];
    pall = [];
    
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifgcor_bothsess);
%         
%         pall(1,t_win) = p(1,2);

        [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, diff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        ball(1,t_win) = b(2,1);

        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(2), subplot(2,3,source);
    
    hold on
    
    figure(2);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    %plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); %legendflex({'MEM','PLA'},'fontsize',14)
    %boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    
    %Positive assoc
%     plot(p(ball>0),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    %Negative
%     plot(p(ball<0),'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)


    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) 'PLA MRSIFGgabacor Assoc in MMN Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    % Mean calculations
    
    m1difftemp=squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    %m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %X = [ExtTable_gabaifg_bothsess - mean(ExtTable_gabaifg_bothsess)];

    [b, stats] = robustfit(ExtTable_gabaifg_bothsess, m1difftemp);
    
    
    MMNrall(source,1) = b(2,1);
    
    MMNpall(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifg_bothsess, m1difftemp);
    
    MMNpall_coeff(source,1) = p(1,2);

    
    figure(3), subplot(2,3,source)
    
    
    hold on
    
    
    figure(3), scatter(ExtTable_gabaifg_bothsess, m1difftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    
    set(gca, 'Fontsize',12); xlabel('Mean MMN PLA'); ylabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    % !! Repeat for corrected calculations
    
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifgcor_bothsess);
%     
%     
%     MMNrall_cor(source,1) = r(1,2);
%     
%     MMNpall_cor(source,1) = p(1,2);
%     


    [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, m1difftemp);
    
    
    MMNrall_cor(source,1) = b(2,1);
    
    MMNpall_cor(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifgcor_bothsess, m1difftemp);
    
    MMNpall_cor_coeff(source,1) = p(1,2);
    

    figure(4), subplot(2,3,source)
    
    
    hold on
    
    
    figure(4), scatter(ExtTable_gabaifgcor_bothsess, m1difftemp, 'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
    
    x = ExtTable_gabaifgcor_bothsess;
    
    plot(x,b(1)+b(2)*x,'k', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN PLA'); xlabel('gaba rIFG');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
end



%% Rep3

for source = 1:nsources
   
    
    %Placebo
    
    diff_pat = avgstd3_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
    
    m1diff_pat=squeeze(nanmean(diff_pat,1));
    s1diff_pat=squeeze(nanstd(diff_pat)./sqrt(size(diff_pat,1)));
    
   
        
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifg_bothsess);
%         
%         pall(1,t_win) = p(1,2);


        [b, stats] = robustfit(ExtTable_gabaifg_bothsess, diff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(5), subplot(2,3,source);
    
    hold on
    
    figure(5);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    %plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); %legendflex({'MEM','PLA'},'fontsize',14)
    %boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) 'PLA MRSIFGgaba Assoc in MMN3 Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
   
    
    
    % Corrected values
    
    ball = [];
    pall = [];
    
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifgcor_bothsess);
%         
%         pall(1,t_win) = p(1,2);

        [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, diff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        ball(1,t_win) = b(2,1);

        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(6), subplot(2,3,source);
    
    hold on
    
    figure(6);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    %plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); %legendflex({'MEM','PLA'},'fontsize',14)
    %boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    
    %Positive assoc
%     plot(p(ball>0),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    %Negative
%     plot(p(ball<0),'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)


    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) 'PLA MRSIFGgabacor Assoc in MMN3 Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    % Mean calculations
    
    
    %m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
    m1difftemp=squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));

    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %X = [ExtTable_gabaifg_bothsess - mean(ExtTable_gabaifg_bothsess)];

    [b, stats] = robustfit(ExtTable_gabaifg_bothsess, m1difftemp);
    
    
    MMNrall_rep3(source,1) = b(2,1);
    
    MMNpall_rep3(source,1) = stats.p(2,1);
    

    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifg_bothsess, m1difftemp);
    
    MMNpall_rep3_coeff(source,1) = p(1,2);
    
    
    figure(7), subplot(2,3,source)
    
    
    hold on
    
    
    figure(7), scatter(ExtTable_gabaifg_bothsess, m1difftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    
    set(gca, 'Fontsize',12); xlabel('PLA Mean MMN3'); ylabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    % !! Repeat for corrected calculations
    
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifgcor_bothsess);
%     
%     
%     MMNrall_cor(source,1) = r(1,2);
%     
%     MMNpall_cor(source,1) = p(1,2);
%     


    [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, m1difftemp);
    
    
    MMNrall_cor_rep3(source,1) = b(2,1);
    
    MMNpall_cor_rep3(source,1) = stats.p(2,1);
    
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifgcor_bothsess, m1difftemp);
    
    MMNpall_cor_rep3_coeff(source,1) = p(1,2);
    
    

    figure(8), subplot(2,3,source)
    
    
    hold on
    
    
    figure(8), scatter(ExtTable_gabaifgcor_bothsess, m1difftemp, 'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
    
    x = ExtTable_gabaifgcor_bothsess;
    
    plot(x,b(1)+b(2)*x,'k', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('PLA Mean MMN3'); xlabel('gaba rIFG');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
end




% Save output

%Tables

%Rep6
[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor;
MMNstats.r = MMNrall_cor;

MMNstats.corrp = MMNpall_cor_coeff;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN_PLA_Pats_fullsample_stats.mat'], 'MMNstats');


%Rep3

[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor_rep3);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor_rep3;
MMNstats.r = MMNrall_cor_rep3;

MMNstats.corrp = MMNpall_cor_rep3_coeff;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN3_PLA_Pats_fullsample_stats.mat'], 'MMNstats');

   

%Save Figures


saveas(figure(1), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'PLA_Pats' '_fullsample.tif']);


saveas(figure(2), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'PLA_Pats' '_fullsample.tif']);


saveas(figure(3), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'meanMMN_PLA_Pats' '_fullsample.tif']);


saveas(figure(4), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_PLA_Pats' '_fullsample.tif']);


savefig(figure(4), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_PLA_Pats' '_fullsample.fig']);


saveas(figure(5), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'PLA3_Pats' '_fullsample.tif']);


saveas(figure(6), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'PLA3_Pats' '_fullsample.tif']);


saveas(figure(7), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'meanMMN3_PLA_Pats' '_fullsample.tif']);


saveas(figure(8), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_PLA_Pats' '_fullsample.tif']);


savefig(figure(8), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_PLA_Pats' '_fullsample.fig']);


%%=========================================================================
%% Now MRS assoc with drug difference
%%=========================================================================


MMNrall = [];

MMNpall = [];


MMNrall_cor = [];

MMNpall_cor = [];


%% First rep6

for source = 1:nsources
   
    
    %Placebo
    
    diff_pat = avgstd_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
    
    m1diff_pat=squeeze(nanmean(diff_pat,1));
    s1diff_pat=squeeze(nanstd(diff_pat)./sqrt(size(diff_pat,1)));
    
    
    %Drug
    
    diff_d_pat = avgstd_d_pat_col_keep(:,:,source)-avgdev_d_pat_col_keep(:,:,source);   
    
    m1diff_d_pat = squeeze(nanmean(diff_d_pat,1));
    s1diff_d_pat = squeeze(nanstd(diff_d_pat)./sqrt(size(diff_d_pat,1)));
    
    
    %PLA - MEM
    
    diff_drugdiff_pat = diff_pat - diff_d_pat;
        
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_drugdiff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifg_bothsess);
%         
%         pall(1,t_win) = p(1,2);


        [b, stats] = robustfit(ExtTable_gabaifg_bothsess, diff_drugdiff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(9), subplot(2,3,source);
    
    hold on
    
    figure(9);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); legendflex({'MEM','PLA'},'fontsize',14)
    boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) ' MRSIFGgaba Assoc in MMN Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
   
    
    
    % Corrected values
    
    ball = [];
    pall = [];
    
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_drugdiff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifgcor_bothsess);
%         
%         pall(1,t_win) = p(1,2);

        [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, diff_drugdiff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        ball(1,t_win) = b(2,1);

        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(10), subplot(2,3,source);
    
    hold on
    
    figure(10);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); legendflex({'MEM','PLA'},'fontsize',14)
    boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    
    %Positive assoc
%     plot(p(ball>0),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    %Negative
%     plot(p(ball<0),'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)


    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) ' MRSIFGgabacor Assoc in MMN Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    % Mean calculations
    
    
    m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %X = [ExtTable_gabaifg_bothsess - mean(ExtTable_gabaifg_bothsess)];

    [b, stats] = robustfit(ExtTable_gabaifg_bothsess, m1drugdifftemp);
    
    
    MMNrall(source,1) = b(2,1);
    
    MMNpall(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifg_bothsess, m1drugdifftemp);
    
    MMNpall_coeff(source,1) = p(1,2);

    
    figure(11), subplot(2,3,source)
    
    
    hold on
    
    
    figure(11), scatter(ExtTable_gabaifg_bothsess, m1drugdifftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    
    set(gca, 'Fontsize',12); xlabel('Mean Drug Diff MMN'); ylabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    % !! Repeat for corrected calculations
    
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifgcor_bothsess);
%     
%     
%     MMNrall_cor(source,1) = r(1,2);
%     
%     MMNpall_cor(source,1) = p(1,2);
%     


    [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, m1drugdifftemp);
    
    
    MMNrall_cor(source,1) = b(2,1);
    
    MMNpall_cor(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [r,p,rlo,rup] = corrcoef(ExtTable_gabaifgcor_bothsess, m1drugdifftemp);
    
    MMNpall_cor_coeff(source,1) = p(1,2);
    
    MMNpall_cor_coeff_r2(source,1) = r(1,2)*r(1,2);
    
    MMNpall_cor_coeff_CI(source,1:2) = [rlo(1,2); rup(1,2)];
    

    figure(12), subplot(2,3,source)
    
    
    hold on
    
    
    figure(12), scatter(ExtTable_gabaifgcor_bothsess, m1drugdifftemp, 'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
    
    x = ExtTable_gabaifgcor_bothsess;
    
    plot(x,b(1)+b(2)*x,'k', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean Drug Diff MMN (PLA-MEM)'); xlabel('gaba rIFG');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    %Save output table for LMM 
    
    lmm_druggabaint = [];
    
    lmm_druggabaint(1:length(ExtTable_gabaifgcor_bothsess), 1) = 1:length(ExtTable_gabaifgcor_bothsess);
    
    lmm_druggabaint(:,2) = ExtTable_gabaifgcor_bothsess;
    
    lmm_druggabaint(:,3) = squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,4) = squeeze(nanmean(diff_d_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,5) = m1drugdifftemp;
    
    lmm_druggabaint(:,6) = cell2mat(extdiaginfo);
    
    
    lmm_druggabaint_table = table(lmm_druggabaint(:,1), lmm_druggabaint(:,2), lmm_druggabaint(:,3), lmm_druggabaint(:,4), lmm_druggabaint(:,5), lmm_druggabaint(:,6));
    
    
    lmm_druggabaint_table.Properties.VariableNames = {'Subject', 'GABA', 'PLA', 'MEM', 'Diff', 'Diag'};
    
    
    writetable(lmm_druggabaint_table, [FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN_DrugDiff_Pats_fullsample_stats_' 'LMMtableforR_' grandavgdev_pat.label{source} '.txt'], 'Delimiter','tab');
    
    
end



%% Rep3

for source = 1:nsources
   
    
    %Placebo
    
    diff_pat = avgstd3_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
    
    m1diff_pat=squeeze(nanmean(diff_pat,1));
    s1diff_pat=squeeze(nanstd(diff_pat)./sqrt(size(diff_pat,1)));
    
    
    %Drug
    
    diff_d_pat = avgstd3_d_pat_col_keep(:,:,source)-avgdev_d_pat_col_keep(:,:,source);   
    
    m1diff_d_pat = squeeze(nanmean(diff_d_pat,1));
    s1diff_d_pat = squeeze(nanstd(diff_d_pat)./sqrt(size(diff_d_pat,1)));
    
    
    %PLA - MEM
    
    diff_drugdiff_pat = diff_pat - diff_d_pat;
        
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_drugdiff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifg_bothsess);
%         
%         pall(1,t_win) = p(1,2);


        [b, stats] = robustfit(ExtTable_gabaifg_bothsess, diff_drugdiff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(13), subplot(2,3,source);
    
    hold on
    
    figure(13);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); legendflex({'MEM','PLA'},'fontsize',14)
    boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) ' MRSIFGgaba Assoc in MMN3 Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
   
    
    
    % Corrected values
    
    ball = [];
    pall = [];
    
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_drugdiff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifgcor_bothsess);
%         
%         pall(1,t_win) = p(1,2);

        [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, diff_drugdiff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        ball(1,t_win) = b(2,1);

        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(14), subplot(2,3,source);
    
    hold on
    
    figure(14);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); legendflex({'MEM','PLA'},'fontsize',14)
    boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    
    %Positive assoc
%     plot(p(ball>0),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    %Negative
%     plot(p(ball<0),'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)


    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) ' MRSIFGgabacor Assoc in MMN3 Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    % Mean calculations
    
    
    m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %X = [ExtTable_gabaifg_bothsess - mean(ExtTable_gabaifg_bothsess)];

    [b, stats] = robustfit(ExtTable_gabaifg_bothsess, m1drugdifftemp);
    
    
    MMNrall_rep3(source,1) = b(2,1);
    
    MMNpall_rep3(source,1) = stats.p(2,1);
    

    %Norm corr
    
    [r,p,rlo,rup] = corrcoef(ExtTable_gabaifg_bothsess, m1drugdifftemp);
    
    MMNpall_rep3_coeff(source,1) = p(1,2);
    
    MMNpall_rep3_coeff_r2(source,1) = r(1,2)*r(1,2);
    
    MMNpall_rep3_coeff_CI(source,1:2) = [rlo(1,2); rup(1,2)];
    
    
    
    figure(15), subplot(2,3,source)
    
    
    hold on
    
    
    figure(15), scatter(ExtTable_gabaifg_bothsess, m1drugdifftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    
    set(gca, 'Fontsize',12); xlabel('Mean Drug Diff MMN3'); ylabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    %Save output table for LMM
    
    lmm_druggabaint = [];
    
    lmm_druggabaint(1:length(ExtTable_gabaifg_bothsess), 1) = 1:length(ExtTable_gabaifg_bothsess);
    
    lmm_druggabaint(:,2) = ExtTable_gabaifg_bothsess;
    
    lmm_druggabaint(:,3) = squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,4) = squeeze(nanmean(diff_d_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,5) = m1drugdifftemp;
    
    lmm_druggabaint(:,6) = cell2mat(extdiaginfo);
    
    lmm_druggabaint(:,7) = ExtTable_Dems_Age;
    
    
    lmm_druggabaint_table = table(lmm_druggabaint(:,1), lmm_druggabaint(:,2), lmm_druggabaint(:,3), lmm_druggabaint(:,4), lmm_druggabaint(:,5), lmm_druggabaint(:,6), lmm_druggabaint(:,7));
    
    
    lmm_druggabaint_table.Properties.VariableNames = {'Subject', 'GABA', 'PLA', 'MEM', 'Diff', 'Diag', 'Age'};
    
    
    writetable(lmm_druggabaint_table, [FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabaassoc_' 'meanMMN3_DrugDiff_Pats_fullsample_stats_' 'LMMtableforR_' grandavgdev_pat.label{source} '.txt'], 'Delimiter','tab');
   
    
    
    % !! Repeat for corrected calculations
    
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifgcor_bothsess);
%     
%     
%     MMNrall_cor(source,1) = r(1,2);
%     
%     MMNpall_cor(source,1) = p(1,2);
%     


    [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, m1drugdifftemp);
    
    
    MMNrall_cor_rep3(source,1) = b(2,1);
    
    MMNpall_cor_rep3(source,1) = stats.p(2,1);
    
    
    [r,p,rlo,rup] = corrcoef(ExtTable_gabaifgcor_bothsess, m1drugdifftemp);
    
    MMNpall_cor_rep3_coeff(source,1) = p(1,2);
    
    MMNpall_cor_rep3_coeff_r2(source,1) = r(1,2)*r(1,2);
    
    MMNpall_cor_rep3_coeff_CI(source,1:2) = [rlo(1,2); rup(1,2)];


    figure(16), subplot(2,3,source)
    
    
    hold on
    
    
    figure(16), scatter(ExtTable_gabaifgcor_bothsess, m1drugdifftemp, 'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
    
    x = ExtTable_gabaifgcor_bothsess;
    
    plot(x,b(1)+b(2)*x,'k', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean Drug Diff MMN3 (PLA-MEM)'); xlabel('gaba rIFG');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    %Save output table for LMM
    
    lmm_druggabaint = [];
    
    lmm_druggabaint(1:length(ExtTable_gabaifgcor_bothsess), 1) = 1:length(ExtTable_gabaifg_bothsess);
    
    lmm_druggabaint(:,2) = ExtTable_gabaifgcor_bothsess;
    
    lmm_druggabaint(:,3) = squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,4) = squeeze(nanmean(diff_d_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,5) = m1drugdifftemp;
    
    lmm_druggabaint(:,6) = cell2mat(extdiaginfo);
    
    lmm_druggabaint(:,7) = ExtTable_Dems_Age;
    
    
    lmm_druggabaint_table = table(lmm_druggabaint(:,1), lmm_druggabaint(:,2), lmm_druggabaint(:,3), lmm_druggabaint(:,4), lmm_druggabaint(:,5), lmm_druggabaint(:,6), lmm_druggabaint(:,7));
    
    
    lmm_druggabaint_table.Properties.VariableNames = {'Subject', 'GABA', 'PLA', 'MEM', 'Diff', 'Diag', 'Age'};
    
    
    writetable(lmm_druggabaint_table, [FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN3_DrugDiff_Pats_fullsample_stats_' 'LMMtableforR_' grandavgdev_pat.label{source} '.txt'], 'Delimiter','tab');
   
    
    
    %Cooks distance
    
    temptbl = table(ExtTable_gabaifgcor_bothsess, m1drugdifftemp);
    
    lm = fitlm(temptbl,'m1drugdifftemp~ExtTable_gabaifgcor_bothsess');
    
    findout = find((lm.Diagnostics.CooksDistance)>4*mean(lm.Diagnostics.CooksDistance));
  
    MMNpall_cor_rep3_cookd{1, source} = findout;
    
    
   
    
end




% Save output

%Tables

%Rep6
[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor;
MMNstats.b = MMNrall_cor;

MMNstats.corrp = MMNpall_cor_coeff;

MMNstats.r2 = MMNpall_cor_coeff_r2;

MMNstats.CI = MMNpall_cor_coeff_CI;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN_DrugDiff_Pats_fullsample_stats.mat'], 'MMNstats');


%Rep3

%Not-corrected 

[h, ~, ~, adj_p]=fdr_bh(MMNpall_rep3);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_rep3;

MMNstats.b = MMNrall_rep3;

MMNstats.corrp = MMNpall_rep3_coeff;

MMNstats.r2 = MMNpall_rep3_coeff_r2;

MMNstats.CI = MMNpall_rep3_coeff_CI;



save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabaassoc_' 'meanMMN3_DrugDiff_Pats_fullsample_stats.mat'], 'MMNstats');



%Corrected

[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor_rep3);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor_rep3;

MMNstats.b = MMNrall_cor_rep3;

MMNstats.corrp = MMNpall_cor_rep3_coeff;

MMNstats.r2 = MMNpall_cor_rep3_coeff_r2;

MMNstats.CI = MMNpall_cor_rep3_coeff_CI;

%And now save cooks d

MMNstats.cookd = MMNpall_cor_rep3_cookd;

save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN3_DrugDiff_Pats_fullsample_stats.mat'], 'MMNstats');



%Save Figures


saveas(figure(9), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'DrugDiff_Pats' '_fullsample.tif']);


saveas(figure(10), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'DrugDiff_Pats' '_fullsample.tif']);


saveas(figure(11), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'meanMMN_DrugDiff_Pats' '_fullsample.tif']);


saveas(figure(12), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_DrugDiff_Pats' '_fullsample.tif']);


savefig(figure(12), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_DrugDiff_Pats' '_fullsample.fig']);


saveas(figure(13), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'DrugDiff3_Pats' '_fullsample.tif']);


saveas(figure(14), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'DrugDiff3_Pats' '_fullsample.tif']);


saveas(figure(15), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'meanMMN3_DrugDiff_Pats' '_fullsample.tif']);


saveas(figure(16), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_DrugDiff_Pats' '_fullsample.tif']);


savefig(figure(16), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_DrugDiff_Pats' '_fullsample.fig']);


%%=========================================================================
%% Add Controls
%%=========================================================================


Con_P = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Con_P.txt';

Con_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Con_D.txt';



fid = fopen(Con_P);
PlacData = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
PlacData = PlacData{1,1};

fid = fopen(Con_D);
DrugData = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
DrugData = DrugData{1,1};


%%=========================================================================
%% Load Data
%%========================================================================


submatcount = 0;



for plasubj = 1:length(PlacData(:,1))
    
    
    filen = [preprocdir '/' PlacData{plasubj} '/' LFPprefix '_' preprocsteps '_' PlacData{plasubj} '_' 's1' '_' LFPpostfix];
    
    
    filen_d = [preprocdir '/' PlacData{plasubj} '/' LFPprefix '_' preprocsteps '_' PlacData{plasubj} '_' 's2' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount = submatcount + 1;
        
        
        allcons_bothsess{submatcount,1} = PlacData{plasubj};
        
        
        
        D=spm_eeg_load(filen);
        
        
        D2=spm_eeg_load(filen_d);
        
        
        
        
        
        %Standards
        
%         dipmatch = strcmp(plasubj, PlacToFlip(:,1));
        
        
        
%         if nnz(dipmatch)==1
%             
%             nodetoflip = cell2mat(PlacToFlip(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 flipndtrl = -(D(nodetoflip, :, stdtrl));
%                 
%                 D(nodetoflip, :, stdtrl) = flipndtrl;
%                 
%                 %D(nodetoflip, :, devtrl) = -(D(nodetoflip, :, devtrl));
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     flipndtrl = -(D(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     D(nodetoflip(nodeflip), :, stdtrl) = flipndtrl;
%                     
%                     %D(nodetoflip(nodeflip), :, devtrl) = -(D(nodetoflip(nodeflip), :, devtrl));
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
%         
%         
%         
%         %Deviants
%         
%         dipmatch = strcmp(plasubj, PlacToFlip_Dev(:,1));
%         
%         
%         
%         if nnz(dipmatch)==1
%             
%             
%             nodetoflip = cell2mat(PlacToFlip_Dev(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 %D(nodetoflip, :, stdtrl) = -(D(nodetoflip, :, stdtrl));
%                 
%                 flipndtrl = -(D(nodetoflip, :, devtrl));
%                 
%                 D(nodetoflip, :, devtrl) = flipndtrl;
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     %D(nodetoflip(nodeflip), :, stdtrl) = -(D(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     flipndtrl = -(D(nodetoflip(nodeflip), :, devtrl));                    
%                     
%                     D(nodetoflip(nodeflip), :, devtrl) = flipndtrl;
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
        
        
        %Other processing, including baseline correction and smoothing
        
        timelock=D.fttimelock;
        
        
        %Baseline correct
        cfg = [];
        cfg.baseline = [-0.100 -0.001];
        [baselinecorr] = ft_timelockbaseline(cfg, timelock);
        
        
        %Std
        cfg = [];
        cfg.trials = timelock.trialinfo == 7;
        data_std = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Std3
        cfg = [];
        cfg.trials = timelock.trialinfo == 4;
        data_std3 = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20,'moving');
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20,'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_col(submatcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_col(submatcount,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_col(submatcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev{submatcount,1}=data_dev;
        
        avgstd{submatcount,1}=data_std;
        
        avgstd3{submatcount,1}=data_std3;
        
        
        
        % Now with automatic correction
        
%         timelockac = timelock;
%         
%         for source = 1:nsources
%             
%             if timelockac.trial(devtrl, source, 85) > 0
%                 
%                 timelockac.trial(devtrl, source, :) = -(timelockac.trial(devtrl, source, :));
%                 
%             end
%             
%             
%             if timelockac.trial(stdtrl, source, 85) > 0
%                 
%                 timelockac.trial(stdtrl, source, :) = -(timelockac.trial(stdtrl, source, :));
%                 
%             end
%             
%         end
%         
%         
%         %Baseline correct
%         cfg = [];
%         cfg.baseline = [-0.100 -0.001];
%         [baselinecorr_ac] = ft_timelockbaseline(cfg, timelockac);
%         
%         
%         %Std
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 7;
%         data_std_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         %Dev
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 1;
%         data_dev_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         
%         %Smoothing
%         
%         for source = 1:nsources
%             
%             %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             %Dev
%             
%             data_dev_ac.avg(source,:) = smooth(data_dev_ac.avg(source,:),20,'moving');
%             
%         end
%         
%         
%         %Combine into 3D matrix
%         
%         avgstd_pat_col_ac(submatcount_pat,:, :) = permute(data_std_ac.avg, [3 2 1]);
%         
%         avgdev_pat_col_ac(submatcount_pat,:, :) = permute(data_dev_ac.avg, [3 2 1]);
        

        
        
        %And repeat for drug session
        

        %Other processing, including baseline correction and smoothing
        
        timelock=D2.fttimelock;
        
        
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
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20,'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20,'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_d_col(submatcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_d_col(submatcount,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_d_col(submatcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev_d{submatcount,1}=data_dev;
        
        avgstd_d{submatcount,1}=data_std;
        
        avgstd3_d{submatcount,1}=data_std3;
        
        
        
        % Auto-correction
%         
%         timelockac = D2.fttimelock;
%         
%         for source = 1:nsources
%             
%             if timelockac.trial(devtrl, source, 85) > 0
%                 
%                 timelockac.trial(devtrl, source, :) = -(timelockac.trial(devtrl, source, :));
%                 
%             end
%             
%             
%             if timelockac.trial(stdtrl, source, 85) > 0
%                 
%                 timelockac.trial(stdtrl, source, :) = -(timelockac.trial(stdtrl, source, :));
%                 
%             end
%             
%         end
%         
%         
%         %Baseline correct
%         cfg = [];
%         cfg.baseline = [-0.100 -0.001];
%         [baselinecorr_ac] = ft_timelockbaseline(cfg, timelockac);
%         
%         
%         %Std
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 7;
%         data_std_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         %Dev
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 1;
%         data_dev_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         
%         %Smoothing
%         
%         for source = 1:nsources
%             
%             %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             %Dev
%             
%             data_dev_ac.avg(source,:) = smooth(data_dev_ac.avg(source,:),20,'moving');
%             
%         end
%         
%         
%         %Combine into 3D matrix
%         
%         avgstd_d_pat_col_ac(submatcount_pat,:, :) = permute(data_std_ac.avg, [3 2 1]);
%         
%         avgdev_d_pat_col_ac(submatcount_pat,:, :) = permute(data_dev_ac.avg, [3 2 1]);
        
        
    end
    
    
end



for drugsubj = 1:length(DrugData(:,1))
    
    
    filen = [preprocdir '/' DrugData{drugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData{drugsubj} '_' 's2' '_' LFPpostfix];
    
    filen_d = [preprocdir '/' DrugData{drugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData{drugsubj} '_' 's1' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount = submatcount + 1;
        
        
        allcons_bothsess{submatcount,1} = DrugData{drugsubj};
        
        
        
        D=spm_eeg_load(filen);
        
        D2=spm_eeg_load(filen_d);
        
        
        
        drugsubj = DrugData{drugsubj};
        
        
%         dipmatch = strcmp(drugsubj, PlacToFlip(:,1));
%         
%         
%         
%         %Standards
%         
%         if nnz(dipmatch)==1
%             
%             
%             nodetoflip = cell2mat(PlacToFlip(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 flipndtrl = -(D(nodetoflip, :, stdtrl));
%                 
%                 D(nodetoflip, :, stdtrl) = flipndtrl;
%                 
%                 %D(nodetoflip, :, devtrl) = -(D(nodetoflip, :, devtrl));
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     flipndtrl = -(D(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     D(nodetoflip(nodeflip), :, stdtrl) = flipndtrl;
%                     
%                     %D(nodetoflip(nodeflip), :, devtrl) = -(D(nodetoflip(nodeflip), :, devtrl));
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
%         
%         
%         
%         %Deviants
%         
%         
%         dipmatch = strcmp(drugsubj, PlacToFlip_Dev(:,1));
%         
%         
%         if nnz(dipmatch)==1
%             
%             
%             nodetoflip = cell2mat(PlacToFlip_Dev(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 %D(nodetoflip, :, stdtrl) = -(D(nodetoflip, :, stdtrl));
%                 
%                 flipndtrl = -(D(nodetoflip, :, devtrl));
%                 
%                 D(nodetoflip, :, devtrl) = flipndtrl;
%                 
%             
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     %D(nodetoflip(nodeflip), :, stdtrl) = -(D(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     flipndtrl = -(D(nodetoflip(nodeflip), :, devtrl));
%                     
%                     D(nodetoflip(nodeflip), :, devtrl) = flipndtrl;
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
        
        
        %Other processing, including baseline correction and smoothing
        
        timelock=D.fttimelock;
        
        %Baseline correct
        cfg = [];
        cfg.baseline = [-0.100 -0.001];
        [baselinecorr] = ft_timelockbaseline(cfg, timelock);
        
        
        %Std
        cfg = [];
        cfg.trials = timelock.trialinfo == 7;
        data_std = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Std3
        cfg = [];
        cfg.trials = timelock.trialinfo == 4;
        data_std3 = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20,'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20,'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_col(submatcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_col(submatcount,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_col(submatcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        %Combine into a group cell
        
        avgdev{submatcount,1}=data_dev;
        
        avgstd{submatcount,1}=data_std;
        
        avgstd3{submatcount,1}=data_std3;
        
        
        %Autocorrect

        
%         timelockac = timelock;
%         
%         for source = 1:nsources
%             
%             if timelockac.trial(devtrl, source, 85) > 0
%                 
%                 timelockac.trial(devtrl, source, :) = -(timelockac.trial(devtrl, source, :));
%                 
%             end
%             
%             
%             if timelockac.trial(stdtrl, source, 85) > 0
%                 
%                 timelockac.trial(stdtrl, source, :) = -(timelockac.trial(stdtrl, source, :));
%                 
%             end
%             
%         end
%         
%         
%         %Baseline correct
%         cfg = [];
%         cfg.baseline = [-0.100 -0.001];
%         [baselinecorr_ac] = ft_timelockbaseline(cfg, timelockac);
%         
%         
%         %Std
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 7;
%         data_std_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         %Dev
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 1;
%         data_dev_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         
%         %Smoothing
%         
%         for source = 1:nsources
%             
%             %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             
%                         %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             %Dev
%             
%             data_dev_ac.avg(source,:) = smooth(data_dev_ac.avg(source,:),20,'moving');
%             
%         end
%         
%         
%         %Combine into 3D matrix
%         
%         avgstd_pat_col_ac(submatcount_pat,:, :) = permute(data_std_ac.avg, [3 2 1]);
%         
%         avgdev_pat_col_ac(submatcount_pat,:, :) = permute(data_dev_ac.avg, [3 2 1]);
        

        
        
        %And repeat for drug session
        
        
%         dipmatch = strcmp(drugsubj, DrugToFlip(:,1));
%         
%         
%         if nnz(dipmatch)==1
%             
%             nodetoflip = cell2mat(DrugToFlip(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 flipndtrl = -(D2(nodetoflip, :, stdtrl));
%                 
%                 D2(nodetoflip, :, stdtrl) = flipndtrl;
%                 
%                 %D2(nodetoflip, :, devtrl) = -(D2(nodetoflip, :, devtrl));
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     flipndtrl = -(D2(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     D2(nodetoflip(nodeflip), :, stdtrl) = flipndtrl;
%                     
%                     %D2(nodetoflip(nodeflip), :, devtrl) = -(D2(nodetoflip(nodeflip), :, devtrl));
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
%         
%         
%         %Deviants
%         
%         dipmatch = strcmp(drugsubj, DrugToFlip_Dev(:,1));
%         
%         
%         if nnz(dipmatch)==1
%             
%             nodetoflip = cell2mat(DrugToFlip_Dev(dipmatch==1, 2));
%             
%             
%             if numel(nodetoflip)==1
%                 
%                 
%                 %D2(nodetoflip, :, stdtrl) = -(D2(nodetoflip, :, stdtrl));
%                 
%                 flipndtrl = -(D2(nodetoflip, :, devtrl));
%                 
%                 D2(nodetoflip, :, devtrl) = flipndtrl;
%                 
%                 
%             else
%                 
%                 
%                 for nodeflip = 1:numel(nodetoflip)
%                     
%                     
%                     %D2(nodetoflip(nodeflip), :, stdtrl) = -(D2(nodetoflip(nodeflip), :, stdtrl));
%                     
%                     flipndtrl = -(D2(nodetoflip(nodeflip), :, devtrl));
%                     
%                     D2(nodetoflip(nodeflip), :, devtrl) = flipndtrl;
%                     
%                     
%                 end
%                 
%                 
%             end
%             
%             
%         end
        
        
        
        %Other processing, including baseline correction and smoothing
        
        timelock=D2.fttimelock;
        
        
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
        data_std = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_d_col(submatcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_d_col(submatcount,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_d_col(submatcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev_d{submatcount_pat,1}=data_dev;
        
        avgstd_d{submatcount_pat,1}=data_std;
        
        avgstd3_d{submatcount_pat,1}=data_std3;
        
        
        
        % Auto-correction
        
%         timelockac = D2.fttimelock;
%         
%         for source = 1:nsources
%         
%         if timelockac.trial(devtrl, source, 85) > 0
%             
%             timelockac.trial(devtrl, source, :) = -(timelockac.trial(devtrl, source, :));
%             
%         end
%         
%         
%         if timelockac.trial(stdtrl, source, 85) > 0
%             
%             timelockac.trial(stdtrl, source, :) = -(timelockac.trial(stdtrl, source, :));
%             
%         end
%         
%         end
%         
%         %Baseline correct
%         cfg = [];
%         cfg.baseline = [-0.100 -0.001];
%         [baselinecorr_ac] = ft_timelockbaseline(cfg, timelockac);
%         
%         
%         %Std
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 7;
%         data_std_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         %Dev
%         cfg = [];
%         cfg.trials = timelock.trialinfo == 1;
%         data_dev_ac = ft_redefinetrial(cfg, baselinecorr_ac);
%         
%         
%         %Smoothing
%         
%         for source = 1:nsources
%             
%             %Std
%             
%             data_std_ac.avg(source,:) = smooth(data_std_ac.avg(source,:), 20,'moving');
%             
%             %Dev
%             
%             data_dev_ac.avg(source,:) = smooth(data_dev_ac.avg(source,:),20,'moving');
%             
%         end
%         
%         
%         %Combine into 3D matrix
%         
%         avgstd_d_pat_col_ac(submatcount_pat,:, :) = permute(data_std_ac.avg, [3 2 1]);
%         
%         avgdev_d_pat_col_ac(submatcount_pat,:, :) = permute(data_dev_ac.avg, [3 2 1]);
%         
        
        
    end
    
    
end


%%=========================================================================
%% Find matching data in MRS table
%%=========================================================================

% Remember to include only those with both MEG and MRS - likely some will
% not have MRS


MRStable_cons = MRStable(1:20,:);

MRStable_conids_lc = lower(MRStable_cons.Participant);



%Now those w/ both sess

match = 0;

ncons_bothsess = size(avgstd,1);

cons_bothsess_keep = zeros(1, ncons_bothsess);



for isub = 1:ncons_bothsess
    
    
    % subj
    
    subj = allcons_bothsess{isub};
    
    SubjMatch=strcmp(MRStable_conids_lc, [subj]);
    

    if ~isempty(SubjMatch) && ~isnan(MRStable_cons.gaba_ifg_corr(SubjMatch))
        
        
    cons_bothsess_keep(isub) = 1;
    
    match = match + 1;
    
    
    ExtTable_gabaifg_bothsess_cons(match,:) = MRStable_cons.gaba_ifg(SubjMatch==1);
    
    ExtTable_gabaifgcor_bothsess_cons(match,:) = MRStable_cons.gaba_ifg_corr(SubjMatch==1);
    
    
    end
    
    
end


cons_bothsess_keep = cons_bothsess_keep';



% Remove subjects without MRS

%Grand averages

cfg=[];


%Plac

% grandavgdev_pat  = ft_timelockgrandaverage(cfg, avgdev_pat{pats_bothsess_keep==1});
% 
% grandavgstd_pat  = ft_timelockgrandaverage(cfg, avgstd_pat{pats_bothsess_keep==1});
% 
% 
% %Drug
% 
% grandavgdev_d_pat  = ft_timelockgrandaverage(cfg, avgdev_d_pat{pats_bothsess_keep==1});
% 
% grandavgstd_d_pat  = ft_timelockgrandaverage(cfg, avgstd_d_pat{pats_bothsess_keep==1});
% 



% And collapsed data

%Non-corrected

avgstd_col_keep = avgstd_col(cons_bothsess_keep==1,:,:);

avgstd3_col_keep = avgstd3_col(cons_bothsess_keep==1,:,:);

avgdev_col_keep = avgdev_col(cons_bothsess_keep==1,:,:);


avgstd_d_col_keep = avgstd_d_col(cons_bothsess_keep==1,:,:);

avgstd3_d_col_keep = avgstd3_d_col(cons_bothsess_keep==1,:,:);

avgdev_d_col_keep = avgdev_d_col(cons_bothsess_keep==1,:,:);



%%=========================================================================
%% First rep6
%%=========================================================================

% Just need the scatters I think.. no int as well

MMNpall_cor = [];

MMNrall_cor = [];

MMNpall_cor_coeff = [];

MMNpall_lm_cor_int = [];


for source = 1:nsources
   
    
    %Placebo Cons
    
    diff = avgstd_col_keep(:,:,source)-avgdev_col_keep(:,:,source);
    
    
    %Placebo - Pats
    
    diff_pat = avgstd_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
    
    
       
    
    % Mean calculations
    
    m1difftemp=squeeze(nanmean(diff(:,MMNs:MMNf),2));
    
    m1difftemp_pat=squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    %m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %X = [ExtTable_gabaifg_bothsess - mean(ExtTable_gabaifg_bothsess)];

    [b, stats] = robustfit(ExtTable_gabaifg_bothsess_cons, m1difftemp);
    
    [b_pat, stats_pat] = robustfit(ExtTable_gabaifg_bothsess, m1difftemp_pat);
    
    
    MMNrall(source,1) = b(2,1);
    
    MMNpall(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifg_bothsess_cons, m1difftemp);
    
    MMNpall_coeff(source,1) = p(1,2);

    
    
    figure(17), subplot(2,3,source)
    
    
    hold on
    
    
    figure(17), scatter(ExtTable_gabaifg_bothsess_cons, m1difftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    figure(17), scatter(ExtTable_gabaifg_bothsess, m1difftemp_pat, 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifg_bothsess_cons;
    
    plot(x,b(1)+b(2)*x,'b', 'LineWidth', 2)
    
    
    x_pat = ExtTable_gabaifg_bothsess;
    
    plot(x_pat,b_pat(1)+b_pat(2)*x_pat,'r', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN PLA'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    
    % !! Repeat for corrected calculations
    
    

    [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess_cons, m1difftemp);

    
    [b_pat, stats_pat] = robustfit(ExtTable_gabaifgcor_bothsess, m1difftemp_pat);

    
    MMNrall_cor(source,1) = b(2,1);
    
    MMNpall_cor(source,1) = stats.p(2,1);
    
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifgcor_bothsess_cons, m1difftemp);
    
    MMNpall_cor_coeff(source,1) = p(1,2);
    
    
    %And lets do their interation..
    
    diaggrp = zeros(length(ExtTable_gabaifgcor_bothsess)+length(ExtTable_gabaifgcor_bothsess_cons),1);
    
    diaggrp(1:length(ExtTable_gabaifgcor_bothsess),1) = -1;
    
    diaggrp(length(ExtTable_gabaifgcor_bothsess)+1:end,1) = 1;
    
    
    temptbl = table(cat(1, ExtTable_gabaifgcor_bothsess, ExtTable_gabaifgcor_bothsess_cons), cat(1, m1difftemp_pat, m1difftemp), categorical(diaggrp), 'VariableNames',{'GABA','Diff','Group'}); 
    
    lm = fitlm(temptbl, 'Diff~GABA+Group'); 
    
    
    intterms = 'GABA*Group';
    lmint = addTerms(lm,intterms);
    
    
    MMNpall_lm_cor_int(source,1) = lmint.Coefficients.pValue(4);
    

    
    figure(18), subplot(2,3,source)
    
    
    hold on
    
    
    figure(18), scatter(ExtTable_gabaifgcor_bothsess_cons, m1difftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    figure(18), scatter(ExtTable_gabaifgcor_bothsess, m1difftemp_pat, 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifgcor_bothsess_cons;
    
    plot(x,b(1)+b(2)*x,'b', 'LineWidth', 2)
    
    
    x_pat = ExtTable_gabaifgcor_bothsess;
    
    plot(x_pat,b_pat(1)+b_pat(2)*x_pat,'r', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN PLA'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
end



%% Rep3

MMNpall_cor_rep3 = [];

MMNpall_cor_rep3_coeff = [];

MMNrall_cor_rep3 = [];


for source = 1:nsources
   
    
    %Placebo - Cons
    
    diff = avgstd3_col_keep(:,:,source)-avgdev_col_keep(:,:,source);
    
    
    %Placebo - Pats
    
    diff_pat = avgstd3_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
 
    
    
    % Mean calculations
    
    m1difftemp=squeeze(nanmean(diff(:,MMNs:MMNf),2));
    
    m1difftemp_pat=squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    
    %m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %X = [ExtTable_gabaifg_bothsess - mean(ExtTable_gabaifg_bothsess)];

    [b, stats] = robustfit(ExtTable_gabaifg_bothsess_cons, m1difftemp);
    
    [b_pat, stats_pat] = robustfit(ExtTable_gabaifg_bothsess, m1difftemp_pat);
    
    
    MMNrall_rep3(source,1) = b(2,1);
    
    MMNpall_rep3(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifg_bothsess_cons, m1difftemp);
    
    MMNpall_rep3_coeff(source,1) = p(1,2);

    
    
    figure(19), subplot(2,3,source)
    
    
    hold on
    
    
    figure(19), scatter(ExtTable_gabaifg_bothsess_cons, m1difftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    figure(19), scatter(ExtTable_gabaifg_bothsess, m1difftemp_pat, 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifg_bothsess_cons;
    
    plot(x,b(1)+b(2)*x,'b', 'LineWidth', 2)
    
    
    x_pat = ExtTable_gabaifg_bothsess;
    
    plot(x_pat,b_pat(1)+b_pat(2)*x_pat,'r', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN3 PLA'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    
    % !! Repeat for corrected calculations
    
    

    [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess_cons, m1difftemp);
    
    [b_pat, stats_pat] = robustfit(ExtTable_gabaifgcor_bothsess, m1difftemp_pat);
    
    
    MMNrall_cor_rep3(source,1) = b(2,1);
    
    MMNpall_cor_rep3(source,1) = stats.p(2,1);
    
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifgcor_bothsess_cons, m1difftemp);
    
    MMNpall_cor_rep3_coeff(source,1) = p(1,2);
    

    
    figure(20), subplot(2,3,source)
    
    
    hold on
    
    
    figure(20), scatter(ExtTable_gabaifgcor_bothsess_cons, m1difftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    figure(20), scatter(ExtTable_gabaifgcor_bothsess, m1difftemp_pat, 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifgcor_bothsess_cons;
    
    plot(x,b(1)+b(2)*x,'b', 'LineWidth', 2)
    
    
    x_pat = ExtTable_gabaifgcor_bothsess;
    
    plot(x_pat,b_pat(1)+b_pat(2)*x_pat,'r', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN3 PLA'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
   
    
end




% Save output

%Tables

%Rep6
[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor;
MMNstats.r = MMNrall_cor;

MMNstats.corrp = MMNpall_cor_coeff;

MMNstats.lmint = MMNpall_lm_cor_int;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN_PLA_PatsandCons_fullsample_stats.mat'], 'MMNstats');


%Rep3

[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor_rep3);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor_rep3;
MMNstats.r = MMNrall_cor_rep3;

MMNstats.corrp = MMNpall_cor_rep3_coeff;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN3_PLA_Cons_fullsample_stats.mat'], 'MMNstats');

   

%Save Figures



saveas(figure(17), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'meanMMN_PLA_PatsandCon' '_fullsample.tif']);


saveas(figure(17), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_PLA_PatsandCon' '_fullsample.tif']);


savefig(figure(18), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_PLA_PatsandCon' '_fullsample.fig']);


saveas(figure(19), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'meanMMN3_PLA_PatsandCon' '_fullsample.tif']);


saveas(figure(20), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_PLA_PatsandCon' '_fullsample.tif']);


savefig(figure(20), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_PLA_PatsandCon' '_fullsample.fig']);



%%=========================================================================
%% Now MRS assoc with drug difference
%%=========================================================================


MMNrall = [];

MMNpall = [];

MMNpall_coeff = [];


MMNrall_cor = [];

MMNpall_cor = [];

MMNpall_cor_coeff = [];



%% First rep6

for source = 1:nsources
   
    
    %Controls
    %Placebo
    
    diff = avgstd_col_keep(:,:,source)-avgdev_col_keep(:,:,source);
    
    
    %Drug
    
    diff_d = avgstd_d_col_keep(:,:,source)-avgdev_d_col_keep(:,:,source);   
       
    
    %PLA - MEM
    
    diff_drugdiff = diff - diff_d;
    
    
    %Patients
    %Placebo
    
    diff_pat = avgstd_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
   
    
    
    %Drug
    
    diff_d_pat = avgstd_d_pat_col_keep(:,:,source)-avgdev_d_pat_col_keep(:,:,source);   
    
    
    
    %PLA - MEM
    
    diff_drugdiff_pat = diff_pat - diff_d_pat;
        
    
    
    
    % Mean calculations
    
    m1drugdifftemp=squeeze(nanmean(diff_drugdiff(:,MMNs:MMNf),2));
    
    m1drugdifftemp_pat=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
    %m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %X = [ExtTable_gabaifg_bothsess - mean(ExtTable_gabaifg_bothsess)];

    [b, stats] = robustfit(ExtTable_gabaifg_bothsess_cons, m1drugdifftemp);
    
    [b_pat, stats_pat] = robustfit(ExtTable_gabaifg_bothsess, m1drugdifftemp_pat);
    
    
    MMNrall(source,1) = b(2,1);
    
    MMNpall(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifg_bothsess_cons, m1drugdifftemp);
    
    MMNpall_coeff(source,1) = p(1,2);

    
    
    figure(21), subplot(2,3,source)
    
    
    hold on
    
    
    figure(21), scatter(ExtTable_gabaifg_bothsess_cons, m1drugdifftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    figure(21), scatter(ExtTable_gabaifg_bothsess, m1drugdifftemp_pat, 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifg_bothsess_cons;
    
    plot(x,b(1)+b(2)*x,'b', 'LineWidth', 2)
    
    
    x_pat = ExtTable_gabaifg_bothsess;
    
    plot(x_pat,b_pat(1)+b_pat(2)*x_pat,'r', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN DrugDiff'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    
    % !! Repeat for corrected calculations
    
    

    [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess_cons, m1drugdifftemp);

    
    [b_pat, stats_pat] = robustfit(ExtTable_gabaifgcor_bothsess, m1drugdifftemp_pat);

    
    MMNrall_cor(source,1) = b(2,1);
    
    MMNpall_cor(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [r,p,rlo,rup] = corrcoef(ExtTable_gabaifgcor_bothsess_cons, m1drugdifftemp);
    
    MMNpall_cor_coeff(source,1) = p(1,2);

    

    
    figure(22), subplot(2,3,source)
    
    
    hold on
    
    
    figure(22), scatter(ExtTable_gabaifgcor_bothsess_cons, m1drugdifftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    figure(22), scatter(ExtTable_gabaifgcor_bothsess, m1drugdifftemp_pat, 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifgcor_bothsess_cons;
    
    plot(x,b(1)+b(2)*x,'b', 'LineWidth', 2)
    
    
    x_pat = ExtTable_gabaifgcor_bothsess;
    
    plot(x_pat,b_pat(1)+b_pat(2)*x_pat,'r', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean DrugDiff MMN'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
end



%% Rep3

MMNpall_rep3 = [];

MMNpall_rep3_coeff = [];

MMNrall_rep3 = [];


MMNpall_cor_rep3 = [];

MMNpall_cor_rep3_coeff = [];

MMNrall_cor_rep3 = [];



for source = 1:nsources
   
    %Controls
    %Placebo
    
    diff = avgstd3_col_keep(:,:,source)-avgdev_col_keep(:,:,source);
    
    
    %Drug
    
    diff_d = avgstd3_d_col_keep(:,:,source)-avgdev_d_col_keep(:,:,source);   
    
    
    %PLA - MEM
    
    diff_drugdiff = diff - diff_d;
    
    
    %Patients
    %Placebo
    
    diff_pat = avgstd3_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
   
    
    %Drug
    
    diff_d_pat = avgstd3_d_pat_col_keep(:,:,source)-avgdev_d_pat_col_keep(:,:,source);   
      
    
    %PLA - MEM
    
    diff_drugdiff_pat = diff_pat - diff_d_pat;
        
   
    
    % Mean calculations
    
    m1drugdifftemp=squeeze(nanmean(diff_drugdiff(:,MMNs:MMNf),2));
    
    m1drugdifftemp_pat=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
    %m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %X = [ExtTable_gabaifg_bothsess - mean(ExtTable_gabaifg_bothsess)];

    [b, stats] = robustfit(ExtTable_gabaifg_bothsess_cons, m1drugdifftemp);
    
    [b_pat, stats_pat] = robustfit(ExtTable_gabaifg_bothsess, m1drugdifftemp_pat);
    
    
    MMNrall_rep3(source,1) = b(2,1);
    
    MMNpall_rep3(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifg_bothsess_cons, m1drugdifftemp);
    
    MMNpall_rep3_coeff(source,1) = p(1,2);

    
    
    figure(23), subplot(2,3,source)
    
    
    hold on
    
    
    figure(23), scatter(ExtTable_gabaifg_bothsess_cons, m1drugdifftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    figure(23), scatter(ExtTable_gabaifg_bothsess, m1drugdifftemp_pat, 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifg_bothsess_cons;
    
    plot(x,b(1)+b(2)*x,'b', 'LineWidth', 2)
    
    
    x_pat = ExtTable_gabaifg_bothsess;
    
    plot(x_pat,b_pat(1)+b_pat(2)*x_pat,'r', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN3 DrugDiff'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    
    % !! Repeat for corrected calculations
    
    

    [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess_cons, m1drugdifftemp);

    
    [b_pat, stats_pat] = robustfit(ExtTable_gabaifgcor_bothsess, m1drugdifftemp_pat);

    
    MMNrall_cor_rep3(source,1) = b(2,1);
    
    MMNpall_cor_rep3(source,1) = stats.p(2,1);
    
    
    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifgcor_bothsess_cons, m1drugdifftemp);
    
    MMNpall_cor_rep3_coeff(source,1) = p(1,2);
    
    
    
    %And lets do their interation..
    
    diaggrp = zeros(length(ExtTable_gabaifgcor_bothsess)+length(ExtTable_gabaifgcor_bothsess_cons),1);
    
    diaggrp(1:length(ExtTable_gabaifgcor_bothsess),1) = -1;
    
    diaggrp(length(ExtTable_gabaifgcor_bothsess)+1:end,1) = 1;
    
    
    temptbl = table(cat(1, ExtTable_gabaifgcor_bothsess, ExtTable_gabaifgcor_bothsess_cons), cat(1, m1drugdifftemp_pat, m1drugdifftemp), categorical(diaggrp), 'VariableNames',{'GABA','DrugDiff','Group'}); 
    
    lm = fitlm(temptbl, 'DrugDiff~GABA+Group'); 
    
    
    intterms = 'GABA*Group';
    lmint = addTerms(lm,intterms);
    
    
    MMNpall_lm_rep3_int(source,1) = lmint.Coefficients.pValue(4);
    
    
    figure(24), subplot(2,3,source)
    
    
    hold on
    
    
    figure(24), scatter(ExtTable_gabaifgcor_bothsess_cons, m1drugdifftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    figure(24), scatter(ExtTable_gabaifgcor_bothsess, m1drugdifftemp_pat, 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifgcor_bothsess_cons;
    
    plot(x,b(1)+b(2)*x,'b', 'LineWidth', 2)
    
    
    x_pat = ExtTable_gabaifgcor_bothsess;
    
    plot(x_pat,b_pat(1)+b_pat(2)*x_pat,'r', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean DrugDiff MMN3'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
end



% Save output

%Tables

%Rep6
[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor;
MMNstats.r = MMNrall_cor;

MMNstats.corrp = MMNpall_cor_coeff;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN_DrugDiff_PatsandCons_fullsample_stats.mat'], 'MMNstats');


%Rep3

[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor_rep3);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor_rep3;
MMNstats.r = MMNrall_cor_rep3;

MMNstats.corrp = MMNpall_cor_rep3_coeff;


MMNstats.lmintp = MMNpall_lm_rep3_int;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN3_DrugDiff_PatsandCons_fullsample_stats.mat'], 'MMNstats');

   

%Save Figures


saveas(figure(21), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'meanMMN_DrugDiff_PatsandCons' '_fullsample.tif']);


saveas(figure(22), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_DrugDiff_PatsandCons' '_fullsample.tif']);


savefig(figure(22), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_DrugDiff_PatsandCons' '_fullsample.fig']);


saveas(figure(23), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'meanMMN3_DrugDiff_PatsandCons' '_fullsample.tif']);


saveas(figure(24), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_DrugDiff_PatsandCons' '_fullsample.tif']);


savefig(figure(24), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_DrugDiff_PatsandCons' '_fullsample.fig']);


%%=========================================================================
%% Diagnostic group information
%%=========================================================================


%Need to do placebo and drug differences
%For both rep3 and rep6


%% Placebo


% Just need the scatters I think..


MMNpall_cor = [];

MMNrall_cor = [];

MMNpall_lm_cor = [];


X = cat(2, ExtTable_gabaifgcor_bothsess, cell2mat(extdiaginfo(:,1)));



% First rep6

for source = 1:nsources
   
    

    %Placebo - Pats
    
    diff_pat = avgstd_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
    
    
 
    % Mean calculations
    
    
    m1difftemp_pat=squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    
    %Do together controlling for diag info
    
    
    % !! Straight to corrected calculations
        

    [b, stats] = robustfit(X, m1difftemp_pat);


    
    MMNrall_cor(source,1) = b(2,1);
    
    MMNpall_cor(source,1) = stats.p(2,1);
    
    
        
    %And also with norm reg model
    
    temptbl = table(ExtTable_gabaifgcor_bothsess, m1difftemp_pat, categorical(cell2mat(extdiaginfo)), 'VariableNames',{'GABA','Diff','Group'}); 
    
    lm = fitlm(temptbl, 'Diff~GABA+Group');
    
    
    MMNpall_lm_cor(source,1) = lm.Coefficients.pValue(2);

    
    
    %Now separately for plotting purposes
    
    %bv == 1
    [b_bv, ~] = robustfit(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==1), m1difftemp_pat(cell2mat(extdiaginfo(:,1))==1));
    
    %psp == 2
    [b_psp, ~] = robustfit(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==2), m1difftemp_pat(cell2mat(extdiaginfo(:,1))==2));

    
    figure(25), subplot(2,3,source)
    
    
    hold on
    
    
    %bv = blue
    figure(25), scatter(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==1), m1difftemp_pat(cell2mat(extdiaginfo(:,1))==1), 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    %psp = red
    figure(25), scatter(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==2), m1difftemp_pat(cell2mat(extdiaginfo(:,1))==2), 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifgcor_bothsess;
    
    
    %bv
    
    plot(x,b_bv(1)+b_bv(2)*x,'b', 'LineWidth', 2)
    
    
    %psp
    
    plot(x,b_psp(1)+b_psp(2)*x,'r', 'LineWidth', 2)
    
    
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN PLA'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source) ' ' 'bv/PSP'], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
end



%% Rep3

MMNpall_cor_rep3 = [];

MMNrall_cor_rep3 = [];


for source = 1:nsources
   
    

    %Placebo - Pats
    
    diff_pat = avgstd3_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
    
    
 
    % Mean calculations
    
    
    m1difftemp_pat=squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    
    %Do together controlling for diag info
    
    
    % !! Straight to corrected calculations
    
   
    %Set up DM
    
    X = cat(2, ExtTable_gabaifgcor_bothsess, cell2mat(extdiaginfo));

    [b, stats] = robustfit(X, m1difftemp_pat);

   
    
    MMNrall_cor_rep3(source,1) = b(2,1);
    
    MMNpall_cor_rep3(source,1) = stats.p(2,1);
    
    

    
    

    %Now separately for plotting purposes
    
    %bv == 1
    [b_bv, ~] = robustfit(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==1), m1difftemp_pat(cell2mat(extdiaginfo(:,1))==1));
    
    %psp == 2
    [b_psp, ~] = robustfit(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==2), m1difftemp_pat(cell2mat(extdiaginfo(:,1))==2));

    
    figure(26), subplot(2,3,source)
    
    
    hold on
    
    
    %bv = blue
    figure(26), scatter(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==1), m1difftemp_pat(cell2mat(extdiaginfo(:,1))==1), 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    %psp = red
    figure(26), scatter(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==2), m1difftemp_pat(cell2mat(extdiaginfo(:,1))==2), 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifgcor_bothsess;
    
    
    %bv
    
    plot(x,b_bv(1)+b_bv(2)*x,'b', 'LineWidth', 2)
    
    
    %psp
    
    plot(x,b_psp(1)+b_psp(2)*x,'r', 'LineWidth', 2)
    
    
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN3 PLA'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source) ' ' 'bv/PSP'], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
end


% Save output

%Tables

%Rep6
[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor;
MMNstats.r = MMNrall_cor;

MMNstats.plm = MMNpall_lm_cor;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN_PLA_Pats_wDIAG_fullsample_stats.mat'], 'MMNstats');


%Rep3

[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor_rep3);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor_rep3;
MMNstats.r = MMNrall_cor_rep3;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN3_PLA_Pats_wDIAG_fullsample_stats.mat'], 'MMNstats');

   

%Save Figures


saveas(figure(25), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_PLA_Pats_wDIAG' '_fullsample.tif']);


savefig(figure(25), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_PLA_Pats_wDIAG' '_fullsample.fig']);


saveas(figure(26), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_PLA_Pats_wDIAG' '_fullsample.tif']);


savefig(figure(26), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_PLA_Pats_wDIAG' '_fullsample.fig']);



%%=========================================================================
%% Now MRS assoc with drug difference
%%=========================================================================


%% First rep6

MMNrall_cor = [];

MMNpall_cor = [];


for source = 1:nsources
   
    
   
    %Patients
    %Placebo
    
    diff_pat = avgstd_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
   
    
    
    %Drug
    
    diff_d_pat = avgstd_d_pat_col_keep(:,:,source)-avgdev_d_pat_col_keep(:,:,source);   
    
    
    
    %PLA - MEM
    
    diff_drugdiff_pat = diff_pat - diff_d_pat;
        
    
    
    
    % Mean calculations
    
    
    m1drugdifftemp_pat=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
    
    
    % !! Straight to corrected calculations
    
    
    %All together 

    [b, stats] = robustfit(X, m1drugdifftemp_pat);

   
    
    MMNrall_cor(source,1) = b(2,1);
    
    MMNpall_cor(source,1) = stats.p(2,1);
    
    
    
    %Now separately for plotting purposes
    
    %bv == 1
    [b_bv, ~] = robustfit(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==1), m1drugdifftemp_pat(cell2mat(extdiaginfo(:,1))==1));
    
    %psp == 2
    [b_psp, ~] = robustfit(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==2), m1drugdifftemp_pat(cell2mat(extdiaginfo(:,1))==2));

    
    figure(27), subplot(2,3,source)
    
    
    hold on
    
    
    %bv = blue
    figure(27), scatter(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==1), m1drugdifftemp_pat(cell2mat(extdiaginfo(:,1))==1), 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    %psp = red
    figure(27), scatter(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==2), m1drugdifftemp_pat(cell2mat(extdiaginfo(:,1))==2), 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifgcor_bothsess;
    
 
    
    %bv
    
    plot(x,b_bv(1)+b_bv(2)*x,'b', 'LineWidth', 2)
    
    
    %psp
    
    plot(x,b_psp(1)+b_psp(2)*x,'r', 'LineWidth', 2)
    
    
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN DrugDiff'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source) ' ' 'bv/PSP'], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
end



%% Rep3



MMNpall_cor_rep3 = [];

MMNrall_cor_rep3 = [];

MMNpall_lm_cor_rep3 = [];


for source = 1:nsources
   
    
   
    %Patients
    %Placebo
    
    diff_pat = avgstd3_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
   
    
    
    %Drug
    
    diff_d_pat = avgstd3_d_pat_col_keep(:,:,source)-avgdev_d_pat_col_keep(:,:,source);   
    
    
    
    %PLA - MEM
    
    diff_drugdiff_pat = diff_pat - diff_d_pat;
        
    
    
    
    % Mean calculations
    
    
    m1drugdifftemp_pat=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
    
    
    % !! Straight to corrected calculations
    
    
    %All together 

    X = cat(2, ExtTable_gabaifgcor_bothsess, cell2mat(extdiaginfo));
    
    [b, stats] = robustfit(X, m1drugdifftemp_pat);

   
    
    MMNrall_cor_rep3(source,1) = b(2,1);
    
    MMNpall_cor_rep3(source,1) = stats.p(2,1);
    
    
    %And also with norm reg model
    
    temptbl = table(ExtTable_gabaifgcor_bothsess, m1drugdifftemp_pat, categorical(cell2mat(extdiaginfo)), 'VariableNames',{'GABA','DrugDiff','Group'}); 
    
    lm = fitlm(temptbl, 'DrugDiff~GABA+Group');
    
    
    MMNpall_lm_cor_rep3(source,1) = lm.Coefficients.pValue(2);
    
    
    
    %Now separately for plotting purposes
    
    %bv == 1
    [b_bv, ~] = robustfit(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==1), m1drugdifftemp_pat(cell2mat(extdiaginfo(:,1))==1));
    
    %psp == 2
    [b_psp, ~] = robustfit(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==2), m1drugdifftemp_pat(cell2mat(extdiaginfo(:,1))==2));

    
    figure(28), subplot(2,3,source)
    
    
    hold on
    
    
    %bv = blue
    figure(28), scatter(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==1), m1drugdifftemp_pat(cell2mat(extdiaginfo(:,1))==1), 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    %psp = red
    figure(28), scatter(ExtTable_gabaifgcor_bothsess(cell2mat(extdiaginfo(:,1))==2), m1drugdifftemp_pat(cell2mat(extdiaginfo(:,1))==2), 'filled','MarkerFaceColor','r','MarkerEdgeColor','k')
    
    
    x = ExtTable_gabaifgcor_bothsess;
    
 
    
    %bv
    
    plot(x,b_bv(1)+b_bv(2)*x,'b', 'LineWidth', 2)
    
    
    %psp
    
    plot(x,b_psp(1)+b_psp(2)*x,'r', 'LineWidth', 2)
    
    
    
    set(gca, 'Fontsize',12); ylabel('Mean MMN3 DrugDiff'); xlabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source) ' ' 'bv/PSP'], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
end


% Save output

%Tables

%Rep6
[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor;
MMNstats.r = MMNrall_cor;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN_DrugDiff_Pats_wDIAG_fullsample_stats.mat'], 'MMNstats');



%Rep3

[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor_rep3);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor_rep3;
MMNstats.r = MMNrall_cor_rep3;

MMNstats.plm = MMNpall_lm_cor_rep3;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN3_DrugDiff_Pats_wDIAG_fullsample_stats.mat'], 'MMNstats');

   

%Save Figures


saveas(figure(27), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_DrugDiff_Pats_wDIAG' '_fullsample.tif']);


savefig(figure(27), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN_DrugDiff_Pats_wDIAG' '_fullsample.fig']);


saveas(figure(28), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_DrugDiff_Pats_wDIAG' '_fullsample.tif']);


savefig(figure(28), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN3_DrugDifg_Pats_wDIAG' '_fullsample.fig']);



%% Not finished yet - further control analyses


%% Rep1

for source = 1:nsources
   
    
    %Placebo
    
    diff_pat = avgstd1_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
    
    m1diff_pat=squeeze(nanmean(diff_pat,1));
    s1diff_pat=squeeze(nanstd(diff_pat)./sqrt(size(diff_pat,1)));
    
    
    %Drug
    
    diff_d_pat = avgstd1_d_pat_col_keep(:,:,source)-avgdev_d_pat_col_keep(:,:,source);   
    
    m1diff_d_pat = squeeze(nanmean(diff_d_pat,1));
    s1diff_d_pat = squeeze(nanstd(diff_d_pat)./sqrt(size(diff_d_pat,1)));
    
    
    %PLA - MEM
    
    diff_drugdiff_pat = diff_pat - diff_d_pat;
        
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_drugdiff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifg_bothsess);
%         
%         pall(1,t_win) = p(1,2);


        [b, stats] = robustfit(ExtTable_gabaifg_bothsess, diff_drugdiff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(29), subplot(2,3,source);
    
    hold on
    
    figure(29);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); legendflex({'MEM','PLA'},'fontsize',14)
    boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) ' MRSIFGgaba Assoc in MMN1 Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
   
    
    
    % Corrected values
    
    ball = [];
    pall = [];
    
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_drugdiff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifgcor_bothsess);
%         
%         pall(1,t_win) = p(1,2);

        [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, diff_drugdiff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        ball(1,t_win) = b(2,1);

        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(30), subplot(2,3,source);
    
    hold on
    
    figure(30);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); legendflex({'MEM','PLA'},'fontsize',14)
    boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    
    %Positive assoc
%     plot(p(ball>0),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    %Negative
%     plot(p(ball<0),'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)


    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) ' MRSIFGgabacor Assoc in MMN1 Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    % Mean calculations
    
    
    m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %X = [ExtTable_gabaifg_bothsess - mean(ExtTable_gabaifg_bothsess)];

    [b, stats] = robustfit(ExtTable_gabaifg_bothsess, m1drugdifftemp);
    
    
    MMNrall_rep1(source,1) = b(2,1);
    
    MMNpall_rep1(source,1) = stats.p(2,1);
    

    %Norm corr
    
    [~,p,~,~] = corrcoef(ExtTable_gabaifg_bothsess, m1drugdifftemp);
    
    MMNpall_rep1_coeff(source,1) = p(1,2);
    
    
    figure(31), subplot(2,3,source)
    
    
    hold on
    
    
    figure(31), scatter(ExtTable_gabaifg_bothsess, m1drugdifftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    
    set(gca, 'Fontsize',12); xlabel('Mean Drug Diff MMN1'); ylabel('gabaifg');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    % !! Repeat for corrected calculations
    
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifgcor_bothsess);
%     
%     
%     MMNrall_cor(source,1) = r(1,2);
%     
%     MMNpall_cor(source,1) = p(1,2);
%     


    [b, stats] = robustfit(ExtTable_gabaifgcor_bothsess, m1drugdifftemp);
    
    
    MMNrall_cor_rep1(source,1) = b(2,1);
    
    MMNpall_cor_rep1(source,1) = stats.p(2,1);
    
    
    [r,p,rlo,rup] = corrcoef(ExtTable_gabaifgcor_bothsess, m1drugdifftemp);
    
    
    MMNpall_cor_rep1_coeff(source,1) = p(1,2);
    
    MMNpall_cor_rep1_coeff_r2(source,1) = r(1,2)*r(1,2);
    
    MMNpall_cor_rep1_coeff_CI(source,1:2) = [rlo(1,2); rup(1,2)];
    
    

    figure(32), subplot(2,3,source)
    
    
    hold on
    
    
    figure(32), scatter(ExtTable_gabaifgcor_bothsess, m1drugdifftemp, 'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
    
    x = ExtTable_gabaifgcor_bothsess;
    
    plot(x,b(1)+b(2)*x,'k', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean Drug Diff MMN1 (PLA-MEM)'); xlabel('gaba rIFG');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
    
    %Save output table for LMM - do it because to prove sig is not just
    %rep3..
    
    lmm_druggabaint = [];
    
    lmm_druggabaint(1:length(ExtTable_gabaifgcor_bothsess), 1) = 1:length(ExtTable_gabaifgcor_bothsess);
    
    lmm_druggabaint(:,2) = ExtTable_gabaifgcor_bothsess;
    
    lmm_druggabaint(:,3) = squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,4) = squeeze(nanmean(diff_d_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,5) = m1drugdifftemp;
    
    lmm_druggabaint(:,6) = cell2mat(extdiaginfo);
    
    
    lmm_druggabaint_table = table(lmm_druggabaint(:,1), lmm_druggabaint(:,2), lmm_druggabaint(:,3), lmm_druggabaint(:,4), lmm_druggabaint(:,5), lmm_druggabaint(:,6));
    
    
    lmm_druggabaint_table.Properties.VariableNames = {'Subject', 'GABA', 'PLA', 'MEM', 'Diff', 'Diag'};
    
    
    writetable(lmm_druggabaint_table, [FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN1_DrugDiff_Pats_fullsample_stats_' 'LMMtableforR_' grandavgdev_pat.label{source} '.txt'], 'Delimiter','tab');
    
    
    
end


%Save output one last time


%Table

%Rep1

[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor_rep1);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor_rep1;

MMNstats.b = MMNrall_cor_rep1;

MMNstats.corrp = MMNpall_cor_rep1_coeff;

MMNstats.r2 = MMNpall_cor_rep1_coeff_r2;

MMNstats.CI = MMNpall_cor_rep1_coeff_CI;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSIFGgabacorassoc_' 'meanMMN1_DrugDiff1_Pats' '_fullsample_stats.mat'], 'MMNstats');


%Figures 


saveas(figure(29), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'DrugDiff1_Pats' '_fullsample.tif']);


saveas(figure(30), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'DrugDiff1_Pats' '_fullsample.tif']);


saveas(figure(31), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabaassoc' '_' 'meanMMN1_DrugDiff_Pats' '_fullsample.tif']);


saveas(figure(32), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN1_DrugDiff_Pats' '_fullsample.tif']);


savefig(figure(32), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSIFGgabacorassoc' '_' 'meanMMN1_DrugDiff_Pats' '_fullsample.fig']);



%% Rep3 - OCC GABA

MMNpall_cor_rep3 = [];

MMNrall_cor_rep3 = [];


for source = 1:nsources
   
    
    %Placebo
    
    diff_pat = avgstd3_pat_col_keep(:,:,source)-avgdev_pat_col_keep(:,:,source);
    
    m1diff_pat=squeeze(nanmean(diff_pat,1));
    s1diff_pat=squeeze(nanstd(diff_pat)./sqrt(size(diff_pat,1)));
    
    
    %Drug
    
    diff_d_pat = avgstd3_d_pat_col_keep(:,:,source)-avgdev_d_pat_col_keep(:,:,source);   
    
    m1diff_d_pat = squeeze(nanmean(diff_d_pat,1));
    s1diff_d_pat = squeeze(nanstd(diff_d_pat)./sqrt(size(diff_d_pat,1)));
    
    
    %PLA - MEM
    
    diff_drugdiff_pat = diff_pat - diff_d_pat;
        
    
    %Note, correlation analysis now for every time point..
    
    
    % Corrected values
    
    ball = [];
    pall = [];
    
    
    %Note, correlation analysis now for every time point..
    
    for t_win = 1:length(diff_drugdiff_pat(1,:))
        
%         [r, p, ~, ~] = corrcoef(diff_drugdiff_pat(:,t_win),ExtTable_gabaifgcor_bothsess);
%         
%         pall(1,t_win) = p(1,2);

        [b, stats] = robustfit(ExtTable_gabaocccor_bothsess, diff_drugdiff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        ball(1,t_win) = b(2,1);

        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(33), subplot(2,3,source);
    
    hold on
    
    figure(33);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); legendflex({'MEM','PLA'},'fontsize',14)
    boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    
    %Positive assoc
%     plot(p(ball>0),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    %Negative
%     plot(p(ball<0),'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)


    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([grandavgdev_pat.label(source) ' MRSOCCgabacor Assoc in MMN3 Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    % Mean calculations
    
    
    m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
    
    
    % !! Repeat for corrected calculations
    
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifgcor_bothsess);
%     
%     
%     MMNrall_cor(source,1) = r(1,2);
%     
%     MMNpall_cor(source,1) = p(1,2);
%     


    [b, stats] = robustfit(ExtTable_gabaocccor_bothsess, m1drugdifftemp);
    
    
    MMNrall_cor_rep3(source,1) = b(2,1);
    
    MMNpall_cor_rep3(source,1) = stats.p(2,1);
    
    
    [~,p,~,~] = corrcoef(ExtTable_gabaocccor_bothsess, m1drugdifftemp);
    
    MMNpall_cor_rep3_coeff(source,1) = p(1,2);
    
    

    figure(34), subplot(2,3,source)
    
    
    hold on
    
    
    figure(34), scatter(ExtTable_gabaocccor_bothsess, m1drugdifftemp, 'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
    
    x = ExtTable_gabaocccor_bothsess;
    
    plot(x,b(1)+b(2)*x,'k', 'LineWidth', 2)
    
    set(gca, 'Fontsize',12); ylabel('Mean Drug Diff MMN3 (PLA-MEM)'); xlabel('gaba OCC');
    
    
    title([grandavgdev_pat.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;
    
   
    %Save output table for LMM
    
    lmm_druggabaint = [];
    
    lmm_druggabaint(1:length(ExtTable_gabaocccor_bothsess), 1) = 1:length(ExtTable_gabaifg_bothsess);
    
    lmm_druggabaint(:,2) = ExtTable_gabaocccor_bothsess;
    
    lmm_druggabaint(:,3) = squeeze(nanmean(diff_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,4) = squeeze(nanmean(diff_d_pat(:,MMNs:MMNf),2));
    
    lmm_druggabaint(:,5) = m1drugdifftemp;
    
    
    lmm_druggabaint_table = table(lmm_druggabaint(:,1), lmm_druggabaint(:,2), lmm_druggabaint(:,3), lmm_druggabaint(:,4), lmm_druggabaint(:,5));
    
    
    lmm_druggabaint_table.Properties.VariableNames = {'Subject', 'GABAocc', 'PLA', 'MEM', 'Diff'};
    
    
    writetable(lmm_druggabaint_table, [FigOutDir '/' figoutprefix '_LFPs_' 'MRSOCCgabacorassoc_' 'meanMMN3_DrugDiff_Pats_fullsample_stats_' 'LMMtableforR_' grandavgdev_pat.label{source} '.txt'], 'Delimiter','tab');
    
    
end


%Save output


%Rep3

[h, ~, ~, adj_p]=fdr_bh(MMNpall_cor_rep3);

MMNstats = struct;
MMNstats.adj_p = adj_p;

MMNstats.praw = MMNpall_cor_rep3;
MMNstats.r = MMNrall_cor_rep3;

MMNstats.corrp = MMNpall_cor_rep3_coeff;


save([FigOutDir '/' figoutprefix '_LFPs_' 'MRSOCCgabacorassoc_' 'meanMMN_DrugDiff3_Pats' '_fullsample_stats.mat'], 'MMNstats');



saveas(figure(33), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSOCCgabacorassoc' '_' 'DrugDiff3_Pats' '_fullsample.tif']);


saveas(figure(34), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSOCCgabacorassoc' '_' 'meanMMN3_DrugDiff_Pats' '_fullsample.tif']);


savefig(figure(34), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'MRSOCCgabacorassoc' '_' 'meanMMN3_DrugDiff_Pats' '_fullsample.fig']);



%% Wow, Finished

sprintf('finito, output saved to %s \n', FigOutDir)

   
end