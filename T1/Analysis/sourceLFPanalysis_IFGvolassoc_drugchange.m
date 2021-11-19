function sourceLFPanalysis_IFGvolassoc_drugchange


%% Setup

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

figoutprefix = 'adv_ssst_newmaxf_fixICA_wfids_nodipcor_IFGarea_rawscores_newT1preproc';


%When removing that outlier
%figoutprefix = 'adv_ssst_newmaxf_fixICA_wfids_nodipcor_remo';



%% IFG info

ifgtable = '/imaging/rowe/Michelle/AlistairIFG/B_Data_preproc_0001/c1_Wmap_4849OP8_GMV_tGMVtable_wgender.csv';
IFGTable = readtable(ifgtable);




%% Other setup
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


%And when removing outliers

Pat_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Pat_D_nop10p19.txt';



fid = fopen(Pat_P);
PlacData_pat = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
PlacData_pat = PlacData_pat{1,1};

fid = fopen(Pat_D);
DrugData_pat = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
DrugData_pat = DrugData_pat{1,1};


% Can call it directly as it wont change
DiagData = readcell('/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/PatDiagInfo.txt');



%% Patients


submatcount_pat = 0;



for pplasubj = 1:length(PlacData_pat(:,1))
    
    
    filen = [preprocdir '/' PlacData_pat{pplasubj} '/' LFPprefix '_' preprocsteps '_' PlacData_pat{pplasubj} '_' 's1' '_' LFPpostfix];
    
    
    filen_d = [preprocdir '/' PlacData_pat{pplasubj} '/' LFPprefix '_' preprocsteps '_' PlacData_pat{pplasubj} '_' 's2' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount_pat = submatcount_pat + 1;
        
        
        allpats_bothsess{submatcount_pat,1} = PlacData_pat{pplasubj};
        
        
        %Also find diag info
        
        PatMatch=strcmp(DiagData(:,1), PlacData_pat{pplasubj});

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
        
        avgstd_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
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
        
        avgstd_d_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_d_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
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
        
        PatMatch=strcmp(DiagData(:,1), DrugData_pat{pdrugsubj});

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
        
        avgstd_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
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
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_d_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_d_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
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


%% Find matching data in IFG table
% Remember to include only those with both MEG and IFG data
% Ensures right order as well..


IFGtable_patids_lc = lower(IFGTable.ID);



%Now those w/ both sess

match = 0;

npats_bothsess = size(avgstd_pat,1);



for isub = 1:npats_bothsess
    
    
    % subj
    
    subj = allpats_bothsess{isub};
    
    SubjMatch=strcmp(IFGtable_patids_lc, [subj]);
    
    
    if nnz(SubjMatch)
        
        match = match + 1;
        
        pats_bothsess_keep(isub) = 1;
        
        ExtTable_IFG(match,:) = IFGTable.tGMV(SubjMatch==1);
        
        ExtTable_Age(match,:) = IFGTable.Age(SubjMatch==1);
        
        ExtTable_TIV(match,:) = IFGTable.TIV(SubjMatch==1);
        
    end
    
    
end


pats_bothsess_keep = pats_bothsess_keep';



% And collapsed data

avgstd_pat_col_keep = avgstd_pat_col(pats_bothsess_keep==1,:,:);

avgstd3_pat_col_keep = avgstd3_pat_col(pats_bothsess_keep==1,:,:);

avgdev_pat_col_keep = avgdev_pat_col(pats_bothsess_keep==1,:,:);


avgstd_d_pat_col_keep = avgstd_d_pat_col(pats_bothsess_keep==1,:,:);

avgstd3_d_pat_col_keep = avgstd3_d_pat_col(pats_bothsess_keep==1,:,:);

avgdev_d_pat_col_keep = avgdev_d_pat_col(pats_bothsess_keep==1,:,:);



% And lastly diagnostic info

extdiaginfo = extdiaginfo(pats_bothsess_keep==1,:);



%% Basic plotting and analyses - association with clin scores


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
%% Now IFG assoc with drug difference
%%=========================================================================


MMNrall = [];

MMNpall = [];


MMNrall_cor = [];

MMNpall_cor = [];


%% Rep3
%Straight to it

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


        [b, stats] = robustfit(ExtTable_IFG, diff_drugdiff_pat(:,t_win));
        
        pall(1,t_win) = stats.p(2,1);
        
        
    end

    k=max(m1diff_pat)+max(s1diff_pat);k=k*1.1;
    p=k*(pall<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    
    figure(1), subplot(2,3,source);
    
    hold on
    
    figure(1);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1diff_d_pat,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[0/255 0/255 255/255]); legendflex({'MEM','PLA'},'fontsize',14)
    boundedline([1:size(m1diff_d_pat,2)],m1diff_d_pat(1,:),s1diff_d_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12); xlabel('Time (ms)'); ylabel('Std - Dev')
    title([data_std.label(source) ' IFG Assoc in MMN3 Pat'], 'Fontsize',12); xlim([0 251]); %ylim([0 1]);%axis square
   
    
    
  
    % Mean calculations
    %And adjust for co-variates
    
    
    m1drugdifftemp=squeeze(nanmean(diff_drugdiff_pat(:,MMNs:MMNf),2));
    
%     [r, p, ~, ~] = corrcoef(m1drugdifftemp, ExtTable_gabaifg_bothsess);
%     
%     
%     MMNrall(source,1) = r(1,2);
%     
%     MMNpall(source,1) = p(1,2);


    %Put into table for lm
    
    %First residuals
    
    
    lmm_drugifgint_restable = table(ExtTable_Age, ExtTable_IFG, ExtTable_TIV);
    
    lmm_drugifgint_restable.Properties.VariableNames = {'Age', 'IFGvol', 'TIV'};

    
    lmmodel = fitlm(lmm_drugifgint_restable, 'IFGvol ~ TIV + Age');
    
    IFGres = lmmodel.Residuals.Raw;
    
    
    %Now actual analysis
    
    lmm_drugifgint_table = table(ExtTable_Age, ExtTable_IFG, ExtTable_TIV, m1drugdifftemp, IFGres, cell2mat(extdiaginfo));
  
    lmm_drugifgint_table.Properties.VariableNames = {'Age', 'IFGvol', 'TIV', 'Diff', 'IFGres', 'DiagGroup'};
    
    
    lmmodel = fitlm(lmm_drugifgint_table, 'Diff ~ IFGvol + TIV + Age');
    
    
    
    figure(2), subplot(2,3,source)
    
    
    hold on
    
    
    figure(2), scatter(IFGres, m1drugdifftemp, 'filled','MarkerFaceColor','b','MarkerEdgeColor','k')
    
    
    set(gca, 'Fontsize',12); xlabel('IFG'); ylabel('Mean Drug Diff MMN3');
    
    
    title([data_std.label(source)], 'FontSize', 14)
    
    
    ax=gca;
    
    ax.YColor = 'black';
    ax.XColor = 'black';
    ax.FontWeight = 'bold';
    ax.FontSize = 10;


    if source == 6 %Just contain it to RAUD for space purposes

        
        writetable(lmm_drugifgint_table, [FigOutDir '/' figoutprefix '_LFPs_' 'IFGvolassoc_' 'meanMMN3_DrugDiff_Pats_fullsample_stats_' 'LMMtableforR_' data_std.label{source} '.txt'], 'Delimiter','tab');
 
        
    end
    
    
end


%Save output


%Figures

saveas(figure(1), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'IFGvolassoc' '_' 'DrugDiff3_Pats' '_fullsample.tif']);


saveas(figure(2), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'IFGvolassoc' '_' 'meanMMN3_DrugDiff_Pats' '_fullsample.tif']);


savefig(figure(2), [FigOutDir '/' figoutprefix '_' 'LFPs_'  'IFGvolassoc' '_' 'meanMMN3_DrugDiff_Pats' '_fullsample.fig']);




%% Finished

sprintf('finito, output saved to %s \n', FigOutDir)
   

end