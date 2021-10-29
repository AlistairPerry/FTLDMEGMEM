function sourceLFPanalysis_basicdruggroup

% mutualsubjscon='Con_BothSess_subjs.txt';
% mutualsubjspat='Pat_BothSess_subjs.txt';
% LFPplaclistcon, LFPdruglistcon, LFPplaclistpat, LFPdruglistpat)
% [MMNpall] = sourceLFPanalysis_basiccond_druggroup_rmanov('Con_BothSess_subjs.txt', 'Pat_BothSess_subjs.txt', 'Con_AllP_Sess.txt', 'Con_AllD_Sess.txt', 'Pat_P_fnames.txt', 'Pat_D_fnames.txt')

%Open up SPM dirs (and field trip), and others

addpath('/imaging/rowe/users/ap09/Toolbox')

%for legend flex
addpath('/imaging/rowe/users/ap09/Toolbox/ekMEG/ntadscripts')

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest')

addpath('/imaging/rowe/users/ap09/Toolbox/fieldtrip')

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/boundedline-pkg'))

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/VBA-toolbox'))


ft_defaults()

FigOutDir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/analysis/SourceSpace/LFPs_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots/DrugAnalysis';

figoutprefix = 'adv_ssst_newmaxf_fixICA_wfids_nodipcor_nop10p19';


%Other setup

nsources = 6;

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

LFPpostfix = 'rov123_sss.mat'; %for now


Con_P = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Con_P.txt';

Con_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Con_D.txt';

Pat_P = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Pat_P.txt';

Pat_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Pat_D_nop10p19.txt';


fid = fopen(Con_P);
PlacData = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
PlacData = PlacData{1,1};

fid = fopen(Con_D);
DrugData = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
DrugData = DrugData{1,1};

fid = fopen(Pat_P);
PlacData_pat = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
PlacData_pat = PlacData_pat{1,1};

fid = fopen(Pat_D);
DrugData_pat = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
DrugData_pat = DrugData_pat{1,1};


%Diagnosis info

% Can call it directly as it will never change

DiagData = readcell('/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/PatDiagInfo.txt');




%% Controls


submatcount = 0;

placsubcount = 0;

drugsubcount = 0;


for cplacsubj = 1:length(PlacData(:,1))
    
    
    filen = [preprocdir '/' PlacData{cplacsubj} '/' LFPprefix '_' preprocsteps '_' PlacData{cplacsubj} '_' 's1' '_' LFPpostfix];
    
    filen_d = [preprocdir '/' PlacData{cplacsubj} '/' LFPprefix '_' preprocsteps '_' PlacData{cplacsubj} '_' 's2' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount = submatcount + 1;
        
        allcons_bothsess{submatcount,1} = PlacData{cplacsubj};
        
        
        D=spm_eeg_load(filen);
        
        D2=spm_eeg_load(filen_d);
        
        
       
        
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
        
        
        %Smooth data
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            
            %Std3
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix 
        
        avgstd_col(submatcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_col(submatcount,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_col(submatcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        avg_col_bothtrl(submatcount,:, :,1) = permute(data_dev.avg, [3 2 1]);
        
        avg_col_bothtrl(submatcount,:, :,2) = permute(data_std3.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev{submatcount,1}=data_dev;
        
        avgstd{submatcount,1}=data_std;
        
        avgstd3{submatcount,1}=data_std3;
        
        
        
        %And repeat for their drug session
        

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
            
            
            %Std
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into a 3D matrix
        
        avgstd_d_col(submatcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_d_col(submatcount,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_d_col(submatcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        avg_d_col_bothtrl(submatcount,:, :,1) = permute(data_dev.avg, [3 2 1]);
        
        avg_d_col_bothtrl(submatcount,:, :,2) = permute(data_std3.avg, [3 2 1]);
        
        
        
        
        %Combine into a group cell
        
        avgdev_d{submatcount,1}=data_dev;
        
        avgstd_d{submatcount,1}=data_std;
        
        avgstd3_d{submatcount,1}=data_std3;
        
    end
    
    
    %Pull out even if individual session
    
    
    if isfile(filen)
        
        
        placsubcount = placsubcount + 1;
        
        
        D=spm_eeg_load(filen);
        
        
       
        
        
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
        
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_all_col(placsubcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgdev_all_col(placsubcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        %Combine into a group cell
        
        avgdev_all{placsubcount,1}=data_dev;
        
        avgstd_all{placsubcount,1}=data_std;
        
        
    end
    
    
    if isfile(filen_d)
        
        drugsubcount = drugsubcount + 1;
        
        
        D=spm_eeg_load(filen_d);
        

        
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
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_d_all_col(drugsubcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgdev_d_all_col(drugsubcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        %Combine into a group cell
        
        avgdev_d{drugsubcount,1}=data_dev;
        
        avgstd_d{drugsubcount,1}=data_std;
        
        
    end
    
    
end



%Repeat for individuals with drug session first


for cdrugsubj = 1:length(DrugData(:,1))
    
    
    filen = [preprocdir '/' DrugData{cdrugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData{cdrugsubj} '_' 's2' '_' LFPpostfix];
    
    filen_d = [preprocdir '/' DrugData{cdrugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData{cdrugsubj} '_' 's1' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount = submatcount + 1;
        
        allcons_bothsess{submatcount,1} = DrugData{cdrugsubj};
        
        
        D=spm_eeg_load(filen);
        
        D2=spm_eeg_load(filen_d);
        
        
        
        
        
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
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            
            %Std
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into a 3D matrix
        
        avgstd_col(submatcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_col(submatcount,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_col(submatcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        avg_col_bothtrl(submatcount,:, :,1) = permute(data_dev.avg, [3 2 1]);
        
        avg_col_bothtrl(submatcount,:, :,2) = permute(data_std3.avg, [3 2 1]);
        
        
        
        
        %Combine into a group cell
        
        avgdev{submatcount,1}=data_dev;
        
        avgstd{submatcount,1}=data_std;
        
        avgstd3{submatcount,1}=data_std3;
        
        
        
        
        %And repeat for drug session
        
   
        
        %Other processing, including baseline correction and smoothing
        
        timelock=D2.fttimelock;
        
        
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
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matric
        
        avgstd_d_col(submatcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_d_col(submatcount,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_d_col(submatcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        avg_d_col_bothtrl(submatcount,:, :,1) = permute(data_dev.avg, [3 2 1]);
        
        avg_d_col_bothtrl(submatcount,:, :,2) = permute(data_std3.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev_d{submatcount,1}=data_dev;
        
        avgstd_d{submatcount,1}=data_std;
        
        avgstd3_d{submatcount,1}=data_std3;
        
    end
    
    
    %Now for all their individual sessions
    
    
    if isfile(filen)
        
        
        placsubcount = placsubcount + 1;
        
        
        D=spm_eeg_load(filen);
        
        

        
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
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_all_col(placsubcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgdev_all_col(placsubcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        

        
        
        %Combine into a group cell
        
        avgdev_all{placsubcount,1}=data_dev;
        
        avgstd_all{placsubcount,1}=data_std;
        
        
    end
    
    
    if isfile(filen_d)
        
        drugsubcount = drugsubcount + 1;
        
        
        D=spm_eeg_load(filen_d);
        
        
       
        
        
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
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_d_all_col(drugsubcount,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgdev_d_all_col(drugsubcount,:, :) = permute(data_dev.avg, [3 2 1]);
        
        

        
        
        %Combine into a group cell
        
        avgdev_d_all{drugsubcount,1}=data_dev;
        
        avgstd_d_all{drugsubcount,1}=data_std;
        
        
    end
    
    
end


%% Patients


submatcount_pat = 0;

placsubcount_pat = 0;

drugsubcount_pat = 0;


for pplacsubj = 1:length(PlacData_pat(:,1))
    
    
    filen = [preprocdir '/' PlacData_pat{pplacsubj} '/' LFPprefix '_' preprocsteps '_' PlacData_pat{pplacsubj} '_' 's1' '_' LFPpostfix];
    
    filen_d = [preprocdir '/' PlacData_pat{pplacsubj} '/' LFPprefix '_' preprocsteps '_' PlacData_pat{pplacsubj} '_' 's2' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount_pat = submatcount_pat + 1;
        
        allpats_bothsess{submatcount_pat,1} = PlacData_pat{pplacsubj};
        
        
        %find also diag info
        
        PatMatch=strcmp(DiagData(:,1), PlacData_pat{pplacsubj});

        extdiaginfo(submatcount_pat,1) = DiagData(PatMatch==1,2);
        
        
        D=spm_eeg_load(filen);
        
        D2=spm_eeg_load(filen_d);
        
        
       
        
        
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
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_pat_col(submatcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        avg_pat_col_bothtrl(submatcount_pat,:, :,1) = permute(data_dev.avg, [3 2 1]);
        
        avg_pat_col_bothtrl(submatcount_pat,:, :,2) = permute(data_std3.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev_pat{submatcount_pat,1}=data_dev;
        
        avgstd_pat{submatcount_pat,1}=data_std;
        
        avgstd3_pat{submatcount_pat,1}=data_std3;
        
        
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
        
        
        avg_d_pat_col_bothtrl(submatcount_pat,:, :,1) = permute(data_dev.avg, [3 2 1]);
        
        avg_d_pat_col_bothtrl(submatcount_pat,:, :,2) = permute(data_std3.avg, [3 2 1]);
        
        
        
        %Combine into a group cell
        
        avgdev_d_pat{submatcount_pat,1}=data_dev;
        
        avgstd_d_pat{submatcount_pat,1}=data_std;
        
        avgstd3_d_pat{submatcount_pat,1}=data_std3;
        
        
    end
    
    
    if isfile(filen)
        
        
        placsubcount_pat = placsubcount_pat + 1;
        
        
        D=spm_eeg_load(filen);
        
         
        
        
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
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_all_pat_col(placsubcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgdev_all_pat_col(placsubcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        

        
        
        %Combine into a group cell
        
        avgdev_all_pat{placsubcount_pat,1}=data_dev;
        
        avgstd_all_pat{placsubcount_pat,1}=data_std;
        
        
    end
    
    
    if isfile(filen_d)
        
        
        drugsubcount_pat = drugsubcount_pat + 1;
        
        
        D=spm_eeg_load(filen_d);
        
        
      
        
        
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
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Smoothing
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_d_all_pat_col(drugsubcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgdev_d_all_pat_col(drugsubcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
 
        
        
        %Combine into a group cell
        
        avgdev_d_all_pat{drugsubcount_pat,1}=data_dev;
        
        avgstd_d_all_pat{drugsubcount_pat,1}=data_std;
        
        
    end
    
    
end


for pdrugsubj = 1:length(DrugData_pat(:,1))
    
    
    filen = [preprocdir '/' DrugData_pat{pdrugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData_pat{pdrugsubj} '_' 's2' '_' LFPpostfix];
    
    filen_d = [preprocdir '/' DrugData_pat{pdrugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData_pat{pdrugsubj} '_' 's1' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount_pat = submatcount_pat + 1;
        
        allpats_bothsess{submatcount_pat,1} = DrugData_pat{pdrugsubj};
        
        
        
        %find also diag info
        
        PatMatch=strcmp(DiagData(:,1), DrugData_pat{pdrugsubj});

        extdiaginfo(submatcount_pat,1) = DiagData(PatMatch==1,2);
        
        
        D=spm_eeg_load(filen);
        
        D2=spm_eeg_load(filen_d);
        
        
       
        
        
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
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            
            %Std3
            
            data_std3.avg(source,:) = smooth(data_std3.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_pat_col(submatcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgstd3_pat_col(submatcount_pat,:, :) = permute(data_std3.avg, [3 2 1]);
        
        avgdev_pat_col(submatcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        
        avg_pat_col_bothtrl(submatcount_pat,:, :,1) = permute(data_dev.avg, [3 2 1]);
        
        avg_pat_col_bothtrl(submatcount_pat,:, :,2) = permute(data_std3.avg, [3 2 1]);
        
   
        
        %Combine into a group cell
        
        avgdev_pat{submatcount_pat,1}=data_dev;
        
        avgstd_pat{submatcount_pat,1}=data_std;
        
        avgstd3_pat{submatcount_pat,1}=data_std3;
        
        
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

        
        avg_d_pat_col_bothtrl(submatcount_pat,:, :,1) = permute(data_dev.avg, [3 2 1]);
        
        avg_d_pat_col_bothtrl(submatcount_pat,:, :,2) = permute(data_std3.avg, [3 2 1]);
        
        
        
        
        %Combine into a group cell
        
        avgdev_d_pat{submatcount_pat,1}=data_dev;
        
        avgstd_d_pat{submatcount_pat,1}=data_std;
        
        avgstd3_d_pat{submatcount_pat,1}=data_std3;
        
    end
    
    
    %And for their indiviudal sessions
    
    
    if isfile(filen)
        
        
        placsubcount_pat = placsubcount_pat + 1;
        
        
        D=spm_eeg_load(filen);
        
        
        
        
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
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        %Smoothing
        
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        end
        
        
        %Combine into 3D matrix
        
        avgstd_all_pat_col(placsubcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgdev_all_pat_col(placsubcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        

        
        
        %Combine into a group cell
        
        avgdev_all_pat{placsubcount_pat,1}=data_dev;
        
        avgstd_all_pat{placsubcount_pat,1}=data_std;
        
        
    end
    
    
    if isfile(filen_d)
        
        drugsubcount_pat = drugsubcount_pat + 1;
        
        
        D=spm_eeg_load(filen_d);
        
        
      
        
        
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
        
        %Dev
        cfg = [];
        cfg.trials = timelock.trialinfo == 1;
        data_dev = ft_redefinetrial(cfg, baselinecorr);
        
        
        for source = 1:nsources
            
            %Std
            
            data_std.avg(source,:) = smooth(data_std.avg(source,:), 20, 'moving');
            
            %Dev
            
            data_dev.avg(source,:) = smooth(data_dev.avg(source,:),20,'moving');
            
        
        
      
        
        end
        
        
        %Combine into 3D matrix
        
        avgstd_d_all_pat_col(drugsubcount_pat,:, :) = permute(data_std.avg, [3 2 1]);
        
        avgdev_d_all_pat_col(drugsubcount_pat,:, :) = permute(data_dev.avg, [3 2 1]);
        
        

        
        
        
        %Combine into a group cell
        
        avgdev_d_all_pat{drugsubcount_pat,1}=data_dev;
        
        avgstd_d_all_pat{drugsubcount_pat,1}=data_std;
        
        
    end
    
    
    
    
end


%Save figure metadata

FigMetaData = struct;

FigMetaData.npats = length(allpats_bothsess);

FigMetaData.ncons = length(allcons_bothsess);


FigMetaData.PatIDs = allpats_bothsess;

FigMetaData.ConIDs = allcons_bothsess;

FigMetaData.DiagInfo = cell2mat(extdiaginfo);


save([FigOutDir '/' figoutprefix '_' 'LFP_'  'ControlsandPats_' date '_wbothsess_FigMetaData.mat'], 'FigMetaData');



%% Grand averages

cfg=[];

%Plac

grandavgdev  = ft_timelockgrandaverage(cfg, avgdev{:});

grandavgstd  = ft_timelockgrandaverage(cfg, avgstd{:});

grandavgdev_all = ft_timelockgrandaverage(cfg, avgdev_all{:});

grandavgstd_all = ft_timelockgrandaverage(cfg, avgstd_all{:});


%Drug

grandavgdev_d_pat  = ft_timelockgrandaverage(cfg, avgdev_d_pat{:});

grandavgstd_d_pat  = ft_timelockgrandaverage(cfg, avgstd_d_pat{:});

grandavgdev_d_all_pat  = ft_timelockgrandaverage(cfg, avgdev_d_all_pat{:});

grandavgstd_d_all_pat  = ft_timelockgrandaverage(cfg, avgstd_d_all_pat{:});


%Patient Plac

grandavgdev_pat  = ft_timelockgrandaverage(cfg, avgdev_pat{:});

grandavgstd_pat  = ft_timelockgrandaverage(cfg, avgstd_pat{:});

grandavgdev_all_pat  = ft_timelockgrandaverage(cfg, avgdev_all_pat{:});

grandavgstd_all_pat  = ft_timelockgrandaverage(cfg, avgstd_all_pat{:});


%Patient Drug

grandavgdev_d_pat  = ft_timelockgrandaverage(cfg, avgdev_d_pat{:});

grandavgstd_d_pat  = ft_timelockgrandaverage(cfg, avgstd_d_pat{:});

grandavgdev_d_all_pat  = ft_timelockgrandaverage(cfg, avgdev_d_all_pat{:});

grandavgstd_d_all_pat  = ft_timelockgrandaverage(cfg, avgstd_d_all_pat{:});



%% Basic plotting and analyses - all subjects and sessions


%Plotting as function of session


%Standards

for source = 1:nsources
    

    %Start with controls
    
    
    data = avgstd_all_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgstd_d_all_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    k=max(max([m1(1,:); m1_d(1,:)]))+max(max([s1(1,:);s1_d(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(1), subplot(2,3,source);
    
    figure(1);set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Normalised signal', 'FontSize', 10)
    title([grandavgdev.label(source) ' STD in Cons'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Patients in here as well
 

    %Start with controls
    
    
    data = avgstd_all_pat_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgstd_d_all_pat_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    k=max(max([m1(1,:); m1_d(1,:)]))+max(max([s1(1,:);s1_d(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(2), subplot(2,3,source);
    
    figure(2); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Normalised signal', 'FontSize', 10)
    title([grandavgdev.label(source) ' STD in Pats'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    

    
end
    

   %Save figure
    saveas(figure(1), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Controls' '_Std_wDrug_fullsample.tif']);
    
    
        %Save figure
    saveas(figure(2), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Pats' '_Std_wDrug_fullsample.tif']);
    
    
 
    %Deviants

    figcount = 2;
    
 for source = 1:nsources
   

    %Start with controls
    
    
    data = avgdev_all_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgdev_d_all_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    k=max(max([m1(1,:); m1_d(1,:)]))+max(max([s1(1,:);s1_d(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+1), subplot(2,3,source);
    
    figure(figcount+1); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Normalised signal', 'FontSize', 10)
    title([grandavgdev.label(source) ' Dev in Cons'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Patients in here as well


    %Start with controls
    
    
    data = avgdev_all_pat_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgdev_d_all_pat_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    m1(1,:),s1(1,:)
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+2), subplot(2,3,source);
    
    figure(figcount+2); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Normalised signal', 'FontSize', 10)
    title([grandavgdev.label(source) ' Dev in Pats'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    

    
end
    

   %Save figure
    saveas(figure(figcount+1), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Controls' '_Dev_wDrug_fullsample.tif']);
    
    
%Save figure
    saveas(figure(figcount+2), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Pats' '_Dev_wDrug_fullsample.tif']);
    
 
    
    
 %Standards (rep3) and deviants
 
 
 figcount = 4;
 
 
  for source = 1:nsources
    
      %Controls
      
    data = avg_col_bothtrl(:,:,source,:);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avg_d_col_bothtrl(:,:,source,:);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    %Patients
    
    data = avg_pat_col_bothtrl(:,:,source,:);
    
    m1_pat=squeeze(nanmean(data,1));
    s1_pat=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avg_d_pat_col_bothtrl(:,:,source,:);
    
    m1_d_pat=squeeze(nanmean(data_d,1));
    s1_d_pat=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    
    %Placebo first
    
    figure(figcount+1), subplot(2,3,source);
    
    figure(figcount+1); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    plot(m1(:,2),'--','LineWidth',2,'Color',[0/255 0/255 255/255]); hold on
    plot(m1(:,1),':','LineWidth',2,'Color',[0/255 0/255 255/255]);
    plot(m1_pat(:,2),'--','LineWidth',2,'Color',[255/255 0/255 0/255]);
    plot(m1_pat(:,1),':','LineWidth',2,'Color',[255/255 0/255 0/255]); %legendflex({'rep3 CON','dev CON','rep3 FTLD','dev FTLD'},'fontsize',12, 'fontweight', 'bold')
    boundedline([1:size(m1,1)],m1(:,2),s1(:,2),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,1)],m1(:,1),s1(:,1),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,1)],m1_pat(:,2),s1_pat(:,2),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,1)],m1_pat(:,1),s1_pat(:,1),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial'); ylabel({'fT/mm'}, 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    %title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
        box off

    set(gca, 'TickDir', 'Out')
    
    xlim([0 251]);
    
     
    
    %Now drug
    
       figure(figcount+2), subplot(2,3,source);
    
    figure(figcount+2); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    plot(m1_d(:,2),'--','LineWidth',2,'Color',[0/255 0/255 255/255]); hold on
    plot(m1_d(:,1),':','LineWidth',2,'Color',[0/255 0/255 255/255]);
    plot(m1_d_pat(:,2),'--','LineWidth',2,'Color',[255/255 0/255 0/255]);
    plot(m1_d_pat(:,2),':','LineWidth',2,'Color',[255/255 0/255 0/255]); %legendflex({'rep3 CON','dev CON','rep3 FTLD','dev FTLD'},'fontsize',12, 'fontweight', 'bold')
    boundedline([1:size(m1_d,1)],m1_d(:,2),s1_d(:,2),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_d,1)],m1_d(:,1),s1_d(:,1),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_d_pat,1)],m1_d_pat(:,2),s1_d_pat(:,2),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_d_pat,1)],m1_d_pat(:,1),s1_d_pat(:,1),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial'); ylabel({'Amplitude (a.u)'}, 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    %title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
        box off

    set(gca, 'TickDir', 'Out')
    
    xlim([0 251]);
      
  end
 
 
     %Save figure
    saveas(figure(figcount+1), [FigOutDir '/' figoutprefix '_' 'LFP_'  'ConandPats' '_Rep3Dev_wboth_PlacSess_fullsample.tif']);
    
    
%Save figure
    saveas(figure(figcount+2), [FigOutDir '/' figoutprefix '_' 'LFP_'  'ConandPats' '_Rep3Dev_wboth_DrugSess_fullsample.tif']);
    

 
 
 figcount = 6;
    
    
%Difference
    
for source = 1:nsources
    

    %Start with controls
    
    
    data = avgstd_all_col(:,:,source)-avgdev_all_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgstd_d_all_col(:,:,source)-avgdev_d_all_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    k=max(max([m1(1,:); m1_d(1,:)]))+max(max([s1(1,:);s1_d(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+1), subplot(2,3,source);
    
    figure(figcount+1); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Std - Dev', 'FontSize', 10)
    title([grandavgdev.label(source) ' Diff in Cons'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Patients in here as well
 

    %Start with controls
    
    
    data = avgstd_all_pat_col(:,:,source)-avgdev_all_pat_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgstd_d_all_pat_col(:,:,source)-avgdev_d_all_pat_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    k=max(max([m1(1,:); m1_d(1,:)]))+max(max([s1(1,:);s1_d(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+2), subplot(2,3,source);
    
    figure(figcount+2); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Std - Dev', 'FontSize', 10)
    title([grandavgdev.label(source) ' Diff in Pats'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    

    
end    

saveas(figure(figcount+1), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Controls_' 'Diff_wDrug_fullsample.tif']);

saveas(figure(figcount+2), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Pats_' '_Diff_wDrug_fullsample.tif']);


    
%% Now difference of difference for subjects with both sessions
%Rep6


meanMMNcol_all_allnodes = [];

meanP3col_all_allnodes = [];

peaklat_all_allnodes = [];


time = D.time;

bltime = find(time == 0);


%Setup first for ANOVA

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


%Design matrix for ANOVA..

nfilesc = size(avgstd_col,1);

nfilesp = size(avgstd_pat_col,1);


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



%Start rep6

figcount = 8;

for source = 1:nsources
   

    
    %Start with controls
    
    
    data = avgstd_col(:,:,source)-avgdev_col(:,:,source);
    
    data_d = avgstd_d_col(:,:,source)-avgdev_d_col(:,:,source);
    
    
    data_diff = data - data_d;
    
    
    
    %And now patients
    
    data_pat = avgstd_pat_col(:,:,source)-avgdev_pat_col(:,:,source);
    
    data_d_pat = avgstd_d_pat_col(:,:,source)-avgdev_d_pat_col(:,:,source);
    
    
    data_diff_pat = data_pat - data_d_pat;
    
    
    %Calculate means
    
    m1=squeeze(nanmean(data_diff,1));
    
    s1=squeeze(nanstd(data_diff)./sqrt(size(data_diff,1)));
    
    
    m1_pat = squeeze(nanmean(data_diff_pat,1));

    s1_pat = squeeze(nanstd(data_diff_pat)./sqrt(size(data_diff_pat,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data_diff),squeeze(data_diff_pat));
    k=max(max([m1(1,:); m1_pat(1,:)]))+max(max([s1(1,:);s1_pat(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+1), subplot(2,3,source);
    
    figure(figcount+1); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_pat(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'Pat','Con'},'fontsize',10)
    boundedline([1:size(m1_pat,2)],m1_pat(1,:),s1_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN PLA - MEM')
    title([grandavgdev.label(source) ' MMN DrugDiff'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
   
   %Amplitude
   %Since I have them loaded in can't I do the other analyses now?
    
    %MMN
    
    m1_mmn = nanmean(data(:,MMNs:MMNf),2);
    
    m1_mmn_d = nanmean(data_d(:,MMNs:MMNf),2);
    
    
    m1_mmn_pat = nanmean(data_pat(:,MMNs:MMNf),2);
    
    m1_mmn_d_pat = nanmean(data_d_pat(:,MMNs:MMNf),2);
    
    
    %And their interaction significance
    
    meanMMNcol_all = cat(1, m1_mmn, m1_mmn_d, m1_mmn_pat, m1_mmn_d_pat);
    
    meanMMNcol_all_allnodes = cat(2, meanMMNcol_all_allnodes, meanMMNcol_all);
    
    X = cat(2, meanMMNcol_all, F1all, F2all, Sall);
    
    
    [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
    
    
    pgrp = Ps(1);
    pcond = Ps(3);
    pint = Ps(4);
    
    MMNpall(source,1)=pgrp;
    MMNpall(source,2)=pcond;
    MMNpall(source,3)=pint;
    
    
    %And now P3 
    
    m1_p3 = nanmean(data(:,p3s:p3f),2);
    
    m1_p3_d = nanmean(data_d(:,p3s:p3f),2);
    
    
    m1_p3_pat = nanmean(data_pat(:,p3s:p3f),2);
    
    m1_p3_d_pat = nanmean(data_d_pat(:,p3s:p3f),2);
    
    
    
    %Calculating interaction signficance
    
    meanP3col_all = cat(1, m1_p3, m1_p3_d, m1_p3_pat, m1_p3_d_pat);
    
    meanP3col_all_allnodes = cat(2, meanP3col_all_allnodes, meanP3col_all);
    
    
    X = cat(2, meanP3col_all, F1all, F2all, Sall);
    
    
    [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
    
    
    pgrp = Ps(1);
    pcond = Ps(3);
    pint = Ps(4);
    
    P3pall(source,1)=pgrp;
    P3pall(source,2)=pcond;
    P3pall(source,3)=pint;
    
    
    
    %Write out data for JASP Bayes RM ANOVA - would be too complex to write it out
    %containing multiple nodes
    
    meanMMNcol_lw = cat(2, m1_mmn, m1_mmn_d);
    
    meanMMNcol_lw_pat = cat(2, m1_mmn_pat, m1_mmn_d_pat);
    
    meanMMNcol_lw_all = cat(1, meanMMNcol_lw, meanMMNcol_lw_pat);

    meanMMNcol_lw_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));

    meanMMNcol_lw_subgrp = cat(1, repmat(3, [nfilesc 1]), cell2mat(extdiaginfo));
    
    
    MMNtable_lw = table(meanMMNcol_lw_all(:,1), meanMMNcol_lw_all(:,2), meanMMNcol_lw_grp, meanMMNcol_lw_subgrp);
    
    
    MMNtable_lw.Properties.VariableNames = {'meanMMNcol_PLA', 'meanMMNcol_MEM', 'Group', 'PatSubGrp'};

    writetable(MMNtable_lw, [FigOutDir '/' figoutprefix '_' 'LFP_' grandavgdev.label{source} '_MMNmean' '_' 'ConsandPats' '_DrugInt_' 'fullsample_forJASP.txt'], 'Delimiter','tab')
    
   
    
    %Peak latency
    
    subjthm_wid = [];
    subjthm_ups = [];
    subjthm_wid_d = [];
    subjthm_ups_d = [];
    
    subjthm_wid_pat = [];
    subjthm_ups_pat = [];
    subjthm_wid_d_pat = [];
    subjthm_ups_d_pat = [];
        
        
    for csubj = 1:nfilesc
        
        
        %Placebo
        
        subjbl = data(csubj, bltime);
        
        
        %100 to 150 equals 100 to 200ms
        subjpeak = min(data(csubj,100:150));
        
        subjipeak = find(data(csubj,:)==subjpeak);
        
        subjtpeak = time(find(data(csubj,:) == subjpeak));
            
            
            subjtpeak_all(csubj,1) = subjtpeak;
            
            
%             if subjtpeak > 0.250
%                 
%                 
%             subjthm_ups(csubj,1) = 0;
%             
%             subjthm_wid(csubj,1) = 0;
%                 
%                 
%             else
%             
%         
%             subjhm = (subjpeak - subjbl)./2;
%             
%             
%             %Find them using closest to the point (half max) through subtraction
%             
%             hmminusdiff = subjhm-diff(csubj,:);
%             
%             
%             %To left of peak
%             
%             lefthmind = find(abs(hmminusdiff(1:subjipeak))==min(abs(hmminusdiff(1:subjipeak))));
%             
%             
%             %And right
%             
%             righthmind = find(abs(hmminusdiff(1,subjipeak:188))==min(abs(hmminusdiff(1,subjipeak:188))));
%             
%             
%             subj_ltimehm = time(lefthmind);
%             
%             subj_rtimehm = time(subjipeak + righthmind);
%             
%             
%             subjthm_ups(csubj,1) = subj_rtimehm;
%             
%             subjthm_wid(csubj,1) = subj_rtimehm-subj_ltimehm;
            
            
            
            
            %Drug
            
            subjbl = data_d(csubj, bltime);
            
            subjpeak = min(data_d(csubj,100:150));
            
            subjipeak = find(data_d(csubj,:)==subjpeak);
            
            subjtpeak = time(find(data_d(csubj,:) == subjpeak));
            
            subjtpeak_d_all(csubj,1) = subjtpeak;
            
            
%             if subjtpeak > 0.250
%                 
%                 
%                 subjthm_ups_d(csubj,1) = 0;
%                 
%                 subjthm_wid_d(csubj,1) = 0;
%                 
%                 
%             else
%                 
%                 
%                 subjhm = (subjpeak - subjbl)./2;
%                 
%                 
%                 %Find them using closest to the point (half max) through subtraction
%                 
%                 hmminusdiff_d = subjhm-diff_d(csubj,:);
%                 
%                 
%                 %To left of peak
%                 
%                 lefthmind = find(abs(hmminusdiff_d(1:subjipeak))==min(abs(hmminusdiff_d(1:subjipeak))));
%                 
%                 
%                 %And right
%                 
%                 righthmind = find(abs(hmminusdiff_d(1,subjipeak:188))==min(abs(hmminusdiff_d(1,subjipeak:188))));
%                 
%                 
%                 subj_ltimehm = time(lefthmind);
%                 
%                 subj_rtimehm = time(subjipeak + righthmind);
%                 
%                 
%                 subjthm_ups_d(csubj,1) = subj_rtimehm;
%                 
%                 subjthm_wid_d(csubj,1) = subj_rtimehm-subj_ltimehm;
%             
%             
%             end
        
  
    end
        
        
        %Repeat for patients
        
        
        for psubj = 1:nfilesp
            
            
            %Placebo
            
            subjbl = data_pat(psubj, bltime);
            
            subjpeak = min(data_pat(psubj,100:150));
            
            subjipeak = find(data_pat(psubj,:)==subjpeak);
            
            subjtpeak = time(find(data_pat(psubj,:) == subjpeak));
            
            subjtpeak_pat_all(psubj,1) = subjtpeak;
            
            
%             if subjtpeak > 0.250
%                 
%                 
%             subjthm_ups_pat(psubj,1) = 0;
%             
%             subjthm_wid_pat(psubj,1) = 0;
%                 
%                 
%             else
%             
%             subjhm = (subjpeak - subjbl)./2;
%             
%             
%             %Find them using closest to the point (half max) through subtraction
%             
%             hmminusdiff_pat = subjhm-diff_pat(psubj,:);
%             
%             
%             %To left of peak
%             
%             lefthmind = find(abs(hmminusdiff_pat(1:subjipeak))==min(abs(hmminusdiff_pat(1:subjipeak))));
%             
%             
%             %And right
%             
%             righthmind = find(abs(hmminusdiff_pat(1,subjipeak:188))==min(abs(hmminusdiff_pat(1,subjipeak:188))));
%             
%             
%             subj_ltimehm = time(lefthmind);
%             
%             subj_rtimehm = time(subjipeak + righthmind);
%             
%             
%             subjthm_ups_pat(psubj,1) = subj_rtimehm;
%             
%             subjthm_wid_pat(psubj,1) = subj_rtimehm-subj_ltimehm;
%             
%             
%             end
            
            
            %Drug
            
            subjbl = data_d_pat(psubj, bltime);
            
            subjpeak = min(data_d_pat(psubj,100:150));
            
            subjipeak = find(data_d_pat(psubj,:)==subjpeak);
            
            subjtpeak = time(find(data_d_pat(psubj,:) == subjpeak));
            
            
            subjtpeak_d_pat_all(psubj,1) = subjtpeak;
            
            
%             if subjtpeak > 0.250
%                 
%                 
%             subjthm_ups_d_pat(psubj,1) = 0;
%             
%             subjthm_wid_d_pat(psubj,1) = 0;
%                 
%                 
%             else
%             
%             
%             subjhm = (subjpeak - subjbl)./2;
%             
%             
%             %Find them using closest to the point (half max) through subtraction
%             
%             hmminusdiff_d_pat = subjhm-diff_d_pat(psubj,:);
%             
%             
%             %To left of peak
%             
%             lefthmind = find(abs(hmminusdiff_d_pat(1:subjipeak))==min(abs(hmminusdiff_d_pat(1:subjipeak))));
%             
%             
%             %And right
%             
%             righthmind = find(abs(hmminusdiff_d_pat(1,subjipeak:188))==min(abs(hmminusdiff_d_pat(1,subjipeak:188))));
%             
%             
%             subj_ltimehm = time(lefthmind);
%             
%             subj_rtimehm = time(subjipeak + righthmind);
%             
%             
%             subjthm_ups_d_pat(psubj,1) = subj_rtimehm;
%             
%             subjthm_wid_d_pat(psubj,1) = subj_rtimehm-subj_ltimehm;
%             
%             
%             end
            
        end
   
        
%     %Calculate FWHM interaction significance    
%         
%         
%     FWHM_all = cat(1,  subjthm_wid, subjthm_wid_d, subjthm_wid_pat, subjthm_wid_d_pat);
%     
%     
%     X = cat(2, FWHM_all, F1all, F2all, Sall);
%     
%      
%     zerovals = find(FWHM_all(:,1)==0);
%     
%     X(zerovals(:,1),:) = [];
%     
%     
%     %[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
% 
%     
%     Y = X(:,1);
%     
%     Group{1,1} = X(:,2);
%     
%     Group{1,2} = X(:,3);
%     
%     
%     [p, ~, ~ , ~] = anovan(Y, Group, 'model', 'interaction', 'varnames', {'Group', 'Drug'}, 'display', 'off');
% 
%     
%     FWHMpall(i-3,:) = p;
%     
%     
%     %Can push out into rainclouds for plotting in R :)
%     
%     FWHMtable = table(X(:,1), X(:,2), X(:,3));
%     
%     FWHMtable.Properties.VariableNames = {'FWHMwid', 'Group', 'Drug'};
%     
%     writetable(FWHMtable, [FigOutDir '/' OUTpre '_' 'ERF_'  'FWHMwid' '_' 'ConsandPats_DrugInt_' varnm{i} '_fullsample.txt'], 'Delimiter','tab')
    
    
    %MMN peak latency interaction significance
    
    peaklat_all = cat(1, subjtpeak_all, subjtpeak_d_all, subjtpeak_pat_all, subjtpeak_d_pat_all);
    
    peaklat_all_allnodes = cat(2, peaklat_all_allnodes, peaklat_all);
    
    
    X = cat(2, peaklat_all, F1all, F2all, Sall);
    
    
    [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
    
    
    pgrp = Ps(1);
    pcond = Ps(3);
    pint = Ps(4);
    
    peaklatpall(source,1)=pgrp;
    peaklatpall(source,2)=pcond;
    peaklatpall(source,3)=pint;

    
    %Write out data for JASP Bayes RM ANOVA - would be too complex to write it out
    %containing multiple nodes
    
    peaklatcol_lw = cat(2, subjtpeak_all, subjtpeak_d_all);
    
    peaklatcol_lw_pat = cat(2, subjtpeak_pat_all, subjtpeak_d_pat_all);
    
    peaklatcol_lw_all = cat(1, peaklatcol_lw, peaklatcol_lw_pat);

    peaklatcol_lw_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));
    
    peaklattable_lw = table(peaklatcol_lw_all(:,1), peaklatcol_lw_all(:,2), peaklatcol_lw_grp);
    
    
    peaklattable_lw.Properties.VariableNames = {'peaklatcol_PLA', 'peaklatcol_MEM', 'Group'};

    writetable(peaklattable_lw, [FigOutDir '/' figoutprefix '_' 'LFP_' grandavgdev.label{source} '_peaklat' '_' 'ConsandPats' '_DrugInt_' 'fullsample_forJASP.txt'], 'Delimiter','tab')
    
  

end
    

%And full data tables


%Can push out into rainclouds for plotting in R :)

% MMN

MMNtable = table(meanMMNcol_all_allnodes(:,1), meanMMNcol_all_allnodes(:,2), meanMMNcol_all_allnodes(:,3), meanMMNcol_all_allnodes(:,4), meanMMNcol_all_allnodes(:,5), meanMMNcol_all_allnodes(:,6), F1all, F2all, Sall);

MMNtable.Properties.VariableNames = {'meanMMNcol_LIFG', 'meanMMNcol_LSTG', 'meanMMNcol_LAUD', 'meanMMNcol_RIFG', 'meanMMNcol_RSTG', 'meanMMNcol_RAUD', 'Group', 'Drug', 'ID'};

writetable(MMNtable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'MMNmean' '_' 'ConsandPats' '_DrugInt_' 'fullsample.txt'], 'Delimiter','tab')


%P3

P3table = table(meanP3col_all_allnodes(:,1), meanP3col_all_allnodes(:,2), meanP3col_all_allnodes(:,3), meanP3col_all_allnodes(:,4), meanP3col_all_allnodes(:,5), meanP3col_all_allnodes(:,6), F1all, F2all, Sall);

P3table.Properties.VariableNames = {'meanP3col_LIFG', 'meanP3col_LSTG', 'meanP3col_LAUD', 'meanP3col_RIFG', 'meanP3col_RSTG', 'meanP3col_RAUD', 'Group', 'Drug', 'ID'};

writetable(P3table, [FigOutDir '/' figoutprefix '_' 'LFP_'  'P3mean' '_' 'ConsandPats' '_DrugInt_' 'fullsample.txt'], 'Delimiter','tab')


%Peak latency

peaklattable = table(peaklat_all_allnodes(:,1), peaklat_all_allnodes(:,2), peaklat_all_allnodes(:,3), peaklat_all_allnodes(:,4), peaklat_all_allnodes(:,5), peaklat_all_allnodes(:,6), F1all, F2all, Sall);

peaklattable.Properties.VariableNames = {'peaklatcol_LIFG', 'peaklatcol_LSTG', 'peaklatcol_LAUD', 'peaklatcol_RIFG', 'peaklatcol_RSTG', 'peaklatcol_RAUD', 'Group', 'Drug', 'ID'};

writetable(peaklattable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'peaklat' '_' 'ConsandPats' '_DrugInt_' 'fullsample.txt'], 'Delimiter','tab')
    


%Save figure
saveas(figure(figcount+1), [FigOutDir '/' figoutprefix '_' 'LFP_'  'ControlsandPats' '_Diff_DrugDiff_fullsample.tif']);

    
%Save output tables of ANOVA
%MMN 

summrestable = table(MMNpall(:,1), MMNpall(:,2), MMNpall(:,3));
summrestable.Properties.VariableNames = {'Group' 'Drug' 'Int'};

writetable(summrestable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'MMNmean' '_restable_2wayANOVA_ConsandPats_DrugInt.txt']);


%And P3

summrestable = table(P3pall(:,1), P3pall(:,2), P3pall(:,3));
summrestable.Properties.VariableNames = {'Group' 'Drug' 'Int'};

writetable(summrestable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'P3mean' '_restable_2wayANOVA_ConsandPats_DrugInt.txt']);


%And lastly peak latency

summrestable = table(peaklatpall(:,1), peaklatpall(:,2), peaklatpall(:,3));
summrestable.Properties.VariableNames = {'Group' 'Drug' 'Int'};

writetable(summrestable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'peaklat' '_restable_2wayANOVA_ConsandPats_DrugInt.txt']);

    
    %
    %     figure(2), subplot(2,3,source)
    %
    %     hold on
    %
    %
    %     for isub = 1:nsubjs
    %
    %         %plot individuals across each group and their connectivity trajectories
    %         %first group one (plac)
    %
    %         figure(2), plot([1 2], [meanMMNcol([isub-1]*4+1) meanMMNcol([isub-1]*4+2)], '-b','LineWidth',1,'Marker','o','MarkerFaceColor','b','MarkerEdgeColor','k','MarkerSize',6)
    %
    %         %now group two (drug)
    %
    %         figure(2), plot([1 2], [meanMMNcol([isub-1]*4+3) meanMMNcol([isub-1]*4+4)], '-r','LineWidth',1,'Marker','o','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',6)
    %
    %     end
    %
    %     title([data_std_plac.label{source, 1}])
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
    %     ax.XTickLabel = {'Std' 'Dev'};
    
    %
    %     [h, ~, ~, adj_p]=fdr_bh(pall);
    %
    %
    %     if nnz(h)~=0
    %
    %         sigpout(h==1)=adj_p(h==1);
    %
    %     end
    
    
    
    %    statoutput(:, source)=sigpout;
    
%% Rep 3


%Plotting as function of session

figcount = 9;

for source = 1:nsources
    

    %Start with controls
    
    
    data = avgstd3_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgstd3_d_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    k=max(max([m1(1,:); m1_d(1,:)]))+max(max([s1(1,:);s1_d(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    
    figure(figcount+1), subplot(2,3,source);
    
    figure(figcount+1);set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(:,1),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(:,1),s1(:,1),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Normalised signal', 'FontSize', 10)
    title([grandavgdev.label(source) ' STD3 in Cons'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    %Patients in here as well    
    
    data = avgstd3_pat_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgstd3_d_pat_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    k=max(max([m1(1,:); m1_d(1,:)]))+max(max([s1(1,:);s1_d(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+2), subplot(2,3,source);
    
    figure(figcount+2); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Normalised signal', 'FontSize', 10)
    title([grandavgdev.label(source) ' STD3 in Pats'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square

    
end
    

%Save figure
saveas(figure(figcount+1), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Controls' '_Std3_wDrug_fullsample.tif']);


%Save figure
saveas(figure(figcount+2), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Pats' '_Std3_wDrug_fullsample.tif']);

    

%Difference
    
figcount = 11;


for source = 1:nsources
    

    %Start with controls
    
    
    data = avgstd3_col(:,:,source)-avgdev_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgstd3_d_col(:,:,source)-avgdev_d_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    k=max(max([m1(1,:); m1_d(1,:)]))+max(max([s1(1,:);s1_d(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+1), subplot(2,3,source);
    
    figure(figcount+1); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Std - Dev', 'FontSize', 10)
    title([grandavgdev.label(source) ' MMN3 in Cons'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Patients in here as well
 

    %Start with controls
    
    
    data = avgstd3_pat_col(:,:,source)-avgdev_pat_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgstd3_pat_col(:,:,source)-avgdev_d_pat_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data),squeeze(data_d));
    k=max(max([m1(1,:); m1_d(1,:)]))+max(max([s1(1,:);s1_d(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+2), subplot(2,3,source);
    
    figure(figcount+2); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_d(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'MEM','PLA'},'fontsize',10)
    boundedline([1:size(m1_d,2)],m1_d(1,:),s1_d(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('Std - Dev', 'FontSize', 10)
    title([grandavgdev.label(source) ' MMN3 in Pats'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    

    
end    


saveas(figure(figcount+1), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Controls_' 'Diff3_wDrug_fullsample.tif']);

saveas(figure(figcount+2), [FigOutDir '/' figoutprefix '_' 'LFP_'  'Pats_' '_Diff3_wDrug_fullsample.tif']);



%Difference of differences

figcount = 13;

for source = 1:nsources
   
    
    %Start with controls
    
    
    data = avgstd3_col(:,:,source)-avgdev_col(:,:,source);
    
    data_d = avgstd3_d_col(:,:,source)-avgdev_d_col(:,:,source);
    
    
    data_diff = data - data_d;
    
    
    
    %And now patients
    
    data_pat = avgstd3_pat_col(:,:,source)-avgdev_pat_col(:,:,source);
    
    data_d_pat = avgstd3_d_pat_col(:,:,source)-avgdev_d_pat_col(:,:,source);
    
    
    data_diff_pat = data_pat - data_d_pat;
    
    
    %Calculate means
    
    m1=squeeze(nanmean(data_diff,1));
    
    s1=squeeze(nanstd(data_diff)./sqrt(size(data_diff,1)));
    
    
    m1_pat = squeeze(nanmean(data_diff_pat,1));

    s1_pat = squeeze(nanstd(data_diff_pat)./sqrt(size(data_diff_pat,1)));
    
    
    [h,p,ci,stats]=ttest2(squeeze(data_diff),squeeze(data_diff_pat));
    k=max(max([m1(1,:); m1_pat(1,:)]))+max(max([s1(1,:);s1_pat(1,:)]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+1), subplot(2,3,source);
    
    figure(figcount+1); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_pat(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'Pat','Con'},'fontsize',10)
    boundedline([1:size(m1_pat,2)],m1_pat(1,:),s1_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN PLA - MEM')
    title([grandavgdev.label(source) ' MMN3 DrugDiff'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
   
   %Amplitude
   %Since I have them loaded in can't I do the other analyses now?
    
    %MMN
    
    m1_mmn = nanmean(data(:,MMNs:MMNf),2);
    
    m1_mmn_d = nanmean(data_d(:,MMNs:MMNf),2);
    
    
    m1_mmn_pat = nanmean(data_pat(:,MMNs:MMNf),2);
    
    m1_mmn_d_pat = nanmean(data_d_pat(:,MMNs:MMNf),2);
    
    
    %And their interaction significance
    
    meanMMNcol_all = cat(1, m1_mmn, m1_mmn_d, m1_mmn_pat, m1_mmn_d_pat);
    
    meanMMNcol_all_allnodes = cat(2, meanMMNcol_all_allnodes, meanMMNcol_all);
    
    X = cat(2, meanMMNcol_all, F1all, F2all, Sall);
    
    
    [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
    
    
    pgrp = Ps(1);
    pcond = Ps(3);
    pint = Ps(4);
    
    MMNpall(source,1)=pgrp;
    MMNpall(source,2)=pcond;
    MMNpall(source,3)=pint;
    
    
    %And now P3 
    
    m1_p3 = nanmean(data(:,p3s:p3f),2);
    
    m1_p3_d = nanmean(data_d(:,p3s:p3f),2);
    
    
    m1_p3_pat = nanmean(data_pat(:,p3s:p3f),2);
    
    m1_p3_d_pat = nanmean(data_d_pat(:,p3s:p3f),2);
    
    
    
    %Calculating interaction signficance
    
    meanP3col_all = cat(1, m1_p3, m1_p3_d, m1_p3_pat, m1_p3_d_pat);
    
    meanP3col_all_allnodes = cat(2, meanP3col_all_allnodes, meanP3col_all);
    
    
    X = cat(2, meanP3col_all, F1all, F2all, Sall);
    
    
    [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
    
    
    pgrp = Ps(1);
    pcond = Ps(3);
    pint = Ps(4);
    
    P3pall(source,1)=pgrp;
    P3pall(source,2)=pcond;
    P3pall(source,3)=pint;
    
    
    
    %Write out data for JASP Bayes RM ANOVA - would be too complex to write it out
    %containing multiple nodes
    
    meanMMNcol_lw = cat(2, m1_mmn, m1_mmn_d);
    
    meanMMNcol_lw_pat = cat(2, m1_mmn_pat, m1_mmn_d_pat);
    
    meanMMNcol_lw_all = cat(1, meanMMNcol_lw, meanMMNcol_lw_pat);

    meanMMNcol_lw_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));
    
    meanMMNcol_lw_subgrp = cat(1, repmat(3, [nfilesc 1]), cell2mat(extdiaginfo));
    
    
    MMNtable_lw = table(meanMMNcol_lw_all(:,1), meanMMNcol_lw_all(:,2), meanMMNcol_lw_grp, meanMMNcol_lw_subgrp);
    
    
    MMNtable_lw.Properties.VariableNames = {'meanMMNcol_PLA', 'meanMMNcol_MEM', 'Group', 'PatSubGrp'};


   writetable(MMNtable_lw, [FigOutDir '/' figoutprefix '_' 'LFP_' grandavgdev.label{source} '_MMNmean' '_rep3_' 'ConsandPats' '_DrugInt_' 'fullsample_forJASP.txt'], 'Delimiter','tab')
    


end
    

%And full data tables


%Can push out into rainclouds for plotting in R :)

% MMN

MMNtable = table(meanMMNcol_all_allnodes(:,1), meanMMNcol_all_allnodes(:,2), meanMMNcol_all_allnodes(:,3), meanMMNcol_all_allnodes(:,4), meanMMNcol_all_allnodes(:,5), meanMMNcol_all_allnodes(:,6), F1all, F2all, Sall);

MMNtable.Properties.VariableNames = {'meanMMNcol_LIFG', 'meanMMNcol_LSTG', 'meanMMNcol_LAUD', 'meanMMNcol_RIFG', 'meanMMNcol_RSTG', 'meanMMNcol_RAUD', 'Group', 'Drug', 'ID'};

writetable(MMNtable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'MMN3mean' '_' 'ConsandPats' '_DrugInt_' 'fullsample.txt'], 'Delimiter','tab')


%P3

P3table = table(meanP3col_all_allnodes(:,1), meanP3col_all_allnodes(:,2), meanP3col_all_allnodes(:,3), meanP3col_all_allnodes(:,4), meanP3col_all_allnodes(:,5), meanP3col_all_allnodes(:,6), F1all, F2all, Sall);

P3table.Properties.VariableNames = {'meanP3col_LIFG', 'meanP3col_LSTG', 'meanP3col_LAUD', 'meanP3col_RIFG', 'meanP3col_RSTG', 'meanP3col_RAUD', 'Group', 'Drug', 'ID'};

writetable(P3table, [FigOutDir '/' figoutprefix '_' 'LFP_'  'P3mean' '_rep3_' 'ConsandPats' '_DrugInt_' 'fullsample.txt'], 'Delimiter','tab')




%Save figure
saveas(figure(figcount+1), [FigOutDir '/' figoutprefix '_' 'LFP_' 'rep3_' 'ControlsandPats' '_Diff_DrugDiff_fullsample.tif']);

    
%Save output tables of ANOVA
%MMN 

summrestable = table(MMNpall(:,1), MMNpall(:,2), MMNpall(:,3));
summrestable.Properties.VariableNames = {'Group' 'Drug' 'Int'};

writetable(summrestable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'MMNmean' '_rep3_restable_2wayANOVA_ConsandPats_DrugInt.txt']);


%And P3

summrestable = table(P3pall(:,1), P3pall(:,2), P3pall(:,3));
summrestable.Properties.VariableNames = {'Group' 'Drug' 'Int'};

writetable(summrestable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'P3mean' '_rep3_restable_2wayANOVA_ConsandPats_DrugInt.txt']);



    
%%=========================================================================
%% RFT
%%=========================================================================

% figcount = 13;
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
% for source = 4:nsources
%    
%     
%     %Start with controls
%     
%     
%     data = avgstd3_col(:,:,source)-avgdev_col(:,:,source);
%     
%     data_d = avgstd3_d_col(:,:,source)-avgdev_d_col(:,:,source);
%     
%     
%     data_diff = data - data_d;
%     
%     
%     
%     %And now patients
%     
%     data_pat = avgstd3_pat_col(:,:,source)-avgdev_pat_col(:,:,source);
%     
%     data_d_pat = avgstd3_d_pat_col(:,:,source)-avgdev_d_pat_col(:,:,source);
%     
%     
%     data_diff_pat = data_pat - data_d_pat;
%    
%     
%     
%     figcount = figcount + 1;
%         
%         
%     mmndiffall = cat(1, data_diff, data_diff_pat);
%     
%     %If want to crop baseline period
%     %y = mmndiffall(:,51:end);
%     
%     y = mmndiffall;
%     
%     
%     [stat,out] = RFT_GLM_contrast(DM,y,c,'F',1,1);
%     
%     
%     %Save figure
%     
%     saveas(figure(figcount), [FigOutDir '/' figoutprefix '_' 'LFP_' 'rep3_' 'ControlsandPats_' grandavgdev.label{source} '_Diff_DrugDiff_fullsample_RFT.tif']);
%     
%     
%     %Saveoutput
%     
%     save([FigOutDir '/' figoutprefix '_' 'LFP_' 'rep3_' 'ControlsandPats_' grandavgdev.label{source} '_Diff_DrugDiff_fullsample_RFT.mat'], 'stat');
%     
%     
% end
    


%%=========================================================================
%% DiagnosisGroup
%%=========================================================================


%First determine number of bv/psp

extdiaginfo = cell2mat(extdiaginfo);

proc_pats_bv = find(extdiaginfo==1);
proc_pats_psp = find(extdiaginfo==2);

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



%Rep 6

figcount = 14;


for source = 1:nsources
   

    
    %Start with controls
    
    
    data = avgstd_col(:,:,source)-avgdev_col(:,:,source);
    
    data_d = avgstd_d_col(:,:,source)-avgdev_d_col(:,:,source);
    
    
    data_diff = data - data_d;
    
    
    
    %And now patients
    %bv
    
    data_pat = avgstd_pat_col(proc_pats_bv,:,source)-avgdev_pat_col(proc_pats_bv,:,source);
    
    data_d_pat = avgstd_d_pat_col(proc_pats_bv,:,source)-avgdev_d_pat_col(proc_pats_bv,:,source);
    
    
    data_diff_pat_bv = data_pat - data_d_pat;
    
    
    %PSP
    
    data_pat = avgstd_pat_col(proc_pats_psp,:,source)-avgdev_pat_col(proc_pats_psp,:,source);
    
    data_d_pat = avgstd_d_pat_col(proc_pats_psp,:,source)-avgdev_d_pat_col(proc_pats_psp,:,source);
    
    
    data_diff_pat_psp = data_pat - data_d_pat;
    
    
    %Calculate means
    
    m1=squeeze(nanmean(data_diff,1));
    
    s1=squeeze(nanstd(data_diff)./sqrt(size(data_diff,1)));
    
    
    m1_pat_bv = squeeze(nanmean(data_diff_pat_bv,1));

    s1_pat_bv = squeeze(nanstd(data_diff_pat_bv)./sqrt(size(data_diff_pat_bv,1)));
    
    
    m1_pat_psp = squeeze(nanmean(data_diff_pat_psp,1));

    s1_pat_psp = squeeze(nanstd(data_diff_pat_psp)./sqrt(size(data_diff_pat_psp,1)));
    
    
    
    %Start with subgroup plotting
    
    [h,p,ci,stats]=ttest2(squeeze(data_diff_pat_bv),squeeze(data_diff_pat_psp));
    k=max(max([m1_pat_bv(1,:); m1_pat_psp(1,:)]))+max(max([s1_pat_bv(1,:);s1_pat_psp(1,:)])); k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+1), subplot(2,3,source);
    
    figure(figcount+1); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_pat_bv(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1_pat_psp(1,:),'LineWidth',4,'Color',[0 255 0]./255); legendflex({'bvFTD','PSP'},'fontsize',10)
    boundedline([1:size(m1_pat_bv,2)],m1_pat_bv(1,:),s1_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat_psp,2)],m1_pat_psp(1,:),s1_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN PLA - MEM')
    title([grandavgdev.label(source) ' MMN DrugDiff'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    %Now including three groups (i.e. add controls)
    %Can go straight into analysis
        
    drugdiff_allgrps = cat(1, data_diff, data_diff_pat_bv, data_diff_pat_psp);
    
    grpst_allgrps = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp_bv 1]), repmat(3, [nfilesp_psp 1]));
    
    
    %Note, one-way anova now
    for t_win = 1:length(drugdiff_allgrps(1,:))
        p(1,t_win) = anova1(drugdiff_allgrps(:,t_win),grpst_allgrps, 'off');
    end
    
    k=max(max([m1; m1_pat_bv; m1_pat_psp]))+max(max([s1; s1_pat_bv; s1_pat_psp])); k=k*1.1;
    
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %Start figure
    
    figure(figcount+2), subplot(2,3,source);
    
    
    figure(figcount+2);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1,'LineWidth',4,'Color',[0 0 255]./255); hold on
    plot(m1_pat_bv,'LineWidth',4,'Color',[255 0 0]./255); 
    plot(m1_pat_psp,'LineWidth',4,'Color',[0/255 255/255 0/255]); legendflex({'CON', 'bvFTD','PSP'},'fontsize',10)
    
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat_bv,2)],m1_pat_bv(1,:),s1_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat_psp,2)],m1_pat_psp(1,:),s1_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN PLA - MEM')
    title([grandavgdev.label(source) ' MMN DrugDiff'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
     
end
    

%Save figures

saveas(figure(figcount+1), [FigOutDir '/' figoutprefix '_' 'LFP_'  'PatsSubgrp' '_Diff_DrugDiff_fullsample.tif']);

saveas(figure(figcount+2), [FigOutDir '/' figoutprefix '_' 'LFP_'  'ConssAndPatsSubgrp' '_Diff_DrugDiff_fullsample.tif']);
    


%Now rep3


MMNpall = [];


figcount = 16;


for source = 1:nsources
    
    
    
    %Start with controls
    
    
    data = avgstd3_col(:,:,source)-avgdev_col(:,:,source);
    
    data_d = avgstd3_d_col(:,:,source)-avgdev_d_col(:,:,source);
    
    
    data_diff = data - data_d;
    
    
    
    %And now patients
    %bv
    
    data_pat_bv = avgstd3_pat_col(proc_pats_bv,:,source)-avgdev_pat_col(proc_pats_bv,:,source);
    
    data_d_pat_bv = avgstd3_d_pat_col(proc_pats_bv,:,source)-avgdev_d_pat_col(proc_pats_bv,:,source);
    
    
    data_diff_pat_bv = data_pat_bv - data_d_pat_bv;
    
    
    %PSP
    
    data_pat_psp = avgstd3_pat_col(proc_pats_psp,:,source)-avgdev_pat_col(proc_pats_psp,:,source);
    
    data_d_pat_psp = avgstd3_d_pat_col(proc_pats_psp,:,source)-avgdev_d_pat_col(proc_pats_psp,:,source);
    
    
    data_diff_pat_psp = data_pat_psp - data_d_pat_psp;
    
    

    %Calculate means
    
    m1=squeeze(nanmean(data_diff,1));
    
    s1=squeeze(nanstd(data_diff)./sqrt(size(data_diff,1)));
    
    
    m1_pat_bv = squeeze(nanmean(data_diff_pat_bv,1));
    
    s1_pat_bv = squeeze(nanstd(data_diff_pat_bv)./sqrt(size(data_diff_pat_bv,1)));
    
    
    m1_pat_psp = squeeze(nanmean(data_diff_pat_psp,1));
    
    s1_pat_psp = squeeze(nanstd(data_diff_pat_psp)./sqrt(size(data_diff_pat_psp,1)));
    
    
    
    %Start with subgroup plotting
    
    [h,p,ci,stats]=ttest2(squeeze(data_diff_pat_bv),squeeze(data_diff_pat_psp));
    k=max(max([m1_pat_bv(1,:); m1_pat_psp(1,:)]))+max(max([s1_pat_bv(1,:);s1_pat_psp(1,:)])); k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %   Don't need to do the below
    %     temp_var = strcat( 'figure',num2str(source) );
    %
    %     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount+1), subplot(2,3,source);
    
    figure(figcount+1); set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    plot(m1_pat_bv(1,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1_pat_psp(1,:),'LineWidth',4,'Color',[0 255 0]./255); legendflex({'bvFTD','PSP'},'fontsize',10)
    boundedline([1:size(m1_pat_bv,2)],m1_pat_bv(1,:),s1_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat_psp,2)],m1_pat_psp(1,:),s1_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN3 PLA - MEM')
    title([grandavgdev.label(source) ' MMN3 DrugDiff'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    
    %Now including three groups (i.e. add controls)
    %Can go straight into analysis
    
    drugdiff_allgrps = cat(1, data_diff, data_diff_pat_bv, data_diff_pat_psp);
    
    grpst_allgrps = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp_bv 1]), repmat(3, [nfilesp_psp 1]));
    
    
    %Note, one-way anova now
    for t_win = 1:length(drugdiff_allgrps(1,:))
        p(1,t_win) = anova1(drugdiff_allgrps(:,t_win),grpst_allgrps, 'off');
    end
    
    k=max(max([m1; m1_pat_bv; m1_pat_psp]))+max(max([s1; s1_pat_bv; s1_pat_psp])); k=k*1.1;
    
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %Start figure
    
    figure(figcount+2), subplot(2,3,source);
    
    
    figure(figcount+2);set(gcf, 'Units','normalized', 'Position', [0 0 1 1]);
    plot(m1,'LineWidth',4,'Color',[0 0 255]./255); hold on
    plot(m1_pat_bv,'LineWidth',4,'Color',[255 0 0]./255);
    plot(m1_pat_psp,'LineWidth',4,'Color',[0/255 255/255 0/255]); legendflex({'CON', 'bvFTD','PSP'},'fontsize',10)
    
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat_bv,2)],m1_pat_bv(1,:),s1_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat_psp,2)],m1_pat_psp(1,:),s1_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'fontsize',10); xlabel('Time (ms)'); ylabel('MMN3 PLA - MEM')
    title([grandavgdev.label(source) ' MMN3 DrugDiff'], 'Fontsize',10); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %ANOVA
    
    %Calculate means across drug session
    
    
    m1_mmn = nanmean(data(:,MMNs:MMNf),2);
    
    m1_mmn_d = nanmean(data_d(:,MMNs:MMNf),2);
    
    
    m1_mmn_pat_bv = nanmean(data_pat_bv(:,MMNs:MMNf),2);
    
    m1_mmn_d_pat_bv = nanmean(data_d_pat_bv(:,MMNs:MMNf),2);
    
    
    m1_mmn_pat_psp = nanmean(data_pat_psp(:,MMNs:MMNf),2);
    
    m1_mmn_d_pat_psp = nanmean(data_d_pat_psp(:,MMNs:MMNf),2);
    
    
    
    %And their interaction significance
    
    meanMMNcol_all = cat(1, m1_mmn, m1_mmn_d, m1_mmn_pat_bv, m1_mmn_d_pat_bv, m1_mmn_pat_psp, m1_mmn_d_pat_psp);
    
    X = cat(2, meanMMNcol_all, F1all, F2all, Sall);
    
    
    [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
    
    
    pgrp = Ps(1);
    pcond = Ps(3);
    pint = Ps(4);
    
    MMNpall(source,1)=cell2mat(pgrp);
    MMNpall(source,2)=cell2mat(pcond);
    MMNpall(source,3)=cell2mat(pint);
    
    
end


%Save output tables of ANOVA
%MMN 

summrestable = table(MMNpall(:,1), MMNpall(:,2), MMNpall(:,3));
summrestable.Properties.VariableNames = {'Group' 'Drug' 'Int'};

writetable(summrestable, [FigOutDir '/' figoutprefix '_' 'LFP_'  'MMNmean' '_rep3_restable_2wayANOVA_ConsandPats_DrugInt_wDIAG.txt']);



%Save figures

saveas(figure(figcount+1), [FigOutDir '/' figoutprefix '_' 'LFP_' 'rep3_'  'PatsSubgrp' '_Diff_DrugDiff_fullsample.tif']);

saveas(figure(figcount+2), [FigOutDir '/' figoutprefix '_' 'LFP_' 'rep3_' 'ConssAndPatsSubgrp' '_Diff_DrugDiff_fullsample.tif']);
    


%% Wow, Finished

sprintf('finito, output saved to %s \n', FigOutDir)



end