function sourceLFPanalysis_plotstddevtrls_druggroup


%Open up SPM dirs (and field trip), and others

addpath('/imaging/rowe/users/ap09/Toolbox')

%for legend flex
addpath('/imaging/rowe/users/ap09/Toolbox/ekMEG/ntadscripts')

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest')

addpath('/imaging/rowe/users/ap09/Toolbox/fieldtrip')

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/boundedline-pkg'))

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/VBA-toolbox'))


ft_defaults()

FigOutDir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/analysis/SourceSpace/LFPs_icafixes/LFPs_COH_wrad_allMEGch_wfids/C_Plots';

figoutprefix = 'adv_ssst_newmaxf_fixICA_wfids_noc22_nodipcor_nop10p19';


%Other setup

nsources = 6;

sources={'RIFG', 'RSTG', 'RAUD'};


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

%2806 - remove error subj
%Con_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Con_D.txt';
Con_D = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/Con_D_noc22.txt';

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
    
    
   
    
    
end


for pdrugsubj = 1:length(DrugData_pat(:,1))
    
    
    filen = [preprocdir '/' DrugData_pat{pdrugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData_pat{pdrugsubj} '_' 's2' '_' LFPpostfix];
    
    filen_d = [preprocdir '/' DrugData_pat{pdrugsubj} '/' LFPprefix '_' preprocsteps '_' DrugData_pat{pdrugsubj} '_' 's1' '_' LFPpostfix];
    
    
    
    if isfile(filen) && isfile(filen_d)
        
        
        submatcount_pat = submatcount_pat + 1;
        
        allpats_bothsess{submatcount_pat,1} = DrugData_pat{pdrugsubj};
        
       
        
        
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
    
    

end


figcount = 0;

%% Start plots

%Only RH nodes

%Placebo first

for source = 4:nsources
    
    figcount = figcount + 1;
    
    
    %Start with controls   
    
    data = avgstd3_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgdev_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    %And patients
        
    data_pat = avgstd3_pat_col(:,:,source);
    
    m1_pat=squeeze(nanmean(data_pat,1));
    s1_pat=squeeze(nanstd(data_pat)./sqrt(size(data_pat,1)));
    
    data_d_pat = avgdev_pat_col(:,:,source);
    
    m1_d_pat=squeeze(nanmean(data_d_pat,1));
    s1_d_pat=squeeze(nanstd(data_d_pat)./sqrt(size(data_d_pat,1)));
    
    
    %Start figure
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 0.33 0.33]);
    plot(m1,'--','LineWidth',2.5,'Color',[0/255 0/255 255/255]); hold on
    plot(m1_d,':','LineWidth',2.5,'Color',[0/255 0/255 255/255]);
    plot(m1_pat,'--','LineWidth',2.5,'Color',[255/255 0/255 0/255]);
    plot(m1_d_pat,':','LineWidth',2.5,'Color',[255/255 0/255 0/255]); %legendflex({'rep3 CON','dev CON','rep3 FTLD','dev FTLD'},'fontsize',14)
    boundedline([1:size(m1,2)],m1,s1,'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1_d,s1_d,'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat,s1_pat,'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_d_pat,s1_d_pat,'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    %plot(p,'-o','LineWidth',2.5,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial'); ylabel({'Amplitude (a.u)'}, 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Arial');
    %title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
    box off
    
    set(gca, 'TickDir', 'Out')
    
    xlim([0 251]);
    
    
    %Save figure
    saveas(figure(figcount), ([FigOutDir '/' figoutprefix '_' 'LFP_'  'CONandPATS_' 'rep3anddevtrls_' 'PLA_' sources{figcount} '.tif']));    

    
end



%Now Drug

for source = 4:nsources
    
    figcount = figcount + 1;
    
    
    %Start with controls   
    
    data = avgstd3_d_col(:,:,source);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    data_d = avgdev_d_col(:,:,source);
    
    m1_d=squeeze(nanmean(data_d,1));
    s1_d=squeeze(nanstd(data_d)./sqrt(size(data_d,1)));
    
    
    %And patients
        
    data_pat = avgstd3_d_pat_col(:,:,source);
    
    m1_pat=squeeze(nanmean(data_pat,1));
    s1_pat=squeeze(nanstd(data_pat)./sqrt(size(data_pat,1)));
    
    data_d_pat = avgdev_d_pat_col(:,:,source);
    
    m1_d_pat=squeeze(nanmean(data_d_pat,1));
    s1_d_pat=squeeze(nanstd(data_d_pat)./sqrt(size(data_d_pat,1)));
    
    
    %Start figure
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 0.33 0.33]);
    plot(m1,'--','LineWidth',2.5,'Color',[0/255 0/255 255/255]); hold on
    plot(m1_d,':','LineWidth',2.5,'Color',[0/255 0/255 255/255]);
    plot(m1_pat,'--','LineWidth',2.5,'Color',[255/255 0/255 0/255]);
    plot(m1_d_pat,':','LineWidth',2.5,'Color',[255/255 0/255 0/255]); %legendflex({'rep3 CON','dev CON','rep3 FTLD','dev FTLD'},'fontsize',14)
    boundedline([1:size(m1,2)],m1,s1,'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1_d,s1_d,'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat,s1_pat,'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_d_pat,s1_d_pat,'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    %plot(p,'-o','LineWidth',2.5,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial'); ylabel({'Amplitude (a.u)'}, 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Arial');
    %title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
    box off
    
    set(gca, 'TickDir', 'Out')
    
    xlim([0 251]);
    
    
    %Save figure
    saveas(figure(figcount), ([FigOutDir '/' figoutprefix '_' 'LFP_'  'CONandPATS_' 'rep3anddevtrls_' 'MEM_' sources{figcount-3} '.tif']));    

    
end



%Now just fig legend

    figure(figcount+1);set(gcf, 'Units','normalized', 'Position', [0 0 0.33 0.33]);
    plot(m1,'--','LineWidth',2.5,'Color',[0/255 0/255 255/255]); hold on
    plot(m1_d,':','LineWidth',2.5,'Color',[0/255 0/255 255/255]);
    plot(m1_pat,'--','LineWidth',2.5,'Color',[255/255 0/255 0/255]);
    plot(m1_d_pat,':','LineWidth',2.5,'Color',[255/255 0/255 0/255]); legendflex({'rep3 CON','dev CON','rep3 FTLD','dev FTLD'},'fontsize',14)
    boundedline([1:size(m1,2)],m1,s1,'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1_d,s1_d,'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat,s1_pat,'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_d_pat,s1_d_pat,'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    %plot(p,'-o','LineWidth',2.5,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',14, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 14, 'FontName', 'Arial'); ylabel({'Amplitude (a.u)'}, 'FontWeight', 'bold', 'FontSize', 16, 'FontName', 'Arial');
    %title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
    box off
    
    set(gca, 'TickDir', 'Out')
    
    xlim([0 251]);
    
    
    %Save figure
    saveas(figure(figcount+1), ([FigOutDir '/' figoutprefix '_' 'LFP_'  'CONandPATS_' 'rep3anddevtrls_' 'MEM_' 'figleg' '.tif']));    




end