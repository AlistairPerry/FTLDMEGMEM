function preprocall_MMN_sourcespace
%% Preprocessing script for Roving MMN
% Rik Henson (2017), Modifications by Ece K (2020).
% Modifications indicated by EK

%% Processing paths

%Note, wont be able to find all dirs
addpath /imaging/local/meg_misc/
addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/;
addpath /hpc-software/matlab/cbu/;

%addpath /imaging/projects/cbu/ntad/scripts
%addpath (genpath('/imaging/projects/cbu/ntad/scripts/functions'))

% spm eeg;

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/fieldtrip/'); %EK

ft_defaults()

addpath(genpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/mne/')); %EK

%addpath(genpath('/imaging/projects/cbu/ntad/scripts/functions/FastICA_25/'))

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/osl-core'));
addpath(genpath('/imaging/rowe/users/ap09/Toolbox/MATLAB/FastICA_25'))

%And EK's MEG scrips
addpath('/imaging/rowe/users/ap09/Toolbox/ekMEG')


clear all; clc;


%% Setup
maxfiltdir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/preproc/Raw/B_Data';

origfname = 'rov123_sss.fif';

%Subjects
%Start with one for now
%subjects = {'c10'}; %debug epoching
%subjects = {'p22'}; %to determine if artifact rejection has worked..

%All
%subjects = {'c1','c10','c11','c12','c13','c14','c15','c16','c2','c21','c22','c23','c24','c3','c4','c5','c6','c7','c8','c9','p1','p10','p11','p12','p13','p14','p15','p16','p17','p18','p19','p2','p21','p22','p23','p24','p25','p3','p4','p5','p6','p7','p8','p9'};

%1404 - take out p3
subjects = {'c1','c10','c11','c12','c13','c14','c15','c16','c2','c21','c22','c23','c24','c3','c4','c5','c6','c7','c8','c9','p1','p10','p11','p12','p13','p14','p15','p16','p17','p18','p19','p2','p21','p22','p23','p24','p25','p3','p4','p5','p6','p7','p8','p9'};
% subjects = {'p3'}; %prob subject - something wrong with p3 (session 2)
% subjects = {'p4','p5','p6','p7','p8','p9'}; %ran again because of aforementioned subject session

%290622 - take out c22
subjects = {'c1','c10','c11','c12','c13','c14','c15','c16','c2','c21','c23','c24','c3','c4','c5','c6','c7','c8','c9','p1','p10','p11','p12','p13','p14','p15','p16','p17','p18','p19','p2','p21','p22','p23','p24','p25','p3','p4','p5','p6','p7','p8','p9'};

sessions = {'s1', 's2'};

% 110221/missing subject
%subjects = {'p10'};
%sessions = {'s1'};

preprocdir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/preproc/preproc/Source/LFPs/B_Data/Adv_fix_icafixes';


%Ensure that preproc directory exitsts
mkdir(preprocdir);


% Directories
% ana_dir = '/imaging/projects/cbu/ntad/meg_data/'; %AP - change to preprocdir
% task='mmn';
% session = 'TW';


% Processing switches
doconvert = 1; % Convert to SPM
dodown = 1; % Downsample
dofilter = 1; % Filter
donotch = 1; % Notch filter
doartefact = 1; % Artefact rejection,
doica = 1; % Run ICA
doepoch = 1; % Epoch
doartefact2 = 1; % Second artefact rejection (not sure about this, will need to ask again) - actually this is bad trials being removed :)
doaverage = 0; % Average trials and low pass
docontrast = 0; % Compute contrasts
dofilter2 = 0; % Filter to analysis-specific frequency band
dorename = 0; % Rename files to something sensible

% Parameters
basefile = 'transdef'; %Maxfilter output to preprocess; %AP - ok need to be careful what fname its picking up
pads = [0 0]; % extra padding to allow lowpass filtering, for visualisation only [-100 100]
epochWin = [-100 400]; % epoch window in ms

%con_values = [100 103 106]; % trigger values that correspond to your conditions of interest %AP to check
%AP! - need to fix

offset = [64 64 64]; % stimulus delivery delay for each condition (NI card & Presentation; Audio: 64 ms; %AP to check Vis: 35 ms)

%con_labels = {'DEV','REP3','REP6'}; % condition names in the same order as
%trigger values - AP to change
con_labels = {'DEV','REP5'}; % condition names in the same order as trigger values

filtwindow = [0.1 125]; % Filter window
filtwindow2 = [0.1 45]; % Analysis specific frequency range
drate = 500; % Downsampling rate
PCdim = 60; % Number of PCs for ICA
sampletime=10; % time segment to visualise in IC time series
samprate=1000; %sampling rate in Hz

%ref_chans = {'VEOG', 'HEOG', 'ECG'}; % external references for ICA %AP - to change (below as well)
ref_chans = {'VEOG', 'HEOG'}; % external references for ICA %AP - to change (below as well)

mods = {'MEGMAG','MEGPLANAR','EEG'}; % modalities

%ses_nam=session; AP remove

%%=========================================================================
%% Setting up Matlab Pool
%%=========================================================================

%parpool(15)

%It has been turned off now, apparently, but can be easily switched on

%Related prob is to be able to run all subjects, when they dont have a
%sess..

%%=========================================================================
%% Convert Data to SPM
%%=========================================================================

if doconvert
    chanfile = fullfile('chan_select_MEG_EEG_STI101.mat'); %EK
    %AP - unsure what to do here
    % Channels to read: MEG+3EEG+EOG+ECG+STI101 (excludes MISC for subjects with no eye-tracking)
    chan = load(chanfile);
    
    %ref_old_names = {'EOG061','EOG062','ECG063'};
    %ref_new_names = {'HEOG';'VEOG';'ECG'}; %AP - take out ECG comp
    ref_old_names = {'EOG061','EOG062'};
    ref_new_names = {'HEOG';'VEOG'}; %AP - take out ECG comp
    
    %     parfor sub=1:length(subjects)
    for sub = 1:length(subjects)
        
        %ensure subj processing dir exists
        mkdir([preprocdir '/' subjects{sub}])
        
        for sess = 1:length(sessions)
            
            mkdir([preprocdir '/' subjects{sub}], sessions{sess});
            
            
            if isfile(fullfile([maxfiltdir,'/',subjects{sub},'/', sessions{sess}, '/', subjects{sub} '_' sessions{sess} '_' origfname]))
                
                %sub_dir = fullfile(ana_dir,subjects{sub},session,task);
                sub_dir = fullfile(maxfiltdir,subjects{sub},sessions{sess}); %AP - both sessions located in one folder
                
                all_trl{sub} = {}; all_con{sub} = {}; %AP! - revisit this
                
                %ses_out_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/'];
                ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
                
                mkdir(ses_out_dir)
                
                %outfile = [ses_dir basefile '.mat'];
                basefile = [subjects{sub} '_' sessions{sess} '_rov123_sss'];
                outfile = [ses_out_dir '/' 'S_' basefile '.mat'];
                
                
                S = [];
                S.outfile = outfile;
                
                %S.dataset = fullfile(ses_dir,sprintf('%s.fif',basefile));
                S.dataset = fullfile(sub_dir,sprintf('%s.fif',basefile));
                
                
                S.save = 0;  % saved below anywapadsy (if S.save=1, then prompts for new filename!)
                S.reviewtrials = 0;
                S.channels = chan.label;
                S.continuous = 1;pads
                S.checkboundary = 0;
                
                D = spm_eeg_convert(S);
                
                pads
                tc = find(strcmp(D.chanlabels,'STI101'));
                D(tc,:) = D(tc,:)/1e6;  % 1e6 new scaling introduced by SPM12!
                D = units(D,tc,'V');pads
                
                for c = 1:length(ref_old_names)
                    ch = indchannel(D,ref_old_names{c});
                    D = chanlabels(D,ch,ref_new_names{c});
                end
                
                D.save;
                
            end
            
        end
        
    end
    
end

%%=========================================================================
%% Downsample EK
%%=========================================================================

%AP -ok, straight forward

if dodown
    %     parfor sub=1:length(subjects)
    for sub=1:length(subjects)
        
        for sess = 1:length(sessions)
            
        ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];

            
         if isfile(fullfile([ses_out_dir '/' 'S_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
             
                
                basefile = ['S_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                
                
                %ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task,'/')
                
                S = [];
                
                %S.D = [ses_dir '/' basefile '.mat'];
                
                S.D = [ses_out_dir '/' basefile '.mat'];
                
                S.method = 'downsample';
                S.fsample_new = drate;
                
                D = spm_eeg_downsample(S);
                
            end
            
        end
        
    end
    
end

%%=========================================================================
%% Filter
%%=========================================================================

if dofilter
    
    low=filtwindow(2); high=filtwindow(1);
    
    %     parfor sub=1:length(subjects)
    for sub=1:length(subjects)
        
        for sess = 1:length(sessions)
            
             
            ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];

            
            if isfile(fullfile([ses_out_dir '/' 'dS_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
                
                
                %ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task)
                
                basefile = ['dS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                
                
                S = [];
                
                S.D = [ses_out_dir '/' basefile '.mat'];
                S.type = 'butterworth';
                S.order = 5;
                S.band = 'low';
                S.freq = low;
                D = spm_eeg_filter(S);
                
                S= [];
                S.D = D;
                S.type = 'butterworth';
                S.order = 5;
                S.band = 'high';
                S.freq = high;
                D = spm_eeg_filter(S);
                
            end
            
        end
        
    end
    
end


%%=========================================================================
%% Notch filter EK
%%=========================================================================

if donotch
    
    %    parfor sub=1:length(subjects)
    for sub=1:length(subjects)
        
        for sess = 1:length(sessions)
            
        ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];

            
            if isfile(fullfile([ses_out_dir '/' 'dS_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
                
                
                %ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task)
                
                basefile = ['ffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                
                S = [];
                
                %S.D = [ses_dir '/ffd'  basefile '.mat'];
                S.D = [ses_out_dir '/'  basefile '.mat'];
                
                S.type = 'butterworth';
                S.order = 5;
                S.band = 'stop';
                S.freq = [45 55];
                
                D = spm_eeg_filter(S);
                
                S = [];
                S.D = D;
                S.type = 'butterworth';
                S.order = 5;
                S.band = 'stop';
                S.freq = [95 105];
                
                D = spm_eeg_filter(S);
                
            end
            
        end
        
    end
    
end

%%=========================================================================
%% Artefact rejection EK
%%=========================================================================

if doartefact % uses Oxford's OSL, works much better than SPM's.
    
    for sub=1:length(subjects)
        
        for sess = 1:length(sessions)
            
                            ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];

            
            if isfile(fullfile([ses_out_dir '/' 'ffffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
                
                
                %ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
                
                basefile = ['ffffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                
                %D = spm_eeg_load([ses_dir '/ffffd' basefile '.mat']); %EK
                D = spm_eeg_load([ses_out_dir '/' basefile '.mat']);
                
                D=D.badchannels(1:D.nchannels,0);
                D2=osl_detect_artefacts(D);D=D2;
                
                %save([ses_dir '/affffd' basefile '.mat'],'D','-v7.3')
                save([ses_out_dir '/' 'a' basefile '.mat'],'D','-v7.3')
                
            end
            
        end
        
    end
    
end

%%=========================================================================
%% ICA %EK
%%=========================================================================

submaticacount = 0;

if doica
    
    %     ref_old_names = {'EOG061','EOG062','ECG063'};
    %     ref_new_names = {'HEOG';'VEOG';'ECG'}; %AP -again take out
    ref_old_names = {'EOG061','EOG062'};
    ref_new_names = {'HEOG';'VEOG'};
    
    ICA = [];
    ICA.PCA_dim = 60; % Number of PCs for ICA
    ICA.Nperm = 0; %round((4*sqrt(CorPval)/CorPval)^2); % If want accurate p-values
    ICA.Rseed = 1; % to make reproducible
    all_ica_remove = {}; all_ica_TraMat = {}; all_ica_out = {};
    
    %arttopos = load('/imaging/rh01/Methods/MEGEEG70ArtifactTemplateTopographies');
    arttopos = load('/imaging/rowe/users/ek01/collab/Alistair_FTD_MEM/MEGEEG70ArtifactTemplateTopographies');

    
    sparefs = {};
    mind = find(ismember({'MEGMAG', 'MEGPLANAR', 'EEG'},mods));
    
    for m=1:length(mods)
        sparefs{m} = {};
        for r = 1:length(ref_chans)
            tmp = getfield(arttopos,upper(ref_chans{r}));  % If doing blinks only
            sparefs{m}{r} = tmp{mind(m)}';
        end
    end
    
    %     parfor sub=1:length(subjects)
    for sub=1:length(subjects)
        
        for sess = 1:length(sessions)
            
          
            ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
            
            
            if isfile(fullfile([ses_out_dir '/' 'affffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
                
                submaticacount = submaticacount+1; %AP add with icafix
             
                
                basefile = ['affffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                
                
                %ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task);
                cd(ses_out_dir);
                
                %D=spm_eeg_load(['affffd' basefile '.mat']);
                D=spm_eeg_load([basefile '.mat']);
                
                disp([subjects{sub} ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++'])
                
                %outfile = [ses_dir '/affffd' basefile '.mat'];
                outfile = [ses_out_dir '/' basefile '.mat']; %bit weird to specify it as output, but lets see
                
                if size(D,1)>330
                    mod={'MEGMAG', 'MEGPLANAR', 'EEG'};
                else
                    mod={'MEGMAG', 'MEGPLANAR'};
                end;
                
                %icafile = fullfile(ses_dir,sprintf('%s.mat',['Maffffd' basefile]));
                icafile = fullfile(ses_out_dir,sprintf('%s.mat',['M' basefile]));
                
                
                %AP -take out ECG stuff
                S = ICA;
                
                S.refs.tem = {}; % Reference signal for correlating with ICs
                if ~isempty(find(ismember(D.chanlabels,'VEOG'))) %EK
                    S.refs.tem{1}=D(find(ismember(D.chanlabels,'VEOG')),:); %VEOG
                    S.refs.tem{2}=D(find(ismember(D.chanlabels,'HEOG')),:); %VEOG
                    %S.refs.tem{3}=D(find(ismember(D.chanlabels,'ECG')),:); %VEOG
                else
                    S.refs.tem{1}=D(find(ismember(D.chanlabels,'EOG062')),:); %VEOG
                    S.refs.tem{2}=D(find(ismember(D.chanlabels,'EOG061')),:); %VEOG
                    %S.refs.tem{3}=D(find(ismember(D.chanlabels,'ECG063')),:); %VEOG
                end
                
                chans = {}; remove = {}; weights = {}; temcor = {}; spacor = {}; TraMat = {};
                
                for m = 1:length(mod)
                    
                    S.refs.spa = sparefs{m};
                    chans = indchantype(D,mod{m});
                    S.d  = D(chans,:);
                    [Out,ICs] = detect_ICA_artefacts_ek(S);
                    Out.chans = chans;
                    %ses_nam = session;
                    icafile = fullfile(ses_out_dir,sprintf('%s_ica.mat',['affffdS_' mod{m}]));
                    all_subs_ica_out{submaticacount}{m}=Out;
                    
                end
                
            end
            
        end
        
    end
    
    % Save ICA output
    %ses_nam=session; AP - I have two sessions
    
    
    submaticacount = 0; %AP add with icafix
    
    for sub = 1:length(subjects)
        
        for sess = 1:length(sessions)
            
            ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
            
            
            if isfile(fullfile([ses_out_dir '/' 'affffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
                
                
                submaticacount = submaticacount + 1; %AP add with icafix
                
                     
                basefile = ['affffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                
                %ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task,'/'); cd(ses_dir)
                %outfile = [ses_dir '/affffd' basefile '.mat'];
                
                outfile = [ses_out_dir '/' basefile '.mat'];
                
                D = spm_eeg_load(outfile);
                
                if size(D,1)>330 %EK
                    mod={'MEGMAG', 'MEGPLANAR', 'EEG'};
                else
                    mod={'MEGMAG', 'MEGPLANAR'};
                end
                
                for m = 1:length(mod)
                    
                    %icafile = fullfile(ses_Outdir,sprintf('%s_ica.mat',['affffd' basefile '_' mod{m}]));
                    icafile = fullfile(ses_out_dir,sprintf('%s_ica.mat',[basefile '_' mod{m}]));
                    Out=all_subs_ica_out{submaticacount}{m};
                    
                    %parsave(icafile,Out);
                    %save(icafile,Out); AP - apparently only likes ''
                    
                    save(icafile,'Out');
                    
                end
                
            end
            
        end
        
    end
    
    %     parfor sub=1:length(subjects)
    for sub=1:length(subjects)
        
        for sess = 1:length(sessions)
            
            
            ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
            
            
            if isfile(fullfile([ses_out_dir '/' 'affffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
                    
                ses_nam = sessions{sess};
                
                %ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task);
                basefile = ['affffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                
                TraMat = {}; chans = [];
                %D=spm_eeg_load([ses_dir '/affffd' basefile '.mat']); %original file - AP indent out
                D=spm_eeg_load([ses_out_dir '/' basefile '.mat']); %original file
                
                if size(D,1)>330
                    mods={'MEGMAG', 'MEGPLANAR', 'EEG'};
                else
                    mods={'MEGMAG', 'MEGPLANAR'};
                end
                
                for m = 1:length(mods)
                    %icafile = fullfile(ses_dir,sprintf('%s_ica.mat',['affffd' basefile '_' mods{m}]));
                    icafile = fullfile(ses_out_dir,sprintf('%s_ica.mat',[basefile '_' mods{m}]));
                    
                    ica = load(icafile); %ica output
                    TraMat{m} = ica.Out.TraMat;
                    chans = [chans ica.Out.chans];
                end
                
                if sum(TraMat{1,1}(:))~=size(TraMat{1,1},1)
                    S = []; S.D = D;
                    for c = 1:length(chans)
                        S.montage.labelorg{c} = D.chanlabels{chans(c)};
                    end
                    S.montage.labelnew = S.montage.labelorg;
                    S.montage.tra      = blkdiag(TraMat{:});
                    S.keepothers       = 1;
                    D = spm_eeg_montage(S); % becomes 'Mfftransdef.mat' GOODIES
                    
                else
                    
                    S=[]; S.D = D;
                    %S.outfile = [ses_dir '/Maffffd' basefile '.mat'];
                    S.outfile = [ses_out_dir '/' 'M' basefile '.mat'];
                    
                    D = spm_eeg_copy(S);
                    
                end
                
            end
            
        end
        
    end
    
    % Plot removed ICs EK
    for sub=1:length(subjects)
        
        for sess = 1:length(sessions)
            
            ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
            
            
            if isfile(fullfile([ses_out_dir '/' 'MaffffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
                
                
                ses_nam = sessions{sess};
                
                %ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task);
                basefile = ['affffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                
                
                %D=spm_eeg_load([ses_dir '/affffd' basefile '.mat']); %original file
                D=spm_eeg_load([ses_out_dir '/' basefile '.mat']); %original file
                
                if size(D,1)>330; mod={'MEGMAG','MEGPLANAR','EEG'}; else; mod={'MEGMAG','MEGPLANAR'};end
                
                for m=1:length(mod)
                    
                    %load([ses_dir '/affffd' basefile '_' mod{m} '_ica.mat']);
                    
                    load([ses_out_dir '/' 'affffdS' '_' subjects{sub} '_' sessions{sess} '_rov123_sss_' mod{m} '_ica.mat']); %bit confused here as its pointing to pre-ICA file
                    w=pinv(Out.weights)'; tcorrs=Out.temcor;scorrs=Out.spacor;
                    idx=find(strcmp(D.chantype,mod{m}));
                    in.noButtons=1;in.type=mod{m};in.plotpos=0;
                    %w = weights(idx,:)';
                    close all;
                    rem=unique([Out.temrem, Out.sparem]);
                    
                    for com=1:length(rem)
                        bc=rem(com);
                        figure
                        set(gca,'FontSize',9)
                        [ZI,f] =spm_eeg_plotScalpData_ek2(w(bc,:)',D.coor2D(idx),D.chanlabels(idx),in);
                        close(figure(2));
                        imagesc(flipdim(ZI,1));
                        colordata=colormap;cc=colordata(1,:);
                        colordata(4:67,:)=colormap;colordata(2:3,:)=repmat(cc,2,1);colordata(1,:)=[1 1 1];
                        colormap(colordata); clear colordata cc
                        axis off;
                        
                        %AP - again take out below stuff
                        if tcorrs(1,bc)>=0.3 || scorrs(1,bc)>=0.7
                            txt=['HEOG_tem_' num2str(tcorrs(1,bc),'%.2f') '_spa_' num2str(scorrs(1,bc),'%.2f')];
                            
                        else tcorrs(2,bc)>=0.3 || scorrs(2,bc)>=0.7
                            txt=['VEOG_tem_' num2str(tcorrs(2,bc),'%.2f') '_spa_' num2str(scorrs(2,bc),'%.2f')];
                            %                     elseif tcorrs(2,bc)>=0.3 || scorrs(2,bc)>=0.7
                            %                         txt=['VEOG_tem_' num2str(tcorrs(2,bc),'%.2f') '_spa_' num2str(scorrs(2,bc),'%.2f')];
                            %                     elseif tcorrs(3,bc)>=0.3 || scorrs(3,bc)>=0.7
                            %                         txt=['ECG_tem_' num2str(tcorrs(3,bc),'%.2f') '_spa_' num2str(scorrs(3,bc),'%.2f')];
                        end
                        
                        title({[ 'IC' num2str(bc) ' ' strrep(txt,'_',' ') ' ' mod{m}]});
                        %print(figure(1),'-dtiff', [ses_dir '/' mod{m} '_ic' num2str(bc) '_' subjects{sub} '_' txt '.tiff'],'-r200');
                        
                        print(figure(1),'-dtiff', [ses_out_dir '/' mod{m} '_ic' num2str(bc) '_' subjects{sub} '_' txt '.tiff'],'-r200');
                        
                        close(figure(1));
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

%%=========================================================================
%% Read triggers & Epoch
%%=========================================================================
%Could be most difficult bit - trying to read in triggers etc and epoching
%to trial conditions
%AP maybe I could just insert scripts from previous scripts

if doepoch
    
    %Indent all out so you can subsitute your own epoching script
    %     all_trl = {}; all_con = {};
    %     padding = pads;
    %     ses_nam = session;
    %
    % %     parfor sub=1:length(subjects)
    %     for sub=1:length(subjects)
    %         if isfile(fullfile(ana_dir,subjects{sub},ses_nam,'/',task),'dir')
    %
    %             ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/'];
    %             if ~isfile([ses_dir 'trial_info.txt'],'file')
    %                 filelog=fopen([ses_dir 'trial_info.txt'],'w');
    %                 fclose(filelog);%EK
    %             end
    %
    %             filelog=fopen([ses_dir 'trial_info.txt'],'a');%EK
    %             D = spm_eeg_load([ses_dir 'dtransdef.mat']); %change dtsss to dtransdef
    %             d = D(indchannel(D,'STI101'),:);
    %             if max(d(:))>5000; d=d/1000000; end; d=round(d); %EK
    %             if isequal(unique(d),0); d = D(indchannel(D,'STI101'),:); end
    %             fprintf('unique trigger codes: %s\n',num2str(unique(d)))  % triggers found
    %             fprintf(filelog, 'unique trigger codes: %s\n',num2str(unique(d)));
    %
    %             dd = [0 diff(d)]; % onset of triggers
    %             ons = find(ismember(dd,con_values));
    %             trl = []; conditionlabels = {};
    %             for o = 1:length(ons)
    %                 con = find(con_values == dd(ons(o)));
    %                 conditionlabels{end+1} = con_labels{con};
    %                 trl(end+1,1) = ons(o) + round((epochWin(1) + padding(1) + offset(con))*(D.fsample/1000));
    %                 trl(end,2)   = ons(o) + round((epochWin(2) + padding(2) + offset(con))*(D.fsample/1000));
    %                 trl(end,3)   = round((epochWin(1) + padding(1))*(D.fsample/1000));
    %             end
    %
    %             fprintf('%s\n', subjects{sub});
    %             for con = 1:length(con_values)
    %                 fprintf('%d \ttriggers of type %s\n',length(find(strcmp(conditionlabels,con_labels{con}))),con_labels{con});
    %                 fprintf(filelog, '%d \ttriggers of type %s\n',length(find(strcmp(conditionlabels,con_labels{con}))),con_labels{con}); %EK
    %             end
    %
    %             fclose(filelog);%EK
    %             all_trl{sub} = trl; all_con{sub} = conditionlabels;
    %         end
    %     end
    %
    %     for sub = 1:length(subjects)
    %         if isfile(fullfile(ana_dir,subjects{sub},ses_nam,task))
    %             ses_nam = session;
    %             ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task);
    %             trl = all_trl{sub}; conditionlabels = all_con{sub};
    %             trlfile = fullfile(ses_dir,sprintf('%s_trl.mat','tsss')); %change tsss to transdef
    %             save(trlfile,'trl','conditionlabels');
    %         end
    %     end
    %
    %     bline=epochWin(1);
    % %         parfor sub=1:length(subjects)
    %     for sub=1:length(subjects)
    %
    %         ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/'];
    %         trl = load(fullfile(ses_dir,sprintf('%s_trl.mat',basefile)));
    %         %trl = load(fullfile(ses_dir,sprintf('%s_trl.mat',basefile(1:4))));
    %         D = spm_eeg_load([ses_dir '/Maffffd' basefile '.mat']); %EK
    %         S = []; S.D = D;
    %         S.trl = trl.trl;
    %         S.conditionlabels = trl.conditionlabels;
    %         S.eventpadding = 0;  % Already included above!
    %         S.bc = 0;            % No need to baseline correct - done below...
    %         D = spm_eeg_epochs(S);
    %
    %         % Baseline correct (ie ignoring padding)
    %         S = []; S.D = D;
    %         S.timewin = [bline 0];
    %         S.save = 1;
    %         S.prefix='';
    %         D = spm_eeg_bc(S);
    %
    %     end
    
    for sub=1:length(subjects)
        
        if strcmp(subjects{sub}, 'p3') %do this because p3 has s2 has epoching issues
            
            session='s1';
            
            ses_out_dir = [preprocdir '/' subjects{sub} '/' session];
            
            
            if isfile(fullfile([ses_out_dir '/' 'MaffffdS_' subjects{sub} '_' session '_rov123_sss.mat']))
                
                
                basefile = ['MaffffdS_' subjects{sub} '_' session '_rov123_sss'];
                
                
                D = spm_eeg_load([ses_out_dir '/' basefile '.mat']);
                
                i = find(strcmp(D.chanlabels,'STI101'));
                vals = [10, 30, 50, 70, 90, 110, 130, 150, 170]; % with +1 for reps % each of these are different frequencies - so first dev = 10, the second = 11, etc, for the first frequency
                D(i,:,:) = D(i,:,:)/10^7;
                
                % prepare for spm trial define:
                S.D        = D;
                S.timewin  = [-100 400];
                S.save     = 0;
                S.reviewtrials = 0;
                
                nn     = 0;
                for n  = 1:10
                    if n == 1
                        nm = 'Dev';
                    else nm = sprintf('rep%d',(n-1));
                    end
                    for v = 1:length(vals)
                        nn = nn +1;
                        S.trialdef(nn).conditionlabel = nm;
                        S.trialdef(nn).eventtype      = 'STI101_up';
                        S.trialdef(nn).eventvalue     = (vals(v)+(n-1));
                    end
                end
                
                [trl,conditionlabels,S2] = spm_eeg_definetrial(S);
                
                S = [];
                S.D = D;
                S.bc = 1;
                S.trl = trl;
                S.conditionlabels = conditionlabels;
                
                D2 = spm_eeg_epochs(S);
                
                % save D2:
                D2.save;
                clear D;
                D = D2;
                
            end
            
        else
            
            for sess = 1:length(sessions)
                
                
                ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
                
                
                if isfile(fullfile([ses_out_dir '/' 'MaffffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
                    
                    
                    basefile = ['MaffffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                    
                    
                    D = spm_eeg_load([ses_out_dir '/' basefile '.mat']);
                    
                    i = find(strcmp(D.chanlabels,'STI101'));
                    vals = [10, 30, 50, 70, 90, 110, 130, 150, 170]; % with +1 for reps % each of these are different frequencies - so first dev = 10, the second = 11, etc, for the first frequency
                    D(i,:,:) = D(i,:,:)/10^7;
                    
                    % prepare for spm trial define:
                    S.D        = D;
                    S.timewin  = [-100 400];
                    S.save     = 0;
                    S.reviewtrials = 0;
                    
                    nn     = 0;
                    for n  = 1:10
                        if n == 1
                            nm = 'Dev';
                        else nm = sprintf('rep%d',(n-1));
                        end
                        for v = 1:length(vals)
                            nn = nn +1;
                            S.trialdef(nn).conditionlabel = nm;
                            S.trialdef(nn).eventtype      = 'STI101_up';
                            S.trialdef(nn).eventvalue     = (vals(v)+(n-1));
                        end
                    end
                    
                    [trl,conditionlabels,S2] = spm_eeg_definetrial(S);
                    
                    S = [];
                    S.D = D;
                    S.bc = 1;
                    S.trl = trl;
                    S.conditionlabels = conditionlabels;
                    
                    D2 = spm_eeg_epochs(S);
                    
                    % save D2:
                    D2.save;
                    clear D;
                    D = D2;
                    
                end
                
            end
            
        end
        
    end
    
end

%%=========================================================================
%% Artefact rejection and removal
%%=========================================================================

if doartefact2 % EK
    
    for sub=1:length(subjects)
        
        for sess=1:length(sessions)
  
            ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
            
            
            if isfile(fullfile([ses_out_dir '/' 'eMaffffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss.mat']))
                
                %         ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
                %         D = spm_eeg_load([ses_dir '/eMaffffd' basefile '.mat']); %EK
                
                basefile = ['eMaffffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
                
                ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
                
                D = spm_eeg_load([ses_out_dir '/' basefile '.mat']); %AP
                
                if ~isequal(max(ismember(D.chantype,'MEGMAG')),1)
                    %load([ses_dir '/eMaffffd' basefile '.mat']); AP
                    load([ses_out_dir '/' basefile '.mat'])
                    
                    %AP -remove ECG below
                    %chans={'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211';'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';'MEG0331';'MEG0341';'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531';'MEG0541';'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721';'MEG0731';'MEG0741';'MEG0811';'MEG0821';'MEG0911';'MEG0921';'MEG0931';'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121';'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311';'MEG1321';'MEG1331';'MEG1341';'MEG1411';'MEG1421';'MEG1431';'MEG1441';'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';'MEG1631';'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821';'MEG1831';'MEG1841';'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011';'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';'MEG2131';'MEG2141';'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331';'MEG2341';'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521';'MEG2531';'MEG2541';'MEG2611';'MEG2621';'MEG2631';'MEG2641';'MEG0113';'MEG0112';'MEG0122';'MEG0123';'MEG0132';'MEG0133';'MEG0143';'MEG0142';'MEG0213';'MEG0212';'MEG0222';'MEG0223';'MEG0232';'MEG0233';'MEG0243';'MEG0242';'MEG0313';'MEG0312';'MEG0322';'MEG0323';'MEG0333';'MEG0332';'MEG0343';'MEG0342';'MEG0413';'MEG0412';'MEG0422';'MEG0423';'MEG0432';'MEG0433';'MEG0443';'MEG0442';'MEG0513';'MEG0512';'MEG0523';'MEG0522';'MEG0532';'MEG0533';'MEG0542';'MEG0543';'MEG0613';'MEG0612';'MEG0622';'MEG0623';'MEG0633';'MEG0632';'MEG0642';'MEG0643';'MEG0713';'MEG0712';'MEG0723';'MEG0722';'MEG0733';'MEG0732';'MEG0743';'MEG0742';'MEG0813';'MEG0812';'MEG0822';'MEG0823';'MEG0913';'MEG0912';'MEG0923';'MEG0922';'MEG0932';'MEG0933';'MEG0942';'MEG0943';'MEG1013';'MEG1012';'MEG1023';'MEG1022';'MEG1032';'MEG1033';'MEG1043';'MEG1042';'MEG1112';'MEG1113';'MEG1123';'MEG1122';'MEG1133';'MEG1132';'MEG1142';'MEG1143';'MEG1213';'MEG1212';'MEG1223';'MEG1222';'MEG1232';'MEG1233';'MEG1243';'MEG1242';'MEG1312';'MEG1313';'MEG1323';'MEG1322';'MEG1333';'MEG1332';'MEG1342';'MEG1343';'MEG1412';'MEG1413';'MEG1423';'MEG1422';'MEG1433';'MEG1432';'MEG1442';'MEG1443';'MEG1512';'MEG1513';'MEG1522';'MEG1523';'MEG1533';'MEG1532';'MEG1543';'MEG1542';'MEG1613';'MEG1612';'MEG1622';'MEG1623';'MEG1632';'MEG1633';'MEG1643';'MEG1642';'MEG1713';'MEG1712';'MEG1722';'MEG1723';'MEG1732';'MEG1733';'MEG1743';'MEG1742';'MEG1813';'MEG1812';'MEG1822';'MEG1823';'MEG1832';'MEG1833';'MEG1843';'MEG1842';'MEG1912';'MEG1913';'MEG1923';'MEG1922';'MEG1932';'MEG1933';'MEG1943';'MEG1942';'MEG2013';'MEG2012';'MEG2023';'MEG2022';'MEG2032';'MEG2033';'MEG2042';'MEG2043';'MEG2113';'MEG2112';'MEG2122';'MEG2123';'MEG2133';'MEG2132';'MEG2143';'MEG2142';'MEG2212';'MEG2213';'MEG2223';'MEG2222';'MEG2233';'MEG2232';'MEG2242';'MEG2243';'MEG2312';'MEG2313';'MEG2323';'MEG2322';'MEG2332';'MEG2333';'MEG2343';'MEG2342';'MEG2412';'MEG2413';'MEG2423';'MEG2422';'MEG2433';'MEG2432';'MEG2442';'MEG2443';'MEG2512';'MEG2513';'MEG2522';'MEG2523';'MEG2533';'MEG2532';'MEG2543';'MEG2542';'MEG2612';'MEG2613';'MEG2623';'MEG2622';'MEG2633';'MEG2632';'MEG2642';'MEG2643';'EEG001';'EEG002';'EEG003';'EEG004';'EEG005';'EEG006';'EEG007';'EEG008';'EEG009';'EEG010';'EEG011';'EEG012';'EEG013';'EEG014';'EEG015';'EEG016';'EEG017';'EEG018';'EEG019';'EEG020';'EEG021';'EEG022';'EEG023';'EEG024';'EEG025';'EEG026';'EEG027';'EEG028';'EEG029';'EEG030';'EEG031';'EEG032';'EEG033';'EEG034';'EEG035';'EEG036';'EEG037';'EEG038';'EEG039';'EEG040';'EEG041';'EEG042';'EEG043';'EEG044';'EEG045';'EEG046';'EEG047';'EEG048';'EEG049';'EEG050';'EEG051';'EEG052';'EEG053';'EEG054';'EEG055';'EEG056';'EEG057';'EEG058';'EEG059';'EEG060';'EEG065';'EEG066';'EEG067';'EEG068';'EEG069';'EEG070';'EEG071';'EEG072';'EEG073';'EEG074';'HEOG';'VEOG';'ECG';'STI101'};
                    %types={'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EOG';'EOG';'ECG';'Other'};
                    
                    chans={'MEG0111';'MEG0121';'MEG0131';'MEG0141';'MEG0211';'MEG0221';'MEG0231';'MEG0241';'MEG0311';'MEG0321';'MEG0331';'MEG0341';'MEG0411';'MEG0421';'MEG0431';'MEG0441';'MEG0511';'MEG0521';'MEG0531';'MEG0541';'MEG0611';'MEG0621';'MEG0631';'MEG0641';'MEG0711';'MEG0721';'MEG0731';'MEG0741';'MEG0811';'MEG0821';'MEG0911';'MEG0921';'MEG0931';'MEG0941';'MEG1011';'MEG1021';'MEG1031';'MEG1041';'MEG1111';'MEG1121';'MEG1131';'MEG1141';'MEG1211';'MEG1221';'MEG1231';'MEG1241';'MEG1311';'MEG1321';'MEG1331';'MEG1341';'MEG1411';'MEG1421';'MEG1431';'MEG1441';'MEG1511';'MEG1521';'MEG1531';'MEG1541';'MEG1611';'MEG1621';'MEG1631';'MEG1641';'MEG1711';'MEG1721';'MEG1731';'MEG1741';'MEG1811';'MEG1821';'MEG1831';'MEG1841';'MEG1911';'MEG1921';'MEG1931';'MEG1941';'MEG2011';'MEG2021';'MEG2031';'MEG2041';'MEG2111';'MEG2121';'MEG2131';'MEG2141';'MEG2211';'MEG2221';'MEG2231';'MEG2241';'MEG2311';'MEG2321';'MEG2331';'MEG2341';'MEG2411';'MEG2421';'MEG2431';'MEG2441';'MEG2511';'MEG2521';'MEG2531';'MEG2541';'MEG2611';'MEG2621';'MEG2631';'MEG2641';'MEG0113';'MEG0112';'MEG0122';'MEG0123';'MEG0132';'MEG0133';'MEG0143';'MEG0142';'MEG0213';'MEG0212';'MEG0222';'MEG0223';'MEG0232';'MEG0233';'MEG0243';'MEG0242';'MEG0313';'MEG0312';'MEG0322';'MEG0323';'MEG0333';'MEG0332';'MEG0343';'MEG0342';'MEG0413';'MEG0412';'MEG0422';'MEG0423';'MEG0432';'MEG0433';'MEG0443';'MEG0442';'MEG0513';'MEG0512';'MEG0523';'MEG0522';'MEG0532';'MEG0533';'MEG0542';'MEG0543';'MEG0613';'MEG0612';'MEG0622';'MEG0623';'MEG0633';'MEG0632';'MEG0642';'MEG0643';'MEG0713';'MEG0712';'MEG0723';'MEG0722';'MEG0733';'MEG0732';'MEG0743';'MEG0742';'MEG0813';'MEG0812';'MEG0822';'MEG0823';'MEG0913';'MEG0912';'MEG0923';'MEG0922';'MEG0932';'MEG0933';'MEG0942';'MEG0943';'MEG1013';'MEG1012';'MEG1023';'MEG1022';'MEG1032';'MEG1033';'MEG1043';'MEG1042';'MEG1112';'MEG1113';'MEG1123';'MEG1122';'MEG1133';'MEG1132';'MEG1142';'MEG1143';'MEG1213';'MEG1212';'MEG1223';'MEG1222';'MEG1232';'MEG1233';'MEG1243';'MEG1242';'MEG1312';'MEG1313';'MEG1323';'MEG1322';'MEG1333';'MEG1332';'MEG1342';'MEG1343';'MEG1412';'MEG1413';'MEG1423';'MEG1422';'MEG1433';'MEG1432';'MEG1442';'MEG1443';'MEG1512';'MEG1513';'MEG1522';'MEG1523';'MEG1533';'MEG1532';'MEG1543';'MEG1542';'MEG1613';'MEG1612';'MEG1622';'MEG1623';'MEG1632';'MEG1633';'MEG1643';'MEG1642';'MEG1713';'MEG1712';'MEG1722';'MEG1723';'MEG1732';'MEG1733';'MEG1743';'MEG1742';'MEG1813';'MEG1812';'MEG1822';'MEG1823';'MEG1832';'MEG1833';'MEG1843';'MEG1842';'MEG1912';'MEG1913';'MEG1923';'MEG1922';'MEG1932';'MEG1933';'MEG1943';'MEG1942';'MEG2013';'MEG2012';'MEG2023';'MEG2022';'MEG2032';'MEG2033';'MEG2042';'MEG2043';'MEG2113';'MEG2112';'MEG2122';'MEG2123';'MEG2133';'MEG2132';'MEG2143';'MEG2142';'MEG2212';'MEG2213';'MEG2223';'MEG2222';'MEG2233';'MEG2232';'MEG2242';'MEG2243';'MEG2312';'MEG2313';'MEG2323';'MEG2322';'MEG2332';'MEG2333';'MEG2343';'MEG2342';'MEG2412';'MEG2413';'MEG2423';'MEG2422';'MEG2433';'MEG2432';'MEG2442';'MEG2443';'MEG2512';'MEG2513';'MEG2522';'MEG2523';'MEG2533';'MEG2532';'MEG2543';'MEG2542';'MEG2612';'MEG2613';'MEG2623';'MEG2622';'MEG2633';'MEG2632';'MEG2642';'MEG2643';'EEG001';'EEG002';'EEG003';'EEG004';'EEG005';'EEG006';'EEG007';'EEG008';'EEG009';'EEG010';'EEG011';'EEG012';'EEG013';'EEG014';'EEG015';'EEG016';'EEG017';'EEG018';'EEG019';'EEG020';'EEG021';'EEG022';'EEG023';'EEG024';'EEG025';'EEG026';'EEG027';'EEG028';'EEG029';'EEG030';'EEG031';'EEG032';'EEG033';'EEG034';'EEG035';'EEG036';'EEG037';'EEG038';'EEG039';'EEG040';'EEG041';'EEG042';'EEG043';'EEG044';'EEG045';'EEG046';'EEG047';'EEG048';'EEG049';'EEG050';'EEG051';'EEG052';'EEG053';'EEG054';'EEG055';'EEG056';'EEG057';'EEG058';'EEG059';'EEG060';'EEG065';'EEG066';'EEG067';'EEG068';'EEG069';'EEG070';'EEG071';'EEG072';'EEG073';'EEG074';'HEOG';'VEOG';'STI101'};
                    types={'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGMAG';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'MEGPLANAR';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EEG';'EOG';'EOG';'Other'};
                    
                    
                    for ch=1:size(D.data,1)
                        D.channels(ch).label=chans{ch};
                        D.channels(ch).type=types{ch};
                        %end; save([ses_dir '/eMaffffd' basefile '.mat'],'D','-v7.3')
                    end; save([ses_out_dir '/' basefile '.mat'],'D','-v7.3')
                    
                    %clear D; D = spm_eeg_load([ses_dir '/eMaffffd' basefile '.mat']); %EK
                    
                    clear D; D = spm_eeg_load([ses_out_dir '/' basefile '.mat']); %AP
                    
                end
                
                D=D.badchannels(1:D.nchannels,0);
                
                
                D2=osl_detect_artefacts(D, 'modalities', {'MEGMAG', 'MEGPLANAR'}); D=D2;
                
                %             save([ses_dir '/aeMaffffd' basefile '.mat'],'D','-v7.3')
                %             D = spm_eeg_load([ses_dir '/aeMaffffd' basefile '.mat']); %EK
                
                save([ses_out_dir '/' 'a' basefile '.mat'],'D','-v7.3')
                
                D = spm_eeg_load([ses_out_dir '/' 'a' basefile '.mat']); %AP
                
                S = []; S.D = D;
                
                D = spm_eeg_remove_bad_trials(S);
                
            end
            
        end
        
    end
    
end

%%=========================================================================
%% Robust averaging & repeat lowpass
%%=========================================================================

if doaverage %EK
    %     parfor sub=1:length(subjects)
    parfor sub=1:length(subjects)
        
        for sess=1:length(sessions)
            
            %         ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
            %         D = spm_eeg_load([ses_dir '/raeMaffffd' basefile '.mat']); AP indent out
            
            
            %basefile = ['eMaffffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
            
            ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
            
            
            if isfile(fullfile([ses_out_dir,'/', 'raeMaffffdS' '_' subjects{sub} '_' sessions{sess} '_' 'rov123_sss.mat']))
                
                D = spm_eeg_load([ses_out_dir,'/', 'raeMaffffdS' '_' subjects{sub} '_' sessions{sess} '_' 'rov123_sss.mat']);
                
                S   = []; S.D = D;
                S.robust.savew = 0;
                S.robust.bycondition = 1;
                S.robust.ks = 3;
                S.robust.removebad = 0; %won't like this, as bad trials
                %aren't marked
                D = spm_eeg_average(S);
                
            end
            
        end
        
    end
    
    low=filtwindow(2);
    %     parfor sub=1:length(subjects)
    
    %  for sub=1:length(subjects); AP indent
    
    for sub=1:length(subjects)
        
        
        for sess=1:length(sessions)
            
            %         ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
            %         D = spm_eeg_load([ses_dir '/raeMaffffd' basefile '.mat']); AP indent out
            
            
            %basefile = ['eMaffffdS_' subjects{sub} '_' sessions{sess} '_rov123_sss'];
            
            ses_out_dir = [preprocdir '/' subjects{sub} '/' sessions{sess}];
            
            
            if isfile(fullfile([ses_out_dir,'/', 'mraeMaffffdS' '_' subjects{sub} '_' sessions{sess} '_' 'rov123_sss.mat']))
                
                %ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task); AP indent
                %out
                
                S = [];
                %S.D = [ses_dir '/mraeMaffffd' basefile '.mat']; AP indent out
                
                S.D = [ses_out_dir,'/', 'mraeMaffffdS' '_' subjects{sub} '_' sessions{sess} '_' 'rov123_sss.mat'];
                S.type = 'butterworth';
                S.order = 5;
                S.band = 'low';
                S.freq = low;
                D = spm_eeg_filter(S);
                
            end
            
        end
    end
    
end

%AP -this seems to be the end of it
%Its kept in sensor space, will I need to integrate it w/ co-registration and LFP extraction..? 

%%=========================================================================
%% Contrast
%%=========================================================================

if docontrast 
       
%     parfor sub=1:length(subjects)
   parfor sub=1:length(subjects)
        
        ses_dir = [ana_dir subjects{sub} '/' ses_nam '/' task '/']
        D = spm_eeg_load([ses_dir '/fmraeMaffffd' basefile '.mat']);
        
        S = [];
        S.D = D;
        S.label = {'DEV','REP3','REP6','DEV-REP3','DEV-REP6'};
        S.c = zeros(5,3);
        i=find(strcmp(D.conditions,'DEV')) %EK
        S.c(1,i) = 1;
        j=find(strcmp(D.conditions,'REP3'))
        S.c(2,j) = 1;
        k=find(strcmp(D.conditions,'REP6'))
        S.c(3,k) = 1;
        S.c(4,i) = 1; S.c(4,j) = -1;
        S.c(5,i) = 1; S.c(5,k) = -1;
        
        S.weighted = 0;
        D = spm_eeg_contrast(S);
        
    end
end


%%=========================================================================
%% Analysis specific low pass filter
%%=========================================================================

if dofilter2
    low=filtwindow2(2);
%     parfor sub=1:length(subjects)
    parfor sub=1:length(subjects)
        ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task)
        S = [];
        S.D = [ses_dir '/wfmraeMaffffd' basefile '.mat'];
        S.type = 'butterworth';
        S.order = 5;
        S.band = 'low';
        S.freq = low;
        D = spm_eeg_filter(S);
    end
end

%%=========================================================================
%% Rename files
%%=========================================================================
 if dorename
     
    parfor sub=1:length(subjects)
        ses_dir = fullfile(ana_dir,subjects{sub},ses_nam,task)
        S = [];
        S.D = [ses_dir '/fwfmraeMaffffd' basefile '.mat'];
        S.outfile = [ses_dir '/' subjects{sub} '_' ses_nam '_' task '.mat'];
        D = spm_eeg_copy(S);
    end
     
 end

end