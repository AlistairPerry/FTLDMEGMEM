function maxfilter_FTDMMN_AP
%% Maxfilter
% Rik Henson, 2017; Modifications by AP, following LH and Ece K 2017

% Script detects and interpolates bad channels, performs signal space
% separation and shifts data to a default space


addpath /imaging/local/meg_misc/
addpath /neuro/meg_pd_1.2/

addpath /imaging/ap09/Toolbox/ekMEG/ntadscripts

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/'); %EK
addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/fieldtrip/'); %EK

% Directories

raw_dir = '/megdata/cbu/ftdrug/'; % Where MEG data lives (be careful you don't want to overwrite anything)
preproc_dir = '/imaging/ap09/Projects/FTD-MEG-MEM/newmaxfilter/preproc/Raw/B_Data/'; % Where you want to put the pre-processed MEG data


basefname = 'rov123';


%% Subjects

addpath('/imaging/ap09/Toolbox/extDCM/extDCMpreproc_scripts');

%Con_Sub = Mem_Control_Data();

%NA's in 17:20
%Con_Sub(17:20) = [];
% 
% %Broke at subject 12 (26.10.20) for some some reason
%Con_Sub(1:11) = [];


%Have to create these again for subjects that didn't run initially (cap
%letters in fname..)

Con_Sub = Mem_Control_Data();

%Pat_Sub = Mem_Patient_Data_redo();

%Pat_Sub(20) = [];

Pat_Sub = Mem_Patient_Data();

%Do I just pull out placebo sessions first?

%All_Subs = Con_Sub;

All_Subs = cat(2, Con_Sub, Pat_Sub);

nsubjs = length(All_Subs);


%I say, yes

for subj = 1:nsubjs
    
    subjects{subj,1} = All_Subs{1, subj}.Name{1, 3};
    
    
    if All_Subs{1, subj}.DrugSess == 1
        
        sessions{subj,1} = 's1';
        
        %Note some will not have scans
        scans{subj,1} = All_Subs{1, subj}.Name{1, 1};
        
    else
        
        sessions{subj,1} = 's2';
        
        %Note some will not have scans
        scans{subj,1} = All_Subs{1, subj}.Name{1, 2};
        
    end
    
    
end




%% Other setup

%which MF processes to run
run_autobad = 1;
run_tsss = 1;
run_trans_default = 1;
OverWrite = 1;


%settings
trans_offset_maxf = [0 -13 +6];

clear cmdb; clear cmd_badgo; clear cmd; clear badchstr; clear badch
clear filename; clear headpoints; clear fit



%% Locate files and create subject directories



for subj = 1:length(subjects)
    
    sess = sessions{subj};
    
    
    if strmatch(sessions{subj}, 's1')
        
      sessint = 1;
        
    else
        
      sessint = 2;
        
    end
    
    subj_lc = lower(subjects{subj});
    
    
    %make a dir to put the new MF sss files in
    
    mkdir([preproc_dir '/' subj_lc '/'], sess);
    
    
    % *** locate raw data in cbu dir **
    
    rawfiles{subj}.rp1 = strcat(raw_dir, All_Subs{subj}.Name{sessint}, '/'); %Parent folder
    rawfiles{subj}.rp2 = dir(rawfiles{subj}.rp1); %Scan folder
    
    rawfiles{subj}.rawpath = strcat(rawfiles{subj}.rp1, rawfiles{subj}.rp2(3).name, '/'); %Complete scan folder
    
    rawfiles{subj}.allraw = spm_select('List', rawfiles{subj}.rawpath, '.*raw.*fif'); %% %% All raw files
    
%     
%     %Fix strings now, as spm_select is returning fnames with '   '
%     
%     for nnn = 1:size(rawfiles{subj}.allraw)
%         
%         A = rawfiles{subj}.allraw(nnn,:);
%         A = A(~isspace(A));
%         
%         rawfiles{subj}.allraw(nnn,:) = A;
%         
%     end
       
    rawfiles{subj}.allraw1 = lower(rawfiles{subj}.allraw); %change to lower case
    
    rawfiles{subj}.allraw_rov = [];
    
    
    for nnn = 1:size(rawfiles{subj}.allraw1, 1)
        if strfind(rawfiles{subj}.allraw1(nnn,:), 'rov123')
            rawfiles{subj}.allraw_rov = [rawfiles{subj}.allraw_rov; rawfiles{subj}.allraw(nnn,:)];
        end
    end
    
    
    %Remove empty strings just in case
    
    if strmatch(subj_lc, 'p3')
        
    rawfiles{1, subj}.allraw_rov = rawfiles{1, subj}.allraw_rov(2,:);   
        
    end
    
    
    A = rawfiles{subj}.allraw_rov;
    
    A = A(~isspace(A));
    
    rawfiles{subj}.allraw_rov = A;
    
    rawfiles{subj}.allraw = rawfiles{subj}.allraw_rov;
    
    
    rawfiles{subj}.Numoffiles = size(rawfiles{subj}.allraw,1); %number of files, tho for here should be 1
    
    
end


%%=========================================================================
%% Begin
%%=========================================================================

basestr = ' -ctc /neuro/databases/ctc/ct_sparse.fif -cal /neuro/databases/sss/sss_cal.dat';
basestr = [basestr ' -linefreq 50 -hpisubt amp'];
basestr = [basestr ' -force'];
maxfstr = '!/neuro/bin/util/x86_64-pc-linux-gnu/maxfilter-2.2.12 ';
tStsssr = ' -st 10 -corr 0.98'; %'tStsssr = '';
dsstr = ''; % dsstr = ' -ds 4';  % downsample to 250Hz - can cause MF crash?
incEEG = 0;
TransTarget = '';


for subj = 1:length(subjects)
    %parfor sub = 1:length(subjects)
    
    sess = sessions{subj};
    
    
    if ~isempty(rawfiles{subj}.allraw)
        
        
        subj_lc = lower(subjects{subj});
        
        sub_dir = fullfile(preproc_dir,'/',subj_lc,'/',sess);
        
        movfile = fullfile(sub_dir, [subj_lc '_' sess '_' basefname '_' 'trans_move.txt']); % This file will record translations between runs
        
        rik_eval(sprintf('!touch %s',movfile));
        rik_eval(sprintf('!echo ''%s'' >> %s',datestr(now),movfile));
        
        
        % Fit sphere (since better than MaxFilter does)
        [co ki] = hpipoints([rawfiles{subj}.rawpath '' rawfiles{subj}.allraw]);
        
        cd(sub_dir)
        
        tmppoints = co(:,ki>1)';  % don't include the fiducial points
        
        nloc = find(~(tmppoints(:,3) < 0 & tmppoints(:,2) > 0));  % Remove the nose points if you have any
        
        headpoints = tmppoints(nloc,:);
        
        save([sub_dir '/hpi.txt'],'-ASCII','headpoints') % save headpoints
        cmd_fit = ['/imaging/local/software/neuromag/bin/util/fit_sphere_to_points hpi.txt'];
        
        
        if exist([sub_dir '/fittmp.txt'],'file') %EK
            s=dir('fittmp.txt');
            while s.bytes==0
                try delete('fittmp.txt'), end
                eval(['! ' cmd_fit ' > fittmp.txt']);
                s=dir('fittmp.txt');
            end
        else
            s=[];
            s.bytes=0;
            while s.bytes==0
                try delete('fittmp.txt'), end
                eval(['! ' cmd_fit ' > fittmp.txt']);
                s=dir('fittmp.txt');
            end
        end
        
        spherefit = dlmread('fittmp.txt');
        %orig = str2num(spherefit)*1000;  % m to mm;
        
        if length(spherefit)<5;
            error('Spherefit failed for a different reason!')
        end
        
        orig=spherefit*1000;
        origstr = sprintf(' -origin %d %d %d -frame head',orig(1),orig(2),orig(3))
        
        
        %%=================================================================
        %% 1. Autobad
        %%=================================================================
        % (this email says important if doing tSSS later
        % https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=NEUROMEG;d3f363f3.1205)
        
        if run_autobad
            outfile = fullfile(sub_dir, '/' , [subj_lc '_' sess '_' basefname '_bad']); 
            badfile = fullfile(sub_dir, '/' , [subj_lc '_' sess '_' basefname '_bad.txt']); 
            logfile = fullfile(sub_dir, '/' , [subj_lc '_' sess '_' basefname '_bad.log']);
            %badstr  = sprintf(' -autobad %d -badlimit %d',1800,7); % 1800s is 30mins - ie enough for all do_sessions?
            badstr  = sprintf(' -autobad %d -badlimit %d',20,7); % 1800s is 30mins - ie enough for all do_sessions?
            
            
            if ~exist(logfile,'file') | OverWrite
                filestr = sprintf(' -f %s -o %s.fif', [rawfiles{subj}.rawpath '' rawfiles{subj}.allraw], outfile);
                posfile = fullfile(sub_dir, '/', 'headpos.txt');
                compstr = sprintf(' -headpos -hpistep 10 -hp %s',posfile);
                finstr = [maxfstr filestr origstr basestr badstr compstr sprintf(' -v | tee %s.log',outfile)]
                rik_eval(finstr);
                %delete(sprintf('%s.fif',outfile));
            end
            
            % Pull out bad channels from logfile:
            delete(badfile); %ask Ek about this
            rik_eval(sprintf('!echo '' '' > %s', badfile));
            rik_eval(sprintf('!cat %s.log | sed -n -e ''/Detected/p'' -e ''/Static/p'' | cut -f 5- -d '' '' >> %s',outfile,badfile));
            tmp=dlmread(badfile,' '); % Nbuf = size(tmp,1)-1; Doesn't work if subset of channels bad for short time
            tmp=reshape(tmp,1,prod(size(tmp)));
            tmp=tmp(tmp>0); % Omit zeros (padded by dlmread):
            
            % Mark bad based on threshold (currently ~5% of buffers (assuming 500 buffers)):
            [frq,allbad] = hist(tmp,unique(tmp));
            badchans = allbad(frq>0.05*max(frq)); %save badchans.mat allbad badchans; %EK assuming that there is always at least one sensor that's bad throughout eg 813
            if isempty(badchans) badstr = '';
                
            else
                badstr = sprintf(' -bad %s',num2str(badchans))
            end
            
        end
        
        
        %%=================================================================
        %% 2. tSSS
        %%=================================================================
        % (ie, align within subject if multiple sessions)
        
        if run_tsss
            if ~exist(sprintf('%s.fif',outfile),'file') | ~exist([outfile '.log'],'file') | OverWrite
                % waitbar(0.66, h, ['Running tSSS and movecomp on Sub:' subjects{sub}(7:10) ' Task:' do_sessions{ses} '...']); %EK
                transfstfile = [outfile '.fif'];
                %if ~isempty(TransTarget) &
                %~strcmp(do_tasks{ses},TransTarget) - do_tasks is
                %problematic here
                
                %Note I am changing from tsss to sss for compatibility
                %with later scripts
                
                if ~isempty(TransTarget)
                    outfile = fullfile(sub_dir,sprintf('tsss_trans_%s',TransTarget))
                    transtr = sprintf(' -trans %s/tsss.fif',fullfile(sub_dir,TransTarget));
                else
                    transtr = '';
                    outfile = fullfile(sub_dir, '/' , [subj_lc '_' sess '_' basefname '_' sprintf('sss')]);
                end
                
                compstr = sprintf(' -movecomp inter'); %compstr = sprintf(' -movecomp');
                filestr = sprintf(' -f %s -o %s.fif', [rawfiles{subj}.rawpath '' rawfiles{subj}.allraw], outfile);
                finstr = [maxfstr filestr basestr badstr tStsssr compstr origstr transtr dsstr sprintf(' -v | tee %s.log',outfile)]
                rik_eval(finstr);
            
            end
            
            
            %if ~isempty(TransTarget) & ~strcmp(do_tasks{ses},TransTarget)
            %- again is problematic - ask Ece about this
                %fprintf(fp,'\nTransfirst %s to %s: ',do_sessions{ses},TransTarget  
                %rik_eval(sprintf('!echo ''Trans %s to %s'' >> %s',do_tasks{ses},TransTarget,movfile));
                if ~isempty(TransTarget)
                rik_eval(sprintf('!cat %s.log | sed -n ''/Position change/p'' | cut -f 7- -d '' '' >> %s',outfile,movfile));
            end
        end
        
        %%=================================================================
        %% 3. Trans default
        %%=================================================================
        %Note changing to ssst for compatibility with later scripts
        
        if run_trans_default
            transdeffile = fullfile(sub_dir, [subj_lc '_' sess '_' basefname '_' sprintf('ssst')]);
            
            if ~exist([transdeffile '.fif'],'file') | ~exist([transdeffile '.log'],'file') | OverWrite
                outfile = fullfile(sub_dir, '/' , [subj_lc '_' sess '_' basefname '_' sprintf('sss')]);
                %waitbar(0.99, h, ['Running trans default on Sub:' subjects{sub}(7:10) ' Task:' do_sessions{ses} '...']); %EK
                transtr = [' -trans default -origin ' num2str(orig(1)+0) ' ' num2str(orig(2)-13) ' ' num2str(orig(3)+6) ' -frame head -force'];
                filestr = sprintf(' -f %s.fif -o %s.fif',outfile,transdeffile);
                finstr = [maxfstr filestr transtr sprintf(' -v | tee %s.log',transdeffile)]
                rik_eval(finstr);
            end
            
            %rik_eval(sprintf('!echo ''Transdef for %s'' >> %s',do_tasks{ses},movfile));
            rik_eval(sprintf('!echo ''Transdef for %s'' >> %s',subj_lc,movfile));
            
            %fprintf(fp,'\nTransDef %s: ',do_sessions{ses});
            rik_eval(sprintf('!cat %s.log | sed -n ''/Position change/p'' | cut -f 7- -d '' '' >> %s',transdeffile,movfile));
            % delete(h); %EK
            %fclose(fp);
        end
        
    end
    
end

end