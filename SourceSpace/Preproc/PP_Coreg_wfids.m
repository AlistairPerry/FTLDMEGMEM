function D = PP_Coreg_wfids(D, subj, struct)

%For now let it run not including subjects w/out T1.. (should change)

od = cd(fileparts(D));
spm_jobman('initcfg'); 

% specify:    
m{1}.spm.meeg.source.headmodel.D = {D};
m{1}.spm.meeg.source.headmodel.val = 1;
m{1}.spm.meeg.source.headmodel.comment = '';

if strcmp('',struct)
    m{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1; 
    warning('!!! --- USING TEMPLATE STRUCTURAL --- !!!')
else
    m{1}.spm.meeg.source.headmodel.meshing.meshes.mri = {struct};
end

m{1}.spm.meeg.source.headmodel.meshing.meshres = 2;


% Pick up ficudials

%MEM

FIDSDIR = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/preproc/FIDS';

%TGB

%FIDSDIR = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/preproc/FIDS_TGB';


%New check to determine indivs w/out fids (but has T1)
%Nasion file

nasfile = [FIDSDIR '/' subj '_' 'nas.txt'];


if strcmp('',struct) || ~exist(nasfile)
% warning('!!! --- NOFIDS then--- !!!')

m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'Nasion';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.select = 'nas';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'LPA';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.select = 'lpa';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.select = 'rpa';

else 

%Now read them    

nasdata = dlmread(nasfile);

m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).fidname = 'Nasion';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(1).specification.type = [nasdata];


%LPA

lpafile = [FIDSDIR '/' subj '_' 'lpa.txt'];

lpadata = dlmread(lpafile);

m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).fidname = 'LPA';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(2).specification.type = [lpadata];


%RPA

rpafile = [FIDSDIR '/' subj '_' 'rpa.txt'];

rpadata = dlmread(rpafile);

m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).fidname = 'RPA';
m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.fiducial(3).specification.type = [rpadata];


end


m{1}.spm.meeg.source.headmodel.coregistration.coregspecify.useheadshape = 1;


m{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
m{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';

spm('defaults','EEG');
spm_jobman('run',m);

% [f1,f2,f3] = fileparts(D);
% load(D)
% D.fname = ['C' f2 f3];
% D.data.fname = [f1 filesep 'C' f2 '.dat'];
% save([f1 filesep 'C' f2 f3],'D')
% delete([f1 filesep f2 f3])
% movefile([f1 filesep f2 '.dat'],[f1 filesep 'C' f2 '.dat'])

cd(od)


end