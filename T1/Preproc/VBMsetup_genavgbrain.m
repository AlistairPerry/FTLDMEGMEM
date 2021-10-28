function VBMsetup_genavgbrain(mrilist)

addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/;


fid = fopen(mrilist);
Data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
mrilist = Data{1,1};

nrun = length(mrilist);
jobfile = {'/imaging/rowe/Michelle/FTLD_7T_CPSP/scripts/VBM_batch_normalise_unmodulated_unsmoothed.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

%core_imagepaths = [group1_mrilist; group2_mrilist(1:length(group1_mrilist))];
split_stem = regexp(mrilist, '/', 'split');
%split_stem_template = regexp(core_imagepaths, '/', 'split');
path_to_template_6 = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/Template_6.nii']);
%path_to_template_6 = cellstr(['/' fullfile(split_stem_template{1}{1:end-1}) '/Template_6.nii']);

for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}(1:end-4) '_Template.nii']); %Uses transformation to template space
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
    %Put in new dir so files hopefully get processed there
    inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/' 'avgbrain' '/' 'bext_' '' split_stem{crun}{end}]); %Native file
end


normaliseforaverageworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        normaliseforaverageworkedcorrectly(crun) = 1;
    catch
        normaliseforaverageworkedcorrectly(crun) = 0;
    end
end

nrun = 1;
jobfile = {'/imaging/rowe/Michelle/FTLD_7T_CPSP/scripts/VBM_batch_imcalc_average.m'};
inputs = cell(2, nrun);

% split_stem = regexp(core_imagepaths, '/', 'split');
% inputs{1,1} = cell(length(core_imagepaths),1);

split_stem = regexp(mrilist, '/', 'split');
inputs{1,1} = cell(length(mrilist),1);


% for i = 1:length(core_imagepaths)
%     inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '//w' split_stem{i}{end}]);
% end

%This should pick up all norm files in avg brain dir
for i = 1:length(mrilist)
    inputs{1,1}(i) = cellstr(['/' fullfile(split_stem{i}{1:end-1}) '/avgbrain/' 'wbext_' '' split_stem{i}{end}]);
end

inputs{2,1} = ['average_matched_control_patient_T1head'];

imcalcworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun %Still submit as a parfor to avoid overloading a login node
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        imcalcworkedcorrectly(crun) = 1;
    catch
        imcalcworkedcorrectly(crun) = 0;
    end
end


end