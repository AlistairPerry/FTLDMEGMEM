function normalizeGM_forIFGvol(mrilist)


%% Now normalise all scans

%Orig mri list - for pulling out warp files

fid = fopen(mrilist);
mrilist = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
mrilist = mrilist{1,1};

nrun = length(mrilist);
jobfile = {'/imaging/rowe/Michelle/AlistairIFG/A_Scripts/VBM_batch_normalise_modulated_unsmoothed.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);

split_stem = regexp(mrilist, '/', 'split');
%split_stem_template = regexp(core_imagepaths, '/', 'split');
path_to_template_6 = cellstr(['/' fullfile(split_stem{1}{1:end-1}) '/Template_6.nii']);

for crun = 1:nrun
    inputs{1, crun} = path_to_template_6; % Normalise to MNI Space: Dartel Template - cfg_files
    inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}(1:end-4) '_Template.nii']);
    if ~exist(char(inputs{2,crun}),'file')
        inputs{2, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/u_rc1' split_stem{crun}{end}]); %Deal with the difference in naming convention depending if part of the dartel templating or not
    end
    %inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/c1' split_stem{crun}{end}]);
end


for crun = 1:nrun
    inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/' 'c1' split_stem{crun}{end}]);
end

normaliseworkedcorrectly = zeros(1,nrun);
jobs = repmat(jobfile, 1, 1);


save('/imaging/rowe/Michelle/AlistairIFG/A_Scripts/inputs_norm_GMIFGvol.mat','inputs');


parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    try
        spm_jobman('run', jobs, inputs{:,crun});
        normaliseworkedcorrectly(crun) = 1;
    catch
        normaliseworkedcorrectly(crun) = 0;
    end
end


end