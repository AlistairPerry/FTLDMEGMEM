function VBMsetup_genbrainmask(mrilist)

addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/;


fid = fopen(mrilist);
Data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
mrilist = Data{1,1};

nrun = length(mrilist);

jobfile = {'/imaging/rowe/Michelle/AlistairIFG/A_Scripts/spmbrainext_batch_job.m'};

split_stem = regexp(mrilist, '/', 'split');


jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);


for crun = 1:nrun

inputs{1, crun} = [cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/c1' split_stem{crun}{end}]); cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/c2' split_stem{crun}{end}]); cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/' split_stem{crun}{end}])];

inputs{2, crun} = ['bext_' '' split_stem{crun}{end}(1:end-4)];

inputs{3, crun} = cellstr(['/' fullfile(split_stem{crun}{1:end-1}) '/' 'avgbrain' '/']);

end


jobs = repmat(jobfile, 1, 1);

parfor crun = 1:nrun
    spm('defaults', 'PET');
    spm_jobman('initcfg')
    spm_jobman('run', jobs, inputs{:,crun}); 
end


end