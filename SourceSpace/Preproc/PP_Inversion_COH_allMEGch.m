function D = PP_Inversion_COH_allMEGch(D,varargin)

    
if isempty(varargin), varargin{1} = 1; end
if isempty(varargin) || length(varargin)<2, varargin{2} = [0 180]; end
if isempty(varargin) || length(varargin)<3, varargin{3} = {'Dev' 'rep1' 'rep2' 'rep3' 'rep4' 'rep5' 'rep6' 'rep7' 'rep8' 'rep9'}; end
%Meth = {'IID' 'MSP'};

m{1}.spm.meeg.source.invert.D = {D};
m{1}.spm.meeg.source.invert.val = varargin{1};
m{1}.spm.meeg.source.invert.whatconditions.condlabel = varargin{3};
m{1}.spm.meeg.source.invert.isstandard.custom.invtype = 'COH'; % 'IID'; % %Meth{varargin{1}};
m{1}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
m{1}.spm.meeg.source.invert.isstandard.custom.foi = varargin{2};
m{1}.spm.meeg.source.invert.isstandard.custom.hanning = 0;
m{1}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
m{1}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
m{1}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
m{1}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
m{1}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
m{1}.spm.meeg.source.invert.modality = {
                                        'MEG'
                                        'MEGPLANAR'
                                        }';

spm('defaults','EEG');
spm_jobman('run',m);


end

