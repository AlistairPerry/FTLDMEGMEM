function D = PP_Extract(D,RouteDir,varargin)


Lpos = [-46    20     8;  % These from Garrido 2008 / Phillips 2015
        -61   -32     8;
        -42   -22     7;
        46    20     8;
        59   -25     8;
        46   -14     8;
        ]; 
    
    Sname = {'LIFG';
        'LSTG';
        'LAud';
        'RIFG';
        'RSTG';
        'RAud';
        };
    
[f1,f2,f3] = fileparts(D);
if isempty(varargin), varargin{1} = 1; end
fn = [RouteDir 'LFP' num2str(length(Sname)) 'inv'  num2str(varargin{1}) '_' f2 f3];

% specify:
m{1}.spm.meeg.source.extract.D(1) = {D};
m{1}.spm.meeg.source.extract.val = varargin{1};
m{1}.spm.meeg.source.extract.source(1).label = Sname{1};
m{1}.spm.meeg.source.extract.source(1).xyz = Lpos(1,:);
m{1}.spm.meeg.source.extract.source(2).label = Sname{2};
m{1}.spm.meeg.source.extract.source(2).xyz = Lpos(2,:);
m{1}.spm.meeg.source.extract.source(3).label = Sname{3};
m{1}.spm.meeg.source.extract.source(3).xyz = Lpos(3,:);
m{1}.spm.meeg.source.extract.source(4).label = Sname{4};
m{1}.spm.meeg.source.extract.source(4).xyz = Lpos(4,:);
m{1}.spm.meeg.source.extract.source(5).label = Sname{5};
m{1}.spm.meeg.source.extract.source(5).xyz = Lpos(5,:);
m{1}.spm.meeg.source.extract.source(6).label = Sname{6};
m{1}.spm.meeg.source.extract.source(6).xyz = Lpos(6,:);
m{1}.spm.meeg.source.extract.rad = 7;
m{1}.spm.meeg.source.extract.type = 'trials';
m{1}.spm.meeg.source.extract.fname = fn;

spm('defaults','EEG');
spm_jobman('run',m);

clear D
load(fn);
x = num2cell(Lpos(:,1));
y = num2cell(Lpos(:,2));
[D.channels.X_plot2D] = deal(x{:});
[D.channels.Y_plot2D] = deal(y{:});
save(fn,'D')


end

