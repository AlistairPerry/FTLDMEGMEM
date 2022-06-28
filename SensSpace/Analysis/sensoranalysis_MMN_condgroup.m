function sensoranalysis_MMN_groupdrug


%% NTAD Sensor level RMS analysis for Roving (Ece K, 2019)

%restoredefaultpath


%Add paths

addpath('/imaging/rowe/users/ap09/Toolbox/')

addpath('/imaging/rowe/users/ap09/Toolbox/ekMEG/ntadscripts')

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/boundedline-pkg'))

addpath(genpath('/imaging/rowe/users/ap09/Toolbox/VBA-toolbox'))


% spm eeg;

addpath /imaging/local/software/spm_cbu_svn/releases/spm12_latest/;

addpath('/imaging/local/software/spm_cbu_svn/releases/spm12_latest/external/fieldtrip/'); %EK

ft_defaults()


%% Other setup

subjs = {'c2','c3','c7','c8','c9','c11','c13','c14','c18','c20','c23','c24','c1','c4','c5','c6','c10','c12','c15','c16','c17','c19','c21','c22', 'p1','p3','p6','p7','p9','p12','p13','p14','p17','p18','p21','p24','p2','p4','p5','p8','p10','p11','p15','p16','p19','p20','p22','p23','p25'};

sessions = {'s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s1','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2','s2'};



%Final subjects - those w/both sessions



subjs_new = {'c2';'c3';'c7';'c8';'c9';'c11';'c13';'c14';'c23';'c24';'c1';'c4';'c5';'c6';'c10';'c12';'c15';'c16';'c21';'c22'; 'p1';'p6';'p7';'p9';'p12';'p13';'p14';'p17';'p18';'p21';'p24';'p2';'p4';'p5';'p8';'p11';'p16';'p22';'p23';'p25'};


%Now find their sessions

for i = 1:length(subjs_new)

match = strcmp(subjs, subjs_new{i});


if nnz(match)

sessions_new(i,1) = sessions(match);

end
    
end


sessions = sessions_new;

subjs = subjs_new;


preprocdir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/preproc/preproc/SensSpace/ERF/B_Data/Adv_fix_7';

preprocsteps = 'fmraeMaffffdS';

basefname = 'rov123_ssst.mat'; %for now


%OUTpre = 'adv_ssst_newmaxf_fixICA_normALLtrls_automeg_nonorm';

%Remove outlier

%OUTpre = 'adv_ssst_newmaxf_fixICA_normALLtrls_automeg_nonorm_remo';


%Both sessions no outliers

OUTpre = 'adv_ssst_newmaxf_fixICA_normALLtrls_automeg_nonorm_remo_bothsess_newsensors';


FigOutDir = '/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/newmaxfilter/analysis/SensSpace/ERF_icafixes/C_Plots';


xticks=[1 50 100 150 200 251];
xlabels={'-100' '0' '100' '200' '300' '400'};
conditions={'DEV', 'REP6'};


trialsint = [1;4;7];

devtrl = 1;
stdtrl = 3;
std3trl = 2;



%% Disease group info
% Can call it directly as it will never change

DiagData = readcell('/imaging/rowe/users/ap09/Projects/FTD-MEG-MEM/ERP/analysis/LFPs/A_Scripts/PatDiagInfo.txt');



%% Channel selection

%Old
%fr_mag={'MEG0111','MEG0121','MEG0341','MEG0131','MEG1411','MEG1421','MEG1221','MEG1441'};
%fr_grad={'MEG0112','MEG0113','MEG0123','MEG0122','MEG0342','MEG0343','MEG0133','MEG0132',...
    %'MEG1412','MEG1413','MEG1422','MEG1423','MEG1222','MEG1223','MEG1442','MEG1443'};

fr_eeg={'EEG002','EEG004','EEG005','EEG006','EEG007','EEG008','EEG009','EEG010','EEG011','EEG012','EEG013','EEG014','EEG015',...
    'EEG016','EEG017','EEG019','EEG020','EEG021','EEG022','EEG023','EEG024','EEG025','EEG026','EEG027'};

%tem_mag={'MEG0241','MEG0231','MEG0441','MEG1611','MEG1621','MEG1811','MEG1821','MEG1641','MEG1631','MEG1841',...
    %'MEG1131','MEG1341','MEG1331','MEG2211','MEG2221','MEG2411','MEG2431','MEG2441'};
%tem_grad={'MEG0242','MEG0243','MEG0233','MEG0232','MEG0442','MEG0443','MEG1612','MEG1613','MEG1623','MEG1622','MEG1812',...
    %'MEG1813','MEG1823','MEG1822','MEG1642','MEG1643','MEG1633','MEG1632','MEG1842','MEG1843',...
    %'MEG1132','MEG1133','MEG1342','MEG1343','MEG1332','MEG1333','MEG2212','MEG2213','MEG2222','MEG2223','MEG2412','MEG2413',...
    %'MEG2432','MEG2433','MEG2442','MEG2443'};
tem_eeg={'EEG030','EEG031','EEG032','EEG033','EEG034','EEG035','EEG036','EEG037','EEG041','EEG042','EEG043','EEG044','EEG045',...
    'EEG046','EEG047','EEG048'};


%New sensors
fr_mag={'MEG0111','MEG0121','MEG0341','MEG0321','MEG0131','MEG0211','MEG0221',...
    'MEG1231', 'MEG1221','MEG1411','MEG1421','MEG1311','MEG1321','MEG1441'};
fr_grad={'MEG0112','MEG0113','MEG0123','MEG0122','MEG0342','MEG0343','MEG0322','MEG0323','MEG0133','MEG0132','MEG0212','MEG0213','MEG0222','MEG0223',...
    'MEG1232','MEG1233','MEG1222','MEG1223','MEG1412','MEG1413','MEG1422','MEG1423','MEG1312','MEG1313','MEG1322','MEG1323','MEG1442','MEG1443'};
tem_mag={'MEG0241','MEG0231','MEG0441','MEG0431','MEG1521','MEG1611','MEG1621','MEG1811','MEG1821','MEG1721','MEG1641','MEG1631','MEG1841',...
    'MEG1141','MEG1131','MEG1341','MEG1331','MEG2211','MEG2221','MEG2411','MEG2421','MEG2641','MEG2231','MEG2441','MEG2431','MEG2521'};
tem_grad={'MEG0242','MEG0243','MEG0233','MEG0232','MEG0442','MEG0443','MEG0432','MEG0433','MEG1522','MEG1523','MEG1612','MEG1613','MEG1623','MEG1622','MEG1812','MEG1813',...
    'MEG1823','MEG1822','MEG1722','MEG1723','MEG1642','MEG1643','MEG1633','MEG1632','MEG1842','MEG1843','MEG1142','MEG1143','MEG1132','MEG1133','MEG1342','MEG1343','MEG1332','MEG1333','MEG2212','MEG2213','MEG2222','MEG2223','MEG2412','MEG2413','MEG2422','MEG2423',...
    'MEG2642','MEG2643','MEG2232','MEG2233','MEG2442','MEG2443','MEG2432','MEG2433','MEG2522','MEG2523'};



%% Gather data - altogether

c_ind=find(contains(subjs,'c'));

subcount = 0;

submatcount_pat = 0;


for sub=1:length(subjs)
    
    sess = sessions{sub};
    
    
    filen=[preprocdir '/' subjs{sub} '/' sess '/' preprocsteps '_' subjs{sub} '_' sess '_' basefname];
    
    
    if isfile(fullfile([preprocdir '/' subjs{sub} '/' sess '/' preprocsteps '_' subjs{sub} '_' sess '_' basefname]))
    
    subcount = subcount + 1;
    
    
    allsubjs{subcount,1} = subjs{sub};
    
    
    cmatch = contains(subjs{sub},'c');
    
    
    if cmatch == 1
    
        grpst(subcount,1) = 1;
        
        patgrpst(subcount,1) = 3;
        
    else
        
        grpst(subcount,1) = -1;
        
        
        %A dumb way to do it but oh well
        
        submatcount_pat = submatcount_pat + 1;
        
        
        diagmatch = strcmp(subjs{sub}, DiagData(:,1));
        
        patgrpst(subcount,1) = DiagData{diagmatch==1,2};
        
        
    end
        
    
    D=spm_eeg_load(filen);
    
%     origdata(subcount,1:size(D,1),:,:)=D(:,:,:);
    
    %     if size(D,1)>320
    %         coormeg(subcount,:,1:380)=D.coor2D(1:380);
    %     else
    %         coormeg(subcount,:,1:306)=D.coor2D;
    %     end
    
    data{1}(subcount,:,:,:)=D(D.selectchannels('MEGMAG'),:,trialsint);
    data{2}(subcount,:,:,:)=D(D.selectchannels(fr_mag),:,trialsint);
    data{3}(subcount,:,:,:)=D(D.selectchannels(tem_mag),:,trialsint);
    data{4}(subcount,:,:,:)=D(D.selectchannels('MEGPLANAR'),:,trialsint);
    data{5}(subcount,:,:,:)=D(D.selectchannels(fr_grad),:,trialsint);
    data{6}(subcount,:,:,:)=D(D.selectchannels(tem_grad),:,trialsint);
    
    if size(D,1)>350
        data{7}(subcount,:,:,:)=D([307:366 370:379],:,trialsint);
        data{8}(subcount,:,:,:)=D(D.selectchannels(fr_eeg),:,trialsint);
        data{9}(subcount,:,:,:)=D(D.selectchannels(tem_eeg),:,trialsint);
    else
        data{7}(subcount,:,:,:)=NaN;
        data{8}(subcount,:,:,:)=NaN;
        data{9}(subcount,:,:,:)=NaN;
    end
    
    end

end



for i=1:9
    for s=1:size(data{i},1)
        for tr=1:size(data{i},4)
            for ch=1:size(data{i},2)
                if data{i}(s,ch,85,tr)>0
                    data{i}(s,ch,:,tr)=-(data{i}(s,ch,:,tr)); % N100 is negative
                    %data{i}(s,ch,:,tr)=NaN(1,1,501);
                end
            end
        end
    end
end


for i=1:9
    %datarms{i}(:,:,:,:)=permute(data{i}(:,:,:,:),[2 1 4 3]);
    %datarms{i}=squeeze(rms_ek(datarms{i},1,[1 100]));
    d=squeeze(nanmean(data{i}(:,:,:,:),2));
    d=permute(d,[1 3 2]);
    datarms{i}=d;
    
    for j=1:size(datarms{i},1)
        for k=1:size(datarms{i},2)
            d=squeeze(datarms{i}(j,k,:));
            bl=nanmean(d(1:51)); d=d-bl; % baseline correct again just in case
            datarms{i}(j,k,:)=smooth(d,20,'moving');
        end
    end
    
     ind=find(isnan(datarms{i}));
%     for s=1:size(datarms{i},1)
%         datarms{i}(s,:,:)=mat2gray(datarms{i}(s,:,:));
%     end
    
    datarms{i}(ind)=NaN;
    datarms{i}=datarms{i}(:,:,1:251);
    
end


%% Plot Means
clear data
varnm={'MAG','Frontal MAG','Temporal MAG','GRAD','Frontal GRAD','Temporal GRAD','EEG','Frontal EEG','Temporal EEG'};
%

proc_ctrs = find(grpst==1);
for i=1:9
    
    data=datarms{i}(proc_ctrs,:,:);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    [h,p,ci,stats]=ttest(squeeze(data(:,devtrl,:)),squeeze(data(:,stdtrl,:)));
    k=max(max([m1(2,:); m1(1,:)]))+max(max([s1(2,:);s1(1,:)]));k=k*1.1;
    p=0.9.*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
 
    
%   Don't need to do the below 
%     temp_var = strcat( 'figure',num2str(source) );
%     
%     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(i);set(gcf, 'Units','normalized', 'Position', [0 0 .25 .4]);
    plot(m1(2,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'REP6','DEV'},'fontsize',14)
    boundedline([1:size(m1,2)],m1(2,:),s1(2,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('Normalised signal')
    title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Save figure
    saveas(figure(i), [FigOutDir '/' OUTpre '_' 'ERF_'  'Controls' '_' varnm{i} '_fullsample.tif']);
    
end


figcount = i;


%Repeat with patients

proc_pats = find(grpst==-1);
for i=1:9
    
    figcount = figcount + 1;
    
    data=datarms{i}(proc_pats,:,:);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    [h,p,ci,stats]=ttest(squeeze(data(:,devtrl,:)),squeeze(data(:,stdtrl,:)));
    k=max(max([m1(2,:); m1(1,:)]))+max(max([s1(2,:);s1(1,:)]));k=k*1.1;
    p=0.9.*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
 
    
%   Don't need to do the below 
%     temp_var = strcat( 'figure',num2str(source) );
%     
%     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 .25 .4]);
    plot(m1(2,:),'LineWidth',4,'Color',[255/255 0/255 0/255]); hold on
    plot(m1(1,:),'LineWidth',4,'Color',[0 0 255]./255); legendflex({'REP6','DEV'},'fontsize',14)
    boundedline([1:size(m1,2)],m1(2,:),s1(2,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0 0 255]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('Normalised signal')
    title([varnm{i} ' mean in patients']); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'Patients' '_' varnm{i} '_fullsample.tif']);
    
end


%% Difference calculations - standard 6


%First, set up design matrix for RFT

nfilesc = length(proc_ctrs);
nfilesp = length(proc_pats);

DM = zeros(nfilesc+nfilesp,2);

DM(1:nfilesc,1) = repmat(1, [nfilesc 1]);

DM(nfilesc+1:end,2) = repmat(1, [nfilesp 1]);


c = [1; -1];


for i=1:9

    figcount = figcount + 1;
    
    %Load in again (sigh)
    data=datarms{i}(:,:,:);
    
    
    %Std minus diff - note rep6 will now be third
    mmndiff = squeeze(data(:,stdtrl,:)-data(:,devtrl,:));
    
    
    %Separate for controls and patients
    
    mmndiffcon = mmndiff(proc_ctrs,:);
    mmndiffpat = mmndiff(proc_pats,:);
    
    
    %Means and stds of difference
    m1diff_con=squeeze(nanmean(mmndiffcon,1));
    s1diff_con=squeeze(nanstd(mmndiffcon)./sqrt(size(mmndiffcon,1)));
    
    
    m1diff_pat=squeeze(nanmean(mmndiffpat,1));
    s1diff_pat=squeeze(nanstd(mmndiffpat)./sqrt(size(mmndiffpat,1)));

    
    %Note, independent ttest now
    [h,p,ci,stats]=ttest2(mmndiffcon,mmndiffpat);
    k=max(max([m1diff_con; m1diff_pat]))+max(max([s1diff_con;s1diff_pat]));k=k*1.1;
    
    
    %FDR-corrected
    [~, ~, ~, adj_p]=fdr_bh(p(52:end));
    adj_p_wbl = cat(2, NaN(1, 51), adj_p);
    
    %Uncorrected
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
 
    
%   Don't need to do the below 
%     temp_var = strcat( 'figure',num2str(source) );
%     
%     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 0.33 0.25]);   
    plot(m1diff_con,'LineWidth',4,'Color',[0/255 0/255 255/255]); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[255 0 0]./255); %legendflex({'CON','FTLD'},'fontsize',14, 'FontWeight', 'bold')
    boundedline([1:size(m1diff_con,2)],m1diff_con(1,:),s1diff_con(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[255 0 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    %plot(p(adj_p_wbl<0.05),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',13, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial'); ylabel({'Mismatch response'; '(rep6-dev)'}, 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    title([varnm{i} ' MMN in CON and FTD']); xlim([0 251]); %ylim([0 1]);%axis square
    
    box off

    set(gca, 'TickDir', 'Out')
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff_' 'rep6' '_' 'ConsandPats_' varnm{i} '_fullsample.tif']);
    
   
    % RFT for gradients only
    
    if i > 3 && i < 7
        
        
        figcount = figcount + 1;
        
        
        mmndiffall = cat(1, mmndiffcon, mmndiffpat);
        
        y = mmndiffall;
        
        
        [stat,out] = RFT_GLM_contrast(DM,y,c,'F',1,1);
        
        
        %Save figure
        
        saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff' '_' 'RFT' '_' 'ConsandPats_' varnm{i} '_fullsample.tif']);

        
        %Saveoutput 
        
        save([FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff' '_' 'RFTstats' '_' 'ConsandPats_' varnm{i} '_fullsample.mat'], 'stat');
        
        
    end
    
    
end



% Difference calculations - standard rep 3

peakneg_all = [];


for i=4:6
%for i=4:6

    figcount = figcount + 1;
    
    %Load in again (sigh)
    data=datarms{i}(:,:,:);
    
    
    %Std minus diff - note rep3 will now be second
    mmndiff = squeeze(data(:,std3trl,:)-data(:,devtrl,:));
    
    
    %Separate for controls and patients
    
    mmndiffcon = mmndiff(proc_ctrs,:);
    mmndiffpat = mmndiff(proc_pats,:);
    
    
    %Means and stds of difference
    m1diff_con=squeeze(nanmean(mmndiffcon,1));
    s1diff_con=squeeze(nanstd(mmndiffcon)./sqrt(size(mmndiffcon,1)));
    
    
    m1diff_pat=squeeze(nanmean(mmndiffpat,1));
    s1diff_pat=squeeze(nanstd(mmndiffpat)./sqrt(size(mmndiffpat,1)));

    
    %Find their average peak negative deflection

    conpeakneg = [find(m1diff_con(51:251)==min(m1diff_con(51:251)))]*2;
    
    patpeakneg = [find(m1diff_pat(51:251)==min(m1diff_pat(51:251)))]*2;


    peakneg_all(i-3,1) = conpeakneg;

    peakneg_all(i-3,2) = patpeakneg;
    

    %Note, independent ttest now
    [h,p,ci,stats]=ttest2(mmndiffcon,mmndiffpat);
    k=max(max([m1diff_con; m1diff_pat]))+max(max([s1diff_con;s1diff_pat]));k=k*1.1;
    
    
    %FDR-corrected
    [~, ~, ~, adj_p]=fdr_bh(p(52:end));
    adj_p_wbl = cat(2, NaN(1, 51), adj_p);
    
    %Uncorrected
    p_all=p;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
 
    
%   Don't need to do the below 
%     temp_var = strcat( 'figure',num2str(source) );
%     
%     eval(sprintf('%s = %s',temp_var, ['figure(' num2str(source) ')']))
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 0.33 0.25]);    
    plot(m1diff_con,'LineWidth',4,'Color',[0/255 0/255 255/255]); hold on
    plot(m1diff_pat,'LineWidth',4,'Color',[255 0 0]./255); %legendflex({'CON','FTLD'},'fontsize',14, 'FontWeight', 'bold')
    boundedline([1:size(m1diff_con,2)],m1diff_con(1,:),s1diff_con(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat,2)],m1diff_pat(1,:),s1diff_pat(1,:),'cmap',[255 0 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    %plot(p(adj_p_wbl<0.05),'-o','LineWidth',2,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0],'MarkerSize',8)
    
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',13, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial'); ylabel({'Mismatch response'; '(rep3-dev)'}, 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    title([varnm{i} ' MMN3 in CON and FTD']); xlim([0 251]); %ylim([0 1]);%axis square
    
    box off

    set(gca, 'TickDir', 'Out')
    


    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff_' 'rep3' '_' 'ConsandPats_' varnm{i} '_fullsample.tif']);
    

    figcount = figcount + 1;
    
    
    mmndiffall = cat(1, mmndiffcon, mmndiffpat);
    
    y = mmndiffall;
    
    
    [stat,out] = RFT_GLM_contrast(DM,y,c,'F',1,1);
    
    
    %Save figure
    
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff_' 'rep3' '_' 'RFT' '_' 'ConsandPats_' varnm{i} '_fullsample.tif']);
    
    
    %Saveoutput
    
    save([FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff_' 'rep3' '_' 'RFTstats' '_' 'ConsandPats_' varnm{i} '_fullsample.mat'], 'stat');
    
    
end


%Save average latency
dlmwrite([FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff_' 'rep3' '_' 'ConsandPats_' 'avgnegpeak' '_fullsample.txt'], peakneg_all, 'delimiter', '\t');


%%=========================================================================
%% Add in mean amplitude stuff
%%=========================================================================
%Setup mean MMN

% Calculate means for MMN ANOVA

ds = 1;

if ds == 1
    
    %Determine MMN start/finish
    
    MMNs = 226./2;
    MMNf = 276./2;
    
    
    %And also p3
    
    p3s = floor(375./2);
    p3f = floor(425./2);
    
else
    
    MMNs = 226;
    MMNf = 276;
    
    p3s = 275./2;
    p3f = 325./2;
    
end


%Design matrix for ANOVA..

nfilesc = length(proc_ctrs);

nfilesp = length(proc_pats);


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



%% Run mean amplitude analysis
% First MMN for rep6

for i=1:9
    
%Load in data (guess what, again)

data=datarms{i}(:,:,:);
    

%Mean amplitude for std and dev

constd = data(proc_ctrs, stdtrl, :);

condev = data(proc_ctrs, devtrl, :);


patstd = data(proc_pats, stdtrl, :);

patdev = data(proc_pats, devtrl, :);


% m1_constdmmn = nanmean(constd(:,MMNs:MMNf),2);
% 
% m1_condevmmn = nanmean(condev(:,MMNs:MMNf),2);
% 
% 
% m1_patstdmmn = nanmean(patstd(:,MMNs:MMNf),2);
% 
% m1_patdevmmn = nanmean(patdev(:,MMNs:MMNf),2);


%And their interaction significance

% meanMMNcol_all = cat(1, m1_constdmmn, m1_condevmmn, m1_patstdmmn, m1_patdevmmn);
% 
% X = cat(2, meanMMNcol_all, F1all, F2all, Sall);
% 
% 
% [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
% 
% 
% pgrp = Ps(1);
% pcond = Ps(3);
% pint = Ps(4);
% 
% MMNpall(i,1)=pgrp;
% MMNpall(i,2)=pcond;
% MMNpall(i,3)=pint;


%Can push out into rainclouds for plotting in R :)
% MMNtable = table(meanMMNcol_all, F1all, F2all, Sall);
% 
% MMNtable.Properties.VariableNames = {'meanMMNcol', 'Group', 'Trial', 'ID'};
% 
% writetable(MMNtable, [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNmean' '_' 'ConsandPats_' varnm{i} '_fullsample.txt'], 'Delimiter','tab')
%   


%Add difference score for use in JASP

conmmndiff = nanmean(constd(:,MMNs:MMNf),2)-nanmean(condev(:,MMNs:MMNf),2);

patmmndiff = nanmean(patstd(:,MMNs:MMNf),2)-nanmean(patdev(:,MMNs:MMNf),2);



[~,p,~,stats] = ttest2(conmmndiff,patmmndiff,'tail','both');

MMNpall(i,1) = p;
MMNtall(i,1) = stats.tstat;



%Collate for output table

meanMMNdiffcol_all = cat(1, conmmndiff, patmmndiff);

meanMMNdiffcol_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1])); 


MMNdifftable = table(meanMMNdiffcol_all, meanMMNdiffcol_grp, patgrpst);

MMNdifftable.Properties.VariableNames = {'meanMMNcol', 'Group', 'PatSubgrp'};

writetable(MMNdifftable, [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiffmean' '_' 'ConsandPats_' varnm{i} '_fullsample.txt'], 'Delimiter','tab')


%Below not needed

% meanMMNcol_lw = cat(2, m1_constdmmn, m1_condevmmn);
% 
% meanMMNcol_lw_pat = cat(2, m1_patstdmmn, m1_patdevmmn);
% 
% meanMMNcol_lw_all = cat(1, meanMMNcol_lw, meanMMNcol_lw_pat);
% 
% meanMMNcol_lw_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));
% 
% MMNtable_lw = table(meanMMNcol_lw_all(:,1), meanMMNcol_lw_all(:,2), meanMMNcol_lw_grp);
% 
% 
% MMNtable_lw.Properties.VariableNames = {'meanMMNcol_STD', 'meanMMNcol_Dev', 'Group'};
% 
% 
% writetable(MMNtable_lw, [FigOutDir '/' OUTpre '_' 'ERF_'  'STDDEV' '_' 'ConsandPats_' varnm{i} '_fullsample_forJASP.txt'], 'Delimiter','tab')


end


%Save output table of mean MMN

MMNstats = struct;

MMNstats.t = MMNtall;

MMNstats.p = MMNpall;


save([FigOutDir '/' OUTpre '_' 'ERF_'  'MMNmean' '_restable_ConsandPats_ttest.mat'], 'MMNstats');




%% Now p300 
for i=1:9
    
%Load in data (guess what, again)

data=datarms{i}(:,:,:);
    

%Mean amplitude for std and dev

constd = data(proc_ctrs, stdtrl, :);

condev = data(proc_ctrs, devtrl, :);


patstd = data(proc_pats, stdtrl, :);

patdev = data(proc_pats, devtrl, :);


m1_constd_p3 = nanmean(constd(:,p3s:p3f),2);

m1_condev_p3 = nanmean(condev(:,p3s:p3f),2);


m1_patstd_p3 = nanmean(patstd(:,p3s:p3f),2);

m1_patdev_p3 = nanmean(patdev(:,p3s:p3f),2);


%And their interaction significance

meanp3col_all = cat(1, m1_constd_p3, m1_condev_p3, m1_patstd_p3, m1_patdev_p3);

X = cat(2, meanp3col_all, F1all, F2all, Sall);


[SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);


pgrp = Ps(1);
pcond = Ps(3);
pint = Ps(4);

p3pall(i,1)=pgrp;
p3pall(i,2)=pcond;
p3pall(i,3)=pint;


%Can push out into rainclouds for plotting in R :)
p3table = table(meanp3col_all, F1all, F2all, Sall);

p3table.Properties.VariableNames = {'meanp3col', 'Group', 'Trial', 'ID'};

writetable(p3table, [FigOutDir '/' OUTpre '_' 'ERP_'  'p3mean' '_' 'ConsandPats_' varnm{i} '_fullsample.txt'], 'Delimiter','tab')
    


%Add difference score for use in JASP

conp3diff = nanmean(constd(:,p3s:p3f),2)-nanmean(condev(:,p3s:p3f),2);

patp3diff = nanmean(patstd(:,p3s:p3f),2)-nanmean(patdev(:,p3s:p3f),2);


meanp3diffcol_all = cat(1, conp3diff, patp3diff);

meanp3diffcol_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));


p3difftable = table(meanp3diffcol_all, meanp3diffcol_grp);

p3difftable.Properties.VariableNames = {'meanp3col', 'Group'};

writetable(p3difftable, [FigOutDir '/' OUTpre '_' 'ERF_'  'p3diffmean' '_' 'ConsandPats_' varnm{i} '_fullsample.txt'], 'Delimiter','tab')


end


%Save output table of mean MMN

summrestable = table(p3pall(:,1), p3pall(:,2), p3pall(:,3));
summrestable.Properties.VariableNames = {'Group' 'Cond' 'Int'};

writetable(summrestable, [FigOutDir '/' OUTpre '_' 'ERP_'  'p3mean' '_restable_2wayANOVA_CondGrps_wInt.txt']);



%% And do again for rep3 (yikes)

MMNstats = [];

for i=1:9
    
    %Load in data (guess what, again)
    
    data=datarms{i}(:,:,:);
    
    
    %Mean amplitude for std and dev
    
    constd = data(proc_ctrs, std3trl, :);
    
    condev = data(proc_ctrs, devtrl, :);
    
    
    patstd = data(proc_pats, std3trl, :);
    
    patdev = data(proc_pats, devtrl, :);
    
    
    m1_constdmmn = nanmean(constd(:,MMNs:MMNf),2);
    
    m1_condevmmn = nanmean(condev(:,MMNs:MMNf),2);
    
    
    m1_patstdmmn = nanmean(patstd(:,MMNs:MMNf),2);
    
    m1_patdevmmn = nanmean(patdev(:,MMNs:MMNf),2);
    
    
    %And their interaction significance
    
%     meanMMNcol_all = cat(1, m1_constdmmn, m1_condevmmn, m1_patstdmmn, m1_patdevmmn);
%     
%     X = cat(2, meanMMNcol_all, F1all, F2all, Sall);
%     
%     
%     [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
%     
%     
%     pgrp = Ps(1);
%     pcond = Ps(3);
%     pint = Ps(4);
%     
%     MMNpall(i,1)=pgrp;
%     MMNpall(i,2)=pcond;
%     MMNpall(i,3)=pint;
    
    
    %Can push out into rainclouds for plotting in R :)
%     MMNtable = table(meanMMNcol_all, F1all, F2all, Sall);
%     
%     MMNtable.Properties.VariableNames = {'meanMMNcol', 'Group', 'Trial', 'ID'};
%     
%     writetable(MMNtable, [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNmean' '_rep3_' 'ConsandPats_' varnm{i} '_fullsample.txt'], 'Delimiter','tab')
%     
    
    
%Add difference score for use in JASP

conmmndiff = nanmean(constd(:,MMNs:MMNf),2)-nanmean(condev(:,MMNs:MMNf),2);

patmmndiff = nanmean(patstd(:,MMNs:MMNf),2)-nanmean(patdev(:,MMNs:MMNf),2);



[~,p,~,stats] = ttest2(conmmndiff,patmmndiff,'tail','both');

MMNpall(i,1) = p;
MMNtall(i,1) = stats.tstat;



%Collate for output table

meanMMNdiffcol_all = cat(1, conmmndiff, patmmndiff);

meanMMNdiffcol_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1])); 


MMNdifftable = table(meanMMNdiffcol_all, meanMMNdiffcol_grp, patgrpst);

MMNdifftable.Properties.VariableNames = {'meanMMNcol', 'Group', 'PatSubgrp'};

writetable(MMNdifftable, [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiffmean' '_rep3_' 'ConsandPats_' varnm{i} '_fullsample.txt'], 'Delimiter','tab')



    %Write out data for JASP Bayes RM ANOVA - would be too complex to write it out
    %containing multiple nodes
    
%     meanMMNcol_lw = cat(2, m1_constdmmn, m1_condevmmn);
%     
%     meanMMNcol_lw_pat = cat(2, m1_patstdmmn, m1_patdevmmn);
%     
%     meanMMNcol_lw_all = cat(1, meanMMNcol_lw, meanMMNcol_lw_pat);
% 
%     meanMMNcol_lw_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));
%     
%     MMNtable_lw = table(meanMMNcol_lw_all(:,1), meanMMNcol_lw_all(:,2), meanMMNcol_lw_grp);
%     
%     
%     MMNtable_lw.Properties.VariableNames = {'meanMMNcol_STD', 'meanMMNcol_Dev', 'Group'};
% 
%     
%     writetable(MMNtable_lw, [FigOutDir '/' OUTpre '_' 'ERF_'  'STDDEV' '_rep3_' 'ConsandPats_' varnm{i} '_fullsample_forJASP.txt'], 'Delimiter','tab')
    
   
end

%Save output table of mean MMN

MMNstats = struct;

MMNstats.t = MMNtall;

MMNstats.p = MMNpall;


save([FigOutDir '/' OUTpre '_' 'ERF_'  'MMNmean' '_rep3_' 'restable_ConsandPats_ttest.mat'], 'MMNstats');



%Save output table of mean MMN

% summrestable = table(MMNpall(:,1), MMNpall(:,2), MMNpall(:,3));
% summrestable.Properties.VariableNames = {'Group' 'Cond' 'Int'};
% 
% writetable(summrestable, [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNmean' '_rep3_' 'restable_2wayANOVA_CondGrps_wInt.txt']);




%% Now p300
for i=1:9
    
    %Load in data (guess what, again)
    
    data=datarms{i}(:,:,:);
    
    
    %Mean amplitude for std and dev
    
    constd = data(proc_ctrs, std3trl, :);
    
    condev = data(proc_ctrs, devtrl, :);
    
    
    patstd = data(proc_pats, std3trl, :);
    
    patdev = data(proc_pats, devtrl, :);
    
    
    m1_constd_p3 = nanmean(constd(:,p3s:p3f),2);
    
    m1_condev_p3 = nanmean(condev(:,p3s:p3f),2);
    
    
    m1_patstd_p3 = nanmean(patstd(:,p3s:p3f),2);
    
    m1_patdev_p3 = nanmean(patdev(:,p3s:p3f),2);
    
    
    %And their interaction significance
    
    meanp3col_all = cat(1, m1_constd_p3, m1_condev_p3, m1_patstd_p3, m1_patdev_p3);
    
    X = cat(2, meanp3col_all, F1all, F2all, Sall);
    
    
    [SSQs, DFs, MSQs, Fs, Ps] = mixed_between_within_anova(X);
    
    
    pgrp = Ps(1);
    pcond = Ps(3);
    pint = Ps(4);
    
    p3pall(i,1)=pgrp;
    p3pall(i,2)=pcond;
    p3pall(i,3)=pint;
    
    
    %Can push out into rainclouds for plotting in R :)
    p3table = table(meanp3col_all, F1all, F2all, Sall);
    
    p3table.Properties.VariableNames = {'meanp3col', 'Group', 'Trial', 'ID'};
    
    writetable(p3table, [FigOutDir '/' OUTpre '_' 'ERP_'  'p3mean' '_rep3_' 'ConsandPats_' varnm{i} '_fullsample.txt'], 'Delimiter','tab')
    
    
    
    %Add difference score for use in JASP
    
    conp3diff = nanmean(constd(:,p3s:p3f),2)-nanmean(condev(:,p3s:p3f),2);
    
    patp3diff = nanmean(patstd(:,p3s:p3f),2)-nanmean(patdev(:,p3s:p3f),2);
    
    
    meanp3diffcol_all = cat(1, conp3diff, patp3diff);
    
    meanp3diffcol_grp = cat(1, repmat(1, [nfilesc 1]), repmat(2, [nfilesp 1]));
    
    
    p3difftable = table(meanp3diffcol_all, meanp3diffcol_grp);
    
    p3difftable.Properties.VariableNames = {'meanp3col', 'Group'};
    
    writetable(p3difftable, [FigOutDir '/' OUTpre '_' 'ERF_'  'p3diffmean' '_rep3_' 'ConsandPats_' varnm{i} '_fullsample.txt'], 'Delimiter','tab')
    
    
end


%Save output table of mean MMN

summrestable = table(p3pall(:,1), p3pall(:,2), p3pall(:,3));
summrestable.Properties.VariableNames = {'Group' 'Cond' 'Int'};

writetable(summrestable, [FigOutDir '/' OUTpre '_' 'ERP_'  'p3mean' '_rep3_' '_restable_2wayANOVA_CondGrps_wInt.txt']);




%%=========================================================================
%% Disease Group
%%=========================================================================

%For simplicity define them here first

proc_pats_bv = find(patgrpst==1);
proc_pats_psp = find(patgrpst==2);
proc_cons = find(patgrpst==3);


%Rep 6

for i = 1:6
    
    
    figcount = figcount + 1;
    
    
    %Load in data (guess what, again)
    
    data=datarms{i}(proc_cons,:,:);
    
    
    data_pat_bv=datarms{i}(proc_pats_bv,:,:);
    
    
    data_pat_psp=datarms{i}(proc_pats_psp,:,:);
    
    
    
    %Std minus diff in placebo
    
    diff = squeeze(data(:,stdtrl,:)-data(:,devtrl,:));
    
    diff_pat_bv = squeeze(data_pat_bv(:,stdtrl,:)-data_pat_bv(:,devtrl,:));
    
    diff_pat_psp = squeeze(data_pat_psp(:,stdtrl,:)-data_pat_psp(:,devtrl,:));
    
    
    
    
    %Means and stds of difference (differences)
    
    m1diff=squeeze(nanmean(diff,1));
    s1diff=squeeze(nanstd(diff)./sqrt(size(diff,1)));
    
    
    m1diff_pat_bv=squeeze(nanmean(diff_pat_bv,1));
    s1diff_pat_bv=squeeze(nanstd(diff_pat_bv)./sqrt(size(diff_pat_bv,1)));
    
    
    m1diff_pat_psp=squeeze(nanmean(diff_pat_psp,1));
    s1diff_pat_psp=squeeze(nanstd(diff_pat_psp)./sqrt(size(diff_pat_psp,1)));
    
    
    %Note, independent ttest now just for pat subgroups
    [h,p,ci,stats]=ttest2(diff_pat_bv,diff_pat_psp);
    k=max(max([m1diff_pat_bv; m1diff_pat_psp]))+max(max([s1diff_pat_bv;s1diff_pat_psp]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 .25 .4]);
    plot(m1diff_pat_bv,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat_psp,'LineWidth',4,'Color',[0/255 255/255 0/255]); legendflex({'bvFTD','PSP'},'fontsize',14)
    boundedline([1:size(m1diff_pat_bv,2)],m1diff_pat_bv(1,:),s1diff_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_psp,2)],m1diff_pat_psp(1,:),s1diff_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('Std-Dev')
    title([varnm{i} ' MMN in bv/PSP']); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff' '_' 'PatsSubGrp_' varnm{i} '_fullsample.tif']);
    
    
    
    %Now including three groups (i.e. add controls)
    %Can go straight into analysis
    
    diff_allgrps = cat(1, diff, diff_pat_bv, diff_pat_psp);
    
    
    %Note, one-way anova now
    for t_win = 1:length(diff_allgrps(1,:))
        p(1,t_win) = anova1(diff_allgrps(:,t_win),patgrpst, 'off');
    end
    
    k=max(max([m1diff; m1diff_pat_bv; m1diff_pat_psp]))+max(max([s1diff; s1diff_pat_bv; s1diff_pat_psp])); k=k*1.1;
    
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %start figure
    
    figcount = figcount + 1;
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 .25 .4]);
    plot(m1diff,'LineWidth',4,'Color',[0 0 255]./255); hold on
    plot(m1diff_pat_bv,'LineWidth',4,'Color',[255 0 0]./255);
    plot(m1diff_pat_psp,'LineWidth',4,'Color',[0/255 255/255 0/255]); legendflex({'CON', 'bvFTD','PSP'},'fontsize',14)
    
    boundedline([1:size(m1diff,2)],m1diff(1,:),s1diff(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_bv,2)],m1diff_pat_bv(1,:),s1diff_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_psp,2)],m1diff_pat_psp(1,:),s1diff_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('Std-Dev')
    title([varnm{i} ' MMN']); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff' '_' 'ConANDPatsSubGrp_' varnm{i} '_fullsample.tif']);
    
    
end


%Rep3


for i = 1:6
    
    
    figcount = figcount + 1;
    
    
    %Load in data (guess what, again)
    
    data=datarms{i}(proc_cons,:,:);
    
    
    data_pat_bv=datarms{i}(proc_pats_bv,:,:);
    
    
    data_pat_psp=datarms{i}(proc_pats_psp,:,:);
    
    
    
    %Std minus diff in placebo
    
    diff = squeeze(data(:,std3trl,:)-data(:,devtrl,:));
    
    diff_pat_bv = squeeze(data_pat_bv(:,std3trl,:)-data_pat_bv(:,devtrl,:));
    
    diff_pat_psp = squeeze(data_pat_psp(:,std3trl,:)-data_pat_psp(:,devtrl,:));
    
    
    
    
    %Means and stds of difference (differences)
    
    m1diff=squeeze(nanmean(diff,1));
    s1diff=squeeze(nanstd(diff)./sqrt(size(diff,1)));
    
    
    m1diff_pat_bv=squeeze(nanmean(diff_pat_bv,1));
    s1diff_pat_bv=squeeze(nanstd(diff_pat_bv)./sqrt(size(diff_pat_bv,1)));
    
    
    m1diff_pat_psp=squeeze(nanmean(diff_pat_psp,1));
    s1diff_pat_psp=squeeze(nanstd(diff_pat_psp)./sqrt(size(diff_pat_psp,1)));
    
    
    %Note, independent ttest now just for pat subgroups
    [h,p,ci,stats]=ttest2(diff_pat_bv,diff_pat_psp);
    k=max(max([m1diff_pat_bv; m1diff_pat_psp]))+max(max([s1diff_pat_bv;s1diff_pat_psp]));k=k*1.1;
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 .25 .4]);
    plot(m1diff_pat_bv,'LineWidth',4,'Color',[255 0 0]./255); hold on
    plot(m1diff_pat_psp,'LineWidth',4,'Color',[0/255 255/255 0/255]); legendflex({'bvFTD','PSP'},'fontsize',14)
    boundedline([1:size(m1diff_pat_bv,2)],m1diff_pat_bv(1,:),s1diff_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_psp,2)],m1diff_pat_psp(1,:),s1diff_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('Std-Dev')
    title([varnm{i} ' MMN3 in bv/PSP']); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff' '_rep3_' '_' 'PatsSubGrp_' varnm{i} '_fullsample.tif']);
    
    
    
    %Now including three groups (i.e. add controls)
    %Can go straight into analysis
    
    diff_allgrps = cat(1, diff, diff_pat_bv, diff_pat_psp);
    
    
    %Note, one-way anova now
    for t_win = 1:length(diff_allgrps(1,:))
        p(1,t_win) = anova1(diff_allgrps(:,t_win),patgrpst, 'off');
    end
    
    k=max(max([m1diff; m1diff_pat_bv; m1diff_pat_psp]))+max(max([s1diff; s1diff_pat_bv; s1diff_pat_psp])); k=k*1.1;
    
    p=k*(p<0.05); p(find(p==0))=NaN; p(1:51)=NaN;
    
    
    %start figure
    
    figcount = figcount + 1;
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 .25 .4]);
    plot(m1diff,'LineWidth',4,'Color',[0 0 255]./255); hold on
    plot(m1diff_pat_bv,'LineWidth',4,'Color',[255 0 0]./255);
    plot(m1diff_pat_psp,'LineWidth',4,'Color',[0/255 255/255 0/255]); legendflex({'CON', 'bvFTD','PSP'},'fontsize',14)
    
    boundedline([1:size(m1diff,2)],m1diff(1,:),s1diff(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_bv,2)],m1diff_pat_bv(1,:),s1diff_pat_bv(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1diff_pat_psp,2)],m1diff_pat_psp(1,:),s1diff_pat_psp(1,:),'cmap',[0 255 0]./255, 'transparency', 0.3, 'orientation', 'vert', 'alpha');
    plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',16); xlabel('Time (ms)'); ylabel('Std-Dev')
    title([varnm{i} ' MMN3']); xlim([0 251]); %ylim([0 1]);%axis square
    
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'MMNdiff' '_rep3_' '_' 'ConANDPatsSubGrp_' varnm{i} '_fullsample.tif']);
    
    
end


%Write out metadata

writecell(allsubjs, [FigOutDir '/' OUTpre '_' 'ERF_' 'PlacGroupDifferenceResults_' date '_subjlist.txt']); 



%% Plot trial waveforms of interest

i = 4; %Gradiometers
    

    data=datarms{i}(proc_ctrs,:,:);
    
    m1=squeeze(nanmean(data,1));
    s1=squeeze(nanstd(data)./sqrt(size(data,1)));
    
    
    data_pat=datarms{i}(proc_pats,:,:);
    
    m1_pat=squeeze(nanmean(data_pat,1));
    s1_pat=squeeze(nanstd(data_pat)./sqrt(size(data_pat,1)));

    
    %Deviants and rep3
    
    figcount = figcount+1;
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 0.25 0.25]);
    plot(m1(2,:),'--','LineWidth',2,'Color',[0/255 0/255 255/255]); hold on
    plot(m1(1,:),':','LineWidth',2,'Color',[0/255 0/255 255/255]);
    plot(m1_pat(2,:),'--','LineWidth',2,'Color',[255/255 0/255 0/255]);
    plot(m1_pat(1,:),':','LineWidth',2,'Color',[255/255 0/255 0/255]); %legendflex({'rep3 CON','dev CON','rep3 FTLD','dev FTLD'},'fontsize',14)
    boundedline([1:size(m1,2)],m1(2,:),s1(2,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat(2,:),s1_pat(2,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat(1,:),s1_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial'); ylabel({'fT/mm'}, 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    %title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
        box off

    set(gca, 'TickDir', 'Out')
    
    xlim([0 251]);
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'ControlsandPats' '_' varnm{i} 'rep3anddev_fullsample.tif']);

    
    %Deviants and rep6
    
    figcount = figcount + 1;
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 0.25 0.25]);
    plot(m1(3,:),'--','LineWidth',2,'Color',[0/255 0/255 255/255]); hold on
    plot(m1(1,:),':','LineWidth',2,'Color',[0/255 0/255 255/255]);
    plot(m1_pat(3,:),'--','LineWidth',2,'Color',[255/255 0/255 0/255]);
    plot(m1_pat(1,:),':','LineWidth',2,'Color',[255/255 0/255 0/255]); %legendflex({'rep6 CON','dev CON','rep6 FTLD','dev FTLD'},'fontsize',14)
    boundedline([1:size(m1,2)],m1(3,:),s1(3,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat(3,:),s1_pat(3,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat(1,:),s1_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial'); ylabel({'fT/mm'}, 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    %title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
        box off

    set(gca, 'TickDir', 'Out')
    
    xlim([0 251]);
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'ControlsandPats' '_' varnm{i} 'rep6anddev_fullsample.tif']);

    
    %Add in legends sep
    
    %Deviants and rep3
    
    figcount = figcount+1;
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 0.25 0.25]);
    plot(m1(2,:),'--','LineWidth',2,'Color',[0/255 0/255 255/255]); hold on
    plot(m1(1,:),':','LineWidth',2,'Color',[0/255 0/255 255/255]);
    plot(m1_pat(2,:),'--','LineWidth',2,'Color',[255/255 0/255 0/255]);
    plot(m1_pat(1,:),':','LineWidth',2,'Color',[255/255 0/255 0/255]); legendflex({'rep3 CON','dev CON','rep3 FTLD','dev FTLD'},'fontsize',12, 'fontweight', 'bold')
    boundedline([1:size(m1,2)],m1(2,:),s1(2,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat(2,:),s1_pat(2,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat(1,:),s1_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial'); ylabel({'fT/mm'}, 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    %title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
        box off

    set(gca, 'TickDir', 'Out')
    
    xlim([0 251]);
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'ControlsandPats' '_' varnm{i} 'rep3anddev_fullsample_wleg.tif']);

    
    %Deviants and rep6
    
    figcount = figcount + 1;
    
    figure(figcount);set(gcf, 'Units','normalized', 'Position', [0 0 0.25 0.25]);
    plot(m1(3,:),'--','LineWidth',2,'Color',[0/255 0/255 255/255]); hold on
    plot(m1(1,:),':','LineWidth',2,'Color',[0/255 0/255 255/255]);
    plot(m1_pat(3,:),'--','LineWidth',2,'Color',[255/255 0/255 0/255]);
    plot(m1_pat(1,:),':','LineWidth',2,'Color',[255/255 0/255 0/255]); legendflex({'rep6 CON','dev CON','rep6 FTLD','dev FTLD'},'fontsize',12, 'fontweight', 'bold')
    boundedline([1:size(m1,2)],m1(3,:),s1(3,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1,2)],m1(1,:),s1(1,:),'cmap',[0/255 0/255 255/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat(3,:),s1_pat(3,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    boundedline([1:size(m1_pat,2)],m1_pat(1,:),s1_pat(1,:),'cmap',[255/255 0/255 0/255], 'transparency', 0.1, 'orientation', 'vert', 'alpha');
    %plot(p,'-o','LineWidth',2,'MarkerEdgeColor',[0.6 0.6 0.6],'MarkerFaceColor',[0.6 0.6 0.6],'MarkerSize',8)
    set(gca, 'XTick',xticks, 'XTickLabel',xlabels,'Fontsize',12, 'FontWeight', 'bold', 'FontName', 'Arial'); xlabel('Time (ms)', 'FontWeight', 'bold', 'FontSize', 12, 'FontName', 'Arial'); ylabel({'fT/mm'}, 'FontWeight', 'bold', 'FontSize', 13, 'FontName', 'Arial');
    %title([varnm{i} ' mean in controls']); xlim([0 251]); %ylim([0 1]);%axis square
    
        box off

    set(gca, 'TickDir', 'Out')
    
    xlim([0 251]);
    
    %Save figure
    saveas(figure(figcount), [FigOutDir '/' OUTpre '_' 'ERF_'  'ControlsandPats' '_' varnm{i} 'rep6anddev_fullsample_wleg.tif']);

    
    
end