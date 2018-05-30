% Example code to compute spatial modes, and run statistical analyses of
% mode variance. The code is by no means optimized, for the sake of
% interpretability, but should run relatively quickly nonetheless. I wrote
% this code in MATLAB 2012a, and it should be forward compatible. 
%
% The code is written so as to produce modes that are more strongly
% expressed in the atomoxetine condition than in the placebo condition, and
% compare the variance of the first mode between conditions using
% cross-validation (both the average and with ROC analysis). 
%
% Ruud van den Brink, 2017
%
% If you use any of this this code for a publication, please cite: 
% van den Brink, Nieuwenhuis, & Donner (2018) Amplification and Suppression
% of Distinct Brain-wide Activity Patterns by Catecholamines  

%% clear contents
clear 
close all
warning('off','all')
clc

%% load the data
load M.mat %this matrix contains Z-scored BOLD time series, and is of size: participants X condition (1=placebo, 2=drug) X volume X brain region

%% the number of iterations for permutation testing
npermutes = 10000;

%% compute covariance matrices (eq. 1 in the article)

%initialize covariance matrices. Size: participants X condition X brain region X brain region
C  = zeros(size(M,1),size(M,2),size(M,4),size(M,4)); %covariance across all the volumes
C1 = zeros(size(M,1),size(M,2),size(M,4),size(M,4)); %covariance across the first half of volumes (for statistics later on)
C2 = zeros(size(M,1),size(M,2),size(M,4),size(M,4)); %covariance across the second half of volumes (for statistics later on)  

for subi = 1:size(M,1) %loop over participants
    for condi = 1:size(M,2) %loop over conditions        
        C(subi,condi,:,:)  = cov(squeeze(M(subi,condi,:,:))); %covariance across all the volumes
        C1(subi,condi,:,:) = cov(squeeze(M(subi,condi,1:size(M,3)/2,:))); %covariance across the first half of volumes (for statistics later on)
        C2(subi,condi,:,:) = cov(squeeze(M(subi,condi,ceil(size(M,3)/2):size(M,3),:))); %covariance across the second half of volumes (for statistics later on)        
    end
end

%covariance for all the data
Cp = squeeze(mean(C(:,1,:,:),1)); %covariance placebo, averaged across participants
Ca = squeeze(mean(C(:,2,:,:),1)); %covariance atomoxetine, averaged across participants

%covariance for first half of data (for cross-validation)
Cp1 = squeeze(mean(C1(:,1,:,:),1)); %covariance placebo, averaged across participants
Ca1 = squeeze(mean(C1(:,2,:,:),1)); %covariance atomoxetine, averaged across participants

%covariance for second half of data (for cross-validation)
Cp2 = squeeze(mean(C2(:,1,:,:),1)); %covariance placebo, averaged across participants
Ca2 = squeeze(mean(C2(:,2,:,:),1)); %covariance atomoxetine, averaged across participants
%% decompose based on the full dataset (visualization only)

%this is where the decomposition happens
[V,lambda] = eig(Ca,Cp); %decompose (atomoxetine > placebo)
% [V,lambda] = eig(Cp,Ca); %decompose (placebo > atomoxetine) (eq. 4 in the article)

[~, I] = sort(diag(lambda),'descend'); %sort in decending order
lambda = lambda(I,I); %sort in decending order
V = V(:,I); %sort V by lambda
l_sum = sum(diag(lambda)); %get the sum of eigenvalues

%Compute mode time series
%These you would use as regressors onto the voxel-level data to produce
%mode spatial maps. 
ti = zeros(size(M));
for subi = 1:size(M,1) %loop over participants
    for condi = 1:2 %loop over conditions        
        for modei = 1:size(V,2) %loop over modes
            %project the mode onto the data
            p  = V(:,modei); %spatial mode (ROI by 1)
            mi = squeeze(M(subi,condi,:,:)); %data for this participant and this run
            %time-course correspinding to spatial mode, size: participants X condition X volume X mode
            ti(subi,condi,:,modei) = mi*p; %(eq. 5 in the article)
        end
    end
end

%% make some plots (covariance matrices, V, lambda, and a bar plot of lambda)
figure
subplot(2,2,1)
imagesc(Cp)
title('C_P')
set(gca,'xtick',[],'ytick',[],'clim',[-1 1])
axis square

subplot(2,2,2)
imagesc(Ca)
title('C_A')
set(gca,'xtick',[],'ytick',[],'clim',[-1 1])
axis square

subplot(2,2,3)
imagesc(V)
title('Eigenvectors (V)')
set(gca,'xtick',[],'ytick',[],'clim',[-1 1])
axis square

subplot(2,2,4)
imagesc(lambda)
title('Eigenvalues (\lambda)')
set(gca,'xtick',[],'ytick',[],'clim',[-2.3 2.3])
axis square

colormap(cmap);

figure
box off
bar(diag(lambda),'k')
xlabel('Mode number')
ylabel('Eigenvalue (a.u.)')
xlim([0 size(lambda,1)])
set(gca,'tickdir','out')
box off
%% define modes on half the data and test for differences in variance in the other half of the data

%The code below generates modes that are more strongly expressed in the
%atomoxetine condition than in the placebo condition (i.e. the mode
%captures a pattern of brain regions in which activity cofluctuates more
%strongly following atomoxetine). To run the decomposition in the reverse
%direction (i.e. find patterns in which activity cofluctuates less
%strongly), swap 'Ca' and Cp' in the relevant lines below. Also, change
%the tail of the ptest at the end of the cell. You would also need to make
%the neccesary changes to the code for ROC analysis (i.e. swap the inputs
%in the decomposition lines, and swap the inputs in the ROC analysis
%itself)
%
% Modes here are computed based on covariance in half of the data, and then
% projected onto the other (independent) half. This will tell us if the
% mode captures meaningful changes in cofluctuation strength rather than
% explain more variance in one condition than the other by definition. 

%decompose based on first half of data (eq. 4 in the article)
[V1,lambda1] = eig(Ca1, Cp1); % (atomoxetine > placebo)
% [V1,lambda1] = eig(Cp1, Ca1); % (placebo > atomoxetine)  
[~, I] = sort(diag(lambda1),'descend'); %sort in decending order
lambda1 = lambda1(I,I); %sort in decending order
V1 = V1(:,I); %sort V by lambda
l_sum1 = sum(diag(lambda1)); %get the sum of eigenvalues

%decompose based on second half of data (eq. 4 in the article)
[V2,lambda2] = eig(Ca2, Cp2); % (atomoxetine > placebo)
% [V2,lambda2] = eig(Cp2, Ca2); % (placebo > atomoxetine)  
[~, I] = sort(diag(lambda2),'descend'); %sort in decending order
lambda2 = lambda2(I,I); %sort in decending order
V2 = V2(:,I); %sort V by lambda
l_sum2 = sum(diag(lambda2)); %get the sum of eigenvalues

%get first and second half of data to project modes onto
M1 = M(:,:,1:size(M,3)/2,:);
M2 = M(:,:,ceil(size(M,3)/2):size(M,3),:);

%initialize variables
ti1 = zeros(size(M,1),2,size(M2,3)); %variance explained in second half of data
ti2 = zeros(size(M,1),2,size(M1,3)); %variance explained in first half of data
s   = zeros(size(M,1),2); %average across the above two

%now compute the percentage of variance explained by the first mode
for subi = 1:size(M,1) %loop over participants
    for condi = 1:2
        m1 = squeeze(M1(subi,condi,:,:)); %first half of data for this participant and condition (ROI by TR)
        m2 = squeeze(M2(subi,condi,:,:)); %scond half of data for this participant and condition (ROI by nTRs)
        
        %project the mode (only the first) onto the data, separately for
        %each half of the data
        p  = V1(:,1); %spatial mode (ROI by 1)
        ti1(subi,condi,:) = m2*p; %time-course correspinding to spatial mode (nTRs by 1) (eq. 5 in the article)
        s1 = (squeeze(ti1(subi,condi,:))' * squeeze(ti1(subi,condi,:))) / l_sum1; %(eq. 6 in the article)
        
        p  = V2(:,1); %spatial mode (ROI by 1)
        ti2(subi,condi,:) = m1*p; %time-course correspinding to spatial mode (nTRs by 1) (eq. 5 in the article)
        s2 = (squeeze(ti2(subi,condi,:))' * squeeze(ti2(subi,condi,:))) / l_sum2; %(eq. 6 in the article)
        
        %average elements of s across first and second half of the data so
        %that we have one cross-validated value to test
        s(subi,condi) = mean([s1 s2]);        
    end
end

%test one-tailed because we expect variance to be greater in the
%atomoxetine condition (change tail when decomposing in the reverse 
%direction)
% [~, p] = permtest(squeeze(s(:,2)),squeeze(s(:,1)),npermutes,0.05,'left'); 
[~, p] = permtest(squeeze(s(:,2)),squeeze(s(:,1)),npermutes,0.05,'right');
disp(['Comparison of mode variance between drug and placebo: p=' num2str(p)])

%% bar plot of variance explained in the individual conditions

figure
bar(squeeze(mean(s)),'facecolor',[0.8 0.8 0.8],'edgecolor','none');
hold on

%add error bars
sems  = squeeze(std(s)./sqrt(size(M,1)))';
means = squeeze(mean(s))';
plot([1; 1], [means(1,:)-sems(1,:); means(1,:)+sems(1,:) ],'k','linewidth',2);
plot([2; 2], [means(2,:)-sems(2,:); means(2,:)+sems(2,:) ],'k','linewidth',2);

%formatting
xlim([0 3])
ylabel('Variance explained (%)')
title('Comparison of mode 1 variance')
set(gca,'tickdir','out','xtick',[1 2],'xticklabel',{'Placebo' 'Atomoxetine'})
box off; axis square

%indicate if the comparison was significant
if p < 0.05 && p > 0.01; text2plot = '*'; elseif p < 0.01 && p > 0.001; text2plot = '**'; elseif p < 0.001; text2plot = '***'; else; text2plot = 'n.s.'; end
text(1.5,max(means)+max(sems),text2plot,'HorizontalAlignment','center','fontsize',30,'VerticalAlignment','middle')

%% prepare the data for ROC analysis

%In thic cell, modes are computed on 4 different bins of data. The
%remaining part of the data is then segmented into 20 parts, and the modes
%are projected onto each of those parts. For each of those 20 segments, we
%compute mode variance. The resulting distribution of variances (4 x 20) is
%submitted to ROC analysis in the next cell, for each participant. 

nbins  = 4;  %how many bins of data to define modes on
nsecs  = 20; %how many sections to cut the remainder of TRs in, to project modes onto

%initialize variables
nTRs   = size(M,3); %how many volumes are in the data
binidx = 1:floor(nTRs/(nbins)):nTRs; %indices of the different bins
C_r = zeros(size(M,1),size(M,2),nbins,size(M,4),size(M,4)); %covariance matrices to define modes on, size: participat X condition X bin X brain region X brain region
M_r = zeros(size(M,1),size(M,2),nbins,nTRs-(binidx(2)-binidx(1)),size(M,4)); %data to project modes onto, size: participat X condition X bin X volume X brain region

%compute covariance matrices and save the volumes on which covariance was
%not computed
for subi = 1:size(M,1) %loop over participants
    %get covariance matrices of part of the data, to define modes on
    for condi = 1:2 %loop over conditions
        for bini = 1:nbins %loop over bins
            %for each bin, get the covariance of the data in that bin
            C_r(subi,condi,bini,:,:) = cov(zscore(squeeze(M(subi,condi,binidx(bini)+1:binidx(bini+1),:))));
            %then, get the data of all time-points outside of the current
            %bin
            tmptcs = squeeze(M(subi,condi,:,:));
            tmptcs(binidx(bini)+1:binidx(bini+1),:) = [];
            M_r(subi,condi,bini,:,:) = zscore(tmptcs);
        end
    end
end

%set up segmentation of the remaining volumes (the part of the data that was not
%included in defining the mode is cut up into pieces, and variance
%explained [si] is computed for each of those pieces)
rTRs = size(M_r,4); %the number of remaining volumes 
tridx = round(linspace(0 ,rTRs, nsecs+1)); %indices of the segments of the remainging volumes
ti = zeros(size(M,1),size(M,2),nbins,size(M_r,4)); %initialize
ss = zeros(size(M,1),size(M,2),nbins,nsecs); %his will contain mode variance for each bin and segment of the remaining data
%loop over bins, decompose for each bin, and compute mode variance for the
%20 segments of remaining TRs 
for bini = 1:nbins
     %decompose
    [V, lambda] = eig(squeeze(mean(C_r(:,2,bini,:,:),1)),squeeze(mean(C_r(:,1,bini,:,:),1))); % (atomoxetine > placebo)
%     [V, lambda] = eig(squeeze(mean(C_r(:,1,bini,:,:),1)),squeeze(mean(C_r(:,2,bini,:,:),1))); % (placebo > atomoxetine)
    [~, I] = sort(diag(lambda),'descend'); %sort in decending order
    lambda = lambda(I,I); %sort in decending order
    V = V(:,I); %sort V by lambda
    l_sum = sum(diag(lambda)); %get the sum of eigenvalues
    
    for subi = 1:size(M,1) %loop over participants
        for condi = 1:2 %loop over conditions
            mi = squeeze(M_r(subi,condi,bini,:,:)); %data time-courses (ROI by nTRs)
            p  = V(:,1); %spatial mode (ROI by 1)
            ti(subi,condi,bini,:) = mi*p; %time-course correspinding to spatial mode (nTRs by 1)
            
            %compute variance explained for each little section of the
            %remaining volumes
            for seci = 1:length(tridx)-1
                s = squeeze(ti(subi,condi,bini,tridx(seci)+1:tridx(seci+1)))' * squeeze(ti(subi,condi,bini,tridx(seci)+1:tridx(seci+1)));
                ss(subi,condi,bini,seci) = s; %collect variance for all the segments
            end %end seci
        end %end condi
    end %end subi
end %end bini


%% run ROC analysis

auc = zeros(size(M,1),nbins); %will contain the ROC index
for subi = 1:size(M,1)
    for bini = 1:nbins  
        %in the ROC analysis values of 0.5+ will indicate better
        %classification if mean(A) > mean(B), so when running the
        %decomposition in the reverse direction, swap A and B below
        A = squeeze(ss(subi,2,bini,:)); % (atomoxetine > placebo)
        B = squeeze(ss(subi,1,bini,:)); % (atomoxetine > placebo)
%         A = squeeze(ss(subi,1,bini,:)); % (placebo > atomoxetine)
%         B = squeeze(ss(subi,2,bini,:)); % (placebo > atomoxetine)  
        
        %run ROC analysis
        [auc(subi,bini), tp, fp] = sroc(A,B);
        
        %initialize
        if subi == 1 && bini == 1; tproc = zeros(size(M,1),nbins,length(tp)); fproc = zeros(size(M,1),nbins,length(fp)); end
        tproc(subi,bini,:) = tp;
        fproc(subi,bini,:) = fp;        
    end
end

mauc = squeeze(mean(mean(auc,1),2)); %mean area under the curve
sauc = std(squeeze(mean(auc,2)))./sqrt(size(M,1)); %SEM across participants

%% plot the results of the ROC analysis

figure

%plot the ROC index
subplot(1,2,2), hold on
bar(mauc,'facecolor',[0.5 0.5 0.5],'edgecolor','none')
axis square
ylim([0.4 0.8])
hold on
plot([0 2],[0.5 0.5],'k--')
xlim([0 2])
set (gca,'xtick',[],'tickdir','out')
ylabel('Area under curve (a.u.)')
hold on
plot([1; 1],[mauc'-sauc; mauc'+sauc],'k','linewidth',2)
[j, p] = permtest(squeeze(mean(auc,2)),0.5,npermutes,0.05,'right');

if p < 0.05 && p > 0.01; text2plot = '*'; elseif p < 0.01 && p > 0.001; text2plot = '**'; elseif p < 0.001; text2plot = '***'; else; text2plot = 'n.s.'; end
text(1,mauc+sauc+0.05,text2plot,'HorizontalAlignment','center','fontsize',30,'VerticalAlignment','middle')
disp(['ROC index: ' num2str(mauc) ', p=' num2str(p)])

title('AUC discriminating atomoxetine and placebo')

%plot the ROC curve itself
mtproc = squeeze(mean(tproc,2)); %mean true positive rate
mfproc = squeeze(mean(fproc,2)); %mean false positive rate
subplot(1,2,1), hold on, L=line([0 1],[0 1]); set(L,'LineStyle','--','color','k')
plot(squeeze(mean(mfproc(:,:),1)),squeeze(mean(mtproc(:,:),1)) ,'k' );
axis square
xlabel('1-Specificity (False positive rate)')
ylabel('Sensitivity (True positive rate)')
set (gca,'xtick',[0 .5 1],'ytick',[0 .5 1],'tickdir','out')
title('ROC curve')
