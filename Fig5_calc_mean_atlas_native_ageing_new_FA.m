clear all

close all

data_dir=' ';

metadata=importdata('metadata.txt');

rats={};

d1=dir(data_dir);

count=0;

param={'FA'}; %;'FR'

for i=1:size(metadata.textdata,1)

    rats{i}=metadata.textdata{i,1};

end

gm=1:70;

wm=[1 4 5 19 20 29 32 42 46 59 61];

PFC=[55 56 67 68];

for i=wm

    gm(gm==i)=[];

end

template=load_nii('AtlasV4.0.nii');atlas1=template.img;

atlas1(:,:,1:40)=0; ant_gm=unique(atlas1);

atlas1=template.img; atlas1(:,:,30:end)=0; post_gm=unique(atlas1);

atlas1=template.img; atlas1(:,:,[1:29 39:64])=0; med_gm=unique(atlas1);

rois=med_gm(2:end)'; % pick rois
%
rois=PFC;

for i=1:size(wm,2)
    rois(rois==wm(i))=[];
end

rois=wm;


for k=1:length(rats)

    for p=1:length(param)

        try

            atlas=load_nii([data_dir 'atlas_FA_' rats{k} '.nii' ] ); atlas=atlas.img;

            V=load_nii([data_dir  param{p} '_' rats{k} '.nii']); V=V.img;

            for i=1:max(atlas(:))

                tabl(k,p,i)=mean(V(atlas==i));


            end

        catch

            for i=1:max(atlas(:))

                tabl(k,p,i)=nan;

            end
        end

    end

end

exper=metadata.data(:,3);

age=metadata.data(:,4);

pnd=metadata.data(:,5);

subj=metadata.data(:,6);

count=0;


for r=rois ; 

    count=count+1;

    for p=1:length(param)

        count1=0;

        for i=[120 180 400 520]

            count1=count1+1;

            ag=i;

            t=squeeze(tabl(:,p,r));

        end

        yy=t(exper==1&age<400);

        xx=pnd(exper==1&age<400);

        ss=subj(exper==1&age<400);

        cc=0;

        for i=unique(ss)'

            cc=cc+1;

            if size(xx(ss==i),1)>1

                p=polyfit(xx(ss==i),yy(ss==i),1);

                sl(cc)=p(1);

            else

                sl(cc)=nan;

            end

        end

        s1_1(count,:)=sl;

        yy=t(exper==1&age>180);

        xx=pnd(exper==1&age>180);

        ss=subj(exper==1&age>180);

        cc=0; sl=[];

        for i=unique(ss)'

            cc=cc+1;

            if size(xx(ss==i),1)>1

                p=polyfit(xx(ss==i),yy(ss==i),1);

                sl(cc)=p(1);

            else

                sl(cc)=nan;

            end

        end

        s1_2(count,:)=sl;

        yy=t(exper==2&age<400);

        xx=pnd(exper==2&age<400);

        ss=subj(exper==2&age<400);

        cc=0; sl=[];

        for i=unique(ss)'

            cc=cc+1;

            if size(xx(ss==i),1)>1

                p=polyfit(xx(ss==i),yy(ss==i),1);

                sl(cc)=p(1);

            else

                sl(cc)=nan;

            end

        end

        s2_1(count,:)=sl;

        yy=t(exper==2&age>180);

        xx=pnd(exper==2&age>180);

        ss=subj(exper==2&age>180);

        cc=0; sl=[];

        for i=unique(ss)'

            cc=cc+1;

            if size(xx(ss==i),1)>1

                p=polyfit(xx(ss==i),yy(ss==i),1);

                sl(cc)=p(1);

            else

                sl(cc)=nan;

            end

        end

        s2_2(count,:)=sl;

        yy=t(exper==3&age<400);

        xx=pnd(exper==3&age<400);

        ss=subj(exper==3&age<400);

        cc=0; sl=[];

        for i=unique(ss)'

            cc=cc+1;

            if size(xx(ss==i),1)>1

                p=polyfit(xx(ss==i),yy(ss==i),1);

                sl(cc)=p(1);

            else

                sl(cc)=nan;

            end

        end

        s3_1(count,:)=sl;

        yy=t(exper==3&age>180);

        xx=pnd(exper==3&age>180);

        ss=subj(exper==3&age>180);

        cc=0; sl=[];

        for i=unique(ss)'

            cc=cc+1;

            if size(xx(ss==i),1)>1

                p=polyfit(xx(ss==i),yy(ss==i),1);

                sl(cc)=p(1);

            else

                sl(cc)=nan;

            end

        end

        s3_2(count,:)=sl;


    end

 
end

%-------------------
% ANOVA timepoint 1
%-------------------

% Dimensions
[nROI,n1] = size(s1_1);
[~,n2] = size(s2_1);
[~,n3] = size(s3_1);

% Data
slope = [s1_1(:);s2_1(:);s3_1(:)];

% Group
Group = [repmat("G1",nROI*n1,1);repmat("G2",nROI*n2,1);repmat("G3",nROI*n3,1)];

% ROI
ROI = [repmat((1:nROI)',n1,1);repmat((1:nROI)',n2,1);repmat((1:nROI)',n3,1)];

% Subject ID
Subject = [repelem((1:n1)',nROI);repelem((n1+1:n1+n2)',nROI);repelem((n1+n2+1:n1+n2+n3)',nROI)];

% Table
T = table(slope, categorical(Group), categorical(ROI), ...
          categorical(Subject), ...
          'VariableNames',{'Slope','Group','ROI','Subject'});

% Linear mixed-effects model
lme = fitlme(T,'Slope ~ Group*ROI + (1|Subject)');

% ANOVA
anova(lme)

%figure
% Mean slope across ROIs for each subject
meanSlope1 = mean(s1_1,1,'omitnan');
meanSlope2 = mean(s2_1,1,'omitnan');
meanSlope3 = mean(s3_1,1,'omitnan');

data = [meanSlope1(:); meanSlope2(:); meanSlope3(:)];

groups = [ones(5,1);
          2*ones(5,1);
          3*ones(5,1)];

figure;
hold on;

% Individual subjects
x = [ones(5,1); 2*ones(5,1); 3*ones(5,1)];
xjitter = x + (rand(size(x))-0.5)*0.12;

scatter(xjitter,data,50,'k','filled');

% Mean ± SEM
means = [mean(meanSlope1,'omitnan'), ...
         mean(meanSlope2,'omitnan'), ...
         mean(meanSlope3,'omitnan')];

sem = [std(meanSlope1,'omitnan')/sqrt(sum(~isnan(meanSlope1))), ...
       std(meanSlope2,'omitnan')/sqrt(sum(~isnan(meanSlope2))), ...
       std(meanSlope3,'omitnan')/sqrt(sum(~isnan(meanSlope3)))];

errorbar(1:3,means,sem,'k','LineWidth',1.5,...
         'LineStyle','none','CapSize',10);

% Formatting
xlim([0.5 3.5]);
xticks(1:3);
xticklabels({'Group 1','Group 2','Group 3'});
ylabel('Mean slope across ROIs');

box off;
set(gca,'FontSize',12);

for r = 1:nROI
    
    data = [s1_1(r,:)'; s2_1(r,:)'; s3_1(r,:)'];
    
    group = [repmat({'G1'},5,1);
             repmat({'G2'},5,1);
             repmat({'G3'},5,1)];
    
    [p(r),~,stats] = anova1(data,group,'off');
    
    if p(r) < 0.05
        c = multcompare(stats,'Display','off');
        posthoc{r} = c;
    end
end

% FDR correction across the 11 ROIs
pFDR = mafdr(p,'BHFDR',true);

% if no roi effect
Tg = table(T.Slope, T.Group, T.Subject, ...
    'VariableNames',{'Y','Group','Subject'});

lme_g = fitlme(Tg,'Y ~ Group + (1|Subject)');
anova(lme_g)
[p,~,stats] = anova1(Tg.Y, Tg.Group, 'off');
c = multcompare(stats,'CType','tukey-kramer');

%-------------------
% ANOVA timepoint 2
%-------------------

% Dimensions
[nROI,n1] = size(s1_2);
[~,n2] = size(s2_2);
[~,n3] = size(s3_2);

% Data
slope = [s1_2(:);s2_2(:);s3_2(:)];

% Group
Group = [repmat("G1",nROI*n1,1);repmat("G2",nROI*n2,1);repmat("G3",nROI*n3,1)];

% ROI
ROI = [repmat((1:nROI)',n1,1);repmat((1:nROI)',n2,1);repmat((1:nROI)',n3,1)];

% Subject ID
Subject = [repelem((1:n1)',nROI);repelem((n1+1:n1+n2)',nROI);repelem((n1+n2+1:n1+n2+n3)',nROI)];

% Table
T = table(slope, categorical(Group), categorical(ROI), ...
          categorical(Subject), ...
          'VariableNames',{'Slope','Group','ROI','Subject'});

% Linear mixed-effects model
lme = fitlme(T,'Slope ~ Group*ROI + (1|Subject)');

% ANOVA
anova(lme)

for r = 1:nROI
    
    data = [s1_2(r,:)'; s2_2(r,:)'; s3_2(r,:)'];
    
    group = [repmat({'G1'},5,1);
             repmat({'G2'},5,1);
             repmat({'G3'},5,1)];
    
    [p(r),~,stats] = anova1(data,group,'off');
    
    if p(r) < 0.05
        c = multcompare(stats,'Display','off');
        posthoc{r} = c;
    end
end

% FDR correction across the 11 ROIs
pFDR = mafdr(p,'BHFDR',true);

% Names of fixed effects
names = lme.Coefficients.Name;

nCoef = numel(names);
nROI = 11;

p12 = nan(nROI,1);   % G1 vs G2
p13 = nan(nROI,1);   % G1 vs G3
p23 = nan(nROI,1);   % G2 vs G3

for r = 1:nROI

    L12 = zeros(1,nCoef);
    L13 = zeros(1,nCoef);
    L23 = zeros(1,nCoef);

    % ---- ROI 1 ----
    if r == 1

        % G1 vs G2
        L12(strcmp(names,'Group_G2')) = 1;

        % G1 vs G3
        L13(strcmp(names,'Group_G3')) = 1;

        % G2 vs G3
        L23(strcmp(names,'Group_G2')) = 1;
        L23(strcmp(names,'Group_G3')) = -1;

    % ---- ROI 2-11 ----
    else

        % interaction terms
        intG2 = sprintf('Group_G2:ROI_%d',r);
        intG3 = sprintf('Group_G3:ROI_%d',r);
        roiName = sprintf('ROI_%d',r);

        % G1 vs G2
        L12(strcmp(names,'Group_G2')) = 1;
        L12(strcmp(names,intG2)) = 1;

        % G1 vs G3
        L13(strcmp(names,'Group_G3')) = 1;
        L13(strcmp(names,intG3)) = 1;

        % G2 vs G3
        L23(strcmp(names,'Group_G2')) = 1;
        L23(strcmp(names,intG2)) = 1;
        L23(strcmp(names,'Group_G3')) = -1;
        L23(strcmp(names,intG3)) = -1;

    end

    % Test contrasts
    p12(r) = coefTest(lme,L12);
    p13(r) = coefTest(lme,L13);
    p23(r) = coefTest(lme,L23);

end

q12 = mafdr(p12,'BHFDR',true);
q13 = mafdr(p13,'BHFDR',true);
q23 = mafdr(p23,'BHFDR',true);
sig12 = q12 < 0.05
sig13 = q13 < 0.05
sig23 = q23 < 0.05


