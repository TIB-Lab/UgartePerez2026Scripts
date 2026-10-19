clear all

close all

stat_dir=' '; % mask or stat

data_dir=' '; % FC maps

file_list=importdata(' '); % list of file names

n = length(file_list);

Animal = cell(n,1);
ScanCode = cell(n,1);

for i = 1:n
    parts = strsplit(strtrim(file_list{i}));
    
    Animal{i}   = parts{1}; % AM3
    ScanCode{i} = parts{2}; % A050521
end

tstat=load_nii([stat_dir 'IC1_DMN_mask.nii']);tstat_p=tstat.img; % anterior

tstat_p(:,1:26,:)=0;

tstat.img=tstat_p;
save_nii(tstat,'fMRI_ant_DMN.nii')

tstat=load_nii([stat_dir 'IC1_DMN_mask.nii']);tstat_n=tstat.img; % posterior
tstat_n(:,16:end,:)=0;

tstat.img=tstat_n;
save_nii(tstat,'fMRI_post_DMN.nii')

d1=dir(data_dir);

count=0;

sex_364=[1 2 1 2 2 1 1 1 2 2 2]; % 1 male 2 female

sex_200=[1 1 1 1 1 1 2 2 2 2 2 2 2 1]; % 1 male 2 female

sex_546=[2 2 2 1 1 1 2];


for i=1:length(d1)

    if strfind(d1(i,:).name,'546_')

        count=count+1;

        V=load_nii([data_dir d1(i,:).name]);V=V.img;

        val_546_p(count)=mean(V(tstat_p>0.95));

        val_546_n(count)=mean(V(tstat_n>0.95));

        idx = find(strcmp(ScanCode, d1(i).name(5:11)));

        if isempty(idx)
            continue
        end

        animal_546{count} = Animal{idx,1};

    end

end

count=0;

for i=1:length(d1)

    if strfind(d1(i,:).name,'36')

        count=count+1;

        V=load_nii([data_dir d1(i,:).name]);V=V.img;

        val_364_p(count)=mean(V(tstat_p>0.95));

        val_364_n(count)=mean(V(tstat_n>0.95));

        idx = find(strcmp(ScanCode, d1(i).name(5:11)));

        if isempty(idx)
            continue
        end

        animal_364{count} = Animal{idx,1};

    end

end


count=0;

for i=1:length(d1)

    if strfind(d1(i,:).name,'200_')

        count=count+1;

        V=load_nii([data_dir d1(i,:).name]);V=V.img;

        val_200_p(count)=mean(V(tstat_p>0.95));

        val_200_n(count)=mean(V(tstat_n>0.95));

        idx = find(strcmp(ScanCode, d1(i).name(5:11)));

        if isempty(idx)
            continue
        end

        animal_200{count} = Animal{idx,1};

    end

end



figure(1);
%male blue
bar([1 2 3 4]-0.1,[mean(val_200_p(sex_200==1)) mean(val_546_p(sex_546==1)) mean(val_200_n(sex_200==1)) mean(val_546_n(sex_546==1))],0.2)
hold on
plot([ones(size(val_200_p(sex_200==1))) 2*ones(size(val_546_p(sex_546==1))) 3*ones(size(val_200_n(sex_200==1))) 4*ones(size(val_546_n(sex_546==1)))]-0.1,[val_200_p(sex_200==1) val_546_p(sex_546==1) val_200_n(sex_200==1) val_546_n(sex_546==1)],'ob')

hold on
%female yellow, ANT 200 500 POST 200 500
bar([1 2 3 4]+0.1,[mean(val_200_p(sex_200==2)) mean(val_546_p(sex_546==2)) mean(val_200_n(sex_200==2)) mean(val_546_n(sex_546==2))],0.2)

plot([ones(size(val_200_p(sex_200==2))) 2*ones(size(val_546_p(sex_546==2))) 3*ones(size(val_200_n(sex_200==2))) 4*ones(size(val_546_n(sex_546==2)))]+0.1,[val_200_p(sex_200==2) val_546_p(sex_546==2) val_200_n(sex_200==2) val_546_n(sex_546==2)],'or')


xlim([0 5])


% linear model P - ANTERIOR

Y = [val_200_p(:); val_546_p(:)];

Age = [repmat(200, length(val_200_p),1); ...
       repmat(546, length(val_546_p),1)];

Sex = [sex_200(:); sex_546(:)];

Animal = [animal_200(:); animal_546(:)];

tbl = table(Y, categorical(Age), categorical(Sex), categorical(Animal), ...
    'VariableNames', {'Y','Age','Sex','Animal'});

lme = fitlme(tbl, 'Y ~ Age * Sex + (1|Animal)','DummyVarCoding', 'effects');
anova(lme)

% males
tbl_m = tbl(tbl.Sex==categorical(1),:);
lme_m = fitlme(tbl_m, 'Y ~ Age + (1|Animal)','DummyVarCoding', 'effects');

statm=anova(lme_m);

% females
tbl_f = tbl(tbl.Sex==categorical(2),:);
lme_f = fitlme(tbl_f, 'Y ~ Age + (1|Animal)','DummyVarCoding', 'effects');

statf=anova(lme_f);

[p_corr, h] = bonferroni_holm_corrections([statm.pValue(2) statf.pValue(2)])

[statm.FStat(2) statf.FStat(2)]

% linear model N - posterior

Y = [val_200_n(:); val_546_n(:)];

Age = [repmat(200, length(val_200_n),1); ...
       repmat(546, length(val_546_n),1)];

Sex = [sex_200(:); sex_546(:)];

Animal = [animal_200(:); animal_546(:)];

tbl = table(Y, categorical(Age), categorical(Sex), categorical(Animal), ...
    'VariableNames', {'Y','Age','Sex','Animal'});

lme = fitlme(tbl, 'Y ~ Age * Sex + (1|Animal)','DummyVarCoding', 'effects');
anova(lme)

% males
tbl_m = tbl(tbl.Sex==categorical(1),:);
lme_m = fitlme(tbl_m, 'Y ~ Age + (1|Animal)','DummyVarCoding', 'effects');

statm=anova(lme_m);

% females
tbl_f = tbl(tbl.Sex==categorical(2),:);
lme_f = fitlme(tbl_f, 'Y ~ Age + (1|Animal)','DummyVarCoding', 'effects');

statf=anova(lme_f);

[p_corr, h] = bonferroni_holm_corrections([statm.pValue(2) statf.pValue(2)])

[statm.FStat(2) statf.FStat(2)]
