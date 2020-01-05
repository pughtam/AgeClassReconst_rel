% Make Figure 3 (regional plots of human-driven changes over 1900-2015 from LUH2)
%
% Reads csv data created using age_class_reconstruction.m
%
% T. Pugh
% 05.01.20

%---
%Fig. 2a - tropical forest

nyout=13;

%Load data - tropical forest
trop_luh2woodharv_year=NaN(nyout,1);
trop_luh2woodharv_baseline=NaN(nyout,15);
trop_luh2woodharv_baseline_old=NaN(nyout,15);
trop_luh2woodharv_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2woodharv_unmasked_sens_region_tropical.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    %Get the year
    data=textscan(fid,'%d\n',1);
    trop_luh2woodharv_year(nn)=data{1};
    %Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    %Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2woodharv_baseline(nn,:)=cell2mat(data);
    %Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2woodharv_baseline_old(nn,:)=cell2mat(data);
    %Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2woodharv_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

trop_luh2_year=NaN(nyout,1);
trop_luh2_baseline=NaN(nyout,15);
trop_luh2_baseline_old=NaN(nyout,15);
trop_luh2_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2_unmasked_sens_region_tropical.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    %Get the year
    data=textscan(fid,'%d\n',1);
    trop_luh2_year(nn)=data{1};
    %Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    %Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2_baseline(nn,:)=cell2mat(data);
    %Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2_baseline_old(nn,:)=cell2mat(data);
    %Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

%Load data - temperate and Mediterranean forests
temp_luh2woodharv_year=NaN(nyout,1);
temp_luh2woodharv_baseline=NaN(nyout,15);
temp_luh2woodharv_baseline_old=NaN(nyout,15);
temp_luh2woodharv_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2woodharv_unmasked_sens_region_temperate.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    %Get the year
    data=textscan(fid,'%d\n',1);
    temp_luh2woodharv_year(nn)=data{1};
    %Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    %Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2woodharv_baseline(nn,:)=cell2mat(data);
    %Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2woodharv_baseline_old(nn,:)=cell2mat(data);
    %Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2woodharv_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

temp_luh2_year=NaN(nyout,1);
temp_luh2_baseline=NaN(nyout,15);
temp_luh2_baseline_old=NaN(nyout,15);
temp_luh2_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2_unmasked_sens_region_temperate.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    %Get the year
    data=textscan(fid,'%d\n',1);
    temp_luh2_year(nn)=data{1};
    %Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    %Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2_baseline(nn,:)=cell2mat(data);
    %Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2_baseline_old(nn,:)=cell2mat(data);
    %Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

%Load data - boreal forest
bor_luh2woodharv_year=NaN(nyout,1);
bor_luh2woodharv_baseline=NaN(nyout,15);
bor_luh2woodharv_baseline_old=NaN(nyout,15);
bor_luh2woodharv_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2woodharv_unmasked_sens_region_boreal.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    %Get the year
    data=textscan(fid,'%d\n',1);
    bor_luh2woodharv_year(nn)=data{1};
    %Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    %Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2woodharv_baseline(nn,:)=cell2mat(data);
    %Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2woodharv_baseline_old(nn,:)=cell2mat(data);
    %Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2woodharv_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

bor_luh2_year=NaN(nyout,1);
bor_luh2_baseline=NaN(nyout,15);
bor_luh2_baseline_old=NaN(nyout,15);
bor_luh2_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2_unmasked_sens_region_boreal.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    %Get the year
    data=textscan(fid,'%d\n',1);
    bor_luh2_year(nn)=data{1};
    %Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    %Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2_baseline(nn,:)=cell2mat(data);
    %Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2_baseline_old(nn,:)=cell2mat(data);
    %Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

if length(trop_luh2woodharv_year)~=length(temp_luh2woodharv_year) || length(trop_luh2woodharv_year)~=length(bor_luh2woodharv_year)
    error('1: Different number of years in input files!')
end
if length(trop_luh2_year)~=length(temp_luh2_year) || length(trop_luh2_year)~=length(bor_luh2_year)
    error('2: Different number of years in input files!')
end
if length(trop_luh2woodharv_year)~=length(trop_luh2_year)
    error('3: Different number of years in input files!')
end

trop_luh2woodharv_baseline_regrowth=sum(trop_luh2woodharv_baseline(:,1:14),2);
trop_luh2woodharv_baseline_old_regrowth=sum(trop_luh2woodharv_baseline_old(:,1:14),2);
trop_luh2woodharv_baseline_young_regrowth=sum(trop_luh2woodharv_baseline_young(:,1:14),2);
trop_luh2woodharv_baseline_oldgrowth=trop_luh2woodharv_baseline(:,15);
trop_luh2woodharv_baseline_old_oldgrowth=trop_luh2woodharv_baseline_old(:,15);
trop_luh2woodharv_baseline_young_oldgrowth=trop_luh2woodharv_baseline_young(:,15);

trop_luh2_baseline_regrowth=sum(trop_luh2_baseline(:,1:14),2);
trop_luh2_baseline_old_regrowth=sum(trop_luh2_baseline_old(:,1:14),2);
trop_luh2_baseline_young_regrowth=sum(trop_luh2_baseline_young(:,1:14),2);
trop_luh2_baseline_oldgrowth=trop_luh2_baseline(:,15);
trop_luh2_baseline_old_oldgrowth=trop_luh2_baseline_old(:,15);
trop_luh2_baseline_young_oldgrowth=trop_luh2_baseline_young(:,15);

temp_luh2woodharv_baseline_regrowth=sum(temp_luh2woodharv_baseline(:,1:14),2);
temp_luh2woodharv_baseline_old_regrowth=sum(temp_luh2woodharv_baseline_old(:,1:14),2);
temp_luh2woodharv_baseline_young_regrowth=sum(temp_luh2woodharv_baseline_young(:,1:14),2);
temp_luh2woodharv_baseline_oldgrowth=temp_luh2woodharv_baseline(:,15);
temp_luh2woodharv_baseline_old_oldgrowth=temp_luh2woodharv_baseline_old(:,15);
temp_luh2woodharv_baseline_young_oldgrowth=temp_luh2woodharv_baseline_young(:,15);

temp_luh2_baseline_regrowth=sum(temp_luh2_baseline(:,1:14),2);
temp_luh2_baseline_old_regrowth=sum(temp_luh2_baseline_old(:,1:14),2);
temp_luh2_baseline_young_regrowth=sum(temp_luh2_baseline_young(:,1:14),2);
temp_luh2_baseline_oldgrowth=temp_luh2_baseline(:,15);
temp_luh2_baseline_old_oldgrowth=temp_luh2_baseline_old(:,15);
temp_luh2_baseline_young_oldgrowth=temp_luh2_baseline_young(:,15);

bor_luh2woodharv_baseline_regrowth=sum(bor_luh2woodharv_baseline(:,1:14),2);
bor_luh2woodharv_baseline_old_regrowth=sum(bor_luh2woodharv_baseline_old(:,1:14),2);
bor_luh2woodharv_baseline_young_regrowth=sum(bor_luh2woodharv_baseline_young(:,1:14),2);
bor_luh2woodharv_baseline_oldgrowth=bor_luh2woodharv_baseline(:,15);
bor_luh2woodharv_baseline_old_oldgrowth=bor_luh2woodharv_baseline_old(:,15);
bor_luh2woodharv_baseline_young_oldgrowth=bor_luh2woodharv_baseline_young(:,15);

bor_luh2_baseline_regrowth=sum(bor_luh2_baseline(:,1:14),2);
bor_luh2_baseline_old_regrowth=sum(bor_luh2_baseline_old(:,1:14),2);
bor_luh2_baseline_young_regrowth=sum(bor_luh2_baseline_young(:,1:14),2);
bor_luh2_baseline_oldgrowth=bor_luh2_baseline(:,15);
bor_luh2_baseline_old_oldgrowth=bor_luh2_baseline_old(:,15);
bor_luh2_baseline_young_oldgrowth=bor_luh2_baseline_young(:,15);

%Make the plot
figure
s1=subplot(3,1,1);
hold on
e1=errorbar(trop_luh2woodharv_year,trop_luh2woodharv_baseline_regrowth,abs(trop_luh2woodharv_baseline_regrowth-trop_luh2woodharv_baseline_old_regrowth),...
    abs(trop_luh2woodharv_baseline_regrowth-trop_luh2woodharv_baseline_young_regrowth),'color','b');
e2=errorbar(trop_luh2woodharv_year,trop_luh2woodharv_baseline_oldgrowth,abs(trop_luh2woodharv_baseline_oldgrowth-trop_luh2woodharv_baseline_old_oldgrowth),...
    abs(trop_luh2woodharv_baseline_oldgrowth-trop_luh2woodharv_baseline_young_oldgrowth),'color','k');
ylabel('Forest area (million km^2)')
set(gca,'XTickLabel','')
title('(a) Tropical')

s2=subplot(3,1,2);
hold on
e3=errorbar(temp_luh2woodharv_year,temp_luh2woodharv_baseline_regrowth,abs(temp_luh2woodharv_baseline_regrowth-temp_luh2woodharv_baseline_old_regrowth),...
    abs(temp_luh2woodharv_baseline_regrowth-temp_luh2woodharv_baseline_young_regrowth),'color','b');
e4=errorbar(temp_luh2woodharv_year,temp_luh2woodharv_baseline_oldgrowth,abs(temp_luh2woodharv_baseline_oldgrowth-temp_luh2woodharv_baseline_old_oldgrowth),...
    abs(temp_luh2woodharv_baseline_oldgrowth-temp_luh2woodharv_baseline_young_oldgrowth),'color','k');
ylabel('Forest area (million km^2)')
set(gca,'XTickLabel','')
title('(b) Temperate and Mediterranean')

s3=subplot(3,1,3);
hold on
e5=errorbar(bor_luh2woodharv_year,bor_luh2woodharv_baseline_regrowth,abs(bor_luh2woodharv_baseline_regrowth-bor_luh2woodharv_baseline_old_regrowth),...
    abs(bor_luh2woodharv_baseline_regrowth-bor_luh2woodharv_baseline_young_regrowth),'color','b');
e6=errorbar(bor_luh2woodharv_year,bor_luh2woodharv_baseline_oldgrowth,abs(bor_luh2woodharv_baseline_oldgrowth-bor_luh2woodharv_baseline_old_oldgrowth),...
    abs(bor_luh2woodharv_baseline_oldgrowth-bor_luh2woodharv_baseline_young_oldgrowth),'color','k');
xlabel('Years')
ylabel('Forest area (million km^2)')
title('(c) Boreal')

%Add on the lines for LUC only
e7=errorbar(s1,trop_luh2_year,trop_luh2_baseline_regrowth,abs(trop_luh2_baseline_regrowth-trop_luh2_baseline_old_regrowth),...
    abs(trop_luh2_baseline_regrowth-trop_luh2_baseline_young_regrowth),'color','b','linestyle',':');
e8=errorbar(s1,trop_luh2_year,trop_luh2_baseline_oldgrowth,abs(trop_luh2_baseline_oldgrowth-trop_luh2_baseline_old_oldgrowth),...
    abs(trop_luh2_baseline_oldgrowth-trop_luh2_baseline_young_oldgrowth),'color','k','linestyle',':');

e9=errorbar(s2,temp_luh2_year,temp_luh2_baseline_regrowth,abs(temp_luh2_baseline_regrowth-temp_luh2_baseline_old_regrowth),...
    abs(temp_luh2_baseline_regrowth-temp_luh2_baseline_young_regrowth),'color','b','linestyle',':');
e10=errorbar(s2,temp_luh2_year,temp_luh2_baseline_oldgrowth,abs(temp_luh2_baseline_oldgrowth-temp_luh2_baseline_old_oldgrowth),...
    abs(temp_luh2_baseline_oldgrowth-temp_luh2_baseline_young_oldgrowth),'color','k','linestyle',':');

e11=errorbar(s3,bor_luh2_year,bor_luh2_baseline_regrowth,abs(bor_luh2_baseline_regrowth-bor_luh2_baseline_old_regrowth),...
    abs(bor_luh2_baseline_regrowth-bor_luh2_baseline_young_regrowth),'color','b','linestyle',':');
e12=errorbar(s3,bor_luh2_year,bor_luh2_baseline_oldgrowth,abs(bor_luh2_baseline_oldgrowth-bor_luh2_baseline_old_oldgrowth),...
    abs(bor_luh2_baseline_oldgrowth-bor_luh2_baseline_young_oldgrowth),'color','k','linestyle',':');

legend(s1,'<140 y; LUC+WH','\geq140 y; LUC+WH','<140 y; LUC','\geq140 y; LUC')

set(s1,'Position',[0.13 0.68 0.775 0.26])
set(s2,'Position',[0.13 0.38 0.775 0.26])
set(s3,'Position',[0.13 0.08 0.775 0.26])
