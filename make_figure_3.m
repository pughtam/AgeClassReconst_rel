% Make Figure 3 (regional plots of human-driven changes over 1900-2015 from LUH2)
%
% Reads csv data created using age_class_reconstruction.m
%
% Dependencies:
% - shadedErrorBar (Rob Campbell (2020). raacampbell/shadedErrorBar (https://www.github.com/raacampbell/shadedErrorBar),
% GitHub. Retrieved January 5, 2020.)
%
% T. Pugh
% 05.01.20

%--- Fig. 2a - tropical forest ---

nyout=13;

% Load data - tropical forest
trop_luh2woodharv_year=NaN(nyout,1);
trop_luh2woodharv_baseline=NaN(nyout,15);
trop_luh2woodharv_baseline_old=NaN(nyout,15);
trop_luh2woodharv_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2woodharv_unmasked_sens_region_tropical.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    % Get the year
    data=textscan(fid,'%d\n',1);
    trop_luh2woodharv_year(nn)=data{1};
    % Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    % Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2woodharv_baseline(nn,:)=cell2mat(data);
    % Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2woodharv_baseline_old(nn,:)=cell2mat(data);
    % Get the baseline young-skew simulation data
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
    % Get the year
    data=textscan(fid,'%d\n',1);
    trop_luh2_year(nn)=data{1};
    % Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    % Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2_baseline(nn,:)=cell2mat(data);
    % Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2_baseline_old(nn,:)=cell2mat(data);
    % Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    trop_luh2_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

% Load data - temperate and Mediterranean forests
temp_luh2woodharv_year=NaN(nyout,1);
temp_luh2woodharv_baseline=NaN(nyout,15);
temp_luh2woodharv_baseline_old=NaN(nyout,15);
temp_luh2woodharv_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2woodharv_unmasked_sens_region_temperate.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    % Get the year
    data=textscan(fid,'%d\n',1);
    temp_luh2woodharv_year(nn)=data{1};
    % Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    % Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2woodharv_baseline(nn,:)=cell2mat(data);
    % Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2woodharv_baseline_old(nn,:)=cell2mat(data);
    % Get the baseline young-skew simulation data
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
    % Get the year
    data=textscan(fid,'%d\n',1);
    temp_luh2_year(nn)=data{1};
    % Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    % Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2_baseline(nn,:)=cell2mat(data);
    % Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2_baseline_old(nn,:)=cell2mat(data);
    % Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    temp_luh2_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

% Load data - boreal forest
bor_luh2woodharv_year=NaN(nyout,1);
bor_luh2woodharv_baseline=NaN(nyout,15);
bor_luh2woodharv_baseline_old=NaN(nyout,15);
bor_luh2woodharv_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2woodharv_unmasked_sens_region_boreal.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    % Get the year
    data=textscan(fid,'%d\n',1);
    bor_luh2woodharv_year(nn)=data{1};
    % Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    % Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2woodharv_baseline(nn,:)=cell2mat(data);
    % Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2woodharv_baseline_old(nn,:)=cell2mat(data);
    % Get the baseline young-skew simulation data
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
    % Get the year
    data=textscan(fid,'%d\n',1);
    bor_luh2_year(nn)=data{1};
    % Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    % Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2_baseline(nn,:)=cell2mat(data);
    % Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    bor_luh2_baseline_old(nn,:)=cell2mat(data);
    % Get the baseline young-skew simulation data
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

% Make the plot
figure
s1=subplot(3,1,1);
hold on
e1=shadedErrorBar(trop_luh2woodharv_year,trop_luh2woodharv_baseline_regrowth,cat(2,abs(trop_luh2woodharv_baseline_regrowth-trop_luh2woodharv_baseline_old_regrowth),...
    abs(trop_luh2woodharv_baseline_regrowth-trop_luh2woodharv_baseline_young_regrowth)),'lineProps','b');
e2=shadedErrorBar(trop_luh2woodharv_year,trop_luh2woodharv_baseline_oldgrowth,cat(2,abs(trop_luh2woodharv_baseline_oldgrowth-trop_luh2woodharv_baseline_young_oldgrowth),...
    abs(trop_luh2woodharv_baseline_oldgrowth-trop_luh2woodharv_baseline_old_oldgrowth)),'lineProps','k');
ylabel('Forest area (million km^2)')
set(gca,'XTickLabel','')
title('(a) Tropical')
set(e1.mainLine,'linewidth',2)
set(e2.mainLine,'linewidth',2)
set(e1.edge,'linestyle','none')
set(e2.edge,'linestyle','none')

s2=subplot(3,1,2);
hold on
e3=shadedErrorBar(temp_luh2woodharv_year,temp_luh2woodharv_baseline_regrowth,cat(2,abs(temp_luh2woodharv_baseline_regrowth-temp_luh2woodharv_baseline_old_regrowth),...
    abs(temp_luh2woodharv_baseline_regrowth-temp_luh2woodharv_baseline_young_regrowth)),'lineProps','b');
e4=shadedErrorBar(temp_luh2woodharv_year,temp_luh2woodharv_baseline_oldgrowth,cat(2,abs(temp_luh2woodharv_baseline_oldgrowth-temp_luh2woodharv_baseline_young_oldgrowth),...
    abs(temp_luh2woodharv_baseline_oldgrowth-temp_luh2woodharv_baseline_old_oldgrowth)),'lineProps','k');
ylabel('Forest area (million km^2)')
set(gca,'XTickLabel','')
title('(b) Temperate and Mediterranean')
set(e3.mainLine,'linewidth',2)
set(e4.mainLine,'linewidth',2)
set(e3.edge,'linestyle','none')
set(e4.edge,'linestyle','none')

s3=subplot(3,1,3);
hold on
e5=shadedErrorBar(bor_luh2woodharv_year,bor_luh2woodharv_baseline_regrowth,cat(2,abs(bor_luh2woodharv_baseline_regrowth-bor_luh2woodharv_baseline_old_regrowth),...
    abs(bor_luh2woodharv_baseline_regrowth-bor_luh2woodharv_baseline_young_regrowth)),'lineProps','b');
e6=shadedErrorBar(bor_luh2woodharv_year,bor_luh2woodharv_baseline_oldgrowth,cat(2,abs(bor_luh2woodharv_baseline_oldgrowth-bor_luh2woodharv_baseline_young_oldgrowth),...
    abs(bor_luh2woodharv_baseline_oldgrowth-bor_luh2woodharv_baseline_old_oldgrowth)),'lineProps','k');
xlabel('Years')
ylabel('Forest area (million km^2)')
title('(c) Boreal')
set(e5.mainLine,'linewidth',2)
set(e6.mainLine,'linewidth',2)
set(e5.edge,'linestyle','none')
set(e6.edge,'linestyle','none')

% Add on the lines for land-use change only
subplot(3,1,1);
e7=shadedErrorBar(trop_luh2_year,trop_luh2_baseline_regrowth,cat(2,abs(trop_luh2_baseline_regrowth-trop_luh2_baseline_old_regrowth),...
    abs(trop_luh2_baseline_regrowth-trop_luh2_baseline_young_regrowth)),'lineProps','b--');
e8=shadedErrorBar(trop_luh2_year,trop_luh2_baseline_oldgrowth,cat(2,abs(trop_luh2_baseline_oldgrowth-trop_luh2_baseline_young_oldgrowth),...
    abs(trop_luh2_baseline_oldgrowth-trop_luh2_baseline_old_oldgrowth)),'lineProps','k--');
set(e7.mainLine,'linewidth',2)
set(e8.mainLine,'linewidth',2)
set(e7.edge,'linestyle','none')
set(e8.edge,'linestyle','none')
box on

subplot(3,1,2);
e9=shadedErrorBar(temp_luh2_year,temp_luh2_baseline_regrowth,cat(2,abs(temp_luh2_baseline_regrowth-temp_luh2_baseline_old_regrowth),...
    abs(temp_luh2_baseline_regrowth-temp_luh2_baseline_young_regrowth)),'lineProps','b--');
e10=shadedErrorBar(temp_luh2_year,temp_luh2_baseline_oldgrowth,cat(2,abs(temp_luh2_baseline_oldgrowth-temp_luh2_baseline_young_oldgrowth),...
    abs(temp_luh2_baseline_oldgrowth-temp_luh2_baseline_old_oldgrowth)),'lineProps','k--');
set(e9.mainLine,'linewidth',2)
set(e10.mainLine,'linewidth',2)
set(e9.edge,'linestyle','none')
set(e10.edge,'linestyle','none')
box on

subplot(3,1,3);
e11=shadedErrorBar(bor_luh2_year,bor_luh2_baseline_regrowth,cat(2,abs(bor_luh2_baseline_regrowth-bor_luh2_baseline_old_regrowth),...
    abs(bor_luh2_baseline_regrowth-bor_luh2_baseline_young_regrowth)),'lineProps','b--');
e12=shadedErrorBar(bor_luh2_year,bor_luh2_baseline_oldgrowth,cat(2,abs(bor_luh2_baseline_oldgrowth-bor_luh2_baseline_young_oldgrowth),...
    abs(bor_luh2_baseline_oldgrowth-bor_luh2_baseline_old_oldgrowth)),'lineProps','k--');
set(e11.mainLine,'linewidth',2)
set(e12.mainLine,'linewidth',2)
set(e11.edge,'linestyle','none')
set(e12.edge,'linestyle','none')
box on

legend(s1,'<140 y; LUC+WH','\geq140 y; LUC+WH','<140 y; LUC','\geq140 y; LUC')

set(gcf,'color','w')

set(s1,'XLim',[1900 2015])
set(s2,'XLim',[1900 2015])
set(s3,'XLim',[1900 2015])

set(s1,'Position',[0.13 0.68 0.775 0.26])
set(s2,'Position',[0.13 0.38 0.775 0.26])
set(s3,'Position',[0.13 0.08 0.775 0.26])
