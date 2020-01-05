%Make Figure 2 (global plots)

%---
%Fig. 2a - human-driven changes over 1900-2015 from LUH2

nyout=13;

%Load data

luh2woodharv_year=NaN(nyout,1);
luh2woodharv_baseline=NaN(nyout,15);
luh2woodharv_baseline_old=NaN(nyout,15);
luh2woodharv_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2woodharv_unmasked_sens_global.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    %Get the year
    data=textscan(fid,'%d\n',1);
    luh2woodharv_year(nn)=data{1};
    %Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    %Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    luh2woodharv_baseline(nn,:)=cell2mat(data);
    %Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    luh2woodharv_baseline_old(nn,:)=cell2mat(data);
    %Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    luh2woodharv_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

luh2_year=NaN(nyout,1);
luh2_baseline=NaN(nyout,15);
luh2_baseline_old=NaN(nyout,15);
luh2_baseline_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2_unmasked_sens_global.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    %Get the year
    data=textscan(fid,'%d\n',1);
    luh2_year(nn)=data{1};
    %Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    %Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    luh2_baseline(nn,:)=cell2mat(data);
    %Get the baseline old-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    luh2_baseline_old(nn,:)=cell2mat(data);
    %Get the baseline young-skew simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    luh2_baseline_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

if length(luh2woodharv_year)~=length(luh2_year)
    error('Different number of years in each input file!')
end

luh2woodharv_baseline_regrowth=sum(luh2woodharv_baseline(:,1:14),2);
luh2woodharv_baseline_old_regrowth=sum(luh2woodharv_baseline_old(:,1:14),2);
luh2woodharv_baseline_young_regrowth=sum(luh2woodharv_baseline_young(:,1:14),2);
luh2woodharv_baseline_oldgrowth=luh2woodharv_baseline(:,15);
luh2woodharv_baseline_old_oldgrowth=luh2woodharv_baseline_old(:,15);
luh2woodharv_baseline_young_oldgrowth=luh2woodharv_baseline_young(:,15);

luh2_baseline_regrowth=sum(luh2_baseline(:,1:14),2);
luh2_baseline_old_regrowth=sum(luh2_baseline_old(:,1:14),2);
luh2_baseline_young_regrowth=sum(luh2_baseline_young(:,1:14),2);
luh2_baseline_oldgrowth=luh2_baseline(:,15);
luh2_baseline_old_oldgrowth=luh2_baseline_old(:,15);
luh2_baseline_young_oldgrowth=luh2_baseline_young(:,15);

%Make the plot
figure
hold on
e1=errorbar(luh2woodharv_year,luh2woodharv_baseline_regrowth,abs(luh2woodharv_baseline_regrowth-luh2woodharv_baseline_old_regrowth),...
    abs(luh2woodharv_baseline_regrowth-luh2woodharv_baseline_young_regrowth),'color','b');
e2=errorbar(luh2woodharv_year,luh2woodharv_baseline_oldgrowth,abs(luh2woodharv_baseline_oldgrowth-luh2woodharv_baseline_old_oldgrowth),...
    abs(luh2woodharv_baseline_oldgrowth-luh2woodharv_baseline_young_oldgrowth),'color','k');
e3=errorbar(luh2_year,luh2_baseline_regrowth,abs(luh2_baseline_regrowth-luh2_baseline_old_regrowth),...
    abs(luh2_baseline_regrowth-luh2_baseline_young_regrowth),'color','b','linestyle','--');
e4=errorbar(luh2_year,luh2_baseline_oldgrowth,abs(luh2_baseline_oldgrowth-luh2_baseline_old_oldgrowth),...
    abs(luh2_baseline_oldgrowth-luh2_baseline_young_oldgrowth),'color','k','linestyle','--');
xlabel('Years')
ylabel('Forest area (million km^2)')
legend('<140 y; LUC+WH','\geq140 y; LUC+WH','<140 y; LUC','\geq140 y; LUC')

clear all

%---
%Fig. 2b - changes from all sources in closed-canopy forests over 1900-2100, including a disturbance scenario

nyout=22;

%Load data

luh2dist_year=NaN(nyout,1);
luh2dist_baseline=NaN(nyout,15);
%luh2dist_baseline_old=NaN(nyout,15);
%luh2dist_baseline_young=NaN(nyout,15);
luh2dist_sens=NaN(nyout,15);
%luh2dist_sens_old=NaN(nyout,15);
%luh2dist_sens_young=NaN(nyout,15);
fid=fopen('age_reconstruction_luh2dist_masked_fut_scen_global.csv');
dump=textscan(fid,'%s\n','delimiter',',','delimiter','\n');
for nn=1:nyout
    %Get the year
    data=textscan(fid,'%d\n',1);
    luh2dist_year(nn)=data{1};
    %Ignore the column headings
    dump=textscan(fid,'%s\n',1,'delimiter','\n');
    %Get the baseline simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    luh2dist_baseline(nn,:)=cell2mat(data);
    %Get the baseline old-skew simulation data
    %data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    %luh2dist_baseline_old(nn,:)=cell2mat(data);
    %Get the baseline young-skew simulation data
    %data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    %luh2dist_baseline_young(nn,:)=cell2mat(data);
    %Get the sensitivity simulation data
    data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    luh2dist_sens(nn,:)=cell2mat(data);
    %Get the sensitivity old-skew simulation data
    %data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    %luh2dist_sens_old(nn,:)=cell2mat(data);
    %Get the sensitivity young-skew simulation data
    %data=textscan(fid,'%*s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',1,'delimiter',',');
    %luh2dist_sens_young(nn,:)=cell2mat(data);
end
clear nn data dump
fclose(fid);

%Remove the 2015 data point
aa=find(luh2dist_year==2015);
luh2dist_year(aa)=[];
luh2dist_baseline(aa,:)=[];
%luh2dist_baseline_old(aa,:)=[];
%luh2dist_baseline_young(aa,:)=[];
luh2dist_sens(aa,:)=[];
%luh2dist_sens_old(aa,:)=[];
%luh2dist_sens_young(aa,:)=[];
clear aa

luh2dist_baseline_regrowth=sum(luh2dist_baseline(:,1:14),2);
%luh2dist_baseline_old_regrowth=sum(luh2dist_baseline_old(:,1:14),2);
%luh2dist_baseline_young_regrowth=sum(luh2dist_baseline_young(:,1:14),2);
luh2dist_baseline_old_regrowth=luh2dist_baseline_regrowth;
luh2dist_baseline_young_regrowth=luh2dist_baseline_regrowth;
luh2dist_baseline_oldgrowth=luh2dist_baseline(:,15);
%luh2dist_baseline_old_oldgrowth=luh2dist_baseline_old(:,15);
%luh2dist_baseline_young_oldgrowth=luh2dist_baseline_young(:,15);
luh2dist_baseline_old_oldgrowth=luh2dist_baseline_oldgrowth;
luh2dist_baseline_young_oldgrowth=luh2dist_baseline_oldgrowth;

luh2dist_sens_regrowth=sum(luh2dist_sens(:,1:14),2);
%luh2dist_sens_old_regrowth=sum(luh2dist_sens_old(:,1:14),2);
%luh2dist_sens_young_regrowth=sum(luh2dist_sens_young(:,1:14),2);
luh2dist_sens_old_regrowth=luh2dist_sens_regrowth;
luh2dist_sens_young_regrowth=luh2dist_sens_regrowth;
luh2dist_sens_oldgrowth=luh2dist_sens(:,15);
%luh2dist_sens_old_oldgrowth=luh2dist_sens_old(:,15);
%luh2dist_sens_young_oldgrowth=luh2dist_sens_young(:,15);
luh2dist_sens_old_oldgrowth=luh2dist_sens_oldgrowth;
luh2dist_sens_young_oldgrowth=luh2dist_sens_oldgrowth;

%Make the plot
figure
hold on
e1=errorbar(luh2dist_year,luh2dist_baseline_regrowth,abs(luh2dist_baseline_regrowth-luh2dist_baseline_old_regrowth),...
    abs(luh2dist_baseline_regrowth-luh2dist_baseline_young_regrowth),'color','b');
e2=errorbar(luh2dist_year,luh2dist_baseline_oldgrowth,abs(luh2dist_baseline_oldgrowth-luh2dist_baseline_old_oldgrowth),...
    abs(luh2dist_baseline_oldgrowth-luh2dist_baseline_young_oldgrowth),'color','k');
e3=errorbar(luh2dist_year,luh2dist_sens_regrowth,abs(luh2dist_sens_regrowth-luh2dist_sens_old_regrowth),...
    abs(luh2dist_sens_regrowth-luh2dist_sens_young_regrowth),'color','b','linestyle','--');
e4=errorbar(luh2dist_year,luh2dist_sens_oldgrowth,abs(luh2dist_sens_oldgrowth-luh2dist_sens_old_oldgrowth),...
    abs(luh2dist_sens_oldgrowth-luh2dist_sens_young_oldgrowth),'color','k','linestyle','--');
xlabel('Years')
ylabel('Closed-canopy forest area (million km^2)')
legend('<140 y; Baseline','\geq140 y; Baseline','<140 y; Inc. dist','\geq140 y; Inc. dist')

%---
%Fig. 2c - Same data as for Fig. 2b, but showing age distributions

ages=5:10:150;

ind1=find(luh2dist_year==2000);
ind2=find(luh2dist_year==2050);
ind3=find(luh2dist_year==2100);

figure
[p1 h1 h2]=plotyy(ages(1:14),luh2dist_baseline(ind1,1:14),ages(15),luh2dist_baseline(ind1,15));
hold(p1(1)); hold(p1(2));
set(h1,'marker','.','markersize',10,'color','k')
set(h2,'marker','.','markersize',15,'color','k','linestyle','none')
set(p1(1),'ycolor','k'); set(p1(2),'ycolor','k');
%Other years for baseline simulation
plot(p1(1),ages(1:14),luh2dist_baseline(ind2,1:14),'k-.','markersize',5,'marker','^')
plot(p1(1),ages(1:14),luh2dist_baseline(ind3,1:14),'k-.','markersize',5,'marker','square')
plot(p1(2),ages(15),luh2dist_baseline(ind2,15),'k.','markersize',5,'marker','^')
plot(p1(2),ages(15),luh2dist_baseline(ind3,15),'k.','markersize',5,'marker','square')
%All years for sensitivity simulation
plot(p1(1),ages(1:14),luh2dist_sens(ind2,1:14),'r-.','markersize',5,'marker','^')
plot(p1(1),ages(1:14),luh2dist_sens(ind3,1:14),'r--','markersize',5,'marker','square')
plot(p1(2),ages(15),luh2dist_sens(ind2,15),'r.','markersize',5,'marker','^')
plot(p1(2),ages(15),luh2dist_sens(ind3,15),'r.','markersize',5,'marker','square')

set(p1(1),'XLim',[0 155],'YLim',[0 Inf],'Box','off', 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
set(p1(2),'XLim',[0 155],'YLim',[0 Inf], 'YTickMode', 'auto', 'YTickLabelMode', 'auto')
set(p1(1),'XTick',10:10:150,'XTickLabel',{'1-10','11-20','21-30','31-40','41-50','51-60',...
    '61-70','71-80','81-90','91-100','101-110','111-120','121-130','131-140','OG'})
set(p1(1),'XTickLabelRotation',300)
ylabel(p1(1),'Closed-canopy young forest area (million km^{-2})')
ylabel(p1(2),'Closed-canopy OG forest area (million km^{-2})')
set(get(p1(2),'Ylabel'),'Rotation',270,'VerticalAlignment','bottom')
legend('Baseline; 2000','Baseline; 2050','Baseline; 2100','Inc. dist; 2050','Inc. dist; 2100')
