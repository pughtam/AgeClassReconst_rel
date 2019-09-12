function age_class_reconstruction_regionalsel_plot_func(fage_out_decade_globe,fage_out_decade_reg,regsel,...
    ages,yname,nyout,regions)
%Make a figure with the global age class distribution in the main plot and 3 regional
%subplots.
%
%T. Pugh
%12.09.19

xticks={'1-10','11-20','21-30','31-40','41-50','51-60',...
    '61-70','71-80','81-90','91-100','101-110','111-120','121-130','131-140','OG'};

ndata=length(fage_out_decade_globe);

figure
ycols={'k','b','r'};
s1=subplot(2,2,1);
hold on
for yy=1:nyout
    plot(ages,fage_out_decade_globe(:,yy),'.-','markersize',15,'color',ycols{yy})
end
legend('1900','1950','2015')
ylabel(yname)
set(gca,'XTick',10:10:140,'XTickLabel',xticks(1:ndata))
set(gca,'XTickLabelRotation',300)
title('Global')
set(gca,'XLim',[0 ndata*10])

s2=subplot(2,2,2);
hold on
for yy=1:nyout
    plot(ages,fage_out_decade_reg(regsel(1),:,yy),'.-','markersize',15,'color',ycols{yy})
end
ylabel(yname)
set(gca,'XTick',10:10:140,'XTickLabel',xticks(1:ndata))
set(gca,'XTickLabelRotation',300)
title(regions{regsel(1)})
set(gca,'XLim',[0 ndata*10])

s3=subplot(2,2,3);
hold on
for yy=1:nyout
    plot(ages,fage_out_decade_reg(regsel(2),:,yy),'.-','markersize',15,'color',ycols{yy})
end
set(gca,'XTick',10:10:140,'XTickLabel',xticks(1:ndata))
set(gca,'XTickLabelRotation',300)
title(regions{regsel(2)})
set(gca,'XLim',[0 ndata*10])
xlabel('Age class')

s4=subplot(2,2,4);
hold on
for yy=1:nyout
    plot(ages,fage_out_decade_reg(regsel(3),:,yy),'.-','markersize',15,'color',ycols{yy})
end
set(gca,'XTick',10:10:140,'XTickLabel',xticks(1:ndata))
set(gca,'XTickLabelRotation',300)
title(regions{regsel(3)})
set(gca,'XLim',[0 ndata*10])

set(s1,'Position',[0.1 0.6 0.85 0.3])
set(s2,'Position',[0.1 0.15 0.25 0.3])
set(s3,'Position',[0.4 0.15 0.25 0.3])
set(s4,'Position',[0.7 0.15 0.25 0.3])