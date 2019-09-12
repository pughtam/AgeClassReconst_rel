function age_class_reconstruction_9region_gfad_plot_func(fage_out_decade_reg,...
    gfad_fage_area_sum_reg_stan,gfad_fage_area_sum_reg_lower,gfad_fage_area_sum_reg_upper,...
    ages,regions,nregion,plotyear,plot_gfad)
%Make a figure with regional age class distributions for both the
%simulations and GFAD.
%
%T. Pugh
%12.09.19

yy=plotyear;
figure
for nn=1:nregion
    ss(nn)=subplot(3,3,nn);
    hold on
    plot(ages(1:14),fage_out_decade_reg(nn,1:14,yy,1),'r.-','markersize',15)
    plot(ages(1:14),fage_out_decade_reg(nn,1:14,yy,2),'r-.')
    plot(ages(1:14),fage_out_decade_reg(nn,1:14,yy,3),'r--')
    plot(ages(15),fage_out_decade_reg(nn,15,yy,1),'r.','markersize',15)
    plot(ages(15),fage_out_decade_reg(nn,15,yy,2),'r.','markersize',10,'marker','^')
    plot(ages(15),fage_out_decade_reg(nn,15,yy,3),'r.','markersize',10,'marker','v')
    if plot_gfad
        plot(ages(1:14),gfad_fage_area_sum_reg_upper(nn,1:14),'k-','markersize',15)
        plot(ages(1:14),gfad_fage_area_sum_reg_stan(nn,1:14),'-','markersize',15,'color',[0.7 0.7 0.7])
        plot(ages(1:14),gfad_fage_area_sum_reg_lower(nn,1:14),'-','markersize',15,'color',[0.7 0.7 0.7])
        plot(ages(15),gfad_fage_area_sum_reg_upper(nn,15),'k.','markersize',15)
        plot(ages(15),gfad_fage_area_sum_reg_stan(nn,15),'.','markersize',15,'color',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7])
        plot(ages(15),gfad_fage_area_sum_reg_lower(nn,15),'.','markersize',15,'color',[0.7 0.7 0.7],'markerfacecolor',[0.7 0.7 0.7])
    end
    title(regions{nn})
    set(gca,'XLim',[0 155])
    if nn==7 || nn==8 || nn==9
        set(gca,'XTick',10:10:150,'XTickLabel',{'1-10','11-20','21-30','31-40','41-50','51-60',...
            '61-70','71-80','81-90','91-100','101-110','111-120','121-130','131-140','OG'})
        set(gca,'XTickLabelRotation',300)
    end
    if nn<=6
        set(gca,'XTick','','XTickLabel','')
    end
    if nn==1 || nn==4 || nn==7
        ylabel('million km^{-2}')
    end
end
clear nn
% if gfad_comp
%     legend('LUH2+dist','GFAD upper','GFAD','GFAD lower')
% end
set(ss(1),'Position',[0.1 0.7 0.25 0.25])
set(ss(2),'Position',[0.4 0.7 0.25 0.25])
set(ss(3),'Position',[0.7 0.7 0.25 0.25])
set(ss(4),'Position',[0.1 0.4 0.25 0.25])
set(ss(5),'Position',[0.4 0.4 0.25 0.25])
set(ss(6),'Position',[0.7 0.4 0.25 0.25])
set(ss(7),'Position',[0.1 0.1 0.25 0.25])
set(ss(8),'Position',[0.4 0.1 0.25 0.25])
set(ss(9),'Position',[0.7 0.1 0.25 0.25])