%Script to to construct the forest age distribution based on land-use change from LUH2 and on random
%disturbances as defined by Pugh et al. (2019).
%
%Note that this scripts needs upwards of 20Gb RAM to run for current settings
%
%T. Pugh
%10.08.19

%---
%Input options
year1=1900; %First calendar year of simulation for which output is required (i.e. after spin-up)
nyear=115; %Total number of years of simulation for which output is required
inc_woodharv=false; %Include the wood harvest transitions? Standard assumption is inc_woodharv=false
inc_dist=true; %Include background disturbance (true) or just LUH2 transitions (false)

use_dist_scen=true; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=1; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=2; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

output_crosscheck_plots=0; %Make diagnostic cross-check plots
output_years=[1900 1950 2015]; %Years for which to provide outputs
gfad_comp=false; %Include GFAD in the output plots

%---
%Settings
nages=150; %Number of age classes (1 class = 1 year)

%LUH2 input
luh2_filepath_trans='/data/LUH2/transitions.nc';
luh2_filepath_states='/data/LUH2/states.nc';
firstluh2year=850; %First calendar year in LUH2 dataset

%Background disturbance return time input
distfile='/home/adf/pughtam/data/Disturbance/netcdfs_for_deposition/tauO/tauO_standard_forest-area_LUcorrected.nc';

%GFAD data (only used for comparison in output plots
gfad_filepath_stan='/data/Poulter_forest_age/global_forestAgeClasses_2011b_CORRECTED.nc';
gfad_filepath_lower='/data/Poulter_forest_age/Uncertainties/global_forestAgeClasses_2011_minus40perc_tropmin_corrected.nc';
gfad_filepath_upper='/data/Poulter_forest_age/Uncertainties/global_forestAgeClasses_2011_plus40perc_tropmax_corrected.nc';

%---
%Read in necessary data

%Read initial forest state
nspinup=nages; %Start nages years before first output year to allow for spin-up of age classes
year1ind=year1-nspinup-firstluh2year+1;
[primf_init,secdf_init]=luh2_forstates_read(luh2_filepath_states,year1ind);

%Read forest loss and gain transitions from LUH2
nyearluh2=nyear+nspinup;
[luh2_forlu_loss_prim_1deg,luh2_forlu_loss_sec_1deg]=luh2_forloss_read(luh2_filepath_trans,year1ind,nyearluh2,inc_woodharv);
luh2_forlu_gain_1deg=luh2_forgain_read(luh2_filepath_trans,year1ind,nyearluh2,inc_woodharv);

%Read in disturbance return period
distint=ncread(distfile,'tauO');
distint=fliplr(distint);

%---
%Calculations section
%Track primary and secondary forest from LUH2 separately and only merge together in the final output array.

if use_dist_scen
    %If using a disturbance scenario, initialise the annual multiplicative array here
    dist_scen=NaN(nyear,1);
    for nn=1:nyear
        df=dist_scen_end-dist_scen_start;
        dist_scen(nn)=dist_scen_start+(nn*(df/nyear));
    end
    clear nn
end

%Initialise forest fractions as the maximum age
fage_sec=zeros(360,180,nages);
fage_sec(:,:,nages)=secdf_init;
fage_prim=zeros(360,180,nages);
fage_prim(:,:,nages)=primf_init;

%Now start the loop for forest age calculation
fage=zeros(360,180,nages,nyear+nspinup);
for yy=1:nyear+nspinup
    fprintf('%d\n',yy)

    %First increment the ages of forest
    %Secondary forest
    fage_sec(:,:,nages)=fage_sec(:,:,nages)+fage_sec(:,:,nages-1);
    fage_sec(:,:,2:nages-1)=fage_sec(:,:,1:nages-2);
    fage_sec(:,:,1)=zeros(360,180);
    %Primary forest
    fage_prim(:,:,nages)=fage_prim(:,:,nages)+fage_prim(:,:,nages-1);
    fage_prim(:,:,2:nages-1)=fage_prim(:,:,1:nages-2);
    fage_prim(:,:,1)=zeros(360,180);

    %Subtract forest loss
    %Remove losses from random age class until all losses are allocated
    %Primary forest
    for ii=1:360
        for jj=1:180
            if luh2_forlu_loss_prim_1deg(ii,jj,yy)>0
                to_lose=luh2_forlu_loss_prim_1deg(ii,jj,yy); %Running total of lost forest fraction still to be allocated
                while to_lose>0.00000001
                    hf=find(fage_prim(ii,jj,:)>0);
                    if ~isempty(hf)
                        rr=randi(length(hf)); %Randomly choose age class
                        if fage_prim(ii,jj,hf(rr)) > to_lose
                            fage_prim(ii,jj,hf(rr))=fage_prim(ii,jj,hf(rr))-to_lose;
                            break
                        else
                            to_lose=to_lose-fage_prim(ii,jj,hf(rr));
                            fage_prim(ii,jj,hf(rr))=0;
                        end
                    else
                        %If there is no primary forest to lose then ignore this loss
                        break
                    end
                end
                clear to_lose rr
            end
        end
        clear jj
    end
    clear ii
    %Secondary forest
    for ii=1:360
        for jj=1:180
            if luh2_forlu_loss_sec_1deg(ii,jj,yy)>0
                to_lose=luh2_forlu_loss_sec_1deg(ii,jj,yy); %Running total of lost forest fraction still to be allocated
                while to_lose>0.00000001
                    hf=find(fage_sec(ii,jj,:)>0);
                    if ~isempty(hf)
                        rr=randi(length(hf)); %Randomly choose age class
                        if fage_sec(ii,jj,hf(rr)) > to_lose
                            fage_sec(ii,jj,hf(rr))=fage_sec(ii,jj,hf(rr))-to_lose;
                            break
                        else
                            to_lose=to_lose-fage_sec(ii,jj,hf(rr));
                            fage_sec(ii,jj,hf(rr))=0;
                        end
                    else
                        %If there is no secondary forest to lose then ignore this loss
                        break
                    end
                end
                clear to_lose rr
            end
        end
        clear jj
    end
    clear ii

    %Add forest gain to youngest secondary age class
    fage_sec(:,:,1)=fage_sec(:,:,1)+luh2_forlu_gain_1deg(:,:,yy);

    if inc_dist
        %Carry out random disturbance (equal probability across all ages of secondary and primary forest) and add
        %to youngest age class of secondary forest.
        %Disturb a fixed fraction per year defined by 1/distint

        if use_dist_scen
            %Modify the disturbance rate
            if yy<=nspinup
                distrate=(1./distint)*dist_scen(1);
            else
                distrate=(1./distint)*dist_scen(yy-nspinup);
            end
        else
            distrate=(1./distint);
        end

        %Secondary forest
        frac_dist_sec=fage_sec(:,:,2:nages).*repmat(distrate,[1 1 nages-1]);
        fage_sec(:,:,2:nages)=fage_sec(:,:,2:nages)-frac_dist_sec;
        fage_sec(:,:,1)=fage_sec(:,:,1)+sum(frac_dist_sec,3);
        clear frac_dist_sec
        %Primary forest
        frac_dist_prim=fage_prim(:,:,2:nages).*repmat(distrate,[1 1 nages-1]);
        fage_prim(:,:,2:nages)=fage_prim(:,:,2:nages)-frac_dist_prim;
        fage_prim(:,:,1)=fage_prim(:,:,1)+sum(frac_dist_prim,3);
        clear frac_dist_prim
    end

    %Sum both primary and secondary forest in output array
    fage(:,:,:,yy)=fage_sec+fage_prim;
end

if output_crosscheck_plots
    %Cross-check total forest fraction with states file for end of the calculation
    %NOTE: These should only agree if wood harvest is included in the transitions (i.e. inc_woodharv=1)
    [primf_end,secdf_end]=luh2_forstates_read(luh2_filepath_states,year1ind+nspinup+nyear);
    res_check=sum(fage(:,:,:,nspinup+nyear),3)-primf_end-secdf_end;
    
    figure
    p1=pcolor(flipud(res_check'));
    set(p1,'linestyle','none')
    colorbar
end

%---
%Postprocessing of fage

%Select desired years
nyout=length(output_years);
output_years_ind=output_years-year1+nspinup;
fage_out=NaN(360,180,nages,nyout);
for nn=1:nyout
    fage_out(:,:,:,nn)=fage(:,:,:,output_years_ind(nn));
end
clear nn output_years_ind

%Group to decadal timesteps
nages_dec=nages/10;
fage_out_decade=NaN(360,180,nages_dec,nyout);
for nn=1:nyout
    dd=0;
    for aa=1:10:nages
        dd=dd+1;
        fage_out_decade(:,:,dd,nn)=sum(fage_out(:,:,aa:aa+9,nn),3);
    end
    clear dd aa
end
clear nn

%---
%Plot age distributions for a variety of regions

%Read in forest area
fmask=fliplr(ncread('/data/Disturbance/netcdfs_for_deposition/forestmask/hansen_forested_frac_1deg.nc4','forested_50_percent'));
fmask=double(fmask)./100;

%Calculate grid-cell area
garea=global_grid_area_1deg()';

%Load region mask (currently ESA)
addpath('/data/ESA_landcover/')
[rmask,regions,nregion]=esa_forest_9regions_new_1deg_func(false);
rmask=fliplr(rmask');

%Convert fractions to areas
fage_out_decade_totfor=squeeze(sum(fage_out_decade,3)); %Total forest fraction
fage_out_decade_frac=NaN(size(fage_out_decade));
fage_out_decade_area=NaN(size(fage_out_decade));
for nn=1:nyout
    %Standardise by total forest fraction (to allow to use forest fraction from another database)
    fage_out_decade_frac(:,:,:,nn)=fage_out_decade(:,:,:,nn)./repmat(fage_out_decade_totfor(:,:,nn),[1 1 nages_dec]);
    %Convert to areas using provided forest fraction file
    fage_out_decade_area(:,:,:,nn)=fage_out_decade_frac(:,:,:,nn).*repmat(fmask.*garea,[1 1 nages_dec]);
end
clear nn

%Aggregate age distributions over regions
fage_out_decade_reg=NaN(nregion,nages_dec,nyout);
for nn=1:nyout
    for rr=1:nregion
        for aa=1:nages_dec
            fage_out_decade_area_sel=squeeze(fage_out_decade_area(:,:,aa,nn));
            fage_out_decade_reg(rr,aa,nn)=squeeze(nansum(fage_out_decade_area_sel(rmask==rr)))/1e12;
        end
    end
end
clear nn aa fage_out_decade_area_sel

if gfad_comp
    %Read in the GFAD data
    [gfad_fage_area_sum_reg_stan,gfad_fage_area_sum_reg_lower,gfad_fage_area_sum_reg_upper]=...
            gfad_region_read(fmask,garea,rmask,nregion,gfad_filepath_stan,gfad_filepath_lower,gfad_filepath_upper);
end

%Make plot
ages=5:10:nages;
figure
for nn=1:nregion
    ss(nn)=subplot(3,3,nn);
    hold on
    plot(ages(1:14),fage_out_decade_reg(nn,1:14),'r.-','markersize',15)
    plot(ages(15),fage_out_decade_reg(nn,15),'r.','markersize',15)
    if gfad_comp
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




%ALLOW DISTURBANCE RATE TO VARY IN TIME FROM YEAR1


