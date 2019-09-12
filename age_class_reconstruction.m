%Script to to construct the forest age distribution based on land-use change from LUH2 and on random
%disturbances as defined by Pugh et al. (2019).
%
%T. Pugh
%12.09.19

%---
%Input options
year1=1900; %First calendar year of simulation for which output is required (i.e. after spin-up)
nyear=115; %Total number of years of simulation for which output is required
inc_luh2=true; %Include the LUH2 transitions?
inc_woodharv=false; %Include the wood harvest transitions? Standard assumption is inc_woodharv=false (only functional if inc_luh2=true)
inc_dist=true; %Include background disturbance (true) or just LUH2 transitions (false)

output_crosscheck_plots=0; %Make diagnostic cross-check plots
output_years=[1900 1950 2015]; %Years for which to provide outputs
gfad_comp=true; %Include GFAD in the output plots
hansenmask=false; %Mask results according to year 2000 closed-canopy forest cover?

regmask='RECCAP'; %Region mask to use (ESA or RECCAP)

loadinputdata=true; %Do not load input data if equal to false (for rapid reruns when input data is already in memory)

outputcsv=true;
csvname='age_reconstruction_luh2dist_global_unmasked_multisens_v3.csv';

%---
%Settings
nages=150; %Number of age classes (1 class = 1 year)

%LUH2 input
luh2_filepath_trans='/data/LUH2/transitions.nc';
luh2_filepath_states='/data/LUH2/states.nc';
firstluh2year=850; %First calendar year in LUH2 dataset

%Background disturbance return time input
distfile='/Users/pughtam/Documents/GAP_work/Disturbance/netcdfs_for_deposition/tauO/tauO_standard_forest-area_LUcorrected.nc';

%Forest mask file
formaskfile='/Users/pughtam/Documents/GAP_work/Disturbance/netcdfs_for_deposition/forestmask/hansen_forested_frac_1deg.nc4';

%GFAD data (only used for comparison in output plots
gfad_filepath_stan='/data/GFAD_V1-1/GFAD_V1-1.nc';
gfad_filepath_lower='/data/GFAD_V1-1/GFAD_V1-1_lowerbound.nc';
gfad_filepath_upper='/data/GFAD_V1-1/GFAD_V1-1_upperbound.nc';

nyout=length(output_years);
nages_dec=nages/10;

%---
if loadinputdata
    %Read in necessary data
    
    %Read initial forest state
    nspinup=nages; %Start nages years before first output year to allow for spin-up of age classes (outputs are insensitive to longer spinup, even for disturbance)
    year1ind=year1-nspinup-firstluh2year+1;
    [primf_init,secdf_init]=luh2_forstates_read(luh2_filepath_states,year1ind);
    
    %Read forest loss and gain transitions from LUH2
    nyearluh2=nyear+nspinup;
    [luh2_forlu_loss_prim_1deg,luh2_forlu_loss_sec_1deg]=luh2_forloss_read(luh2_filepath_trans,year1ind,nyearluh2,inc_woodharv);
    luh2_forlu_gain_1deg=luh2_forgain_read(luh2_filepath_trans,year1ind,nyearluh2,inc_woodharv);
    
    %save luh_array_read.mat
    
    %Read in disturbance return period
    distint=ncread(distfile,'tauO');
    distint=fliplr(distint);
    
    %Assign disturbance rates for areas outside of the closed-canopy mask.
    %Do this using medians for ESA forest types
    
    addpath('/data/ESA_landcover/')
    [rmask,~,nregion]=esa_forest_9regions_new_1deg_func(false);
    rmask=fliplr(rmask');
    
    distint_fill=ones(size(distint))*100; %Background value of 100 years for any oulying areas
    for nn=1:nregion
        distmed=nanmedian(distint(rmask==nn));
        distint_fill(rmask==nn)=distmed;
    end
    clear nn distmed
    distint_fill(isnan(distint)==0)=distint(isnan(distint)==0);
    clear rmask nregion
    distint=distint_fill;
    clear distint_fill
    
end

%---
%Calculations section
%Track primary and secondary forest from LUH2 separately and only merge together in the final output array.

simulation_names={'Baseline','Baseline old-skew','Baseline young-skew',...
    '2to1','2to1 old-skew','2to1 young-skew',...
    '0p5to1','0p5to1 old-skew','0p5to1 young-skew'};

%Baseline simulation with uniform weighting
LUC_age_weights=0; %Whether to have equal likelihood of LUC conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=0; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=false; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=1; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=1; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

fage_out_decade_base=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec);

%Baseline simulation with weighting towards older age classes
LUC_age_weights=1; %Whether to have equal likelihood of LUC conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=1; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=false; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=1; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=1; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

fage_out_decade_base_oldweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec);

%Baseline simulation with weighting towards younger age classes
LUC_age_weights=2; %Whether to have equal likelihood of LUC conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=2; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=false; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=1; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=1; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

fage_out_decade_base_youngweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec);


%2to1 simulation with uniform weighting
LUC_age_weights=0; %Whether to have equal likelihood of LUC conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=0; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=true; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=2; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=1; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

fage_out_decade_2to1=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec);

%2to1 simulation with weighting towards older age classes
LUC_age_weights=1; %Whether to have equal likelihood of LUC conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=1; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=true; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=2; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=1; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

fage_out_decade_2to1_oldweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec);

%2to1 simulation with weighting towards younger age classes
LUC_age_weights=2; %Whether to have equal likelihood of LUC conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=2; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=true; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=2; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=1; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

fage_out_decade_2to1_youngweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec);


%0p5to1 simulation with uniform weighting
LUC_age_weights=0; %Whether to have equal likelihood of LUC conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=0; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=true; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=0.5; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=1; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

fage_out_decade_0p5to1=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec);

%0p5to1 simulation with weighting towards older age classes
LUC_age_weights=1; %Whether to have equal likelihood of LUC conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=1; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=true; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=0.5; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=1; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

fage_out_decade_0p5to1_oldweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec);

%0p5to1 simulation with weighting towards younger age classes
LUC_age_weights=2; %Whether to have equal likelihood of LUC conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=2; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=true; %Modify the background disturbance rate by a multiplicative scenario
dist_scen_start=0.5; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
dist_scen_end=1; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)

fage_out_decade_0p5to1_youngweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    dist_scen_start,dist_scen_end,nyout,nages_dec);

%Make a combined array for postprocessing
fage_out_decade=cat(5,fage_out_decade_base,fage_out_decade_base_oldweight,fage_out_decade_base_youngweight,...
    fage_out_decade_2to1,fage_out_decade_2to1_oldweight,fage_out_decade_2to1_youngweight,...
    fage_out_decade_0p5to1,fage_out_decade_0p5to1_oldweight,fage_out_decade_0p5to1_youngweight);
%fage_out_decade=cat(5,fage_out_decade_base,fage_out_decade_base_oldweight);

size_fage_out_decade=size(fage_out_decade);
nsens=size_fage_out_decade(5);
clear size_fage_out_decade

%---
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
%Calculate regional and global age distributions

%Read in forest area
if hansenmask
    fmask=fliplr(ncread(formaskfile,'forested_50_percent'));
    fmask=double(fmask)./100;
else
    fmask=ones(360,180);
end

%Calculate unmasked grid-cell area
garea=global_grid_area_1deg()';

%Load region mask
if strcmp(regmask,'ESA')
    addpath('/data/ESA_landcover/')
    [rmask,regions,nregion]=esa_forest_9regions_new_1deg_func(false);
    rmask=fliplr(rmask');
elseif strcmp(regmask,'RECCAP')
    rmask=fliplr(ncread('/data/Masks_etc/RECCAP_MASK11_Mask.nc','Region_Map'));
    regions={'North America','South America','Europe','Africa','Russia','Middle East',...
        'China/Japan','South-East Asia','Australasia'};
    nregion=length(regions);
else
    error('Region mask setting is invalid')
end

%Convert fractions to areas
fage_out_decade_totfor=squeeze(sum(fage_out_decade_base,3)); %Total forest fraction
fage_out_decade_frac=NaN(size(fage_out_decade));
fage_out_decade_area=NaN(size(fage_out_decade));
for nn=1:nyout
    if hansenmask
        %Standardise by total forest fraction (i.e. convert to fraction of gridcell to allow to use forest fraction from another database)
        fage_out_decade_frac(:,:,:,nn,:)=fage_out_decade(:,:,:,nn,:)./repmat(fage_out_decade_totfor(:,:,nn),[1 1 nages_dec 1 nsens]);
    else
        fage_out_decade_frac(:,:,:,nn,:)=fage_out_decade(:,:,:,nn,:);
    end
    %Convert to areas using provided forest fraction file
    fage_out_decade_area(:,:,:,nn,:)=fage_out_decade_frac(:,:,:,nn,:).*repmat(fmask.*garea,[1 1 nages_dec 1 nsens]);
end
clear nn

%Aggregate age distributions over regions
fage_out_decade_reg=NaN(nregion,nages_dec,nyout,nsens);
for ss=1:nsens
    for nn=1:nyout
        for rr=1:nregion
            for aa=1:nages_dec
                fage_out_decade_area_sel=squeeze(fage_out_decade_area(:,:,aa,nn,ss));
                fage_out_decade_reg(rr,aa,nn,ss)=squeeze(nansum(fage_out_decade_area_sel(rmask==rr)))/1e12;
            end
        end
    end
end
clear ss nn aa rr fage_out_decade_area_sel

%Aggregate age distributions over the globe
fage_out_decade_globe=NaN(nages_dec,nyout,nsens);
for ss=1:nsens
    for nn=1:nyout
        for aa=1:nages_dec
            fage_out_decade_area_sel=squeeze(fage_out_decade_area(:,:,aa,nn,ss));
            fage_out_decade_globe(aa,nn,ss)=squeeze(nansum(fage_out_decade_area_sel(:)))/1e12;
        end
    end
end
clear ss nn aa fage_out_decade_area_sel

%Calculate percentages
fage_out_decade_reg_perc=(fage_out_decade_reg./repmat(sum(fage_out_decade_reg,2),[1 nages_dec 1]))*100;
fage_out_decade_globe_perc=(fage_out_decade_globe./repmat(sum(fage_out_decade_globe,1),[nages_dec 1]))*100;

%---
%Make plot by global, and three selected regions
ages=5:10:nages;

if strcmp(regmask,'ESA')
    regsel=[1 5 6];
elseif strcmp(regmask,'RECCAP')
    regsel=[1 3 4];
    %regsel=[1 3 5];
end

%Percentages of forest area for baseline simulation
yname='% forest area';
age_class_reconstruction_regionalsel_plot_func(fage_out_decade_globe_perc(1:14,:,1),fage_out_decade_reg_perc(:,1:14,:,1),regsel,...
    ages(1:14),yname,nyout,regions);

%Absolute values for baseline simulation
yname='Forest area (Mkm^2)';
age_class_reconstruction_regionalsel_plot_func(fage_out_decade_globe(:,:,1),fage_out_decade_reg(:,:,:,1),regsel,...
    ages,yname,nyout,regions);

%---
if outputcsv
    %Output csv file with global outputs for each simulation
    fid=fopen(csvname,'w');
    fprintf(fid,'Units: million km2 closed-canopy forest area\n');
    fprintf(fid,'%d\n',output_years(1));
    fprintf(fid,'Simulation,1-10,11-20,21-30,31-40,41-50,51-60,61-70,71-80,81-90,91-100,101-110,111-120,121-130,131-140,OG\n');
    for ss=1:nsens
        fprintf(fid,'%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',...
            simulation_names{ss},fage_out_decade_globe(:,1,ss));
    end
    clear ss
    
    fprintf(fid,'%d\n',output_years(2));
    fprintf(fid,'Simulation,1-10,11-20,21-30,31-40,41-50,51-60,61-70,71-80,81-90,91-100,101-110,111-120,121-130,131-140,OG\n');
    for ss=1:nsens
        fprintf(fid,'%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',...
            simulation_names{ss},fage_out_decade_globe(:,2,ss));
    end
    clear ss
    
    fprintf(fid,'%d\n',output_years(3));
    fprintf(fid,'Simulation,1-10,11-20,21-30,31-40,41-50,51-60,61-70,71-80,81-90,91-100,101-110,111-120,121-130,131-140,OG\n');
    for ss=1:nsens
        fprintf(fid,'%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',...
            simulation_names{ss},fage_out_decade_globe(:,3,ss));
    end
    clear ss
    fclose(fid);
end

%Simple plot based on global outputs
figure
subplot(2,1,1)
hold on
plot(output_years,sum(fage_out_decade_globe(1:14,:,1),1),'b-o')
plot(output_years,fage_out_decade_globe(15,:,1),'k-o')
plot(output_years,sum(fage_out_decade_globe(1:14,:,2),1),'b^')
plot(output_years,sum(fage_out_decade_globe(1:14,:,3),1),'bv')
plot(output_years,fage_out_decade_globe(15,:,2),'k^')
plot(output_years,fage_out_decade_globe(15,:,3),'kv')
ylabel('Forest area (M km^2)')
legend('Regrowth','Old-growth')

subplot(2,1,2)
hold on
plot(output_years,sum(fage_out_decade_globe(1:14,:,1),1)./sum(fage_out_decade_globe(1:15,:,1),1),'k-o')
plot(output_years,sum(fage_out_decade_globe(1:14,:,2),1)./sum(fage_out_decade_globe(1:15,:,3),1),'k^')
plot(output_years,sum(fage_out_decade_globe(1:14,:,3),1)./sum(fage_out_decade_globe(1:15,:,2),1),'kv')
ylabel('Ratio young:old forests')

%---
if gfad_comp
    %Make plot by all regions including GFAD
    %NOTE: If the number/type of simulations change, the arrays passed to
    %the plotting function will need editing.
    
    %Read in the GFAD data
    if hansenmask
        fmask_for_gfad=fmask;
    else
        %Use the LUH2 mask for the most recent output year
        fmask_for_gfad=nansum(fage_out_decade(:,:,:,nyout),3);
    end
    [gfad_fage_area_sum_reg_stan,gfad_fage_area_sum_reg_lower,gfad_fage_area_sum_reg_upper,...
        fage_area_sum_globe_stan,fage_area_sum_globe_lower,fage_area_sum_globe_upper]=...
        gfad_region_read(fmask_for_gfad,garea,rmask,nregion,gfad_filepath_stan,gfad_filepath_lower,gfad_filepath_upper);
    
    %Baseline simulations
    plotyear=3; %Which output year to plot for
    plot_gfad=true;
    age_class_reconstruction_9region_gfad_plot_func(fage_out_decade_reg(:,:,:,1:3),...
        gfad_fage_area_sum_reg_stan,gfad_fage_area_sum_reg_lower,gfad_fage_area_sum_reg_upper,...
        ages,regions,nregion,plotyear,plot_gfad)
    
    %2to1 simulations
    plotyear=3; %Which output year to plot for
    plot_gfad=true;
    age_class_reconstruction_9region_gfad_plot_func(fage_out_decade_reg(:,:,:,4:6),...
        gfad_fage_area_sum_reg_stan,gfad_fage_area_sum_reg_lower,gfad_fage_area_sum_reg_upper,...
        ages,regions,nregion,plotyear,plot_gfad)
    
    %0p5to1 simulations
    plotyear=3; %Which output year to plot for
    plot_gfad=true;
    age_class_reconstruction_9region_gfad_plot_func(fage_out_decade_reg(:,:,:,7:9),...
        gfad_fage_area_sum_reg_stan,gfad_fage_area_sum_reg_lower,gfad_fage_area_sum_reg_upper,...
        ages,regions,nregion,plotyear,plot_gfad)
    
end
