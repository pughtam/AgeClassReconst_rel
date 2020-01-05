% Script to construct the forest stand age distribution based on gridded land-use change data and wood harvest from LUHv2h 
% and/or on a gridded disturbance rate as defined by Pugh et al. (2019, Nature Geoscience 12, 730-735). Outputs csv files of
% age distributions at regional and global aggregations.
%
% Age is categorised in annual bins up to 140 years old, after which it is designated as old-growth. Output is then grouped
% into 10-year bins for display purposes. Note that simulation must therefore start 140 years before the first year for which
% an age distribution is required.
%
% An implicit assumption when using the LUHv2h data is that all land-use change or wood harvest actions in that dataset are
% enacted through transitioning patches of land at least 0.1 ha in size.
%
% Scenarios can be constructed, multiplying non-LUC or wood harvest disturbance rates by an arbitrary factor. Different
% assumptions for the likelihood of disturbance with stand age are also tested.
%
% All calculations are made on a 1 x 1 degree global grid, aggregating the LUHv2h and forest mask data to this resolution.
%
% LUHv2h data from:
% https://luh.umd.edu/data.shtml
%
% Gridded disturbance rate data (based on 2000-2014 observations) from:
% https://dataguru.lu.se/doi:10.18161/disturbance_tauo.201905
%
% Year 2000 closed canopy forest cover from:
% https://dataguru.lu.se/doi:10.18161/disturbance_forestmask.201905
% (described in Pugh et al. 2019, Nature Geoscience 12, 730-735; based on data from Hansen et al. 2013, Science 342, 850-853)
%
% Region mask based on Olson et al. (2001, Bioscience 51(11), 933-938) and converted to raster using olson_biom_to_raster.m
%
% Dependencies:
% - luh2_forstates_read.m
% - luh2_forloss_read.m
% - luh2_forgain_read.m
% - global_grid_area_1deg.m
% - age_class_reconstruction_func.m
%
% T. A .M. Pugh
% t.a.m.pugh@bham.ac.uk
% 12.09.19

%---
% Input options
year1=1900; %First calendar year of simulation for which output is required (i.e. after spin-up)
nyear=115; %Total number of years of simulation for which output is required (if a future scenario is used, this is added to this variable later)
inc_luh2=true; %Include the LUH2 transitions?
inc_luh2futscen=true; %Include an LUH2 future scenario to 2100?
inc_woodharv=false; %Include the wood harvest transitions? Standard assumption is inc_woodharv=false (only functional if inc_luh2=true)
inc_dist=true; %Include background disturbance (true) or just LUH2 transitions (false)

inc_agesens=false; %Include sensitivity tests on preferrential age of conversion
inc_distsens=true; %Include disturbance sensitivity scenarios (settings for these are in the Calculations section, below)

output_crosscheck_plots=0; %Make diagnostic cross-check plots between states and transitions
ccanopymask=true; %Mask results according to year 2000 closed-canopy forest cover?

loadinputdata=true; %Do not load input data if equal to false (for rapid reruns when input data is already in memory)

outputcsv=true; %Whether to output csv files with summary data
csvname_stub='age_reconstruction_luh2dist_masked_fut_scen';


%---
% Settings and filepaths

nages=150; %Number of age classes (1 class = 1 year)

% LUH2 historical input
luh2_filepath_trans='/data/LUH2/transitions.nc';
luh2_filepath_states='/data/LUH2/states.nc';
firstluh2year=850; %First calendar year in LUH2 dataset

% LUH2 future input
luh2_filepath_trans_fut='/data/LUH2/multiple-transitions_input4MIPs_landState_ScenarioMIP_UofMD-GCAM-ssp434-2-1-f_gn_2015-2100.nc';

% Background disturbance return time input
distfile='/Users/pughtam/Documents/GAP_and_other_work/Disturbance/netcdfs_for_deposition/tauO/tauO_standard_forest-area_LUcorrected.nc';

% Closed-canopy forest mask file 
formaskfile='/Users/pughtam/Documents/GAP_and_other_work/Disturbance/netcdfs_for_deposition/forestmask/hansen_forested_frac_0p5deg.nc4';

%Region mask file
regmaskfile='/data/Olson_biome_mask/Olson_biomes_1deg.nc';

% Years for which to provide outputs (holding for all years is too memory heavy)
if inc_luh2futscen
    output_years=[1900:10:2010 2015 2020:10:2100];
else
    output_years=[1900:10:2010 2015];
end

nyout=length(output_years);
nages_dec=nages/10;


%---
% Basic checks

if inc_woodharv && inc_dist
    error('Using the background disturbance datasets and wood harvest from LUH2 double-counts wood harvest. Set either inc_woodharv or inc_dist to false')
end

if inc_dist && ~ccanopymask
    ccanopymask=true;
    fprintf('WARNING: Setting ccanopymask=true because inc_dist==true and disturbance dataset is only valid for closed-canopy forest')
end
    

%---
% Now beginning processing

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
    
    %Add future scenario data to luh2 arrays if necessary
    if inc_luh2futscen
        if year1==1900 && nyear==115
            nyear=nyear+85;
        else
            error('If inc_luh2futscen==true, year1 must be 1900 and nyear must be 115')
        end
        [luh2_forlu_loss_fut_prim_1deg,luh2_forlu_loss_fut_sec_1deg]=luh2_forloss_read(luh2_filepath_trans_fut,1,85,inc_woodharv);
        luh2_forlu_gain_fut_1deg=luh2_forgain_read(luh2_filepath_trans_fut,1,85,inc_woodharv);

        luh2_forlu_loss_prim_1deg=cat(3,luh2_forlu_loss_prim_1deg,luh2_forlu_loss_fut_prim_1deg);
        luh2_forlu_loss_sec_1deg=cat(3,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_fut_sec_1deg);
        luh2_forlu_gain_1deg=cat(3,luh2_forlu_gain_1deg,luh2_forlu_gain_fut_1deg);
        clear luh2_forlu_loss_fut_prim_1deg luh2_forlu_loss_fut_sec_1deg luh2_forlu_gain_fut_1deg
    end
    
    %Read in forest area
    if ccanopymask
        fmask_0p5deg=fliplr(ncread(formaskfile,'forested_50_percent'));
        fmask_0p5deg=double(fmask_0p5deg)./100;
        %Aggregate to 1 degree
        fmask=NaN(360,180);
        for xx=1:360
            for yy=1:180
                xx_s=xx*2-1;
                xx_e=xx*2;
                yy_s=yy*2-1;
                yy_e=yy*2;
                temp=fmask_0p5deg(xx_s:xx_e,yy_s:yy_e);
                fmask(xx,yy)=nansum(temp(:))/4;
            end
        end
        clear xx yy xx_s xx_e yy_s yy_e
    else
        fmask=ones(360,180);
    end
    
    %Read in disturbance return period
    if inc_dist
        distint_in=ncread(distfile,'tauO');
        distint_in=fliplr(distint_in);
        
        %Assign disturbance rates for areas outside of the mask used in the disturbance return period dataset using the nearest
        %neighbour rule in order to have values for all gridcells containing closed-canopy forest (as in Pugh et al., 2019,
        %Nature Geoscience 12, 730-735).
        
        distint=distint_in;
        nansleft=length(find(isnan(distint)==1));
        while nansleft>1075 %1075 is number of cells around the edge of the domain which we don't fill
            temp=distint;
            
            for ii=2:359
                for jj=2:179
                    if isnan(distint(ii,jj))==1
                        temp_ext=temp(ii-1:ii+1,jj-1:jj+1);
                        distint(ii,jj)=nanmean(temp_ext(:));
                    end
                end
            end
            nansleft=length(find(isnan(distint)==1));
        end
        clear ii jj nansleft temp temp_ext
        
        %Do not extrapolate to gridcells without at least 5% closed-canopy forest cover to avoid overextrapolation.
        %For these gridcells assign the median global disturbance return period from the original dataset.
        distint(fmask<0.05)=nanmedian(distint_in(:));
    else
        distint=zeros(360,180);
    end

end

%---
% Calculations section

if inc_agesens && inc_distsens
    simulation_names={'Baseline','Baseline old-skew','Baseline young-skew',...
        'Sensitivity','Sensitivity old-skew','Sensitivity young-skew'};
elseif inc_agesens
    simulation_names={'Baseline','Baseline old-skew','Baseline young-skew'};
elseif inc_distsens
    simulation_names={'Baseline',...
        'Sensitivity'};
else
    simulation_names={'Baseline'};
end

% Baseline simulation with uniform weighting
LUC_age_weights=0; %Whether to have equal likelihood of LUC or wood harvest conversion for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
dist_age_weights=0; %Whether to have equal likelihood of disturbance for all ages classes (0), higher likelihood for older classes (1) or higher likelihood for younger classes (2)
use_dist_scen=false; %Modify the background disturbance rate by a multiplicative scenario?

fage_out_decade_base=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
    nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
    0,0,nyout,nages_dec,0,0);

if inc_agesens
    %Baseline simulation with weighting towards older age classes
    LUC_age_weights=1; 
    dist_age_weights=1;
    use_dist_scen=false;
    
    fage_out_decade_base_oldweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
        nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
        0,0,nyout,nages_dec,0,0);
    
    %Baseline simulation with weighting towards younger age classes
    LUC_age_weights=2;
    dist_age_weights=2;
    use_dist_scen=false;
    
    fage_out_decade_base_youngweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
        nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
        0,0,nyout,nages_dec,0,0);
end

if inc_distsens
    %Sensitivity simulation with uniform weighting
    LUC_age_weights=0;
    dist_age_weights=0;
    use_dist_scen=true;
    dist_scen_start=1; %Multiplier for background disturbance rate at year1 and during spin-up (if use_dist_scen=true)
    dist_scen_end=2; %Multiplier for background disturbance rate at end of simulation (if use_dist_scen=true)
    firstscenyear=2015; %Year to start disturbance scenario (if use_dist_scen=true)
    lastscenyear=2050; %Year to start disturbance scenario and fix rates thereafter (if use_dist_scen=true)
    
    fage_out_decade_sens=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
        nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
        dist_scen_start,dist_scen_end,nyout,nages_dec,firstscenyear,lastscenyear);
    
    if inc_agesens
        %Sensitivity simulation with weighting towards older age classes
        LUC_age_weights=1;
        dist_age_weights=1;
        use_dist_scen=true;
        dist_scen_start=1;
        dist_scen_end=2;
        
        fage_out_decade_sens_oldweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
            nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
            dist_scen_start,dist_scen_end,nyout,nages_dec,firstscenyear,lastscenyear);
        
        %Sensitivity simulation with weighting towards younger age classes
        LUC_age_weights=2;
        dist_age_weights=2;
        use_dist_scen=true;
        dist_scen_start=1;
        dist_scen_end=2;
        
        fage_out_decade_sens_youngweight=age_class_reconstruction_func(use_dist_scen,LUC_age_weights,inc_dist,inc_luh2,dist_age_weights,output_years,...
            nages,secdf_init,primf_init,luh2_forlu_loss_sec_1deg,luh2_forlu_loss_prim_1deg,luh2_forlu_gain_1deg,distint,nyear,nspinup,year1,...
            dist_scen_start,dist_scen_end,nyout,nages_dec,firstscenyear,lastscenyear);
    end
    
end


% Make a combined array for postprocessing
if inc_distsens
    if inc_agesens
        fage_out_decade=cat(5,fage_out_decade_base,fage_out_decade_base_oldweight,fage_out_decade_base_youngweight,...
            fage_out_decade_sens,fage_out_decade_sens_oldweight,fage_out_decade_sens_youngweight);
    else
        fage_out_decade=cat(5,fage_out_decade_base,...
            fage_out_decade_sens);
    end
else
    if inc_agesens
        fage_out_decade=cat(5,fage_out_decade_base,fage_out_decade_base_oldweight,fage_out_decade_base_youngweight);
    else
        fage_out_decade=fage_out_decade_base;
    end
end

size_fage_out_decade=size(fage_out_decade);
if length(size_fage_out_decade)==4
    nsens=1; %Total number of senstivity studies
else
    nsens=size_fage_out_decade(5);
end
clear size_fage_out_decade

%---
if output_crosscheck_plots
    %Cross-check total forest fraction with states file for end of the calculation. Difference should be no more than
    %rounding errors (influenced by minimum for "to_lose" in age_class_reconstruction_func.m).
    %NOTE: These should only agree if wood harvest is included in the transitions (i.e. inc_woodharv=1)
    [primf_end,secdf_end]=luh2_forstates_read(luh2_filepath_states,year1ind+nspinup+nyear);
    res_check=sum(fage_out_decade(:,:,:,3),3)-primf_end-secdf_end; %Assumes that year 2015 is 3rd index in array
    
    figure
    p1=pcolor(flipud(res_check'));
    set(p1,'linestyle','none')
    colorbar
end

%---
% Calculate regional and global age distributions

% Calculate unmasked grid-cell area
garea=global_grid_area_1deg()';

% Load region mask for tropical, temperate and boreal forests
rmask_raw=ncread(regmaskfile,'Olson_biomes');
rmask_raw=fliplr(rmask_raw);
rmask=zeros(size(rmask_raw));
rmask(rmask_raw==1 | rmask_raw==2 | rmask_raw==3 | rmask_raw==14)=1; %Tropical forest
rmask(rmask_raw==4 | rmask_raw==5 | rmask_raw==12)=2; %Temperate and Mediterranean forest
rmask(rmask_raw==6)=3; %Boreal forest and taiga
clear rmask_raw
nregion=3;
regions={'Tropical forest','Temperate and Mediterranean forest','Boreal forest and taiga'};
regions_short={'tropical','temperate','boreal'};

% Convert fractions to areas
fage_out_decade_totfor=squeeze(nansum(fage_out_decade,3)); %Total forest fraction
fage_out_decade_frac=NaN(size(fage_out_decade));
fage_out_decade_area=NaN(size(fage_out_decade));
for nn=1:nyout
    if ccanopymask
        %Standardise by total forest fraction (i.e. convert to fraction of gridcell to allow to use forest fraction from another database)
        fage_out_decade_frac(:,:,:,nn,:)=fage_out_decade(:,:,:,nn,:)./permute(repmat(fage_out_decade_totfor(:,:,nn,:),[1 1 1 1 nages_dec]),[1 2 5 3 4]);
    else
        fage_out_decade_frac(:,:,:,nn,:)=fage_out_decade(:,:,:,nn,:);
    end
    %Convert to areas using provided forest fraction file
    fage_out_decade_area(:,:,:,nn,:)=fage_out_decade_frac(:,:,:,nn,:).*repmat(fmask.*garea,[1 1 nages_dec 1 nsens]);
end
clear nn

% Aggregate age distributions over regions
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

% Aggregate age distributions over the globe
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

% Calculate percentages
fage_out_decade_reg_perc=(fage_out_decade_reg./repmat(sum(fage_out_decade_reg,2),[1 nages_dec 1]))*100;
fage_out_decade_globe_perc=(fage_out_decade_globe./repmat(sum(fage_out_decade_globe,1),[nages_dec 1]))*100;

%---
if outputcsv
    %Output csv file with global outputs for each simulation
    fid=fopen([csvname_stub,'_global.csv'],'w');
    if ccanopymask
        fprintf(fid,'Units: million km2 closed-canopy forest area\n');
    else
        fprintf(fid,'Units: million km2 forest area\n');
    end
    for nn=1:nyout
        fprintf(fid,'%d\n',output_years(nn));
        fprintf(fid,'Simulation,1-10,11-20,21-30,31-40,41-50,51-60,61-70,71-80,81-90,91-100,101-110,111-120,121-130,131-140,OG\n');
        for ss=1:nsens
            fprintf(fid,'%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',...
                simulation_names{ss},fage_out_decade_globe(:,nn,ss));
        end
        clear ss
    end
    clear nn
    fclose(fid);
    
    %Output csv files with regional outputs for each simulation
    for rr=1:nregion
        fid=fopen([csvname_stub,'_region_',regions_short{rr},'.csv'],'w');
        if ccanopymask
            fprintf(fid,'Units: million km2 closed-canopy forest area\n');
        else
            fprintf(fid,'Units: million km2 forest area\n');
        end
        for nn=1:nyout
            fprintf(fid,'%d\n',output_years(nn));
            fprintf(fid,'Simulation,1-10,11-20,21-30,31-40,41-50,51-60,61-70,71-80,81-90,91-100,101-110,111-120,121-130,131-140,OG\n');
            for ss=1:nsens
                fprintf(fid,'%s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',...
                    simulation_names{ss},fage_out_decade_reg(rr,:,nn,ss));
            end
            clear ss
        end
        clear nn
        fclose(fid);
    end
end

%Plot based on global numbers
figure
hold on
plot(output_years,sum(fage_out_decade_globe(1:14,:,1),1),'b-o')
plot(output_years,fage_out_decade_globe(15,:,1),'k-o')
if inc_agesens
    plot(output_years,sum(fage_out_decade_globe(1:14,:,2),1),'b^')
    plot(output_years,sum(fage_out_decade_globe(1:14,:,3),1),'bv')
    plot(output_years,fage_out_decade_globe(15,:,2),'k^')
    plot(output_years,fage_out_decade_globe(15,:,3),'kv')
end
ylabel('Forest area (M km^2)')
legend('Regrowth','Old-growth')
