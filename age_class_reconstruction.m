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
inc_woodharv=0; %Include the wood harvest transitions? Standard assumption is inc_woodharv=0
output_crosscheck_plots=0; %Make diagnostic cross-check plots
output_years=[1900 1950 2015]; %Years for which to provide outputs

%---
%Settings
nages=150; %Number of age classes (1 class = 1 year)

%LUH2 input
luh2_filepath_trans='/data/LUH2/transitions.nc';
luh2_filepath_states='/data/LUH2/states.nc';
firstluh2year=850; %First calendar year in LUH2 dataset

%Background disturbance return time input
distfile='/home/adf/pughtam/data/Disturbance/netcdfs_for_deposition/tauO/tauO_standard_forest-area_LUcorrected.nc';

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

    %Carry out random disturbance (equal probability across all ages of secondary and primary forest) and add
    %to youngest age class of secondary forest.
    %Disturb a fixed fraction per year defined by 1/distint
    %Secondary forest
    frac_dist_sec=fage_sec(:,:,2:nages).*repmat(1./distint,[1 1 nages-1]);
    fage_sec(:,:,2:nages)=fage_sec(:,:,2:nages)-frac_dist_sec;
    fage_sec(:,:,1)=fage_sec(:,:,1)+sum(frac_dist_sec,3);
    clear frac_dist_sec
    %Primary forest
    frac_dist_prim=fage_prim(:,:,2:nages).*repmat(1./distint,[1 1 nages-1]);
    fage_prim(:,:,2:nages)=fage_prim(:,:,2:nages)-frac_dist_prim;
    fage_prim(:,:,1)=fage_prim(:,:,1)+sum(frac_dist_prim,3);
    clear frac_dist_prim

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
fage_out_decade=NaN(360,180,nages/10,nyout);
for nn=1:nyout
    dd=0;
    for aa=1:10:nages
        dd=dd+1;
        fage_out_decade(:,:,dd,nn)=sum(fage_out(:,:,aa:aa+9,nn),3);
    end
    clear dd aa
end
clear nn

%Plot age distributions for a variety of regions


%Plot an age distribution example NOTE: PROBABLY NEED TO CORRECT BY AREA IF DOING FOR REGIONS
fage_europe=fage(160:220,20:60,:,265);
fage_europe_mean=NaN(nages,1);
for aa=1:nages
    temp=fage_europe(:,:,aa);
    fage_europe_mean(aa)=nanmean(temp(:));
end
clear aa temp
plot(1:nages-1,fage_europe_mean(1:nages-1))

nn=0;
fage_europe_mean_decade=NaN(nages/10,1);
for aa=1:10:nages
    nn=nn+1;
    fage_europe_mean_decade(nn)=sum(fage_europe_mean(aa:aa+9));
end
clear nn


%ALLOW DISTURBANCE RATE TO VARY IN TIME FROM YEAR1

%PLOT AGE DISTRIBUTIONS FOR A VARIETY OF STANDARDISED REGIONS (E.G. TRANSCOM)


%Remove spin-up info
%fage(:,:,:,1:nages)=[];
