function [fage_area_sum_reg_stan,fage_area_sum_reg_lower,fage_area_sum_reg_upper]=...
            gfad_region_read(fmask,garea,rmask,nregion,gfad_filepath_stan,gfad_filepath_lower,gfad_filepath_upper)

maxclass=15;

%Read in the GFAD files
fage_stan=ncread(gfad_filepath_stan,'age'); %Units are fraction of the grid-cell occupied by that PFT and age class
fage_stan=squeeze(sum(fage_stan,3)); %Merge all the PFTs together

fage_lower_in=ncread(gfad_filepath_lower,'age');
fage_lower_in=squeeze(sum(fage_lower_in,3));
fage_lower=zeros(size(fage_stan)); %Make the classes in the lower bounda array consistent with the other arrays
fage_lower(:,:,1:9)=fage_lower_in(:,:,1:9);
fage_lower(:,:,15)=fage_lower_in(:,:,10);
clear fage_lower_in

fage_upper=ncread(gfad_filepath_upper,'age');
fage_upper=squeeze(sum(fage_upper,3));

%Aggregate to 1 degree
fage_stan_1deg=NaN(360,180,maxclass);
fage_lower_1deg=NaN(360,180,maxclass);
fage_upper_1deg=NaN(360,180,maxclass);
for xx=1:360
    for yy=1:180
        ind_x=(xx*2)-1;
        ind_y=(yy*2)-1;
        for nn=1:maxclass
            temp_stan=fage_stan(ind_x:ind_x+1,ind_y:ind_y+1,nn);
            fage_stan_1deg(xx,yy,nn)=nanmean(temp_stan(:));
            temp_lower=fage_lower(ind_x:ind_x+1,ind_y:ind_y+1,nn);
            fage_lower_1deg(xx,yy,nn)=nanmean(temp_lower(:));
            temp_upper=fage_upper(ind_x:ind_x+1,ind_y:ind_y+1,nn);
            fage_upper_1deg(xx,yy,nn)=nanmean(temp_upper(:));
        end
        clear temp_stan temp_lower temp_upper
        clear ind_x ind_y nn
    end
    clear yy
end
clear xx
clear fage_stan fage_lower fage_upper

%Convert fractions to areas
poul_totfor_stan=sum(fage_stan_1deg,3); %Total forest fraction
fage_frac_stan=fage_stan_1deg./repmat(poul_totfor_stan,[1 1 maxclass]); %Standardise by total forest fraction (to allow to use forest fraction from another database)
fage_area_stan=fage_frac_stan.*repmat(fmask.*garea,[1 1 maxclass]); %Convert to areas using provided forest fraction file
%fage_area_sum_stan=squeeze(sum(nansum(fage_area_stan,2),1))/1e12;

poul_totfor_lower=sum(fage_lower_1deg,3);
fage_frac_lower=fage_lower_1deg./repmat(poul_totfor_lower,[1 1 maxclass]);
fage_area_lower=fage_frac_lower.*repmat(fmask.*garea,[1 1 maxclass]);
%fage_area_sum_lower=squeeze(sum(nansum(fage_area_lower,2),1))/1e12;

poul_totfor_upper=sum(fage_upper_1deg,3);
fage_frac_upper=fage_upper_1deg./repmat(poul_totfor_upper,[1 1 maxclass]);
fage_area_upper=fage_frac_upper.*repmat(fmask.*garea,[1 1 maxclass]);
%fage_area_sum_upper=squeeze(sum(nansum(fage_area_upper,2),1))/1e12;

%ages=10:10:150;

fage_area_sum_reg_stan=NaN(nregion,maxclass);
fage_area_sum_reg_lower=NaN(nregion,maxclass);
fage_area_sum_reg_upper=NaN(nregion,maxclass);
for nn=1:nregion
    for aa=1:maxclass
        fage_area_sel=squeeze(fage_area_stan(:,:,aa));
        fage_area_sum_reg_stan(nn,aa)=squeeze(nansum(fage_area_sel(rmask==nn)))/1e12;
        fage_area_sel=squeeze(fage_area_lower(:,:,aa));
        fage_area_sum_reg_lower(nn,aa)=squeeze(nansum(fage_area_sel(rmask==nn)))/1e12;
        fage_area_sel=squeeze(fage_area_upper(:,:,aa));
        fage_area_sum_reg_upper(nn,aa)=squeeze(nansum(fage_area_sel(rmask==nn)))/1e12;
    end
end
clear nn aa fage_area_sel