function luh2_forlu_gain_1deg=luh2_forgain_read(filepath,year1ind,nyear,inc_woodharv)
%Load in the LUH2 transitions data for a specified period
%Sum up the forest gain transitions fractions
%Aggregate to 1 degree resolution
%
%filepath should be the path to the transitions.nc file as downloaded from http://luh.umd.edu/data.shtml
%
%T. Pugh
%10.08.19

%---
%Do all processing in latitude bands to minimise memory requirements
luh2_forlu_gain_1deg=NaN(360,180,nyear);

llint=80; %index range for latitude bands

for lls=1:llint:720
    
    primn_to_secdf=ncread(filepath,'primn_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    primn_to_secdf(primn_to_secdf>1e19)=NaN;
    luh2_forlu_gain=primn_to_secdf;
    clear primn_to_secdf
    
    secdn_to_secdf=ncread(filepath,'secdn_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    secdn_to_secdf(secdn_to_secdf>1e19)=NaN;
    luh2_forlu_gain=luh2_forlu_gain+secdn_to_secdf;
    clear secdn_to_secdf
    
    urban_to_secdf=ncread(filepath,'urban_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    urban_to_secdf(urban_to_secdf>1e19)=NaN;
    luh2_forlu_gain=luh2_forlu_gain+urban_to_secdf;
    clear urban_to_secdf
    
    c3ann_to_secdf=ncread(filepath,'c3ann_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    c3ann_to_secdf(c3ann_to_secdf>1e19)=NaN;
    luh2_forlu_gain=luh2_forlu_gain+c3ann_to_secdf;
    clear c3ann_to_secdf
    
    c4ann_to_secdf=ncread(filepath,'c4ann_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    c4ann_to_secdf(c4ann_to_secdf>1e19)=NaN;
    luh2_forlu_gain=luh2_forlu_gain+c4ann_to_secdf;
    clear c4ann_to_secdf
    
    c3per_to_secdf=ncread(filepath,'c3per_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    c3per_to_secdf(c3per_to_secdf>1e19)=NaN;
    luh2_forlu_gain=luh2_forlu_gain+c3per_to_secdf;
    clear c3per_to_secdf
    
    c4per_to_secdf=ncread(filepath,'c4per_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    c4per_to_secdf(c4per_to_secdf>1e19)=NaN;
    luh2_forlu_gain=luh2_forlu_gain+c4per_to_secdf;
    clear c4per_to_secdf
    
    c3nfx_to_secdf=ncread(filepath,'c3nfx_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    c3nfx_to_secdf(c3nfx_to_secdf>1e19)=NaN;
    luh2_forlu_gain=luh2_forlu_gain+c3nfx_to_secdf;
    clear c3nfx_to_secdf
    
    pastr_to_secdf=ncread(filepath,'pastr_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    pastr_to_secdf(pastr_to_secdf>1e19)=NaN;
    luh2_forlu_gain=luh2_forlu_gain+pastr_to_secdf;
    clear pastr_to_secdf
    
    range_to_secdf=ncread(filepath,'range_to_secdf',[1 lls year1ind],[Inf llint nyear]);
    range_to_secdf(range_to_secdf>1e19)=NaN;
    luh2_forlu_gain=luh2_forlu_gain+range_to_secdf;
    clear range_to_secdf
    
    if inc_woodharv
        primf_harv=ncread(filepath,'primf_harv',[1 lls year1ind],[Inf llint nyear]);
        primf_harv(primf_harv>1e19)=NaN;
        luh2_forlu_gain=luh2_forlu_gain+primf_harv;
        clear primf_harv
        
        secmf_harv=ncread(filepath,'secmf_harv',[1 lls year1ind],[Inf llint nyear]);
        secmf_harv(secmf_harv>1e19)=NaN;
        luh2_forlu_gain=luh2_forlu_gain+secmf_harv;
        clear secmf_harv
        
        secyf_harv=ncread(filepath,'secyf_harv',[1 lls year1ind],[Inf llint nyear]);
        secyf_harv(secyf_harv>1e19)=NaN;
        luh2_forlu_gain=luh2_forlu_gain+secyf_harv;
        clear secyf_harv
    end
    
    %Aggregate to 1 x 1 degree
    indaggs=ceil(lls/4);
    indagge=ceil((lls+llint-1)/4);
    for xx=1:360
        for yy=indaggs:indagge
            ind_x=(xx*4)-3;
            %ind_y=(yy*4)-3;
            ind_y=(yy*4)-3-lls+1;
            temp=luh2_forlu_gain(ind_x:ind_x+3,ind_y:ind_y+3,:);
            luh2_forlu_gain_1deg(xx,yy,:)=nansum(nansum(temp,2),1)/16; %Calculate mean transition fraction, preserving gridcell area
            clear temp
            clear ind_x ind_y
        end
        clear yy
    end
    clear xx
    
end
