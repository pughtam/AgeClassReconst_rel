function [luh2_forlu_loss_prim_1deg,luh2_forlu_loss_sec_1deg]=luh2_forloss_read(filepath,year1ind,nyear,inc_woodharv)
%Load in the LUH2 transitions data for a specified period
%Sum up the forest loss transitions fractions
%Aggregate to 1 degree resolution
%
%filepath should be the path to the transitions.nc file as downloaded from http://luh.umd.edu/data.shtml
%
%T. Pugh
%10.08.19

primf_to_secdn=ncread(filepath,'primf_to_secdn',[1 1 year1ind],[Inf Inf nyear]);
primf_to_secdn(primf_to_secdn>1e19)=NaN;
luh2_forlu_loss_prim=primf_to_secdn;
clear primf_to_secdn

primf_to_urban=ncread(filepath,'primf_to_urban',[1 1 year1ind],[Inf Inf nyear]);
primf_to_urban(primf_to_urban>1e19)=NaN;
luh2_forlu_loss_prim=luh2_forlu_loss_prim+primf_to_urban;
clear primf_to_urban

primf_to_c3ann=ncread(filepath,'primf_to_c3ann',[1 1 year1ind],[Inf Inf nyear]);
primf_to_c3ann(primf_to_c3ann>1e19)=NaN;
luh2_forlu_loss_prim=luh2_forlu_loss_prim+primf_to_c3ann;
clear primf_to_c3ann

primf_to_c4ann=ncread(filepath,'primf_to_c4ann',[1 1 year1ind],[Inf Inf nyear]);
primf_to_c4ann(primf_to_c4ann>1e19)=NaN;
luh2_forlu_loss_prim=luh2_forlu_loss_prim+primf_to_c4ann;
clear primf_to_c4ann

primf_to_c3per=ncread(filepath,'primf_to_c3per',[1 1 year1ind],[Inf Inf nyear]);
primf_to_c3per(primf_to_c3per>1e19)=NaN;
luh2_forlu_loss_prim=luh2_forlu_loss_prim+primf_to_c3per;
clear primf_to_c3per

primf_to_c4per=ncread(filepath,'primf_to_c4per',[1 1 year1ind],[Inf Inf nyear]);
primf_to_c4per(primf_to_c4per>1e19)=NaN;
luh2_forlu_loss_prim=luh2_forlu_loss_prim+primf_to_c4per;
clear primf_to_c4per

primf_to_c3nfx=ncread(filepath,'primf_to_c3nfx',[1 1 year1ind],[Inf Inf nyear]);
primf_to_c3nfx(primf_to_c3nfx>1e19)=NaN;
luh2_forlu_loss_prim=luh2_forlu_loss_prim+primf_to_c3nfx;
clear primf_to_c3nfx

primf_to_pastr=ncread(filepath,'primf_to_pastr',[1 1 year1ind],[Inf Inf nyear]);
primf_to_pastr(primf_to_pastr>1e19)=NaN;
luh2_forlu_loss_prim=luh2_forlu_loss_prim+primf_to_pastr;
clear primf_to_pastr

primf_to_range=ncread(filepath,'primf_to_range',[1 1 year1ind],[Inf Inf nyear]);
primf_to_range(primf_to_range>1e19)=NaN;
luh2_forlu_loss_prim=luh2_forlu_loss_prim+primf_to_range;
clear primf_to_range

if inc_woodharv
    primf_harv=ncread(filepath,'primf_harv',[1 1 year1ind],[Inf Inf nyear]);
    primf_harv(primf_harv>1e19)=NaN;
    luh2_forlu_loss_prim=luh2_forlu_loss_prim+primf_harv;
    clear primf_harv
end



secdf_to_secdn=ncread(filepath,'secdf_to_secdn',[1 1 year1ind],[Inf Inf nyear]);
secdf_to_secdn(secdf_to_secdn>1e19)=NaN;
luh2_forlu_loss_sec=secdf_to_secdn;
clear secdf_to_secdn

secdf_to_urban=ncread(filepath,'secdf_to_urban',[1 1 year1ind],[Inf Inf nyear]);
secdf_to_urban(secdf_to_urban>1e19)=NaN;
luh2_forlu_loss_sec=luh2_forlu_loss_sec+secdf_to_urban;
clear secdf_to_urban

secdf_to_c3ann=ncread(filepath,'secdf_to_c3ann',[1 1 year1ind],[Inf Inf nyear]);
secdf_to_c3ann(secdf_to_c3ann>1e19)=NaN;
luh2_forlu_loss_sec=luh2_forlu_loss_sec+secdf_to_c3ann;
clear secdf_to_c3ann

secdf_to_c4ann=ncread(filepath,'secdf_to_c4ann',[1 1 year1ind],[Inf Inf nyear]);
secdf_to_c4ann(secdf_to_c4ann>1e19)=NaN;
luh2_forlu_loss_sec=luh2_forlu_loss_sec+secdf_to_c4ann;
clear secdf_to_c4ann

secdf_to_c3per=ncread(filepath,'secdf_to_c3per',[1 1 year1ind],[Inf Inf nyear]);
secdf_to_c3per(secdf_to_c3per>1e19)=NaN;
luh2_forlu_loss_sec=luh2_forlu_loss_sec+secdf_to_c3per;
clear secdf_to_c3per

secdf_to_c4per=ncread(filepath,'secdf_to_c4per',[1 1 year1ind],[Inf Inf nyear]);
secdf_to_c4per(secdf_to_c4per>1e19)=NaN;
luh2_forlu_loss_sec=luh2_forlu_loss_sec+secdf_to_c4per;
clear secdf_to_c4per

secdf_to_c3nfx=ncread(filepath,'secdf_to_c3nfx',[1 1 year1ind],[Inf Inf nyear]);
secdf_to_c3nfx(secdf_to_c3nfx>1e19)=NaN;
luh2_forlu_loss_sec=luh2_forlu_loss_sec+secdf_to_c3nfx;
clear secdf_to_c3nfx

secdf_to_pastr=ncread(filepath,'secdf_to_pastr',[1 1 year1ind],[Inf Inf nyear]);
secdf_to_pastr(secdf_to_pastr>1e19)=NaN;
luh2_forlu_loss_sec=luh2_forlu_loss_sec+secdf_to_pastr;
clear secdf_to_pastr

secdf_to_range=ncread(filepath,'secdf_to_range',[1 1 year1ind],[Inf Inf nyear]);
secdf_to_range(secdf_to_range>1e19)=NaN;
luh2_forlu_loss_sec=luh2_forlu_loss_sec+secdf_to_range;
clear secdf_to_range

%Regrid to 1 x 1 degree
luh2_forlu_loss_prim_1deg=zeros(360,180,nyear);
luh2_forlu_loss_sec_1deg=zeros(360,180,nyear);
for xx=1:360
    for yy=1:180
        ind_x=(xx*4)-3;
        ind_y=(yy*4)-3;
        temp_prim=luh2_forlu_loss_prim(ind_x:ind_x+3,ind_y:ind_y+3,:);
        temp_sec=luh2_forlu_loss_sec(ind_x:ind_x+3,ind_y:ind_y+3,:);
        luh2_forlu_loss_prim_1deg(xx,yy,:)=nansum(nansum(temp_prim,2),1)/16; %Calculate mean transition fraction, preserving gridcell area
        luh2_forlu_loss_sec_1deg(xx,yy,:)=nansum(nansum(temp_sec,2),1)/16;
        clear temp_sec temp_prim
        clear ind_x ind_y
    end
    clear yy
end
clear xx
