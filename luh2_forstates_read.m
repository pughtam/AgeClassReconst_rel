function [primf_init_1deg,secdf_init_1deg]=luh2_forstates_read(filepath,year1ind)
%Load in the LUH2 states data for a specified period
%Aggregate to 1 degree resolution
%
%filepath should be the path to the states.nc file as downloaded from http://luh.umd.edu/data.shtml
%
%T. Pugh
%10.08.19

primf_init=ncread(filepath,'primf',[1 1 year1ind],[Inf Inf 1]);
primf_init(primf_init>1e19)=NaN;

secdf_init=ncread(filepath,'secdf',[1 1 year1ind],[Inf Inf 1]);
secdf_init(secdf_init>1e19)=NaN;

%Aggregate to 1 x 1 degree
primf_init_1deg=NaN(360,180);
secdf_init_1deg=NaN(360,180);
for xx=1:360
    for yy=1:180
        ind_x=(xx*4)-3;
        ind_y=(yy*4)-3;
        temp_prim=primf_init(ind_x:ind_x+3,ind_y:ind_y+3);
        temp_sec=secdf_init(ind_x:ind_x+3,ind_y:ind_y+3);
        primf_init_1deg(xx,yy)=nansum(temp_prim(:))/16; %Calculate mean transition fraction, preserving gridcell area
        secdf_init_1deg(xx,yy)=nansum(temp_sec(:))/16;
        clear temp_sec temp_prim
        clear ind_x ind_y
    end
    clear yy
end
clear xx