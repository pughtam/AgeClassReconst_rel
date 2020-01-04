% Script to make a raster mask of the Olson biomes based on the shapefile provided at
% https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
% Output is written to netcdf files at 1 degree and at 12 minutes
%
% T. Pugh
% 04.01.20

inputfiles='/Users/pughtam/data/Olson_biome_mask/wwf_terr_ecos.shp';

%Set domain boundary
latlim=[-90 90];
lonlim=[-180 180];

%Read in the shapefile data and convert to raster
biome_raster=NaN(180*5,360*5);
for bb=1:14
    fprintf('Processing biome %d\n',bb)
    
    Sb = shaperead(inputfiles,'UseGeoCoords',true,'Selector',...
           {@(v1) (v1 == bb),'BIOME'});

    nSb=length(Sb);
    
    for nn=1:nSb
        lat=Sb(nn).Lat;
        lon=Sb(nn).Lon;
        
        [Z, R] = vec2mtx(lat, lon, 5,latlim,lonlim, 'filled');
        
        biome_raster(Z<2)=bb;
        if mod(nn,100)==0
            fprintf('Processed %5.1f%% for biome %d\n',(nn/nSb)*100,bb)
        end
    end
end
clear bb

%Aggregate to 1 x 1 degree
biome_raster_1deg=NaN(180,360);
for xx=1:360
    for yy=1:180
        xx_s=(xx*5)-5+1;
        xx_e=xx*5;
        yy_s=(yy*5)-5+1;
        yy_e=yy*5;
        
        temp=biome_raster(yy_s:yy_e,xx_s:xx_e);
        biome_raster_1deg(yy,xx)=mode(temp(:));
    end
end
clear xx yy xx_s xx_e yy_s yy_e


%---
%Write out to netcdf files

%1 degree
outfile='Olson_biomes_1deg.nc';
nccreate(outfile,'Olson_biomes','Dimensions',{'longitude',360,'latitude',180})
nccreate(outfile,'latitude','Dimensions',{'latitude',180})
nccreate(outfile,'longitude','Dimensions',{'longitude',360})
ncwrite(outfile,'Olson_biomes',biome_raster_1deg');
ncwrite(outfile,'latitude',-89.5:1:89.5);
ncwrite(outfile,'longitude',-179.5:1:179.5);
ncwriteatt(outfile,'Olson_biomes','1','Tropical & Subtropical Moist Broadleaf Forests');
ncwriteatt(outfile,'Olson_biomes','2','Tropical & Subtropical Dry Broadleaf Forests');
ncwriteatt(outfile,'Olson_biomes','3','Tropical & Subtropical Coniferous Forests');
ncwriteatt(outfile,'Olson_biomes','4','Temperate Broadleaf & Mixed Forests');
ncwriteatt(outfile,'Olson_biomes','5','Temperate Conifer Forests');
ncwriteatt(outfile,'Olson_biomes','6','Boreal Forests/Taiga');
ncwriteatt(outfile,'Olson_biomes','7','Tropical & Subtropical Grasslands, Savannas & Shrublands');
ncwriteatt(outfile,'Olson_biomes','8','Temperate Grasslands, Savannas & Shrublands');
ncwriteatt(outfile,'Olson_biomes','9','Flooded Grasslands & Savannas');
ncwriteatt(outfile,'Olson_biomes','10','Montane Grasslands & Shrublands');
ncwriteatt(outfile,'Olson_biomes','11','Tundra');
ncwriteatt(outfile,'Olson_biomes','12','Mediterranean Forests, Woodlands & Scrub');
ncwriteatt(outfile,'Olson_biomes','13','Deserts & Xeric Shrublands');
ncwriteatt(outfile,'Olson_biomes','14','Mangroves');

ncwriteatt(outfile,'longitude','Units','degrees_east')
ncwriteatt(outfile,'latitude','Units','degrees_north')
ncwriteatt(outfile,'/','Comment','Created based on https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world using olson_biome_to_raster.m')
ncwriteatt(outfile,'/','Resolution','1 arc degree')
ncwriteatt(outfile,'/','Version','1')
ncwriteatt(outfile,'/','Contact','T.A.M. Pugh, t.a.m.pugh@bham.ac.uk')

%12 minute
outfile='Olson_biomes_12min.nc';
nccreate(outfile,'Olson_biomes','Dimensions',{'longitude',360*5,'latitude',180*5})
nccreate(outfile,'latitude','Dimensions',{'latitude',180*5})
nccreate(outfile,'longitude','Dimensions',{'longitude',360*5})
ncwrite(outfile,'Olson_biomes',biome_raster');
ncwrite(outfile,'latitude',-89.9:0.2:89.9);
ncwrite(outfile,'longitude',-179.9:0.2:179.9);
ncwriteatt(outfile,'Olson_biomes','1','Tropical & Subtropical Moist Broadleaf Forests');
ncwriteatt(outfile,'Olson_biomes','2','Tropical & Subtropical Dry Broadleaf Forests');
ncwriteatt(outfile,'Olson_biomes','3','Tropical & Subtropical Coniferous Forests');
ncwriteatt(outfile,'Olson_biomes','4','Temperate Broadleaf & Mixed Forests');
ncwriteatt(outfile,'Olson_biomes','5','Temperate Conifer Forests');
ncwriteatt(outfile,'Olson_biomes','6','Boreal Forests/Taiga');
ncwriteatt(outfile,'Olson_biomes','7','Tropical & Subtropical Grasslands, Savannas & Shrublands');
ncwriteatt(outfile,'Olson_biomes','8','Temperate Grasslands, Savannas & Shrublands');
ncwriteatt(outfile,'Olson_biomes','9','Flooded Grasslands & Savannas');
ncwriteatt(outfile,'Olson_biomes','10','Montane Grasslands & Shrublands');
ncwriteatt(outfile,'Olson_biomes','11','Tundra');
ncwriteatt(outfile,'Olson_biomes','12','Mediterranean Forests, Woodlands & Scrub');
ncwriteatt(outfile,'Olson_biomes','13','Deserts & Xeric Shrublands');
ncwriteatt(outfile,'Olson_biomes','14','Mangroves');

ncwriteatt(outfile,'longitude','Units','degrees_east')
ncwriteatt(outfile,'latitude','Units','degrees_north')
ncwriteatt(outfile,'/','Comment','Created based on https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world using olson_biome_to_raster.m')
ncwriteatt(outfile,'/','Resolution','12 arc minutes')
ncwriteatt(outfile,'/','Version','1')
ncwriteatt(outfile,'/','Contact','T.A.M. Pugh, t.a.m.pugh@bham.ac.uk')
        