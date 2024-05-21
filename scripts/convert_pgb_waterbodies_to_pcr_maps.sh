
#~ (pcrglobwb_python3) edwinaha@tcn1179.local.snellius.surf.nl:/scratch-shared/edwinaha/dynqual_waterbodies$ cdo showvar waterBodies5ArcMin_2010.nc
 #~ fracWaterInp waterBodyIds waterBodyTyp resMaxCapInp resSfAreaInp
#~ cdo    showname: Processed 5 variables [0.14s 80MB].

cdo selvar,fracWaterInp waterBodies5ArcMin_2010.nc fracWaterInp_waterBodies5ArcMin_2010.nc
cdo selvar,waterBodyIds waterBodies5ArcMin_2010.nc waterBodyIds_waterBodies5ArcMin_2010.nc
cdo selvar,waterBodyTyp waterBodies5ArcMin_2010.nc waterBodyTyp_waterBodies5ArcMin_2010.nc
cdo selvar,resMaxCapInp waterBodies5ArcMin_2010.nc resMaxCapInp_waterBodies5ArcMin_2010.nc
cdo selvar,resSfAreaInp waterBodies5ArcMin_2010.nc resSfAreaInp_waterBodies5ArcMin_2010.nc

gdal_translate -of PCRaster fracWaterInp_waterBodies5ArcMin_2010.nc fracWaterInp_waterBodies5ArcMin_2010.map
gdal_translate -of PCRaster waterBodyIds_waterBodies5ArcMin_2010.nc waterBodyIds_waterBodies5ArcMin_2010.map
gdal_translate -of PCRaster waterBodyTyp_waterBodies5ArcMin_2010.nc waterBodyTyp_waterBodies5ArcMin_2010.map
gdal_translate -of PCRaster resMaxCapInp_waterBodies5ArcMin_2010.nc resMaxCapInp_waterBodies5ArcMin_2010.map
gdal_translate -of PCRaster resSfAreaInp_waterBodies5ArcMin_2010.nc resSfAreaInp_waterBodies5ArcMin_2010.map

rm *.xml
