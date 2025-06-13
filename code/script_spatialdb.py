#---
# puma grass spatial database
# building the database
# crs: Albers Equal Area, datum SIRGAS 2000 (GRS80)
#
# bernardo niebuhr
# 2020-09-17
#---

python

# modules
# import modules
import os
import grass.script as grass
from grass.pygrass.modules.shortcuts import general as g
from grass.pygrass.modules.shortcuts import vector as v
from grass.pygrass.modules.shortcuts import raster as r

#---------------------------------------
# Setup

# Load maps in mapset PERMANENT

set GRASS_COMPRESSOR=ZLIB

#---------------------------------------
# Load data

# folder
root = r"D:\bernardo\00_academico\07_projetos\04_oncas_on_the_move\_Pardas_Tiete\products\ms_Niebuhr_puma_dispersal_HR\spatial_data"
os.chdir(root)

#-------------------
# Load vectors

#-----
# buffer for SP

# folder
pa = root+"/br_states"
os.chdir(pa)

v.in_ogr(input = "SP_buff100km_albers_sirgas.gpkg", output = "SP_buff100km_albers_sirgas",
	overwrite = True)

#-----
# roads

# folder
pa = root+"/roads"
os.chdir(pa)

v.in_ogr(input = "roads_dnit2014_forestgis_editado_albers_sirgas.gpkg", 
    output = "roads_dnit2014_forestgis_editado_albers_sirgas")

grass.run_command("v.import", input = "roads_dnit2014_forestgis_editado_albers_sirgas.gpkg", 
    output = "roads_dnit2014_forestgis_editado_albers_sad69",
    overwrite = True)

#-------------------
# Load rasters

# folder
pa = root+"/use"
os.chdir(pa)

#-----
# sugarcane CanaSat

# folder of the other location
dbase = r"D:\bernardo\00_academico\07_projetos\04_oncas_on_the_move\_Pardas_Tiete\grassDB"

# reproject from other location
g.region(raster = "FBDS_30m_map_all_states_albers_sirgas_compressed", flags = "ap")
r.proj(location = "location5m", mapset = "landscapemetrics", input = "canasat_2013_albers_sirgas2000",
    output = "canasat_2013_albers_sad69", dbase = dbase)

# export it to have it saved
r.out_gdal(input = "canasat_2013_albers_sad69", output = "canasat_2013_albers_sad69.tif", 
    createopt = "TFW=YES,COMPRESS=DEFLATE")

# reimport
r.in_gdal(input = "canasat_2013_albers_sirgas2000_ok.tif", output = "canasat_2013_albers_sirgas")

# map to align
map_to_align = "canasat_2013_albers_sirgas"

#-----
# land use FBDS
r.in_gdal(input = "FBDS_30m_map_all_states_albers_sirgas_compressed.tif", output = "FBDS_30m_map_all_states_albers_sirgas_compressed")

# align map
g.region(raster = "FBDS_30m_map_all_states_albers_sirgas_compressed", align = map_to_align, flags = "p")
r.mapcalc("FBDS_30m_map_all_states = FBDS_30m_map_all_states_albers_sirgas_compressed")

#-----
# pasture Lapig

# import pasture map for the Atlantic Forest
r.in_gdal(input = "mata_atlantica_pasture_2013_albers_sirgas.tif", 
    output = "mata_atlantica_pasture_2013_albers_sirgas", overwrite = True)

# read vector
v.in_ogr(input = "mata_atlantica_pasture_2013_albers_sirgas.gpkg", output = "mata_atlantica_pasture_2013_albers_sirgas")

# region
g.region(vector = "mata_atlantica_pasture_2013_albers_sirgas", 
    align = map_to_align, flags = "p")
# rasterize
v.to_rast(input = "mata_atlantica_pasture_2013_albers_sirgas", output = "mata_atlantica_pasture_2013_albers_sirgas_rast",
    use = "val", value = 8, overwrite = True)

# import pasture map for the Cerrado
r.in_gdal(input = "cerrado_pasture_2013_albers_sirgas.tif", 
    output = "cerrado_pasture_2015_albers_sirgas", overwrite = True)

# read vector
v.in_ogr(input = "cerrado_pasture_2013_albers_sirgas.gpkg", output = "cerrado_pasture_2013_albers_sirgas")

# region
g.region(vector = "cerrado_pasture_2013_albers_sirgas", 
    align = map_to_align, flags = "p")
# rasterize
v.to_rast(input = "cerrado_pasture_2013_albers_sirgas", output = "cerrado_pasture_2013_albers_sirgas_rast",
    use = "val", value = 8, overwrite = True)

#-----
# mosaic pasture maps for Mata Atlantica and Cerrado
pasture_maps = ("mata_atlantica_pasture_2013_albers_sirgas_rast", 
    "cerrado_pasture_2013_albers_sirgas_rast")
g.region(raster = pasture_maps, align = map_to_align, flags = "p")

r.patch(input = pasture_maps, output = "pasture_2015_lapig_cerrado_mata_atlantica_rast", overwrite = True)

#-----
# water from FBDS - until now, only within SP
# and with big water bodies mixed with large and small water courses

# import water map from FBDS for SP - 5m
r.in_gdal(input = "FBDS_hidro_all_SP_sirgas_albers.tif", 
    output = "FBDS_hidro_all_SP_sirgas_albers_5m", overwrite = True)

# import water map from FBDS for SP - 30m
r.in_gdal(input = "FBDS_hidro_all_SP_sirgas_albers_30m.tif", 
    output = "FBDS_hidro_all_SP_sirgas_albers_30m", overwrite = True)

#---------------------------------------
# Process data

#-----
# land use

# get pasture that is within antropic class from FBDS and is not sugarcane from CanaSat
g.region(raster = "pasture_2015_lapig_cerrado_mata_atlantica_rast", flags = "p")
expression = "pasture_lapig_2015_antropic_nosugarcane = if(FBDS_30m_map_all_states_albers_sirgas_compressed == 2 && isnull(canasat_2013_albers_sirgas), \
pasture_2015_lapig_cerrado_mata_atlantica_rast, null())"
r.mapcalc(expression, overwrite = True)

# get sugarcane that is within antropic class from FBDS
g.region(raster = "canasat_2013_albers_sirgas", flags = "p")
expression = "canasat_2013_albers_sirgas2000_cut_FBDS_antropic = if(FBDS_30m_map_all_states_albers_sirgas_compressed == 2 && !isnull(canasat_2013_albers_sirgas), \
7, null())"
r.mapcalc(expression, overwrite = True)

# merge all land uses
landuse_maps = ("canasat_2013_albers_sirgas2000_cut_FBDS_antropic", "pasture_lapig_2015_antropic_nosugarcane", 
    "FBDS_30m_map_all_states_albers_sirgas_compressed")
g.region(raster = landuse_maps, align = map_to_align, flags = "p")

r.patch(input = landuse_maps, output = "FBDS_30m_map_all_states_sugarcane7_pasture8", overwrite = True)

#-----
# export the whole map

# folder
pa = root+"/use"
os.chdir(pa)

# region
g.region(raster = "FBDS_30m_map_all_states_sugarcane7_pasture8", flags = "p")

# export
r.out_gdal(input = "FBDS_30m_map_all_states_sugarcane7_pasture8", output = "FBDS_30m_map_all_states_sugarcane7_pasture8.tif",
    createopt = "TFW=YES,COMPRESS=DEFLATE")

# export only for the study area

# region
g.region(vector = "SP_buff100km_albers_sirgas", align = map_to_align, flags = "ap")
# mask
r.mask(vector = "SP_buff100km_albers_sirgas")
# export
r.out_gdal(input = "FBDS_30m_map_all_states_sugarcane7_pasture8", output = "FBDS_30m_map_all_states_sugarcane7_pasture8_SP.tif",
    createopt = "TFW=YES,COMPRESS=DEFLATE")
# remove mask
r.mask(flags = "r")

#-------------
# Calculate distances and densities

# mask for SP buffer
# region
g.region(vector = "SP_buff100km_albers_sirgas", align = map_to_align, flags = "ap")
# mask
r.mask(vector = "SP_buff100km_albers_sirgas", overwrite = True)

#--------
# cities
# get cities from FBDS as 1/0
expression = "cities_FBDS = if(FBDS_30m_map_all_states_sugarcane7_pasture8 == 3, 1, 0)"
r.mapcalc(expression, overwrite = True)

# get cities from FBDS as 1/null
expression = "cities_FBDS_null = if(FBDS_30m_map_all_states_sugarcane7_pasture8 == 3, 1, null())"
r.mapcalc(expression, overwrite = True)

# distance from urban areas
r.grow_distance(input = "cities_FBDS_null", distance = "dist_cities_FBDS_30m")

# density of urban areas in a 10km radius
size_m = 10000 # in meters
res = 30 # resolution in meters
size_pixels = round(size_m/res)
r.neighbors(input = "cities_FBDS", output = "urban_density_n_pixels_per10km2", 
    method = "average", size = size_pixels, flags = "c")
    
# density of urban areas using multiple radii
sizes_m = [100, 250, 500, 1000, 2500, 5000, 10000] # in meters
for i in sizes_m:
    r.resamp_filter(input = "cities_FBDS", output = "urban_density_scale_"+str(i), 
        filter = "bartlett", radius = i)

#--------
# roads

# rasterize
v.to_rast(input = "roads_dnit2014_forestgis_editado_albers_sirgas", output = "roads_dnit2014",
    use = "val", value = 1, overwrite = True)

# roads as 1/0
r.mapcalc("roads_dnit2014_1_0 = roads_dnit2014")
r.null(map = "roads_dnit2014_1_0", null = 0)

# distance from roads
r.grow_distance(input = "roads_dnit2014", distance = "dist_roads_dnit2014")

# density of road pixels areas in a 10km radius
size_m = 10000 # in meters
res = 30 # resolution in meters
size_pixels = round(size_m/res)
r.neighbors(input = "roads_dnit2014_1_0", output = "road_density_n_pixels_per10km2", 
    method = "sum", size = size_pixels, flags = "c")

# density of roads in km of roads per 100km2
expression = "road_density_km_per100km2 = 30.0 / 1000 * road_density_n_pixels_per10km2"
r.mapcalc(expression, overwrite = True)

# density of roads using multiple radii
sizes_m = [100, 250, 500, 1000, 2500, 5000, 10000] # in meters
for i in sizes_m:
    r.resamp_filter(input = "roads_dnit2014_1_0", output = "road_density_scale_"+str(i), 
        filter = "bartlett", radius = i)

#--------
# rivers/water
g.region(raster = "FBDS_hidro_all_SP_sirgas_albers_5m", flags = "ap")

# distance from water
r.grow_distance(input = "FBDS_hidro_all_SP_sirgas_albers_5m", distance = "dist_hidro_5m")

# mask for SP buffer
# region
g.region(vector = "SP_buff100km_albers_sirgas", align = map_to_align, flags = "ap")
# mask
r.mask(vector = "SP_buff100km_albers_sirgas", overwrite = True)

# rivers as 1/0
r.mapcalc("FBDS_hidro_1_0 = FBDS_hidro_all_SP_sirgas_albers_5m")
r.null(map = "FBDS_hidro_1_0", null = 0)

# density of water using multiple radii
sizes_m = [100, 250, 500, 1000, 2500, 5000, 10000] # in meters
for i in sizes_m:
    r.resamp_filter(input = "FBDS_hidro_1_0", output = "water_density_scale_"+str(i), 
        filter = "bartlett", radius = i)

#--------
# forest

# binary forest 1/0
expression = "forest_1_0 = if(FBDS_30m_map_all_states_sugarcane7_pasture8 == 4, 1, 0)"
r.mapcalc(expression, overwrite = True)

# forest as 1/null
expression = "forest_1_null = if(FBDS_30m_map_all_states_sugarcane7_pasture8 == 4, 1, null())"
r.mapcalc(expression, overwrite = True)

# distance from forest
r.grow_distance(input = "forest_1_null", distance = "dist_forest_FBDS_30m")

# density of forest using multiple radii
sizes_m = [100, 250, 500, 1000, 2500, 5000, 10000] # in meters
for i in sizes_m:
    r.resamp_filter(input = "forest_1_0", output = "forest_density_scale_"+str(i), 
        filter = "bartlett", radius = i)
        
#--------
# pasture

# binary pasture 1/0
expression = "pasture_1_0 = if(pasture_lapig_2015_antropic_nosugarcane == 8, 1, 0)"
r.mapcalc(expression, overwrite = True)
r.null(map = "pasture_1_0", null = 0)

# pasture as 1/null
expression = "pasture_1_null = if(pasture_lapig_2015_antropic_nosugarcane == 8, 1, null())"
r.mapcalc(expression, overwrite = True)

# distance from pasture
r.grow_distance(input = "pasture_1_null", distance = "dist_pasture_FBDS_30m")

# density of pasture using multiple radii
sizes_m = [100, 250, 500, 1000, 2500, 5000, 10000] # in meters
for i in sizes_m:
    r.resamp_filter(input = "pasture_1_0", output = "pasture_density_scale_"+str(i), 
        filter = "bartlett", radius = i, overwrite = True)

#--------
# sugarcane

# binary sugarcane 1/0
expression = "sugarcane_1_0 = if(canasat_2013_albers_sirgas2000_cut_FBDS_antropic == 7, 1, 0)"
r.mapcalc(expression, overwrite = True)
r.null(map = "sugarcane_1_0", null = 0)

# sugarcane as 1/null
expression = "sugarcane_1_null = if(canasat_2013_albers_sirgas2000_cut_FBDS_antropic == 7, 1, null())"
r.mapcalc(expression, overwrite = True)

# distance from sugarcane
r.grow_distance(input = "sugarcane_1_null", distance = "dist_sugarcane_FBDS_30m")

# density of sugarcane using multiple radii
sizes_m = [100, 250, 500, 1000, 2500, 5000, 10000] # in meters
for i in sizes_m:
    r.resamp_filter(input = "sugarcane_1_0", output = "sugarcane_density_scale_"+str(i), 
        filter = "bartlett", radius = i, overwrite = True)


#--------
# forestry

# binary forestry 1/0
expression = "forestry_1_0 = if(FBDS_30m_map_all_states_sugarcane7_pasture8 == 6, 1, 0)"
r.mapcalc(expression, overwrite = True)

# forestry as 1/null
expression = "forestry_1_null = if(FBDS_30m_map_all_states_sugarcane7_pasture8 == 6, 1, null())"
r.mapcalc(expression, overwrite = True)

# distance from forestry
r.grow_distance(input = "forestry_1_null", distance = "dist_forestry_FBDS_30m")

# density of forestry using multiple radii
sizes_m = [1000, 2500, 5000, 10000]#100, 250, 500, 1000, 2500, 5000, 10000] # in meters
for i in sizes_m:
    r.resamp_filter(input = "forestry_1_0", output = "forestry_density_scale_"+str(i), 
        filter = "bartlett", radius = i, overwrite = True)

#---------
# get elevation


#-----------------
# export maps

# folder
pa = root+"/analysis"
os.chdir(pa)

# maps to export
maps_to_export = ["FBDS_30m_map_all_states_sugarcane7_pasture8", 
"road_density_km_per100km2", "dist_roads_dnit2014", 
"urban_density_n_pixels_per10km2", "dist_cities_FBDS_30m", "dist_hidro_5m"]

# new simple names
names = ["land_use_FBDS_8_30m", "road_density_km_per100km2", "dist_road_m",
"urban_density_n_pixels_per100km2", "dist_urban_m", "dist_hidro_m"]

map_to_align = "FBDS_30m_map_all_states_sugarcane7_pasture8"

# export
for i in range(len(maps_to_export)):
    if "cities" in maps_to_export[i]:
        # region
        g.region(vector = "SP_buff100km_albers_sirgas", res = 30, 
            align = map_to_align, flags = "ap")
        # export
        r.out_gdal(input = maps_to_export[i], output = names[i]+".tif", 
            createopt = "TFW=YES,COMPRESS=DEFLATE", overwrite = True)

