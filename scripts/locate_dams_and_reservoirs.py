#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil

import pcraster as pcr
import virtualOS as vos

# table/column input file - this should be sorted based on reservoir capacities (the largest the lower ids)
column_input_file = "../aha_hydropowers/existing_reservoirs_v2024-04-10_selected_no-headers.csv"

# ~ # - for testing with 74 and 85 (multiple pixel problem)
# ~ column_input_file = "../aha_hydropowers/existing_reservoirs_v2024-04-10_selected_no-headers-test74and85.csv"


# set the clone (your map grid system and extent)
clone_map_file = "../pcrglobwb_maps/lddsound_05min_version_20210330.map"
pcr.setclone(clone_map_file)

# ldd map (river network)
ldd_map_file = "../pcrglobwb_maps/lddsound_05min_version_20210330.map"
ldd_map = pcr.readmap(ldd_map_file)

# cell area (unit: m2)
cell_area_file = "../pcrglobwb_maps/cdo_gridarea_clone_global_05min_correct_lats.nc.map"
cell_area = pcr.readmap(cell_area_file)

# hydrolakes
hydrolakes_file = "../pcrglobwb_maps/areaIDs.map"
hydrolakes_ids = pcr.readmap(hydrolakes_file)

# calculate catchment area (in km2) based on pcrglobwb ldd
catchment_area_km2 = pcr.catchmenttotal(cell_area, ldd_map) / (1000.*1000.)

# calculate rough estimation of hydrolakes surface area based on pcrglobwb cell area
hydrolakes_pcrglobwb_area = pcr.areatotal(cell_area, hydrolakes_ids)

# convert table/column to a pcraster map of dam ids
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 1 -S -s ',' " + column_input_file + " dam_ids.map" 
print(cmd); os.system(cmd)
# - read dam ids as a variable
dam_ids = pcr.readmap("dam_ids.map") 

# convert table/column to a pcraster map of catchment areas based on AHA
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 2 -S -s ',' " + column_input_file + " aha_catchment_area_km2.map"
print(cmd); os.system(cmd)
# - read aha_catchment_area_km2 as a variable
aha_catchment_area_km2 = pcr.readmap("aha_catchment_area_km2.map") 

# convert table/column to a pcraster map of the latitude coordinates based on AHA
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 3 -S -s ',' " + column_input_file + " aha_latitudes.map"
print(cmd); os.system(cmd)
# - read as a variable
aha_latitudes = pcr.readmap("aha_latitudes.map") 

# convert table/column to a pcraster map of the latitude coordinates based on AHA
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 4 -S -s ',' " + column_input_file + " aha_longitudes.map"
print(cmd); os.system(cmd)
# - read as a variable
aha_longitudes = pcr.readmap("aha_longitudes.map") 

# convert table/column to a pcraster map of the reservoir surface area based on AHA (unit: km2)
cmd = "col2map --clone " + clone_map_file + " -M -x 4 -y 3 -v 6 -S -s ',' " + column_input_file + " aha_surface_area_km2.map"
print(cmd); os.system(cmd)
# - read as a variable
aha_surface_area_km2 = pcr.readmap("aha_surface_area_km2.map") 

# get the pcrglobwb catchment area for every dam id
dam_ids_pcrglobwb_catchment_area_km2 = pcr.ifthen(pcr.defined(dam_ids), catchment_area_km2)

# calculate the relative difference between two catchment areas
rel_dif_catchment_area = pcr.abs(aha_catchment_area_km2 - dam_ids_pcrglobwb_catchment_area_km2) / dam_ids_pcrglobwb_catchment_area_km2

# loop through all dams (from the largest to the smallest), if rel_dif_catchment_area > threshold, we have to reposition it
threshold = 0.1
number_of_dams = 131

# ~ # for testing
# ~ number_of_dams = 10
# ~ number_of_dams = 2


def get_distances_to_a_reference_latlon_coodinate_in_meter(lon_pcrmap, lat_pcrmap, lon_ref, lat_ref):
    pcr_distance_map = pcr.scalar(pcr.acos(pcr.sin(lat_pcrmap)*pcr.sin(lat_ref) + pcr.cos(lat_pcrmap)*pcr.cos(lat_ref)*pcr.cos(lon_ref - lon_pcrmap)))/360.*6371000. 
    return pcr_distance_map


for dam_id in range(1, number_of_dams + 1):
    
    print(dam_id)
    
    # make a point map of this dam
    this_dam_point = pcr.ifthen(dam_ids == dam_id, pcr.boolean(1.0))
    
    # evaluate relative difference for this particular dam
    rel_dif_catchment_area_this_dam = pcr.ifthen(this_dam_point, rel_dif_catchment_area)
    # - get its cell value
    rel_dif_catchment_area_this_dam_cell_value = pcr.cellvalue(pcr.mapmaximum(rel_dif_catchment_area_this_dam),1)[0]
    
    # correcting the location of dam
    if rel_dif_catchment_area_this_dam_cell_value < threshold:

        location_corrected_dam_id = pcr.ifthen(this_dam_point, pcr.nominal(dam_id))
        
    else:
    
        # expanding the point to its neighbours
        search_window = pcr.windowmajority(this_dam_point, 5./60. * 3.)
        # - note using the window_length = 2 to avoid 'too large' window size
        
        # get the pcrglobwb catchment area within this search_window
        catchment_area_within_search_window = pcr.ifthen(pcr.defined(search_window), catchment_area_km2)
        
        # compare the above to aha catchment_area
        aha_catchment_area_km2_this_dam = pcr.ifthen(this_dam_point, aha_catchment_area_km2)
        aha_catchment_area_km2_this_dam = pcr.windowmaximum(aha_catchment_area_km2_this_dam, 5./60 * 3.)

        # calculate the absolute difference
        difference_catch_area = pcr.abs(catchment_area_within_search_window - aha_catchment_area_km2_this_dam)
        
        location_corrected_dam_id = pcr.ifthen(difference_catch_area == pcr.mapminimum(difference_catch_area), pcr.nominal(dam_id))
        

        # make sure that there is only one pixel in location_corrected_dam_id; if this is NOT the case, we just find the closest one to the original coordinates 

        # - get the latitude and longitude coordinates based on AHA, and transfer them to single values maps (using mapmaximum)
        lat_aha_coordinate = pcr.mapmaximum(pcr.ifthen(dam_ids == dam_id, aha_latitudes))
        lon_aha_coordinate = pcr.mapmaximum(pcr.ifthen(dam_ids == dam_id, aha_longitudes))
        
        # - calculate the distance of every candidate point to the AHA coordinate
        distance_map = get_distances_to_a_reference_latlon_coodinate_in_meter(lon_pcrmap = pcr.xcoordinate(pcr.defined(location_corrected_dam_id)),\
                                                                              lat_pcrmap = pcr.ycoordinate(pcr.defined(location_corrected_dam_id)),\
                                                                              lon_ref = lon_aha_coordinate,\
                                                                              lat_ref = lat_aha_coordinate)
        area_order = pcr.areaorder(distance_map, location_corrected_dam_id)
        # - choose the one with shortest distance (order/rank = 1)
        location_corrected_dam_id = pcr.ifthen(area_order == 1, pcr.nominal(dam_id))
        
    if dam_id == 1:    
        all_location_corrected_dam_ids = location_corrected_dam_id
    else:
        all_location_corrected_dam_ids = pcr.cover(all_location_corrected_dam_ids, location_corrected_dam_id)
    
    # get the reservoir surface area based on AHA (unit: m2)
    aha_surface_area_m2 = pcr.ifthen(this_dam_point, aha_surface_area_km2) * 1000.*1000.
    aha_surface_area_m2_this_dam_cell_value = pcr.cellvalue(pcr.mapmaximum(aha_surface_area),1)[0]
    
    # get the cell area for this dam (unit: m2)
    cell_area_m2_this_dam_cell_value = pcr.cellvalue(pcr.mapmaximum(pcr.ifthen(this_dam_point, cell_area)))
    
    # obtaining the reservoir extent (including estimating surface area)
    if aha_surface_area_m2_this_dam_cell_value < cell_area_m2_this_dam_cell_value:
        
        # this means the reservoir extent will be within a cell
        reservoir_surface_area           = pcr.scalar(aha_surface_area_m2_this_dam_cell_value)
        reservoir_surface_area_per_cell  = reservoir_surface_area
        
        # fraction of surface water within a cell
        reservoir_fraction_water = reservoir_surface_area_per_cell/cell_area

    else:

        # find the hydrolakes ID and its extent
        hydrolakes_id_selected = ???

        # expanding the point to its neighbours
        search_window = pcr.windowmajority(this_dam_point, 5./60. * 3.)
        # - note using the window_length = 2 to avoid 'too large' window size
        
        # find the hydrolakes ids within this search window 
        hydrolakes_ids_within_search_window = pcr.ifthen(pcr.defined(search_window), hydro_lakes)
        
        # check whether there are more than one hydrolakes_ids_within_search_window
        number_of_hydrolakes_ids_within_search_window = ???

        if number_of_hydrolakes_ids_within_search_window > 1:
            
            # find the one that has the most similar surface area to the estimate on hydrolakes_pcrglobwb_area
            absolute_difference_surface_area = pcr.abs(aha_surface_area_km2 - hydrolakes_pcrglobwb_area)

            area_order = pcr.areaorder(absolute_difference_surface_area, location_corrected_dam_id)
            
            # - choose the one with shortest distance (order/rank = 1)
            location_corrected_dam_id = pcr.ifthen(area_order == 1, pcr.nominal(dam_id))

            
            pass
        
        else:

            hydrolakes_id_selected = hydrolakes_ids_within_search_window

        # identify the number of cells based on hydrolakes
        number_of_cells_according_to_hydrolakes

        # obtaining reservoir surface area
        reservoir_surface_area          = number_of_cells_according_to_hydrolakes * reservoir_surface_area_per_cell
        reservoir_surface_area_per_cell = reservoir_surface_area / number_of_cells_according_to_hydrolakes
        
        # fraction of surface water within a cell
        reservoir_fraction_water = reservoir_surface_area_per_cell / cell_area
        
    
    
    # compare aha_surface_area with individual cell_area
    
    
    

    
        

# save the all_location_corrected_dam_ids to a pcraster map
pcr.report(all_location_corrected_dam_ids, "corrected_dam_ids.map")

# obtain the catchment areas of all_location_corrected_dam_ids
corrected_dam_catchment_area_km2 = pcr.ifthen(pcr.defined(all_location_corrected_dam_ids), catchment_area_km2)
pcr.report(corrected_dam_catchment_area_km2, "corrected_dam_catchment_area_km2.map")

# obtain a table/column format for all_location_corrected_dam_ids and corrected_dam_catchment_area_km2
cmd = "map2col corrected_dam_ids.map corrected_dam_catchment_area_km2.map corrected_dams.txt"
print(cmd); os.system(cmd)        


