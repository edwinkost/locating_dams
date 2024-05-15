#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import shutil

import pcraster as pcr
import virtualOS as vos

# ~ # table/column input file - this should be sorted based on reservoir capacities (the largest the lower ids)
# ~ column_input_file = "../aha_hydropowers/existing_reservoirs_v2024-04-10_selected_no-headers.csv"

# column file with dam locations that have been corrected (adjusted to PCR-GLOBWB) 
column_input_file = "../aha_hydropowers/corrected_dams_2024-04-26_no-headers.csv"

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
# ~ hydrolakes_file = "../pcrglobwb_maps/areaIDs.map"
# - using pcrglobwb lake and reservoir extent
hydrolakes_file = "../pcrglobwb_maps/waterBodyIds_waterBodies5ArcMin_2010.map"
hydrolakes_ids = pcr.nominal(pcr.readmap(hydrolakes_file))
hydrolakes_ids = pcr.ifthen(pcr.scalar(hydrolakes_ids) > 0, hydrolakes_ids)
# ~ pcr.aguila(hydrolakes_ids)

# fraction water
hydrolakes_fraction_water_file = "../pcrglobwb_maps/fracWaterInp_waterBodies5ArcMin_2010.map" 
hydrolakes_fraction_water = pcr.readmap(hydrolakes_fraction_water_file)

# calculate catchment area (in km2) based on pcrglobwb ldd
catchment_area_km2 = pcr.catchmenttotal(cell_area, ldd_map) / (1000.*1000.)

# calculate rough estimation of hydrolakes surface area based on pcrglobwb cell area (unit: m2)
hydrolakes_pcrglobwb_area_m2 = pcr.areatotal(cell_area, hydrolakes_ids)

# convert table/column to a pcraster map of dam ids
cmd = "col2map --clone " + clone_map_file + " -M -x 1 -y 2 -v 3 -S -s ',' " + column_input_file + " dam_ids.map" 
print(cmd); os.system(cmd)
# - read dam ids as a variable
dam_ids = pcr.readmap("dam_ids.map") 

# convert table/column to a pcraster map of catchment areas based on AHA
cmd = "col2map --clone " + clone_map_file + " -M -x 1 -y 2 -v 4 -S -s ',' " + column_input_file + " aha_catchment_area_km2.map"
print(cmd); os.system(cmd)
# - read aha_catchment_area_km2 as a variable
aha_catchment_area_km2 = pcr.readmap("aha_catchment_area_km2.map") 

# convert table/column to a pcraster map of the latitude coordinates based on AHA
cmd = "col2map --clone " + clone_map_file + " -M -x 1 -y 2 -v 2 -S -s ',' " + column_input_file + " aha_latitudes.map"
print(cmd); os.system(cmd)
# - read as a variable
aha_latitudes = pcr.readmap("aha_latitudes.map") 

# convert table/column to a pcraster map of the latitude coordinates based on AHA
cmd = "col2map --clone " + clone_map_file + " -M -x 1 -y 2 -v 1 -S -s ',' " + column_input_file + " aha_longitudes.map"
print(cmd); os.system(cmd)
# - read as a variable
aha_longitudes = pcr.readmap("aha_longitudes.map") 

# convert table/column to a pcraster map of the reservoir surface area based on AHA (unit: km2)
cmd = "col2map --clone " + clone_map_file + " -M -x 1 -y 2 -v 6 -S -s ',' " + column_input_file + " aha_surface_area_km2.map"
print(cmd); os.system(cmd)
# - read as a variable
aha_surface_area_km2 = pcr.readmap("aha_surface_area_km2.map") 

# get the pcrglobwb catchment area for every dam id
dam_ids_pcrglobwb_catchment_area_km2 = pcr.ifthen(pcr.defined(dam_ids), catchment_area_km2)

# calculate the relative difference between two catchment areas
rel_dif_catchment_area = pcr.abs(aha_catchment_area_km2 - dam_ids_pcrglobwb_catchment_area_km2) / dam_ids_pcrglobwb_catchment_area_km2

# loop through all dams
number_of_dams = 131

number_of_dams = 10

number_of_dams = 2


for dam_id in range(1, number_of_dams + 1):
    
    print(dam_id)
    
    # make a point map of this dam
    this_dam_point = pcr.ifthen(dam_ids == dam_id, pcr.boolean(1.0))
   
    # get the reservoir surface area based on AHA (unit: m2)
    aha_surface_area_m2 = pcr.ifthen(this_dam_point, aha_surface_area_km2) * 1000.*1000.
    aha_surface_area_m2_this_dam_cell_value = pcr.cellvalue(pcr.mapmaximum(aha_surface_area_m2),1)[0]
    
    # get the cell area for this dam (unit: m2)
    cell_area_m2_this_dam_cell_value = pcr.cellvalue(pcr.mapmaximum(pcr.ifthen(this_dam_point, cell_area)),1)[0]
    
    # obtaining the reservoir extent (including estimating surface area)
    if aha_surface_area_m2_this_dam_cell_value < cell_area_m2_this_dam_cell_value:
        
        # this means the reservoir extent will be within a cell
        reservoir_surface_area           = pcr.ifthen(this_dam_point, pcr.spatial(pcr.scalar(aha_surface_area_m2_this_dam_cell_value)))
        reservoir_surface_area_per_cell  = reservoir_surface_area
        
        # fraction of surface water within a cell
        reservoir_fraction_water = reservoir_surface_area_per_cell/cell_area

        # reservoir extent
        reservoir_extent = this_dam_point

    else:

        # this is for the case if reservoirs covering multiple cells
        
        # make a search window to find the nearest hydrolakes - expanding the point to its neighbours
        search_window = pcr.windowmajority(this_dam_point, 5./60. * 3.)
        # - note using the window_length = 2 to avoid 'too large' window size
        
        # find the hydrolakes ids within this search window 
        hydrolakes_ids_within_search_window = pcr.ifthen(pcr.defined(search_window), hydrolakes_ids)
        hydrolakes_ids_within_search_window = pcr.areamajority(hydrolakes_ids_within_search_window, hydrolakes_ids)
        
        # ~ pcr.aguila(hydrolakes_ids_within_search_window)
        
        # identify the outlets - based on pcrglobwb catchment areas
        hydrolakes_ids_within_search_window_catchment_areas_order = pcr.areaorder( pcr.ifthen(pcr.defined(hydrolakes_ids_within_search_window), catchment_area_km2 * -1.0), hydrolakes_ids)
        hydrolakes_ids_within_search_window_outlets = pcr.ifthen(hydrolakes_ids_within_search_window_catchment_areas_order == 1, pcr.nominal(hydrolakes_ids_within_search_window))
        
        # check whether there are more than one hydrolakes_ids_within_search_window
        number_of_hydrolakes_ids_within_search_window = pcr.cellvalue(pcr.mapmaximum(pcr.scalar(pcr.clump(hydrolakes_ids_within_search_window_outlets))),1)[0]

        print(number_of_hydrolakes_ids_within_search_window)
        
        if number_of_hydrolakes_ids_within_search_window > 0:
            
            hydrolakes_ids_within_search_window_area_m2 = pcr.areatotal(cell_area * hydrolakes_fraction_water, hydrolakes_ids_within_search_window)
            
            pcr.aguila(hydrolakes_ids_within_search_window_area_m2)
            
            # find the one that has the most similar surface area to the estimate on hydrolakes_pcrglobwb_area
            absolute_difference_surface_area = pcr.abs(hydrolakes_ids_within_search_window_area_m2 - aha_surface_area_m2_this_dam_cell_value)/aha_surface_area_m2_this_dam_cell_value
            
            print(aha_surface_area_m2_this_dam_cell_value)
            
            pcr.aguila(absolute_difference_surface_area)
            
            absolute_difference_surface_area = pcr.ifthen(pcr.defined(hydrolakes_ids_within_search_window_outlets), absolute_difference_surface_area)
            
            class_map_boolean = pcr.defined(absolute_difference_surface_area)
            class_map_nominal = pcr.ifthen(class_map_boolean, pcr.nominal(1.0))
            
            area_order = pcr.areaorder(absolute_difference_surface_area, class_map_nominal)
            
            # - choose the one with order/rank = 1 (minimum difference in surface area)
            hydrolakes_id_selected = pcr.mapmaximum(pcr.scalar(pcr.ifthen(area_order == 1, pcr.nominal(hydrolakes_ids_within_search_window))))
            
            # - the chosen hydro lakes id and assign in to the dam point
            hydrolakes_id_for_this_dam_point = pcr.ifthen(this_dam_point, pcr.nominal(hydrolakes_id_selected))
            
            # make the sub catchment until the dam point and its ldd
            sub_catchment = pcr.subcatchment(ldd_map, this_dam_point)
            ldd_above_the_dam_point = pcr.lddmask(ldd_map, sub_catchment)
            
            # define the extent based on the hydrolakes
            hydrolakes_id_selected_extent = pcr.ifthen(pcr.scalar(hydrolakes_ids) == hydrolakes_id_selected, pcr.boolean(1.0))
            
            # - expand the extent until the dam point until the entire reservoirs
            hydrolakes_id_selected_extent = pcr.path(ldd_above_the_dam_point, hydrolakes_id_selected_extent)
            
            # - provide an id to this reservoir
            hydrolakes_id_for_this_dam_point = pcr.ifthen(hydrolakes_id_selected_extent, pcr.spatial(pcr.scalar(dam_id)))
            
            # identify the number of cells based on the hydrolakes
            number_of_cells_according_to_hydrolakes = pcr.areatotal(pcr.spatial(pcr.scalar(1.0)), pcr.nominal(hydrolakes_id_for_this_dam_point))
		    
            # obtaining reservoir surface area
            reservoir_surface_area_per_cell = aha_surface_area_m2_this_dam_cell_value / number_of_cells_according_to_hydrolakes
            
            # fraction of surface water within a cell
            reservoir_fraction_water = reservoir_surface_area_per_cell / cell_area

            # reservoir extent
            reservoir_extent = pcr.defined(hydrolakes_id_for_this_dam_point)
            reservoir_extent = pcr.ifthen(reservoir_extent, reservoir_extent)

            reservoir_surface_area             = pcr.ifthen(reservoir_extent, pcr.spatial(pcr.scalar(aha_surface_area_m2_this_dam_cell_value)))
            reservoir_surface_area_per_cell    = pcr.ifthen(reservoir_extent, reservoir_surface_area_per_cell)
            reservoir_fraction_water           = pcr.ifthen(reservoir_extent, reservoir_fraction_water)
            
            # ~ pcr.aguila(reservoir_extent)

        else:

            # if not identified in the hydrolakes

            # start with assuming everything from zero
            reservoir_surface_area             = 0.0
            reservoir_surface_area_per_cell    = 0.0
            reservoir_fraction_water           = 0.0
            number_of_cells_for_this_reservoir = 0.0

            # assign the current cell as the reservoir 
            reservoir_extent                   = this_dam_point
            reservoir_surface_area             = pcr.ifthen(reservoir_extent, cell_area)
            reservoir_surface_area_per_cell    = pcr.ifthen(reservoir_extent, cell_area)
            reservoir_fraction_water           = pcr.ifthen(reservoir_extent, pcr.spatial(pcr.scalar(1.0)))
            number_of_cells_for_this_reservoir = 1.0
            
            # calculate the remaining surface area that needs to be covered
            remaining_area = aha_surface_area_m2_this_dam_cell_value - pcr.cellvalue(pcr.maptotal(reservoir_surface_area_per_cell),1)[0]
            
            while remaining_area > 0:
            
                # identify the upstream cells
                upstream_cells_of_reservoirs = pcr.boolean(pcr.downstream(ldd_map, pcr.cover(pcr.scalar(reservoir_extent), 0.0)))

                # number of upstream cells
                num_of_upstream_cells = pcr.maptotal(pcr.scalar(upstream_cells_of_reservoirs))
                
                # reservoir surface area at the upstream cells
                additional_reservoir_surface_area_per_cell = pcr.ifthen(upstream_cells_of_reservoirs, pcr.min(remaining_area / num_of_upstream_cells, cell_area))
                
                # the updated reservoir surface area per cell
                reservoir_surface_area_per_cell = pcr.cover(reservoir_surface_area_per_cell, additional_reservoir_surface_area_per_cell)
                
                # calculate the remaining surface area that needs to be covered
                remaining_area = aha_surface_area_m2_this_dam_cell_value - pcr.cellvalue(pcr.maptotal(reservoir_surface_area_per_cell))
            
            # update the reservoir with upstream cells
            reservoir_extent         = pcr.cover(reservoir_extent, pcr.ifthen(reservoir_surface_area_per_cell > 0.0, pcr.boolean(1.0)))
            reservoir_surface_area   = pcr.ifthen(reservoir_extent, pcr.maptotal(reservoir_surface_area_per_cell))
            reservoir_fraction_water = reservoir_surface_area_per_cell / cell_area
        

    # we summarize all reservoirs into single variables 
    if dam_id == 1:    
        all_reservoir_extent_ids     = pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.spatial(pcr.scalar(dam_id))))
        all_reservoir_surface_area   = pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(reservoir_surface_area)))
        all_reservoir_fraction_water = pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(reservoir_fraction_water))) 
    else:
        all_reservoir_extent_ids     = pcr.cover(all_reservoir_extent_ids,     pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.spatial(pcr.scalar(dam_id)))))
        all_reservoir_surface_area   = pcr.cover(all_reservoir_surface_area,   pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(reservoir_surface_area))))
        all_reservoir_fraction_water = pcr.cover(all_reservoir_fraction_water, pcr.ifthen(reservoir_extent, pcr.mapmaximum(pcr.scalar(reservoir_fraction_water)))) 


# ~ # save all_location_corrected_dam_ids to a pcraster map - this will be the location of all dams/outlets
# ~ pcr.report(all_location_corrected_dam_ids, "corrected_dam_ids.map")

# ~ # obtain the catchment areas of all_location_corrected_dam_ids
# ~ corrected_dam_catchment_area_km2 = pcr.ifthen(pcr.defined(all_location_corrected_dam_ids), catchment_area_km2)
# ~ pcr.report(corrected_dam_catchment_area_km2, "corrected_dam_catchment_area_km2.map")

# ~ # obtain a table/column format for all_location_corrected_dam_ids and corrected_dam_catchment_area_km2
# ~ cmd = "map2col corrected_dam_ids.map corrected_dam_catchment_area_km2.map corrected_dams.txt"
# ~ print(cmd); os.system(cmd)        

all_reservoir_extent_ids_masked     = pcr.ifthen(pcr.defined(ldd_map), pcr.cover(pcr.nominal(all_reservoir_extent_ids), pcr.nominal(0.0)))
all_reservoir_surface_area_masked   = pcr.ifthen(pcr.defined(ldd_map), pcr.cover(pcr.scalar(all_reservoir_surface_area), 0.0))
all_reservoir_fraction_water_masked = pcr.ifthen(pcr.defined(ldd_map), pcr.cover(pcr.scalar(all_reservoir_fraction_water), 0.0))

# save also the following variables for pcrglobwb input:
pcr.report(all_reservoir_extent_ids_masked,     "reservoir_extent_ids.map")
pcr.report(all_reservoir_surface_area_masked,   "reservoir_surface_area_ids.map")
pcr.report(all_reservoir_fraction_water_masked, "reservoir_fraction_water_ids.map")
