import os
import sys
import shutil

import pcraster as pcr
import virtualOS as vos

# pcrglobwb/dynqual files
# ~ (pcrglobwb_python3) edwinaha@tcn1179.local.snellius.surf.nl:/scratch-shared/edwinaha/dynqual_waterbodies$ ls -lah *.map
# ~ -rw-r--r--. 1 edwinaha edwinaha 36M May 21 13:15 fracWaterInp_waterBodies5ArcMin_2010.map
# ~ -rw-r--r--. 1 edwinaha edwinaha 36M May 21 13:15 resMaxCapInp_waterBodies5ArcMin_2010.map
# ~ -rw-r--r--. 1 edwinaha edwinaha 36M May 21 13:15 resSfAreaInp_waterBodies5ArcMin_2010.map
# ~ -rw-r--r--. 1 edwinaha edwinaha 36M May 21 13:15 waterBodyIds_waterBodies5ArcMin_2010.map
# ~ -rw-r--r--. 1 edwinaha edwinaha 36M May 21 13:15 waterBodyTyp_waterBodies5ArcMin_2010.map

# AHA files
# ~ (pcrglobwb_python3) edwinaha@tcn1179.local.snellius.surf.nl:/home/edwinaha/github/edwinkost/locating_dams/scripts$ ls -lah existing*
# ~ -rw-r--r--. 1 edwinaha edwinaha 18K May 21 12:55 existing_locate_dams_and_reservoirs_2nd_part.py
# ~ -rw-r--r--. 1 edwinaha edwinaha 36M May 21 13:06 existing_reservoir_capacity_ids_million_m3.map
# ~ -rw-r--r--. 1 edwinaha edwinaha 36M May 21 13:06 existing_reservoir_extent_ids.map
# ~ -rw-r--r--. 1 edwinaha edwinaha 36M May 21 13:06 existing_reservoir_fraction_water_ids.map
# ~ -rw-r--r--. 1 edwinaha edwinaha 36M May 21 13:06 existing_reservoir_surface_area_ids_m2.map

# IDs based on pcrglowbb
pgb_ids_file = "/scratch-shared/edwinaha/dynqual_waterbodies/waterBodyIds_waterBodies5ArcMin_2010.map"
pgb_ids = pcr.nominal(pcr.readmap(pgb_ids_file))
pgb_ids = pcr.ifthen(pcr.scalar(pgb_ids) > 0.0, pgb_ids)

# IDs based on AHA
aha_ids_file = "/home/edwinaha/github/edwinkost/locating_dams/scripts/existing_reservoir_extent_ids.map"
aha_ids = pcr.nominal(pcr.readmap(aha_ids_file))
aha_ids = pcr.ifthen(pcr.scalar(aha_ids) > 0.0, aha_ids)
aha_ids_new_ids_scalar = pcr.mapmaximum(pcr.scalar(pgb_ids)) * 10.0 + pcr.scalar(aha_ids)
aha_ids = pcr.nominal(aha_ids_new_ids_scalar)

# exclude pgb_ids that are in aha_ids
pgb_ids_ids_exclusion = pcr.ifthenelse(pcr.defined(aha_ids), pcr.nominal(0), pgb_ids)
pgb_ids = pcr.nominal(pcr.areaminimum(pcr.scalar(pgb_ids_ids_exclusion), pgb_ids))
pgb_ids = pcr.ifthen(pcr.scalar(pgb_ids) > 0.0, pgb_ids)


# read all attribute of pcrglobwb water bodies, and exclude the ones identified in AHA
pgb_fracWaterInp_file = "/scratch-shared/edwinaha/dynqual_waterbodies/fracWaterInp_waterBodies5ArcMin_2010.map"
pgb_resMaxCapInp_file = "/scratch-shared/edwinaha/dynqual_waterbodies/resMaxCapInp_waterBodies5ArcMin_2010.map"
pgb_resSfAreaInp_file = "/scratch-shared/edwinaha/dynqual_waterbodies/resSfAreaInp_waterBodies5ArcMin_2010.map"
pgb_waterBodyIds_file = "/scratch-shared/edwinaha/dynqual_waterbodies/waterBodyIds_waterBodies5ArcMin_2010.map"
pgb_waterBodyTyp_file = "/scratch-shared/edwinaha/dynqual_waterbodies/waterBodyTyp_waterBodies5ArcMin_2010.map"
pgb_fracWaterInp = pcr.ifthen(pcr.defined(pgb_ids), pcr.readmap(pgb_fracWaterInp_file))
pgb_resMaxCapInp = pcr.ifthen(pcr.defined(pgb_ids), pcr.readmap(pgb_resMaxCapInp_file))
pgb_resSfAreaInp = pcr.ifthen(pcr.defined(pgb_ids), pcr.readmap(pgb_resSfAreaInp_file))
pgb_waterBodyIds = pcr.ifthen(pcr.defined(pgb_ids), pcr.readmap(pgb_waterBodyIds_file))
pgb_waterBodyTyp = pcr.ifthen(pcr.defined(pgb_ids), pcr.readmap(pgb_waterBodyTyp_file))


# read/set AHA properties and merge them to pcrglobwb
# - ids and extent
aha_waterBodyIds    = aha_ids
merged_waterBodyIds = pcr.cover(aha_waterBodyIds, pgb_ids)
# - fraction water 
aha_fracWaterInp_file = "/home/edwinaha/github/edwinkost/locating_dams/scripts/existing_reservoir_fraction_water_ids.map"
aha_fracWaterInp      = pcr.readmap(aha_fracWaterInp_file)
merged_fracWaterInp   = pcr.cover(aha_fracWaterInp, pgb_fracWaterInp)
# - reservoir capacity
aha_resMaxCapInp_file = "/home/edwinaha/github/edwinkost/locating_dams/scripts/existing_reservoir_capacity_ids_million_m3.map"
aha_resMaxCapInp      = pcr.readmap(aha_resMaxCapInp_file)
merged_resMaxCapInp   = pcr.cover(aha_resMaxCapInp, pgb_resMaxCapInp)
# - reservoir surface area
aha_resSfAreaInp_file   = "/home/edwinaha/github/edwinkost/locating_dams/scripts/existing_reservoir_surface_area_ids_m2.map"
aha_resSfAreaInp        = pcr.readmap(aha_resSfAreaInp_file)
merged_aha_resSfAreaInp = pcr.cover(aha_resSfAreaInp, aha_resSfAreaInp)
# - reservoir type
aha_waterBodyTyp      = pcr.ifthen(pcr.defined(aha_ids), pcr.scalar(2.0))
merged_waterBodyTyp   = pcr.cover(aha_waterBodyTyp, pgb_waterBodyTyp)

# save all reservoir types
pcr.report(merged_waterBodyIds, "merged_waterBodyIds_existing.map")
pcr.report(merged_fracWaterInp, "merged_fracWaterInp_existing.map")
pcr.report(merged_resMaxCapInp, "merged_resMaxCapInp_existing.map")
pcr.report(merged_resMaxCapInp, "merged_resMaxCapInp_existing.map")
pcr.report(merged_waterBodyTyp, "merged_waterBodyTyp_existing.map")

