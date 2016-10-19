############ POPULATION GRIDDING AFRICA #################
# Combines SSP population projections and               #
# RCP urban fractions to grid population for the future #
#                                                       #
# Author: Niklas Boke-Olen                              #
# niklas.boke-olen@nateko.lu.se                         #
#                                                       #
#########################################################


R - Files need to be run in the following order


1. RCPurbanfrac_stacks.R   /// To create the urbanfraction stacked rasters needed later
2. create_unique_pop_data.R // create the unique dataset for year 2000
3. population_gridding.R  /// main population gridding -- runs in parallell specified by cores parameter
