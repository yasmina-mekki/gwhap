if (!require(devtools)){install.packages("devtools")}
remove.packages('gwhap')
devtools::install_github("yasmina-mekki/gwhap", force = TRUE)
library(gwhap)



#_____________________________________________________________
# Part 1
# Update the genetic map
#_____________________________________________________________

# Update the genetic map

sprintf('Time de lancement : %s', Sys.time())
snp_physical_positions = "/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb_hap_"
rutgers_map = "/neurospin/brainomics/bio_resources/rutgers_map_v3"
output = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map'
create_augmented_genetic_map(snp_physical_positions=snp_physical_positions, genetic_map_dir=rutgers_map, save_genetic_map=TRUE, output=output)
sprintf('Time de fin : %s', Sys.time())

# Read the genetic map -not mine-
load("/neurospin/brainomics/2019_gwhap/data/genMap_all.Rdata")

# read the augmented genetic map
# localy
augmented_genetic_map_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map'
augmented_map_df = get_augmented_genetic_map(augmented_genetic_map_dir, chromosomes=1:23)


#_____________________________________________________________
# Part 1
# Debugging
#_____________________________________________________________
sprintf('Time de lancement : %s', Sys.time())
snp_physical_positions = "/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb_hap_"
rutgers_map = "/neurospin/brainomics/bio_resources/rutgers_map_v3"
output = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map_approx_function'
create_augmented_genetic_map(snp_physical_positions=snp_physical_positions, genetic_map_dir=rutgers_map, save_genetic_map=TRUE, output=output)
sprintf('Time de fin : %s', Sys.time())

# Read the genetic map -not mine-
load("/neurospin/brainomics/2019_gwhap/data/genMap_all.Rdata")

# read the augmented genetic map
# localy
augmented_genetic_map_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map_approx_function'
augmented_map_df = get_augmented_genetic_map(augmented_genetic_map_dir, chromosomes=1:23)
