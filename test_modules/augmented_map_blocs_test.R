library(tools)
library(rbgen)
library(DBI)
library(data.table)
library(readr)
library(R.utils)


#_____________________________________________________________
# Part 1
# Update the genetic map using the rutgers map
#_____________________________________________________________

# Update the genetic map

sprintf('Start time : %s', Sys.time())
snp_physical_positions = "/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb_hap_"
rutgers_map = "/neurospin/brainomics/bio_resources/rutgers_map_v3"
output = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map'
create_augmented_genetic_map(snp_physical_positions=snp_physical_positions, genetic_map_dir=rutgers_map, map_name='rutgers', save_genetic_map=TRUE, output=output)
sprintf('End time : %s', Sys.time())

# Read the genetic map -not mine-
load("/neurospin/brainomics/2019_gwhap/data/genMap_all.Rdata")

# read the augmented genetic map
# localy
augmented_genetic_map_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map'
augmented_map_df = get_augmented_genetic_map(augmented_genetic_map_dir, chromosomes=1:23)


#_______________________________________________________________________________________
# Part 2
# Update the genetic map using the rutgers map using 1000 genome interpolation map
#_______________________________________________________________________________________

# Update the genetic map

sprintf('Start time : %s', Sys.time())
snp_physical_positions = "/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb_hap_"
m1000_genome_map = "/neurospin/brainomics/bio_resources/hg19_maps/interpolated_maps/1000-genomes-genetic-maps-master/interpolated_from_hapmap"
output = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map_1000_genome_interpolated'
create_augmented_genetic_map(snp_physical_positions=snp_physical_positions, genetic_map_dir=m1000_genome_map, map_name='1000_genome_interpolated', save_genetic_map=TRUE, output=output)
sprintf('End time : %s', Sys.time())

# Read the genetic map -not mine-
load("/neurospin/brainomics/2019_gwhap/data/genMap_all.Rdata")

# read the augmented genetic map
# localy
augmented_genetic_map_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map_1000_genome_interpolated'
augmented_map_df = get_augmented_genetic_map(augmented_genetic_map_dir, chromosomes=1:23)

# debugging
df_merged_by_rsid = merge(augmented_map_df, genMap_all, by.x = 'rsid', by.y = 'rsid')
same <- df_merged_by_rsid[df_merged_by_rsid$cM.x==df_merged_by_rsid$cM.y,]
setdiff(df_merged_by_rsid$cM.x, df_merged_by_rsid$cM.y)

