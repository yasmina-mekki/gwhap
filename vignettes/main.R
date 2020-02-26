if (!require(devtools)){install.packages("devtools")}
devtools::install_github("yasmina-mekki/gwhap", force = TRUE)
library(gwhap)

# should not load dplyr from main !
library(dplyr)

# results folder
output = './results/'
RESOURCE_ROOT='/neurospin/brainomics/bio_resources'

#_____________________________________________________________
# Part 1
# Update the genetic map, create the blocs and visualise them
#_____________________________________________________________

# Update the genetic map
to_be_updated_genetic_map_path = "./data/small_region_imputed.bim"
rutgers_map = file.path(RESOURCE_ROOT,"rutgers_map_v3")
gen_map_updated = genetic_map_update(genetic_map_path=rutgers_map, snp_physical_positions=to_be_updated_genetic_map_path)

# Create blocs
df_blocks = create_blocks(gen_map_updated, delta=1e-3)

# Create blocs with multiple values of delta
deltas = c(1e-3, 2.5e-3, 5e-3, 7.5e-3)
dfs_blocks = data.frame()
for(delta in deltas){
  df_blocks_tmp = create_blocks(genetics_map=gen_map_updated, delta=delta)
  dfs_blocks <- rbind(dfs_blocks, df_blocks_tmp)
}

# Block distribution plot

#_________________________________________________
# from here now on, winDF_all is going to be used
# winDF_all : blocks created using ukb_hap dataset
#             for delta from 0.001 to 0.025 cM
#_________________________________________________

# get winDF_all
dfs_blocks = get(load("./data/winDF_all.Rdata"))
colnames(dfs_blocks) = c("start_w", "end_w", "win_size", "from_bp", "to_bp", "from_cm", "to_cm", "chr", "delta")


# All Chromosomes
haplotype_block_distribution_per_delta = haplotype_block_distribution(dfs_blocks, xlim=c(1.590790e-10, 1), ylim=c(0, 1.5))
print(haplotype_block_distribution_per_delta)
save_plot(haplotype_block_distribution_per_delta, paste0(output, 'haplotype_block_distribution_per_delta.pdf'))

# Per chromosome
haplotype_block_distribution_per_chr = haplotype_block_distribution(dfs_blocks, xlim=c(1.590790e-10, 1), ylim=c(0, 1.5), per_chromosome=TRUE)
print(haplotype_block_distribution_per_chr)
save_plot(haplotype_block_distribution_per_chr, paste0(output, 'haplotype_block_distribution_per_chr.pdf'))

# Karyotype
df_blocs_delta_3 = dfs_blocks[dfs_blocks$delta==0.001,]
karyotype_plot_obj = karyotype_plot(df_blocs_delta_3)
save_plot(karyotype_plot_obj, paste0(output, 'karyotype_plot.pdf'))


#______________________________________________________
# Part 2
# determine the haplotypes
#______________________________________________________



#______________________________________________________
# Part x
# Phenotypes distribution plots
#______________________________________________________

# Phenotypes distribution plots
phenotypes_path = Sys.glob(file.path("./data/ln_partial_correlation_eth_white/" , '*'))
phenotype_distribution_plot = phenotype_distribution(phenotypes_path = phenotypes_path)
print(phenotype_distribution_plot)

# Save phenotypes distribution plot
save_plot(phenotype_distribution_plot, paste0(output, 'phenotype_distribution_plot.pdf'))
