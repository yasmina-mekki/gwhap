#_____________________________________________________________
# Part 2
# visualize blocs
#_____________________________________________________________
#blocs_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/blocs_using_gen_map_all'
#df_blocs = get_blocs(blocs_dir)
library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
library(karyoploteR)
library(grDevices)

load('/home/ym257837/workspace/R_workspace/khra/gwhap/ressource/winDF_all.Rdata')
colnames(winDF) <- c("startW", "endW", "winsz", "from_bp", "to_bp", "from_cm", "to_cm", "chr", "delta")

# for one delta
df_blocs = winDF[winDF$delta == 1e-3, ]
haplotype_block_distribution_per_delta = haplotype_block_distribution(df_blocs)
print(haplotype_block_distribution_per_delta)

# for multiple delta
dfs_blocs = winDF
colnames(winDF) <- c("startW", "endW", "winsz", "fromBp", "toBp", "from_cm", "to_cm", "chr", "delta")
haplotype_block_distribution_per_deltas = haplotype_block_distribution(dfs_blocs, xlim=c(1.590790e-10, 1e-02))
print(haplotype_block_distribution_per_deltas)

# for multiple delta and per chromosome
haplotype_block_distribution_per_chr = haplotype_block_distribution(dfs_blocs, xlim=c(1.590790e-10, 1), ylim=c(0, 1.2), per_chromosome=TRUE)
print(haplotype_block_distribution_per_chr)

# karyoplote
karyotype_plot_obj = karyotype_plot(df_blocs)











