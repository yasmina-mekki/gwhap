#_____________________________________________________________
# Part 2
# visualize blocs
#_____________________________________________________________
blocs_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/blocs'
df_blocs = get_blocs(blocs_dir)

haplotype_block_distribution_per_delta = haplotype_block_distribution(df_blocs)
print(haplotype_block_distribution_per_delta)

haplotype_block_distribution_per_chr = haplotype_block_distribution(df_blocs, xlim=c(1.590790e-10, 1), ylim=c(0, 1.5), per_chromosome=TRUE)
print(haplotype_block_distribution_per_chr)

karyotype_plot_obj = karyotype_plot(df_blocs)
