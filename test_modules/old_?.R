if (!require(devtools)){install.packages("devtools")}
remove.packages('gwhap')
devtools::install_github("yasmina-mekki/gwhap", force = TRUE)
library(gwhap)



#_____________________________________________________________
# Part 1
# Update the genetic map
#_____________________________________________________________

# Update the genetic map


snp_physical_positions = "~/R_workspace/data/inputs/bgi_files/ukb_hap_"
rutgers_map = "~/workspace/r_projects/gwhap_test/data/rutgers_map"
create_augmented_genetic_map(snp_physical_positions=snp_physical_positions, genetic_map_dir=rutgers_map, save_genetic_map=TRUE, output='~/R_workspace/data')
sprintf('Time de fin : %s', Sys.time())

sprintf('Time de lancement : %s', Sys.time())
snp_physical_positions = "/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb_hap_"
rutgers_map = "/neurospin/brainomics/bio_resources/rutgers_map_v3"
output = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map'
create_augmented_genetic_map(snp_physical_positions=snp_physical_positions, genetic_map_dir=rutgers_map, save_genetic_map=TRUE, output=output)
sprintf('Time de fin : %s', Sys.time())

augmented_map_tmp = c()
for (chr in 1:22){
  augmented_chr = sprintf('~/R_workspace/data/augmented_map_chr%s.txt', chr)
  if(file.exists(augmented_chr)){
    augmented_map_tmp = rbind(augmented_map_tmp, read_delim(augmented_chr, delim='\t'))
  }
}

load("/neurospin/brainomics/2019_gwhap/data/genMap_all.Rdata")

#_____________________________________________________________
# Part 2
# Blocs
#_____________________________________________________________

#_____________________________________________________________
# Part 2.1
# Create blocs
#_____________________________________________________________

# read the augmented genetic map
# localy
augmented_genetic_map_dir = '~/R_workspace/data/outputs/augmented_genetic_map'
augmented_map_df = get_augmented_genetic_map(augmented_genetic_map_dir, chromosomes=1:23)
# mine
create_blocks(augmented_map_df, delta=1e-3, save_blocs=TRUE, output='~/R_workspace/data/outputs/blocs_yme')
# note mine
create_blocks(genMap_all, delta=1e-3, save_blocs=TRUE, output='~/R_workspace/data/outputs/blocs_cc')


# load the blocs

# mine
aa = get_blocs('~/R_workspace/data/outputs/blocs_yme')
# not mine
kk = get_blocs('~/R_workspace/data/outputs/blocs_cc')

# not mine original
load('~/R_workspace/gwhap_ressource/winDF_all.Rdata')
pp = winDF[winDF$min_cM == 1e-3]

get_blocs <- function(blocs_dir, chromosomes=1:23){
  blocs_df = c()
  for (chr in 1:22){
    blocs_chr = sprintf('%s/blocs_chr%s.txt', blocs_dir, chr)
    if(file.exists(blocs_chr)){blocs_df = rbind(blocs_df, read_delim(blocs_chr, delim='\t'))}
  }
  return(blocs_df)
}


#genome wide
sprintf('Start time lecture : %s', Sys.time())
augmented_genetic_map_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map'
augmented_map_df = get_augmented_genetic_map(augmented_genetic_map_dir, chromosomes=1:23)
sprintf('End time lecture : %s', Sys.time())

# Create blocs
sprintf('Start time creation blocs : %s', Sys.time())
output = '/neurospin/tmp/ymekki/new_repo/gwhap/blocs'
df_blocks = create_blocks(augmented_map_df, delta=1e-3, save_blocs=TRUE, output=output)
sprintf('End time creation blocs : %s', Sys.time())


#_____________________________________________________________
# Part 2.2
# visualize blocs
#_____________________________________________________________
blocs_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/blocs'
df_blocs = get_blocs(blocs_dir)

haplotype_block_distribution_per_delta = haplotype_block_distribution(df_blocs)
print(haplotype_block_distribution_per_delta)

haplotype_block_distribution_per_chr = haplotype_block_distribution(df_blocs, xlim=c(1.590790e-10, 1), ylim=c(0, 1.5), per_chromosome=TRUE)
print(haplotype_block_distribution_per_chr)

karyotype_plot_obj = karyotype_plot(df_blocs)


#_____________________________________________________________
# Part 3
# Determine haplotypes
#_____________________________________________________________

load('/neurospin/brainomics/2019_gwhap/data/winDF_all.Rdata')
bgen_file = '/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb_hap'
sample_semi_path = '/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb25251_hap'
winDF_1e_3 = winDF[winDF$min_cM==1e-3, ]

# read .sample file

for(chr in unique(winDF_1e_3$chr)){
  start = head(winDF_1e_3$fromBp[winDF_1e_3$chr == chr], 7)
  end   = head(winDF_1e_3$toBp[winDF_1e_3$chr == chr], 7)
  bgnfile  = sprintf('%s_%s_v2.bgen', bgen_file, chr)

  # index
  #sample_index_file = read_delim(sprintf('%s_%s_v2_s487395.sample', sample_semi_path, chr), delim='\t')
  #samples_index = 2:nrow(sample_index_file)
  samples_index = 2:10000

  list_haplotypes = mapply(FUN = determine_haplotypes, chr=chr, start=start, end=end, bgnfile=bgnfile, samples_index=list(samples_index))

  # binds all list elements by column
  haplotype_combined = list.cbind(list_haplotypes)

  # save the haplotypes
  # Use the sparse matrix ...
  haplotype_combined_path = sprintf('%s/haplotypes_chr%s.txt', output, chr)
  save(haplotype_combined, file=haplotype_combined_path, compress=T)
  break
}

# load phenotypes
library(readr)
library(data.table)
phenotype_path = '/neurospin/brainomics/2019_ln_YME/GWAS/UKB_v2/imputed_data/non_filtred/D/Putamen_d_IFGorb_d/Putamen_d_IFGorb_d.txt'
phenotype = read_delim(phenotype_path, delim='\t')

chr = 'chr10'
sample_semi_path = '/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb25251_hap'
sample_index_file = read_delim(sprintf('%s_%s_v2_s487395.sample', sample_semi_path, chr), delim='  ')
setDT(sample_index_file, keep.rownames = TRUE)[]

caca = as.numeric(sample_index_file[sample_index_file$ID_2 %in% phenotype$IID, ]$rn)


start = winDF_1e_3$fromBp[winDF_1e_3$chr == chr][1]
end   = winDF_1e_3$toBp[winDF_1e_3$chr == chr][1]
samples_index = 2:10000
a = determine_haplotypes("chr10", start, end, "/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb_hap_chr10_v2.bgen", samples_index)

# filtre sur les IDs
# conlationner les matrices et sauvegarder
options(max.print=1000000)
haplotype_combined = list.cbind(b)

mat <- matrix(haplotype_combined, ncol=ncol(haplotype_combined))
print(object.size(mat),units="auto")
mat_sparse <- Matrix(mat, sparse=TRUE)

sp_matrix <- sparseMatrix(i = i,
                          j = j,
                          x = x,
                          dims = list(nrow(haplotype_combined),ncol(haplotype_combined)))

writeMM(obj = haplotype_combined, file="~/Bureau/sparse_matrix.mtx")


A <- as(mat, "sparseMatrix")       # see also `vignette("Intro2Matrix")`
B <- Matrix(mat, sparse = TRUE)    # Thanks to Aaron for pointing this out

identical(A, B)

load("~/Bureau/test.Rbin")



