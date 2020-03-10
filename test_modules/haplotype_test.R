remove.packages('gwhap')
devtools::install_github("yasmina-mekki/gwhap", force = TRUE)
library(gwhap)


#_____________________________________________________________
# Part 3
# Determine haplotypes
#_____________________________________________________________
library(data.table)
library(readr)
library(rlist)
library(DBI)
library(rbgen)
library(dummies)
library(parallel)
#install.packages("vcfR")
#library(vcfR)

# load blocs
load('/neurospin/brainomics/2019_gwhap/data/winDF_all.Rdata')
winDF_1e_3 = winDF[winDF$min_cM==1e-3, ]

# bgen file
bgen_file = '/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb_hap'

# Resting state phenotypes
#phenotype_path = '/neurospin/brainomics/2019_ln_YME/GWAS/UKB_v2/imputed_data/non_filtred/D/Putamen_d_IFGorb_d/Putamen_d_IFGorb_d.txt'
#phenotype = read_delim(phenotype_path, delim='\t')

# Structural phenotypes
phenotype_path = '/neurospin/ukb/workspace/pheno_boxcox/all_pheno_res.tsv'
phenotype = read_delim(phenotype_path, delim='\t')

# index path
sample_semi_path = '/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb25251_hap'

# output
output = '/neurospin/tmp/ymekki/new_repo/gwhap/T1_subjects/haplotypes'

sprintf('Start time determine haplotypes : %s', Sys.time())
for(chr in unique(winDF_1e_3$chr)){
  sprintf('Start time determine haplotypes chr %s : %s', chr, Sys.time())
  start   = winDF_1e_3$fromBp[winDF_1e_3$chr == chr]
  end     = winDF_1e_3$toBp[winDF_1e_3$chr == chr]
  bgnfile = sprintf('%s_%s_v2.bgen', bgen_file, chr)

  # read index file
  sample_index_file = read_delim(sprintf('%s_%s_v2_s487395.sample', sample_semi_path, chr),
                                 delim='  ')
  # get the index columns
  setDT(sample_index_file, keep.rownames = TRUE)[]
  # filtre on the phenotype IID and get the index
  samples_index = as.numeric(sample_index_file[sample_index_file$ID_2 %in% phenotype$IID, ]$rn)
  
  bgifile = sprintf("%s.bgi", bgnfile)
  allVar = get_bgi_file(file_path = bgifile)
  dummybg = bgen.load(
    filename = bgnfile,
    data.frame(
      chromosome = '',
      start = allVar$position[1],
      end = allVar$position[1]
    ),
    max_entries_per_sample = 4
  )
  
  # filter the partcipant ID using their index
  mysamples = dummybg$samples[samples_index]
  # determine the haplotypes and binds all list elements by column
  haplotype_combined <- do.call(cbind, mcmapply(FUN     = determine_haplotypes_per_bloc,
                                               chr     = chr,
                                               start   = start,
                                               end     = end,
                                               bgnfile = bgnfile,
                                               samples_index = list(samples_index),
                                               mysamples = list(mysamples),
                                               mc.cores = 32))

  # save the haplotypes
  
  haplotype_combined_path = sprintf('%s/haplotypes_%s.RData', output, chr)
  save(haplotype_combined, file=haplotype_combined_path, compress=T)
  sprintf('End time determine haplotypes chr %s : %s', chr, Sys.time())
}
sprintf('End time determine haplotypes : %s', Sys.time())


#______________________________________________________________________________________________________________
# read the haplotypes (reconstruction des l'haplotypes pour un bloc)
library(stringr)
load('/neurospin/tmp/ymekki/new_repo/gwhap/haplotypes/haplotypes_chr20.RData')
head(haplotype_combined)
chr = "chr10"

colnames(haplotype_combined) = gsub("NA", chr, colnames(haplotype_combined))
start = unique(vapply(strsplit(colnames(haplotype_combined),"_"), `[`, 2, FUN.VALUE=character(1)))
end = unique(vapply(strsplit(colnames(haplotype_combined),"_"), `[`, 3, FUN.VALUE=character(1)))

blocs = sprintf('%s_%s', start, end)

columns_blocs = str_subset(colnames(haplotype_combined), blocs[1])
X = haplotype_combined[, columns_blocs]

#lapply(str_contains, x = colnames(haplotype_combined), pattern = blocs)















