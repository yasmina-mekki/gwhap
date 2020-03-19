# deviner le format des inputs ...
# Y uniquement la valeur de tes phenotype. avec une column name.

#__________________________________________________________________________________________________
# set the parameters

# read the haplotypes (reconstruction des l'haplotypes pour un bloc)
library(stringr)
library(readr)
library(data.table)
library(dplyr)
library(parallel)


#_________________________________________________________________________________________________________
#???
#


load('/neurospin/tmp/ymekki/new_repo/gwhap/RS_subjects/haplotypes/haplotypes_chr10.RData')
#haplotype_combined
chr = "chr10"

# replace NA added by cbind by the chromosome code
colnames(haplotype_combined) = gsub("NA", chr, colnames(haplotype_combined))
# get the start and end bloc's position
start = unique(vapply(strsplit(colnames(haplotype_combined),"_"), `[`, 2, FUN.VALUE=character(1)))
end = unique(vapply(strsplit(colnames(haplotype_combined),"_"), `[`, 3, FUN.VALUE=character(1)))
# concatenate the start and end bloc position
blocs = sprintf('%s_%s', start, end)

# filtre sur le premier bloc ...
columns_blocs = str_subset(colnames(haplotype_combined), blocs[1])
haplotype_one_bloc = haplotype_combined[, columns_blocs]

#_________________________________________________________________________________________________________
# get the IID
#_________________________________________________________________________________________________________
chr = "chr5"

# phenotypes
phenotype_path = '/neurospin/brainomics/2019_ln_YME/GWAS/UKB_v2/imputed_data/non_filtred/D/Putamen_d_IFGorb_d/Putamen_d_IFGorb_d.txt'
phenotype = read_delim(phenotype_path, delim='\t')

# filtre sur les IDs common phenotype and genotype
# read index file
sample_semi_path = '/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2/ukb25251_hap'
sample_index_file = read_delim(sprintf('%s_%s_v2_s487395.sample', sample_semi_path, chr), delim='  ')

# get the index columns
setDT(sample_index_file, keep.rownames = TRUE)[]

# get common phenotype IID and genotype IID
phenotype_filtred = phenotype[phenotype$IID %in% sample_index_file$ID_2, ]

#___________________________________________________________________________________________________
# résidualisation
#___________________________________________________________________________________________________
covars_path = '/neurospin/brainomics/2019_ln_YME/covars/UKB/genetics/covar_GCTA/covar_sex_modified_GCTA_SexArray.cov'
qcovar_path = '/neurospin/brainomics/2019_ln_YME/covars/UKB/genetics/covar_GCTA/qcovar_GCTA_Age_10PCS.cov'
covars  = read_delim(covars_path, delim='\t')
qcovars = read_delim(qcovar_path, delim='\t')

cov_UKB = merge(covars, qcovars, by="IID")
cov_pheno_UKB = merge(cov_UKB, phenotype_filtred, by='IID')

X_covar = c("Array", "sex_bis", "Age", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")
Y_pheno = c('connectivity')
X_col_name = paste(X_covar, collapse = '+')
Y_col_name = paste(Y_pheno, collapse = ',')

# linear regression
full_Lm = lapply(sprintf("cbind( %s ) ~ %s ", X_col_name, Y_col_name), lm, data = cov_pheno_UKB)
sum_Lm = lapply(full_Lm, summary)

# get the residu
Y_residualisee = data.frame(sum_Lm[[1]]$residuals)
colnames(Y_residualisee) = c('connectivity_residualisee')

Y_cov_residualisee = cbind(Y_residualisee, cov_pheno_UKB)
Y_filtred_residualisee = Y_cov_residualisee[, c('IID', 'connectivity_residualisee')]

#___________________________________________________________________________________________________
# run the test part
#___________________________________________________________________________________________________



# format Y : IID, phenotype value

phenotypes_sulci = fread("/neurospin/ukb/workspace/pheno_boxcox/all_sulci_boxcox.tsv")
Y_filtred = phenotypes_sulci[, c('IID', 'FCMpost_left')]

sprintf('Start time : %s', Sys.time())
#output = '/neurospin/tmp/ymekki/new_repo/gwhap/RS_subjects/lm_test'
output = '/neurospin/tmp/ymekki/new_repo/gwhap/T1_subjects/lm_test'

for (i in 1:22){
  #haplotypes_path = '/neurospin/tmp/ymekki/new_repo/gwhap/RS_subjects/haplotypes/haplotypes_'
  haplotypes_path = '/neurospin/tmp/ymekki/new_repo/gwhap/T1_subjects/haplotypes/haplotypes_'
  
  load(sprintf('%schr%s.RData', haplotypes_path, i))
  #haplotype_combined
  chr = sprintf('chr%s', i)
  
  # separate the blocs
  # replace NA added by cbind by the chromosome code
  colnames(haplotype_combined) = gsub("NA", chr, colnames(haplotype_combined))
  # get the start and end bloc's position
  start = unique(vapply(strsplit(colnames(haplotype_combined),"_"), `[`, 2, FUN.VALUE=character(1)))
  end = unique(vapply(strsplit(colnames(haplotype_combined),"_"), `[`, 3, FUN.VALUE=character(1)))
  # concatenate the start and end bloc position
  blocs = sprintf('%s_%s', start, end)
  
  results_chr = mcmapply(FUN = lm_test_haplotypes_per_bloc,
                         blocs = blocs,
                         haplotype_combined = list(haplotype_combined),
                         sample_index_file = list(sample_index_file),
                         Y_filtred = list(Y_filtred),
                         mc.cores = 32)
  
  # results_chr is represented by matrix : rows representing the different test and columns representing each blocs
  # row bind the results test separetely
  bloc_test_results     = data.frame()
  complete_test_results = data.frame()
  single_test_results   = data.frame()
  for(col in colnames(results_chr)){
    bloc_test_results     = rbind(bloc_test_results, results_chr[1, col][[1]])
    complete_test_results = rbind(complete_test_results, results_chr[2, col][[1]])
    single_test_results   = rbind(single_test_results, results_chr[3, col][[1]])
  } 
  
  
  # save the results
  dir.create(sprintf('%s/chr_%s', output, i))
  write.table(bloc_test_results, sprintf('%s/chr_%s/bloc_test_results.tsv', output, i), sep="\t", row.names=FALSE)
  write.table(complete_test_results, sprintf('%s/chr_%s/complete_test_results.tsv', output, i), sep="\t", row.names=FALSE)
  write.table(single_test_results, sprintf('%s/chr_%s/single_test_results.tsv', output, i), sep="\t", row.names=FALSE)
}
sprintf('End time : %s', Sys.time())


lm_test_haplotypes_per_bloc = function(blocs, haplotype_combined, sample_index_file, Y_filtred){
  
  # change Y_filtred column name
  colnames(Y_filtred) = c('IID', 'phenotype')
  
  # get only the haplotypes corresponding to bloc got as input
  columns_blocs = str_subset(colnames(haplotype_combined), blocs)
  haplotype_one_bloc = haplotype_combined[, columns_blocs]
  
  # This is done in order to be sure that the subject are well alligned
  # get the index columns
  setDT(haplotype_one_bloc, keep.rownames = TRUE)[]
  # merge the index file and the haplotypes using the index
  df_tmp_merged = merge(sample_index_file, haplotype_one_bloc, by="rn")
  # merge the df obtained above with the phenotype using the ID
  df_merged = merge(df_tmp_merged, Y_filtred, by.x = 'ID_2', by.y = 'IID')
  
  # set X, Y and run the haplotypique test
  X = df_merged[, ..columns_blocs]
  Y = df_merged[, 'phenotype']
  results = lm_test_haplotypes(X, Y, kind='all')
  
  return(results)
}






