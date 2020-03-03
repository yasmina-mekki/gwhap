# deviner le format des inputs ...
# Y uniquement la valeur de tes phenotype. avec une column name.

#__________________________________________________________________________________________________
# set the parameters

# haplotypes
output = '/neurospin/tmp/ymekki/new_repo/gwhap/haplotypes_7_bloc_per_chromosome'
chr = 'chr1'
load(sprintf('%s/haplotypes_%s.RData', output, chr))
X = haplotype_combined

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
Y = phenotype[phenotype$IID %in% sample_index_file$ID_2, ]

#___________________________________________________________________________________________________

# get the index columns
setDT(haplotype_combined, keep.rownames = TRUE)[]
df_tmp_merged = merge(sample_index_file, haplotype_combined, by="rn")
df_merged = merge(df_tmp_merged, phenotype, by.x = 'ID_2', by.y = 'IID')
df_merged$ID_2 <- NULL
df_merged$ID_1 <- NULL
df_merged$rn <- NULL
df_merged$missing <- NULL
df_merged$FID <- NULL
df_merged$sex <- NULL


Y$FID <- NULL
Y$IID <- NULL
#___________________________________________________________________________________________________

# get the columns name
col_orgY = colnames(Y)
col_orgX = colnames(X)

# appliquer une somme sur les colonnes de X
# nombre d'individus ayant cet haplotype ? que signifi les nombre 1, 2, 0 pour les haplotypes (dummification) 0, 1 pour présence ou pas mais le 2 ?
allelcount = unlist(apply(X, 2, sum))
varL = colnames(X)

# get the index of the max of allel count
# get haplotypes columns without the common one.
if (length(varL) > 1) {varL = colnames(X)[setdiff(1:ncol(X), which.max(allelcount))]}


form_Y = paste(col_orgY, collapse = ',')
form_X = paste(varL, collapse = '+')


full_Lm = lapply(sprintf("cbind( %s ) ~ %s ", col_orgY, form_X),
                 lm,
                 data = data.frame(cbind(Y, X)))


sum_Lm = lapply(full_Lm,summary)
names(sum_Lm) = col_orgY # rownames(full_Lm$coefficients)
