
# Begin of the script
#########################################
data_root='../data'

bgen_fpath = file.path(data_root, 'ukb_hap_chr21_v2.bgen')
blocs_dir = file.path(data_root)
# index path
sample_fpath = file.path(data_root, 'ukb25251_hap_chr21_v2_s487395.sample')
# phenotype
phenotype_fpath = file.path(data_root, 'all_pheno_res.tsv')
# output
output_dir = file.path(data_root, 'output')


# chromosome
chromosome=c(21)
# get blocs
df_blocs  = get_blocs(blocs_dir, chromosomes=chromosome)
# read pheno
phenotype = read_delim(phenotype_fpath, delim='\t')


start     = df_blocs$from_bp[df_blocs$chr == chromosome]
end       = df_blocs$to_bp[df_blocs$chr == chromosome]


# read the index file & remove the first line
sample_index = read_delim(sample_fpath, delim='  ')
sample_index = sample_index[-1, ]

# get the index as columns
setDT(sample_index, keep.rownames = TRUE)[]

# filtre on the phenotype IID and get the index as well as the participant IID ordered
filtered_sample_index = as.numeric(sample_index[sample_index$ID_2 %in% phenotype$IID, ]$rn)
filtered_sample_index_ID = as.numeric(sample_index[sample_index$ID_2 %in% phenotype$IID, ]$ID_2)

# Get a data_loader of the phased_data
phased_dl = phased_data_loader.bgen(bgen_fpath)

# get snps position using the bgi file
annot = getAnnotVariants(phased_dl)
# filter the partcipant bgen ID using their index
internal_samples = getInternalIID(phased_dl)
#internal_samples = internal_samples[filtered_sample_index]

# determine haplotypes for all blocs of one chromosome
print(class(phased_dl))
phased_dl[["annot_internalIID"]]= phased_dl[["annot_internalIID"]][1:3]
phased_dl[["annot_variants"]]=phased_dl[["annot_variants"]][1:2,]

# in this example consider blocs #10, #11, #12
haplotype_combined = determine_haplotypes_per_chromosome(
                            phased_dl,
                            chromosome=chromosome,
                            start[10:12],
                            end[10:12],
                            sample_iid = filtered_sample_index_ID,
                            sample_bgen_iid_code = internal_samples[filtered_sample_index],
                            nb_core=detectCores()-1,
                            verbose=FALSE
                           )

# save the haplotypes in tsv file (very heavy, you should consider to use an other way to save the haplotypes)
#save_haplotypes(haplotype_combined, chromosome=chr, output=output_dir)

# save the haplotypes in RData object
save(haplotype_combined, file=sprintf('%s/haplotypes_%s.RData', output_dir, chromosome), compress=T)
