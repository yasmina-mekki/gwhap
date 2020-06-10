########################################################################
# Brainomics - Copyright (C) CEA, 2020
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html
# for details.
#
# Author : V Frouin, S Karkar, Y Mekki
########################################################################

########################################################################
##  libraries

library(gwhap)
library(parallel)
library(optparse)
library(readr)
library(data.table)

########################################################################
##  root pathes

data_root='/neurospin/ukb'
# genetic map root 
genmap_root = paste(sep='',
      '/neurospin/brainomics/bio_resources',
      '/hg19_maps',
      '/interpolated_maps',
      '/1000-genomes-genetic-maps-master/interpolated_from_hapmap'
      )
# ukb root to the first pass of imputation - so-called phased data
#                                          - aka  haplotyped data
ukb_phased_root = paste(sep='', data_root,
     '/genetic/GENETIC_DATA_500k/HAPLOTYPES/v2')

# output name
output_root = '/neurospin/tmp/gwhap_ejhg'

########################################################################
##  option parsing

pipe = c(
    "Load genetic map",                                         #1
    "Load bucket of SNP",                                       #2
    "Interpolation from the snp bucket into the interp. map",   #3
    "Block splitting",                                          #4
    "Read phenotype file",                                      #5
    "Adjust and align genetic/phenotype data",                  #6
    "Haplotype calling",                                        #7
    "Haplotype testing"                                         #8
)
mask_pipe = pipe

########################################################################
##  option parsing

option_list <- list(
  make_option(c("-d", "--startrun"),
              action="store_true", default=FALSE,
              help="Start run to fix path and resources [default]"),
  make_option(c("-s", "--save"),
              action="store_true", default=FALSE,
              help="Save intermediate data [default]"),
  make_option(c("", "--burnin"),
              action="store_true", default=FALSE,
              help="Prepare haplotype block variant for a given delta [default]"),
  make_option(c("", "--testsonly"),
              action="store_true", default=FALSE,
              help="Read the prepared haplotype block and run tests [default]"),
  make_option(c("", "--phenotypeall"),
              action="store_true", default=FALSE,
              help="The process is applied to all phenotypes available. [default]"),
  make_option(c("-n", "--nb_core"),
              type="integer", default=1,
              help="Number of cores [default %default]"),
  make_option(c("-c", "--chromosome"),
              type="integer", default=21,
              help="Chromosome number [default %default]"),
  make_option(c("-p", "--phenotype"),
              type="character", default="FCLa_right",
              help="Phenotype name [default %default]"),
  make_option(c("-k", "--kind"),
              type="character", default="block",
              help="Kind of test to consider (comma separeted list c('single','block','complete'))[default %default]"),
  make_option(c("-o", "--outputdir"),
              type="character", default=output_root,
              help="Output dir [default %default]")
)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser, positional_arguments = 1)
opt <- args$options # opt accessible like opt$outputdir
                    # args$args will be sued as pheno full path


chromosome <- opt$chromosome
save       <- opt$save
startrun   <- opt$startrun
save       <- opt$save
outputdir  <- opt$outputdir
nb_core    <- opt$nb_core
list_col_phenotype <- unlist(strsplit(opt$phenotype, split=',')) # name of the specifi sulcus in pheno file
phenotype_fpath <- args$args
burnin     <- opt$burnin
testsonly  <- opt$testsonly
phenotypeall <- opt$phenotypeall
if (all(sapply(opt$kind,  function(kk) kk %in% c('single', 'complete', 'block')))) {
    kind <-opt$kind
} else
{
    print('Kind should be in c(single, complete, complete)')
    exit(1)
}

#~ chromosome <- 21
#~ save <- TRUE
#~ startrun <- T
#~ outputdir <- '/neurospin/tmp/gwhap_ejhg'
#~ nb_core <- 8
#~ list_col_phenotype <- c("FCLa_right", "FCLa_left")
#~ phenotype_fpath = file.path(data_root, '/reports-internal/Haplo-paper-2020-/all_pheno_res.tsv')
#~ burnin     <- F
#~ testsonly  <- T
#~ phenotypeall <- F
#~ kind <- c('block', 'complete', 'single')
#~ kind <- c('block')
                    
# Do some pre-interpretations
if (burnin) {mask_pipe=pipe[1:7]}
if (testsonly) {mask_pipe=pipe[c(5, 6, 8)]}
if (phenotypeall) {cat(sprint("--phenotypeall set. All phenotypes will be considered"))}



if (TRUE){
    cat(sprintf("Parameters: chrom=%d\n            pheno=%s\nSetup:\tsave=%s\tstartrun=%s\tnb_core=%d\tkind=%s\n\n",
        chromosome, phenotype_fpath, save, startrun, nb_core, kind)
        )
    }


########################################################################
##  full file pathes

# haplotype or phased data
# these are the SNPs from the UK believe chip that pass the first step
# of the imputation process
bgen_fpath.bgi = file.path(ukb_phased_root, 
                    sprintf('ukb_hap_chr%d_v2.bgen.bgi',chromosome))
bgen_fpath.bgen = file.path(ukb_phased_root, 
                    sprintf('ukb_hap_chr%d_v2.bgen', chromosome))
bgen_fpath.samples = file.path(ukb_phased_root,
                    sprintf('ukb25251_hap_chr%d_v2_s487395.sample', chromosome))



# outputs
output_blocks_dir = file.path(output_root, 'blocks')
if (save) {dir.create(output_blocks_dir, showWarnings=F, recursive=T)}

output_haplo_dir = file.path(output_root, 'haps')
if ("Haplotype calling" %in% mask_pipe) {dir.create(output_haplo_dir, showWarnings=F, recursive=T)}

output_test_dir = file.path(output_root, 'tests')
invisible(lapply(kind, function(i) {
    if (save) { dir.create(file.path(output_test_dir,i), showWarnings=F, recursive=T)}
}))
########################################################################
##  a few checks
file.exists(bgen_fpath.bgen)
file.exists(bgen_fpath.bgi)
file.exists(bgen_fpath.samples)


########################################################################
# Pipeline start here
########################################################################

########################################################################
##  1) Load the reference genetic map
##     instanciate a Genetic_Map object
##     a) Two parameters have to be build filepaths and encodings (see  after)
##
##     b) First create the object 
##     c) and then read it

pipe.step = "Load genetic map"
if (pipe.step %in% mask_pipe) {

start.time <- Sys.time()
## the a) part
f1 = file.path(genmap_root, 
                sprintf('chr%d.interpolated_genetic_map.gz', chromosome))
# build a NAMED list mapping the chrom information
chr = list(chromosome)
names(chr) = c(f1)
# now build the parameter for Genetic_Map object. Inspect the file to
# assign the right column name for cM and position.
# Here the chrom info is not read from the file it is read from the
# NAMED list chr
filepaths = c(f1)
encodings = list("cM"="V3", "position"="V2","chr"=chr, "format"="")

## the b) part get an instance of Genetic_Map
genetic_map = Genetic_Map(filepaths=filepaths, encodings=encodings)
# print(genetic_map)

## the c) part finally read
genetic_map = readData(genetic_map)
# print(genetic_map)

## summary for this step
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(sprintf("Step %s in %s sec.\n",pipe.step, time.taken))
}

########################################################################
##  2) Load the bucket of SNP on which the blocks will be determined
##     instanciate a Snp_Bucket object (get snp and physical position)
##     a) Two parameters have to be build filepaths and encodings (see  after)
##
##     b) First create the object 
##     c) and then read it

pipe.step = "Load bucket of SNP"
if (pipe.step %in% mask_pipe) {

start.time <- Sys.time()
##  Part a)
fsnpb = bgen_fpath.bgi
filepaths = c(fsnpb)
chr = list(chromosome)
names(chr) = c(fsnpb)
encodings = list('snp'='snp', 'position'='position', 'chr'=chr, "format"="bgen")
## Part b)
snp_bucket = Snp_Bucket(filepaths=filepaths, encodings=encodings)
## Part c)
snp_bucket = readData(snp_bucket)

## summary for this step
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(sprintf("Step %s in %s sec.\n",pipe.step, time.taken))
}

########################################################################
##  3) Perform the interpolation in each snp_bucket snp (based on the
##     reference genetic map genetic_map obtained in 1)

pipe.step = "Interpolation from the snp bucket into the interp. map"
if (pipe.step %in% mask_pipe) {

start.time <- Sys.time()
interpolated_map=create_augmented_genetic_map(
                        snp_bucket=snp_bucket,
                        genetic_map=genetic_map,
                        save_genetic_map=FALSE)

## summary for this step
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(sprintf("Step %s in %s sec.\n",pipe.step, time.taken))
}

########################################################################
##  4) Split in bloc according to the delta value (see EJHG paper)

pipe.step = "Block splitting"

if (pipe.step %in% mask_pipe) {

start.time <- Sys.time()
df_blocs=create_blocs(genetics_map=interpolated_map, 
                         delta=1e-3,
                         save_blocs=FALSE,
                         output=output_dir)

# transforme list of df in df
df_blocs=do.call(rbind,df_blocs)

blocs.start = df_blocs$from_bp[df_blocs$chr == sprintf('chr%d', chromosome)]
blocs.end   = df_blocs$to_bp[df_blocs$chr == sprintf('chr%d', chromosome)]


## summary for this step
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(sprintf("Step %s in %s sec.\n",pipe.step, time.taken))
}

########################################################################
##  5) Read phenotype data

pipe.step = "Read phenotype file"

if (pipe.step %in% mask_pipe) {

start.time <- Sys.time()
options(warn=-1)
phenotype = read_delim(phenotype_fpath, delim='\t')
if (phenotypeall) {
    list_col_phenotype <- colnames(phenotype)[3:ncol(phenotype)]
}

## summary for this step
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(sprintf("Step %s in %s sec.\n",pipe.step, time.taken))
}

########################################################################
##  6) Synchronize genetic and phenotype data
##      a) get the genetic samples available 
##      b) load the genetic data as phased data and
##         (not as snp_bucket like in 3).


pipe.step = "Adjust and align genetic/phenotype data"

if (pipe.step %in% mask_pipe) {

start.time <- Sys.time()
## part a) read the index file & remove the first line
sample_index = read_delim(bgen_fpath.samples, delim='  ')
sample_index = sample_index[-1, ]

# get the index as columns
sample_index=setDT(sample_index, keep.rownames = TRUE)[] #avoiding output

# filter on the phenotype IID and get the index as well as the participant IID ordered
filtered_sample_index = as.numeric(sample_index[sample_index$ID_2 %in% phenotype$IID, ]$rn)
filtered_sample_index_ID = as.numeric(sample_index[sample_index$ID_2 %in% phenotype$IID, ]$ID_2)

## part b) Get a data_loader of the phased_data

phased_dl = phased_data_loader.bgen(bgen_fpath.bgen)

# get snps position using the bgi file
annot = getAnnotVariants(phased_dl)

# filter the partcipant bgen ID using their index using the 
# phased_data_loader object methods

internal_samples = getInternalIID(phased_dl)
internal_samples = internal_samples[filtered_sample_index]

#~ filtered_sample_index[1:5]
#~ filtered_sample_index_ID[1:5]
#~ filtered_sample_index_ID[1:5] %in% phenotype$FID
#~ phenotype[phenotype$FID %in% filtered_sample_index_ID[1:5] , ]

## summary for this step
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(sprintf("Step %s in %s sec.\n",pipe.step, time.taken))
}

########################################################################
##  7) Determine the haplotypes
##     This step is modulated by the startrun flag to fix resources

pipe.step = "Haplotype calling"

if (pipe.step %in% mask_pipe) {

start.time <- Sys.time()
if (startrun) {
    dryrun.nb_block = 5
    dryrun.nb_subj = 1000

    hap = determine_haplotypes_per_chromosome(phased_dl,
                                chromosome=chromosome,
                                start=blocs.start[1:dryrun.nb_block],
                                end=blocs.end[1:dryrun.nb_block],
                                sample_iid = filtered_sample_index_ID[1:dryrun.nb_subj],
                                sample_bgen_iid_code = internal_samples[1:dryrun.nb_subj],
                                nb_core=nb_core #detectCores()-1
                                )
} else
{
    hap = determine_haplotypes_per_chromosome(phased_dl,
                                chromosome=chromosome,
                                start=blocs.start,
                                end=blocs.end,
                                sample_iid = filtered_sample_index_ID,
                                sample_bgen_iid_code = internal_samples,
                                nb_core=nb_core #detectCores()-1
                                )
}

if (save) {
    save_haplotypes(hap, chromosome=sprintf('chr%d', chromosome), output=output_haplo_dir)
}

## summary for this step
end.time <- Sys.time()
time.taken <- end.time - start.time
if (startrun) {
    cat(sprintf("Step %s in %s sec. [STARTRUN mode :%d blocks and %d subjects]",
                  pipe.step, time.taken, dryrun.nb_block, dryrun.nb_subj)
         )
} else {
    cat(sprintf("Step %s in %s sec.\n",pipe.step, time.taken))
}

}

## ---------------------------------------------------------------------
# CHEMIN DES HAPLOTYPES calculé dans le papier
#/neurospin/ukb/genetic/GENETIC_DATA_500k/HAPLOTYPES/blocks_19k_gz/all_sulci_boxcox.csv.hapDF_chr21*gz
## ---------------------------------------------------------------------

########################################################################
##  8) Perform the test
##     Three tests : block, single, complete
##      a) 

pipe.step = "Haplotype testing"

if (pipe.step %in% mask_pipe) {

start.time <- Sys.time()
# check wether hap data are available and try to read them from cache
if (!("hap" %in% ls())) {
    tryCatch({
        hap = load_haplotypes(chromosome=sprintf('chr%d', chromosome), output_haplo_dir)
    },
    error = function(error_condition) {
        fp = file.path(output_haplo_dir, sprintf('chr%d', chromosome))
        cat(sprintf("Cannot read procomputed haplotypes from %s", fp))
    }
    )
}


# do the vecorization and mc run
haplo.start = na.omit(unique(vapply(strsplit(colnames(hap),"_"), `[`, 2, FUN.VALUE=character(1))))
haplo.end = na.omit(unique(vapply(strsplit(colnames(hap),"_"), `[`, 3, FUN.VALUE=character (1))))
haplo.blocks = sprintf('chr%d_%s_%s', chromosome, haplo.start, haplo.end)
haplo.names = colnames(hap)



# Now inited from optarg : 
# kind = c("block", "single", "complete")
outkind = lapply(kind, function(p) file.path(output_test_dir, p))
phe_outkind = lapply(outkind, function(p) file.path(p, list_col_phenotype))
names(phe_outkind) = names(outkind) = kind
for (k in kind){names(phe_outkind[[k]])=list_col_phenotype} 

#if (STOP){}

if (startrun) {
    dryrun.nb_subj = 1000
    
    for (k in kind){
        for (col_phenotype in list_col_phenotype) {
            test= do.call(rbind,
                         mcmapply(
                            FUN=lm_test_haplotypes,
                            X=lapply(haplo.blocks, function(b) hap[1:dryrun.nb_subj,startsWith(haplo.names, b)]),
                            MoreArgs=list(
                                Y=as.data.frame(phenotype)[1:dryrun.nb_subj, col_phenotype, drop=FALSE],
                                kind=k),
                           mc.cores=nb_core
                           )
                        )
            rownames(test) = c()
            if (save){
                dir.create(phe_outkind[[k]][[col_phenotype]], showWarnings=F, recursive=T)
                save_tests(test, chromosome, phe_outkind[[k]][[col_phenotype]])
            }
        }
    }
} else {
    for (k in kind){
        for (col_phenotype in list_col_phenotype) {
            test= do.call(rbind,
                         mcmapply(
                            FUN=lm_test_haplotypes,
                            X=lapply(haplo.blocks, function(b) hap[,startsWith(haplo.names, b)]),
                            MoreArgs=list(
                                Y=as.data.frame(phenotype)[, col_phenotype, drop=FALSE],
                                kind=k),
                           mc.cores=nb_core
                           )
                        )
            rownames(test) = c()
            if (save){
                dir.create(phe_outkind[[k]][[col_phenotype]], showWarnings=F, recursive=T)
                save_tests(test, chromosome, phe_outkind[[k]][[col_phenotype]])
            }
        }
    }
}

#~         test= do.call(rbind,
#~                      mcmapply(
#~                         FUN=lm_test_haplotypes,
#~                         X=lapply(haplo.blocks, function(b) hap[,startsWith(haplo.names, b)]),
#~                         MoreArgs=list(
#~                             Y=as.data.frame(phenotype)[, list_col_phenotype, drop=FALSE],
#~                             kind=k),
#~                        mc.cores=nb_core
#~                        )
#~                     )
#~         rownames(test) = c()
#~         if (save){
#~             phe_outkind = lapply(outkind, function(p) file.path(p, list_col_phenotype))
#~             dir.create(phe_outkind[[k]], showWarnings=F, recursive=T)
#~             for (ph in list_col_phenotype){
#~                 save_tests(test[which(test$test==k & test$phname==ph),],
#~                        chromosome, phe_outkind[[k]])
#~                    }
#~         }


## summary for this step
end.time <- Sys.time()
time.taken <- end.time - start.time
cat(sprintf("Step %s in %s sec.\n",pipe.step, time.taken))
}
