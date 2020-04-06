#__________________________________________________________________________________________________
# TO DO :
# merge write_genetic_map and write_blocs functions
# merge get_augmented_genetic_map and get_blocs functions
#__________________________________________________________________________________________________

#' Get rutgers map
#'
#' @description Get rutgers map
#'
#' @param file_path the location of the rutgers map
#'
#' @return the rutgers map into a data table structure
#' @import data.table
#' @export
#'
get_rutgers_map <- function(file_path){
  chr_rutgers_map = suppressMessages(fread(file_path))
  # filter on the desired columns : Sex_averaged_start_map_position, Build37_start_physical_positio, Marker_name and the chromosome code
  return(data.frame(rsid = chr_rutgers_map$Marker_name,
                    position = chr_rutgers_map$Build37_start_physical_position,
                    cM = chr_rutgers_map$Sex_averaged_start_map_position))
}

#' Get 1000 genome interpolated map
#'
#' @description Get 1000 genome interpolated map
#'
#' @param file_path the location of the 1000 genome interpolated map
#'
#' @return the 1000 genome interpolated map into a data table structure
#' @import data.table
#' @export
#'
get_1000_genome_interpolated_map <- function(file_path){
  chr_1000_genome_interpolated_map = suppressMessages(fread(file_path))
  chr_1000_genome_interpolated_map_filtred = data.frame(rsid = chr_1000_genome_interpolated_map[, 1],
                                                        position = chr_1000_genome_interpolated_map[, 2],
                                                        cM = chr_1000_genome_interpolated_map[, 3])
  colnames(chr_1000_genome_interpolated_map_filtred) = c('rsid', 'position', 'cM')
  return(chr_1000_genome_interpolated_map_filtred)
}


#' Get 1000 genome map
#'
#' @description Get 1000 genome map
#'
#' @param file_path the location of the 1000 genome map
#'
#' @return the 1000 genome map into a data table structure
#' @import data.table
#' @export
#'
get_1000_genome_map <- function(file_path){
  chr_1000_genome_map = suppressMessages(fread(file_path))

  # replacing the () and / characters by _ in the columns name
  colnames(get_1000_genome_map) = c('position', 'COMBINED_rate_cM_Mb', 'Genetic_Map_cM')

  return(data.frame(position = chr_1000_genome_map$position,
                    cM = chr_1000_genome_map$Genetic_Map_cM))
}



#' Read the genetic map
#'
#' @description Read the genetic map
#'
#' @param genetic_map_dir A path to a dir containing the maps
#' @param chromosomes list of integer representing the chromosome that one want to read
#'
#' @return List representing the genetic map loaded.
#' @export
#'
get_genetic_map <- function(genetic_map_dir, chromosomes=1:23) {
  gen_map = list()
  for (chr in 1:23){
    # read the rutgers map
    chr_map = get_rutgers_map(sprintf('%s/RUMapv3_B137_chr%s.txt', genetic_map_dir, chr))

    # append to gen_map list
    gen_map[[chr]] = data.frame(cM=chr_map$cM, pos=chr_map$position, rsid=chr_map$rsid, chr=chr)
  }
  return(gen_map)
}


#' Get bim file
#'
#' @description Get bim file
#'
#' @param file_path A path to a file with tablulation as a delimiter
#'
#' @return The .bim file in a tibbles structure (data frame).
#' @import readr
#' @export
#'
get_bim_file <- function(file_path){
  return(read_delim(file_path, delim='\t', col_names=c('chromosome', 'snp', 'bp')))
}


#' @export
#'
getAnnotVariants <- function(obj)
  UseMethod("getAnnotVariants", obj)

#' @export
#'
getAnnotVariants.default <- function(obj) {
  stop("getAnnotVariants not defined on this object")
}

#' @export
#'
getAnnotVariants.vcf <- function(obj) {
  stop("getAnnotVariants not defined on this object")
}


#' @export
#'
getAnnotVariants.bgen <- function(obj) {
    return(obj[["annot_variants"]])
}



#' @export
#'
getInternalIID <- function(obj)
  UseMethod("getInternalIID", obj)

#' @export
#'
getInternalIID.default <- function(obj) {
  stop("getInternalIID not defined on this object")
}

#' @export
#'
getInternalIID.vcf <- function(obj) {
  stop("getInternalIID not defined on this object")
}


#' @export
#'
getInternalIID.bgen <- function(obj) {
    return(obj[["annot_internalIID"]])
}



#' Get bgi file
#'
#' @description Get bgi file
#'
#' @param file_path A path to a .bgi file
#'
#' @return data frame with the following columns: chromosome, position, rsid, number_of_alleles, allele1, allele2, file_start_position size_in_bytes
#' @import DBI
#' @export
#'
get_bgi_file <- function(file_path){
  # create a connection to the data base managment system
  con = (dbConnect(RSQLite::SQLite(), file_path))
  # read the variant table and store it as a data frame structure
  bgi_dataframe = data.frame(dbReadTable(con, "Variant"))
  # close the connection and frees ressources
  dbDisconnect(con)
  return(bgi_dataframe)
}

#' Get bgen file
#' @description Get bgen file
#'
#' @param file_path A path to a .bgen file
#' @param start the start genomic position
#' @param end the end genomic position
#' @param chromosome String. The chromosome code. '' by default
#' @param max_entries_per_sample An integer specifying the maximum number of probabilities expected per variant per sample.
#' This is used to set the third dimension of the data matrix returned. 4 by default.
#' @param samples A character vector specifying the IDs of samples to load data for.
#'
#' @return bgen file loaded in a bgen format
#' @import rbgen
#' @export
#'
get_bgen_file <- function(file_path, start, end, samples=samples, chromosome='', max_entries_per_sample=4){
  return(bgen.load(filename = file_path,
                   data.frame(chromosome=chromosome, start=start, end=end),
                   samples = samples,
                   max_entries_per_sample = max_entries_per_sample))
}


#' Write genetic map
#'
#' @description Write genetic map
#'
#' @param output A dir path where the map is saved
#' @param dataframe dataframe representing the augmented genetic map for one chromosome
#' @import utils
#'
#' @export
#'
write_genetic_map <- function(dataframe, output){
  write.table(dataframe, output, sep="\t", row.names=FALSE, quote=FALSE)
}


#' Get augmented genetic map
#'
#' @description Get augmented genetic map
#'
#' @param augmented_genetic_map_dir A path to the augmeneted genetic map dir
#' @param chromosomes A list of chromosomes that one want to read
#'
#' @return the augmented genetic map into a data table structure
#' @export
#'
get_augmented_genetic_map <- function(augmented_genetic_map_dir, chromosomes=1:22){
  augmented_genetic_map_df = c()
  for (chr in 1:22){
    augmented_chr = sprintf('%s/augmented_map_chr%s.txt', augmented_genetic_map_dir, chr)
    if(file.exists(augmented_chr)){augmented_genetic_map_df = rbind(augmented_genetic_map_df, read_delim(augmented_chr, delim='\t'))}
  }
  return(augmented_genetic_map_df)
}


#' Get blocs
#'
#' @description Get blocs
#'
#' @param blocs_dir A path to the blocs dir
#' @param chromosomes A list of chromosomes that one want to read
#' @import readr
#'
#' @return the blocs concatenated into a data table structure
#' @export
#'
get_blocs <- function(blocs_dir, chromosomes=1:22){
  blocs_df = c()
  for (chr in chromosomes){
    blocs_chr = sprintf('%s/blocs_chr%s.txt', blocs_dir, chr)
    print(blocs_chr)
    if(file.exists(blocs_chr)){blocs_df = rbind(blocs_df, read_delim(blocs_chr, delim='\t'))}
  }
  return(blocs_df)
}


#' Write blocs
#'
#' @description Write blocs
#'
#' @param dataframe dataframe representing the blocs created for one chromosome
#' @param output A dir path where the blocs are saved
#' @import utils
#'
#' @export
#'
write_blocs <- function(dataframe, output){
  write.table(dataframe, output, sep="\t", row.names=FALSE, quote=FALSE)
}


#' Save haplotypes
#'
#' @description Save haplotypes per chromosome. Each rows represent the subject with their IID as index.
#' Each column represent the haplotypes name that basicaly contain the follow information chromosome code, bloc start bp, end bloc bp and the haplotypes code
#' @param haplotype_combined haplotype dataframe. The rows correspond to the subject while the column correspond to the haplotypes name
#' @param chromosome chromosome code
#' @param output A dir path where the haplotypes are saved
#'
#' @return None
#' @import utils
#' @export
#'
save_haplotypes <- function(haplotype_combined, chromosome, output){

  # set the output path
  haplotype_combined_path = sprintf('%s/haplotypes_%s.tsv', output, chromosome)

  # remove NA in the column name added by cbind
  colnames(haplotype_combined) = vapply(strsplit(colnames(haplotype_combined),"[.]"), `[`, 2, FUN.VALUE=character(1))

  # save the haplotype as tsv file
  #write.table(haplotype_combined, haplotype_combined_path, sep="\t", row.names=TRUE, quote=FALSE)

  # save the haplotypes as .RData
  save(haplotype_combined, file=sprintf('%s/haplotypes_%s.RData', output, chromosome), compress=T)
}


#' Summary haplotypes test
#'
#' @description Filter on the results obtained and keep only the significant p values
#' @param test_path A dir path where the tests are saved
#' @param threshold threshold
#' @param verbose silent warning messages. FALSE by default.
#' @import utils
#'
#' @return None
#' @export
#'
summary_haplotypes_test <- function(test_path, threshold = 5e-6, verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

  # init
  test_possible = list('bloc_test_results', 'complete_test_results', 'single_test_results')

  # init the outputs data frames
  bloc_test_results     = data.frame()
  complete_test_results = data.frame()
  single_test_results   = data.frame()

  # read each test and concatenate it into one dataframe for all blocs and chromosomes
  chromosme_test_path = Sys.glob(file.path(test_path, '*'))

  # create a summary dir
  dir.create(sprintf('%s/summary', test_path))

  for(chromosome_path in chromosme_test_path){
    for(test in test_possible){
      unit_test_path = Sys.glob(file.path(sprintf('%s/%s', chromosome_path, test), '*'))
      for(unit_path in unit_test_path){
        if(test=='bloc_test_results'){bloc_test_results <- rbind(bloc_test_results, data.frame(read_tsv(unit_path)))}
        if(test=='complete_test_results'){complete_test_results <- rbind(complete_test_results, data.frame(read_tsv(unit_path)))}
        if(test=='single_test_results'){single_test_results <- rbind(single_test_results, data.frame(read_tsv(unit_path)))}
      }
    }
  }

  # filtre on the significant p values
  bloc_test_results = bloc_test_results[bloc_test_results$p_value < threshold, ]
  complete_test_results = complete_test_results[complete_test_results$p_value < threshold, ]
  single_test_results = single_test_results[single_test_results$p_value < threshold, ]

  # write the summary
  write.table(bloc_test_results, sprintf('%s/summary/bloc_test_results.tsv', test_path), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(complete_test_results, sprintf('%s/summary/complete_test_results.tsv', test_path), sep="\t", row.names=FALSE, quote=FALSE)
  write.table(single_test_results, sprintf('%s/summary/single_test_results.tsv', test_path), sep="\t", row.names=FALSE, quote=FALSE)
}


#' Download rutgers maps
#'
#' @description download rutgers maps using the following url : http://compgen.rutgers.edu/downloads/rutgers_map_v3.zip
#'
#' @return None
#' @export
download_rutgers_map <- function(){
  # dont use linux command
  # use the native R cmd instead

  # download the rutgers map
  system('wget http://compgen.rutgers.edu/downloads/rutgers_map_v3.zip')
  # unzip
  system('unzip rutgers_map_v3.zip')
  # remove the zip file
  system('rm rutgers_map_v3.zip')
}


#' get rutgers map in genMap format [bp,cM]
#'
#' @param genetic_map_path the location of the rutgers map
#' @param chr chromosome code
#' @return genetic_map : the genetetic mapin genMap format
#'
#' @import data.table
#'
#' @export
get_rutgers_genMap_SexAverage <- function(genetic_map_path, chr){
  ru_map=(suppressMessages(fread(sprintf('%s/RUMapv3_B137_chr%s.txt', genetic_map_path, chr))))
  return(data.frame(bp=ru_map$Build37_start_physical_position,
                    cM=ru_map$Sex_averaged_start_map_position))
}


#' Create a S3 object ready to be queried from a haps file
#'
#' @param bgen_filename : full path name to the bgen file of the phased data
#' @return phased_data_loader : the genetetic mapin genMap format
#'
#' @import rbgen
#'
#' @export
phased_data_loader.haps <- function(haps_filename) {
  # check the existence of haps_filename file
  # TODO

  # read 2 flavors of haps file with 1 or 2 cols describing the snps
  sep = " "
  hap_field_num = count.fields(haps_filename, sep=sep)[1]
  phased_data = read_table(haps_filename, col_names=FALSE)
  if ((hap_field_num%%2) == 0){
    phased_data = phased_data[-2]
  }
  samples_num = (length(colnames(phased_data)) - 5)/2
  tmp = sprintf("sample_%d", 0:(samples_num-1))
  new_col_names =  c(c('chrom', 'rsid', 'pos', 'allele_1', 'allele_2'),
                     unlist(lapply(tmp, function(s) sprintf("%s_strand%d", s, 1:2))))
  colnames(phased_data) <- new_col_names

  ret_obj <- list(phased_data=phased_data, is_phased=TRUE, full_fname_haps=haps_filename)
  class(ret_obj) <- c(class(ret_obj), "phased", "haps")

  return(ret_obj)
}




#' Create a S3 object ready to be queried from a bgen file
#'
#' @param bgen_filename : full path name to the bgen file of the phased data
#' @return phased_data_loader : the genetetic mapin genMap format
#'
#' @import rbgen
#'
#' @export
phased_data_loader.bgen <- function(bgen_filename) {
  # silent warning messages
  options(warn=-1)

  # check the existence of bgen.bgi file
  # TODO

  # get the annotation
  full_fname_bgi=sprintf("%s.bgi", bgen_filename)
  annot_variants = get_bgi_file(full_fname_bgi)

  # open and check that data are phased
  data = get_bgen_file(file_path = bgen_filename,
                          start = annot_variants$position[1],
                          end = annot_variants$position[1],
                          samples = c(),
                          chromosome = '',
                          max_entries_per_sample = 4)
  #  print(str(data))
  annot_internalIID <-data$samples
  
  # In ukb chromosome names is not in the bgen/bgi :degenerated FLAG
  chrom_name_degenerated = FALSE
  if (unique(annot_variants$chromosome) == "") {
      chrom_name_degenerated = TRUE
  }
      
  ret_obj = list(full_fname_bgen=bgen_filename,
                 is_phased=TRUE,
                 max_entries=4,
                 annot_internalIID=annot_internalIID,
                 annot_variants=annot_variants,
                 full_fname_bgi=full_fname_bgi,
                 chrom_name_degenerated=chrom_name_degenerated)

  # create S3 object
  class(ret_obj) <- c(class(ret_obj), "phased", "bgen")

  return(ret_obj)
}

