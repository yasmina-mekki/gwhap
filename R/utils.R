#__________________________________________________________________________________________________
# TO DO :
# get_rutgers_map NOT necessary ... have been replaced by get_genetic_map. You should remove it.
# merge write_genetic_map and write_blocs functions
# merge get_augmented_genetic_map and get_blocs functions
#__________________________________________________________________________________________________

#' get rutgers map
#'
#' @param file_path the location of the rutgers map
#'
#' @return the rutgers map in a data.table structure
#' @import readr
#' @export
get_rutgers_map <- function(file_path){
  return(suppressMessages(fread(file_path)))
}

#' read the genetic map
#'
#' @param genetic_map_dir A path to a dir containing the maps
#' @param chromosomes list of integer representing the chromosome that one want to read
#'
#' @return
#' @export
#'
#' @examples
get_genetic_map <- function(genetic_map_dir, chromosomes=1:23) {

  gen_map = list()
  for (chr in 1:23){
    # read the rutgers map
    chr_map = get_rutgers_map(sprintf('%s/RUMapv3_B137_chr%s.txt', genetic_map_dir, chr))#use Sys.glob instead of specifying the path

    # filter on the desired columns : Sex_averaged_start_map_position, Build37_start_physical_positio, Marker_name and the chromosome code
    gen_map[[chr]] = data.frame(cM=chr_map$Sex_averaged_start_map_position,
                                pos=chr_map$Build37_start_physical_position,
                                rsid=chr_map$Marker_name,
                                chr=chr)
  }

  return(gen_map)
}


#' get .bim file
#'
#' @param file_path A path to a file with tablulation as a delimiter
#'
#' @return The .bim file in a tibbles structure (data.frame).
#' @import readr
#' @export
#'
get_bim_file <- function(file_path){
  return(read_delim(file_path, delim='\t', col_names=c('chromosome', 'snp', 'bp')))
}


#' get .bgi file
#'
#' @param file_path A path to a .bgi file
#'
#' @return data.frame with the following columns (chromosome, position, rsid, number_of_alleles, allele1, allele2, file_start_position size_in_bytes)
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

#' get bgen file
#' Need to handle when samples are not specifyed ...
#'
#' @param file_path A path to a .bgen file
#' @param start the start genomic position
#' @param end the end genomic position
#' @param chromosome String. The chromosome code. '' by default
#' @param max_entries_per_sample An integer specifying the maximum number of probabilities expected per variant per sample.
#' This is used to set the third dimension of the data matrix returned. 4 by default.
#' @param samples A character vector specifying the IDs of samples to load data for.
#'
#' @return
#' @import rbgen
#' @export
#'
get_bgen_file <- function(file_path, start, end, samples=samples, chromosome='', max_entries_per_sample=4){
  return(bgen.load(filename=bgnfile,
                   data.frame(chromosome=chromosome, start=start, end=end),
                   samples = samples,
                   max_entries_per_sample=max_entries_per_sample))
}


#' write genetic map
#'
#' @param output A path where a dir is created and that will contain the maps
#' @param dataframe dataframe representing the augmented genetic map for one chromosome
#' @import utils
#'
#' @export
#'
write_genetic_map <- function(dataframe, output){
  write.table(dataframe, output, sep="\t", row.names=FALSE)
}


#' get augmented genetic map
#'
#' @param augmented_genetic_map_dir A path to the augmeneted genetic map dir
#' @param chromosomes A list of chromosomes that one want to read
#'
#' @return Data.table with all the augmented genetic map concatenated
#' @export
#'
get_augmented_genetic_map <- function(augmented_genetic_map_dir, chromosomes=1:23){
  augmented_genetic_map_df = c()
  for (chr in 1:22){
    augmented_chr = sprintf('%s/augmented_map_chr%s.txt', augmented_genetic_map_dir, chr)
    if(file.exists(augmented_chr)){augmented_genetic_map_df = rbind(augmented_genetic_map_df, read_delim(augmented_chr, delim='\t'))}
  }
  return(augmented_genetic_map_df)
}


#' get blocs
#'
#' @param blocs_dir A path to the blocs dir
#' @param chromosomes  A list of chromosomes that one want to read
#'
#' @return Data.table with all the blocs concatenated
#' @export
#'
get_blocs <- function(blocs_dir, chromosomes=1:22){
  blocs_df = c()
  for (chr in 1:22){
    blocs_chr = sprintf('%s/blocs_chr%s.txt', blocs_dir, chr)
    if(file.exists(blocs_chr)){blocs_df = rbind(blocs_df, read_delim(blocs_chr, delim='\t'))}
  }
  return(blocs_df)
}


#' write blocs
#'
#' @param dataframe dataframe representing the blocs created for one chromosome
#' @param output A path where a dir is created and that will contain the blocs
#' @import utils
#'
#' @export
#'
write_blocs <- function(dataframe, output){
  write.table(dataframe, output, sep="\t", row.names=FALSE)
}
