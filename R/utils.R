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
  # the actual name of the column is Genetic_Map(cM). You need to rename chr_1000_genome_map before construcing the dataframe
  return(data.frame(position = chr_1000_genome_map$position,
                    cM = chr_1000_genome_map$Genetic_Map(cM)))
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
  write.table(dataframe, output, sep="\t", row.names=FALSE)
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
  for (chr in 1:22){
    blocs_chr = sprintf('%s/blocs_chr%s.txt', blocs_dir, chr)
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
  write.table(dataframe, output, sep="\t", row.names=FALSE)
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
