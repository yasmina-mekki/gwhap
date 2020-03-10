#' create an augmented genetic map
#'
#' @param snp_physical_positions the snp physical position map. Either a bim or bgi file
#' @param genetic_map_dir i.e. rutgers map
#' @param save_genetic_map Boolean. specify if the augmented genetic map should be saved. FALSE by default.
#' @param map_name the genetic map name. Two maps name are available : rutgers and 1000_genome_interpolated. rutgers by default.
#' @param output A path to a file where the augmented genetic map will be saved.
#'
#' @return if save_genetic_map == TRUE, then the augmented genetic map is saved. Otherwise, it will return a list of dataframe.
#' @import data.table
#' @import readr
#' @import tools
#' @export
create_augmented_genetic_map <- function(snp_physical_positions, genetic_map_dir, map_name='rutgers', save_genetic_map=FALSE, output=''){

  # read the snp physical positions

  # read the one bim file for the whole autosome
  if('bim'==file_ext(snp_physical_positions)){snp_list = get_bim_file(snp_physical_positions)}

  # one bgi file per chromosome. read each one of them and concatenate them into one dataframe
  else{
    snp_list = c()
    for (chr in 1:23){
      bgi_file_path = sprintf('%schr%s_v2.bgen.bgi', snp_physical_positions, chr)
      # check if file path exist
      if(file.exists(bgi_file_path)){
        # get bgi file
        tmp = get_bgi_file(bgi_file_path)
        # add chromosome code
        tmp$chromosome = chr
        # change column name for correspondance with bim file
        colnames(tmp) <- c('chromosome', 'bp', 'snp', 'number_of_alleles', 'allele1', 'allele2', 'file_start_position', 'size_in_bytes')
        # concatenation
        snp_list  = rbind(snp_list, tmp)
      }
    }
  }

  # create the augmented genetic map
  gen_map_updated = list()
  for (chr in unique(snp_list$chromosome)){
    
    # read the genetic map. # use Sys.glob instead of specifying the path
    if(map_name = 'rutgers'){chr_map = get_rutgers_map(sprintf('%s/RUMapv3_B137_chr%s.txt', genetic_map_dir, chr))}
    if(map_name = '1000_genome_interpolated'){chr_map = get_1000_genome_interpolated_map(sprintf('%s/chr%s.interpolated_genetic_map.gz', genetic_map_dir, chr))}
    if(map_name = '1000_genome'){chr_map = get_1000_genome_map(sprintf('%s/genetic_map_chr%s_combined_b37.txt', genetic_map_dir, chr))} 
    
    # get interval of position defined by min-1 and max+1 neighbour variant in the genetic map
    interp_in = chr_map$position > (min(snp_list$bp[snp_list$chromosome == chr])-1) & chr_map$position < (max(snp_list$bp[snp_list$chromosome == chr]+1))

    # perform the linear interpolation
    gen_map_approx.fun <- stats::approxfun(chr_map$position[interp_in],
                                           chr_map$cM[interp_in],
                                           ties="ordered")

    # filter on the position, rsid and snp interpolation for the current chromosome
    snp_interp = gen_map_approx.fun(snp_list$bp[snp_list$chromosome == chr])
    position   = snp_list$bp[snp_list$chromosome == chr]
    rsid       = snp_list$snp[snp_list$chromosome == chr]

    # store the augmented genetic map for each chromosome
    gen_map_updated[[chr]]  = data.frame(cM=snp_interp, pos=position, rsid=rsid, chr=chr)
  }

  # if save_genetic_map, write the augmented genetic map
  if (save_genetic_map) {
    for (chr in unique(snp_list$chromosome)){
      file_path = sprintf('%s/augmented_map_chr%s.txt', output, chr)
      write_genetic_map(gen_map_updated[[chr]], file_path)
    }
    return(0)
  }
  # if not save_genetic_map, return the augmented genetic map as a list of dataframe
  return(gen_map_updated)
}


#' use and interpolate genetic map of a chromosome
#'
#' @param genetic_map  reference genetic map of a chromosome [BP,cM]
#' @param snp_list SNP position list to interpolate [BP]
#'
#' @return gen_map_chr_interp : genetic map for the SNP list : [BP,cM]
#' @export
genetic_map_interp_chr <- function(genetic_map, snp_list){

    #builiding the interpolation model using all reference positions in chromosome chr
    gen_map_approx.fun <- stats::approxfun(genetic_map[,1],
                                           genetic_map[,2],
                                           ties="ordered")
    #snps to interpolate in the chromosome chr
    snp_interp = gen_map_approx.fun(snp_list)

    # update the genetic map
    gen_map_chr_interp = data.frame(pos=snp_list,
                                cM=snp_interp)


  return(gen_map_chr_interp)
}


#' use and interpolate of all chromosomes using a reference map  object
#'
#' @param genetic_map_all  reference genetic map of all chromosomes [chr,BP,cM]
#' @param snp_list SNP position list to interpolate in all chromosomes [chr,BP]
#'
#' @return gen_map_allchr_interp : genetic map for the SNP list : [chr,BP,cM]
#' @export
genetic_map_interp <- function(genetic_map_all, snp_list){
  gen_map_allchr_interp=c()
  #loop over all chr in the snp list
  for (chr in unique(snp_list[,1])){
    map_in=genetic_map_all[,1]==chr
    genetic_map=genetic_map_all[map_in,c(2,3)]
    snps_in=snp_list[,1]==chr
    snps_to_interp= snp_list[snps_in,2]

    # adding interpolation for chomosome chr to resulting map gen_map_allchr_interp
    gen_map_allchr_interp=rbind(gen_map_allchr_interp,
                                data.frame(chr=chr,
                                           genetic_map_interp_chr(genetic_map = genetic_map,
                                                                  snp_list = snps_to_interp)))

  }
  return(gen_map_allchr_interp)

}
