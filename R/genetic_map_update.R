#' update the genetic map
#'
#' @param genetic_map_path i.e. rutgers map
#' @param to_be_updated_genetic_map_path the genetic map that you want to update
#' @param delim  character used to separate fields within a record. By default tab separated values is used
#'
#' @return a genetic map updated
#' @import data.table
#' @import readr
#' @export
read.augmented_genetic_map <- function(genetic_map_dir, snp_physical_positions){

  snp_list = read_delim(snp_physical_positions, delim='\t', col_names=c('chromosome', 'snp', 'bp'))
  gen_map_updated = list()
  for (chr in unique(snp_list$chromosome)){
    gen_map_updated[[chr]] = data.frame()
  }

  for (chr in unique(snp_list[,1])){
    chr_map = get_rutgers_map(sprintf('%s/RUMapv3_B137_chr%s.txt', genetic_map_dir, chr)) #use Sys.glob instead of specifying the path

    interp_in = chr_map$Build37_start_physical_position > (min(snp_list$bp[snp_list$chromosome == chr])-1) &
      chr_map$Build37_start_physical_position < (max(snp_list$bp[snp_list$chromosome == chr]+1))

    gen_map_approx.fun <- stats::approxfun(chr_map$Build37_start_physical_position[interp_in],
                                    chr_map$Sex_averaged_start_map_position[interp_in],
                                    ties="ordered")

    snp_interp = gen_map_approx.fun(snp_list$bp[snp_list$chromosome == chr])
    names(snp_interp) = snp_list$snp

    # update the genetic map
    gen_map_updated[[chr]]  = data.frame(cM=snp_interp, pos=snp_list$bp, rsid=snp_list$snp, chr=chr)
  }


  return(gen_map_updated)
}


read.genetic_map <- function(genetic_map_dir, chromosomes=1:23) {

  gen_map = list()
  for (chr in chromosomes){
    gen_map[[chr]] = data.frame()
  }


  for (chr in 1:23){
    chr_map = get_rutgers_map(sprintf('%s/RUMapv3_B137_chr%s.txt', genetic_map_dir, chr))#use Sys.glob instead of specifying the path

    gen_map[[chr]] = data.frame(cM=chr_map$Sex_averaged_start_map_position,
                                  pos=chr_map$Build37_start_physical_position,
                                  rsid=chr_map$Marker_name,
                                  chr=chr)
  }

  return(gen_map)
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
