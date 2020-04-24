#' Create an augmented genetic map
#'
#' @description  Create an augmented genetic map
#'
#' @param snp_bucket A Snp_Bucket object containing the snp physical
#' positions to interpolate from the reference map. Either a bim or bgi file
#' @param genetic_map_dir i.e. rutgers map
#' @param genetic_map A Genetic_Map object containing the reference 
#' on which the interpolation is based genetic map
#' @param save_genetic_map Boolean. specify if the augmented genetic map should be saved as a txt file. FALSE by default.
#' @param outputfile A file path where the augmented genetic map will be saved as Genetic_Maps instance in rds.
#' @param verbose A Boolean for warning messages. FALSE by default.
#'
#' @return if save_genetic_map == TRUE, the augmented genetic map as Genetic_Maps instance in rds
#' The format of the output would be one txt file per chromosome.
#' Each txt file would have the following information:
#' each line represent a SNP. The columns represent the centimorgan information, snp's position, its rs id and chromosome code.
#' Otherwise, it will return a list of dataframe.
#' @import data.table
#' @import readr
#' @import tools
#' @importFrom stats na.omit
#' @export
create_augmented_genetic_map <- function(snp_bucket, 
                                    genetic_map=NULL,
                                    save_genetic_map=FALSE, 
                                    outputfile='', verbose=FALSE){

    # silent warning messages
    if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

    # create the augmented genetic map
    gen_map_augmented_data = list()

    short_list_chr = intersect(unique(names(snp_bucket@bucketData)), unique(names(genetic_map@gmapData)))
    for (chr in short_list_chr) {
        # get the NAMED list (chr <-> data.frame
        chr_map = genetic_map@gmapData[[chr]]
        chr_snp_list = snp_bucket@bucketData[[chr]]
    
        # get interval of position defined by min-1 and max+1 neighbour variant in the genetic map
        interp_in = (chr_map$position > (min(chr_snp_list$position)-1)) & (chr_map$position < (max(chr_snp_list$position)+1))

        # initialize the linear interpolator
        gen_map_approx.fun <- stats::approxfun(chr_map$position[interp_in],
                                               chr_map$cM[interp_in],
                                               ties="ordered")
        # filter on the position, rsid and snp interpolation for the current chromosome
        snp_interp = gen_map_approx.fun(chr_snp_list$position)
        position   = chr_snp_list$position
        rsid       = chr_snp_list$snp

        # check if NA values are present
        gen_map_augmented_chr = data.frame(cM=snp_interp, pos=position, rsid=rsid, chr=chr)
        if (nrow(gen_map_augmented_chr[is.na(gen_map_augmented_chr), ]) != 0){
          warning(sprintf('Number of NA in chromosome %s : %s. Progress towards removing ...', chr, nrow(gen_map_augmented_chr[is.na(gen_map_augmented_chr), ])))
        }

        # remove the NAs if present and store the augmented genetic map for each chromosome
        gen_map_augmented_data[[chr]]  = na.omit(gen_map_augmented_chr)
        
    }
    
    # Embed the gen_map_augmented_data in an S4 instance of Genetic_Data
    gen_map_augmented = Genetic_Map(
                filepaths="",
                encodings=list("chr"="", "cM"="", "position"="", "format"=""),
                gmapData=gen_map_augmented_data)

    # if save_genetic_map equal to TRUE, write the augmented genetic map on txt files
    if (save_genetic_map) {
        saveData(gen_map_augmented, file=outputfile)
        return(invisible())
    }
    # if save_genetic_map is FALSE, return the augmented genetic map as a Genetic_Map instance
    return(gen_map_augmented)
}


##' Use and interpolate genetic map of a chromosome
##'
##' @param genetic_map  reference genetic map of a chromosome [BP,cM]
##' @param snp_list SNP position list to interpolate [BP]
##'
##' @return gen_map_chr_interp : genetic map for the SNP list : [BP,cM]
##' @export
#genetic_map_interp_chr <- function(genetic_map, snp_list){

#    # builiding the interpolation model using all reference positions in chromosome chr
#    gen_map_approx.fun <- stats::approxfun(genetic_map[,1],
#                                           genetic_map[,2],
#                                           ties="ordered")
#    # snps to interpolate in the chromosome chr
#    snp_interp = gen_map_approx.fun(snp_list)

#    # update the genetic map
#    gen_map_chr_interp = data.frame(pos=snp_list,
#                                cM=snp_interp)


#  return(gen_map_chr_interp)
#}


##' Use and interpolate of all chromosomes using a reference map  object
##'
##' @param genetic_map_all  reference genetic map of all chromosomes [chr,BP,cM]
##' @param snp_list SNP position list to interpolate in all chromosomes [chr,BP]
##'
##' @return gen_map_allchr_interp : genetic map for the SNP list : [chr,BP,cM]
##' @export
#genetic_map_interp <- function(genetic_map_all, snp_list){
#  gen_map_allchr_interp=c()
#  #loop over all chr in the snp list
#  for (chr in unique(snp_list[,1])){
#    map_in=genetic_map_all[,1]==chr
#    genetic_map=genetic_map_all[map_in,c(2,3)]
#    snps_in=snp_list[,1]==chr
#    snps_to_interp= snp_list[snps_in,2]

#    # adding interpolation for chomosome chr to resulting map gen_map_allchr_interp
#    gen_map_allchr_interp=rbind(gen_map_allchr_interp,
#                                data.frame(chr=chr,
#                                           genetic_map_interp_chr(genetic_map = genetic_map,
#                                                                  snp_list = snps_to_interp)))

#  }
#  return(gen_map_allchr_interp)

#}
