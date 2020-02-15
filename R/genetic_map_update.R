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
genetic_map_update <- function(genetic_map_path, to_be_updated_genetic_map_path, delim='\t'){

  snp_list = read_delim(to_be_updated_genetic_map_path, delim='\t', col_names=c('chromosome', 'snp', 'bp'))
  gen_map_updated = data.frame()

  for (chr in unique(snp_list[,1])){
    chr_map = get_rutgers_map(sprintf('%s/RUMapv3_B137_chr%s.txt', genetic_map_path, chr)) #use Sys.glob instead of specifying the path

    interp_in = chr_map$Build37_start_physical_position > (min(snp_list$bp[snp_list$chromosome == chr])-1) &
      chr_map$Build37_start_physical_position < (max(snp_list$bp[snp_list$chromosome == chr]+1))

    gen_map_approx.fun <- approxfun(chr_map$Build37_start_physical_position[interp_in],
                                    chr_map$Sex_averaged_start_map_position[interp_in],
                                    ties="ordered")

    snp_interp = gen_map_approx.fun(snp_list$bp[snp_list$chromosome == chr])
    names(snp_interp) = snp_list$snp

    # update the genetic map
    gen_map_updated = rbind(gen_map_updated, data.frame(cM=snp_interp, pos=snp_list$bp, rsid=snp_list$snp, chr=chr))
  }

  return(gen_map_updated)
}
