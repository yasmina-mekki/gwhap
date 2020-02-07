#' phenotype distribution
#'
#' @param phenotypes_path list of phenotypes path
#'
#' @return ggplot
#' @import ggplot2
#' @import readr
#' @export
phenotype_distribution <- function(phenotypes_path){

  phenotypes_df =  data.frame()
  for(path in phenotypes_path){
    phenotype_file = data.frame(read_csv(path))
    phenotype_file$phenotype_name = strsplit(basename(path), "\\.")[[1]][1]
    phenotypes_df <- rbind(phenotypes_df, phenotype_file)
  }
  return(ggplot(phenotypes_df, aes(x=connectivity, fill=phenotype_name)) + geom_density(alpha = 0.4))
}


#' haplotype block distribution
#'
#' @return ggplot
#' @export
haplotype_block_distribution <- function(){

}


#_________________________
# utils. Not in the right place
#_________________________

#' save plot
#'
#' @param plot the plot to save
#'
#' @return NONE
#' @export
save_plot <- function(plot){
  # ggsave("myplot.pdf")
}
