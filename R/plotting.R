#' Phenotype distribution
#'
#' @description Phenotype distribution
#'
#' @param phenotypes_path list of phenotypes path. The phenotype file should have as column the participant IID and the phenotype values
#' @param alpha transparency param for the entire plot
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return phenotype distribution
#' @import ggplot2
#' @import readr
#' @export
phenotype_distribution <- function(phenotypes_path, alpha=0.4, verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

  phenotypes_df =  data.frame()
  for(path in phenotypes_path){
    phenotype_file = data.frame(read_csv(path))
    phenotype_file$phenotype_name = strsplit(basename(path), "\\.")[[1]][1]
    phenotypes_df <- rbind(phenotypes_df, phenotype_file)
  }
  return(ggplot(phenotypes_df, aes(x=phenotypes_df$connectivity, fill=factor(phenotypes_df$phenotype_name))) +
         geom_density(alpha=alpha) +
         xlab(colnames(phenotypes_df)[2]))
}




#' Haplotype bloc distribution
#'
#' @description Haplotype bloc distribution
#'
#' @param df_blocs data frame structure containing all the blocs created
#' @param xlim range of the x axe
#' @param ylim range of the y axe
#' @param alpha transparency parameter for the entire plot
#' @param per_chromosome if the distribution are grouped by chromosome or not. By default FALSE.
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return haplotype bloc distribution
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @export
haplotype_bloc_distribution <- function(df_blocs, xlim=NULL, ylim=NULL, alpha=0.1, per_chromosome=FALSE, verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

  factor = df_blocs$delta
  if(per_chromosome){factor = df_blocs$chr}

  haplotype_bloc_distribution_plot <- ggplot(df_blocs, aes(x=df_blocs$to_cm-df_blocs$from_cm, fill=factor(factor))) +
                                       geom_density(alpha=alpha) +
                                       scale_x_continuous(trans = 'log10') +
                                       labs(fill = "delta (cM)") +
                                       xlab('lh(cM)')
  if(!is.null(xlim)){haplotype_bloc_distribution_plot <- haplotype_bloc_distribution_plot + xlim(xlim[1], xlim[2])}
  if(!is.null(ylim)){haplotype_bloc_distribution_plot <- haplotype_bloc_distribution_plot + ylim(ylim[1], ylim[2])}
  if(per_chromosome){haplotype_bloc_distribution_plot <- haplotype_bloc_distribution_plot + facet_grid(rows=vars(df_blocs$delta))}

  return(haplotype_bloc_distribution_plot)
}


#_________________________
# utils. Not in the right place
#_________________________

#' Save ggplot object
#'
#' @description Save a ggplot object as a png/pdf file
#'
#' @param plot the ggplot object
#' @param output the path indicating where to save the plot
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return None
#' @export
save_plot <- function(plot, output, verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}
  # save the plot
  ggsave(output)
}
