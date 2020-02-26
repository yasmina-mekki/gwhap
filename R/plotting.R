#' phenotype distribution
#'
#' @param phenotypes_path list of phenotypes path
#' @param alpha alpha
#'
#' @return ggplot
#' @import ggplot2
#' @import readr
#' @export
phenotype_distribution <- function(phenotypes_path, alpha=0.4){

  phenotypes_df =  data.frame()
  for(path in phenotypes_path){
    phenotype_file = data.frame(read_csv(path))
    phenotype_file$phenotype_name = strsplit(basename(path), "\\.")[[1]][1]
    phenotypes_df <- rbind(phenotypes_df, phenotype_file)
  }
  return(ggplot(phenotypes_df, aes(x=phenotypes_df$connectivity, fill=phenotype_file$phenotype_name)) + geom_density(alpha = alpha))
}




#' haplotype block distribution
#'
#' @param df_blocs blocs
#' @param xlim range of the x axe
#' @param ylim range of the y axe
#' @param alpha alpha
#' @param per_chromosome if the distribution are grouped by chromosome or not. By default FALSE.
#'
#' @return ggplot
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @export
haplotype_block_distribution <- function(df_blocs, xlim=NULL, ylim=NULL, alpha=0.1, per_chromosome=FALSE){
  # TO DO
  # to_cm-from_cm same as win_size ?

  factor = df_blocs$delta
  if(per_chromosome){factor = df_blocs$chr}

  haplotype_block_distribution_plot <- ggplot(df_blocs%>% mutate(lh=df_blocs$to_cm-df_blocs$from_cm), aes(x=lh, fill=factor(factor))) +
    geom_density(alpha=alpha) + labs(fill = "delta (cM)")
  if(!is.null(xlim)){haplotype_block_distribution_plot <- haplotype_block_distribution_plot + xlim(xlim[1], xlim[2])}
  if(!is.null(ylim)){haplotype_block_distribution_plot <- haplotype_block_distribution_plot + ylim(ylim[1], ylim[2])}
  if(per_chromosome){haplotype_block_distribution_plot <- haplotype_block_distribution_plot + facet_grid(rows=vars(df_blocs$delta),
                                                                                                         scales="free")}

  return(haplotype_block_distribution_plot)
}


#' Genome coverage of haplotype blocs
#'
#' @param df_blocs data.frame
#' @param colors colors (int)
#'
#' @return karyotype plot object
#' @import karyoploteR
#' @importFrom grDevices rainbow
#' @importFrom dplyr group_by
#' @export
karyotype_plot <- function(df_blocs, colors=NULL){
  # TO DO
  # handle colors params
  # add scale bar

  # params blocs setting
  df_blocs$lbp = df_blocs$to_bp - df_blocs$from_bp
  df_blocs$mb  = round(df_blocs$from_bp/1e6)
  sum_l_bp_chr = df_blocs %>% group_by(df_blocs$chr, df_blocs$mb) %>% summarise(coverage=sum(lbp)/1e6)

  # params plot setting
  pp <- getDefaultPlotParams(plot.type = 1)
  pp$data1height = 200

  # init karyotype object with params plot chosen above
  karyotype_plot_obj <- plotKaryotype(chromosomes=c("autosomal"), plot.params = pp)

  kpHeatmap(karyoplot = karyotype_plot_obj,
            chr       = sum_l_bp_chr$chr,
            x0        = sum_l_bp_chr$mb*1e6,
            x1        = (sum_l_bp_chr$mb+1)*1e6,
            y         = round(sum_l_bp_chr$coverage*100),
            colors    = rainbow(200)[round(sum_l_bp_chr$coverage*100)]) #rainbow(200, start=0.2, end=0.7)
  # add base numbers
  kpAddBaseNumbers(karyotype_plot_obj)

  return(karyotype_plot_obj)
}


#_________________________
# utils. Not in the right place
#_________________________

#' Save plot
#'
#' @param plot the ggplot object
#' @param output the path indicating where to save the plot
#'
#' @return None
#' @export
save_plot <- function(plot, output){
  ggsave(output)
}
