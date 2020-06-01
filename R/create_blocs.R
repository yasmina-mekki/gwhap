#' Create blocs
#'
#' @description Create blocs
#'
#' @param genetics_map genetic map S4 objects with a slot gmapData
#' @param delta delta
#' @param save_blocs Boolean. specify if the created blocs should be saved. FALSE by default.
#' @param output A path to a file where the created blocs will be saved.
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return if save_blocs == TRUE, then the created blocs are saved.
#' The format of the output would be one txt file per chromosome.
#' Each txt file would have the following information:
#' each line represent a bloc. The columns represent the start and end window of a bloc, a window size, the start and end bloc position,
#' the chromosome code and the delta.
#' Otherwise, it will return a list of dataframe representing each chromosome's blocs.
#' @export
create_blocs <- function(genetics_map, delta=1e-3, save_blocs=FALSE, output='', verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

  if (class(interpolated_map) != "Genetic_Map") {  
      return("The first parameter is not a Genetic_Map S4 object")
  }

  genetics_mapData = genetics_map@gmapData
  bloc_df = list()

#  for (chr in unique(genetics_mapData$chr)){
  for (chr in names(genetics_mapData)) {

    idx_chr  = genetics_mapData[[chr]]$chr == chr
    start_cm = genetics_mapData[[chr]]$cM[idx_chr]
    position = genetics_mapData[[chr]]$pos[idx_chr]
    ft       = position[c(1, length(position))]
    from_bp  = ft[1];
    to_bp    = ft[2];
    ft       = position

    cut = (which(diff(start_cm) > delta))

    start_w = c(1, cut+1)
    end_w   = c(cut, length(start_cm))
    size_ok = start_w < end_w
    start_w = start_w[size_ok]
    end_w   = end_w[size_ok]


    from_cm  = start_cm[start_w]
    to_cm    = start_cm[end_w]
    from_bp  = position[start_w]
    to_bp    = position[end_w]
    win_size = apply(data.frame(from_cm, to_cm), 1, diff)

    bloc_df[[chr]]  = data.frame(start_w, end_w, win_size, from_bp, to_bp, from_cm, to_cm, chr, delta)
  }

  if (save_blocs) {
#    for (chr in unique(genetics_mapData[[chr]]$chr)){
    for (chr in names(genetics_mapData)) {
      file_path = sprintf('%s/blocs_chr%s.txt', output, chr)
      write_blocs(bloc_df[[chr]], file_path)
    }
    return(0)
  }
  # if not save_blocs, return the blocs as a list of dataframe for each chromosome
  return(bloc_df)
}
