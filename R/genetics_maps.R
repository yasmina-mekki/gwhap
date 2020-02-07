#_________________________
# Download part
#_________________________

#' download rutgers maps using the following url : http://compgen.rutgers.edu/downloads/rutgers_map_v3.zip
#'
#' @return None
#' @export
download_rutgers_map <- function(){
  # dont use linux command
  # instead use the native R cmd

  # download the rutgers map
  system('wget http://compgen.rutgers.edu/downloads/rutgers_map_v3.zip')
  # unzip
  system('unzip rutgers_map_v3.zip')
  # remove the zip file
  system('rm rutgers_map_v3.zip')
}

#_________________________
# Getteur part
#_________________________

#' get rutgers map
#'
#' @param file_path the location of the rutgers map
#'
#' @return the rutgers map in a data.table structure
#' @export
get_rutgers_map <- function(file_path){
  return(suppressMessages(fread(file_path)))
}
