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
