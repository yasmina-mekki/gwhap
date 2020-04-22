# This class was created followingthe principled description from:
# https://www.cyclismo.org/tutorial/R/s4Classes.html#creating-an-s4-class
# Documentation and NAMESPACE : http://r-pkgs.had.co.nz/man.html#man-classes
# See examples : https://github.com/ropensci/datapack/blob/master/R/DataObject.R

# Create the base Genetic_Map S4 class
#
# This is used to represent the Genetic_Map.

#' An S4 class to represent a Genetic_Map.
#'
#' @name Genetic_Map
#' @rdname Genetic_Map
#' @slot filepaths A filepath character vector
#' @slot encodings A Named list c('cM', 'position', 'chr')
#' @slot gmapData A data.frame for the read genetic map data
#' @export Genetic_Map
Genetic_Map <- setClass(
    # Set the name for the class
    "Genetic_Map",

    # Define the slots
    slots = c(
            filepaths = "character",
            encodings = "list",
            gmapData = "data.frame"
            ),

    # Set the default values for the slots. (optional)
    prototype=list(
            filepaths = c(),
            encodings   = list("chr"=NULL, "cM"=NULL, "position"=NULL),
            gmapData = data.frame()
            ),

    # Make a function that can test to see if the data is consistent.
    # This is not called if you have an initialize function defined!
    validity=function(object)
    {
        # filepaths validity
        if(!is.character(object@filepaths)) {
            return("The filepaths should be a character vector")
        }
        # encodings validity
        if(!is.list(object@encodings)) {
            return("The encodings should be a list")
        }
        if(!all(names(object@encodings) %in% c("chr", "cM", "position"))) {
            print(names(object@encodings))
            return("The encodings list should contain chr, cM and position keys")
        }
        # specific encodings@chr validity checks
        if (!(is.list(object@encodings[["chr"]]) | is.character(object@encodings[["chr"]]))) {
            return("The encodings@chr should be list or a charcacter")
        }
        if (is.list(object@encodings[["chr"]])) {
            if (!all(names(object@encodings[["chr"]]) %in% filepaths)) {
                return("The encodings@chr list does not match wirh filepaths vector")
            }
            if (!all(sapply(encodings[["chr"]], function(i) is.numeric(i)))) {
                return("The encodings@chr should contain only numeric as chromosome reference")
            }
        }
        return(TRUE)
    }
)

# create a method to read the Genetic_Map map
# first define hidden ancillary functions
.read_gmap_no_internal_chrom_info <- function(object) {
    outdf = data.frame()
    for (key in names(object@encodings$chr)) {
        # key in this map (filename <-> #chr)
        chr_num = object@encodings$chr[[key]]
        df = fread(key)
        # get the sepecified column and fix names - robust code to keep
        col_tokeep = names(encodings)[names(encodings) %in% c('cM', 'position')]
        col_toextract = sapply(col_tokeep, function(s) object@encodings[[s]])
        df = df[, ..col_toextract]
        colnames(df) = col_tokeep
        # set the chr number
        df[, 'chr'] = chr_num
        # collate w/ existing
        outdf = rbind(outdf, df)
    }

    return(outdf)
}


# now define the exposet object Genetic_Map methods
#' readData perform the actual reading of the genetic maps data based on the encodings
#' @name readData
#' 
#' @aliases readData
#' @export
setGeneric(name="readData",
    def=function(object)
    {
       standardGeneric("readData")
    }
)

#' @rdname readData
#' @return data for map based on mutiple chromfiles
#' @aliases getData
#' @examples
#' TOBEFIXED
setMethod(f="readData",
    signature="Genetic_Map",
    definition=function(object)
    {
        if (is.list(object@encodings[["chr"]])) {
        # chr information not in the files
            object@gmapData = .read_gmap_no_internal_chrom_info(object)
        }
        else {
        # chr information is in the files
            object@gmapData = data.frame()
        }
        return(object)
    }
)

