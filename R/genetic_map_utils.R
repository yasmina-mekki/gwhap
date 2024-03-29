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
#' @slot gmapData A list of data.frame to containe the genetic map data.
#' This is a named list with key like 'chr1', 'chr2', 'chrX'. Each 
#' data.frame has two columns : 'cM' and 'position'
#' @export Genetic_Map
Genetic_Map <- setClass(
    # Set the name for the class
    "Genetic_Map",

    # Define the slots
    slots = c(
            filepaths = "character",
            encodings = "list",
            gmapData = "list"
            ),

    # Set the default values for the slots. (optional)
    prototype=list(
            filepaths = c(),
            encodings   = list("chr"=NULL, "cM"=NULL, "position"=NULL, "format"=NULL),
            gmapData = list(data.frame())
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
        if(!all(names(object@encodings) %in% c("chr", "cM", "position", "format"))) {
            print(names(object@encodings))
            return("The encodings list should contain chr, cM, position and format keys")
        }
        # specific encodings@format validity checks
        if (!(object@encodings[["format"]] %in% c("", "table", "bgen", "rds"))) {
            return("The encodings$format should be table, rds or bgen")
        }
        # specific encodings@chr validity checks
        if (!(is.list(object@encodings[["chr"]]) | is.character(object@encodings[["chr"]]))) {
            return("The encodings$chr should be list or a charcacter")
        }
        if (is.list(object@encodings[["chr"]])) {
            if (!all(names(object@encodings[["chr"]]) %in% filepaths)) {
                return("The encodings$chr list does not match wirh filepaths vector")
            }
            if (!all(sapply(encodings[["chr"]], function(i) is.numeric(i)))) {
                return("The encodings$chr should contain only numeric as chromosome reference")
            }
        }
        if (!is.list(object@gmapData)) {
            return("The gmapData should be a list of data.frame")
        }
        if (!is.data.frame(object@gmapData[[1]])) {
            return("The gmapData should be a list of data.frame")
        }
        return(TRUE)
    }
)

# create a method to read the Genetic_Map map
# first define hidden ancillary functions
.read_gmap_no_internal_chrom_info <- function(object) {
    outdf_list = list()
    for (key in names(object@encodings$chr)) {
        # key in this map (filename <-> #chr)
        chr_num = object@encodings$chr[[key]]
        df = fread(key)
        # get the sepecified column and fix names - robust code to keep
        col_tokeep = names(object@encodings)[names(object@encodings) %in% c('cM', 'position')]
        col_toextract = sapply(col_tokeep, function(s) object@encodings[[s]])
        df = df[, ..col_toextract]
        colnames(df) = col_tokeep
        # append the dataframe in the existing list w/ named chrom
        outdf_list[[sprintf('chr%d',chr_num)]] = df
    }

    return(outdf_list)
}

.read_gmap_rds <- function(object) {
    # only one filepath is expected in object@filepaths
    # the rds file read is expected to contain one instance of a 
    # instance of Genetic_Map
    # 
    outdf_list = readRDS(object@filepaths)@gmapData

    return( outdf_list )
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
#' @aliases readData
#' @examples
#' TOBEFIXED
setMethod(f="readData",
    signature="Genetic_Map",
    definition=function(object)
    {
        if (is.list(object@encodings[["chr"]])) {
        # chr information NOT IN the files
            object@gmapData = .read_gmap_no_internal_chrom_info(object)
        }
        else {
        # chr information IN the files
            if (object@encodings[["format"]]=="rds") {
                # object@filepaths should be a 1 length vector TODO test
                object@gmapData = .read_gmap_rds(object)
            } else {
                object@gmapData = list(data.frame())
            }
        }
        return(object)
    }
)

# now define the exposet object Genetic_Map methods
#' saveData perform the actual reading of the genetic maps data based on the encodings
#' @name saveData
#' 
#' @aliases saveData
#' @export
setGeneric(name="saveData",
    def=function(object, file)
    {
       standardGeneric("saveData")
    }
)

#' @rdname saveData
#' @return save a genetic map in RDS format
#' @aliases saveData
#' @examples
#' TOBEFIXED
setMethod(f="saveData",
    signature="Genetic_Map",
    definition=function(object, file)
    {
        # not necessary but coherent
        object@encodings$format="rds"
        # now save this instance of Genetic_Map
        saveRDS(object, file=file)

        return(invisible(object))
    }
)
