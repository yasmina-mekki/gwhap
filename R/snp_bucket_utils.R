# This class was created following the principled description from:
# https://www.cyclismo.org/tutorial/R/s4Classes.html#creating-an-s4-class
# Documentation and NAMESPACE : http://r-pkgs.had.co.nz/man.html#man-classes
# See examples : https://github.com/ropensci/datapack/blob/master/R/DataObject.R

# Create the base Snp_Bucket S4 class
#
# This is used to represent the physical position of list of SNPs
#' An S4 class to represent a Snp_Bucket.
#'
#' @name Snp_Bucket
#' @rdname Snp_Bucket
#' @slot filepaths A filepath character vector
#' @slot encodings A Named list c('snp', 'position', 'chr'). The value 
#' corresponding to 'snp' and 'position' are column names or column number
#' in the files of the filepaths list. If the value
#' corresponding to the key 'chr' is a map (named list) with key being the
#' filename in the filepaths list and value corresponding to chromosmome number.
#' 'chr' may also be a numeric or a character and it should be interpreted
#' as a column corresponding to chrom in the files of filepaths.
#' @slot bucketData A list of data.frame to containe the physical position data.
#' This is a named list with key like 'chr1', 'chr2', 'chrX'. Each 
#' data.frame has two columns : 'snp' and 'position'
#' @export Snp_Bucket
Snp_Bucket <- setClass(
    # Set the name for the class
    "Snp_Bucket",

    # Define the slots
    slots = c(
            filepaths = "character",
            encodings = "list",
            bucketData = "list"
            ),

    # Set the default values for the slots. (optional)
    prototype=list(
            filepaths = c(),
            encodings   = list("chr"=NULL, "snp"=NULL, "position"=NULL),
            bucketData = list(data.frame())
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
        if(!all(names(object@encodings) %in% c("chr", "snp", "position", "format"))) {
            print(names(object@encodings))
            return("The encodings list should contain chr, snp, position and format keys")
        }
        # specific encodings@format validity checks
        if (!(object@encodings[["format"]] %in% c("", "table", "bgen"))) {
            return("The encodings$format should be table or bgen")
        }
        # specific encodings@chr validity checks
        if (!(is.list(object@encodings[["chr"]]) | is.character(object@encodings[["chr"]]) | is.numeric(object@encodings[["chr"]]))) {
            return("The encodings@chr should be list or a charcacter/numeric")
        }
        if (is.list(object@encodings[["chr"]])) {
            if (!all(names(object@encodings[["chr"]]) %in% filepaths)) {
                return("The encodings@chr list does not match wirh filepaths vector")
            }
            if (!all(sapply(object@encodings[["chr"]], function(i) is.numeric(i)))) {
                return("The encodings@chr should contain only numeric as chromosome reference")
            }
        }
        if (!is.list(object@bucketData)) {
            return("The gmapData should be a list of data.frame")
        }
        if (!is.data.frame(object@bucketData[[1]])) {
            return("The gmapData should be a list of data.frame")
        }
        return(TRUE)
    }
)

########################################################################
# Define hidden ancillary functions
#
# plink/bim oriented method to read the Snp_Bucket
.read_snp_bucket_with_internal_chrom_info <- function(object) {
    
    outdf = data.frame()
    for (f in object@filepaths) {
        df = fread(f)
        outdf = rbind(outdf, df)
    }

    # test whether col reference are either all numeric or all character
    if (is.numeric(object@encodings$chr)) {
        if (!(is.numeric(object@encodings$position) & is.numeric(object@encodings$snp))) {
            print("encodings$chr and encodings$position/snp should be all the same type (col-index or col-names)")
            return(list(data.frame()))
        }
    }
    
    # select and name the columns based on encodings information
    col_tokeep = names(object@encodings)[names(object@encodings) %in% c('snp', 'position', 'chr')]
    col_toextract = as.vector(sapply(col_tokeep, function(s) object@encodings[[s]]))
    
    #outdf = subset(outdf, select=col_toextract)
    outdf = outdf[, ..col_toextract]
    colnames(outdf) = col_tokeep
    
    # create a list of dataframe with chr key
    outdf_list = list()
    for (i in unique(outdf$chr)) { # we know outdf has a chr column
        outdf_list[[sprintf('chr%d',i)]] = outdf[outdf$chr==i, c('snp', 'position')]
    }

    return(outdf_list)
}

# bgen oriented method to read the Snp_Bucket
.read_snp_bucket_no_internal_chrom_info <- function(object) {

    outdf_list = list()
    for (key in names(object@encodings$chr)) {
        chr_num = object@encodings$chr[[key]]
        if (file.exists(key)) {
            df = get_bgi_file(key)
        }
        # TODO check this hard coded stuff work all the time.
        colnames(df) <- c('chromosome', 'position', 'snp', 'number_of_alleles', 'allele1', 'allele2', 'file_start_position', 'size_in_bytes')
        df = df[, c('position', 'snp')]
        
        # append the dataframe in the existing list w/ named chrom
        outdf_list[[sprintf('chr%d',chr_num)]] = df
    }

    return(outdf_list)
}


# now define the exposet object Genetic_Map methods
#' readData perform the actual reading of the genetic maps data based on the encodings
#' @name readData
#' 
#' @aliases readData
#' @export
# TODO this should factorized elsewhere
#~ setGeneric(name="readData",
#~     def=function(object)
#~     {
#~        standardGeneric("readData")
#~     }
#~ )


#' @import data.table
#' @import DBI
#' @rdname readData
#' @return data for map based on mutiple chromfiles
#' @aliases readData
#' @examples
#' TOBEFIXED
setMethod(f="readData",
    signature="Snp_Bucket",
    definition=function(object)
    {
        if (is.list(object@encodings[["chr"]])) {
        # chr information not in the files
            object@bucketData = .read_snp_bucket_no_internal_chrom_info(object)
        }
        else {
        # chr information is in the files
            object@bucketData = .read_snp_bucket_with_internal_chrom_info(object)
        }

        return(object)
    }
)
