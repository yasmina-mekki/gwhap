# This class was created followingthe principled description from:
# https://www.cyclismo.org/tutorial/R/s4Classes.html#creating-an-s4-class

# Create the base Genetic_Map S4 class
#
# This is used to represent the Genetic_Map.
Genetic_Map <- setClass(
    # Set the name for the class
    "Genetic_Map",

    # Define the slots
    slots = c(
            filepaths = "character",
            encodings = "list",
            gmap = "data.frame"
            ),

    # Set the default values for the slots. (optional)
    prototype=list(
            filepaths = c(),
            encodings   = list("chr"=NULL, "cM"=NULL, "pos"=NULL),
            gmap = data.frame()
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
        if(!all(names(object@encodings) %in% c("chr", "cM", "pos"))) {
            return("The encodings list should contain chr, cM and pos keys")
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
    print(object@filepaths)
    return(data.frame(c("k")))
}

# now define the exposet object Genetic_Map methods
setGeneric(name="readData",
    def=function(object)
    {
       standardGeneric("readData")
    }
)

setMethod(f="readData",
    signature="Genetic_Map",
    definition=function(object)
    {
        if (is.list(object@encodings[["chr"]])) {
        # chr information not in the files
            object@gmap = .read_gmap_no_internal_chrom_info(object)
        }
        else {
        # chr information is in the files
            object@gmap = data.frame()
        }
        return(object)
    }
)

filepaths = c(system.file("extdata", "chr1.interpolated_genetic_map.gz",
                          package="gwhap", mustWork=TRUE),
              system.file("extdata", "chr2.interpolated_genetic_map.gz",
                          package="gwhap", mustWork=TRUE))
chr = list(1, 2)
names(chr) = filepaths
encodings = list("cM"="cM",
                 "pos"="bp",
                 "chr"=chr
                )
a = Genetic_Map(filepaths=filepaths, encodings=encodings)
print(a)
a = readData(a)
print(a)

#~ # example to test:
#~ sas internal
#~ /neurospin//brainomics/bio_resources/hg19_maps/1000GP_Phase3_b37/genetic_map_chr21_combined_b37.txt
#~ /neurospin//brainomics/bio_resources/hg19_maps/interpolated_maps/1000-genomes-genetic-maps-master/interpolated_from_hapmap/chr21.interpolated_genetic_map.gz
#~ /neurospin//brainomics/bio_resources/rutgers_map_v3
# avec internal carte bgen bgi ?
