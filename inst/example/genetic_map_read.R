library(gwhap)

########################################################################
# This R scripts contains example to read genetic map
# object Genetic_Map

########################################################################
# ::TOY:: read 1000 genome interpolated

# give the list of files to consider
f1 = system.file("extdata", "chr1.interpolated_genetic_map.gz",
                          package="gwhap", mustWork=TRUE)
f2 = system.file("extdata", "chr2.interpolated_genetic_map.gz",
                          package="gwhap", mustWork=TRUE)
# build a NAMED list mapping the chrom information
chr = list(1, 2)
names(chr) = c(f1, f2)

# now build the parameter for Genetic_Map object. Inspect the file to
# assign the right column name for cM and position.
# Here the chrom info is not read from the file it is read from the
# NAMED list chr
filepaths = c(f1, f2)
encodings = list("cM"="cM", "position"="bp","chr"=chr)

# get an instance of Genetic_Map
genetic_map = Genetic_Map(filepaths=filepaths, encodings=encodings)
print(genetic_map)

# finally read the corresponding data and have a quick look at it
genetic_map = readData(genetic_map)
print(genetic_map)


########################################################################
# ::ACTUAL DATA:: read 1000 genome interpolated
#~ /neurospin//brainomics/bio_resources/hg19_maps/interpolated_maps/1000-genomes-genetic-maps-master/interpolated_from_hapmap/chr21.interpolated_genetic_map.gz
filepaths
encodings
genetic_map


########################################################################
# ::ACTUAL DATA:: read 1000 genome
#~ /neurospin//brainomics/bio_resources/hg19_maps/1000GP_Phase3_b37/genetic_map_chr21_combined_b37.txt
# filter on the desired columns : position COMBINED_rate(cM/Mb)
#                                 Genetic_Map(cM)

# give the list of files to consider (here file for chr21 and chr22
NS_ROOT = '/neurospin//brainomics/bio_resources/hg19_maps/1000GP_Phase3_b37'
f1 = file.path(NS_ROOT, 'genetic_map_chr21_combined_b37.txt')
f2 = file.path(NS_ROOT, 'genetic_map_chr22_combined_b37.txt')

# build a NAMED list mapping the chrom information
chr = list(21, 22)
names(chr) = c(f1, f2)

# now build the parameter for Genetic_Map object.
# the columns in the file are : 'position', 'COMBINED_rate(cM/Mb)', 'Genetic_Map(cM)'
filepaths = c(f1, f2)
encodings = list("cM"="Genetic_Map(cM)",
                 "position"="position",
                 "chr"=chr)

# get an instance of Genetic_Map
genetic_map = Genetic_Map(filepaths=filepaths, encodings=encodings)
print(genetic_map)

# finally read the corresponding data and have a quick look at it
genetic_map = readData(genetic_map)
print(genetic_map)



########################################################################
# ::ACTUAL DATA:: read Rutgers
# filter on the desired columns : Sex_averaged_start_map_position, 
#                                 Build37_start_physical_positio,
#                                 Marker_name 
#                                 chromosome code
NS_ROOT = '/neurospin//brainomics/bio_resources/rutgers_map_v3'
filelist=sapply(21:22, function(i) file.path(NS_ROOT,sprintf("RUMapv3_B137_chr%d.txt", i)))
chr=lapply(21:22, function(i) i)
names(chr) = filelist

filepaths=filelist
encodings = list("cM"="Sex_averaged_start_map_position",
                 "position"="Build37_start_physical_position",
                 "chr"=chr)

# get an instance of Genetic_Map
genetic_map = Genetic_Map(filepaths=filepaths, encodings=encodings)
print(genetic_map)

# finally read the corresponding data and have a quick look at it
genetic_map = readData(genetic_map)
print(genetic_map)
print(dim(genetic_map@gmapData))
