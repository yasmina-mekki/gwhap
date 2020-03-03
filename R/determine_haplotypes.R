#' haplotypes bloc matrix for bgen file
#' This function is designed for UKB haplotypes standard :
#' one .bgen file, .bgen.bgi variant index file and .sample file containing the participant's ID per chromosome.
#' @param chr string, chromosome code of the desired bloc.
#' @param start An integer specifying the starting position (bp) of the bloc.
#' @param end An integer specifying the ending position (bp) of the bloc.
#' @param bgnfile A path to a .bgen file for the desired chromosome.
#' @param max_entries_per_sample An integer specifying the maximum number of probabilities expected per variant per sample.
#' This is used to set the third dimension of the data matrix returned.
#'
#' @return haplotype count matrix for all haplotypes oserved in the bloc.
#' @import rbgen
#' @import DBI
#' @import dummies
#' @export
#'
determine_haplotypes=function(chr, start, end, bgnfile, samples_index, max_entries_per_sample=4)
{
  # we should use the ID not their index ...
  # for the returned object, we should keep the ID as well

  # get the IDs of the partcipants indexes
  # this part should be done outside of the function
  # use get_bgi_file and get_bgen_file functions
  # use the .sample to get both index and ID ?

  bgifile = sprintf("%s.bgi", bgnfile)
  con = (dbConnect(RSQLite::SQLite(), bgifile))
  allVar = data.frame(dbReadTable(con, "Variant"))
  dbDisconnect(con)
  dummybg = bgen.load(
    filename = bgnfile,
    data.frame(
      chromosome = '',
      start = allVar$position[1],
      end = allVar$position[1]
    ),
    max_entries_per_sample = max_entries_per_sample
  )

  # filter the partcipant ID using the their index
  mysamples = dummybg$samples[samples_index]

  # read the bgen file
  # use the function get_bgen_file ==> mybg = get_bgen_file(bgnfile, start+1, end-1, sample=mysamples, max_entries_per_sample=4)
  mybg = bgen.load(
    filename = bgnfile,
    data.frame(
      chromosome = '',
      start =start-1, # to make sure we have the starting position
      end = end+1 # to make sure we have the ending position
    ),
    samples = mysamples,
    max_entries_per_sample = max_entries_per_sample
  )

  # get genotype probability values, indexed by variants, samples, and by genotype.
  # c(4) and c(2) correspond to the genotype probability fields. they depend on max_entries_per_sample parameters chosen. Justify the value==4.
  allhap_2 = t(mybg$data[, , c(4)])
  allhap_1 = t(mybg$data[, , c(2)])
  remove(mybg)
  colnames(allhap_1) = colnames(allhap_2) = c()

  # concatenate and
  all_hap_1_2 = rbind(allhap_1, allhap_2)
  all_hapscodes = apply(all_hap_1_2, 1, function(x) {paste(unlist(x), collapse = '')})

  # factor
  haprle = rle(sort(all_hapscodes))
  rlefact = factor(haprle$values)
  hapfact = factor(all_hapscodes)

  # dummification
  df = dummy.data.frame(data.frame(levels(rlefact)[as.numeric(hapfact)]))
  my.effect.diplo = df[1:length(mysamples),] + df[(length(mysamples) +1):(2 * length(mysamples)),]
  colNames=colnames(my.effect.diplo)
  haps=unlist(lapply(strsplit(colNames,'_'),function(x) x[length(x)]))

  colnames(my.effect.diplo)=apply(data.frame(chr,start, end,haps),1,paste,collapse='_')
  row.names(my.effect.diplo)=samples_index
  return(my.effect.diplo)
}

