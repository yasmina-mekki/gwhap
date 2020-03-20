#' Determine haplotypes bloc for bgen file
#'
#' @description This function is designed for UKB haplotypes standard:
#' one .bgen file, .bgen.bgi variant index file and .sample file containing the participant's ID.
#'
#' @param chr string, chromosome code of the desired bloc.
#' @param start An integer specifying the starting position (bp) of the bloc.
#' @param end An integer specifying the ending position (bp) of the bloc.
#' @param bgen_file A path to a .bgen file for the desired chromosome.
#' @param samples_index index list corresponding to the participant ID. You can found it in the .sample file
#' @param mysamples correspondance of the partcipant index and ID using .bgi/.bgen file
#' @param max_entries_per_sample An integer specifying the maximum number of probabilities expected per variant per sample.
#' This is used to set the third dimension of the data matrix returned.
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return haplotype count matrix for all haplotypes oserved in the bloc. The haplotypes name is coding as follow: chr_start_end_haps-code
#' @import rbgen
#' @import DBI
#' @import dummies
#' @export
#'
determine_haplotypes_per_bloc = function(chr, start, end, bgen_file, samples_index, mysamples, max_entries_per_sample=4, verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

  # read the bgen file
  mybg = get_bgen_file(file_path = bgen_file,
                       start = start-1, # to make sure we have the starting position
                       end = end+1, # to make sure we have the ending position
                       samples = mysamples,
                       chromosome = '',
                       max_entries_per_sample = max_entries_per_sample)

  # get genotype probability values, indexed by variants, samples, and by genotype.
  # c(4) and c(2) correspond to the genotype probability fields. they depend on max_entries_per_sample parameters chosen. Justify the value==4.
  allhap_2 = t(mybg$data[, , c(4)])
  allhap_1 = t(mybg$data[, , c(2)])
  #remove(mybg)
  colnames(allhap_1) = colnames(allhap_2) = c()

  # concatenate
  all_hap_1_2 = rbind(allhap_1, allhap_2)
  all_hapscodes = apply(all_hap_1_2, 1, function(x) {paste(unlist(x), collapse = '')})

  # factor
  haprle = rle(sort(all_hapscodes))
  rlefact = factor(haprle$values)
  hapfact = factor(all_hapscodes)

  # dummification
  df = dummy.data.frame(data.frame(levels(rlefact)[as.numeric(hapfact)]))
  my.effect.diplo = df[1:length(mysamples), ] + df[(length(mysamples) +1):(2 * length(mysamples)), ]
  colNames = colnames(my.effect.diplo)
  haps = unlist(lapply(strsplit(colNames,'_'), function(x) x[length(x)]))

  # changing haps name format into chr_start_end_haps-code
  haps = vapply(strsplit(haps,"hapfact.."), `[`, 2, FUN.VALUE=character(1))

  colnames(my.effect.diplo) = apply(data.frame(chr, start, end, haps), 1, paste, collapse='_')
  row.names(my.effect.diplo) = samples_index
  return(my.effect.diplo)
}
