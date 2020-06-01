#' Retrieve DiploHaplo matrix from a phased_data_loader instance
#'
#' @description This function retrieve DiploHaplo matrix from a phased_data_loader instance
#'
#' @export
#'
getDiploHaplo <- function(obj, ranges = NULL, samples_selected=NULL, samples_rename=samples_selected)
  UseMethod("getDiploHaplo", obj)


#' @export
#'
getDiploHaplo.default <- function(obj, ranges = NULL, samples_selected=NULL, samples_rename=samples_selected) {
  stop("Area not defined on any object")
}

#' Retrieve DiploHaplo matrix from a phased_data_loader instance
#'
#' @description This function retrieve DiploHaplo matrix from a phased_data_loader instance (vcf implementation)
#'
#' @export
#'
getDiploHaplo.vcf <- function(obj, ranges = NULL, samples_selected=NULL, samples_rename=samples_selected) {
  stop("Area not defined on any object")
}

#' Retrieve DiploHaplo matrix from a phased_data_loader instance
#'
#' @description This function retrieve DiploHaplo matrix from a phased_data_loader instance (haps implementation)
#'
#'
#' @import stringr
#' @export
#'
getDiploHaplo.haps <- function(obj, ranges = NULL, samples_selected=NULL, samples_rename=samples_selected) {

  phased_data.haps = obj[["phased_data"]]
  tmp = colnames(phased_data.haps)
  tmp = unlist(lapply(tmp, function(s) if (!(s %in% c("chrom", "rsid", "pos", "allele_1", "allele_2"))) s ))

  # so called strand1
  tmp1 = unlist(lapply(tmp, function(s) if ("strand1" == str_sub(s, -7)) s))
  first_chrom_seq = t(phased_data.haps[tmp1])
  rownames(first_chrom_seq) = str_sub(rownames(first_chrom_seq), 1, -9)

  # so called strand2
  tmp2 = unlist(lapply(tmp, function(s) if ("strand2" == str_sub(s, -7)) s))
  second_chrom_seq = t(phased_data.haps[tmp2])
  rownames(second_chrom_seq) = str_sub(rownames(second_chrom_seq), 1, -9)

  # concatenate
  all_hap_1_2.haps = rbind(first_chrom_seq, second_chrom_seq)

  return(all_hap_1_2.haps)
}

#' Retrieve DiploHaplo matrix from a phased_data_loader instance
#'
#' @description This function retrieve DiploHaplo matrix from a phased_data_loader instance (bgen implementation)
#'
#' @export
#'
getDiploHaplo.bgen <- function(obj, ranges = NULL, samples_selected=NULL, samples_rename=samples_selected) {

  # get a the bgen file descriptor from the obj passed
  fn_bgen = obj[['full_fname_bgen']]
  max_entries = obj[['max_entries']]
  fn_bgi = obj[['full_fname_bgi']]
  chrom_name_degenerated = obj[['chrom_name_degenerated']]

  # test args
  if (dim(ranges)[1] !=1) {
      # TOFIX raise an error
      print("ERROR ranges should be a 1 row data frame")
  }

  if (chrom_name_degenerated) {
        rangesfix = data.frame(chromosome="",start=ranges$start, end=ranges$end)
  }
  else {
      rangesfix = ranges
  }

  #
  phased_data.bgen = bgen.load(fn_bgen,
                               ranges=rangesfix,
                               samples = samples_selected,
                               max_entries=max_entries,
                               index.filename = fn_bgi)

  # get genotype probability values, indexed by variants, samples, and by genotype.
  # c(4) and c(2) correspond to the genotype probability fields. they depend on max_entries_per_sample parameters chosen. Justify the value==4.
  allhap_2 = t(phased_data.bgen$data[, , c(4)])
  allhap_1 = t(phased_data.bgen$data[, , c(2)])
  #remove(mybg)
  colnames(allhap_1) = colnames(allhap_2) = c()

  # concatenate
  all_hap_1_2.bgen = rbind(allhap_1, allhap_2)

  return(all_hap_1_2.bgen)
}



#' Determine haplotypes bloc
#'
#' @description This function is designed for UKB haplotypes standard:
#' one .bgen file, .bgen.bgi variant index file and .sample file containing the participant's ID.
#'
#' @param chromosome string, chromosome code of the desired bloc.
#' @param start An integer specifying the starting position (bp) of the bloc.
#' @param end An integer specifying the ending position (bp) of the bloc.
#' @param phased_data_loader as obrained by phased_data_loader.bgen
#' @param sample_iid participant ID. It correspond to the ID_2 in the .sample file. NOT CLEAR. TOFIX
#' @param sample_bgen_iid_code a list of index corresponding to the partcipant ID using .bgi/.bgen file. Please notice that, in the bgen format the participant IDs are anonymized.
#' In order to get the correspondance, you need to look at the .sample file.
#' this file can be used as correspondance table between the subject ID and the 'anonimized'ID using the index.
#' Warning: the first row of the .sample file does'nt correspond to a participant ID and should be removed bedore the correspondance process.
#' @param max_entries_per_sample An integer specifying the maximum number of probabilities expected per variant per sample.
#' This is used to set the third dimension of the data matrix returned.
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return haplotype count matrix for all haplotypes oserved in the bloc. The haplotypes name is coding as follow: chr_start_end_haps-code
#'
#' @import rbgen
#' @import DBI
#' @import dummies
#' @export
#'
determine_haplotypes_per_bloc <- function(
                                     chromosome=character(0),
                                     start=integer(0),
                                     end=integer(0),
                                     phased_dl,
                                     sample_iid=NULL,
                                     sample_bgen_iid_code=NULL,
                                     verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

  ranges = data.frame(
    chromosome = chromosome,
    start = start,
    end = end
  )

  samples_selected = sample_bgen_iid_code
  samples_rename = sample_iid

  all_hap_1_2 = getDiploHaplo(phased_dl, ranges=ranges,
                              samples_selected=samples_selected,
                              samples_rename=samples_rename)
                              
  all_hapscodes = apply(all_hap_1_2, 1, function(x) {paste(unlist(x), collapse = '')})

  # factor
  haprle = rle(sort(all_hapscodes))
  rlefact = factor(haprle$values)
  hapfact = factor(all_hapscodes)

  # dummification
  library(dummies)
  df = dummy.data.frame(data.frame(levels(rlefact)[as.numeric(hapfact)]))
  my.effect.diplo = df[1:length(samples_selected), ] + df[(length(samples_selected) +1):(2 * length(samples_selected)), ]
  if (dim(my.effect.diplo)[2]==0) { # happens when all subjects homozyg for this haplo
      return(NULL)
  }
  colNames = colnames(my.effect.diplo)
  haps = unlist(lapply(strsplit(colNames,'_'), function(x) x[length(x)]))

  # changing haps name format into chr_start_end_haps-code
  haps = vapply(strsplit(haps,"hapfact.."), `[`, 2, FUN.VALUE=character(1))

  colnames(my.effect.diplo) = apply(merge(ranges, haps), 1, paste, collapse='_')
  # Fix colnames no type control on chr name but output must be correct
  # chr??_XXX_YYYY_!!!   with symbol to code for
  colnames(my.effect.diplo) = paste("chr", 
            colnames(my.effect.diplo)[!startsWith(colnames(my.effect.diplo), "chr")], sep="")

  row.names(my.effect.diplo) = samples_rename

  return(my.effect.diplo)
}

#' Determine haplotypes per chromosome
#'
#' @description Determine haplotypes per chromosome
#'
#' @param phased_dl A phased_data_loader returned by phased_data_loader.bgen() or phased_data_loader.haps().
#' @param chromosome chromosome code
#' @param start a list of integer specifying the starting position (bp) of the bloc.
#' @param end a list of integer specifying the ending position (bp) of the bloc.
#' @param sample_iid a list of participant ID. It correspond to the ID_2 in the .sample file. NOT CLEAR TOFIX. 
#' @param sample_bgen_iid_code a list of index corresponding to the partcipant ID using .bgi/.bgen file. NOT CLEAR TOFIX. Please notice that, in the bgen format the participant IDs are anonymized.
#' In order to get the correspondance, you need to look at the .sample file.
#' this file can be used as correspondance table between the subject ID and the 'anonimized'ID using the index.
#' Warning: the first row of the .sample file does'nt correspond to a participant ID and should be removed bedore the correspondance process.
#' @param max_entries_per_sample An integer specifying the maximum number of probabilities expected per variant per sample.
#' This is used to set the third dimension of the data matrix returned.
#' @param nb_core number of cores one want to use. The number of core available in the machine -2 are used by default
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return Data frame structure representing the haplotypes determined for all blocs column binded
#' @import parallel
#' @export
#'
determine_haplotypes_per_chromosome <- function(phased_dl,
                                                chromosome,
                                                start,
                                                end,
                                                sample_iid,
                                                sample_bgen_iid_code,
                                                nb_core=detectCores()-1,
                                                verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

  # apply on all blocs of one chromosome
  haplotype_combined <- mcmapply(
                            FUN = determine_haplotypes_per_bloc,
                            chromosome = chromosome,
                            start = start,
                            end   = end,
                            MoreArgs=list(
                                phased_dl=phased_dl,
                                sample_iid = sample_iid,
                                sample_bgen_iid_code = sample_bgen_iid_code,
                                verbose=FALSE),
                            USE.NAMES=TRUE,
                            mc.cores = nb_core)

  # purge NULL items in the list
  haplotype_combined <- do.call(cbind,
               haplotype_combined[!sapply(haplotype_combined, is.null)])

  return(haplotype_combined)
}


#' Determine haplotypes bloc
#'
#' @description This function is designed for UKB haplotypes standard:
#' one .bgen file, .bgen.bgi variant index file and .sample file containing the participant's ID.
#'
#' @param chromosome string, chromosome code of the desired bloc.
#' @param start An integer specifying the starting position (bp) of the bloc.
#' @param end An integer specifying the ending position (bp) of the bloc.
#' @param bgen_file A path to a .bgen file for the desired chromosome.
#' @param sample_iid participant ID. It correspond to the ID_2 in the .sample file
#' @param sample_bgen_iid_code a list of index corresponding to the partcipant ID using .bgi/.bgen file. Please notice that, in the bgen format the participant IDs are anonymized.
#' In order to get the correspondance, you need to look at the .sample file.
#' this file can be used as correspondance table between the subject ID and the 'anonimized'ID using the index.
#' Warning: the first row of the .sample file does'nt correspond to a participant ID and should be removed bedore the correspondance process.
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
legacy_determine_haplotypes_per_bloc <- function(chromosome, start, end, bgen_file, sample_iid, sample_bgen_iid_code, max_entries_per_sample=4, verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

  # read the bgen file
  mybg = get_bgen_file(file_path = bgen_file,
                       start = start-1, # to make sure we have the starting position
                       end = end+1, # to make sure we have the ending position
                       samples = sample_bgen_iid_code,
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
  my.effect.diplo = df[1:length(sample_bgen_iid_code), ] + df[(length(sample_bgen_iid_code) +1):(2 * length(sample_bgen_iid_code)), ]
  colNames = colnames(my.effect.diplo)
  haps = unlist(lapply(strsplit(colNames,'_'), function(x) x[length(x)]))

  # changing haps name format into chr_start_end_haps-code
  haps = vapply(strsplit(haps,"hapfact.."), `[`, 2, FUN.VALUE=character(1))

  colnames(my.effect.diplo) = apply(data.frame(chromosome, start, end, haps), 1, paste, collapse='_')
  row.names(my.effect.diplo) = sample_iid
  return(my.effect.diplo)
}


#' Determine haplotypes per chromosome
#'
#' @description Determine haplotypes per chromosome
#'
#' @param chromosome chromosome code
#' @param start a list of integer specifying the starting position (bp) of the bloc.
#' @param end a list of integer specifying the ending position (bp) of the bloc.
#' @param bgen_file A path to a .bgen file for the desired chromosome.
#' @param sample_iid a list of participant ID. It correspond to the ID_2 in the .sample file
#' @param sample_bgen_iid_code a list of index corresponding to the partcipant ID using .bgi/.bgen file. Please notice that, in the bgen format the participant IDs are anonymized.
#' In order to get the correspondance, you need to look at the .sample file.
#' this file can be used as correspondance table between the subject ID and the 'anonimized'ID using the index.
#' Warning: the first row of the .sample file does'nt correspond to a participant ID and should be removed bedore the correspondance process.
#' @param max_entries_per_sample An integer specifying the maximum number of probabilities expected per variant per sample.
#' This is used to set the third dimension of the data matrix returned.
#' @param nb_core number of cores one want to use. The number of core available in the machine -2 are used by default
#' @param verbose silent warning messages. FALSE by default.
#'
#' @return Data frame structure representing the haplotypes determined for all blocs column binded
#' @import parallel
#' @export
#'
legacy_determine_haplotypes_per_chromosome <- function(chromosome, start, end, bgen_file, sample_iid, sample_bgen_iid_code, max_entries_per_sample=4, nb_core=detectCores()-2, verbose=FALSE){

  # silent warning messages
  if(verbose == TRUE){options(warn=0)} else{options(warn=-1)}

  # apply on all blocs of one chromosome
  haplotype_combined <- do.call(cbind, mcmapply(FUN = legacy_determine_haplotypes_per_bloc,
                                                chromosome = chromosome,
                                                start = start,
                                                end   = end,
                                                bgen_file  = bgen_file,
                                                sample_iid = list(sample_iid),
                                                sample_bgen_iid_code = list(sample_bgen_iid_code),
                                                mc.cores = nb_core))
  return(haplotype_combined)
}
