#' haplotype block matrix
#' depends : rgben, RSQLite,  dummy
#' Assuming ukb_hap standard :
#' one .bgen file / chromosome
#' .bgen.bgi variant index file
#' index 1 of .sample file starts at line 3
#' @param chr : string, chromosome of the block ; ex. "chr1"
#' @param start : int,  starting position (bp) of the block
#' @param end : int, ending position (bp) of the block
#' @param bgnfile : .bgen file for the chromosome
#' @param samples_index : index of participants' ID in the .sample file
#' @return \item{my.effect.diplo}{ haplotype count matrix for all haplotypes oserved in the block }
#' @title  haplotypes block matrix for bgen file
#' @export read_haplotype_bgen
read_haplotype_block_bgen=function(chr,start,end,bgnfile,samples_index)
{

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
    max_entries_per_sample = 4
  )
  mysamples=dummybg$samples[samples_index]
  mybg = bgen.load(
    filename = bgnfile,
    data.frame(
      chromosome = '',
      start =start-1, # to make sure we have the starting position
      end = end+1 # to make sure we have the ending position
    ),
    samples = mysamples,
    max_entries_per_sample = 4
  )


  allhap_2 = t(mybg$data[, , c(4)])
  allhap_1 = t(mybg$data[, , c(2)])
  remove(mybg)


  colnames(allhap_1) = colnames(allhap_2) = c()

  all_hap_1_2 = rbind(allhap_1, allhap_2)
  all_hapscodes = apply(all_hap_1_2, 1, function(x) {
    paste(unlist(x), collapse = '')
  })
  haprle = rle(sort(all_hapscodes))
  rlefact = factor(haprle$values)
  hapfact = factor(all_hapscodes)
  df = dummy(data.frame(levels(rlefact)[as.numeric(hapfact)]), int = T)
  my.effect.diplo = df[1:length(mysamples),] + df[(length(mysamples) +
                                                     1):(2 * length(mysamples)),]
  colNames=colnames(my.effect.diplo)
  haps=unlist(lapply(strsplit(colNames,'_'),function(x) x[length(x)]))

  colnames(my.effect.diplo)=apply(data.frame(chr,start, end,haps),1,paste,collapse='_')
  row.names(my.effect.diplo)=samples_index
  return(my.effect.diplo)


}

