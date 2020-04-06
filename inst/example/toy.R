
# toy bgen
toy_bgen_fn <- system.file("extdata", "haplotypes.bgen",package="gwhap", mustWork=TRUE)
phased_dl.bgen = phased_data_loader.bgen(toy_bgen_fn)
# ranges = data.frame(chromosome = "1", start = 1, end = 4)
# diploHaplo.bgen = getDiploHaplo(phased_dl.bgen, ranges=ranges)
# diploHaplo.bgen
samples_selected = c("sample_0","sample_1","sample_2","sample_3")
haplotypes.bgen = determine_haplotypes_per_bloc(phased_dl.bgen, chromosome="1", start=1, end=4,
                                                   sample_iid=samples_selected,
                                                   sample_bgen_iid_code=samples_selected)
haplotypes.bgen


# toy haps
toy_hap_fn <- system.file("extdata", "haplotypes.haps" ,package="gwhap", mustWork=TRUE)
phased_dl.haps = phased_data_loader.haps(toy_hap_fn)
# diploHaplo.haps = getDiploHaplo(phased_dl.haps)
# diploHaplo.haps
samples_selected = c("sample_0","sample_1","sample_2","sample_3")
haplotypes.haps = determine_haplotypes_per_bloc(phased_dl.haps, chromosome="1", start=1, end=4,
                                                   sample_iid=samples_selected,
                                                   sample_bgen_iid_code=samples_selected)
haplotypes.haps
