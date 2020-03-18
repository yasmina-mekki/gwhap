#_____________________________________________________________
# Part 2
# Blocs
#_____________________________________________________________


# read the augmented genetic map

# rutgers
augmented_genetic_map_dir = '~/R_workspace/data/outputs/augmented_genetic_map'
# 1000 genome interpolated
augmented_genetic_map_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map_1000_genome_interpolated'
augmented_map_df = get_augmented_genetic_map(augmented_genetic_map_dir, chromosomes=1:23)

# create blocs
sprintf('Start time : %s', Sys.time())
output = '/neurospin/tmp/ymekki/new_repo/gwhap/bloc_using_augmented_genetic_map_1000_genome_interpolated'
create_blocks(augmented_map_df, delta=1e-3, save_blocs=TRUE, output=output)
sprintf('End time : %s', Sys.time())


# load the blocs

# mine
blocs_df = get_blocs('/neurospin/tmp/ymekki/new_repo/gwhap/bloc_using_augmented_genetic_map_1000_genome_interpolated')

# note mine
load('/neurospin/brainomics/2019_gwhap/data/winDF_all.Rdata')
# filter on the delta
winDF_1e_3 = winDF[winDF$min_cM==1e-3, ]

# merge
bloc_df_merged = merge(blocs_df, winDF_1e_3, by.x = 'from_bp', by.y = 'fromBp')


