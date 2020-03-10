#_____________________________________________________________
# Part 2
# Blocs
#_____________________________________________________________


# read the augmented genetic map
# localy
augmented_genetic_map_dir = '~/R_workspace/data/outputs/augmented_genetic_map'
augmented_map_df = get_augmented_genetic_map(augmented_genetic_map_dir, chromosomes=1:23)
# mine
create_blocks(augmented_map_df, delta=1e-3, save_blocs=TRUE, output='~/R_workspace/data/outputs/blocs_yme')
# note mine
create_blocks(genMap_all, delta=1e-3, save_blocs=TRUE, output='~/R_workspace/data/outputs/blocs_cc')


# load the blocs

# mine
aa = get_blocs('~/R_workspace/data/outputs/blocs_yme')
# not mine
kk = get_blocs('~/R_workspace/data/outputs/blocs_cc')

# not mine original
load('~/R_workspace/gwhap_ressource/winDF_all.Rdata')
pp = winDF[winDF$min_cM == 1e-3]

get_blocs <- function(blocs_dir, chromosomes=1:23){
  blocs_df = c()
  for (chr in 1:22){
    blocs_chr = sprintf('%s/blocs_chr%s.txt', blocs_dir, chr)
    if(file.exists(blocs_chr)){blocs_df = rbind(blocs_df, read_delim(blocs_chr, delim='\t'))}
  }
  return(blocs_df)
}


#genome wide
sprintf('Start time lecture : %s', Sys.time())
augmented_genetic_map_dir = '/neurospin/tmp/ymekki/new_repo/gwhap/augmented_genetic_map'
augmented_map_df = get_augmented_genetic_map(augmented_genetic_map_dir, chromosomes=1:23)
sprintf('End time lecture : %s', Sys.time())

# Create blocs
sprintf('Start time creation blocs : %s', Sys.time())
output = '/neurospin/tmp/ymekki/new_repo/gwhap/blocs'
df_blocks = create_blocks(augmented_map_df, delta=1e-3, save_blocs=TRUE, output=output)
sprintf('End time creation blocs : %s', Sys.time())


