library(ggtreeExtra)
library(ggtree) 
library(ggplot2) 
library(ggnewscale) 
library(treeio) 
library(tidytree) 
library(dplyr) 
library(RColorBrewer) 
library(ape)

# command line args
args <- commandArgs(trailingOnly = TRUE)
tree_file <- args[1]
tsv_file <- args[2]
output_dir <- args[3]
rxn_name <- args[4]

# data
tree <- read.tree(tree_file)
labs <- c(tree$tip.label)
info <- read.csv(tsv_file)
new <- info[info$hit_id %in% labs,]

new <- new %>%
  select(hit_id, everything())
result <- new
row.names(result) <- new$hit_id

keeptips <- tree$tip.label[tree$tip.label %in% new$hit_id]
p0 <- keep.tip(tree, keeptips)

#trees
p <- ggtree(p0, layout='fan') %<+% result

p1 <- p +
  geom_tippoint(aes(color = Phylum)) +
  theme(legend.position = 'right') 

p2 <-p1 + new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=Genome_type),,
    color=NA,
    size=1,
    offset=0.1,
  ) +
  scale_fill_manual(name = "Genome type", values=c("#B4B9BF","#878D96"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4))


#getting colors ready
colorCount <- n_distinct(info$hit_descript) # number of levels

p3 <-p2 + new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=hit_descript),,
    color=NA,
    size=1,
    offset=0.14,
  ) +
  scale_fill_manual(name = "Prokka function", 
                    values = colorRampPalette(brewer.pal(colorCount, "Pastel1"))(colorCount),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4))


ggsave(paste(output_dir,rxn_name,'.png', sep = ""), p3, width = 12, height = 8)