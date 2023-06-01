
# ggtree example. Full documentation for ggtree can be found https://yulab-smu.top/treedata-book/

library(tidyverse)
library(ggtree)
library(treeio)
library(aplot)

### Tree object ################################################################

# Load a tree. Can either make a random tree or load from a newick file

# rtree comes from ape package but is provided by treeio
random_tree <- rtree(5)
class(random_tree)

newick_file <- "example_tree.newick"
newick_tree <- read.newick(newick_file)

# Then draw with ggtree
ggtree(random_tree)

# The trees are 'phylo' objects from the ape package
class(random_tree)

# To do anything useful with data you can convert it to a 'tbl_tree'
tree_tb <- as_tibble(random_tree)
tree_tb
class(tree_tb)

# Can then add new data
tree_tb$country <- "UK"

# This can be manipulated and converted back to a phylo object.
as.phylo(tree_tb)

# More useful for plotting is to convert to 'treedata' which can keep the associated data
tree_data <- as.treedata(tree_tb)
class(tree_data)


# WARNING: converting with as_tibble twice means it won't convert out so easily
tree_tb2 <- as_tibble(tree_tb)
class(tree_tb2) # note that lack of "tbl_tree" in type

ggtree(as.phylo(tree_tb))
ggtree(as.phylo(tree_tb2)) # see how the branch lengths don't show right
as.treedata(tree_tb)
as.treedata(tree_tb2) # Also no longer includes branch lengths despite them being in tibble


# Overall I like to read newick, convert to tibble, add all metadata and processing,
# Then convert to treedata to draw tree.
# This avoids having to use the %<+% operator as well which felt clunky

### Draw Tree ##################################################################

# Full details https://yulab-smu.top/treedata-book/chapter4.html#basic-tree-visualization

# Basic method
random_tree <- rtree(5)
ggtree(random_tree)

# We can add more data after converting to tibble, and then draw this.
# The most basic additions are tiplab and tippoint, which can take various aesthetics.
# Theme_tree2 adds a scale to the bottom.
tree_tb <- as_tibble(random_tree) %>% 
  mutate(species = case_when(
    is.na(label) ~ 'ancestor',
    node > 2 ~ 'cat',
    T ~ 'dog'
  )) 

tree_data <- tree_tb %>% 
  as.treedata()

tree_data %>% 
  ggtree(aes(color = species), layout='roundrec') +
  geom_tiplab(aes(label = label), size=5, offset=0.1) +
  geom_tippoint(aes(color = species), size = 5) +
  geom_nodelab(aes(label = node), geom='label') +
  theme_tree2()


# lets make a nice function for later:
draw_treedata <- function(tr_data) {
  tr_data %>% 
    ggtree(aes(color = species), layout='roundrec') +
    geom_tiplab(aes(label = label), size=5, offset=0.1) +
    geom_tippoint(aes(color = species), size = 5) +
    geom_nodelab(aes(label = node), geom='label') +
    theme_tree2() %>% 
    return()
}

### Manipulate tree ############################################################

# The tree objects can be manipulated fairly easily.
# Be aware that functions act differently on each tree object type

# Functions are child, parent, offspring, ancestor, sibling, MRCA
# https://yulab-smu.top/treedata-book/chapter2.html#accesor-tidytree


# Note that tbl_tree will return a tbl_tree, whilst phylo/treedata return the nodes

# example to get the cat mrca
cat_nodes <- filter(tree_tb, species=='cat')$node
cat_mrca <- MRCA(tree_tb, cat_nodes)
cat_mrca


### Subset a tree ##############################################################

# treeio has the function tree_subset which takes a node and returns related nodes 
# up to a certain distance bake

sub_tree <- tree_subset(tree_data, 't5', levels_back=2)
draw_treedata(sub_tree)

# Alternatively you can give an internal node:
sub_tree <- tree_subset(tree_data, cat_mrca$node, levels_back=0)
draw_treedata(sub_tree)

# If you want more control you will instead need use drop.tip
tree_data %>% 
  drop.tip('t4') %>% 
  draw_treedata()

# Here's a function to get the subtree made from just the leaves you want to keep
induced_subtree <- function(tree, leaves) {
  # Expects tree to be in tibble format
  if (! is_tibble(tree)){
    print("tree should be tibble")
    return()
  }
  
  all_leaves <- tree %>% 
    select(label) %>% 
    drop_na()
  
  drop_leaves <- all_leaves %>% 
    filter(! label %in% leaves)
  
  tree <- as.treedata(tree) %>% 
    drop.tip(drop_leaves$label) %>% 
    as_tibble()
  
  return(tree)
}

tree_tb %>% 
  induced_subtree(c('t2', 't1', 't5')) %>% 
  as.treedata() %>% 
  draw_treedata()


### Scaling tree and making space ##############################################

# Often you'll have a tree with labels which get cutoff.
# This can be fixed with xlim if you know the exact bounds you want
# Easier is to use hexpand to increase the width by a proportion

rtree(8) %>% 
  as_tibble() %>% 
  mutate(label = sprintf("%s really long words", label)) %>% 
  as.treedata() %>% 
  ggtree() +
  geom_tiplab() +
  theme_tree2() +
  hexpand(0.3, direction=1)

# Many commands will act on the absolute size of the tree.
# As such it can be helpful, when branch length is less important,
# to scale the tree so overall length is 1.
# The following normalise command does this.

get_depth <- function(tree, n) {
  to_search = c(n)
  depths = c(0)
  max_d = 0
  
  while (length(to_search) > 0) {
    m <- tail(to_search, 1)
    to_search <- head(to_search, -1)
    current_depth <- tail(depths, 1)
    depths <- head(depths, -1)
    
    this_branch_length <- tree %>% 
      filter(node == m) %>% 
      .$branch.length %>% 
      sum(na.rm = T)
    
    current_depth <- current_depth + this_branch_length
    max_d <- max(max_d, current_depth)
    children <- child(tree, m)$node
    to_search <- append(to_search, children)
    
    if (length(children) > 0) {
      for (i in 1:length(children)) {
        depths <- append(depths, current_depth)
      }
    }
  }
  
  return(max_d)
}

# Scales tree so that total depth is 1 (or target_length)
normalise_tree <- function(tree_tb, target_length = 1) {
  root <- tree_tb %>% 
    filter(parent == node) %>% 
    .$node
  
  depth <- get_depth(tree_tb, root)
  
  tree_tb %>% 
    mutate(branch.length = branch.length * target_length / depth) %>% 
    return()
}

# which makes it easy to use xlim say
rtree(8) %>% 
  as_tibble() %>% 
  mutate(label = sprintf("%s really long words", label)) %>% 
  normalise_tree() %>% 
  as.treedata() %>% 
  ggtree() +
  geom_tiplab() +
  theme_tree2() +
  xlim(0, 1.2)

### Adding heatmaps to trees with aplot ###################################################

# There are a few ways to do this. 
# The simplest way is to make heatmaps with geom_tile and simply put them next to the tree
# This uses a package called aplot which does the aligning


# Make our random tree
tree_tb <- rtree(10) %>% 
  as_tibble()

gg_tr <- ggtree(as.treedata(tree_tb)) +
  geom_tiplab()

# Let's make some random data
meta <- tree_tb %>% 
  filter(!is.na(label)) %>% 
  mutate(A = runif(n=10), B = runif(n=10), C = runif(n=10))

# We now make our tiles. Note that the aesthetics should be in the ggplot part
tile <- ggplot(meta, aes(x='A value', y=label, fill=A)) +
  geom_tile() +
  geom_text(aes(label = A))
tile

# We then use aplot. Note how it aligns as both have 'label' controlling the y position
tile %>% 
  insert_left(gg_tr)


# WARNING: Seems you must start with another plot and insert the tree.
# The following fails
gg_tr %>% 
  insert_right(tile)

# WARNING: aesthetics for the tile plot must be in the main ggplot
# following fails to align properly
tile2 <- ggplot(meta) +
  geom_tile(aes(x='A value', y=label, fill=A))
tile2 %>% 
  insert_left(gg_tr)


# Can wrap this up to make a bit prettier

make_heatmap <- function(meta, var, name, color_scheme) {
  var <- ensym(var)
  
  p <- ggplot(meta, aes(x = name, y=label, fill=!!var)) +
    geom_tile() +
    xlab(name) +
    theme(axis.title.y = element_blank(),
          axis.text.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(angle = 45, vjust=0.5),
    ) +
    color_scheme
  
  return(p)
}

color_scheme_A <- scale_fill_viridis_c(option="H", name="Value A")
tile_A <- make_heatmap(meta, A, "Value A", color_scheme_A)

color_scheme_B <- scale_fill_viridis_c(option="D", name="Value B")
tile_B <- make_heatmap(meta, B, "Value B", color_scheme_B)

# Can also use width which is always relative to first plot
tile_A %>% 
  insert_left(gg_tr, width=2) %>% 
  insert_right(tile_B, width=0.5)


# Can also add in other kinds of plots
bars <- ggplot(meta, aes(x=label, y=C, fill=A)) +
  geom_col() +
  coord_flip() +
  color_scheme_A
bars

tile_A %>% 
  insert_left(gg_tr, width=2) %>% 
  insert_right(tile_B, width=0.5) %>% 
  insert_right(bars)



### Adding heatmaps with gheatmap ##############################################

# This is a ggtree method that also works with circular layouts.
# WARNING: gheatmap will need the data as a base data.frame where indexes match the leaf labels

# Make our random tree with random data
tree_tb <- rtree(10) %>% 
  as_tibble()

meta <- tree_tb %>% 
  filter(!is.na(label)) %>% 
  mutate(value_A = runif(n=10), value_B = runif(n=10), value_C = runif(n=10))

meta_df <- meta %>% 
  select(value_A, value_B) %>% 
  as.data.frame()
rownames(meta_df) <- meta$label


gg_tr <- ggtree(as.treedata(tree_tb), layout = 'circular') +
  geom_tiplab() 
gg_tr

# now use gheatmap

p <- gheatmap(gg_tr, meta_df, colnames_angle=-60, width=0.5, hjust=0.6, offset=0.5) +
  scale_fill_continuous()
p

# width and offset are all absolute distances.

# To get another heatmap (with different colour scaling etc) we need ggnewscale

library(ggnewscale)

meta_df2 <- meta %>% 
  select(value_C) %>% 
  as.data.frame()
rownames(meta_df2) <- meta$label

p2 <- p + new_scale_fill()
gheatmap(p2, meta_df2, colnames_angle=-60, width=0.25, hjust=0.6, offset=2.0)


# Dealing with all the offsets, width and new scales can get annoying.
# This function wraps up some of that. Works better with non-circular layout

add_heatmap <- function(gg_tr, vars, widths, gaps, color_schemes, meta, font_size=3) {
  vars <- syms(vars)
  meta_df <- meta %>% select(!!!vars) %>% 
    as.data.frame()
  rownames(meta_df) <- meta$label
  
  offset <- 0
  p <- gg_tr + theme_tree2()
  for (i in 1:length(vars)) {
    var <- sym(vars[[i]])
    width <- widths[[i]]
    gap <- gaps[[i]]
    color_scheme <- color_schemes[[i]]
    
    offset <- offset + gap
    p <- gheatmap(p, select(meta_df, !!var), offset = offset - width/2, width=width,
                  font.size = font_size, colnames_angle=-45, hjust=0, color = NA) +
      color_scheme
    p <- p + new_scale_fill()
    offset <- offset + width
  }
  
  return(p)
}

gg_tr <- ggtree(as.treedata(tree_tb)) +
  geom_tiplab() 

heatmap_vars <- c('value_A', 'value_B', 'value_C')
widths <- c(0.25, 0.25, 0.25)
gaps <- c(0.5, 0.25, 0.5)
schemes <- c(
  scale_fill_viridis_c(option="D", name="value A"),
  scale_fill_viridis_c(option="H", name="value B"),
  scale_fill_viridis_c(option="G", name="value C")
)

heatmap <- add_heatmap(gg_tr, heatmap_vars, widths, gaps, schemes, meta)
heatmap


### Adding other kinds of graphs 