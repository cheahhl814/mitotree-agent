#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(getopt))

spec <- matrix(c(
  'tree',     't', 1, "character", "Input Newick tree file (IQ-TREE species tree)",
  'outgroup', 'o', 1, "character", "Outgroup taxon name for re-rooting",
  'output',   'w', 1, "character", "Output Newick file (default: species_tree.nwk)"
), byrow = TRUE, ncol = 5)

opt <- getopt(spec)

if (is.null(opt$tree) || is.null(opt$outgroup)) {
  cat("Usage: reroot_tree.R -t tree.nwk -o outgroup_name [-w output.nwk]\n")
  q(status = 1)
}
if (is.null(opt$output)) opt$output <- "species_tree.nwk"

tree <- read.tree(opt$tree)

if (!opt$outgroup %in% tree$tip.label) {
  cat(sprintf(
    "ERROR: outgroup '%s' not found in tree. Available tip labels:\n  %s\n",
    opt$outgroup, paste(sort(tree$tip.label), collapse = "\n  ")
  ))
  q(status = 1)
}

tree <- root(tree, outgroup = opt$outgroup, resolve.root = TRUE)
write.tree(tree, file = opt$output)
cat(sprintf("Tree re-rooted on '%s' -> %s\n", opt$outgroup, opt$output))
