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

# ape's root(..., resolve.root=TRUE) labels the new root node literally "Root"
# -- this is the only non-numeric internal-node label in an otherwise all-
# numeric-support tree, which breaks toytree's support-vs-name column
# detection downstream (verified live: it silently parsed every support value
# into the "name" column instead, with support left NaN). Blank it so the
# written tree has uniform (numeric-or-empty) internal labels.
if (!is.null(tree$node.label) && length(tree$node.label) > 0 && tree$node.label[1] == "Root") {
  tree$node.label[1] <- ""
}

write.tree(tree, file = opt$output)
cat(sprintf("Tree re-rooted on '%s' -> %s\n", opt$outgroup, opt$output))
