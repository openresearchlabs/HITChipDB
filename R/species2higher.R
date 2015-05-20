species2higher <- function (species.matrix, taxonomy, level, method) {

  # List all species for the given level (L0 / L1 / L2)")
  phylogroups <- levelmap(phylotypes = NULL, from = level, to = "species", taxonomy)

  summarized.matrix <- matrix(NA, nrow = length(phylogroups), ncol = ncol(species.matrix))
  rownames(summarized.matrix) <- sort(names(phylogroups))
  colnames(summarized.matrix) <- colnames(species.matrix)

  # Go through each phylogroup and summarize from species level
  for (pg in names(phylogroups)) {

    specs <- unique(phylogroups[[pg]])
    mat <- matrix(species.matrix[specs,], nrow = length(specs))

    if (method == "ave") { vec <- colMeans(mat) }
    if (method == "sum") { vec <- colSums(mat)  } 
    if (length(grep("rpa", method)) > 0) { vec <- colSums(mat) } # For RPA, use the sum for L1/L2

    summarized.matrix[pg, ] <- vec

  }

  summarized.matrix

}

