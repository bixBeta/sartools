myPCA = function(dds, target_) {
  vsd = varianceStabilizingTransformation(dds, blind = T)

  # calculate the variance for each gene
  rv <- rowVars(assay(vsd))

  # select the ntop genes by variance
  select <- order(rv, decreasing = TRUE)[seq_len(500)]

  pca = prcomp(t(assay(vsd)[select, ]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum(pca$sdev^2)
  pVar.df <- as.data.frame(percentVar)
  pVar.df$x = as.factor(paste0("PC", rownames(pVar.df)))

  pVar.df = pVar.df[, order(names(pVar.df))]
  pVar.df$percentVar = pVar.df$percentVar * 100
  pVar.df$percentVar = round(pVar.df$percentVar, digits = 2)
  group = target$group
  intgroup.df <- as.data.frame(colData(vsd)[, "group", drop = FALSE])
  d <- data.frame(pca$x, name = rownames(pca$x))
  d2 <- left_join(target_, d, by = c("label" = "name"))
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(ggplot2))

  pc1 = ggplot(d2, aes(x = PC1, y = PC2, color = group)) +
    geom_point(size = 2.5) +
    geom_label_repel(
      aes(label = label),
      box.padding = 0.35,
      point.padding = 0.5,
      segment.color = 'grey55',
      show.legend = F
    ) +
    xlab(paste0("PC1: ", pVar.df$percentVar[1], "%")) +
    ylab(paste0("PC2: ", pVar.df$percentVar[2], "%"))

  return(list(pca = pca, plot = pc1))
}

# dds results ====

getDDSresults = function(dds_, combn_) {
  results.l = list()

  for (i in 1:nrow(combn_)) {
    message(paste0("Computing DE for ", combn_$num[i], "_vs_", combn_$denom[i]))
    results.l[[i]] <- results(
      dds_,
      contrast = c("group", combn_$num[i], combn_$denom[i]),
      alpha = 0.05
    ) |>
      data.frame()

    # colnames(results.l[[i]]) <- paste0(
    #   combn_$num[i],
    #   "_vs_",
    #   combn_$denom[i],
    #   ".",
    #   colnames(results.l[[i]])
    # )

    names(results.l)[[i]] <- paste0(
      combn_$num[i],
      "_vs_",
      combn_$denom[i]
    )
  }
  results.l2 = do.call(cbind, results.l)
  return(results.l2)
}
