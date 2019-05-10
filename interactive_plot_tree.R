interactive_plot_tree <- function (physeq, method = "sampledodge", nodelabf = NULL, color = NULL, 
          shape = NULL, size = NULL, min.abundance = Inf, label.tips = NULL, 
          text.size = NULL, sizebase = 5, base.spacing = 0.02, ladderize = FALSE, 
          plot.margin = 0.2, title = NULL, treetheme = NULL, justify = "jagged", tooltip = NULL) 
{
  fix_reserved_vars = function(aesvar) {
    aesvar <- gsub("^abundance[s]{0,}$", "Abundance", aesvar, 
                   ignore.case = TRUE)
    aesvar <- gsub("^OTU[s]{0,}$", "OTU", aesvar, ignore.case = TRUE)
    aesvar <- gsub("^taxa_name[s]{0,}$", "OTU", aesvar, ignore.case = TRUE)
    aesvar <- gsub("^sample[s]{0,}$", "Sample", aesvar, ignore.case = TRUE)
    return(aesvar)
  }
  if (!is.null(label.tips)) {
    label.tips <- fix_reserved_vars(label.tips)
  }
  if (!is.null(color)) {
    color <- fix_reserved_vars(color)
  }
  if (!is.null(shape)) {
    shape <- fix_reserved_vars(shape)
  }
  if (!is.null(size)) {
    size <- fix_reserved_vars(size)
  }
  if (is.null(phy_tree(physeq, FALSE))) {
    stop("There is no phylogenetic tree in the object you have provided.\n", 
         "Try phy_tree(physeq) to see for yourself.")
  }
  if (!inherits(physeq, "phyloseq")) {
    method <- "treeonly"
  }
  treeSegs <- tree_layout(phy_tree(physeq), ladderize = ladderize)
  edgeMap = aes(x = xleft, xend = xright, y = y, yend = y)
  vertMap = aes(x = x, xend = x, y = vmin, yend = vmax)
  p = ggplot(data = treeSegs$edgeDT) + geom_segment(edgeMap) + 
    geom_segment(vertMap, data = treeSegs$vertDT)
  if (is.null(text.size)) {
    text.size <- phyloseq:::manytextsize(ntaxa(physeq))
  }
  if (!is.null(label.tips) & method != "sampledodge") {
    labelDT = treeSegs$edgeDT[!is.na(OTU), ]
    if (!is.null(tax_table(object = physeq, errorIfNULL = FALSE))) {
      taxDT = data.table(tax_table(physeq), OTU = taxa_names(physeq), 
                         key = "OTU")
      labelDT = merge(x = labelDT, y = taxDT, by = "OTU")
    }
    if (justify == "jagged") {
      labelMap <- aes_string(x = "xright", y = "y", label = label.tips, 
                             color = color)
    }
    else {
      labelMap <- aes_string(x = "max(xright, na.rm=TRUE)", 
                             y = "y", label = label.tips, color = color)
    }
    p <- p + geom_text(labelMap, data = labelDT, size = I(text.size), 
                       hjust = -0.1, na.rm = TRUE)
  }
  if (is.null(nodelabf)) {
    nodelabf = phyloseq:::howtolabnodes(physeq)
  }
  p = nodelabf(p, treeSegs$edgeDT[!is.na(label), ])
  p = nodelabf(p, treeSegs$vertDT[!is.na(label), ])
  if (is.null(treetheme)) {
    treetheme <- theme(axis.ticks = element_blank(), axis.title.x = element_blank(), 
                       axis.text.x = element_blank(), axis.title.y = element_blank(), 
                       axis.text.y = element_blank(), panel.background = element_blank(), 
                       panel.grid.minor = element_blank(), panel.grid.major = element_blank())
  }
  if (inherits(treetheme, "theme")) {
    p <- p + treetheme
  }
  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }
  if (method != "sampledodge") {
    return(p)
  }
  dodgeDT = treeSegs$edgeDT[!is.na(OTU), ]
  dodgeDT = merge(x = dodgeDT, y = data.table(psmelt(physeq), 
                                              key = "OTU"), by = "OTU")
  if (justify == "jagged") {
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  if (!is.null(color) | !is.null(shape) | !is.null(size)) {
    setkeyv(dodgeDT, cols = c("OTU", color, shape, size))
  }
  else {
    setkey(dodgeDT, OTU, Sample)
  }
  dodgeDT[, `:=`(h.adj.index, 1:length(xright)), by = OTU]
  if (justify == "jagged") {
    dodgeDT[, `:=`(xdodge, (xright + h.adj.index * base.spacing * 
                              max(xright, na.rm = TRUE)))]
  }
  else {
    dodgeDT[, `:=`(xdodge, max(xright, na.rm = TRUE) + h.adj.index * 
                     base.spacing * max(xright, na.rm = TRUE))]
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  dodgeMap <- aes_string(x = "xdodge", y = "y", color = color, 
                         fill = color, shape = shape, size = size, names=tooltip)
  p <- p + geom_point(dodgeMap, data = dodgeDT, na.rm = TRUE)
  if (!is.null(size)) {
    p <- p + scale_size_continuous(trans = log_trans(sizebase))
  }
  if (any(dodgeDT$Abundance >= min.abundance[1])) {
    pointlabdf = dodgeDT[Abundance >= min.abundance[1], ]
    p <- p + geom_text(mapping = aes(xdodge, y, label = Abundance), 
                       data = pointlabdf, size = text.size, na.rm = TRUE)
  }
  if (!is.null(label.tips)) {
    tiplabDT = dodgeDT
    tiplabDT[, `:=`(xfartiplab, max(xdodge)), by = OTU]
    tiplabDT <- tiplabDT[h.adj.index == 1, .SD, by = OTU]
    if (!is.null(color)) {
      if (color %in% sample_variables(physeq, errorIfNULL = FALSE)) {
        color <- NULL
      }
    }
    labelMap <- NULL
    if (justify == "jagged") {
      labelMap <- aes_string(x = "xfartiplab", y = "y", 
                             label = label.tips, color = color)
    }
    else {
      labelMap <- aes_string(x = "max(xfartiplab, na.rm=TRUE)", 
                             y = "y", label = label.tips, color = color)
    }
    p <- p + geom_text(labelMap, tiplabDT, size = I(text.size), 
                       hjust = -0.1, na.rm = TRUE)
  }
  min.x <- -0.01
  max.x <- dodgeDT[, max(xright, na.rm = TRUE)]
  if ("xdodge" %in% names(dodgeDT)) {
    max.x <- dodgeDT[, max(xright, xdodge, na.rm = TRUE)]
  }
  if (plot.margin > 0) {
    max.x <- max.x * (1 + plot.margin)
  }
  p <- p + scale_x_continuous(limits = c(min.x, max.x))

  return(p)
}
