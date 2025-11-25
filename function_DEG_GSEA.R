library(fgsea)
library(ggplot2)
library(enrichplot)
library(patchwork)
library(clusterProfiler)
library(org.Rn.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(DESeq2)
library(pheatmap)
library(dplyr)
library(grid)
set.seed(1234)


# To run GSEA GO/KEGG/HALLMARK
# 1st: Rank the output object from DEG using function rank_DEG_for_GSEA
# 2nd: Input the ranked object to do GSEA analysis and will produce NES enrichment plot
# 3rd: Make GSEA plot
# 4th: Make heatmap of leading edge genes 


# Generate ranked gene list from DESeq2 output file
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
rank_DEG_for_GSEA <- function(DESEQObject) {
  res_df <- as.data.frame(DESEQObject)
  res_df <- res_df[!is.na(res_df$stat), ]
  res_df$gene <- rownames(res_df)

  res_df <- res_df %>%
    group_by(gene) %>%
    slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
    ungroup()

  gene_ranks <- res_df$stat
  names(gene_ranks) <- res_df$gene
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)
  return(gene_ranks)
}


# Run GO GSEA
# Input ranked gene list
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
run_GSEA_GO <- function (gene_ranks, org = "rat", terms_to_remove = NULL, out_csv = "GO_GSEA_BP.csv", out_pdf = "GO_GSEA_BP_plot.pdf", plot_title = "GO BP GSEA",top_pos_n = 45, top_neg_n = 35) {

  OrgDb <- switch(org,
                  "rat"   = org.Rn.eg.db,
                  "mouse" = org.Mm.eg.db,
                  "human" = org.Hs.eg.db,
                  stop("Org must be 'rat','mouse','human'")
                  )

  gene_df <- bitr(names(gene_ranks),
                  fromType = "SYMBOL",
                  toType = "ENTREZID",
                  OrgDb = OrgDb)

  gene_ranks_entrez <- gene_ranks[gene_df$SYMBOL]
  names(gene_ranks_entrez) <- gene_df$ENTREZID
  gene_ranks_entrez <- sort(gene_ranks_entrez, decreasing = TRUE)

  gse_go <- gseGO( geneList = gene_ranks_entrez,  OrgDb = OrgDb,     ont = "BP",    minGSSize = 20,    maxGSSize = 500,   eps = 0,   pvalueCutoff = 0.05,     verbose = FALSE )

  go_df <- gse_go@result

  top_pos <- go_df %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(top_pos_n)
  top_neg <- go_df %>% filter(NES < 0) %>% arrange(NES) %>% head(top_neg_n)
  top_ids <- c(top_pos$ID, top_neg$ID)

  if (!is.null(terms_to_remove))
    top_ids <- setdiff(top_ids, terms_to_remove)

  go_final <- go_df %>% filter(ID %in% top_ids)
  write.csv(go_final, out_csv, row.names = FALSE, quote = FALSE)

  if (nrow(go_final) == 0) {
    message("No GO terms left after filtering")
    return(go_final)
  }

  plot_df <- go_final %>% mutate(
    GeneRatio = sapply(strsplit(core_enrichment, "/"), length) / setSize,
    logP = -log10(p.adjust)
  ) %>% arrange(NES)

  plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description)

  p <- ggplot(plot_df, aes(NES, Description)) +
    geom_point(aes(size = GeneRatio, color = logP)) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal(base_size = 12) +
    labs(title = plot_title) +
    theme( panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA), axis.text.y = element_text(size = 10),  plot.title = element_text(size = 14, face = "bold") )
  print(p)
  pdf(out_pdf, width = 9, height = 8)
  print(p)
  dev.off()

  return(gse_go)
}


# GSEA plot for a single term
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
plot_gsea_go <- function(gsea_obj, go_term, save_dir = getwd(), show_plot = TRUE) {

  go_row <- gsea_obj@result[gsea_obj@result$Description == go_term, ]
  if (nrow(go_row) == 0)
    stop("GO term not found")

  nes <- round(go_row$NES, 3)
  padj <- signif(go_row$p.adjust, 3)
  plot_title <- paste0(go_term, " (NES=", nes, ", padj=", padj, ")")

  safe_name <- gsub(" ", "_", go_term)
  file_path <- file.path(save_dir, paste0("FigGO_gsea_", safe_name, ".pdf"))

  p <- gseaplot2(gsea_obj, geneSetID = go_row$ID, title = plot_title)

  if (show_plot) print(p)

  pdf(file_path, width = 7, height = 5)
  print(p)
  dev.off()

  message("Saved:", file_path)
}


# Leading-edge heatmap (fixed)
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
plot_leading_edge_heatmap <- function(gseGO_obj, go_term, dds,  top_fraction = 0.9, save_pdf = TRUE,  orgdb = org.Rn.eg.db) {

  go_res_df <- as.data.frame(gseGO_obj@result)
  term_row <- go_res_df[grepl(go_term, go_res_df$Description, TRUE), ][1, ]

  if (is.null(term_row$core_enrichment))
    stop("GO term not found.")

  leading_entrez <- strsplit(term_row$core_enrichment, "/")[[1]]
  leading_symbols <- mapIds(orgdb, keys = leading_entrez, 
                            column = "SYMBOL", keytype = "ENTREZID",
                            multiVals = "first") %>% na.omit()

  out_name <- paste0("genesLEdgeGO_", gsub(" ", "_", go_term), ".csv")
  write.csv(data.frame(GeneSymbol = leading_symbols), out_name, row.names = FALSE)

  norm_counts <- counts(dds, normalized = TRUE)
  common_genes <- intersect(rownames(norm_counts), leading_symbols)
  mat <- norm_counts[common_genes, ]

  if (nrow(mat) < 2)
    stop("Not enough genes for heatmap")

  mat_scaled <- t(scale(t(mat)))

  n_top <- ceiling(nrow(mat_scaled) * top_fraction)
  top_genes <- names(sort(apply(mat_scaled, 1, var), TRUE))[1:n_top]

  sample_anno <- data.frame(Group = colData(dds)$condition)
  rownames(sample_anno) <- colnames(mat_scaled)

  # PRINT heatmap on screen
  pheatmap(mat_scaled[top_genes, ], annotation_col = sample_anno, fontsize_row = 6, main = paste0("Leading-edge: ", term_row$Description))

  if (save_pdf) {
    p_out <- pheatmap(mat_scaled[top_genes, ],
                      annotation_col = sample_anno,
                      fontsize_row = 6,
                      main = paste0("Leading-edge: ", term_row$Description),
                      silent = TRUE)

    pdf(paste0("FigLEdge_", gsub(" ", "_", go_term), ".pdf"),
        width = 6, height = 6)
    grid.newpage()
    grid.draw(p_out$gtable)
    dev.off()
  }
}


# Run KEGG analysis
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

run_kegg_gsea <- function(gene_ranks, orgdb = org.Rn.eg.db, organism_code = "rno", top_n = 20, remove_terms = c(), pdf_file = "FigKEGG_GSEA_plot.pdf") 
{ 
  message("Converting SYMBOL → ENTREZ IDs...")
  gene_df <- bitr(names(gene_ranks), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
  gene_ranks_entrez <- gene_ranks[gene_df$SYMBOL]
  names(gene_ranks_entrez) <- gene_df$ENTREZID
  gene_ranks_entrez <- sort(gene_ranks_entrez, decreasing = TRUE)

  message("Running GSEA (KEGG)...")
  gse_kegg <- gseKEGG(  geneList      = gene_ranks_entrez,    organism      = organism_code,    minGSSize     = 15,    maxGSSize     = 500,    pvalueCutoff  = 0.05,    verbose       = FALSE   )
  kegg_res <- as.data.frame(gse_kegg@result)
  write.csv(kegg_res, "KEGG_GSEA_results.csv", row.names = FALSE)

  message("Filtering top positive/negative NES...")
  top_pos <- kegg_res %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(top_n)
  top_neg <- kegg_res %>% filter(NES < 0) %>% arrange(NES) %>% head(top_n)
  top_ids <- c(top_pos$ID, top_neg$ID)

  message("Removing unwanted KEGG pathways...")
  final_ids <- setdiff(top_ids, remove_terms)
  plot_df <- kegg_res %>%  filter(ID %in% final_ids) %>% mutate( GeneRatio = sapply(strsplit(core_enrichment, "/"), length) / setSize, logP = -log10(p.adjust) ) %>% arrange(NES)

  plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description)

  # ---- BUILD PLOT ----
  p <- ggplot(plot_df, aes(x = NES, y = Description)) +
    geom_point(aes(size = GeneRatio, color = logP)) +
    scale_color_gradient(low = "blue", high = "red") +
    scale_size_continuous(range = c(3, 10)) +
    theme_minimal(base_size = 12) +
    labs( title = "KEGG GSEA", x = "Normalized Enrichment Score (NES)", y = "KEGG Pathway", color = "-log10(p.adjust)", size = "GeneRatio" ) +
    theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text.y = element_text(size = 9), plot.title = element_text(size = 14, face = "bold") )

  # PRINT plot on screen and save
  print(p)
  message("Saving PDF: ", pdf_file)
  pdf(pdf_file, width = 6, height = 8)
  print(p)      # print again INTO the PDF
  dev.off()
  return(list(gsea_object = gse_kegg, results_df = kegg_res, plot_df = plot_df,  plot = p
  ))
}



# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------
plot_kegg_gsea <- function(gsea_result, pathway_term, outdir = ".", save_only = FALSE) {

  # Allow user to pass list returned by run_kegg_gsea()
  if (is.list(gsea_result) && "gsea_object" %in% names(gsea_result)) {
    gsea_result <- gsea_result$gsea_object
  }

  # Ensure output directory exists (current dir usually always exists)
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # Find matching pathway
  match_row <- gsea_result@result[
    grep(pathway_term, gsea_result@result$Description, ignore.case = TRUE),
  ]

  if (nrow(match_row) == 0) {
    message("❌ No matching KEGG pathway found for: ", pathway_term)
    return(NULL)
  }

  # Extract pathway info
  pathway_id <- match_row$ID[1]
  pathway_desc <- match_row$Description[1]
  NES <- round(match_row$NES[1], 2)
  pval <- signif(match_row$p.adjust[1], 3)

  # Generate gseaplot2 combined figure
  gp <- gseaplot2(
    gsea_result,
    geneSetID = pathway_id,
    title = paste0(pathway_desc, "\nNES=", NES, ", adj.P=", pval)
  )

  # Combine components into a single plot
  p <- gp[[1]] + gp[[2]] + gp[[3]] +
    plot_layout(ncol = 1, heights = c(1, 1, 1))

  # Output filename (saved in current directory)
  clean_name <- gsub("[^A-Za-z0-9]", "_", pathway_desc)
  fig_path <- file.path(outdir, paste0("FigKEGG_gsea_", clean_name, ".pdf"))

  # Save PDF
  ggsave(fig_path, plot = p, width = 6, height = 5)

  # Print to screen unless suppressed
  if (!save_only) print(p)

  message("✅ Saved KEGG GSEA plot: ", fig_path)
  invisible(p)
}


plot_leading_edge_KEGGheatmap <- function(gsea_obj, term_name, dds, top_fraction = 0.9, pdf_file = NULL, orgdb = org.Rn.eg.db) {

  # If a list is passed (from run_kegg_gsea), extract the actual gsea object
  if (is.list(gsea_obj) && "gsea_object" %in% names(gsea_obj)) {
    gsea_obj <- gsea_obj$gsea_object
  }

  # Convert GSEA result to data frame
  res_df <- as.data.frame(gsea_obj@result)

  # Find the row for the requested term
  term_row <- res_df[grep(term_name, res_df$Description, ignore.case = TRUE), ][1, ]
  if (is.null(term_row$core_enrichment)) stop("Term not found in GSEA results.")

  # Leading edge genes (Entrez → SYMBOL)
  leading_entrez <- strsplit(term_row$core_enrichment, "/")[[1]]
  leading_symbols <- AnnotationDbi::mapIds(orgdb,
                                           keys = leading_entrez,
                                           column = "SYMBOL",
                                           keytype = "ENTREZID",
                                           multiVals = "first") %>% na.omit()

  # Save leading-edge gene list
  out_csv <- ifelse(is.null(pdf_file),
                    paste0("genesLEdge_", gsub(" ", "_", term_name), ".csv"),
                    sub("\\.pdf$", "_genes.csv", pdf_file))
  write.csv(data.frame(GeneSymbol = leading_symbols), out_csv, row.names = FALSE)
  message("📄 Saved leading-edge genes: ", out_csv)

  # Extract normalized counts from DESeq2
  norm_counts <- counts(dds, normalized = TRUE)
  common_genes <- intersect(rownames(norm_counts), leading_symbols)
  mat <- norm_counts[common_genes, ]

  if (nrow(mat) < 2) {
    message("Not enough leading-edge genes for heatmap.")
    return(invisible(leading_symbols))
  }

  # Scale rows
  mat_scaled <- t(scale(t(mat)))

  # Select top variable genes
  n_top <- ceiling(nrow(mat_scaled) * top_fraction)
  top_genes <- names(sort(apply(mat_scaled, 1, var), decreasing = TRUE))[1:n_top]

  # Sample annotation
  sample_anno <- data.frame(Group = colData(dds)$condition)
  rownames(sample_anno) <- colnames(mat_scaled)

  # Plot on screen
  pheatmap(mat_scaled[top_genes, ],
           annotation_col = sample_anno,
           fontsize_row = 6,
           cluster_cols = TRUE,
           main = paste0("Leading-edge Genes: ", term_row$Description))

  # Save to PDF if requested
  if (!is.null(pdf_file)) {
    p_out <- pheatmap(mat_scaled[top_genes, ],
                      annotation_col = sample_anno,
                      fontsize_row = 6,
                      cluster_cols = TRUE,
                      main = paste0("Leading-edge Genes: ", term_row$Description),
                      silent = TRUE)
    pdf(pdf_file, width = 6, height = 6)
    grid::grid.newpage()
    grid::grid.draw(p_out$gtable)
    dev.off()
    message("✅ Saved heatmap PDF: ", pdf_file)
  }

  return(invisible(leading_symbols))
}




# Run HALLMARK GSEA
# Input is DEG output gene ranks
# Organism: human mouse rat
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------
run_GSEA_HALLMARK <- function ( gene_ranks,    species = "rat", out_csv = "HALLMARK_GSEA_results.csv", out_pdf = "FigHALLMARK_GSEA_plot.pdf",  plot_title = "HALLMARK GSEA",  top_pos_n = 20,    top_neg_n = 20 ) {

  # Map species name for msigdbr
  species_name <- switch(species,
                         "rat"   = "Rattus norvegicus",
                         "mouse" = "Mus musculus",
                         "human" = "Homo sapiens",
                         stop("Species must be 'rat','mouse','human'")
  )

  # Convert SYMBOL → ENTREZID
  gene_df <- bitr(names(gene_ranks),
                  fromType = "SYMBOL",
                  toType   = "ENTREZID",
                  OrgDb    = switch(species,
                                    "rat"   = org.Rn.eg.db,
                                    "mouse" = org.Mm.eg.db,
                                    "human" = org.Hs.eg.db))

  gene_ranks_entrez <- gene_ranks[gene_df$SYMBOL]
  names(gene_ranks_entrez) <- gene_df$ENTREZID
  gene_ranks_entrez <- sort(gene_ranks_entrez, decreasing = TRUE)

  # Load HALLMARK gene sets
  hallmark_sets <- msigdbr(species = species_name, category = "H") %>%
    dplyr::select(gs_name, entrez_gene)

  # Convert to list
  pathways <- split(hallmark_sets$entrez_gene, hallmark_sets$gs_name)

  # Run GSEA
  gsea_hallmark <- GSEA( geneList = gene_ranks_entrez, TERM2GENE = hallmark_sets, minGSSize = 15, maxGSSize = 500, pvalueCutoff = 0.05, verbose = FALSE )
  hall_df <- as.data.frame(gsea_hallmark@result)
  top_pos <- hall_df %>% filter(NES > 0) %>% arrange(desc(NES)) %>% head(top_pos_n)
  top_neg <- hall_df %>% filter(NES < 0) %>% arrange(NES) %>% head(top_neg_n)
  top_ids <- c(top_pos$ID, top_neg$ID)

  hall_final <- hall_df %>% filter(ID %in% top_ids)
  # Save CSV
  write.csv(hall_final, out_csv, row.names = FALSE, quote = FALSE)

  if (nrow(hall_final) == 0) {
    message("No HALLMARK terms left after filtering")
    return(gsea_hallmark)
  }

  # Prepare plot table
  plot_df <- hall_final %>%
    mutate(GeneRatio = sapply(strsplit(core_enrichment, "/"), length) / setSize,
           logP = -log10(p.adjust)) %>% arrange(NES)
  plot_df$Description <- factor(plot_df$Description, levels = plot_df$Description)

  # Plot NES
  p <- ggplot(plot_df, aes(x = NES, y = Description)) +
    geom_point(aes(size = GeneRatio, color = logP)) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal(base_size = 12) +
    labs(title = plot_title, x = "Normalized Enrichment Score (NES)", y = "HALLMARK Pathway", color = "-log10(p.adjust)", size = "GeneRatio") +
    theme( panel.grid = element_blank(), panel.border = element_rect(color = "black", fill = NA), axis.text.y = element_text(size = 10),      plot.title = element_text(size = 14, face = "bold") )

  # Print to screen and save
  print(p)
  pdf(out_pdf, width = 9, height = 8)
  print(p)
  dev.off()

  message("✅ HALLMARK GSEA done. Plot saved to ", out_pdf)

  return(gsea_hallmark)
}


#gene_ranks <- rank_DEG_for_GSEA(res)
#gsea_hallmark <- run_GSEA_HALLMARK( gene_ranks,   species = "rat",   out_csv = "HALLMARK_GSEA_results.csv",  out_pdf = "FigHALLMARK_GSEA_plot.pdf",   plot_title = "HALLMARK GSEA",  top_pos_n = 20,   top_neg_n = 20)






plot_gsea_hallmark <- function(fgsea_res, pathways, gene_ranks, pathway_name, save_dir = "Figures", show_plot = TRUE) {
  # Ensure save directory exists
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # Find the pathway row
  path_row <- fgsea_res[fgsea_res$pathway == pathway_name, ]
  if (nrow(path_row) == 0) {
    message(paste("❌ Pathway not found:", pathway_name))
    return(NULL)
  }

  # Extract stats
  nes <- round(path_row$NES, 3)
  padj <- signif(path_row$padj, 3)
  title_text <- paste0(pathway_name, " (NES=", nes, ", padj=", padj, ")")

  # Create enrichment plot
  plt <- plotEnrichment(pathways[[pathway_name]], gene_ranks) +
    labs(title = title_text, x = "Ranked genes", y = "Enrichment score") +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black"),
      axis.ticks = element_line(color = "black")
    )

  # Print on screen if requested
  if (show_plot) print(plt)

  # Save to PDF
  safe_name <- gsub(" ", "_", pathway_name)
  out_file <- file.path(save_dir, paste0("FigHallmark_GSEA_", safe_name, ".pdf"))
  pdf(out_file, width = 7, height = 5)
  print(plt)
  dev.off()

  message(paste("✅ Saved GSEA plot:", out_file))
}





# ---------------------------------------------------------------------------
# Plot GSEA HALLMARK pathway
# ---------------------------------------------------------------------------
plot_gsea_hallmark <- function(gsea_obj, hallmark_term, pdf_file = NULL, show_plot = TRUE) {

  # Allow user to pass list returned by run_GSEA_HALLMARK()
  if (is.list(gsea_obj) && "gsea_object" %in% names(gsea_obj)) {
    gsea_obj <- gsea_obj$gsea_object
  }

  # Find the row matching the HALLMARK term
  hallmark_row <- gsea_obj@result[grep(hallmark_term, gsea_obj@result$Description, ignore.case = TRUE), ]

  if (nrow(hallmark_row) == 0) {
    message("❌ No matching HALLMARK pathway found for: ", hallmark_term)
    return(NULL)
  }

  # Extract pathway info
  pathway_id <- hallmark_row$ID[1]
  pathway_desc <- hallmark_row$Description[1]
  NES <- round(hallmark_row$NES[1], 3)
  padj <- signif(hallmark_row$p.adjust[1], 3)

  # Generate GSEA plot
  plt <- gseaplot2(
    gsea_obj,
    geneSetID = pathway_id,
    title = paste0(pathway_desc, "\nNES=", NES, ", adj.P=", padj)
  )

  # Patchwork combine if multiple panels
  if (length(plt) > 1) {
    library(patchwork)
    p <- plt[[1]] + plt[[2]] + plt[[3]] + plot_layout(ncol = 1, heights = c(1, 1, 1))
  } else {
    p <- plt
  }

  # Print to screen
  if (show_plot) print(p)

  # Save PDF if requested
  if (!is.null(pdf_file)) {
    ggsave(pdf_file, plot = p, width = 7, height = 5)
    message("✅ Saved HALLMARK GSEA plot to: ", pdf_file)
  }

  invisible(p)
}

# plot_gsea_hallmark(  gsea_obj = gsea_hallmark,   hallmark_term = "HALLMARK_G2M_CHECKPOINT",  pdf_file = "FigHALLMARK_G2M_CHECKPOINT.pdf" )


# ---------------------------------------------------------------------------
# Plot Leading-Edge Heatmap for HALLMARK Pathway
# ---------------------------------------------------------------------------
plot_hallmark_leading_edge_heatmap <- function(
    gsea_obj,
    hallmark_term,
    dds,
    top_fraction = 0.9,
    pdf_file = NULL,
    orgdb = org.Rn.eg.db
) {
  library(pheatmap)
  library(AnnotationDbi)
  
  # Accept list returned by run_GSEA_HALLMARK()
  if (is.list(gsea_obj) && "gsea_object" %in% names(gsea_obj)) {
    gsea_obj <- gsea_obj$gsea_object
  }
  
  # Find the HALLMARK row
  hallmark_row <- gsea_obj@result[grep(hallmark_term, gsea_obj@result$Description, ignore.case = TRUE), ]
  if (nrow(hallmark_row) == 0) {
    message("❌ No matching HALLMARK pathway found for: ", hallmark_term)
    return(NULL)
  }
  
  # Extract leading-edge genes (Entrez IDs)
  leading_entrez <- strsplit(hallmark_row$core_enrichment[1], "/")[[1]]
  
  # Map Entrez → SYMBOL
  leading_symbols <- AnnotationDbi::mapIds(
    orgdb,
    keys = leading_entrez,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  ) %>% na.omit()
  
  if (length(leading_symbols) < 2) {
    message("Not enough leading-edge genes for heatmap.")
    return(NULL)
  }
  
  # Save leading-edge genes
  gene_file <- paste0("genesLEdge_", gsub(" ", "_", hallmark_term), ".csv")
  write.csv(data.frame(GeneSymbol = leading_symbols), gene_file, row.names = FALSE)
  message("📄 Saved leading-edge genes to: ", gene_file)
  
  # Extract normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  common_genes <- intersect(rownames(norm_counts), leading_symbols)
  mat <- norm_counts[common_genes, ]
  
  # Scale rows
  mat_scaled <- t(scale(t(mat)))
  
  # Select top variable genes
  n_top <- ceiling(nrow(mat_scaled) * top_fraction)
  top_genes <- names(sort(apply(mat_scaled, 1, var), decreasing = TRUE))[1:n_top]
  
  # Sample annotation
  sample_anno <- data.frame(Group = colData(dds)$condition)
  rownames(sample_anno) <- colnames(mat_scaled)
  
  # Plot heatmap on screen
  pheatmap(
    mat_scaled[top_genes, ],
    annotation_col = sample_anno,
    fontsize_row = 6,
    cluster_cols = TRUE,
    main = paste0("Leading-edge Genes: ", hallmark_row$Description[1])
  )
  
  # Save PDF if requested
  if (!is.null(pdf_file)) {
    p_out <- pheatmap(
      mat_scaled[top_genes, ],
      annotation_col = sample_anno,
      fontsize_row = 6,
      cluster_cols = TRUE,
      main = paste0("Leading-edge Genes: ", hallmark_row$Description[1]),
      silent = TRUE
    )
    pdf(pdf_file, width = 6, height = 6)
    grid::grid.newpage()
    grid::grid.draw(p_out$gtable)
    dev.off()
    message("✅ Saved heatmap PDF: ", pdf_file)
  }
  
  invisible(leading_symbols)
}
# plot_hallmark_leading_edge_heatmap (  gsea_obj = gsea_hallmark,  hallmark_term = "HALLMARK_G2M_CHECKPOINT",  dds = dds,  top_fraction = 1,  pdf_file = "FigLEdge_HALLMARK_G2M_CHECKPOINT.pdf")






