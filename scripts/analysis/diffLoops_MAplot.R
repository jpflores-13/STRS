## MA plot figure

## Perform differential loops analysis with DESeq2

library(DESeq2)
library(InteractionSet)
library(ggplot2)
library(dplyr)
library(glue)
library(stringr)
library(purrr)
library(apeglm)
library(RColorBrewer)
library(plotgardener)

# Load data ---------------------------------------------------------------
res <- readRDS("data/processed/hic/diffLoops/diffLoops_eGFP-YAP_noDroso_10kb.rds") |> 
  interactions()

# manipulate dataset for custom MA plot -----------------------------------
res_df <- res |>
  as.data.frame() |>
  dplyr::select(baseMean, log2FoldChange, padj) |> 
  mutate(isDE = case_when(
    log2FoldChange > 0 & padj < 0.1 ~ "TRUE - upreg",
    log2FoldChange < 0 & padj < 0.1 ~ "TRUE - downreg",
    padj > 0.1  ~ "FALSE",
    is.character("NA") ~ "FALSE"
  )) |> 
  arrange(isDE)

# custom MA plot ----------------------------------------------------------

(res_gg <- res_df |> 
   ggplot(aes(x = baseMean, y = log2FoldChange, color = isDE)) +
   geom_point(alpha = 1) +
   geom_hline(yintercept = 0,
              linetype = 2,
              color = "grey40") +
   scale_color_manual(values = c("TRUE - upreg" = "#F8766D",
                                 "TRUE - downreg" = "#619CFF",
                                 "FALSE" = "grey80")) +
   ylim(c(-4,4)) +
   scale_x_log10(breaks=c(1,2,5,10,20,50,100)) +
   labs(y = "Hi-C Contact log2FoldChange (sorbitol/untreated)",
        x = "mean of normalized counts") +
   theme_classic() +
   theme(
     legend.position = "NONE",
     # axis.text.x = element_text(size = 10),
     # axis.text.y = element_text(size = 10),
     # axis.title.x = element_text(size = 10),
     # axis.title.y = element_text(size = 10)
   ) +
   annotate(geom = "text",
            label = "Gained",
            x = 100,
            y = 3.5,
            color = "#F8766D") +
   annotate(geom = "text",
            label = paste0("n = ",
                           nrow(subset(res_df,
                                       isDE == "TRUE - upreg"))),
            x = 100,
            y = 3.1,
            color = "#F8766D") +
   annotate(geom = "text",
            label = "Lost",
            x = 100,
            y = -3.5,
            color = "#619CFF") +
   annotate(geom = "text",
            label = paste0("n = ",
                           nrow(subset(res_df,
                                       isDE == "TRUE - downreg"))),
            x = 100,
            y = -3.9,
            color = "#619CFF"))
# annotate(geom = "text",
#          label = "padj < 0.1",
#          x = 0.3,
#          y = 4)

ggsave("plots/diffLoops_MAplot.pdf",
       plot = res_gg,
       device = "pdf",
       width = 8,
       height = 5)

(densityMA <- ggplot(res_df, aes(y = log2FoldChange)) +
    geom_density(
      color = "#619CFF",
      fill = 4,
      alpha = 0.25) +
    geom_hline(yintercept = 0,
               linetype = 2,
               color = "grey40") +
    ylim(c(-4,4)) +
    xlim(c(0, 1.5)) +
    theme_classic() +
    theme(
          legend.position = "NONE",
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.line.x = element_blank())
  )

ggsave("plots/diffLoops_MAplot_density.pdf",
       plot = densityMA,
       device = "pdf",
       width = 3,
       height = 5)