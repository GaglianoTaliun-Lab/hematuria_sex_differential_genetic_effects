## This function was adapted from: https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
# The function highlights in red SNPs that surpass the indicated threshold and draws a horizontal line on the threshold.

# Load Packages --------------------------
library(readr)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Main -----------------------------------
gg.qqplot <- function(df, threshold = -log10(5e-8), title = ""){
  
  # read and format vector of p-values
  colnames(df) <- c("SNP", "P", "TEST")
  hlight <- filter(df, P >= threshold) %>% .[,1]
    
  # compute lambda (genomic control)
  df_subset <- filter(df, TEST == "ADD")
  chisq <- qchisq(1-10^(-(df_subset$P)),1)
  lambda_add <- round(median(chisq)/qchisq(0.5,1),3)
  
  df_subset <- filter(df, TEST == "ADD-SKAT")
  chisq <- qchisq(1-10^(-(df_subset$P)),1)
  lambda_skat <- round(median(chisq)/qchisq(0.5,1),3)
  
  df_subset <- filter(df, TEST == "ADD-SKATO")
  chisq <- qchisq(1-10^(-(df_subset$P)),1)
  lambda_skato <- round(median(chisq)/qchisq(0.5,1),3)
  
  # Add highlight and annotation information
  df.tmp <- df %>%
    mutate( is_highlight=ifelse(P > threshold, "yes", "no")) %>%
    mutate( is_annotate=ifelse(P > threshold, "yes", "no"))
  
  # Main plot
  qqplot <- ggplot(df.tmp, aes(sample = P, color = as.factor(TEST)), alpha=0.5, size=2) +
    # stat_qq() + stat_qq_line() +
    geom_qq(distribution = qexp, dparams=list(rate=log(10))) + geom_qq_line(distribution = qexp, dparams=list(rate=log(10))) +
    scale_color_manual(values = brewer.pal(3, "Set1"), labels = c(stringr::str_c("ADD (lambda = ",lambda_add,")"),
                                                                  stringr::str_c("ADD-SKAT (lambda = ",lambda_skat,")"),
                                                                  stringr::str_c("ADD-SKATO (lambda = ",lambda_skato,")"))) +
    
    # add plot title
    ggtitle(paste0(title)) +
    xlab("Expected -log10(p)") +
    ylab("Observed -log10(p)") +
    
    # Custom the theme:
    theme_bw(base_size = 22) +
    theme( 
      plot.title = element_text(hjust = 0.5)
    ) +
    guides(color=guide_legend(title="Test"))
  
  # Add label using ggrepel and ggplot_build to avoid overlapping
  # geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], aes(label=as.factor(SNP), alpha=0.7), size=5, force=1.3)
  # qqplot.subset <- ggplot(subset(df.tmp, is_annotate=="yes"), aes(sample = P)) + stat_qq() 
  # x.pnts <- ggplot_build(qqplot.subset)$data[[1]]$x
  # y.pnts <- ggplot_build(qqplot.subset)$data[[1]]$y
  
  qqplot # + geom_label_repel(data=df.tmp[df.tmp$is_annotate=="yes",], x = x.pnts, y = y.pnts + 1, aes(label=as.factor(SNP), alpha=0.7), size=5)
}