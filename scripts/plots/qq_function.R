library(qqman)
library(dplyr)

# FUNCTION TO OBTAIN LAMBDA (GENOMIC CONTROL) AND QQPLOTS.
# IT USES ONE FILE (with header = TRUE): 
# file_to_read1 = data frame with p-values only
# IT NEEDS AS INPUT A FILE NAME FOR SAVING THE PLOT (out_name argument)
# The lambda genomic control value will be printed in the screen. The plot will be saved in the directory.

qq_plot <- function(df_input, out_name) {
    
  #read file1:
  colnames(df_input) <- "pvalue"

   #genomic control (lambda):
   chisq <- qchisq(1-df_input$pvalue,1)
   lambda_value <- median(chisq)/qchisq(0.5,1)
   lambda_value <- round(lambda_value, digits = 4)

   qq_filename <- paste("qqplot_", out_name, ".jpeg", sep="")

   #qqplot:
   jpeg(filename = qq_filename, width = 500, height = 500)
   qq(df_input$pvalue, cex.axis = 1.5)
   dev.off()
   
   cat("GC lambda value for ", out_name, "is: ", lambda_value, ".\n")
}
