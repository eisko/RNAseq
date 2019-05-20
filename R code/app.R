#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(devtools)
library(gplots)
library(ggplot2)
library(ggfortify)
library(calibrate)

# Import data
results_genes_sham <- read.csv("~/code/RnaSeq/ballgown/sham/sham_gene_results.tsv", sep = "\t", header = TRUE)
results_genes_aff <- read.csv("~/code/RnaSeq/ballgown/aff/aff_gene_results.tsv", sep = "\t", header = TRUE)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Interactive Volcano Plots"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(
     sidebarPanel(   
      textInput("gois",
                "Genes of Interest in comma seperated list please:",
                value = "Mmaa, Fgf21")
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
          plotOutput("ShamVolcano", hover = "sham_hover", width = "100%"),
          plotOutput("AffVolcano", hover = "aff_hover", width = "100%")
          # verbatimTextOutput("sham_label")
          # verbatimTextOutput("aff_label")
          
         
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
   output$ShamVolcano <- renderPlot({
     
     i <- -1
     file <- "sham"
     q <- 0.1
     fch <- 1
     
     # generate list from input
     goi <- as.vector(unlist(strsplit(input$gois, split = ", ")))
     
     # make volcano plot
     with(results_genes_sham, plot(i*log2(fc), -log10(pval), pch=20, main= "Affecteds vs. Unaffecteds, No Treatment", xlim=c(-10,10), ylim=c(0,11))) 
     
     # Add colored points: orange if q value<q, green if log2FC>1.5, red if log2FC<-1.5, blue if both)
     with(subset(results_genes_sham, qval<q ), points(i*log2(fc), -log10(pval), pch=20, col="orange"))
     with(subset(results_genes_sham, i*log2(fc)> fch), points(i*log2(fc), -log10(pval), pch=20, col="green"))
     with(subset(results_genes_sham, i*log2(fc)< -fch), points(i*log2(fc), -log10(pval), pch=20, col="red"))
     with(subset(results_genes_sham, qval<q & abs(log2(fc))>fch), points(i*log2(fc), -log10(pval), pch=20, col="blue"))
     #
     # plot genes of interest (goi) with one command
     goi_val <- filter(results_genes_sham, gene_name %in% goi)
     with(goi_val, points(i*log2(fc), -log10(pval), pch=20, col="purple"))
     with(goi_val, textxy(i*log2(fc), -log10(pval), labs=gene_name,col="purple",cex=1, font = 2))
     # 
     # Add text saying how many up, down, and total DE genes were found
     text(7, y=9, labels = "Up regulated in Affecteds", cex = 0.9)
     text(-7, y=9, labels = "Down regulated in Affecteds", cex = 0.9)
     # text(0, y=8, labels = paste("Total DE genes: ", nrow(results_genes_sham)), cex = 0.9)
   })
   
   output$AffVolcano <- renderPlot({
     
     i <- 1
     file <- "aff"
     q <- 0.1
     fch <- 1
     
     # generate list from input
     goi <- as.vector(unlist(strsplit(input$gois, split = ", ")))
    
     # make volcano plot
     with(results_genes_aff, plot(i*log2(fc), -log10(pval), pch=20, main= "Vancomycin vs. No Treatment, Affecteds", xlim=c(-10,10), ylim=c(0,11))) 
     
     # Add colored points: orange if q value<q, green if log2FC>1.5, red if log2FC<-1.5, blue if both)
     with(subset(results_genes_aff, qval<q ), points(i*log2(fc), -log10(pval), pch=20, col="orange"))
     with(subset(results_genes_aff, i*log2(fc)> fch), points(i*log2(fc), -log10(pval), pch=20, col="green"))
     with(subset(results_genes_aff, i*log2(fc)< -fch), points(i*log2(fc), -log10(pval), pch=20, col="red"))
     with(subset(results_genes_aff, qval<q & abs(log2(fc))>fch), points(i*log2(fc), -log10(pval), pch=20, col="blue"))
     #
     # plot genes of interest (goi) with one command
     goi_val <- filter(results_genes_aff, gene_name %in% goi)
     with(goi_val, points(i*log2(fc), -log10(pval), pch=20, col="purple"))
     with(goi_val, textxy(i*log2(fc), -log10(pval), labs=gene_name,col="purple",cex=1, font = 2))
     # 
     # Add text saying how many up, down, and total DE genes were found
     text(7, y=9, labels = "Up regulated with Vancomycin", cex = 0.9)
     text(-7, y=9, labels = "Down regulated with Vancomycin", cex = 0.9)
     # text(0, y=8, labels = paste("Total DE genes: ", nrow(results_genes_aff)), cex = 0.9)
   })
   
   # output$sham_label <- renderPrint({
   #   # Import data
   #   x <- input$sham_hover$x
   #   ShamGene <- filter(results_genes_sham, -log2(fc) == x)
   #   ShamGeneName <- ShamGene$gene_name
   #   print(ShamGeneName)
   # })
   
}

# Run the application 
shinyApp(ui = ui, server = server)

