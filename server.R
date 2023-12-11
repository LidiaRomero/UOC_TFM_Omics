packages_to_check <- c(
  "shiny",
  "BiocManager",
  "VariantAnnotation",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "AnnotationDbi",
  "SNPlocs.Hsapiens.dbSNP144.GRCh37",
  "BSgenome.Hsapiens.UCSC.hg19",
  "PolyPhen.Hsapiens.dbSNP131",
  "GO.db",
  "biomaRt",
  "GOstats",
  "topGO",
  "ReactomePA",
  "maftools",
  "clusterProfiler",
  "openxlsx",
  "gridExtra",
  "ggplot2",
  "cowplot"
)
for (package in packages_to_check) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

library(lattice)
library(BiocManager)
library(VariantAnnotation)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(PolyPhen.Hsapiens.dbSNP131)
library(GenomicFeatures)
library(vcfR)
library(GO.db)
library(tidyr)
library(biomaRt)
library(devtools)
library(GOstats)
library(topGO)
library(clusterProfiler)
library(ggplot2)
library(ReactomePA)
library(maftools)
library(openxlsx)
library(gridExtra)
library(ggplot2)
library(cowplot)

options(shiny.maxRequestSize=7000*1024^2)

source("vcfinitialprocessing.R")
source("annotation_1.R")
source("annotation_2.R")
source("annotation_3.R")
source("mafgeneration.R")
source("GOanalysis.R")
source("KEGGanalysis.R")
source("REACTOMEanalysis.R")


function(input, output, session) {

  # Process the data
  processedData <- reactiveValues(matS = NULL)
    observeEvent(input$btn,{
    if(!is.null(input$case$datapath)){
      withProgress(message = "Processing the file...", value = 0, {
        matS <- vcfinitialprocessing(input$case$datapath)
        incProgress(0.25, detail = "Step 1 completed")
        matA <- annotation_1(input$case$datapath)$matA
        variant_location_info <- annotation_1(input$case$datapath)$variant_location_info
        incProgress(0.5, detail = "Step 2 completed")
        matS <- annotation_2(input$case$datapath,matS, matA, variant_location_info)
        incProgress(0.75, detail = "Step 3 completed")
        last_annotation <- annotation_3(matS)
        incProgress(1, detail = "Step 4 completed")
        matS <- last_annotation
        processedData$matS <- matS
        uploadedFilePath(input$case$datapath)
      })
  }
    })
    
  #Filtering
  filterData <- function(matS){
    filtered <- reactive({
      filtered_matS <- matS
      # Gene
      if (input$gene_name != ""){
        filtered_matS <- filtered_matS[grep(input$gene_name, filtered_matS$Gene, ignore.case = TRUE), ]
      }
      #Classification
      if(!is.null(input$classification)){
      classification_map <- c("frameshiftdel" = c("Frame_Shift_Del"),
                                "frameshiftins" = c("Frame_Shift_Ins"),
                                "Missense_Mutation" = c("Missense_Mutation"),
                                "Nonsense_Mutation" = c("Nonsense_Mutation"),
                                "Silent" = c("Silent"),
                                "inframeins" = c("In_Frame_Ins"),
                                "inframdel" = c("In_Frame_Del"))
      filter_values <- unlist(classification_map[input$classification])
      filtered_matS <- filtered_matS[filtered_matS$Variant_Classification %in% filter_values, ]
      }
      #Zygosity
      if (!is.null(input$genot)){
        if ("het" %in% input$genot && "hom" %in% input$genot) {
          all_patterns <- c("0/0", "1/1", "2/2", "0/1", "1/0", "0/2", "2/0", "0/1/2")
          filtered_matS <- filtered_matS[filtered_matS$genotype %in% all_patterns, ]
        } else {
          if ("het" %in% input$genot) {
            heterozygous_patterns <- c("0/1", "1/0", "0/2", "2/0", "0/1/2")
            filtered_matS <- filtered_matS[filtered_matS$genotype %in% heterozygous_patterns, ]
          }
          if ("hom" %in% input$genot) {
            homozygous_patterns <- c("0/0", "1/1", "2/2")
            filtered_matS <- filtered_matS[filtered_matS$genotype %in% homozygous_patterns, ]
          }
        }
      }
      #Disease
      if(input$disease != ""){
        filtered_matS <- filtered_matS[grep(paste0("\\b", input$disease, "\\b"), filtered_matS$DiseaseOMIM, ignore.case = TRUE), ]
      }
      #Polyphen prediction
      if (!is.null(input$polyphen)) {
        polyphen_filters <- c()
        if ("benign" %in% input$polyphen) {
          polyphen_filters <- c(polyphen_filters, "benign")
        }
        if ("probdam" %in% input$polyphen) {
          polyphen_filters <- c(polyphen_filters, "probably damaging")
        }
        if ("possdam" %in% input$polyphen) {
          polyphen_filters <- c(polyphen_filters, "possibly damaging")
        }
        
        filtered_matS <- filtered_matS[filtered_matS$polyphenpredict %in% polyphen_filters, ]
      }
    #MAF
      if (!is.null(input$comp_freq) && !is.null(input$freq) && is.null(input$type) && is.null(input$classification)) {
        input_freq <- as.numeric(gsub(",", ".", input$freq))
        matS_MAF_numeric <- as.numeric(as.character(filtered_matS$MAF))
        
        if (input$comp_freq == ">=") {
          filtered_matS <- filtered_matS[(is.na(matS_MAF_numeric) | matS_MAF_numeric >= input_freq), ]
        } else if (input$comp_freq == "<=") {
          filtered_matS <- filtered_matS[(is.na(matS_MAF_numeric) | matS_MAF_numeric <= input_freq), ]
        }
      }
      return(filtered_matS)
    })
    filtered()
  }
  
  # Download filter data
  applied_filters <- reactive({
    filter_params <- data.frame(
      Parameter = c("Gene Name", "Classification", "Genotype", "Disease", "Polyphen", "MAF"),  
      Value = c(input$gene_name, paste(input$classification, collapse = ", "), paste(input$genot, collapse = ", "), input$disease, paste(input$polyphen, collapse = ", "), input$comp_freq)  
    )
    filter_params
  })
  output$filt <- downloadHandler(
    filename = function() {
      "applied_filters.csv"
    },
    content = function(file) {
      write.csv(applied_filters(), file, row.names = FALSE)
    }
  )  
  
  # Generate results table
  output$resultTable <- renderDataTable({
    filtered_matS <- filterData(processedData$matS)
    filtered_matS
  })
  
  # Clear the processed data
  uploadedFilePath <- reactiveVal(NULL)
  clearUploadedFile <- function() {
    uploadedFilePath(NULL)  
    processedData$matS <- NULL  
  }
  observeEvent(input$clear, {
    clearUploadedFile()  
    session$sendInputMessage("case", list(files = list(NULL)))
  })
  
  # Download the tabular results
  filtered_matS <- reactive({
    filterData(processedData$matS)  
  })
  output$tabres <- downloadHandler(
    filename = function() {
      "filtered_data.xlsx"
    },
    content = function(file) {
      write.xlsx(filtered_matS(), file, rowNames = FALSE)
    }
  )
  
  # Generate graphical results
  filtered_matS_reactive <- reactive({
    filterData(processedData$matS)  
  })
  plot1 <- reactive({
    filtered_data <- filtered_matS_reactive()
    if (!is.null(filtered_data)) {
      maf_object <- mafgeneration(filtered_data)
      plotmafSummary(maf = maf_object, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, fs = 1.1, titleSize = c(1.25,1.25))
    } else {
      ggplot() + ggtitle("No data available for this plot")
    }
  })
  plot2 <- reactive({
    filtered_data <- filtered_matS_reactive()
    if (!is.null(filtered_data)) {
    filtered_non_silent <- filtered_data %>% filter(!is.na(Gene), Variant_Classification != "Silent") %>% group_by(Gene) %>% summarise(Frequency = n()) %>% arrange(desc(Frequency)) %>% slice_head(n=6)
    ggplot(filtered_non_silent, aes(x = reorder(Gene, - Frequency ),y = Frequency, fill = Gene )) +
      geom_bar(stat = "identity") +
      labs(title = "", x = "Gene", y = "Count") + 
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 16))
    }else{
      ggplot() + ggtitle("No data available for this plot")
    }
  })
  plot3 <- reactive({
    filtered_data <- filtered_matS_reactive()
    if (!is.null(filtered_data)) {
    filtered_by_variant_loc <- filtered_data %>% group_by(VariantLoc) %>% summarise(Frequency = n())
    ggplot(filtered_by_variant_loc, aes(x = reorder(VariantLoc, -Frequency), y = Frequency, fill = VariantLoc)) +
      geom_bar(stat = "identity") +
      labs(title = "", x = "Variant Location", y = "Count") +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 16))  
    }else{
      ggplot() + ggtitle("No data available for this plot")
    }
  })
  plot4 <- reactive({
    filtered_data <- filtered_matS_reactive()
    if (!is.null(filtered_data)) {
    filtered_data$start <- as.integer(filtered_data$start)
    filtered_data$end <- as.integer(filtered_data$end)
    ggplot(filtered_data, aes(x=chromosome, y=start, color=Type)) +
      geom_segment(aes(xend=chromosome, yend = 0), linewidth = 2) +
      geom_point(stroke = 5, shape = "|") +
      labs(title = "", x = "Chromosome", y ="Position") + 
      theme_minimal() + 
      coord_flip() +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16)
      )
    }else{
      ggplot() + ggtitle("No data available for this plot")
    }
  }) 
  plot5 <- reactive({
    filtered_data <- filtered_matS_reactive()
    if (!is.null(filtered_data)) {
      params <- GOanalysis(filtered_data)
      if(nrow(params) >0) {
      p1 <- dotplot(params, showCategory=10)
      p2 <- goplot(params, showCategory=4, cex=0.5)
      p3 <- cnetplot(params, node_label = "category",
                     cex.params = list(category_label = 0.8))
      p4 <- cnetplot(params,node_label="gene",
                     cex.params = list(category_label = 0.5))
      plot_grid <- cowplot::plot_grid(p1, p2, p3, p4, ncol = 2, labels = LETTERS[1:4])
      plot_grid <- plot_grid + theme(axis.text = element_text(size = 12),   
                                     plot.margin = unit(rep(1, 4), "cm"))  
      return(plot_grid)
      } else {
        ggplot() +
          annotate("text", x = 0, y = 0, label = "No enriched GO pathways found", size = 6, vjust = 0.5, hjust = 0.5)
        }
    }else{
      ggplot() + ggtitle("No data available for this plot")
    }
  })
  plot6 <- reactive({
    filtered_data <- filtered_matS_reactive()
    if (!is.null(filtered_data)) {
      KEGGparams <- KEGGanalysis(filtered_data)
      if (nrow(KEGGparams) > 0) {
        ggplot(KEGGparams, aes(x = -log10(pvalue), y = ID)) +
          geom_point(stroke = 3) +
          geom_text(aes(label = Description), hjust = -0.1) +
          labs(x = "-log10(P-value)", y = "Pathway ID")
      } else {
        ggplot() +
          annotate("text", x = 0, y = 0, label = "No enriched KEGG pathways found", size = 6, vjust = 0.5, hjust = 0.5)
      }
    } else {
      ggplot() + ggtitle("No data available for this plot")
    }
  })
  plot7 <- reactive({
    filtered_data <- filtered_matS_reactive()
    if (!is.null(filtered_data)) {
      react <- REACTOMEanalysis(filtered_data)
      if (nrow(react) > 0) {
        ggplot(react, aes(x = -log10(pvalue), y = geneID)) +
          geom_point(stroke = 3) +
          geom_text(aes(label = Description), hjust = -0.1) +
          labs(x = "-log10(P-value)", y = "Pathway geneID")
      } else {
        ggplot() +
          annotate("text", x = 0, y = 0, label = "No enriched REACTOME pathways found", size = 6, vjust = 0.5, hjust = 0.5)
      }
    } else {
      ggplot() + ggtitle("No data available for this plot")
    }
    }) 
  
  output$mafgraph <- renderPlot({
    plot1()
  })
  
  
  plots <- list(plot2, plot3, plot4, plot5, plot6, plot7)
  current_plot <- reactiveVal(1)
  plot_titles <- c("Plot 1. Non-silent variants per gene", "Plot 2. Location of variants within the mRNA", "Plot 3. Lollipop plot/s of the variants", "Plot 4. GO enrichment analysis results", "Plot 5. KEGG enrichment analysis results", "Plot 6. REACTOME enrichment analysis results")
  output$grid_plot <- renderPlot({
    plot <- plots[[current_plot()]]()
    title <- plot_titles[current_plot()]
    title_label <- draw_label(title, fontface = 'bold', size = 16)  
    plot_with_title <- cowplot::plot_grid(plot, ggdraw() + title_label, nrow = 2, rel_heights = c(0.9, 0.1))
    plot_with_title
  })
  
  observeEvent(input$prev_btn, {
    if (current_plot() > 1) {
      current_plot(current_plot() - 1)
    } else {
      current_plot(length(plots))
    }
  })
  
  observeEvent(input$next_btn, {
    if (current_plot() < length(plots)) {
      current_plot(current_plot() + 1)
    } else {
      current_plot(1)
    }
  })
  # Download the graphical results
  generateAllPlotsPDF <- function(){
    pdf("analysis_plots.pdf", width = 12, height = 10) 
    for (i in 1:length(plots)) {
      print(plots[[i]]())
    }
    dev.off()
  }
  output$graphres <- downloadHandler(
    filename = function(){
      "analysis.plots.pdf"
    },
    content = function(file) {
      generateAllPlotsPDF()
      file.rename("analysis_plots.pdf", file)
    }
    
  )
}










