library(shiny)

fluidPage(
  tags$head(
    tags$style(
      HTML(
      "
        .front-page {
          text-align: center;
          justify-content: center;
          display: flex;
          align-items: center;
          flex-direction: column;
          height: 100vh; 
          background-image: url('https://drive.google.com/uc?id=1iVw5KQcDtWt_A44nCsoJYNUo11N-LIbb');
          background-size: cover;
          background-position: center;
          border: 2px solid purple;
        }
        .front-page h3 {
          font-size: 1.5em;
          text-align: right;
          margin-bottom: 10px;
        }
        .front-page img {
          max-width: 80%;
          height: auto;
          margin-top: 20px;
        }
        .run-primary{
        font-size: 20px;
        font-weight: bold;
        background-color: #660066;
        color: white;
        }
        .nav-tabs > li > a {
          background-color:   #660066;
          color: white;
          }
        .sidebar{
        background-color: #f5f5fd;
        border: 2px solid purple;
        }
        #clear{
        font-size: 14px;
        background-color: #e0cce0;
        }

        "
    ))
  ),
  tabsetPanel(
    tabPanel(
      "Front Page",
      div(class= "front-page",
          h1(strong("VCF VariantView: an RShiny tool for the functional analysis of variants")),
          h3("UOC - TFM Master in Bioinformatics and Biostatistics"),
          h3("Developed by Lídia Romero Cortadellas"),
          h3("Course 2023-2024"),
          img(src = "https://www.uoc.edu/portal/_resources/common/imatges/sala_de_premsa/noticies/2016/202-nova-marca-uoc.jpg", height = 250, width=250),

      )
    ),
    tabPanel(
      "Data upload and analysis",
      sidebarLayout(
        sidebarPanel(
          class = "sidebar",
          fluidRow(
            h3(strong("File upload")),
          ),
          fluidRow(
            column(12,
            ),
            fileInput("case", label="", multiple=F),
          ),
          
          fluidRow(
            h3(strong("Filtering")),
            
            #Gene 
            column(width=12,
                   h4("Gene"),
                   textInput("gene_name", label = "Introduce the gene",
                             value = ""
                   )),
            
            #Classification
            column(width=12,
                   h4("Classification"),
                   selectInput("classification", label = "Choose the variant classification",
                               choices = c("Frameshift deletion" = "frameshiftdel", "Frameshift insertion" = "frameshiftins", "Missense" = "Missense_Mutation", "Nonsense" = "Nonsense_Mutation", "Silent" = "Silent", "In frame insertion" = "inframeins", "In frame deletion" = "inframdel"),
                               multiple = TRUE),
            ),
            
            #Zygosity
            column(width=12,
                   h4("Zygosity"),
                   selectInput("genot", label = "Choose the zigosity",
                               choices = c("Heterozygous" = "het", "Homozygous" = "hom"),
                               multiple = TRUE)
            ),
            
            #Disease
            column(width=12,
                   h4("Disease"),
                   textInput("disease", label = "Introduce a keyword related to the disease",
                             value = ""
                   )),
            
            #Polyphen prediction
            column(width=12,
                   h4("Polyphen prediction"),
                   selectInput("polyphen", label = "Choose the polyphen prediction of the variant",
                               choices = c("Benign" = "benign", "Probably damaging" = "probdam", "Possibly damaging" = "possdam"),
                               multiple = TRUE),
            ),
            
            #Minor allele frequency
            column(width=12,
                   h4("Minor allele frequency"),
                   selectInput("comp_freq", label="Select the sign",
                               choices = c("Greater or equal" = ">=", "Smaller or equal" = "<=")
                   ),
                   numericInput("freq", label = "Introduce a numeric value", value = 0, min = 0, max = 1, step = 0.001
                   )),
          ),
          fluidRow(
            #Analysis and download search criteria
            h3(strong("Analysis")),
            column(width=12,
                   actionButton("btn", label = "Go!", class = "run-primary"),
                   br(), br(),
                   downloadButton("filt", label ="Save parameters used for this analysis in .csv"),
                   br(), br(),
                   
            ) 
          ),
        ),
        mainPanel(
          tabPanel(
            title = "Description",
            h4(strong("File upload")),
            p("The file must be .vcf and up to 7 GB. Only one file is allowed. Please wait until upload is completed. You can upload a new file at any time. You can also analyse the data without applying any filters."),
            h4(strong("Variant filtering and prioritization")),
            p(strong("Gene"), "Introduce the name of your gene of interest"),
            p(strong("Classification"), "Choose one or more variant classifications among: frameshift, nonsense, missense, silent or in frame insertion/deletion."),
            p(strong("Zygosity"), "Choose the zigosity you expect to see for the variant: homozygous, heterozygous or both."),
            p(strong("Disease"), "Introduce a keyword related to your disease of interest, e.g. anemia."),
            p(strong("Polyphen prediction"), "Choose one or more possibilities of polyphen prediction in pathogenicity for the variant: benign, possibly damaging or probably damaging."),
            p(strong("Minor allele frequency"), "Introduce a MAF numeric value and either you want to search for greater or lower numbers. Decimals are represented by a '.'."),
            p(strong("Other considerations"), "If you want to eliminate one of the chosen options, just select in and once it is marked in blue, press the supress button in your keyboard. You can use either capital or small letters."),
            h4(strong("Variant analysis and storage of parameters")),
            p("Press the Go! button to start the analysis. You can save the selected filter parameters in .csv format. The analysis may require some time (indicated at the bottom right panel), please be patient. Once the processing bar has disappeared, the results will be displayed in the corresponding tabs (tabular results and graphical representations) at the top."),
            h4(strong("Clear analysis")),
            p("Press the following button if you want to clear the analysis once performed."),
            actionButton("clear", "Clear analysis")
            )
        )
      )
    ),
    tabPanel(
      "Variant table",
      dataTableOutput("resultTable"),
      br(), br(),
      downloadButton("tabres", label ="Download the tabular results in .xlsx "),
      
    ),
    tabPanel(
      "Summary of the identified variants",
      plotOutput("mafgraph",height = "80vh", fill=TRUE),
    ),
    tabPanel(
      "Variant analysis results",
      fluidRow(
        column(width = 1, align = "center", actionButton("prev_btn", "◄")),
        column(width = 10, align = "center", plotOutput("grid_plot", height = "80vh", fill=TRUE)),
        column(width = 1, align = "center", actionButton("next_btn", "►")),
      ),
      br(), br(),
      downloadButton("graphres", label ="Download the graphical results in .PDF"),
      p("Generating the PDF may take some time. Please wait until the window for saving the file opens."),
      
    ),
    tabPanel(
      "Help",
      h4(strong("Description of the App")),
      p("This app was developed for variant analysis from a VCF (Variant Call Format) file and coupled with functional analysis and visualization as a powerful tool for researchers and geneticists. The app performs variant annotation and filtering according to the parameters given by the user; and displays the resulting data in both tabular and graphical formats. Several functional analyses are performed if possible and shown in form of graphs. You can download each type of results at the bottom of the table and graphs."),
      h4(strong("Variant table")),
      p("As a result of the analysis, a table with the filtered variants will be displayed. The table can be downloaded in the xlsx format. The following columns will be visible:"),
      p(strong("Variant"), "Complete variant data including chromosome, position and change."),
      p(strong("chromosome"), "Indicates the chromosome where the variant is located."),
      p(strong("start"), "Genomic position in which the variant starts."),
      p(strong("end"), "Genomic position in which the variant ends"),
      p(strong("ref_allele"), "Base/s that locate/s in that position originally."),
      p(strong("alt_allele"), "Base/s by which the original ones have been changed."),
      p(strong("rsID"), "Identifier for the SNP"),
      p(strong("GeneID"), "Entrez Gene identifier"),
      p(strong("Protein_posi"), "Number of amino acid that the variant is affecting."),
      p(strong("ref_AA"), "Amino acid/s that locate/s in that position originally."),
      p(strong("alt_AA"), "Amino acid/s by which the original ones have been changed."),
      p(strong("aaChange"), "Amino acid change in the position."),
      p(strong("SampleID"), "Sample ID."),
      p(strong("Variantloc"), "Variant localization throught the mRNA including: promoter, fiveUTR, coding, intron, threeUTR or intergenic."),
      p(strong("Gene"), "Gene name."),
      p(strong("OMIM"), "OMIM reference number for the gene of interest."),
      p(strong("GO"), "Gene Ontology number for the gene of interest."),
      p(strong("polyphenpredict"), "Polyphen prediction of the amino acid change within the protein: benign, possibly damaging or probably damaging."),
      p(strong("ENSEMBL"), "ENSEMBL numer for the gene of interest."),
      p(strong("GOterm"), "Specifies the branch of GO given in the column GOOntology."),
      p(strong("GOOntology"), "MF: molecular function; describes the elemental activities of a gene product at the molecular level. BP: biological process; defines sets of molecular events or interactions that contribute to the functioning of a cell or organism. CP: cellular component; specifies the locations or complexes where gene products are active."),
      p(strong("genotype"), "Each number represents an allele. 0: reference allele. 1: alternative allele. 2: second alternative allele. 0/0, 1/1 and 2/2: homozygotes for the respectives alleles. 0/1, 1/2, 0/2: heterozygots for the respectives alleles."),
      p(strong("ReadsRef"), "Number of reads of the reference allele."),
      p(strong("ReadsAlt"), "Number of reads of the alternative allele."),
      p(strong("DiseaseOMIM"), "Names of the diseases related to the gene according to OMIM."),
      p(strong("Type"), "Type of base change: Deletion, Insertion os SNP."),
      p(strong("MAF"), "Minor allele frequency for the alterative allele."),
      p(strong("Variant_Classification"), "Type of amino acid mutation: Silent, Missense, Nonsense, Frameshift (specifying either deletion or insertion), In frame changes (deletions and insertions) and Missing Classification (no information available)."),
      h4(strong("Summary of the identified variants")),
      p("Includes several plots that give an overview of the filtered variants."),
      h4(strong("Variant analysis results")),
      p("Includes several plots with an insight in the variant analysis, which can all be downloaded together in pdf format:"),
      p(strong("Plot 1. Non-silent variants per gene."), "Displays the number of non-silent variants per gene."),
      p(strong("Plot 2. Location of variants within the mRNA."), "Identifies the position of variants in the mRNA including promoter, 5'UTR, coding, 3'UTR or intergenic."),
      p(strong("Plot 3. Lollipop plot/s of the variants."), "Generates a lollipop plot per each chromosome, depicting the type of mutation throught them (deletion, insertion or SNP)."),
      p(strong("Plot 4. GO enrichment analysis results."), "Multiple plot with the results from the GO enrichment analysis. A: dotplot del resultat d’enriquiment GO. B: representació gràfica de la ontologia gènica en l’anàlisi d’enriquiment. C i D: Linkage entre gens i termes GO."),
      p(strong("Plot 5. KEGG enrichment analysis results."), "Results from the KEGG enrichment analysis."),
      p(strong("Plot 6. REACTOME enrichment analysis results."), "Results from the REACTOME enrichment analysis."),
      h4(strong("Other considerations")),
      p("Some data are retrieved from online databases and therefore the functionality of the app depends on the proper functioning of such servers and of the internet connection."),
    )
  )
)







