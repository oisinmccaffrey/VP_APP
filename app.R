library(shiny)
library(plotly) # interactive ggplots.. 
library(grid)
library(rmarkdown)
library(magrittr)
library(Cairo)
library(grDevices)
library(ggplot2)
library(vcfR) # For manipulating VCF data
library(dplyr) # Manipulating datasets
library(tidyr)
library(htmlwidgets) # to use saveWidget function
library(shinycssloaders) # for loading spinner while app is loading
library(stringr) # Used for removing pipe separators from VCF files "|" 
library(VariantAnnotation) # exploration and annotation of genetic variants
library(tibble)
library(shinyWidgets) # custom input controls and user interface components for 'Shiny' applications
library(shinydashboard) # For shiny dashboard layout
library(shinydashboardPlus) # Shiny dashboard with 'AdminLTE2' component
library(DT) # to generate interactive datatables
library(plotly) # to generate interactive plots
library(shinyalert) # for error alert
library(shinyjs) # for error message
library(shinyscreenshot) # Capture screenshots in 'Shiny' applications - downloaded as a PNG image
library(webshot)
library(rvest) # Scraping data from websites
library(lubridate) # fast and user friendly parsing of date-time data,

#Importing the data, readVcf from VariantAnnotation package

vcf <- readVcf("/Users/oisinmccaffrey/Desktop/R_Shiny_Summer/SRR.gavin_secondpass.vcf", "GRCh38")
vcf_1 <- as.data.frame(VariantAnnotation::fixed(vcf))
vcf_2 <- as.data.frame(VariantAnnotation::info(vcf))
vcf_3 <- as.data.frame(rowRanges(vcf))


# duplicate info here
vcf_master <- cbind(vcf_3, vcf_1, vcf_2)
# stage an empty DF which we will populate in the 'for' loop.
# If extracting another col, be sure to add it here
collect_ann <- data.frame(Allele=as.character(),
                          Consequence=as.character(),
                          IMPACT=as.character(),
                          Symbol=as.character(),
                          Gene=as.character())

# loop over every row in df
for(i in 1:nrow(vcf_master)){
  # grab the ann column
  ann <- vcf_master$ANN[i]
  # parse ann column
  ann <- str_split(ann, "\\|")
  # convert to DF 
  ann <- as.data.frame(t(unlist(ann)))
  # we want the first 5 cols in this case
  ann <- ann[,1:5]
  # rename to match the collect_ann df
  colnames(ann) <- c("Allele", 
                     "Consequence", 
                     "IMPACT", 
                     "Symbol", 
                     "Gene")
  
  # populate the collect ann df
  collect_ann <- rbind(collect_ann, ann)
}

# append to master df
vcf_master <- cbind(vcf_master, collect_ann)
# Then collect
collect_rlv <- data.frame(Status=as.character(),
                          SRR=as.character(),
                          Prediction=as.character(),
                          CADD_Relevance=as.character())

for(i in 1:nrow(vcf_master)){
  rlv <- vcf_master$RLV[i]
  rlv <- str_split(rlv, "\\|")
  rlv <- as.data.frame(t(unlist(rlv)))
  rlv <- rlv[,c(11,13,15,17)]
  colnames(rlv) <- c("Status", 
                     "SRR", 
                     "Prediction", 
                     "CADD_Relevance")
  
  collect_rlv <- rbind(collect_rlv, rlv)
}
vcf_master <- cbind(vcf_master, collect_rlv)

#Substitute for '|' and 'c' with a blank '' instead
vcf_master$Allele <- gsub('^c\\(|\\)$', '', vcf_master$Allele)

#Substitute numbers for '' instead

vcf_master$Allele <- gsub("[^A-Za-z0-9]", "", vcf_master$Allele)
vcf_master <- vcf_master %>% rownames_to_column(var = "ID_")

# Make the ID 'NA' if it starts with a number, (we want the variants starting with "rs.....")

is.na(vcf_master$ID_) <- 
  startsWith(vcf_master$ID_, "1") |
  startsWith(vcf_master$ID_, "2") |
  startsWith(vcf_master$ID_, "3") |
  startsWith(vcf_master$ID_, "4") |
  startsWith(vcf_master$ID_, "5") |
  startsWith(vcf_master$ID_, "6") |
  startsWith(vcf_master$ID_, "7") |
  startsWith(vcf_master$ID_, "8") |
  startsWith(vcf_master$ID_, "9") |
  startsWith(vcf_master$ID_, "10")|
  startsWith(vcf_master$ID_, "11")|
  startsWith(vcf_master$ID_, "12")|
  startsWith(vcf_master$ID_, "13")|
  startsWith(vcf_master$ID_, "14")|
  startsWith(vcf_master$ID_, "15")|
  startsWith(vcf_master$ID_, "16")|
  startsWith(vcf_master$ID_, "17")|
  startsWith(vcf_master$ID_, "18")|
  startsWith(vcf_master$ID_, "19")|
  startsWith(vcf_master$ID_, "20")|
  startsWith(vcf_master$ID_, "21")|
  startsWith(vcf_master$ID_, "22")|
  startsWith(vcf_master$ID_, "X")

#Rename seqnames
vcf_master <- dplyr::rename(vcf_master, Chr = seqnames)

##N.B. 

# With some of the variants in the file, the CADD score is not under the relevant
#'CADD SCALED Score' column.. therefore it requires extracting from the 'CADD_Relevance'
# Column which provides information on the variant.

#This extraction is performed using the package (stringr) and the function str_extract
#We are taking the numerical value after 'Score of' i.e. CADD score of.

vcf_master$CADD_Extract <- str_extract(vcf_master$CADD_Relevance, '(?i) (?<=score of\\D)\\d+')
vcf_master <- transform(vcf_master, CADD_Extract = as.numeric(CADD_Extract))
vcf_master$CADD_FINAL <- ifelse(is.na(vcf_master$CADD_Extract), vcf_master$CADD_SCALED, vcf_master$CADD_Extract)
vcf_master <- transform(vcf_master, CADD_FINAL = as.numeric(CADD_FINAL))

# transform CADD_SCALED column to be numeric data
vcf_master <- transform(vcf_master, CADD_FINAL = as.numeric(CADD_FINAL))
vcf_master <- transform(vcf_master, EXAC_AF = as.numeric(EXAC_AF))
vcf_master <- transform(vcf_master, ID_ = as.character(ID_))


#Reorder the ID_ column to be rank in terms of CADD_SCALED score
vcf_master$ID_ <- reorder(vcf_master$ID_, vcf_master$CADD_FINAL)


#Log scale ExAC AF
vcf_master$EXAC_AF <- log(vcf_master$EXAC_AF)


#Creating a gene dataframe for the genes table, selecting relevant columns
genes <- vcf_master %>% dplyr::select(c(Symbol, 
                                 Gene, 
                                 Chr, 
                                 start, 
                                 REF, 
                                 Allele, 
                                 Consequence,
                                 EXAC_AF,
                                 IMPACT))

#Renaming the columns to more intuitive names
genes <- dplyr::rename(genes, 
                       From = REF, 
                       To = Allele, 
                       HGNC = Symbol,
                       MAF = EXAC_AF)


#Rounding the CADD Scaled scores to be whole numbers

vcf_master$CADD_FINAL <- round(as.numeric(vcf_master$CADD_FINAL), digits=0)

#Substitute any "_" for blank ""
genes$Consequence <- gsub("_", " ", genes$Consequence)
#Substitute "&" for " and "
genes$Consequence <- gsub("&", " and ", genes$Consequence)

#Rounding the CADD values to nearest whole number using gsub

vcf_master$Consequence <- gsub("_", " ", vcf_master$Consequence)
vcf_master$Consequence <- gsub("&", " and ", vcf_master$Consequence)


#Now using gsub for Status column

vcf_master$Status <- gsub("SRR:", " ", vcf_master$Status)
vcf_master$Status <- gsub("_", " ", vcf_master$Status)

#Quality score dataframe

df<- as.data.frame(vcf_master$QUAL)
df <- rename(df, QUAL = `vcf_master$QUAL`)
df <- transform(df, QUAL = as.numeric(QUAL))
df <- df %>% drop_na(QUAL)

#Read depth dataframe
df_DP<- as.data.frame(vcf_master$DP)
df_DP <- rename(df_DP, DP = `vcf_master$DP`)
df_DP <- transform(df_DP, DP = as.numeric(DP))

#Mapping Quality dataframe

df_MQ<- as.data.frame(vcf_master$MQ)
df_MQ <- rename(df_MQ, MQ = `vcf_master$MQ`)
df_MQ <- transform(df_MQ, MQ = as.numeric(MQ))



### Beginning of Shiny App

# Define UI for application


ui = dashboardPage(controlbar = NULL, footer = NULL,
                   skin = "red",
                   dashboardHeader(title = "Variant Prioritisaion"),
                   
                   # Create side-bar menu with all tab options: 
                   dashboardSidebar(
                     
                     br(),
                     
                     #SidebarMenu to create Home, 
                     #Pathogenicity, Genomic Data, VCF METRICS, About, and Raw Data Tabs
                    
                     
                     sidebarMenu(
                       menuItem(
                         text = "Home",
                         tabName = "home",
                         icon = icon("home")
                       ),
                       menuItem(
                         text = "Pathogenicity", 
                         tabName = "pathogenicity",
                         icon = icon("medkit"),
                         startExpanded = TRUE,
                         menuSubItem("CADD", tabName = "CADD"),
                         menuSubItem("Status", tabName = "Status")),
                       
                       menuItem(
                         text = "Genomic Data", 
                         tabName = "genomicdata",
                         icon = icon("dna"),
                         startExpanded = TRUE,
                         menuSubItem("Gene Table", tabName = "table_gene"),
                         menuSubItem("Gene Overview", tabName = "download_gene")),
                       
                       menuItem(
                         text = "VCF Metrics", 
                         tabName = "metrics",
                         icon = icon("th"),
                         startExpanded = TRUE,
                         menuSubItem("Quality", tabName = "Quality"),
                         menuSubItem("Read Depth", tabName = "DP"),
                         menuSubItem("Mapping Quality", tabName = "MQ")),
                       
                       menuItem(
                         text = "Raw Data", 
                         tabName = "rawdata",
                         icon = icon("code"),
                         startExpanded = TRUE),
                       
                       menuItem(
                         text = "About", 
                         tabName = "About",
                         icon = icon("paper-plane"),
                         startExpanded = TRUE)
                       
                     )),
                   
                   
                   dashboardBody(
                     
    
                     
                     #css the .content-wrapper to change background colour
                     #This Makes the pain page white as oppose to standard grey
                     
                     tags$head(tags$style(HTML('/* body */
                                .content-wrapper, .right-side {
                                background-color: #FFFFFF;
                                }'))),
                     setShadow(class = "dropdown-menu"),
                     setShadow(class = "box"),
                     
                     tabItems(
                       
                       tabItem(tabName = "home",
                               tags$head(
                                 includeCSS("www/styles.css")
                               ),
                               
                               div(id = "page-topright",
                                   div(class = "source_link", 
                                       style="float:right", 
                                       a(href = "https://github.com/oisinmccaffrey/VP_APP", 
                                         "View source code on github", 
                                         icon("github"))),
                               ),
                               
                               
                               #Create Interactive 'More Information' button 
                               #Can be opened or closed
                               
                               HTML("<button type='button' class='btn' 
                               data-toggle='collapse' style='float:left' 
                     data-target='#app_info'><span class='glyphicon 
                     glyphicon-collapse-down'></span> More Information </button>"),
                               

                               br(),
                               br(),
                               
                               div(id = "app_info", class = "collapsible", 
                                   
                                   tags$h2(strong("Welcome to Variant Prioritisation")),
                                   
                                   br(),
                                   
                                   tags$h5(strong("Objective")), 
                                   "The primary objective of this clinical variant prioritisation web application 
                is to improve upon the limited functionality of GAVIN and other popular variant 
                calling/prioritisation tools. The application provides a tiered system ranking variants 
                in terms of their pathogenicity and role in known diseases, as well as providing additional 
                functionality such as summary reports of metrics used in variant classification 
                (e.g., quality, depth and allele frequency). 
                The application provides a clinician/user with the ability to filter by gene panels, 
                and query affected genes via the OMIM API, returning gene-phenotype information supported by literature. 
                As there is a veritable need for improvements to current variant prioritisation and visualisation methods, 
                this project will provide substantial positive progress to the movement 
                for improved sequence variant annotation and prioritisation from NGS projects.",
                                   tags$br(),tags$br(),
                                   
                                   tags$h5(strong("Data")),
                                   
                                   p("This application provides a means of visualising sequencing 
                    data post variant calling. Here we are using data taken 
                    from 1000 Genomes on GRCh38, specifically (Population:
                    British in England and Scotland, European Ancestry) See:"),             
                                   tags$a(href = "https://www.internationalgenome.org/", 
                                          "IGSR"),
                                   p(""),
                                   p("Raw whole exome sequencing (WES) data sourced from IGSR: 
                      The International Genome Sample Resource - "),
                                   tags$a(href = "https://www.internationalgenome.org/data-portal/sample", 
                                          "IGSR Data"),
                                   tags$br(),tags$br(),
                                   tags$img(src = "images/nuig_logo.png", width = "150px", height = "75px")
                               ),
                               
                               
                               
                               br(),  br(),
                               
                       ),
                       
                       #Use Javascript to detect the browser window size 
                       #(initial size and any resize if you make window bigger)
                       
                       tabItem("CADD",
                               fluidRow(
                                 tabPanel("CADD Plot",
                                          fluidPage(
                                            tags$head(tags$script('
                                            
                                            
                        var dimension = [0, 0];
                        $(document).on("shiny:connected", function(e) {
                        dimension[0] = window.innerWidth;
                        dimension[1] = window.innerHeight;
                        Shiny.onInputChange("dimension", dimension);
                        });
                        $(window).resize(function(e) {
                        dimension[0] = window.innerWidth;
                        dimension[1] = window.innerHeight;
                        Shiny.onInputChange("dimension", dimension);
                        });
                         ')), 
                                            
                                        
                                            
                                            plotlyOutput("plot2", width = "auto") %>% withSpinner(color = "red")
                                            

                                          )
                                 ))),
                       
                       #Use Javascript to detect the browser window size 
                       #(initial size and any resize if you make window bigger)
                       
                       tabItem("Status",
                               fluidRow(
                                 tabPanel("Status",
                                          fluidPage(
                                            tags$head(tags$script('
                                            
                        var dimension = [0, 0];
                        $(document).on("shiny:connected", function(e) {
                        dimension[0] = window.innerWidth;
                        dimension[1] = window.innerHeight;
                        Shiny.onInputChange("dimension", dimension);
                        });
                        $(window).resize(function(e) {
                        dimension[0] = window.innerWidth;
                        dimension[1] = window.innerHeight;
                        Shiny.onInputChange("dimension", dimension);
                        });
                         ')),
                                            
                                            plotlyOutput("plot3", width = "auto") %>% withSpinner(color = "red")
                                            
                                            
                                          )
                                 ))),
                       
                       tabItem("table_gene",
                               fluidRow(
                                 box(
                                   title = "Genes Table",
                                   closable = TRUE,
                                   width = 12,
                                   status = "warning",
                                   solidHeader = TRUE,
                                   collapsible = TRUE,
                                   label = boxLabel(
                                     text = 4,
                                     status = "danger",
                                     style = "circle"
                                   ),
                                   tabPanel("Genes",
                                            status = "success",
                                            dataTableOutput("genomic_plot") %>% withSpinner(color = "red"))))),
                       
                       
                       tabItem(
                         tabName = "download_gene",
                         box(
                           title = "Gene Overview",
                           width = 12,
                           status = "warning",
                           collapsible = TRUE,
                           solidHeader = TRUE,
                           type = NULL,
                           src = "shorturl.at/avDHP",
                           color = "aqua-active",
                           useShinyalert(),
                           # Search Bar for user to enter Gene Symbol for patient-level data
                           sidebarSearchForm(textId = "searchvcf", buttonId = "searchButton",
                                             label = "Enter Gene Symbol..."), 
                           br(),
                           boxProfileItem("Gene:", 
                                          span(textOutput("gene_symbol"), 
                                               style='color:black')),
                           br(),
                           boxProfileItem("Chromosome:", 
                                          span(textOutput("gene_chromosome"), 
                                               style='color:black')),
                           br(),
                           boxProfileItem("Reference base:", 
                                          span(textOutput("gene_ref_variant"), 
                                               style='color:black')),
                           br(),
                           boxProfileItem("Alternate base:", 
                                          span(textOutput("gene_alt_variant"), 
                                               style='color:black')),
                           br(),
                           boxProfileItem("Consequence:", 
                                          span(textOutput("gene_consequence"), 
                                               style='color:#0000FF')),
                           br(),
                           boxProfileItem("Minor Allele Frequency: 10^X:", 
                                          span(textOutput("gene_maf"), 
                                               style='color:#0000FF')),
                           br(),
                           boxProfileItem("Impact:", 
                                          span(textOutput("gene_impact"), 
                                               style='color:#ff2400')),
                           br(),
                           boxProfileItem("CADD:", 
                                          span(textOutput("gene_cadd"), 
                                               style='color:#ff2400')),
                           hr())),
                                                 
                       tabItem("rawdata",
                               fluidRow(
                                 tabPanel("Raw Data",
                                          fluidPage(
                                            dataTableOutput("plot1") %>% withSpinner(color = "red")
                                          )
                                 ))),
                       
                       
                       tabItem("About",
                               
                               div(id = "page-topright",
                                   div(class = "source_link", 
                                       style="float:right", 
                                       a(href = "https://github.com/oisinmccaffrey/VP_APP", 
                                         "View source code on github", 
                                         icon("github"))),
                                   
                               ),
                               br(),
                               br(),
                               br(),
                               includeHTML("about.html"),
                               shinyjs::useShinyjs(),
                               tags$head(
                                 tags$link(rel = "stylesheet", 
                                           type = "text/css", 
                                           href = "plugins/carousel.css"),
                                 tags$script(src = "plugins/holder.js")
                               ),
                               
                               
                               tags$style(type="text/css",
                                          ".shiny-output-error { visibility: hidden; }",
                                          ".shiny-output-error:before { visibility: hidden; }"
                               )),
                       
                       
                       tabItem("Quality",
                               
                               fluidRow(
                                 tabPanel("Quality",
                                          fluidPage(
                                            tags$head(
                                              includeCSS("www/styles.css")
                                            ),
                                            
                                            HTML("<button type='button' class='btn' data-toggle='collapse' 
                                                 style='float:left' data-target='#quality_info'
                                                 ><span class='glyphicon glyphicon-collapse-down'
                                                 ></span> More Information</button>"),
                                            
                                            br(),
                                            br(),
                                            
                                            
                                            div(id = "quality_info", class = "collapsible",
                                            
                                            tags$h5(strong("Quality (QUAL)")),
                                             
                                                p("Quality scores or ‘QUAL’ refers to the Phred-scaled probability 
                                                  that a polymorphism actually exists at the specific site 
                                                  (i.e. Ref allele to Alternative allele).
                                                  The data is presented as a histogram of all variants in the file, 
                                                  with the quality (QUAL) of each variant plotted  against  
                                                  the  relative  frequency  of  that  score.  The higher  the  score,  
                                                  the  higher  the  probability  that  thevariant call is correct 
                                                  at that position in the genome.",
                                                  br(),
                                                  tags$h5(strong("Numerical Example:")), 
          
                                                  p("QUAL = 20, there is a 1 % chance that there is no variant at the site"),
                                                  p("QUAL = 50, there is a 1 in 1e5 (0.00001%) chance that there is no variant at the site"),
                                                  tags$a(href = "https://gatk.broadinstitute.org/hc/en-us", "GATK")  
                                                ),
                                                
                                                
                                                
                                                br(),  br(),
                                                
                                            ),
                                            plotlyOutput("qual_hist_plot") %>% withSpinner(color = "red")

                                          )))),
                       
                       
                       tabItem("DP",
                               
                               fluidRow(
                                 tabPanel("DP",
                                          fluidPage(
                                            tags$head(
                                              includeCSS("www/styles.css")
                                            ),
                                            
                                            HTML("<button type='button' class='btn' data-toggle='collapse' 
                                                 style='float:left' data-target='#dp_info'
                                                 ><span class='glyphicon glyphicon-collapse-down'
                                                 ></span> More Information</button>"),
                                            
                                            br(),
                                            br(),
                                            
                                            
                                            div(id = "dp_info", class = "collapsible",
                                                
                                                tags$h5(strong("Read Depth (DP)")),
                                                
                                                p("The read depth or coverage refers to the amount 
                                                  of sequence data available around the site of the variant.
                                                  The read depth can be somewhat considered as 
                                                  reflective of the power one has to make a call on the 
                                                  genotype of the sample at this position. 
                                                  This is distinct from the QUAL calculations -- 
                                                  which filter out reads of poor mapping quality."),
                                                  br(),
                                                  
                                                  tags$a(href = "https://gatk.broadinstitute.org/hc/en-us", "GATK")  
                                                ),
                                            
                                            br(),  br(),
                                            
                                          ),
                                          plotlyOutput("dp_hist_plot") %>% withSpinner(color = "red")
                                 ))),
                       
                       
                       tabItem("MQ",
                               fluidRow(
                                 tabPanel("MQ",
                                          fluidPage(
                                            tags$head(
                                              includeCSS("www/styles.css")
                                            ),
                                            
                                            HTML("<button type='button' class='btn' data-toggle='collapse' 
                                                 style='float:left' data-target='#mq_info'
                                                 ><span class='glyphicon glyphicon-collapse-down'
                                                 ></span> More Information</button>"),
                                            
                                            br(),
                                            br(),
                                          
                                            div(id = "mq_info", class = "collapsible",
                                                
                                                tags$h5(strong("Mapping Quality (MQ)")),
                                                
                                                p("Median mapping quality score calculated taking 
                                                  into account each read that matches 
                                                  the allele from the variant. "),
                                                br(),
                                                
                                                tags$a(href = "https://gatk.broadinstitute.org/hc/en-us", "GATK")  
                                            ),
                                          ),
                                          plotlyOutput("MQ_hist_plot") %>% withSpinner(color = "red")
                                 )
                                 
                               )))))


# Define Server for application


server <- function(input, output, session) {
  
  observeEvent(input$searchButton, {
    # if empty/invalid gene symbol is entered into searchbar and action button is pressed display error message
    if (is.null(input$searchvcf)) {
      shinyalert("Oops!", "Please Enter a Gene Symbol", type = "error") 
    } else if(!(input$searchvcf %in% vcf_master$Symbol)){
      shinyalert("Invalid Gene Symbol Entered", "Please Enter a Valid Gene Symbol", type = "error")
    }
    
  }) 
  
  # Creating reactive functions for the 'Gene Overview' table underneath 
  # the 'Genomic Data' sidebar panel..
  
  vcf_genes <- reactive({
    
    # filter clinical data based on gene symbol entered in searchbar
    vcf_master %>% filter(Symbol == input$searchvcf) 
    
  })
  
  
  output$gene_symbol <- renderText({
    
    # display gene symbol corresponding to gene query in search bar
    vcf_gene <- vcf_genes()
    vcf_gene$Symbol <- as.character(vcf_gene$Symbol)
    
  })
  
  
  output$gene_chromosome <- renderText({
    
    # display chromosome for gene with corresponding gene symbol query in search bar
    vcf_gene <- vcf_genes()
    vcf_gene$Chr
    
  })
  
  
  output$gene_ref_variant <- renderText({
    
    # display reference variant for gene with corresponding gene symbol query in search bar
    vcf_gene <- vcf_genes()
    vcf_gene$REF
    
    
  })
  
  
  output$gene_alt_variant <- renderText({
    
    # display alternative variant for gene with corresponding gene symbol query in search bar
    vcf_gene <- vcf_genes()
    vcf_gene$Allele
    
  })
  
  
  output$gene_consequence <- renderText({
    
    # display the consequence of the variant for gene with corresponding gene symbol query in search bar
    vcf_gene <- vcf_genes()
    vcf_gene$Consequence
    
  })
  
  output$gene_maf <- renderText({
    
    # display the consequence of the variant for gene with corresponding gene symbol query in search bar
    vcf_gene <- vcf_genes()
    vcf_gene$EXAC_AF
    
  })
  
  output$gene_impact <- renderText({
    
    # display the impact for gene with corresponding gene symbol query in search bar
    vcf_gene <- vcf_genes()
    vcf_gene$IMPACT
    
  })
  
  output$gene_cadd <- renderText({
    
    # display CADD score for gene with corresponding gene symbol query in search bar
    vcf_gene <- vcf_genes()
    vcf_gene$CADD_Relevance
    
  })
  
  output$plot1 <- DT::renderDataTable({
    
    vcf_master
    
  })
  
  # Creating interactive plot for CADD vs MAF for each variant
  # The user can filter the interactive plot by the type of variant 
  # e.g. missense/stop gain. Moreover, one can hover over individual points 
  # (variants) to attain the relevant ExAC minor allele frequency, CADD score,  
  # Consequence  (e.g.  missense  variant)  and  
  # SNP ID (e.g. rs200203535) for each variant 
  
  observeEvent(input$dimension,{
    output$plot2 <- renderPlotly({
      
      cadd_id_plot <- ggplot(data = vcf_master, 
                           aes(x=EXAC_AF, y = as.factor(CADD_FINAL), fill = Consequence, ID = ID_)) + 
        geom_point(size = 2) +
        ggtitle("CADD Score vs. Minor Allele Frequency") +
                           
                           theme_minimal() +
                           theme(plot.margin=unit(c(1,1,1.5,1.2),"cm")) +
                           theme(
                             axis.line.y = element_blank(),
                             axis.ticks.y = element_blank(),
                             axis.text.x = element_text(),
                             panel.grid.major.x =element_line(),
                             panel.grid.minor.x = element_blank(),
                             panel.grid.major.y = element_blank(),
                             plot.background = element_rect(fill = "transparent",colour = NA),
                             panel.background=element_blank(),
                             panel.ontop = TRUE, 
                             legend.position= c(0.75, 0.92), legend.direction="horizontal",
                             legend.text = element_text(size = 6), 
                             legend.key.size = unit(0.6, "lines"),
                             legend.title = element_text(size =10))
      
      cadd_id_plot <- cadd_id_plot + ylab("CADD Score") + 
        theme(axis.title.y = element_text(angle = 0)) +
        theme_bw() +
        scale_y_discrete(expand = c(0.025,0.025)) +
        
        scale_x_continuous(limits = c(-12.5, -2.5),
                           expand = c(0,0),
                           breaks = seq(-12.5,-2.5, by = 1),
                           name = "Minor Allele Frequency: 10^X:") +
        
        
        xlab("Minor Allele Frequency: 10^X:") + 
        theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        theme(axis.title.x = element_text(angle = 0))
      
      
      #tooltip used to ensure we can hover each point and obtain relevant information
        
      ggplotly(cadd_id_plot, tooltip = c("x", "y", "ID", "fill"), width = (0.825*as.numeric(input$dimension[1])), 
               height = (0.90*as.numeric(input$dimension[2])))

      
      })
    
    
    # Creating interactive plot for the 'Status' of each variant
    # Informs us on the types of variants on each chromosome in the file.

      output$plot3 <- renderPlotly({
        
       ID_location_plot <- ggplot(vcf_master, aes(x=ID_, y = Chr, fill=Status, CADD = CADD_FINAL, Consequence = Consequence)) + 
          coord_flip() +
          geom_point(size = 2) +
         ggtitle("Genomic location vs. Variant Status") +
          theme(
            axis.title.x = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.ontop = TRUE, 
            panel.grid.major.x = element_line(colour = "grey90", linetype = "dashed", size=0.15),
            panel.grid.minor.x = element_blank(),
            panel.grid.major.y = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=0.5),
            plot.background = element_rect(fill = "transparent",colour = NA),
            plot.margin = unit(c(0.5, 0.5, 2, 0.5), "cm")) +
         theme(axis.title.x = element_text(angle = 0)) +
         theme(axis.text.y = element_blank()) +
         theme(axis.ticks.y = element_blank()) +
          
         scale_x_discrete("CADD Score") +
          scale_y_discrete("Genomic Location") +
         theme(plot.title = element_text(hjust = 0.5, face = "bold"))
         theme_bw() +
           theme(
           legend.position= c(0.75, 0.92), legend.direction="horizontal",
         legend.text = element_text(size = 6), 
         legend.key.size = unit(0.6, "lines"),
         legend.title = element_text(size =10))
                            
       
       ggplotly(ID_location_plot, tooltip = c("x", "y", "fill", "CADD", "Consequence"), width = (0.825*as.numeric(input$dimension[1])), 
                height = (0.90*as.numeric(input$dimension[2])))

      })
      
      
      # The following histograms are featured under 'VCF' metrics.
      # Creating histograms of the quality, read depth, and mapping quality 
      # of the VCF file.
    
    output$qual_hist_plot <- renderPlotly({
      
      qual_plot <- ggplot(df, aes(QUAL)) +
        geom_histogram(bins = 30, color = "black", fill = "#00FFFF") +
        geom_vline(aes(xintercept = mean(QUAL)), 
                   linetype = "dashed", size = 0.6) +
        scale_x_continuous(limits = c(0, 2000), 
                           expand = c(0, 0),
                           breaks = seq(0, 2000, by = 250),
                           name = "QUAL") +
        ggtitle("Quality (QUAL)") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ylab("Frequency") +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) 
      
      ggplotly(qual_plot)
      
    })
    
    output$dp_hist_plot <- renderPlotly({
      
      dp_hist <- ggplot(data=df_DP, aes(DP)) + 
        geom_histogram(bins = 16, color = "black", fill = "#5CD85A") +
        ggtitle("Read Depth (DP)") +
        theme_minimal() +
        geom_vline(aes(xintercept = mean(DP)), 
                   linetype = "dashed", size = 0.6) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ylab("Frequency") +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) 
      
      ggplotly(dp_hist)
      
    })
    
    output$MQ_hist_plot <- renderPlotly({
      
      MQ_hist <- ggplot(data=df_MQ, aes(MQ)) + 
        geom_histogram(bins = 16, color = "black", fill = "#0059b3") +
        ggtitle("Mapping Quality (MQ)") +
        theme_minimal() +
        scale_x_continuous(limits = c(55, 65), 
                           expand = c(0, 0),
                           breaks = seq(55, 65, by = 1),
                           name = "MQ") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
        ylab("Frequency") +
        theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) 
      
      ggplotly(MQ_hist)
      
    })
    
    output$genomic_plot <- DT::renderDataTable({
      genes
      
    })   
    
    
    #For downloading HTML report, not working atm
    output$report <- downloadHandler(
      # For PDF output, change this to "report.pdf"
      filename = "report.html",
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(qual_hist_plot = input$qual_hist_plot)
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv()))
      

    })
    
  })
}
shinyApp(ui = ui, server = server)