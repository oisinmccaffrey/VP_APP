library(shiny)
library(ggplot2)
library(vcfR)
library(dplyr)
library(tidyr)
library(stringr)
library(VariantAnnotation)
library(tibble)
library(shinyWidgets)
library(shinydashboard) # For shiny dashboard layout
library(shinydashboardPlus) # For shiny dashboard layout
library(DT) # to generate interactive datatables
library(plotly) # to generate interactive plots
library(shinycssloaders) # for loading spinner while app is loading
library(randomForest) # to generate RSF classifier models
library(gprofiler2) # for Functional Enrichment analysis of drive genes
library(shinyalert) # for error alert
library(shinyjs) # for error message

#ui.R
library(VariantAnnotation)
vcf <- readVcf("/Users/oisinmccaffrey/Desktop/Masters/SRR/SRR.gavin_secondpass.vcf", "GRCh38")
vcf_1 <- as.data.frame(VariantAnnotation::fixed(vcf))
vcf_2 <- as.data.frame(VariantAnnotation::info(vcf))
vcf_3 <- as.data.frame(rowRanges(vcf))
# duplicate info here, will let you decide whats important
vcf_master <- cbind(vcf_3, vcf_1, vcf_2)
# stage an empty DF which we will populate in the for loop.
# If you want to extract another col, be sure to add it here
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
    colnames(ann) <- c("Allele", "Consequence", "IMPACT", "Symbol", "Gene")
    # populate the collect ann df
    collect_ann <- rbind(collect_ann, ann)
}
# append to master df
vcf_master <- cbind(vcf_master, collect_ann)
# I do not know the official terms for these columns, figure that out
collect_rlv <- data.frame(Status=as.character(),
                          METIN=as.character(),
                          Prediction=as.character(),
                          CADD_Relevance=as.character())
for(i in 1:nrow(vcf_master)){
    rlv <- vcf_master$RLV[i]
    rlv <- str_split(rlv, "\\|")
    rlv <- as.data.frame(t(unlist(rlv)))
    rlv <- rlv[,c(11,13,15,17)]
    colnames(rlv) <- c("Status", "METIN", "Prediction", "CADD_Relevance")
    collect_rlv <- rbind(collect_rlv, rlv)
}
vcf_master <- cbind(vcf_master, collect_rlv)
vcf_master$Allele <- gsub('^c\\(|\\)$', '', vcf_master$Allele)
vcf_master$Allele <- gsub("[^A-Za-z0-9]", "", vcf_master$Allele)
vcf_master <- vcf_master %>% rownames_to_column(var = "ID_")
is.na(vcf_master$ID_) <- startsWith(vcf_master$ID_, "1") | startsWith(vcf_master$ID_, "2") | startsWith(vcf_master$ID_, "3") | startsWith(vcf_master$ID_, "4") | startsWith(vcf_master$ID_, "5") | startsWith(vcf_master$ID_, "6") | startsWith(vcf_master$ID_, "7") | startsWith(vcf_master$ID_, "8") | startsWith(vcf_master$ID_, "9") | startsWith(vcf_master$ID_, "10") | startsWith(vcf_master$ID_, "11") | startsWith(vcf_master$ID_, "12") | startsWith(vcf_master$ID_, "13") | startsWith(vcf_master$ID_, "14") | startsWith(vcf_master$ID_, "15") | startsWith(vcf_master$ID_, "16") | startsWith(vcf_master$ID_, "17") | startsWith(vcf_master$ID_, "18") | startsWith(vcf_master$ID_, "19") | startsWith(vcf_master$ID_, "20") | startsWith(vcf_master$ID_, "21") | startsWith(vcf_master$ID_, "22") | startsWith(vcf_master$ID_, "X")
vcf_master <- vcf_master %>% drop_na(ID_)
vcf_master <- vcf_master %>% drop_na(CADD_SCALED)
vcf_master <- transform(vcf_master, CADD_SCALED = as.numeric(CADD_SCALED))
vcf_master$ID_ <- reorder(vcf_master$ID_, vcf_master$CADD_SCALED)

genes <- vcf_master %>% select(c(IMPACT, Symbol, Gene))




ui = dashboardPage(
    skin = "red",
    dashboardHeader(title = "Variant Prioritisation"),
    # Create side-bar menu with all tab options: 
    dashboardSidebar(
        useShinyalert(),
        sidebarSearchForm(textId = "searchpatient", buttonId = "searchButton",
                          label = "Enter Patient ID..."), # Search Bar for user to enter Patient ID for patient-level data 
        br(),
        sidebarMenu(
            menuItem(
                text = "Home",
                tabName = "home",
                icon = icon("info-circle")
            ),
            menuItem(
                text = "Pathogenicity", 
                tabName = "pathogenicity",
                icon = icon("notes-medical"),
                startExpanded = TRUE,
                menuSubItem("Variant", tabName = "Variant")),
            
            menuItem(
                text = "Genomic Data", 
                tabName = "genomicdata",
                icon = icon("dna"),
                startExpanded = TRUE,
                menuSubItem("Genes", tabName = "rawdata")),

            menuItem(
                text = "At it", 
                tabName = "boxelements",
                icon = icon("th"),
                startExpanded = TRUE,
                menuSubItem("Clinical Prediction", tabName = "clinical_ten")))
    ),
    dashboardBody(
        setShadow(class = "dropdown-menu"),
        setShadow(class = "box"),
        
    tabItems(
            
        tabItem(tabName = "home",
                box(
                    title = h3("Welcome to Variant Prediction"), 
                    closable = FALSE, 
                    width = 12,
                    status = "info", 
                    solidHeader = TRUE, 
                    collapsible = TRUE,
                    h4("An integrated interactive shiny dashboard for visualisation and analysis of the Clinical and Genomic METABRIC Data and generation of ten-year breast cancer survival predictive models"),
                    accordion(
                        accordionItem(
                            id = 1,
                            title = "Clinical Data",
                            color = "info",
                            collapsed = FALSE,
                            h4("Lorem ipsum"),
                            br(),
                            h4("Lorem ipsum"),
                            br(),
                            h4("Lorem ipsum"),
                            hr(),
                            h4("Lorem ipsum"),
                            h4(tags$li("Lorem ipsum")),
                            h4(tags$li("Lorem ipsum", em("Lorem ipsum"))),
                            h4(tags$li("Lorem ipsum", em("Lorem ipsum"))),
                            h4(tags$li("Lorem ipsum", em("Lorem ipsum"))),
                            h4(tags$li("Lorem ipsum", em("Lorem ipsum"))),
                            hr(),
                            
                        accordionItem(
                            id = 2,
                            title = "Access the Code",
                            color = "info",
                            collapsed = FALSE,
                            socialButton(
                                href = "https://github.com",
                                icon = icon("github")))
                            
                        )
                    )
                ),
        ),
        
            tabItem("Variant",
                    fluidRow(
                        tabBox(
                            # Title can include an icon
                            tabPanel("Raw VCF Data",
                                     dataTableOutput("plot1")
                            ),
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
                                         plotlyOutput("plot2", width = "auto")
                                     )
                            )))),
                            
            tabItem("rawdata",
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
                                     dataTableOutput("genomic_plot"))))))))
                       


server <- function(input, output, session) {
    observeEvent(input$searchButton, {
        # if empty/invalid patient ID entered into searchbar and action button is pressed display error message
        if (is.null(input$searchpatient)) {
            shinyalert("Oops!", "Please Enter a Valid Patient ID", type = "error") 
        } else if(!(input$searchpatient %in% vcf_master$ID_)){
            shinyalert("Invalid Patient ID Entered", "Please Enter a Valid Patient ID", type = "error")
        }
    }) 
    output$plot1 <- DT::renderDataTable({
        vcf_master
    })
    observeEvent(input$dimension,{
        output$plot2 <- renderPlotly({
            cadd_id_plot <- ggplot(data = vcf_master, aes(x=ID_, y = CADD_SCALED, fill = Consequence)) + 
                geom_point(size = 2) +
                ggtitle("CADD vs. SNP Accession") +
                scale_y_discrete("CADD Score") +
                theme_minimal() +
                theme(
                    axis.line.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.line.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.text.x = element_text(),
                    axis.title.x = element_blank(),
                    panel.grid.major.x =element_line(),
                    panel.grid.minor.x = element_blank(),
                    panel.grid.major.y = element_blank())
                    
            ggplotly(cadd_id_plot, width = (0.95*as.numeric(input$dimension[1])), height = as.numeric(input$dimension[2]))
        })
        
        
        output$genomic_plot <- DT::renderDataTable({
            genes
            
    })   
        
    })
}
shinyApp(ui = ui, server = server)