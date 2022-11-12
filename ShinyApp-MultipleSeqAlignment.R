
########################################
#..............BCB Diploma ............#
#.............DAV - R Course ..........#
#.............MSA Shiny App............#
#...............29/10/2020.............#
#............Maha Abdelrhman...........#
#.............Yasmine Saied............#
########################################

#sessionInfo()
# R version 4.2.0 (2022-04-22 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19043)
# 
# Matrix products: default
# 
# locale:
# [1] LC_COLLATE=English.utf8  LC_CTYPE=English.utf8    LC_MONETARY=English.utf8
# [4] LC_NUMERIC=C                     LC_TIME=English.utf8    
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] ggmsa_1.2.3          msa_1.28.0           Biostrings_2.64.1    GenomeInfoDb_1.32.4  XVector_0.36.0      
# [6] IRanges_2.30.1       S4Vectors_0.34.0     BiocGenerics_0.42.0  BiocManager_1.30.18  shinydashboard_0.7.2
# [11] rlang_1.0.6          shiny_1.7.2          lme4_1.1-30          Matrix_1.5-1       

# loaded via a namespace (and not attached):
#   [1] compiler_4.2.0 tools_4.2.0   


# Dependencies to be installed first ---
################################################################################
#install.packages("shiny")
#install.packages("rlang")
#install.packages("dplyr")
#install.packages("shinycssloaders")
#install.packages("BiocManager")
#BiocManager::install("ggmsa")
#BiocManager::install("msa")

########################## Load Libraries ######################################


require("shiny")
require("rlang")
require("shinydashboard")
require("BiocManager")
require('ggmsa')
require("ggplot2")
require("ggseqlogo")
library(dplyr)
library(shinycssloaders)
library(ggmsa)
library(msa)
library(ggplot2)
library(ggseqlogo)

# Required external functions
printSplitString <- function(x, width=getOption("width") - 1){
  
  starts <- seq(from=1, to=nchar(x), by=width)
  for (i in 1:length(starts))
    cat(substr(x, starts[i], starts[i] + width - 1), "\n")
}

data("BLOSUM80")

################################################################################

if (interactive()) {
  
  #define UI
ui <- dashboardPage(
  dashboardHeader(title="Multiple Sequence Alignment"),
  ###point 1.
  dashboardSidebar(sidebarMenuOutput("Semi_collapsible_sidebar")),
  ###point 2.
  dashboardBody(tags$head(tags$link(rel = "stylesheet",
                                    type = "text/css", href = "style.css"))),
  
  ###############################################################################
  
  
  body = dashboardBody(

    sidebarLayout(position = "left",
      sidebarPanel(
        conditionalPanel(condition = "input.tabselected==1",
        box(title = "Sequences", width=40 ,status = "primary", solidHeader = TRUE, collapsible = FALSE,
          
        

        fileInput("file1", label = "upload file", width="100%", accept = c(".fa", ".fasta"), 
                                  multiple = FALSE, buttonLabel = "Browse..."),
        
        # Insert action button for posting data  ----
        actionButton("button1", "Submit", class = "btn btn-primary", style = "color: white; background-color: #5c90c4"))),
        
        
        conditionalPanel(condition = "input.tabselected==2",
                         
          box(
            title = "MSA", width=30 ,status = "primary", solidHeader = TRUE, collapsible = FALSE,
              
            selectInput("selectMSA", "select function", width="100%",
                        c(Muscle = "Muscle", ClustalOmega = "ClustalOmega", ClustalW = "ClustalW"), 
                        selected = "ClustalW", multiple = FALSE ),
                        
            checkboxInput("check1", "Complete", value = FALSE),
        
        
            actionButton("button2", label = "Show",  class = "btn-success", 
                          style = "color: white; background-color: #5c90c4"),
            div(style="display:inline-block",downloadButton('downloadData1', 'Download MSA')),
            
            )),
        
        
        
        conditionalPanel(condition = "input.tabselected==3",
                         
        
        box(
          title = "Consensus matrix", width=30, 
          status = "primary", solidHeader = TRUE, collapsible = FALSE,
          
          actionButton("button3", "Compute", icon("paper-plane"))),
        
        
        box(
          title = "Consensus MSA", width=30, 
          status = "primary", solidHeader = TRUE, collapsible = FALSE,
          
          column(6, numericInput("num1", label = h5("Max threshold"), value = 1)),
          
          column(6, numericInput("num2", label = h5("Min threshold"), value = 1)),
          
          actionButton("button4", "Compute", icon("paper-plane")))),
        
        
        
        conditionalPanel(condition = "input.tabselected==4",
                         
          box(
            title = "MSA Plot", width=30, 
            status = "primary", solidHeader = TRUE, collapsible = FALSE,
                           
            column(6, numericInput("num3", label = h5("Start Position"), value = 0)),
                           
            column(6, numericInput("num4", label = h5("End Position"), value = 1)),
                           
          selectInput("selectcolor", "Color residues by: ", width="100%", 
                      c(Chemistry = "Chemistry_AA", Zappo = "Zappo_AA", Letter = "LETTER", 
                        Clustal = "Clustal"), 
                      selected = "Chemistry_AA", multiple = FALSE ), 
          
          actionButton("button6", "Plot", icon("paper-plane")),
                          
          hr(),
          selectInput("annotate", "Annotation modules", width = "100%", 
                      c("Logo", "Consensus"), 
                      selected = "Logo", multiple = FALSE),
                           
          actionButton("button7", "Annotation", icon("paper-plane")),
          
          hr(),
          radioButtons("var1", label = "Select file type:", choices = list("png", "pdf")),
          div(style="display:inline-block",downloadButton('downloadData3', 'Download Plot')),
          
          )),
        
        
        conditionalPanel(condition = "input.tabselected==5",

        box(
          title = "Download", width=30,
          status = "primary", solidHeader = TRUE, collapsible = FALSE,
          div(style="display:inline-block",downloadButton('downloadData2', 'Download Data'))
          
          ))
  
  
      
    ), # closes sidebarPanel
      
      mainPanel( htmlOutput("testHTML"),
                 titlePanel(h3(" Results")),
                 tabsetPanel(type = "tab", id = "tabselected", selected = 1,
                             tabPanel("view", verbatimTextOutput("files"), value = 1), 
                             tabPanel("MSA", verbatimTextOutput("aligns"), value = 2), 
                             tabPanel("Consensus",verbatimTextOutput("cons1"), verbatimTextOutput("cons2"), value = 3), 
                             tabPanel("Plot", plotOutput("ggs1") %>% withSpinner(color="#0dc5c1"), plotOutput("ggs2") %>% withSpinner(color="#c50d11"), value = 4),

                             
                             tabPanel("Info", value = 5,
                                      fluidRow(
                                      conditionalPanel("input.tabs == 'cls'", verbatimTextOutput("text1")),
                                      conditionalPanel("input.tabs == 'com'", verbatimTextOutput("text2")),
                                      conditionalPanel("input.tabs == 'plt'", verbatimTextOutput("text3")),
                                      conditionalPanel("input.tabs == 'readme'", verbatimTextOutput("text4")),
                                      conditionalPanel("input.tabs == 'about'", verbatimTextOutput("text5")),
                                      conditionalPanel("input.tabs == 'structure'", verbatimTextOutput("text6"))
                             
                                      ))
                                                       
    
      
                             
        
       ), #closes tabset panel
    ) #closes main panel   
   ) #closes sidebarlayout
  ) # closes dashborbody  
 ) ## closes dashboarpage
} ### closes user input (UI) section


################################################################################

server <- function(input, output, session) { 
  
  # Captures input data for further processing below in Shiny below ----
  
  readseq <- eventReactive(input$button1, {
    seqdata <- input$file1
    if(is.null(seqdata)){return("upload file")}
    else {
      mySequences <- readAAStringSet(seqdata$datapath)
      return(mySequences)}})
  
  
  
  msaseq <- eventReactive(input$button2,{
    MySequences <- readseq()
    MSAcall <- input$selectMSA
    myAlignment <- msa(MySequences, MSAcall)
    if(is.null(MySequences)){return("No Input")}
    else {
      return(myAlignment)}})
  
  
  down_msa <- reactive({
    msaAlignment <- msaseq()
    sink("msaAlignment.txt", append = TRUE)
    print(msaAlignment, show="complete", halfNrow=NA)
    sink()
  })
  
  conseqmat <- eventReactive(input$button3, {
    myFirstAlignment <- msaseq()
    if(is.null(myFirstAlignment)){return("No Input")}
    else {
      conMat <- consensusMatrix(myFirstAlignment)
      showcons <- msaConsensusSequence(conMat)
      print("Compute consensus sequence using consensus matrix")
      return (showcons)}})
  
  
  
  
  conseqalign <- eventReactive(input$button4, {
    myFirstAlignment <- msaseq()
    maxthres <- input$num1
    minthres <- input$num2
    print("Compute consensus sequence using consensus MSA")
    conalign2 <- printSplitString(msaConsensusSequence(myFirstAlignment))
    
    if(maxthres != 1 || minthres != 1){
      print(paste("Using upperlower type with Threshold =", "(", maxthres, ",", minthres, ")"))
      conalign1 <- printSplitString(msaConsensusSequence(myFirstAlignment, 
                                                         type = "upperlower", thresh = c(maxthres, minthres)))
      return(conalign1)}
    else if((is.null(maxthres)) && (is.null(minthres))){
      return()
    }
  })
  
    
  
  ggmsa_plot1 <- eventReactive(input$button6,{
    start_seq <- input$num3
    end_seq <- input$num4
    color_seq <- input$selectcolor
    myAlignment <- msaseq()
    class(myAlignment) <- "AAMultipleAlignment"
    ggmsa::ggmsa(myAlignment, 
                 start = start_seq, end = end_seq, color = color_seq,
                 
                 char_width = 0.6, font = "helvetical", seq_name = TRUE) })
  
  

  
  
  ggmsa_annotate <- eventReactive(input$button7,{
    gg <- ggmsa_plot1()
    
    if(input$annotate == "Logo"){
      gg + geom_seqlogo() + theme_logo()}
    
    else if(input$annotate == "Consensus"){
      gg + geom_msaBar()}
      
    else if(is.null(input$annotate)){
      return()}
  })
    
  

  
  # Display captured data and processing functions
  
  output$files <- renderPrint({readseq()})
  
  
  output$aligns <- renderPrint({
    if(input$check1 == TRUE){
      showcomp <- 
      return(print(msaseq(), show = "complete"))}
    else{
      return(msaseq())}})
  
  
  output$cons1 <- renderPrint({ conseqmat()})
  
  
  output$cons2 <- renderPrint({ conseqalign()})
  
  
  output$ggs1 <- renderPlot({ggmsa_plot1()})
  
  output$ggs2 <- renderPlot({ggmsa_annotate() })
  
  #download file
  output$downloadData1 <- downloadHandler(
    #specify file name
    filename = "myAlignment.fasta",
    content = function(file){
      down_msa()
      
    }
  )
  

  output$downloadData3 <- downloadHandler(
    filename = function() {paste("color_scheme", input$var1, sep=".") },

    content = function(file) {
      #open the device
      #create the plot
      #close device
      if(input$var1 == "png")
        png(file)
      else
        pdf(file)
      ggmsa_plot1()
      dev.off()

    }
  )

  
  # output$downloadData3 <- downloadHandler(
  #   filename = function() {paste("color_scheme", '.pdf', sep='') },
  #   content = function(file) {
  #     ggsave(file, plot = plotInput(), device = "pdf")
  #   }
  # )
  # 

  
  output$text1 <- renderText({
    paste(
      "The msa function provides a unified interface to the three multiple sequence alignment algorithms in this package:",
      "‘ClustalW’, ‘ClustalOmega’, and ‘MUSCLE’.",
      "",
      "input sequences:  this argument can be a character vector, an object of class XStringSet. ",
      "                  or a single character string and the file must be in FASTA format.",
      "",
      "ClustalW 2.1",
      "     >> cluster:      the default values are Neighbor-Joining.",
      "     >> gapOpening:   the default value for amino acid sequences is 10.0.",
      "     >> gapExtension: the default value for amino acid sequences is 0.2.",
      "",
      "",
      "ClustalOmega 1.2.0",
      "     >> cluster:      the default cluster size is 100.",
      "     >> gapOpening/gapExtension: currently it does not allow to adjust gap penalties.",
      "",
      "",
      "MUSCLE 3.8.31",
      "     >> gapOpening:   the default for amino acid sequences depend on the profile score.",
      "                      *see R documentation for further details*.",
      "     >> gapExtension: the default is 0.",
      sep = "\n" )})
    
  
  output$text2 <- renderText({
    paste(
      "msaConsensusSequence()", 
      "This method computes a consensus sequence from a multiple alignment or a previously computed consensus matrix.",
      "The main purpose of consensus sequences is to get an impression of conservation at individual positions/columns of a multiple alignment.",
      "",
      "     >> input:  an object of class MultipleAlignment or a previously computed consensus matrix",
      "",
      "     >> type:   upperlower... *see R documentation for further details*.",
      "",
      "     >> thresh: a decreasing two-element numeric vector of numbers between 0 and 100 corresponding to the two conservation thresholds.",
      "",
      '     >> ignoreGaps:  a logical (default: FALSE) indicating gaps should be considered.',
      sep = "\n"
    )
  })
  
  output$text3 <- renderText({ 
    paste(
  "This R package (ggmsa, current version: 0.0.6) is avalable via CRAN.",
  "Supports visualizing multiple sequence alignment of DNA and protein sequences using 'ggplot2'.",
  "",
  "  >> msa:   multiple aligned sequence files or objects representing AA sequences.",
  "",
  "  >> start: a numeric vector represent 'Start' position to plot.",
  "",
  "  >> end:   a numeric vector represent 'End' position to plot.",
  "",
  "  >> font:  default is 'helvetical', with 'char_width' = 0.6 ",
  "",
  "  >> color: color schemes for rendering MSA, defaults is 'Chemistry_AA'.",
  "",
  "            Chemistry --> residues are colored according to their side chain chemistry.",
  "            Zappo     --> residues are colored according to their physico-chemical properties.",
  "            LETTER    --> residues are colored according to their letter coding.",
  "            Clustal X  --> residue in the alignment is assigned a colour if the amino acid profile of the alignment at that position meets some minimum criteria specific for the residue type.",
   sep="\n")})
  

  output$text6 <- renderText({
  })
  
  
  output$text4 <- renderText({
    paste(
      "###########################################################################",
      "#....................... Multiple Sequence Alignment .....................#",
      "#..... R Package (current version: 1.29.3) available via Bioconductor.....#",
      "#................ Bioinformatics 31(24):3997--3999 .......................#",
      "#................ DOI: 10.1093/bioinformatics/btv494 .....................#",
      "###########################################################################",
      "",
      
      "Institute of Bioinformatics, Johannes Kepler University Linz.",
      "Authors: U. Bodenhofer, E. Bonatesta, C. Horejs-Kainrath, and S. Hochreiter ", 
      "Contact: msa@bioinf.jku.at",
      sep = "\n")})
  
  
  output$text5 <- renderText({
    paste(
      "########################################",
      "#............. BCB Diploma ............#",
      "#............ DAV - R Course ..........#",
      "#............ MSA Shiny App ...........#",
      "#............. 29/10/2020 .............#",
      "#........... Maha Abdelrhman ..........#",
      "#............ Yasmine Saied ...........#",
      "########################################",
      "",
      "sessionInfo()",
      "R version 4.2.0 (2022-04-22 ucrt)",
      "Platform: x86_64-w64-mingw32/x64 (64-bit)",
      "Running under: Windows 10 x64 (build 19043)",
      
      sep = "\n" )})
  
  
  
  
  ################################# Side-bar Menu ##############################
  
    output$Semi_collapsible_sidebar=renderMenu({
      
        ## add items of menu:
      sidebar <- dashboardSidebar(
        hr(),
      sidebarMenu(id="tabs",
                  
                  menuItem("Codes",  icon = icon("file-text-o"),
                           menuSubItem("MSA type", tabName = "cls", icon = icon("angle-right")),
                           menuSubItem("Consensus", tabName = "com", icon = icon("angle-right")),
                           menuSubItem("ggmsa", tabName = "plt", icon = icon("angle-right")),
                           menuSubItem("geom_logo", tabName = "structure", icon = icon("angle-right"))
                  ),
                  menuItem("ReadMe", tabName = "readme", icon=icon("mortar-board")),
                  menuItem("About", tabName = "about", icon = icon("question")))
                  
            
      
    )## closes sidemenubar
  }) ## closes sidebar output
  
} ### closes user server section

################################################################################

# Run the application 
shinyApp(ui, server)

################################################################################
