
load("ctdf.rda")
load("qtldf.rda")

qtsel <- list(proteins = sort(unique(qtldf$protein)), outcomes = sort(unique(qtldf$outcome)))
ctsel <- list(exposures = sort(unique(ctdf$exposure)), outcomes = sort(unique(ctdf$outcome)))

navbarPage("MR Browser",

  # App title ----
  tabPanel("Complex trait exposures",

           # Sidebar layout with input and output definitions ----
           sidebarLayout(

             # Sidebar panel for inputs ----
             sidebarPanel(

               # Input: Slider for the number of bins ----
               selectizeInput("exposurect", "Exposure", multiple = TRUE,
                              choices = ctsel$exposures),


               selectizeInput("outcomect", "Outcome",
                              choices = ctsel$outcomes)

             ),

             # Main panel for displaying outputs ----
             mainPanel(

               helpText("Choose exposure(s) to show plots. Scroll down for table"),
               # Output: combined snps ----
               h3("Exposure metadata"),
               dataTableOutput(outputId = "metaTab"),
               h3("Combined SNPs"),
               plotOutput(outputId = "mrPlotCT"),
               downloadButton("downloadPlotCT", "Download plot"),
               textInput("wp1", "width (cm)", value = "12", width = "80px"),
               textInput("hp1", "height (cm)", value = "8", width = "80px"),
               h3("individual SNPs"),
               plotOutput(outputId = "mrPlotCT2", height = "800px"),
               downloadButton("downloadPlotCT2", "Download plot"),
               textInput("wp2", "width (cm)", value = "12", width = "80px"),
               textInput("hp2", "height (cm)", value = "8", width = "80px"),
               h3("Numeric results"),
               downloadButton("downloadTabCT", "Download data (xlsx)"),
               dataTableOutput("snpTabCT")

             )
           )
           ),

  tabPanel("Protein QTL exposures",

           sidebarLayout(

             # Sidebar panel for inputs ----
             sidebarPanel(

               # Input: Slider for the number of bins ----
               selectizeInput("exposureqtl", "Protein exposure", multiple = TRUE,
                              choices = qtsel$protein),

               selectInput("cispan", "Test set", c("cis", "pan"), selected = "cis",
                           selectize = FALSE),


               selectizeInput("outcomeqtl", "Outcome",
                              choices = qtsel$outcomes)

             ),

             # Main panel for displaying outputs ----
             mainPanel(

               # Output: Histogram ----
               h3("Combined SNPs"),
               plotOutput(outputId = "mrPlotQTL"),
               downloadButton("downloadPlotQTL", "Download plot"),
               textInput("wp3", "width (cm)", value = "12", width = "80px"),
               textInput("hp3", "height (cm)", value = "8", width = "80px"),
               h3("Individual SNPs"),
               plotOutput(outputId = "mrPlotQTL2", height = "800px"),
               downloadButton("downloadPlotQTL2", "Download plot"),
               textInput("wp4", "width (cm)", value = "12", width = "80px"),
               textInput("hp4", "height (cm)", value = "8", width = "80px"),
               h3("Numeric results for MR analysis"),
               downloadButton("downloadqtlTabCT", "Download data (xlsx)"),
               dataTableOutput(outputId = "qtlTabCT"),
               h3("Numeric results for IVs (SNP -> protein association)"),
               downloadButton("downloadivTabQT", "Download data (xlsx)"),
               dataTableOutput(outputId = "ivTabQT")

             )
           )
           )


)
