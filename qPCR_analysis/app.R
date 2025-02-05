#
#   TO DO
# * Choices: dCT vs exp? SD or SE? consider gaub error propagation form the HK variance? compute SD/SEM in foldchange by antilogging the 2^-dCTs?
# * CHECKBOX TO REMOVE FIELDS
# * - dCT
# * compute geomtric average
# * change the number of columns based on how many samples
# Python code: 


library(shiny)
library(shinythemes)
library(rhandsontable) #editable table to paste data from excel
library(tidyverse)
library(ggsci)
library(shinyjqui) #drag and drop for sample order
library(shinyjs) # hide/show elements

# set as false before deployment
run_local = FALSE

# Define any Python packages needed for the app here:
PYTHON_DEPENDENCIES = c('numpy', "pandas", "xlsxwriter")

if(file.exists("result.xlsx")) {file.remove("result.xlsx")}

empty.df = empty.df <- read.csv("data.csv", header=TRUE)
  
  data.frame(Sample = character(10),
                Gene = character(10),
                CT = double(10))

# Define UI for application that draws a histogram
ui <- fluidPage(
    useShinyjs(),
    theme = shinytheme("flatly"),

    titlePanel("qPCR Analysis"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          h3("Raw Data"),
          helpText("Paste here your row data!", br(),
                    "Don't warry for the space, the size of the table will adapt itself to the content you are pasting"),
          rHandsontableOutput("editable.df")
        ),

        # Show a plot of the generated distribution
        mainPanel(
          textOutput("console"),

          h3("Samples"),
          helpText("Drag and drop to change the order of the samples in the plot"),
          uiOutput("sample_picker"),
          # verbatimTextOutput("print"),
          h3("Loading Control"),
          # Input: Choose House Keeping gene
          uiOutput("selectHK"),
          
          plotOutput("HK.plot"),
          
          radioButtons("value_type", h3("Compute"),
                       choices = list("dCT: Fold changes to HK " = "dCT",
                                      "2^-dCT: Expression Level" = "exp.dCT")),
          
          plotOutput("plot.results"),
          
          h3("Data Preview"),
          tableOutput("data.preview"),
          
          actionButton("run_py", "Generate Excel File"),
          disabled(downloadButton("downloadData", "Download"))
          
        )
    )
)

# initialize to avoid conflicts with python code
#py$df <- data.frame(NULL)
#py$var <- NULL

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # ------------------ App virtualenv setup (Do not edit) ------------------- #
  
  virtualenv_dir = Sys.getenv('VIRTUALENV_NAME')
  python_path = Sys.getenv('PYTHON_PATH')
  
  # Create virtual env and install dependencies
  if (!run_local) {
  reticulate::virtualenv_create(envname = virtualenv_dir, python = python_path)
  reticulate::virtualenv_install(virtualenv_dir, packages = PYTHON_DEPENDENCIES)
  reticulate::use_virtualenv(virtualenv_dir, required = T)
  }
  # --------------------------------------------
  
  # Editable Table for pasting the raw data
  output$editable.df <-renderRHandsontable({rhandsontable(empty.df, allowInvalid = TRUE) %>%
      hot_col(col = "Sample", allowInvalid = TRUE,strict = FALSE) %>%
      hot_col(col = "Gene", allowInvalid = TRUE, strict=FALSE)})
  
  df.unsorted <- reactive({
    hot_to_r(input$editable.df) %>%
      # remove rows with NA values
      filter_all(all_vars(!is.na(.))) %>%
      # drop unuesd levels from previous edits of the table
      droplevels()
    })
      
  
  df <- reactive({
    df.unsorted() %>%
      # for some magical reason, the orderInput automatically appends "_order" 
      # tho the ID of the input so that "input$samples" becames "input$samples_order"
      mutate(Sample = factor(Sample, levels = input$samples_order),
             Gene = factor(Gene))
  })
 
  output$sample_picker <- renderUI({
    orderInput(inputId = 'samples', label='', items = levels(df.unsorted()$Sample), width = "100%")
  })
  
  # output$print <- renderPrint({input$samples_order})
 
  # Reactive input for picking the Loading Contorl Gene
  output$selectHK <- renderUI({
    selectInput("HK","Pick the HK gene", df()$Gene)
  })
  
  
  # plot the CT value of HK genes
  output$HK.plot <- renderPlot({
    df() %>%
      filter(Gene == input$HK) %>%
      ggplot(data = ., aes(x=Sample, y=CT, fill=Sample)) +
        #geom_boxplot() +
        stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.5)+
        geom_jitter(shape=21, position=position_jitter(0), size=5) +
        scale_y_continuous(breaks=c(1:99), expand = c(0.5,0.5))+
        theme_light()+
        theme(
          panel.grid.minor.y = element_blank(),
          legend.position = "none",
          panel.grid.major.x = element_blank())
  })
  
  # compute dCT and 2^-dCT
  computations <- reactive({
    left_join(df() %>% filter(Gene != input$HK),
              df() %>%  #calculate the HK AVG CT for each sample
                filter(Gene == input$HK) %>%
                group_by(Sample) %>%
                summarise(HK_AVG_CT = mean(CT, na.rm = TRUE))) %>%
      mutate(
        dCT = CT - HK_AVG_CT,
        exp.dCT = 2^-(dCT)
      )                  
  })
  
  # summarize
  summary <- reactive({
    computations() %>% 
      # rename the chosen column as key value
      rename(measure_of_int = input$value_type) %>%
      group_by(Sample, Gene) %>%
      # IF else SD or SEM
      summarize(mean = mean(measure_of_int, na.rm = TRUE),
                var = sd(measure_of_int, na.rm=TRUE))
  })
  
  output$plot.results <- renderPlot({
    ggplot(data=summary() %>% filter(Gene != input$HK),
           aes(x=Sample, y=mean, fill=Sample))+
      geom_col(position = "dodge") +
      theme_light() +
      theme(
        legend.position = "none",
        panel.grid.major.x = element_blank())+
      geom_errorbar(aes(ymin = mean - var, ymax = mean + var), position = position_dodge(width = 0.9), width = 0.3) +
      facet_wrap(~ Gene,  scales = "free")
    
  })
  
  # Spread average values
  avg.spread <- reactive({
    summary() %>%
      select(-var) %>%
      spread(key = Sample, value = mean) %>%
      #replace NA with 0 because NA values are not supported by writeXLSX and cause python to carsh
      replace(is.na(.), 0)
  })
  # Spread SD/SEM values
  var.spread <- reactive({
    summary() %>%
      select(-mean) %>%
      spread(key = Sample, value = var) %>%
      replace(is.na(.), 0)
  })
  
  output$data.preview <- renderTable(avg.spread())
  
  
  
  observeEvent(input$run_py, {
    # pass variable from R to Py
    py <- reticulate::import_main()
    py$df <- avg.spread()
    py$var <- var.spread()
    reticulate::py_run_file('py_export_excel.py', local = FALSE, convert = FALSE)
  })
  
  # enable if fileisready is TRUE (delete when app startupp)
  
  output$console <- renderText({
    "test"
  })
  
    
  observe({
    invalidateLater(1000, session)
    # file download is ready
    if(file.exists("result.xlsx")) {enable("downloadData")}
    else {disable("downloadData")}
    })
  
  observeEvent(input$HK, {
    if(input$HK %>% is.null()){
      hide("HK.plot")
    }
  })
  # Downloadable csv of selected dataset ----
  output$downloadData <- downloadHandler(
    filename = "result.xlsx",
    content <- function(file) {
      #wait for the file to be generated
        file.copy("result.xlsx", file)
    }
  )
}

shinyApp(ui = ui, server = server)
