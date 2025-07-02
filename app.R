# ==============
# Adaptation of Making Markov Models Shiny 
# Original authors: Robert Smith & Paul Schneider University of Sheffield
# ==============

## app.R ##

# install.packages("shiny") # necessary if you don't already have the function 'shiny' installed.

# we need the function shiny installed, this loads it from the library.
library(shiny)
library(shinyBS)
library(rmarkdown)
library(bookdown)
library(here)
library(shinyWidgets)
library(shinyjs)
library(shinycssloaders)
library(DT)
library(priceR)
library(scales)
library(mgcv)
library(mgcViz)
library(voi)
library(gt)
library(ggplot2)
library(fontawesome)
library(rsconnect)


# source the wrapper function from directory
source("wrapper.R")

# Define informative discrepancy formulas
ID_formula <- list("\\(\\delta^{U}_{Healthy}\\)","\\(\\delta^{C}_{Healthy}\\)", "\\(\\delta^{U}_{Sick, Treatment}\\)", "\\(\\delta^{C}_{Sick, Treatment}\\)")
ID_spec <- list("$$\\delta^{U}_{Healthy} ~\\sim Beta(\\alpha, \\beta)$$", "$$\\delta^{C}_{Healthy} ~\\sim Gamma(shape, scale)$$", 
                "$$\\delta^{U}_{Sick, Treatment} ~\\sim N(\\mu, \\sigma)$$", "$$\\delta^{C}_{Sick, Treatment} ~\\sim Gamma(shape, scale)$$")
ID_list <- setNames(ID_spec, ID_formula)

#================================================================
#                   Create User Interface
#================================================================

ui <- fluidPage(    # create user interface using fluidpage function
  
  
  tags$head((tags$style(HTML('.modal-lg{
                             width: 95%;
                             }
                             ')))),
  
  titlePanel("Application of dicrepancy terms to Sick Model"),   # title of app
  
  # SIDEBAR
  sidebarLayout(    # indicates layout is going to be a sidebar-layout
    
    sidebarPanel(# open sidebar panel
      useShinyjs(),
      withMathJax(),
      numericInput(inputId = "SI_n_sim",      # id of input, used in server
                   label = "PSA runs",        # label next to numeric input
                   value = 100,              # initial value
                   min = 50,                   # minimum value allowed
                   max = 1000),             # maximum value allowed
      bsPopover(id = "SI_n_sim", 
                title = "", 
                content = paste0("Minimum value of 50. Maximum value of 1,000"), 
                placement = "right", 
                trigger = "hover", 
                options = list(container = "body")),
      autonumericInput(inputId = "WTP", 
                  label = "Willingness to pay threshold",
                  align = "left",
                  value = 20000,
                  step = 1000,
                  minimumValue = 0, 
                  maximumValue = 100000,
                  decimalPlaces = 0,
                  currencySymbol = "£",
                  currencySymbolPlacement = "p",
                  digitGroupSeparator = ","),
                  bsPopover(id = "WTP",
                            title = "",
                            content = paste0("Maximum value of £100,000"),
                            placement = "right",
                            trigger = "hover", 
                            options = list(container = "body")),
      pickerInput(inputId = "switch_model",
                  label = "Model",
                  choices = c("Sick Model", "Sick-Sicker Model"), 
                  selected = "Sick Model"
        ),
      conditionalPanel(
        condition = "input.switch_model == 'Sick Model'",
      pickerInput(inputId = "switch_discrep",
                  label = "Include discrepancy terms",
                  choices = c("Yes", "No"), 
                  selected = "No")), 
      conditionalPanel(
        condition = "input.switch_model == 'Sick Model' && input.switch_discrep == 'Yes'",
      pickerInput(inputId = "switch_spec",
                  label = "Discrepancy term specification",
                  choices = c("Non-informative", "Informative"),
                  selected = "Non-informative")
      ),
      conditionalPanel(
        condition = "input.switch_model == 'Sick Model' && input.switch_discrep == 'Yes' && input.switch_spec == 'Non-informative'",
         # title = "Non-informative discrepancy specification", 
          sliderInput(inputId = "NID_pwcorr", 
                      label = "Pairwise correlation value between discrepancy terms",
                      value = 0,
                      step = 0.1,
                      min = 0,
                      max = 0.9),
          selectInput(inputId = "NID_sd", 
                       label = "Standard deviation (expressed as % of the mean parameter value)",
                       choices = c(0.025, 0.05, 0.1),
                       selected = 0.05)
      ),
      
      conditionalPanel(
        condition = "input.switch_model == 'Sick Model' && input.switch_discrep == 'Yes' && input.switch_spec == 'Informative'",
        selectizeInput(inputId = "ID_utility_select", 
                     label = "Select informative discrepancy term",
                    choices = NULL,
                    options = list(create = TRUE)
                     ),
       uiOutput("ID"),
       
     # hidden(tags$div(id = "ID1", uiOutput("ID1"))),
      
      conditionalPanel(
        condition = paste0("input.ID_utility_select.indexOf('","delta^{U}_{Healthy}","') != -1 "),
        numericInput(inputId = "ID1_param1", label = "\\(\\alpha\\)", value = 3, min = 1, max = 10),
        numericInput(inputId = "ID1_param2", label = "\\(\\beta\\)", value = 15, min = 10, max = 20)),
     conditionalPanel(
        condition = paste0("input.ID_utility_select.indexOf('","delta^{C}_{Healthy}","') != -1 "),
        numericInput(inputId = "ID2_param1", label = "Shape", value = 100, min = 1, max = 500),
        numericInput(inputId = "ID2_param2", label = "Scale", value = 15, min = 1, max = 100)),
     conditionalPanel(
       condition = paste0("input.ID_utility_select.indexOf('","delta^{U}_{Sick, Treatment}","') != -1 "),
       numericInput(inputId = "ID3_param1", label = "\\(\\mu\\)", value = -0.01, min = -0.2, max = 0.2),
       numericInput(inputId = "ID3_param2", label = "\\(\\sigma\\)", value = 0.01, min = 0, max = 0.5)),
     conditionalPanel(
     condition = paste0("input.ID_utility_select.indexOf('","delta^{C}_{Sick, Treatment}","') != -1 "),
       numericInput(inputId = "ID4_param1", label = "Shape", value = 100, min = 1, max = 500),
       numericInput(inputId = "ID4_param2", label = "Scale", value = 20, min = 1, max = 100))
     
      ),

      actionButton(inputId = "run_model",     # id of action button, used in server
                   label   = "Run model"),  # action button label (on button)
      downloadButton("report", "Download report"),
      tags$head(tags$style(HTML(".shiny-notification {
        position: static;
        font-size: 15px;
        }"
        )))
      
    ),  # close sidebarPanel
    

    mainPanel(                                # open main panel
      
      h3("Model results"),                    # heading (results table)                
      
      tableOutput(outputId = "SO_icer_table"),   # tableOutput id = icer_table, from server
      
      conditionalPanel(
        condition = "input.switch_model == 'Sick Model' && input.switch_discrep == 'Yes' && input.switch_spec == 'Non-informative'",

              hidden(tags$div(id = "NID_output",
              h3("EVPPI results"),
              tableOutput(outputId = "evppi_table"),
              htmlOutput(outputId = "evppi_warning"),
              actionButton("diag", "View model diagnostic plots"),
              bsModal(id="diagnostics", title = "Diagnostic Plots", trigger = "diag", size = "large",
                      withSpinner(fluidPage(tabsetPanel(id="diagplots",
                                                        tabPanel(title = "NMB Treatment", plotOutput("gam_check_trt")),
                                                        tabPanel(title = "NMB No Treatment", plotOutput("gam_check_notrt")))))
                      )))
        
        ),
      
      h3("Cost-effectiveness plane"),         # heading (Cost effectiveness plane)
      
      plotOutput(outputId = "CE_plane")       # plotOutput id = CE_plane, from server
      
    ) # close mainpanel    
    
  ) # close sidebarlayout
  
) # close UI fluidpage

#================================================================
#                     Create Server Function
#================================================================

server <- function(input, output, session){   # server = function with two inputs
  
  
  observeEvent(input$run_model,       # when action button pressed ...
               ignoreNULL = F, {
                 
                 
                 # Run model wrapper function with the Shiny inputs and store as data-frame 
                 df_model_res = f_wrapper(n_sim = input$SI_n_sim,
                                          WTP = input$WTP,
                                          model = input$switch_model,
                                          include_discrepancy = input$switch_discrep,
                                          discrepancy_spec = input$switch_spec,
                                          NID_pwcorr = input$NID_pwcorr,
                                          NID_sd = input$NID_sd,
                                          ID1_1 = input$ID1_param1,
                                          ID1_2 = input$ID1_param2, 
                                          ID2_1 = input$ID2_param1,
                                          ID2_2 = input$ID2_param2, 
                                          ID3_1 = input$ID3_param1,
                                          ID3_2 = input$ID3_param2,
                                          ID4_1 = input$ID4_param1,
                                          ID4_2 = input$ID4_param2
                                          )
                 # Populate informative discrepancy term selectize input

                 updateSelectizeInput(session, "ID_utility_select", choices = ID_list, server = TRUE)
                 
                 output$ID <- renderUI({
                   
                   withMathJax(paste0(input$ID_utility_select)) 
                   
                 })
                 
                  #--- CREATE COST EFFECTIVENESS TABLE ---#
                 output$SO_icer_table <- renderTable({ # this continuously updates table
                   
                         df_res_table <- data.frame( # create dataframe
                           Option =  c("Treatment","No Treatment"), 
                           QALYs  =  c(mean(df_model_res$QALY_Trt),mean(df_model_res$QALY_NoTrt)),
                           Costs  =  c(mean(df_model_res$Cost_Trt),mean(df_model_res$Cost_NoTrt)),
                           "Incremental QALYs" = c(mean(df_model_res$QALY_Trt) - mean(df_model_res$QALY_NoTrt),NA),
                           "Incremental Costs" = c(mean(df_model_res$Cost_Trt) - mean(df_model_res$Cost_NoTrt),NA),
                           ICER = c((mean(df_model_res$Cost_Trt) - mean(df_model_res$Cost_NoTrt))/(mean(df_model_res$QALY_Trt) - mean(df_model_res$QALY_NoTrt)), NA), 
                           NMB = c(mean(df_model_res$NMB_Trt), mean(df_model_res$NMB_NoTrt)), 
                           Prob.CE = c(mean(df_model_res$pCE),NA)
                         )
                       # format dataframe to add currency and rounding
                       df_res_table[,c(2,4)] <- round(df_res_table[,c(2,4)],digits = 3)
                       df_res_table[,c(3,5,6,7)] <- format_currency(df_res_table[,c(3,5,6,7)], "£", 0)
                       df_res_table[,8] <- percent(df_res_table[,8])
                       # Remove dots from column headers
                       names(df_res_table) <- gsub("\\.", " ", names(df_res_table))
                       #print the dataframe
                       df_res_table
                 }) # table plot end.

                 # Estimate EVPPI for non-informative discrepancy terms using PSA output and GAMs

                 gam_trt <- gam(NMB_Trt ~ te(NID_utility_trt, NID_cost_trt) + te(NID_utility, NID_cost), data = df_model_res)

                 gam_notrt <- gam(NMB_NoTrt ~ te(NID_utility_trt, NID_cost_trt) +  te(NID_utility, NID_cost), data = df_model_res)

                 # # Diagnostic test to check if specification of number of knots (k) is appropriate
                output$gam_check_trt <- renderPlot({
                  check.gamViz(getViz(gam_trt))
                })
                 
                output$gam_check_notrt <- renderPlot({
                 check.gamViz(getViz(gam_notrt))
                }) 

                 # # Calculate residuals
                 resid_trt <- residuals(gam_trt, type = "deviance")
                 resid_notrt <- residuals(gam_notrt, type = "deviance")

                 # Calculate fitted values
                 g.hat_trt <- gam_trt$fitted
                 g.hat_notrt <- gam_notrt$fitted

                 # # Calculate EVPPI for discrepancy terms & EVPI
                 evppi <- mean(pmax(g.hat_trt, g.hat_notrt)) - max(mean(g.hat_trt), mean(g.hat_notrt))
                 evpi <- evpi(data.frame(df_model_res$NMB_Trt, df_model_res$NMB_NoTrt))
                 evppi_index <- evppi/evpi
                 
                 #--- Create EVPPI Table (for non-informative discrepancy terms) ---#
                 output$evppi_table <- renderTable({
                   df_res_evppi <- data.frame(
                     EVPPI = evppi,
                     "Overall EVPI" = evpi,
                     "EVPPI Index" = evppi_index)
                   
                 # Remove dots from column headers
                   names(df_res_evppi) <- gsub("\\.", " ", names(df_res_evppi))
                 # Format data table to add relevant symbols and rounding
                   df_res_evppi[,c(1,2)] <- format_currency(df_res_evppi[,c(1,2)], "£", 0)
                   df_res_evppi[,3] <- percent(df_res_evppi[,3])
                   # Print data frame
                   df_res_evppi}
                 )
                 
                 warning_text <-
                   ifelse(
                     evppi_index > 1, HTML(paste0('<b>', 'EVPPI Index exceeds 100%', '</b>', '<br>',
                                                  'Check variance parameter, GAM specification and a sufficient number of PSA iterations have been specified.')),
                     ''
                   )

                 output$evppi_warning <- renderText({
                   warning_text
                 })

                 if(input$switch_discrep == "Yes" && input$switch_spec == "Non-informative") {
                   show("NID_output")
                 }
                   else{
                     hide("NID_output")
                   }
                 #---  CREATE COST EFFECTIVENESS PLANE ---#
                 
                 # calculate incremental costs and qalys from results dataframe
                 df_model_res$inc_C <- df_model_res$Cost_Trt - df_model_res$Cost_NoTrt
                 df_model_res$inc_Q <- df_model_res$QALY_Trt - df_model_res$QALY_NoTrt

                 output$CE_plane <- renderPlot({
                   ggplot(df_model_res, aes(x=inc_Q, y=inc_C)) + geom_point() + labs(x = "Incremental QALYs", y = "Incremental costs") +
                     scale_x_continuous(expand = expansion(mult = c(0, 0))) +
                                          scale_y_continuous(expand = expansion(mult = c(0, 0)), labels = label_currency(prefix = "£", big.mark = ",")) +
                                          coord_cartesian(xlim = c(0, 1), ylim = c(0,24000)) +
                                          geom_abline(intercept = 0, slope = input$WTP, colour = "darkred") + theme_classic() +
                                          theme(axis.title = element_text(size = 16), axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14, vjust = 0.25))
                   })
                 
               }) # Observe Event End
  

  # Generate report
  output$report <- downloadHandler(
    filename = function(){
      "report.pdf"},
    content = function(file){
      src <- normalizePath("report.Rmd")
      
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, "report.Rmd", overwrite = TRUE)
      out <- rmarkdown::render("report.Rmd", 
                               params = list(
                                 n = input$slider,
                                 t = input$number
                               ), "pdf_document")
      file.rename(out, file)
    }

    
  )
  
} # Server end


## ----- run app------

shinyApp(ui, server)

