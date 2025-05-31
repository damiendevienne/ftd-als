library(shiny)
library(shinyjs) # For writing messages to console
library(plotly)
library(dplyr)
library(bslib)

df_onset <- readRDS("df_onset.rds")


# Define function to compute the probability of carrying the C9orf72 RE mutation for grandchildren
# equation papa 
pVWH_grandchild <- function(x, pi_t_prime, pi_t) {
  numerator <- (1 - (1 - x) * pi_t) * (1 - (1 - x) * pi_t_prime)
  denominator <- 2 * (2 - (1 - x) * pi_t) - (1 - (1 - x) * pi_t) * (1 - x) * pi_t_prime
  result <- numerator / denominator
  return(result)
}
#equation papa2
pSICK_grandchild <- function(x, pi_t_prime, pi_t) {
  pprime <- (1 - x) * pi_t
  numerator <- (1 - x) * (1 - pi_t_prime) * (1-pprime)
  denominator <- 2 * (2 - pprime) - pi_t_prime * (1 - x) * (1 - pprime) 
  result <- numerator / denominator  
  return(result)
}

      delta <- 30 #age gap
      penetrance_delta <- c(df_onset$Penetrance[(delta+1):length(df_onset$Penetrance)],rep(1, delta))
      penetrances <- cbind(df_onset$Penetrance, penetrance_delta)
      print("ioioioioio")
      print(penetrances)
      print(apply(penetrances, 1, function(w) pVWH_grandchild(w, w[1], w[2]))) #function to compute P(mutation | unaffected) for grandchilds
      print("zzzzzzzzzzz")



# # UI
# ui <- fluidPage(
#   useShinyjs(),  # Initialize shinyjs for console messages
#   # Custom CSS to change font size of slider labels
#   tags$style(HTML("
#     .control-label {
#       font-size: 16px; /* Change font size */
#       font-weight: bold; /* Optional: make the text bold */
#     }
#     .results-panel {
#       background-color: #ffe6e6; /* light red */
#       border: 1px solid #ffcccc;
#       padding: 15px;
#       border-radius: 6px;
#       margin-bottom: 20px;
#     }
#     .result-label {
#       font-weight: bold;
#     }
#   ")),
#   # Change size of the text input box to be smaller
#   tags$script(HTML("
#     $(document).on('shiny:connected', function() {
#       $('#Penetrance').addClass('form-control-sm');
#     });
#   ")),

#   theme = bs_theme(version = 5, bootswatch = "flatly"),
#   fluidRow(
#     column(width = 2),  # Left margin
    
#     column(width = 8,
#       fluidRow(

#           tags$h2("Genetic Counseling (ALS/FTD)", 
#           style = "text-align: center; margin-top: 20px; margin-bottom: 5px;"),
#           tags$h3("Mutation Carrier Probability Estimator of C9orf72 RE", 
#           style = "text-align: center;margin-bottom: 20px;"),
#           tags$p("This application estimates the probability of carrying the C9orf72 RE mutation for unaffected relatives (child or grandchild) of an affected individual. Estimates are based on the formulas presented in de Vienne & de Vienne (In prep), using the age-related penetrance data from Murphy et al. (2017)."),
#           tags$p("If you use this application for publication, please cite:", tags$i("de Vienne D and de Vienne DM. In Prep. Improving Genetic Counseling for C9orf72-ALS/FTD with Age-Based Carrier Estimates"))
#           ),
#       fluidRow(
#         column(width = 6,  # Input panel
#           wellPanel(
#             tags$h4("Options", style = "margin-bottom: 20px;"),

#             tags$div(
#                 style = "display: flex; align-items: baseline;", # Align text and input on the same line
#                 tags$label("Assumed disease penetrance:", style = "margin-right: 10px;"),
#                 numericInput("Penetrance", 
#                             label = NULL,    # Hide label, as we have it above
#                             value = 1, 
#                             min = 0, 
#                             max = 1, 
#                             step = 0.05)
#               ),
# #            sliderInput("Penetrance", "Assumed disease penetrance:", 
#  #                       min = 0, max = 1, value = 1, step = 0.05),
            
#             sliderInput("age", "Patient's current age", value = 40, min = 0, max = 100, step = 1),
            
#             radioButtons("relation", "Family relation to affected individual",
#                          choices = c("Child" = "child", "Grandchild" = "grandchild")),
            
#             conditionalPanel(
#               condition = "input.relation == 'grandchild'",
#               tagList(
#                 sliderInput("parent_age", "Age of the parent (i.e. child of the affected individual)", value = 70, min = 15, max = 100, step = 1),
#                 tags$div(
#                   style = "font-size: 90%; color: #444;",
#                   "Age gap: ",
#                   textOutput("age_gap", inline = TRUE),
#                   " years"
#                 )
#               )
#             )
#           ),
#           br(),
#           wellPanel(
#             class = "results-panel",
#             h4( "Risk Estimates",  style = "margin-bottom: 20px;"),
#             p(
#               span(class = "result-label", "Probability to carry the variant: "),
#               textOutput("carrier_prob_final", inline = TRUE)
#             ),
#             p(
#               span(class = "result-label", "Probability to develop the disease: "),
#               textOutput("sick_prob_final", inline = TRUE)
#             )
#           )
#         ),
        
#         column(width = 6,  # Plot/output panel
#           plotlyOutput("interactive_plot", height = "500px"),
#         )
#       )
#     ),
    
#     column(width = 2)  # Right margin
#   )
# )


# ui <- fluidPage(
#   useShinyjs(),

#   # Styling
#   tags$head(
#     tags$style(HTML("
#       body {
#         background-color: #f7f7f7;
#       }
#       .main-container {
#         background-color: white;
#         max-width: 1200px;
#         margin: 30px auto;
#         padding: 40px;
#         box-shadow: 0 0 15px rgba(0,0,0,0.1);
#         border-radius: 10px;
#         font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
#       }
#       h2, h3 {
#         color: #003366;
#         font-weight: 600;
#       }
#       .results-panel {
#         background-color: #e6f2ff;
#         border-left: 5px solid #3399ff;
#         padding: 15px;
#         border-radius: 6px;
#         margin-top: 20px;
#       }
#       .result-label {
#         font-weight: bold;
#         color: #003366;
#       }
#       .well {
#         background-color: #f9f9f9;
#         border: 1px solid #ccc;
#       }
#     "))
#   ),

#   theme = bs_theme(version = 5, bootswatch = "flatly"),

#   div(class = "main-container",
#     fluidRow(
#       column(
#         width = 12,
#         tags$h2("Genetic Counseling Tool for ALS/FTD", style = "text-align: center;"),
#         tags$h4("C9orf72 RE – Carrier and Risk Estimator", style = "text-align: center; margin-bottom: 25px;"),
#         tags$p("This tool estimates the probability of carrying the C9orf72 repeat expansion and developing ALS/FTD, for unaffected children or grandchildren of affected individuals."),
#         tags$p("Estimates are based on penetrance data from Murphy et al. (2017) and probabilistic models developed by de Vienne & de Vienne (in prep)."),
#         tags$p("If used in publications, please cite: ", tags$i("de Vienne D and de Vienne DM. In Prep. Improving Genetic Counseling for C9orf72-ALS/FTD with Age-Based Carrier Estimates."))
#       )
#     ),

#     br(),

#     fluidRow(
#       column(width = 6,
#         wellPanel(
#           tags$h5("User Inputs"),
#           sliderInput("age", "Age of the person being evaluated:", value = 40, min = 0, max = 100, step = 1),
#           numericInput("Penetrance", "Assumed lifetime penetrance of the mutation:", value = 1, min = 0, max = 1, step = 0.05),
#           radioButtons("relation", "Relationship to the affected individual:",
#             choices = c("Child" = "child", "Grandchild" = "grandchild")
#           ),
#           conditionalPanel(
#             condition = "input.relation == 'grandchild'",
#             tagList(
#               sliderInput("parent_age", "Age of the intermediate parent:", value = 70, min = 15, max = 100, step = 1),
#               tags$div(
#                 style = "font-size: 90%; color: #555;",
#                 "Age gap (parent - person being evaluated): ",
#                 textOutput("age_gap", inline = TRUE),
#                 " years"
#               )
#             )
#           )
#         ),

#         div(class = "results-panel",
#           h5("Estimated Probabilities"),
#           p(span(class = "result-label", "Probability of carrying the mutation: "), textOutput("carrier_prob_final", inline = TRUE)),
#           p(span(class = "result-label", "Probability of developing the disease: "), textOutput("sick_prob_final", inline = TRUE))
#         )
#       ),

#       column(width = 6,
#         h5("Age-dependent probability curves"),
#         plotlyOutput("interactive_plot", height = "500px")
#       )
#     )
#   )
# )

ui <- fluidPage(
  useShinyjs(),

  # Custom CSS
  tags$style(HTML("
    body {
      background-color: #f8f9fa;
      font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
      color: #333;
    }
    .control-label {
      font-size: 16px;
      font-weight: bold;
    }
    .results-panel {
      background-color: #f1f1f1;
      border: 1px solid #ddd;
      padding: 15px;
      border-radius: 6px;
      margin-top: 20px;
    }
    .result-label {
      font-weight: bold;
    }
    .well {
      background-color: #ffffff;
      border: 1px solid #ccc;
    }
  ")),

  # Theme
  theme = bs_theme(version = 5, bootswatch = "flatly"),

  # Main layout
  fluidRow(
    column(width = 2),  # Left margin
    
    column(width = 8,
      # Title and description
      fluidRow(
        column(12,
          tags$h2("Genetic Counseling for ALS/FTD", style = "text-align: center; margin-top: 20px;"),
          tags$h4("Carrier Probability Estimation for C9orf72 Repeat Expansion", style = "text-align: center; margin-bottom: 20px;"),
          tags$p("This application estimates the probability of carrying the C9orf72 repeat expansion mutation for an unaffected first- or second-degree relative of an affected individual."),
          tags$p("Estimates are based on age-specific penetrance data (Murphy et al., 2017) and formulae from de Vienne & de Vienne (in prep)."),
          tags$p("For publication use, please cite:", 
                 tags$i("de Vienne D & de Vienne DM. In prep. Improving Genetic Counseling for C9orf72-ALS/FTD with Age-Based Carrier Estimates."))
        )
      ),

      # Inputs and Outputs
      fluidRow(
        column(width = 6,
          wellPanel(
            tags$h4("Input Parameters", style = "margin-bottom: 15px;"),

            tags$div(style = "margin-bottom: 10px;",
              tags$label("Assumed disease penetrance:", style = "margin-right: 10px;"),
              numericInput("Penetrance", label = NULL, value = 1, min = 0, max = 1, step = 0.05)
            ),

            sliderInput("age", "Current age of the patient:", value = 40, min = 0, max = 100, step = 1),

            radioButtons("relation", "Family relationship to the affected individual:",
                         choices = c("Child" = "child", "Grandchild" = "grandchild")),

            conditionalPanel(
              condition = "input.relation == 'grandchild'",
              sliderInput("parent_age", "Age of the parent (child of the affected individual):", value = 70, min = 15, max = 100, step = 1),
              tags$div(
                style = "font-size: 90%; color: #444;",
                "Age gap: ",
                textOutput("age_gap", inline = TRUE),
                " years"
              )
            )
          ),

          wellPanel(class = "results-panel",
            h4("Estimated Risks", style = "margin-bottom: 15px;"),
            p(span(class = "result-label", "Probability to carry the mutation: "),
              textOutput("carrier_prob_final", inline = TRUE)),
            p(span(class = "result-label", "Probability to develop the disease: "),
              textOutput("sick_prob_final", inline = TRUE))
          )
        ),

        column(width = 6,
          plotlyOutput("interactive_plot", height = "500px")
        )
      )
    ),

    column(width = 2)  # Right margin
  ),

  # Footer
  tags$footer(
    style = "
      background-color: #f0f0f0;
      padding: 15px 20px;
      margin-top: 30px;
      text-align: center;
      font-size: 90%;
      color: #555;
      border-top: 1px solid #ddd;
    ",
    HTML(paste0(
      'Source code on <a href="https://github.com/damiendevienne/ftd-als" target="_blank">GitHub</a> &nbsp;|&nbsp; ',
      'Licensed under the <a href="https://opensource.org/licenses/MIT" target="_blank">MIT License</a> &nbsp;|&nbsp; ',
      'Contact: <a href="mailto:dominique.de-vienne@inrae.fr">dominique.de-vienne@inrae.fr</a>'
    ))
  )
)


# Server
server <- function(input, output, session) {

  # Reactive expression to compute parent age
  age_gap <- reactive({
    req(input$relation == "grandchild", input$age, input$parent_age)
    input$parent_age - input$age
  })

  output$age_gap <- renderText({
    req(input$relation == "grandchild", input$age, input$parent_age)
    input$parent_age - input$age
  })

  # Adjusted dataframe with user-defined penetrance
  # Inside this df we can compute everything we need. 

  adjusted_df <- reactive({ #this reacts to changes in input$Penetrance ONLY for now
    x <- 1-input$Penetrance #proportion of non-penentrant carriers in the population
    print(x)
    # df <- df_onset
    df <- data.frame(Age=df_onset$Age, 
      Penetrance = df_onset$Penetrance
      )
    if (input$relation=="child") {
      df$PvariantWhenHealthy <- (1-(1-x)*df_onset$Penetrance)/(2-(1-x)*df_onset$Penetrance)
      df$PsickInFuture <- ((1-x)*(1-df_onset$Penetrance))/(2-(1-x)*df_onset$Penetrance)
    }
    if (input$relation=="grandchild") {
      #prepare data to compute probabilities
      delta <- input$parent_age - input$age #age gap
      penetrance_delta <- c(df_onset$Penetrance[(delta+1):length(df_onset$Penetrance)],rep(1, delta))
      penetrances <- cbind(df_onset$Penetrance, penetrance_delta)
      df$PvariantWhenHealthy <- apply(penetrances, 1, function(w,z) pVWH_grandchild(z, w[1], w[2]), z=x) #function to compute P(mutation | unaffected) for grandchilds
#      df$PvariantWhenHealthy <- (1-(1-x)*df_onset$Penetrance)/(2-(1-x)*df_onset$Penetrance)
      df$PsickInFuture <- apply(penetrances, 1, function(w,z) pSICK_grandchild(z, w[1], w[2]), z=x)
    }
    df
  })
  
  # Get Carrier Prob and sick prob
  get_carrier_prob <- reactive({
      df <- adjusted_df()
      age <- input$age
      df$PvariantWhenHealthy[df$Age == age]
  })

  get_sick_prob <- reactive({
      df <- adjusted_df()
      age <- input$age
      df$PsickInFuture[df$Age == age]
  })

  output$carrier_prob_final <- renderText({
    round(get_carrier_prob(),2)
  })

  output$sick_prob_final <- renderText({
    round(get_sick_prob(), 2)
  })
  
  output$interactive_plot <- renderPlotly({
    df <- adjusted_df()
    age_input <- input$age
    prob_input1 <- get_carrier_prob()
    prob_input2 <- get_sick_prob()
    
    plt <- plot_ly(
      data = df, x = ~Age, y = ~PvariantWhenHealthy, type = 'scatter', mode = 'lines',
      line = list(color = 'steelblue'),
      hoverinfo = 'x+y',
      name = 'P(mutation | unaffected)'
    ) %>%
      add_trace( 
      y = ~PsickInFuture,
      line = list(color = 'firebrick'),
      name = 'P(disease | unaffected)',
      hoverinfo = 'x+y',
      mode = 'lines'
    ) %>%
      layout(
        title = paste("Probability estimates for a ", input$relation ,sep=""),
        xaxis = list(title = "Age"),
        yaxis = list(title = "Estimated probability"),
        shapes = list(
          list(
            type = "line",
            x0 = age_input, x1 = age_input,
            y0 = 0, y1 = 1,
            line = list(dash = "dot", color = "grey")
          ),
          list(
            type = "line",
            x0 = 0, x1 = 100,
            y0 = prob_input1, y1 = prob_input1,
            line = list(dash = "dot", color = "steelblue")
          ),
          list(
            type = "line",
            x0 = 0, x1 = 100,
            y0 = prob_input2, y1 = prob_input2,
            line = list(dash = "dot", color = "firebrick")
          )
        ), 
        legend = list(
        x = 1,          # horizontal position (0 = left, 1 = right)
        y = 1,          # vertical position (0 = bottom, 1 = top)
        xanchor = "right",  # anchor point of legend box
        yanchor = "top",
        bgcolor = "rgba(255,255,255,0.5)",  # optional: semi-transparent white background
        bordercolor = "gray",
        borderwidth = 1
        )
      ) %>%
      config(
        displayModeBar = FALSE,  # Hides the mode bar (zoom, pan, etc.)
        editable = FALSE        # Disables editing of plot
#        staticPlot = TRUE        # Makes the plot static (no interaction)
      )
    
    
    plt
  })
}

shinyApp(ui, server)
