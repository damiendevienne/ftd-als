require(shiny)
library(shinyjs) # For writing messages to console
library(plotly) # For interactive plots
#library(dplyr)
library(bslib)

#df_onset <- readRDS("df_onset.rds")
df_onset <- read.csv("df_onset.csv", stringsAsFactors = FALSE)
###########################
## FUNCTIONS (JULY 2025) ##
###########################

### Probability of carrying the variant when healthy (pVWH)
# Child
pVWH_child <- function(x, pi_t) {
  pprime_t <- (1 - x) * pi_t
  result <- (1 - pprime_t) / (2 - pprime_t)
  result
}

# Grandchild
pVWH_grandchild <- function(x, pi_t2, pi_t1) {
  pprime_t1 <- (1 - x) * pi_t1
  pprime_t2 <- (1 - x) * pi_t2
  numerator <- (1 - pprime_t2)
  denominator <- 2 * ((2 - pprime_t1) / (1 - pprime_t1)) - pprime_t2
  result <- numerator / denominator
  result
}

### Probability of being affected at age t+n when unaffected at age t (pSICKn)
# Child
pSICKn_child <- function(x, pi_t, pi_t_n) { #pi_t_n is the penetrance at age t+n
  pprime_t <- (1 - x) * pi_t
  pprime_t_n <- (1 - x) * pi_t_n
  result <- (pprime_t_n -pprime_t)/(2-pprime_t)
  result
}

# Grandchild
pSICKn_grandchild <- function(x, pi_t2,pi_t2_n, pi_t1) { #pi_t2_n is the penetrance at age t2+n
  pprime_t1 <- (1 - x) * pi_t1
  pprime_t2 <- (1 - x) * pi_t2
  pprime_t2_n <- (1 - x) * pi_t2_n
  numerator <- (pprime_t2_n - pprime_t2)
  denominator <- 2 * ((2 - pprime_t1) / (1 - pprime_t1)) - pprime_t2
  result <- numerator / denominator
  result
}

# # Ploting function
  make_plot <- function(df, age_input, prob_line, y_col, line_color, plot_title) {
    plot_ly(data = df, x = ~Age, y = as.formula(paste0("~", y_col)), type = 'scatter', mode = 'lines',
            line = list(color = line_color),
            hoverinfo = 'x+y',
            name = plot_title) %>%
      layout(
#        title = list(text = plot_title, font=list(size=16)),
        margin = list(l = 60, r = 10, t = 10, b = 60),
        xaxis = list(title = "Age of the consultand"),
        yaxis = list(title = "Estimated probability"),
        # plot_bgcolor = "rgba(0,0,0,0)",   # transparent plot area
        # paper_bgcolor = "rgba(0,0,0,0)",  # transparent outer area
        shapes = list(
          list(
            type = "line",
            x0 = age_input, x1 = age_input,
            y0 = 0, y1 = 1,
            line = list(dash = "dot", color = line_color)
          ),
          list(
            type = "line",
            x0 = 0, x1 = 100,
            y0 = prob_line, y1 = prob_line,
            line = list(dash = "dot", color = line_color)
          )
        ),
        legend = list(
          x = 1,
          y = 1,
          xanchor = "right",
          yanchor = "top",
          bgcolor = "rgba(255,255,255,0.5)",
          bordercolor = "gray",
          borderwidth = 1
        )
      ) %>%
      config(
        displayModeBar = FALSE,
        editable = FALSE
      )
  }

####################
## USER INTERFACE ##
####################

library(shiny)
library(shinyjs)
library(bslib)

ui <- fluidPage(
  useShinyjs(),

  # Bootstrap Theme
  theme = bs_theme(version = 5, bootswatch = "yeti"),  # Change here for different looks

  # Header

  tags$header(
      class = "bg-primary text-white py-4 mb-4",
      div(class = "container",
        tags$h2("Genetic Counseling for ALS/FTD", class = "mb-2"),
        tags$h4("Carrier Probability Estimation for C9orf72 Repeat Expansion", class = "mb-3"),
        div(
          class = "small",
          HTML(paste0(
            'Source code on <a href="https://github.com/damiendevienne/ftd-als" target="_blank" class="text-white text-decoration-underline">GitHub</a> &nbsp;|&nbsp; ',
            'Licensed under the <a href="https://opensource.org/licenses/MIT" target="_blank" class="text-white text-decoration-underline">MIT License</a>'
          ))
        )
      )
    ),

  # Main layout
  fluidRow(
    column(width = 1),  # Left margin

    column(width = 10,
      # Title and description
      fluidRow(
        column(12,
          tags$p("This application estimates the probability of carrying the C9orf72 repeat expansion mutation for an unaffected first- or second-degree relative of an affected individual. Estimates are based on age-specific penetrance data (Murphy et al., 2017) and methods by de Vienne & de Vienne (in prep)."),
          tags$p("For publication use, please cite:",
                 tags$i("de Vienne D & de Vienne DM. In prep. Improving Genetic Counseling for C9orf72-ALS/FTD with Age-Based Risk Estimates."))
        )
      ),

      # Inputs and Outputs
      fluidRow(
        column(width = 3,
          wellPanel(
            tags$h4("Input Parameters", class = "mb-3"),

            radioButtons("relation", "Family relationship of the consultand to the affected individual:",
                         choices = c("Child" = "child", "Grandchild" = "grandchild")),

            numericInput("Penetrance", label = "Assumed disease penetrance:", value = 1, min = 0, max = 1, step = 0.05),

            sliderInput("age", "Current age of the patient:", value = 40, min = 0, max = 100, step = 1),

            conditionalPanel(
              condition = "input.relation == 'grandchild'",
              sliderInput("parent_age", "Age of the parent (child of the affected individual):", value = 70, min = 15, max = 100, step = 1),
              tags$div(
                class = "small text-muted",
                "Age gap: ",
                textOutput("age_gap", inline = TRUE),
                " years"
              )
            ), 

            numericInput("n_years", 
              label = "Time frame for disease risk estimation (in years):", 
              value = 10, 
              min = 1, 
              max = 100, 
              step = 1)
          )
        ),

        column(width = 9,
          wellPanel(
            tags$h4("Risk estimates visualization"),

            fluidRow(
              column(width = 6,
                tags$h5("Probability of carrying the C9orf72 RE mutation"),
                plotlyOutput("interactive_plot", height = "500px")
              ),
              column(width = 6,
                uiOutput("dynamic_title_nyears"),
                plotlyOutput("interactive_plot2", height = "500px")
              )
            )
          ),

          wellPanel(
            tags$h4("Risk estimates summary", class = "mb-3"),
            tags$p(tags$strong("Probability to carry the mutation: "),
              textOutput("carrier_prob_final", inline = TRUE)),
            tags$p(tags$strong("Probability to develop the disease: "),
              textOutput("sick_prob_final", inline = TRUE))
          )
        )
      )
    ),

    column(width = 1)  # Right margin
  ),

  # Footer
  tags$footer(
    class = "text-center py-3 border-top mt-4",
    HTML(paste0(
      '<img src="logo_lbbe.svg" alt="LBBE logo" height="30" style="vertical-align: middle; margin-right: 10px;">',
      '&nbsp; Contact: <a href="mailto:dominique.de-vienne@inrae.fr">dominique.de-vienne@inrae.fr</a>'
    ))
  )
)


# ui <- fluidPage(
#   useShinyjs(),

#   # Custom CSS
#   tags$style(HTML("
#     body {
#       background-color: #f8f9fa;
#       font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
#       color: #333;
#     }
#     .control-label {
#       font-size: 16px;
#       font-weight: bold;
#     }
#     .results-panel {
#       background-color: #f1f1f1;
#       border: 1px solid #ddd;
#       padding: 15px;
#       border-radius: 6px;
#       margin-top: 20px;
#     }
#     .result-label {
#       font-weight: bold;
#     }
#     .well {
#       background-color: #ffffff;
#       border: 1px solid #ccc;
#     }
#     .plots-panel {
#       background-color: #f0f0f0; 
#     }
#   ")),

#   # Theme
#   theme = bs_theme(version = 5, bootswatch = "darkly"),

#   # Main layout
#   fluidRow(
#     column(width = 1),  # Left margin
    
#     column(width = 10,
#       # Title and description
#       fluidRow(
#         column(12,
#           tags$h2("Genetic Counseling for ALS/FTD", style = "text-align: center; margin-top: 20px;"),
#           tags$h4("Carrier Probability Estimation for C9orf72 Repeat Expansion", style = "text-align: center; margin-bottom: 20px;"),
#           tags$p("This application estimates the probability of carrying the C9orf72 repeat expansion mutation for an unaffected first- or second-degree relative of an affected individual. Estimates are based on age-specific penetrance data (Murphy et al., 2017) and formulae from de Vienne & de Vienne (in prep)."),
#           tags$p("For publication use, please cite:", 
#                  tags$i("de Vienne D & de Vienne DM. In prep. Improving Genetic Counseling for C9orf72-ALS/FTD with Age-Based Risk Estimates."))
#         )
#       ),

#       # Inputs and Outputs
#       fluidRow(
#         column(width = 3,
#           wellPanel(
#             tags$h4("Input Parameters", style = "margin-bottom: 15px;"),

#             radioButtons("relation", "Family relationship of the consultand to the affected individual:",
#                          choices = c("Child" = "child", "Grandchild" = "grandchild")),

#             numericInput("Penetrance", label = "Assumed disease penetrance:", value = 1, min = 0, max = 1, step = 0.05),

#             sliderInput("age", "Current age of the patient:", value = 40, min = 0, max = 100, step = 1),

#             conditionalPanel(
#               condition = "input.relation == 'grandchild'",
#               sliderInput("parent_age", "Age of the parent (child of the affected individual):", value = 70, min = 15, max = 100, step = 1),
#               tags$div(
#                 style = "font-size: 90%; color: #444;",
#                 "Age gap: ",
#                 textOutput("age_gap", inline = TRUE),
#                 " years"
#               )
#             ), 

#             numericInput("n_years", 
#               label = "Time frame for disease risk estimation (in years):", 
#               value = 10, 
#               min = 1, 
#               max = 100, 
#               step = 1)
#           )
#         ),
#         column(width = 9,
#           wellPanel(class= "plots-panel",
#             h4("Results"),
#             fluidRow(
#               column(width = 6,
#                 h5("Probability of carrying the c9orf72 RE mutation"),
#                 plotlyOutput("interactive_plot", height = "500px")
#               ),
#               column(width = 6,
#                  uiOutput("dynamic_title_nyears"),
# #                h5("Probability of developing the disease in the next N years"),
#                 plotlyOutput("interactive_plot2", height = "500px")
#               )
#             )
#           ),
#           wellPanel(class = "results-panel",
#             h4("Risk estimates summary", style = "margin-bottom: 15px;"),
#             p(span(class = "result-label", "Probability to carry the mutation: "),
#               textOutput("carrier_prob_final", inline = TRUE)),
#             p(span(class = "result-label", "Probability to develop the disease: "),
#               textOutput("sick_prob_final", inline = TRUE))
#           )

#         )
#       )
#     ),

#     column(width = 1)  # Right margin
#   ),

#   # Footer
# tags$footer(
#   style = "
#     background-color: #f0f0f0;
#     padding: 15px 20px;
#     margin-top: 30px;
#     text-align: center;
#     font-size: 90%;
#     color: #555;
#     border-top: 1px solid #ddd;
#   ",
#   HTML(paste0(
#     '<img src="logo_lbbe.svg" alt="LBBE logo" height="30" style="vertical-align: middle; margin-right: 10px;"> &nbsp;|&nbsp; ',
#     'Source code on <a href="https://github.com/damiendevienne/ftd-als" target="_blank">GitHub</a> &nbsp;|&nbsp; ',
#     'Licensed under the <a href="https://opensource.org/licenses/MIT" target="_blank">MIT License</a> &nbsp;|&nbsp; ',
#     'Contact: <a href="mailto:dominique.de-vienne@inrae.fr">dominique.de-vienne@inrae.fr</a>'
#   ))
# )
# )




############
## SERVER ##
############

# Server
server <- function(input, output, session) {

  # Reactive expression to compute parent age
  age_gap <- reactive({
    req(input$relation == "grandchild", input$age, input$parent_age)
    input$parent_age - input$age
  })

  output$dynamic_title_nyears <- renderUI({
    h5(paste0("Probability of developing the disease in the next ", input$n_years, " years"))
  })

  output$age_gap <- renderText({
    req(input$relation == "grandchild", input$age, input$parent_age)
    input$parent_age - input$age
  })

  adjusted_df <- reactive({ #this dataframe is updated when variables (n, x, relation) change
    x <- 1-input$Penetrance #proportion of non-penetrant carriers in the population
    n <- input$n_years #time frame for disease risk estimation
    penetrance_n <- c(df_onset$Penetrance[(n+1):length(df_onset$Penetrance)],rep(1, n)) #penetrance at age t+n
    # df <- df_onset
    df <- data.frame(Age=df_onset$Age, 
      Penetrance = df_onset$Penetrance,
      Penetrance_n = penetrance_n #penetrance at age t+n
      )
    if (input$relation=="child") {
      df$PvariantWhenHealthy <- pVWH_child(x, df_onset$Penetrance)
      df$PsickInFuture <- pSICKn_child(x, df_onset$Penetrance, penetrance_n)
    }
    if (input$relation=="grandchild") {
      #prepare data to compute probabilities
      delta <- input$parent_age - input$age #age gap
      penetrance_delta <- c(df_onset$Penetrance[(delta+1):length(df_onset$Penetrance)],rep(1, delta))
      df$PvariantWhenHealthy <- pVWH_grandchild(x, df_onset$Penetrance, penetrance_delta) 
      df$PsickInFuture <- pSICKn_grandchild(x, df_onset$Penetrance, penetrance_n, penetrance_delta)
    }
    df$n <- n
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

  # Server outputs plots using make_plot function
  output$interactive_plot <- renderPlotly({
    df <- adjusted_df()
    make_plot(
      df = df,
      age_input = input$age,
      prob_line = get_carrier_prob(),
      y_col = "PvariantWhenHealthy",
      line_color = "steelblue",
      plot_title = "Probability of carrying the c9orf72 RE mutation"
    )
  })

  output$interactive_plot2 <- renderPlotly({
    df <- adjusted_df()
    make_plot(
      df = df,
      age_input = input$age,
      prob_line = get_sick_prob(),
      y_col = "PsickInFuture",
      line_color = "firebrick",
      plot_title = paste("Probability of developing the disease in the next",df$n[1],"years")
    )
  })  

}

shinyApp(ui, server)
