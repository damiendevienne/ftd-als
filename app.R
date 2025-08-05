require(shiny)
library(shinyjs) # For writing messages to console
library(plotly) # For interactive plots
library(bslib)
library(rmarkdown) #for pdf report
library(knitr) #for pdf report

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


make_combined_plot <- function(df, age_input, carrier_prob, sick_prob, nyears) {
  ymax <- max(carrier_prob, sick_prob)
  
  plot_ly(data = df, type = 'scatter', mode = 'lines') %>%
    add_trace(x = ~Age, y = ~PvariantWhenHealthy,
              name = "Carrier probability",
              line = list(color = "steelblue")) %>%
    add_trace(x = ~Age, y = ~PsickInFuture,
              name = paste("Disease risk (next", nyears, "years)"),
              line = list(color = "firebrick")) %>%
    layout(
      margin = list(l = 60, r = 10, t = 10, b = 60),
      xaxis = list(title = "Age of the consultand"),
      yaxis = list(title = "Estimated probability", range = c(0, 1)),
      shapes = list(
        list(type = "line", x0 = age_input, x1 = age_input, y0 = 0, y1 = ymax,
             line = list(dash = "dot", color = "gray")),
        list(type = "line", x0 = 0, x1 = age_input, y0 = carrier_prob, y1 = carrier_prob,
             line = list(dash = "dot", color = "steelblue")),
        list(type = "line", x0 = 0, x1 = age_input, y0 = sick_prob, y1 = sick_prob,
             line = list(dash = "dot", color = "firebrick"))
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
    config(displayModeBar = FALSE, editable = FALSE)
}


####################
## USER INTERFACE ##
####################

ui <- fluidPage(

  useShinyjs(),

  tags$head(
    tags$style(HTML("
      html, body {
        height: 100%;
        margin: 0;
      }

      .container-fluid {
        display: flex;
        flex-direction: column;
        min-height: 100vh;
      }

      .main-content {
        flex: 1;
      }
      .summary {
        border: 1px solid rgba(0, 0, 0, 0.176);
        padding: 0px 16px 16px 16px;
      } 
      .loader1 {
        visibility: hidden; /* Hidden by default */
      }
      .loader2 {
        visibility: visible; /* visible by default (while loading plot at refresh) */
      }
      .loader {
        position: absolute;
        top: 50%;
        transform: translate(2px, -50%);
        width: 24px;
        height: 24px;
        border: 4px solid lightgrey;
        border-bottom-color: grey;
        border-radius: 50%;
        box-sizing: border-box;
        animation: rotation 1s linear infinite;
      }

      @keyframes rotation {
        0% {
          transform: translate(10px, -50%) rotate(0deg);
        }
        100% {
          transform: translate(10px, -50%) rotate(360deg);
        }
      }

      .help-block {
        font-size: 75%;
      }

      footer {
        margin-top: auto;
      }

    "))
  ),

  # Bootstrap Theme
  theme = bs_theme(version = 5, bootswatch = "yeti"),  # Change here for different looks


  # Header

  tags$header(
      class = "bg-primary text-white py-4 mb-4",
      div(class = "container",
        tags$h2("Genetic Counselling for ALS/FTD", class = "mb-2"),
        tags$h4(HTML("Probability Estimates for <i>C9orf72</i><sup>RE</sup>-related ALS/FTD diseases"), class = "mb-3"),
        div(
          class = "small",
          HTML(paste0(
            'Source code on <a href="https://github.com/damiendevienne/ftd-als" target="_blank" class="text-white text-decoration-underline">GitHub</a> &nbsp;|&nbsp; ',
            'Licensed under the <a href="https://opensource.org/licenses/MIT" target="_blank" class="text-white text-decoration-underline">MIT License</a> &nbsp;|&nbsp; ',
            'Contact <a href="mailto:dominique.de-vienne@universite-paris-saclay.fr" target="_blank" class="text-white text-decoration-underline">Dominique de Vienne</a>'
          ))
        )
      )
    ),

  # Main layout
  div(class="main-content",
    fluidRow(
      column(width = 2),  # Left margin

      column(width = 8,
        # Title and description
        fluidRow(
          column(12,
            tags$p(HTML("This application estimates the probability of carrying the <i>C9orf72</i><sup>RE</sup> mutation for unaffected relatives of a carrier. Estimates are based on the mathematical developments from de Vienne & de Vienne (in prep), applied to the age-specific penetrance data from Murphy <i>et al.</i> (2017). The risk estimates in this report are not definitive predictions. They reflect current knowlegde on age-related penetrance and may be updated as new data become available.")),
            tags$p(HTML("For publication use, please cite: <i>de Vienne & de Vienne. 2025. Age-Based Risk Estimates for </i>C9orf72<i><sup>RE</sup>-related Diseases: Theoretical Developments and Added Value for Genetic Counselling. medRxiv</i>"))
          )
        ),


# autre idée de titre : Improving Genetic Counselling for C9orf72-ALS/FTD with Age-Based Risk Estimates.

        # Inputs and Outputs
        fluidRow(
          div(class = "col-12 col-xl-3 mb-4",
            wellPanel(
              tags$h4("Input Parameters", class = "mb-3"),

              radioButtons("relation", "Relationship of the consultand to the affected carrier:",
                          choices = c("Child (or sibling)" = "child", "Grandchild (or nibling)" = "grandchild")),

              numericInput("Penetrance", label = "Assumed proportion of non-penetrant carriers:", value = 0, min = 0, max = 1, step = 0.05),
              sliderInput("age", "Current age of the consultand:", value = 40, min = 0, max = 100, step = 1),

              conditionalPanel(
                condition = "input.relation == 'grandchild'",
                sliderInput("parent_age", "Age of the parent (child - or sibling - of the affected individual):", value = 70, min = 15, max = 100, step = 1),
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

          div(class = "col-12 col-xl-5 mb-4",
            wellPanel(
              tags$h4("Risk estimates visualization"),
              #withSpinner(plotlyOutput("combined_plot", height = "500px"), type=8, color="grey", size=0.5),
              div(class = "loader loader2"),
              plotlyOutput("combined_plot", height = "500px"),

              div(style = "text-align: right;",
              helpText("Cick and drop to zoom in, double-click to reset zoom")
              )
            )
          ),

          div(class = "col-12 col-xl-4 mb-4",
            wellPanel(
              div(class="summary",
              tags$h4("Risk estimates summary", class = "mb-3"),

              tags$label("Relationship to the affected individual:", textOutput("summary_relation", inline = TRUE)),
              tags$br(),
              tags$label("Current age of the consultand:", textOutput("summary_age", inline = TRUE)),
              tags$br(),
              tags$label("Assumed proportion of non-penetrant carriers:", textOutput("summary_penetrance", inline = TRUE)),

              conditionalPanel(
                condition = "input.relation == 'grandchild'",
                tags$label("Parent's age:", textOutput("summary_parent_age", inline = TRUE)),
                tags$br(),
                tags$label("Age gap (parent - consultand):", textOutput("summary_age_gap", inline = TRUE))
              ),

              #tags$label(tags$strong("Time frame for disease risk estimate:"), textOutput("summary_nyears", inline = TRUE)),

              tags$hr(),

              tags$label(tags$strong("Estimated probability of carrying the mutation:"), 
              span(style = "
                  color: black;
                  font-weight: bold;
                  background-color: #e6f0fa;
                  border: 1px solid steelblue;
                  padding: 4px 8px;
                  border-radius: 4px;
                  display: inline-block;",
                textOutput("carrier_prob_final", inline = TRUE)
              )),
              tags$br(),
              tags$label(tags$strong("Estimated probability of developing the disease in the next ", 
                textOutput("summary_nyears", inline = TRUE),
                ": "),
              span(style = "
                  color: black;
                  font-weight: bold;
                  background-color: #fceaea;
                  border: 1px solid firebrick;
                  padding: 4px 8px;
                  border-radius: 4px;
                  display: inline-block;",
                textOutput("sick_prob_final", inline = TRUE)
              )), 
              div(class="mt-4", style = "text-align: center; position: relative;",
                  downloadButton("download_report", "Download report (pdf)"),
                  span(class="loader loader1")
              )
            ))
          )
        )
      ),
      column(width = 2)  # Right margin
    )
  ),  # End of main content. NOTE this div allows for the footer to be at the bottom of the page

  # Footer
  tags$footer(
    class = "text-center py-3 border-top mt-4",
    HTML(paste0(
      '<img src="logo_lbbe.svg" alt="LBBE logo" height="30" style="vertical-align: middle; margin-right: 10px;">'
    ))
  )
)


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

  #new outputs for summary
  output$summary_relation <- renderText({
    if (input$relation == "child") "Child (or sibling)" else "Grandchild (or nibling)"
  })

  output$summary_age <- renderText({
    paste(input$age, "years")
  })

  output$summary_penetrance <- renderText({
    paste0(input$Penetrance)
  })

  output$summary_parent_age <- renderText({
    paste(input$parent_age, "years")
  })

  output$summary_age_gap <- renderText({
    paste(input$parent_age - input$age, "years")
  })

  output$summary_nyears <- renderText({
    paste(input$n_years, "years")
  })

  adjusted_df <- reactive({ #this dataframe is updated when variables (n, x, relation) change
    x <- input$Penetrance #proportion of non-penetrant carriers in the population
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

  # Server outputs plots using make_combined_plot function
  output$combined_plot <- renderPlotly({
    df <- adjusted_df()   
    pl <- make_combined_plot(
      df = df,
      age_input = input$age,
      carrier_prob = get_carrier_prob(),
      sick_prob = get_sick_prob(),
      nyears = df$n[1]
      )
    htmlwidgets::onRender(pl, "
    function(el, x) {
      $('.loader2').css('visibility', 'hidden');
    }
  ")
   })

output$download_report <- downloadHandler(
  filename = function() {
    paste0("genetic_report_", Sys.Date(), ".pdf")
  },
  content = function(file) {
    # Show spinner
    runjs("$('.loader1').css('visibility', 'visible');")

    tempReport <- file.path(tempdir(), "report-template.Rmd")
    file.copy("report-template.Rmd", tempReport, overwrite = TRUE)
    
    params <- list(
      relation = if (input$relation == "grandchild") "Grandchild (or nibling)" else "Child (or sibling)",
      age = input$age,
      parent_age = if (input$relation == "grandchild") input$parent_age else NA,
      penetrance = input$Penetrance,
      nyears = input$n_years,
      carrier = get_carrier_prob(),
      risk = get_sick_prob(),
      date_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      thedataframe = adjusted_df()
    )

    rmarkdown::render(
      input = tempReport,
      output_file = file,
      params = params,
      envir = new.env(parent = globalenv())
    )
    # Hide spinner after a delay to ensure UI updates
    runjs("$('.loader1').css('visibility', 'hidden');")
  }
  # content = function(file) {
  #   params <- list(
  #     relation = input$relation,
  #     age = input$age,
  #     parent_age = if (input$relation == "grandchild") input$parent_age else NA,
  #     penetrance = input$Penetrance,
  #     nyears = input$n_years,
  #     carrier = get_carrier_prob(),
  #     risk = get_sick_prob(),
  #     date_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), 
  #     thedataframe = adjusted_df()
  #   )

  #   rmarkdown::render(
  #     input = "report-template.Rmd",
  #     output_file = file,
  #     params = params,
  #     envir = new.env(parent = globalenv())
  #   )
  # }
)}

shinyApp(ui, server)
