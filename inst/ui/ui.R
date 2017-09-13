library(shiny)
library(rgl)

shinyUI(fluidPage(
  titlePanel("BIGL"),
  br(),
  sidebarPanel(
    ## Hill slope inputs
    fluidRow(
      column(4,
             numericInput("h1", label = "Hill_1", 1, step = 0.1)),
      column(4,
             numericInput("h2", label = "Hill_2", 1, step = 0.1))
    ),
    ## EC50 inputs
    fluidRow(
      column(4,
             numericInput("e1", label = "EC50_1", 0.1, step = 0.05)),
      column(4,
             numericInput("e2", label = "EC50_2", 0.1, step = 0.05))
    ),
    ## Choice of baseline
    fluidRow(
      column(4,
             numericInput("b", label = "Baseline", 0, step = 0.05))
    ),
    ## Choice of asymptotes
    fluidRow(
      column(4,
             numericInput("m1", label = "Asymptote_1", 1, step = 0.05)),
      column(4,
             numericInput("m2", label = "Asymptote_2", 1, step = 0.05))
    ),
    ## Choice of parameters for null data simulation
    fluidRow(
      column(4,
             numericInput("noise", label = "Noise level", 0, step = 0.05)),
      column(4,
             numericInput("replicates", label = "Replicates", 1, step = 1))
    ),
    ## Whether doses are evenly spaced in log-scale
    checkboxInput("logScale", "Doses are evenly spaced in log-scale", TRUE),
    ## Null model
    radioButtons("null", label = "Null model",
                 choices = c("Generalized Loewe" = "loewe",
                             "Classical Loewe" = "stdloewe",
                             "Highest Single Agent" = "hsa"),
                 selected = "loewe"),
    ## Pre-selected parameters
    strong("Pre-defined examples"),
    br(),
    actionButton("defaultpars", "Agonist / Agonist"),
    br(),
    actionButton("pagonist", "Agonist / Partial agonist"),
    br(),
    actionButton("antagonist", "Agonist / Antagonist")
  ),

  mainPanel(
    tabsetPanel(
      tabPanel("Plots",
               br(),
               strong("Monotherapy coefficients and estimated dose-response curves"),
               br(),
               br(),
               fluidRow(
                 column(3,
                        tableOutput("coefs")),
                 column(9,
                        plotOutput("marginals"))),
               br(),
               strong("Isobologram of the null model"),
               br(),
               br(),
               plotOutput("isobologram"),
               br(),
               br(),
               strong("Expected response surface"),
               rglwidgetOutput("surface", width = "1024px", height = "800px")
               ),

      tabPanel("Table",
               br(),
               DT::dataTableOutput("occuptable")),

      tabPanel("About",
               br(),
               p(strong("BIGL"), "shiny application takes as input parameter values for two ",
                 "4-parameter logistic dose-response curves which share the same baseline."),

               p("In order to test estimation stability, the application allows to add noise to the data ",
                 "generated according to the null model as well as a number of replicates for ",
                 "each dose combination. ", strong("Plots"), " and ", strong("Table"), "tabs use ",
                 "coefficients fitted to this dataset."),

               p("Choice of null model determines which model will be used to construct the ",
                 "monotherapy and response surface plots. This choice has no real impact for ",
                 "the ", strong("Table"), " tab which is always constructed based on Loewe ",
                 "additivity model."),

               p("Doses for both compounds are assumed to be the same and are generated using ",
                 "either ", code("round(seq(0, 3, length.out = 7), 4)"), " or ",
                 code("round(c(0, 3^(-6:0)), 4)"), "commands for both compounds depending on ",
                 "whether evenly spaced logarithmic scale is chosen."),

               br(),
               h4("Plots"),
               p("This tab includes a plot of the monotherapy curves with the estimated monotherapy ",
                 "coefficients. Setting noise level to zero will lead estimated coefficients to be ",
                 "numerically very close, if not identical, to the true coefficients. A table of true ",
                 "and estimated coefficients is provided as well."),
               p("Additionally, 3-dimensional response surface is plotted using the estimated ",
                 "coefficients. Its gradient color represents occupancy values at a given dose ",
                 "combination according to the generalized Loewe model. Transparent color indicates ",
                 "occupancy close to zero, whereas dark blue indicates it being close to one."),
               p("Choice of the null model will be reflected in all components of this tab. If ",
                 "classical Loewe is selected, asymptote estimates are constrained to be the same ",
                 "for both compounds. For generalized Loewe and Highest Single Agent models this ",
                 "restriction does not apply. In the case of 3-dimensional plot, generalized Loewe ",
                 "and classical Loewe models imply that response surface is constructed according ",
                 "to the procedure depicted ",
                 "in the ", strong("Table"), " tab. Expected response for the Highest Single Agent ",
                 "model, on the other hand, is constructed simply by taking either the minimum ",
                 "(if dose-response curves are decreasing) or the maximum (if dose-response curves are ",
                 "increasing) of dose-response values at a particular dose combination."),

               br(),
               h4("Table"),
               p(strong("Table"), " tab contains a detailed summary of how expected response is ",
                 "constructed under the null of either generalized or classical Loewe models given the ",
                 "estimated parameters. HSA model is not represented in this table."),

               p("By definition, expected response is equal to the sum of the baseline and contributions ",
                 "from each compound. Contribution of a compound is obtained by multiplying its monotherapy ",
                 "response at a given dose combination by its weight computed from the occupancy equation.")

               )
    )
  )
))
