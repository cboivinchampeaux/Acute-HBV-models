library(shiny)
# install.packages("shinydashboard")
library(shinydashboard)
# install.packages("shinyWidgets")
library(shinyWidgets)

sidebar <- dashboardSidebar( width=500,
                             hr(),
                             sidebarMenu(id= "tabs",
                                         menuItem("Model", icon = icon("chevron-circle-right"),
                                                  checkboxInput("II", "Two Infected cells populations (Model 5)", TRUE),
                                                  checkboxInput("R", "Refractory cells population (Model 6)", FALSE),
                                                  checkboxInput("A", "Antibody population (Model 7a)", FALSE),
                                                  checkboxInput("AE", "Antibody and Effector cells populations (Model 7b)", FALSE),
                                                  numericInput("num5", label = h4("Duration of Model (Days)"), value = 300)
                                         )),
                             
                             hr()
)


body <- dashboardBody(
    fluidRow(
    column(6,
           tabName = "plot1", 
           box(width = NULL, plotOutput("plot1"), collapsible = TRUE, 
               title = "Simulated Free virus vs. Time Curve - Patient 1", status = "primary", solidHeader = TRUE)
    ),
    column(6,
           tabName = "plot2", 
           box(width = NULL, plotOutput("plot2"), collapsible = TRUE, 
               title = "Simulated Free virus vs. Time Curve - Patient 2", status = "primary", solidHeader = TRUE)
    )
    ),
    fluidRow(
      column(6,
             tabName = "plot3", 
            box(width = NULL, plotOutput("plot3"), collapsible = TRUE, 
                title = "Simulated Free virus vs. Time Curve - Patient 3", status = "primary", solidHeader = TRUE)
      ),  
    column(6,
           tabName = "plot5", 
           box(width = NULL, plotOutput("plot5"), collapsible = TRUE, 
               title = "Simulated Free virus vs. Time Curve - Patient 5", status = "primary", solidHeader = TRUE)
           )
    ),
    fluidRow(
      column(6,
             tabName = "plot6", 
             box(width = NULL, plotOutput("plot6"), collapsible = TRUE, 
                 title = "Simulated Free virus vs. Time Curve - Patient 6", status = "primary", solidHeader = TRUE)
      ),  
      column(6,
             tabName = "plot7", 
             box(width = NULL, plotOutput("plot7"), collapsible = TRUE, 
                 title = "Simulated Free virus vs. Time Curve - Patient 7", status = "primary", solidHeader = TRUE)
      )
     )
)
     

dashboardPage(
  dashboardHeader(title= "Simulation Acute HBV models", titleWidth = 350),
  skin= "purple",
  sidebar,
  body
)