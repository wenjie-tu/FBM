## Group Contribution
## Goliath

# Set the environment language
Sys.setenv(lang="us_en")

# Set the working directory
# setwd("F:/UZH/22Spring/FBM/Session3/R")

# Load libraries
library(dplyr)
library(ggplot2)
library(shiny)

# User interface with navigation bar
ui <- pageWithSidebar(
  ## Header
  headerPanel("Monte Carlo Sampling"),
  
  ## Getting the number of trials
  sidebarPanel(
    sliderInput("N",
                "Number of random points:",
                min = 1000,
                max = 19000,
                step = 1000,
                value = 1000),
  ),
  
  ## Displaying the output
  mainPanel(plotOutput("plot1"), textOutput("text1"))
)


# Server logic
server <- function(input, output) {
  
  ## Monte Carlo simulation
  output$plot1 <- renderPlot({
    ## Set seed for reproducible results
    set.seed(2022)
    trials <- input$N
    
    ## Generate random coordinates within the square
    d.points <- data.frame(
      x = runif(trials, min=-1, max=1),
      y = runif(trials, min=-1, max=1)
    )
    
    ## Adding two columns to the data frame
    d.points <- d.points %>% mutate(
      radius = sqrt(x^2 + y^2),                       # add radius variable
      withinCircle = ifelse(radius <= 1, TRUE, FALSE) # add boolean variable
    )
    
    ## Generate plot
    ggplot(d.points) + 
      geom_point(aes(x=x, y=y, color=withinCircle)) +
      coord_fixed(ratio=1) + 
      theme_minimal()
  })
  
  output$text1 <- renderText({
    
    ## Set seed for reproducible results
    set.seed(2022)
    trials <- input$N
    
    ## Generate random coordinates within the square
    d.points <- data.frame(
      x = runif(trials, min=-1, max=1),
      y = runif(trials, min=-1, max=1)
    )
    
    ## Adding two columns to the data frame
    d.points <- d.points %>% mutate(
      radius = sqrt(x^2 + y^2),                       # add radius variable
      withinCircle = ifelse(radius <= 1, TRUE, FALSE) # add boolean variable
    )
    
    ## Estimate the area
    ratio <- mean(d.points$withinCircle == TRUE) 
    
    ## the estimated area of unit circle
    area <- ratio * 4 
    paste("The area of the unit circle", area)
  })
}

shinyApp(ui = ui, server = server)
