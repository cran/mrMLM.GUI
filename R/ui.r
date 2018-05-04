ui<-(
  fluidPage(
    
    shinyjs::useShinyjs(),
    radioButtons("layout", "",
                 choiceNames = list(
                   strong("mrMLM"),
                   strong("Start")
                   
                 ),
                 choiceValues = list(
                   "mrMLM", "Start"
                 ),inline=TRUE),
    uiOutput("general_ui")
  )
)
