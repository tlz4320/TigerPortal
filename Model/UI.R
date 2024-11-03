ui <- navbarPage(
  title="TIGER portal",
  id="pages",
  tabPanel(title="Predict",
           fluidPage(
             tags$head(
               tags$style(HTML("
                              #tableOutput table {
                                 word-break:break-word;
                               }
                               #tableOutput th:first-child {
                                 width: 50%;
                               }
                               #tableOutput th:nth-child(2) {
                                 width: 10%;
                               }
                               #tableOutput th:nth-child(3) {
                                 width: 10%;
                               }
      "))),
             fluidRow(column(2),
                      column(4,textInput("seq_input", "MutationID(RS)", 
                                         width = "100%", value = "rs397515581")),
                      column(1, textInput("ngg_seq", "PAM", value = "NGG", width = "100%")),
                      column(2,selectInput("database", "CellType:",
                                           c("U2OS" = "U2OS",
                                             "HEK293" = "HEK293"))),
                      height = "10%", align = "center"
             ), 
             fluidRow(column(4),
                      column(2, actionButton("searchSg", "Search gRNA", width = "100%")),
                      column(2,checkboxInput("reverse", "Reverse Complement", value = T)),
                      uiOutput('selAll'),
                      height = "10%", align = "center"
             ),
             fluidRow(column(2), 
                      column(8, 
                             tableOutput('tableOutput')),
                      height = "40%", align = "center"
             ),
             fluidRow(uiOutput('button'),
                      width = "20%", height = '10%', align = "center")
           )
  ),
  tabPanel("Result", fluidPage(
    fluidRow(h1("Please do Predict Fist!", id = "twarn"),
             uiOutput('resultList'),
             width = "20%", height = '10%', align = "center")
  ))
)
