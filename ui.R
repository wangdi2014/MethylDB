library(shiny)

shinyUI(
  fluidPage(
    # Setup the page title
    tagList(tags$head(tags$title("Illumina Methylation Database")), h1(textOutput("title"))),    
    
    fluidRow(
      column(3, 
             wellPanel(includeCSS('boot.css'),
                       includeHTML('header.html'),
                       "A Illumina Methylation Database Application",
                       uiOutput("userPanel"),
                       hr(),
                       textInput("username", "Account Login: ", "xinzhou"),
                       tags$input(id="password",type="password",placeholder="Password", value = "****"),
                       hr(),
                       actionButton("submit", "SignIN"),
                       actionButton("out", "SignOut"),
                       uiOutput("error"),
                       hr(),
                       uiOutput("ui")
             )),
      column(9,
             mainPanel(
               #uiOutput("tabui")
               tabsetPanel(
                 tabPanel("Home",
                          textOutput("text"),
                          imageOutput("icon", height = "220px")),
                 tabPanel("Statistic",
                          HTML("Statisctic of elementary statistic ."),  # draw the histogram of user specific
                          HTML('<table border=0 width="100%"><tr bgcolor="#f5f5f5"><td>'),
                          div(style="width:100%;max-width:600px;",
                              plotOutput("measureHist")),
                          HTML('</td><td>'),
                          plotOutput("measureDiff"),
                          HTML('</td></tr></table>'),
                          HTML('<br>'),
                          HTML("<p style=\"font-size:17px\">Linear Association between Log(Delta(Measurement)) and their probability :</span></p>"),
                          verbatimTextOutput(outputId = "modelPrint")),
                 tabPanel("Diffebumper",
                          HTML("Statisctic of Differential Bumper."),
                          numericInput("topk", label = "Select Top k Significant DMR :",
                                       value = 10, min = 10, max = 1000, step = 10),
                          numericInput("maxGap", label = "Select Clusters' maxGap [500, 2000] :",
                                       value = 500, min = 500, max = 2000, step = 500),
                          plotOutput("nulldist", height = "220px"),
                          HTML('<br>'),
                          downloadButton("downloadPDF", "Automate generate DMR report")),
                 tabPanel("Clustering",
                          HTML("Clustering DMR to find possible DMR Cliques under similar regulation"),
                          HTML("Inter-DMRs correlation of your database"),
                          HTML("<br>"),
                          selectInput("methods", "Cluster measurement",
                                      list("median","area"),
                                      multiple=FALSE, 
                                      selected="median"),
                          selectInput("color", "Heatmap pallete",
                                      list("GreenRed","BlueRed"),
                                      multiple=FALSE, 
                                      selected="GreenRed"),
                          plotOutput("clustermap")),
                 tabPanel("Register",
                          HTML("If you want to create a new methyl database and find out their DMR, please register to create a new DB for youself"),
                          textInput("reguser", label = "username"),
                          HTML('<br>'),
                          HTML("password"),
                          HTML('<br>'),
                          tags$input(id="regpassword", type="password",placeholder="Password", value = ""),
                          HTML('<br>'),
                          HTML("re-password"),
                          HTML('<br>'),
                          tags$input(id="reregpassword", type="password",placeholder="Password", value = ""),
                          HTML('<br>'),
                          textInput("email", label = "Email"),
                          HTML('<br>'),
                          HTML('<hr>'),
                          downloadButton("downloadID", "Download Candidate TCGA ID"),
                          HTML('<br>'),
                          fileInput('metadat', 'MethylDB sample ID from TCGA(please use common to seperate)',
                                    accept = c(
                                      'text/csv',
                                      'text/comma-separated-values',
                                      'text/tab-separated-values',
                                      'text/plain',
                                      '.csv',
                                      '.tsv'
                                    )
                          ),
                          actionButton(inputId = "register", label = "Apply"),
                          textOutput("registText")),
                 tabPanel("Supervisor",
                          uiOutput("super"))
               )
             )
      )
    )
))

