library(shiny)
library(lubridate)
source("BMarray.R")
source("plotHistogram.R")
suppressWarnings(source("DiffMethyl.R"))
library(ggplot2)
library(knitr)
library(dplyr)
library(stringr)
library(heatmap3)

# Load libraries and functions needed to create SQLite databases.
library(RSQLite)
library(RSQLite.extfuns)
saveSQLite <- function(data, name){
  path <- "data/user.sqlite"
  if (!file.exists(path)) {
    message("Caching db at ", path)
    src <- src_sqlite(path, create = TRUE)
    copy_to(src, data, name, temporary = FALSE)
  } else {
    src <- src_sqlite(path)
  }
  return (src)
}

#' sessions
sessions <- data.frame(
  username = c("xinzhou","ruichaohou"),
  password = c("123456", "000000"),
  role     = c("manage", "user"),
  email    = c("xinchoubiology@gmail.com", "tchou@gmail.com"),
  database = c("xinzhou.Rda", "ruichaohou.Rda"),
  tcgaID   = c("TCGA-BH-A1EO,TCGA-BH-A1ES,TCGA-BH-A1ET,TCGA-BH-A1EU,TCGA-BH-A1EV,TCGA-BH-A1EW", "TCGA-BH-A1F0,TCGA-E2-A15I,TCGA-A7-A13G,TCGA-BH-A1EN,TCGA-BH-A1EY,TCGA-BH-A1F2"),
  status   = c("created", "created")
)

src <- saveSQLite(sessions, "sessions")
session_db <- tbl(saveSQLite(sessions, "sessions"), "sessions")

# password match
passwordMatch <- function(user, pwd){
  if(nrow(as.data.frame(filter(session_db, username == user))) == 0){
    return(FALSE)
  }
  #as.data.frame(select(filter(session_db, username == user), password))[1,1] == pwd
  nrow(as.data.frame(filter(filter(session_db, username == user), password == pwd))) == 1
}

# userMatch
userMatch <- function(user){
  dim(filter(session_db, username == user))[1] == 1
}

# role
roleName <- function(user){
  #as.data.frame(select(filter(session_db, username == user), role))[1,1]
  as.data.frame(filter(session_db, username == user))$role
}

# dbConnect
dbGet <- function(user){
  #paste("data", as.data.frame(select(filter(session_db, username == user), database))[1,1], sep = "/")
  paste("data", as.data.frame(filter(session_db, username == user))$database, sep = "/")
}

# bumpdbConnect
bumpdbGet <- function(user, suffix){
  paste("data", paste0(user, suffix, "bump.Rda"), sep = "/")
}

# databaseInsert
dbInsert <- function(userdf){
  sessions = rbind(as.data.frame(session_db), userdf)
  system("rm data/user.sqlite")
  session_db <- tbl(saveSQLite(sessions, "sessions"), "sessions")
  session_db
}

shinyServer(function(input, output) {
  #' Get current user name
  object = NULL
  correlation = NULL
  user <- function(){
    input$submit
    session_db <- tbl(src_sqlite("data/user.sqlite"), "sessions")
    curUser <- input$username
    curPwd  <- input$password
    cat(sprintf("curUser = %s \n", curUser))
    if(curPwd == "" || curUser == ""){
      return(NULL)
    }
    cat(sprintf("curPwd = %s \n", curPwd))
    if(!passwordMatch(curUser, curPwd)){
      return(FALSE)
    }else{
      return(curUser)
    }
  }
  
  isSu <- function(){
    if(is.null(user()) || user() == FALSE){
      return(FALSE)
    }
    role = roleName(user())
    return(role == "manage")
  }
  
  output$userPanel <- renderUI({
    if(input$submit > input$out){
      if(isSu()){
        session_db <- tbl(src_sqlite("data/user.sqlite"), "sessions")
        # The management UI should have a drop-down that allows you to select a 
        # dataset
        tagList(
          HTML(paste0("Logged in as <code>", user(), 
                      "</code> who is a <code>", roleName(user()) ,"</code>.")),
          hr(),
          p("As a manager, you may select any "),
          p("450K TCGA data you wish to view."),
          selectInput("mOrbeta", "Methylation measure:",
                      list("beta","M"),
                      multiple=FALSE, 
                      selected="M")
        )
      }else if(!is.null(user()) && user() != FALSE){
        # It's just a regular user. Just tell them who they are.
        if(file.exists(dbGet(user()))){
          tagList(
            HTML(paste0("Logged in as <code>", user(), "</code> with <code>",
                        "M value measurement</code>.")),
            hr(),
            selectInput("mOrbeta", "Methylation measure:",
                        list("beta","M"),
                        multiple=FALSE, 
                        selected="M"))
        }else{
          tagList(
            hr(),
            HTML(paste0("Your account's database is under processing and building, ")),
            HTML(paste0("please wait for our notice email")),
            hr()
          )
        }
      }else{
        tagList(
          HTML(paste0("Logged in by <code>", "username", "</code> and <code> password </code>")),
          hr()
        )
      }
    }
    else{
      tagList(
        HTML(paste0("Login Methylation 450k Database by <code>", "username", "</code> and <code>", "password", "</code>")),
        hr(),
        HTML(paste0("For example")),
        br(),
        HTML(paste0("manager: xinzhou 123456")),
        br(),
        HTML(paste0("user: ruichaohou 000000"))
      )
    }
  })
  
  output$error <- renderUI({
    if(input$submit > input$out){
      if(is.null(user()) || user() == FALSE){
        tagList(
          br(),
          HTML(paste0("<code>", "pwd incorrect or user inexist", "</code>"))
        )
      }
      else{
        HTML(paste0(""))
      }
    }
  })
  
  output$ui <- renderUI({
    if(input$submit > input$out && !is.null(user()) && user() != FALSE && file.exists(dbGet(user()))){
      tagList(
        sliderInput("breaks", "Breaks of betaOrM histogram", 0, 1, 0.05, 0.01),
        checkboxInput("cond","total v.s tumor-normal"),
        sliderInput("sample", "Sample ratio[500, 5,000]", 500, 5000, 1500, 100),
        sliderInput("nullbreaks", "Breaks of null distribtion histogram", 0, 1, 0.05, 0.01),
        actionButton("draw", "Analysis")
      )
    }else{
      object = NULL
      correlation = NULL
    }
  })
  
  output$measureHist <- renderPlot({
    input$draw
    login = input$draw
    if(is.null(login) || login == 0 || input$submit <= input$out || is.null(user()) || user() == FALSE){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0,1) +ylim(0,1) + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.text =element_blank(), axis.ticks=element_blank())
    }else{
      load(dbGet(user()))
      object = tmp$BMarray
      p <- plotHistogram(object, type = input$mOrbeta, width = input$breaks, cond = input$cond, sampleRatio = input$sample)
      print(p)
    }
  })
  
  output$measureDiff <- renderPlot({
    input$draw
    login = input$draw
    if(is.null(login) || login == 0 || input$submit <= input$out || is.null(user()) || user() == FALSE){
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0,1) +ylim(0,1) + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.text =element_blank(), axis.ticks=element_blank())
    }else{
      if(is.null(object)){
        load(dbGet(user()))
        object = tmp$BMarray
      }
      pl <- plotDiff(object, type = input$mOrbeta, width = input$breaks, sampleRatio = input$sample)
      correlation <- pl$corr
      arrange_ggplot2(pl$p1, pl$p2, ncol = 1)
    }
  })
  
  output$modelPrint <- renderPrint({
    input$draw
    login = input$draw
    if(is.null(login) || login == 0 || input$submit <= input$out || is.null(user()) || user() == FALSE){
      return("Database does not perform M-value packaging")
    }else{
      if(is.null(object)){
        load(dbGet(user()))
        object = tmp$BMarray
      }
      pl <- plotDiff(object, type = input$mOrbeta, width = input$breaks, sampleRatio = input$sample)
      return(summary(pl$corr))
    }
  })
  
  output$icon <- renderImage({
    list(
      src = "images/methyl-logo.png",
      filetype = "image/png",
      alt = "icon"
    )
  }, deleteFile = FALSE)
  
  output$text <- renderText({
    input$draw
    if(is.null(input$draw) || input$submit <= input$out){
      "Welcome to MethylDB Application, MethylDB is an application based on OODB model, it processes Illumina 450k data and find out DMRs"
    }else if(input$draw == 0){
      "MethylDB Database Application For Users, specific different backend database by different username and password"
    }else{
      paste0(user(), "! Welcome to use your Illumina 450k database application MethylDB!")
    }
  })
  
  output$downloadPDF <- downloadHandler(filename = "report.pdf",
                                        content  = function(file){
                                          knit2pdf("report.Rnw")
                                          file.copy("report.pdf", file)
                                          unlink("figure", recursive = TRUE)
                                        },
                                        contentType = "application/pdf")
  
  output$nulldist <- renderPlot({
    input$draw
    login = input$draw
    if(is.null(login) || login == 0 || input$submit <= input$out || is.null(user()) || user() == FALSE){
      #return("Database does not perform M-value packaging")
      df <- data.frame()
      ggplot(df) + geom_point() + xlim(0,1) +ylim(0,1) + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), axis.text =element_blank(),
                                                               axis.ticks=element_blank())
    }else{
      load(bumpdbGet("xinzhou", input$maxGap))
      p <- plotNull(dmrs, width = input$breaks)
      print(p)
    }
  })
  
  ## user information store in data/user.sqlite
  output$registText <- renderText({
    input$register
    apply = input$register
    if(apply < 1 || is.null(input$reguser) || is.null(input$regpassword) || is.null(input$reregpassword) || is.null(input$email) || is.null(input$metadat)){
      ""
    }else if(is.na(str_match(input$email, "@")[1])){
      "Please make sure that your email address is in correct form"
    }else if(userMatch(input$reguser)){
      "username exists"
    }else if(input$regpassword != input$reregpassword){
      "password must be equal to re-password"
    }else{
      metaFile = input$metadat
      tcgaID = as.matrix(read.csv(metaFile$datapath, header = FALSE)[1,])[1,]
      tcgaID = paste(tcgaID, collapse = ",")
      userdf = data.frame(username = input$reguser, password = input$regpassword, 
                          role = "user", email = input$email, database = paste0(input$reguser, ".Rda"), tcgaID = tcgaID, status = "processing")
      dbInsert(userdf)
      system(sprintf("cp %s data/%s.csv", metaFile$datapath, paste0(input$reguser, "meta")), wait = FALSE)
      # modify sessiondb
      ## Rscript DBprocess.R file data/kobemeta.csv username kobe email xinchoubiology@gmail.com
      cat(sprintf("Rscript DBprocess.R file %s username %s email %s \n", paste0("data/", paste0(input$reguser, "meta"), ".csv"), input$reguser, input$email))
      system(sprintf("Rscript DBprocess.R file %s username %s email %s", paste0("data/", paste0(input$reguser, "meta"), ".csv"), input$reguser, input$email), wait = FALSE, ignore.stdout = TRUE)
      HTML(paste("we are preparing for you account's database, and we will notice you via email", input$email, "once it is set"))
    }
  })
  
  output$tabui <- renderUI({
    if(input$submit > input$out && !is.null(user()) && user() != FALSE && file.exists(dbGet(user()))){
      if(isSu()){
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
          tabPanel("CLustering",
                   HTML("Clustering DMR to find possible DMR Cliques under similar regulation"),
                   plotOutput("clustermap")),
          tabPanel("Supervisor",
                   uiOutput("super"))
        )
      }else{
        tagList(
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
                     plotOutput("nulldist"),
                     HTML('<br>'),
                     downloadButton("downloadPDF", "Automate generate DMR report")),
            tabPanel("CLustering",
                     HTML("Clustering DMR to find possible DMR Cliques under similar regulation"),
                     plotOutput("clustermap"))
          )
        )
      }
    }else{
      tagList(
        tabsetPanel(
          tabPanel("Home",
                   imageOutput("icon", height = "220px")),
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
                   textOutput("registText"))
        )
      )
    }
  })
  
  output$super <- renderUI({
    if(input$submit > input$out){
      if(isSu()){
        # The management UI should have a drop-down that allows you to select a 
        # dataset
        tagList(
          HTML(paste0("super manager's database profile")),
          # Create a new Row in the UI for selectInputs
          fluidRow(
            column(4,
                   selectInput("accountid",
                               "username:",
                               c("All", data.frame(session_db)$username)
                   )
            ),
            column(4,
                   selectInput("roles",
                               "role:",
                               c("All", data.frame(session_db)$role)
                   )
            ),
            column(4,
                   selectInput("status",
                               "status:",
                               c("All", c("created", "processing"))
                   )
            )
          ),
          hr(),
          dataTableOutput(outputId="table")
        )
      }else{
        HTML("Welcome to our MethylDB website!")
      }
    }
  })
  
  output$table <- renderDataTable({
    data <- data.frame(session_db)[-2]
    if(input$accountid != "All"){
      data = data[data$username == input$accountid,]
    }
    if(input$roles != "All"){
      data = data[data$role == input$roles,]
    }
    if(input$status != "All"){
      data = data[data$status == input$status,]
    }
    data
  })
})
