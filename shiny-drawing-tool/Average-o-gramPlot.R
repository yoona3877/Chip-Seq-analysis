
library(shiny)
library(gplots)
library(matrixStats)
library(RColorBrewer)
library(tools)


ui <- fluidPage(
  titlePanel("Average-o-gram"),
  
  sidebarPanel(
    fileInput("file", "Choose a file:  ", multiple = FALSE, accept = c("text",".txt")),
    tags$h4("Instruction"),
    tags$p("- Input file must be a txt file."),
    tags$p("- Each average-o-file must be specified with its", tags$strong("absolute path"), "."),
    tags$p("- Dimension of list of average-o-files must be equal."),
    tags$br(),
    tags$h6("* Input file format Example *"),
    tags$em("inputfile.txt"),
    tags$hr(),
    tags$code("/Users/yoona96/average-o-gram1.txt"),
    tags$br(),
    tags$code("/Users/yoona96/average-o-gram2.txt"),
    tags$br(),
    tags$code("/Users/yoona96/average-o-gram3.txt"),
    tags$br(),
    tags$code("/Users/yoona96/average-o-gram4.txt"),
    tags$hr(),
    tags$h6("** Maximum of", tags$strong("8"),"average-o-gram files can be plotted.**")
  ),
  
  mainPanel(
    plotOutput("Plot"),
    downloadButton("save","Save")
  )
)


server<- function(input, output, session){
  Colchoices <- brewer.pal(8, "Set2")
  xvec <- seq(-50,50,1)
  
  data <- reactive({
    inFile <-input$file
    if(is.null(inFile))    return(NULL)
    validate(need(file_ext(inFile$name) == "txt", 
                  "The input file is not a txt file. \n
                  Check the file type again. It must be a txt file.\n
                  Please follow the instruction."))
    read.table(inFile$datapath, header=FALSE)
  })
  
  calculate <- function(dataplotting){
    mean <-colMedians(as.matrix(dataplotting),na.rm=T)
    n <- sqrt(colSums((!is.na(as.matrix(dataplotting)))))
    sd <- colSds(as.matrix(dataplotting),na.rm=T)
    dt <- data.frame(mean,n,sd)
    return(dt)
  }

  gramoPlot <- function(){
    if (is.null(data())) return(NULL)
    
    files <- data()
    fileNames <- list()
    for (i in 1 : nrow(files)){
      fileNames[[i]] = files[i,]
      validate(need(file_ext(basename(as.character(fileNames[[i]]))) == "txt", 
                    "The list of average-o-gram file in the input file is not a txt file. \n
                    Check the input file again.\n All the average-0-files must be a txt file.\n
                    Please follow the instruction."))
    }
    
    # Alert below will trigger if the length of the files is above 8
    validate(
      need(length(fileNames) <=8, 
           "Maximum of 8 average-o-gram files can be plotted\n
           Maximum of 8 average-o-gram files can be plotted\n
           Please try again\nPlease try again")
    )
    
    dataList <- list()
    for (i in 1 : length(fileNames)){
      # Check whether the files exist. If not, error message will appear. 
      validate(
        need(file.exists(as.character(fileNames[[i]])), 
             "Average-o-gram files do not exist!\n
             Check the absolute path of the avarage-o-gram file again")
      )
      dataList[[i]] = read.table(as.character(fileNames[[i]]), header = TRUE)
    }
    
    dataplotting <- dataList[[1]]
    dataplotting[dataplotting == c(-1)] <- NA
    dataplotting$type <- NULL
    
    first <- calculate(dataplotting)
    
    plot(xvec,first$mean,lwd=3,type="l",col=Colchoices[1],ylim=c(0.5,2),ylab="Copy Number",xlab="Distance (in KB)",xlim=c(-50,50))
    arrows(xvec,first$mean-first$sd/first$n,xvec,first$mean+first$sd/first$n,length=0.05, angle=90, code=3, col=Colchoices[1])
    
    if (length(dataList) >1){
      for (i in 2: length(dataList)){
        
        # Check whether the dimension of all the files are equal. If not, error message will appear. 
        validate(
          need(dim(dataList[[1]]) == dim(dataList[[i]]), 
               "Not appropriate average-o-gram files.\n
               Dimension of the list of average-o-files do not match.\n
               Please check again.")
        )
        plotting <- dataList[[i]]
        plotting[plotting == c(-1)] <- NA
        plotting$type <- NULL
      
        others <- calculate(plotting)
      
        lines(xvec,others$mean,lwd=3,col=Colchoices[i])
        arrows(xvec,others$mean-others$sd/others$n,xvec,others$mean+others$sd/others$n,length=0.05, angle=90, code=3, col=Colchoices[i])
        }
    }
    
    legendName <- array()
    for (i in 1 : length(fileNames)){
      legendName[i] = basename(as.character(fileNames[[i]]))
    }
    
    legend("topright",lwd=3,
           legend=legendName, 
           lty=1, col=c(Colchoices[1:8]))
  }
  
  output$Plot <- renderPlot({
    gramoPlot()
  })
  
  output$save <- downloadHandler(
    filename = function() { paste(input$files, '.pdf', sep='') },
    content = function(file) {
      pdf(file)
      gramoPlot()
      dev.off()
    }
  )
}


shinyApp(ui = ui, server = server)