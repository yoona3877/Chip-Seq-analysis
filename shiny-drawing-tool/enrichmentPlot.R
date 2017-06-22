
library(shiny)
library("ggplot2")
library(tools)
options(shiny.maxRequestSize=30*1024^2)


ui <- fluidPage(
  titlePanel("Enrichment Plot"),
  fluidRow(
  
    column(8,wellPanel(
      tags$h4("First input file format"),
      tags$p(" - The file name must ends with ", tags$strong("--chromosome_name.txt")),
      tags$em("** File name example : ",
              tags$code("Brown_AB185_S23.vs.Brown_AB186_S24", tags$strong("--chrV.txt")), " **"),
      tags$p(" - The file must have four column. "),
      tags$h4("Second input file format"),
      tags$p(" - The second input file must be a gff file."))),
    column(4,wellPanel(
           fileInput("inputfile", "Choose the first input file:  ", accept = c("text",".txt")),
           fileInput("inputfile2", "Choose the second input file: ", accept = c("text", ".gff"))
           )),
    column(4, wellPanel(textInput("start", "Start Sequence (MB)", value = 0),
                        textInput("end", "End Sequence (MB)", value = 0))),
    column(4,tags$h4("Select : "),
           checkboxInput(inputId = "ars", label = "ARS", value = TRUE),
           checkboxInput(inputId = "ltr", label = "long_terminal_repeat", value = TRUE),
           checkboxInput(inputId = "centromere", label = "centromere", value = TRUE),
           checkboxInput(inputId = "rr", label = "repeat_region", value = TRUE)),
    column(4,checkboxInput(inputId = "rtt", label = "retrotransposon", value = TRUE),
           checkboxInput(inputId = "rRNA", label = "rRNA", value = TRUE),
           checkboxInput(inputId = "telomere", label = "telomere", value = TRUE),
           checkboxInput(inputId = "teg", label = "transposable_element_gene", value = TRUE),
           checkboxInput(inputId = "tRNA", label = "tRNA", value = TRUE)),
    column(12,mainPanel(
      plotOutput("enrichmentPlot"),
      textOutput("sequence"),
      textOutput("error"),
      downloadButton("save","Save")))
  ))



server<- function(input, output, session){
  offset <- 1;
  
  
  data <- reactive({
    inFile <- input$inputfile
    if (is.null(inFile)) return(NULL)
    
    validate(need(file_ext(input$inputfile$name) == "txt", 
                  "The first input file is not a txt file. \n
                  Check the file type again. It must be a txt file"))
    
    baseFileName <- input$inputfile$name
    chr_name <- substr(baseFileName, regexpr("--",baseFileName)+2, regexpr(".txt",baseFileName)+3)
    
    validate(need(regexpr("--", baseFileName), "The first input file name must ends with '--chromosome_name.txt'\n
                  Please follow the instruction."))
    
    read.table(inFile$datapath)
  })
  
  
  data2 <- reactive({
    inFile2 <- input$inputfile2
    if (is.null(inFile2)) return(NULL)
    validate(need(file_ext(input$inputfile2$name) == "gff", 
                  "The second input file is not a gff file. \n
                  Check the file type again. It must be a gff file"))
    read.table(inFile2$datapath, header = FALSE, fill = T, sep = "\t")
  })
  
  
  plotIndividual <- function(chrData, bd, color, start, fin, s, e){
      for (i in 1 : nrow(chrData)){
        startSequence <- chrData[i, 4]/10^6
        endSequence <- chrData[i,5]/10^6
        if (startSequence >= s && endSequence <= e){
          rect(chrData[i, 4]/10^6, start + offset, chrData[i,5]/10^6, 
               fin + offset, col = color, lwd = 0, border = bd)
        }
      }
  }

  
  plot2 <- function(starting, ending){
    if (is.null(data())) return(NULL)
    
    base <- data()
    
    validate(
      need(ncol(base) == 4, "Not appropriate enrichment file.\n
                            The firt input file does not have column of 4. \n
                            Check the input file again." )
    )
    baseFileName <- input$inputfile$name
    out_file <- substr(baseFileName, 0, regexpr(".txt",baseFileName))
    
    # End of the sequence
    end <- min(max(base[,1]/10^6), as.numeric(ending))
    if (end == 0){
      end <- max(base[,1]/10^6) 
    }
    # Start of the sequence
    start <- max(0, as.numeric(starting))
    
    validate(need(end >= start, "End sequence must be greater or equal to start sequence."))
    
    starts <- start*10^6/50
    ends <- end*10^6/50
    
    enrich = base[,2]/base[,3]
    
    validate(
      need(enrich == base[,4], "Check the file data again.\n
           The second column divided by the third colum mush be equal to the fourth column")
    )
    
    enrich <- log2(enrich)
    enrich[enrich>(2+offset)] <- 2+offset
    enrich[enrich<0] <- 0
    enrich <- enrich[starts: ends]
  
    xaxis <- base[,1]/10^6
    xaxis <- xaxis[starts: ends]
    
    plot(xaxis, enrich, pch=".", ylim=c(0,4+offset),
         xlim=c(start,end*1.18), xlab="Chromosomal Location (MB)",
         ylab="Enrichment Score", xaxs="i", yaxs="i", axes=FALSE, main=out_file)
    axis(side=2, at = c(0,0.5,1,1.5,2,2.5,3))
    axis(side=1, at = seq(start, end, 0.05))
    lines(xaxis,enrich,type="h",col="blue")
    lines(c(start,end),c(2.05+offset,2.05+offset))
    lines(c(end,end),c(0,init[9]+0.1+offset))

    output$sequence <- renderText({
      paste("Sequence : 0 to ", base[nrow(base),1]/10^6, sep = "")
    })
  }
    
  
  plot3 <- function(starting, ending){
    if (is.null(data2())) return(NULL)
    
    if (is.null(data())){
      output$error <- renderText({
        "Data will show up after the first input file is uploaded"
      })
    }
    else{
      output$error <- renderText({
        ""
      })
    }
    init <- c(2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1,2.1)
    legendN <- array()
    color <- array()
    base <-data()
    baseFileName <- input$inputfile$name
    chr_name <- substr(baseFileName, regexpr("--",baseFileName)+2, regexpr(".txt",baseFileName)-1)
    chrData <- data2()
    chrData <- subset(chrData, V1 == chr_name)
    
    e <- min(max(base[,1]/10^6), as.numeric(ending))
    if (e == 0){
      e <- max(base[,1]/10^6) 
    }
    # Start of the sequence
    s <- max(0, as.numeric(starting))
      
    if (input$ars){
      chrData1 <- subset(chrData, V3 == "ARS")
      if (nrow(chrData1) != 0){
        plotIndividual(chrData1, 
                       c("red"), "red", start = init[1], fin = init[1] + 0.2, s , e)
        init[2:length(init)] <- init[2:length(init)] + 0.2
      }
      legendN[length(legendN) + 1] <- "ARS"
      color[length(color) + 1] <- "red"
    }
      
    if (input$ltr){
      chrData2 <- subset(chrData, V3 == "long_terminal_repeat")
      if (nrow(chrData2) != 0){
        plotIndividual(chrData2,
                       c("darkgreen"), "darkgreen", start = init[2], fin = init[2] + 0.2, s, e)
        init[3:length(init)] <- init[3:length(init)] + 0.2
      }
      legendN[length(legendN) + 1] <- "long_terminal_repeat"
      color[length(color) + 1] <- "darkgreen"
    }
      
    if (input$centromere){
      chrData3 <- subset(chrData, V3 == "centromere")
      if (nrow(chrData3) != 0){
        plotIndividual(chrData3,
                       c("grey"), "grey", start = init[3], fin = init[3] + 0.2,s,e)
        init[4:length(init)] <- init[4:length(init)] + 0.2
      }
      legendN[length(legendN) + 1] <- "centromere"
      color[length(color) + 1] <- "grey"
    }
    
    if (input$rr){
      chrData4 <- subset(chrData, V3 == "repeat_region")
      if (nrow(chrData4) != 0){
        plotIndividual(chrData4,
                       c("magenta"), "magenta", start = init[4], fin = init[4] + 0.2,s,e)
        init[5:length(init)] <- init[5:length(init)] + 0.2
      }
      legendN[length(legendN) + 1] <- "repeat_region"
      color[length(color) + 1] <- "magenta"
    }
      
    if (input$rtt){
      chrData5 <- subset(chrData, V3 == "retrotransposon")
      if (nrow(chrData5) != 0){
        plotIndividual(chrData5,
                       c("orange"), "orange", start = init[5], fin = init[5] + 0.2,s,e)
        init[6:length(init)] <- init[6:length(init)] + 0.2
      }
      legendN[length(legendN) + 1] <- "retrotransposon"
      color[length(color) + 1] <- "orange"
    }
      
    if (input$rRNA){
      chrData6 <- subset(chrData, V3 == "rRNA")
      if (nrow(chrData6) != 0){
        plotIndividual(chrData6, 
                       c("pink"), "pink", start = init[6], fin = init[6] + 0.2,s,e)
        init[7:length(init)] <- init[7:length(init)] + 0.2
      }
      legendN[length(legendN) + 1] <- "rRNA"
      color[length(color) + 1] <- "pink"
    }
      
    if (input$telomere){
      chrData7 <- subset(chrData, V3 == "telomere")
      if (nrow(chrData7) != 0){
        plotIndividual(chrData7,
                       c("green"), "green", start = init[7], fin = init[7] + 0.2,s,e )
        init[8:length(init)] <- init[8:length(init)] + 0.2
      }
      legendN[length(legendN) + 1] <- "telomere"
      color[length(color) + 1] <- "green"
    }
      
    if (input$teg){
      chrData8 <- subset(chrData, V3 == "transposable_element_gene")
      if (nrow(chrData8) != 0){
        plotIndividual(chrData8,
                       c("darkblue"), "darkblue", start = init[8], fin = init[8] +0.2,s,e)
        init[length(init)] <- init[length(init)] + 0.2
      }
      legendN[length(legendN) + 1] <- "transposable_element_gene"
      color[length(color) + 1] <- "darkblue"
    }
      
    if (input$tRNA){
      chrData9 <- subset(chrData, V3 == "tRNA")
      if (nrow(chrData9) != 0){
        plotIndividual(chrData9, 
                       c("purple"), "purple", start = init9, fin = init9+ 0.2,s,e)
      }
      legendN[length(legendN) + 1] <- "tRNA"
      color[length(color) + 1] <- "purple"
    }
      
    lines(c(s,e),c(init[9]+0.1+offset,init[9]+0.1+offset))
    lines(c(e,e),c(0,init[9]+0.1+offset))
    lines(c(s,s),c(0,init[9]+0.1+offset))
    
    legend(e,init[9]+0.1+offset,
           legend=legendN,
           col=color,
           pch=19,border="white",lwd=0, box.col=NA)
   }

  
  output$enrichmentPlot <- renderPlot({
    plot2(input$start, input$end)
    plot3(input$start, input$end)
  })
  
  
  
  output$save <- downloadHandler(
    filename = function() { paste(input$inputfile, '.pdf', sep='')},
    content = function(file) {
      pdf(file,width=15,height=8)
      plot2(input$start, input$end)
      plot3(input$start, input$end)
      dev.off()
    }
  )
}


shinyApp(ui = ui, server = server)