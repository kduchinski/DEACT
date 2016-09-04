# server.R definition

server <- function(input, output){
 output$Plot <- renderPlotly({
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    df = as.data.frame(read.csv(inFile$datapath, sep='\t', quote=NULL))
    df$id = c(1:nrow(df))
    df.loss = df[df$lossOrGain == "Loss",]
    df.gain = df[df$lossOrGain == "Gain",]
    df.both = df[df$lossOrGain == "Both",]
    if(input$ChangeChoice == "FPKM Difference"){
      fpkmChangeLoss = as.character(df.both$fpkmdiff.loss)
      fpkmChangeGain = as.character(df.both$fpkmdiff.gain)
    } else{
      fpkmChangeLoss = as.character(df.both$foldChange.loss)
      fpkmChangeGain = as.character(df.both$foldChange.gain)
    }
    assign("df", df, envir = .GlobalEnv)
    assign("df.loss", df.loss, envir = .GlobalEnv)
    assign("df.gain", df.gain, envir = .GlobalEnv)
    assign("df.both", df.both, envir = .GlobalEnv)
    geneNames = as.character(df.both$gene_short_name)
    lossOrGain = as.character(df.both$lossOrGain)
    plot_ly(df.both, 
            x = fpkmChangeLoss, 
            y = fpkmChangeGain, 
            text = geneNames, 
            mode = "markers", 
            source = "subse",
            group = lossOrGain) %>%
      layout(dragmode = "select")
  })
  output$StatsTable = renderDataTable({
    stats.names = c("Gain-of-Function", "Loss-of-Function", "Both Loss and Gain", "Upregulated in Both", "Downregulated in Both", "Up Loss, Down Gain", "Down Loss, Up Gain")
    if(is.null(input$file)){
      stats = c(0) * length(stats.names)
    } else{
      df.both = df[df$lossOrGain == "Both",]
      posBoth = df.both[df.both$fpkmdiff.loss > 0 & df.both$fpkmdiff.gain > 0,]
      negBoth = df.both[df.both$fpkmdiff.loss < 0 & df.both$fpkmdiff.gain < 0,]
      posLossNegGain = df.both[df.both$fpkmdiff.loss > 0 & df.both$fpkmdiff.gain < 0,]
      negLossPosGain = df.both[df.both$fpkmdiff.loss < 0 & df.both$fpkmdiff.gain > 0,]
      stats = c(nrow(df.gain), nrow(df.loss), nrow(df.both), nrow(posBoth), nrow(negBoth), nrow(posLossNegGain), nrow(negLossPosGain))
    }
    statstable = data.frame("Category" = stats.names, "Genes" = stats)
  }, options = c(paging = F, searching = F))
  output$Plot1 <- renderPlotly({
    if (is.null(input$file))
      return(NULL)
    if(input$ChangeChoice == "FPKM Difference"){
      fpkmChangeLoss = as.character(df$fpkmdiff.loss)
      fpkmChangeGain = as.character(df$fpkmdiff.gain)
    } else{
      fpkmChangeLoss = as.character(df$foldChange.loss)
      fpkmChangeGain = as.character(df$foldChange.gain)
    }
    geneNames = as.character(df$gene_short_name)
    lossOrGain = as.character(df$lossOrGain)
    plot_ly(df, 
            x = fpkmChangeLoss, 
            y = fpkmChangeGain, 
            text = geneNames, 
            mode = "markers", 
            source = "subset",
            group = lossOrGain) %>%
      layout(dragmode = "select")
  })
  output$DataTable = renderDataTable({
    if (is.null(input$file))
      return(NULL)
    selection = event_data("plotly_selected", source = "subse")
    quads.df = data.frame("fpkmdiff.loss" = numeric(), "fpkmdiff.gain"=numeric(), "foldChange.loss" = numeric(), "foldChange.gain" = numeric(), "gene_short_name"=character(), "qValueLoss" = numeric(), "qValueGain" = numeric(), "lossOrGain"=character(), "id" = integer())
    if(is.null(selection)){
      selection.df = data.frame("fpkmdiff.loss" = numeric(), "fpkmdiff.gain"=numeric(), "foldChange.loss" = numeric(), "foldChange.gain" = numeric(), "gene_short_name"=character(), "qValueLoss" = numeric(), "qValueGain" = numeric(), "lossOrGain"=character(), "id" = integer())
    } else{
      loss = subset(selection, curveNumber == 0)$pointNumber + 1
      gain = subset(selection, curveNumber == 1)$pointNumber + 1
      both = subset(selection, curveNumber == 2)$pointNumber + 1
      loss.df = df.loss[loss,]
      gain.df = df.gain[gain,]
      both.df = df.both[both,]
      selection.df = rbind(loss.df, gain.df, both.df)
    }
    if(input$quad1Button != 0){
      quads.df = rbind(quads.df, posBoth)
    }
    if(input$quad2Button != 0){
      quads.df = rbind(quads.df, negLossPosGain)
    }
    if(input$quad3Button != 0){
      quads.df = rbind(quads.df, negBoth)
    }
    if(input$quad4Button != 0){
      quads.df = rbind(quads.df, posLossNegGain)
    }
    if(input$addButton != 0){
      quads.df = quads.df[!quads.df$id %in% selection.df$id,]
      newpoints.df = newpoints.df[!newpoints.df$id %in% quads.df$id & !newpoints.df$id %in% selection.df$id,]
      selection.df = rbind(selection.df, quads.df, newpoints.df)
    } else{
      quads.df = quads.df[!quads.df$id %in% selection.df$id,]
      selection.df = rbind(selection.df, quads.df)
    } 
    if(input$ChangeChoice == "log2 Fold Change"){
      if(input$showInf == 'All'){
        df = df
      } else{
        df = df.both
      }
      selection.df = selection.df[!(selection.df$foldChange.gain == Inf | selection.df$foldChange.gain == -Inf | selection.df$foldChange.loss == Inf | selection.df$foldChange.loss == -Inf),]
      if(input$posInfLoss){
        selection.df = rbind(selection.df, df[(!df$id %in% selection.df$id & df$foldChange.loss == Inf),])
      }
      if(input$negInfLoss){
        selection.df = rbind(selection.df, df[(!df$id %in% selection.df$id & df$foldChange.loss == -Inf),])
      }
      if(input$posInfGain){
        selection.df = rbind(selection.df, df[(!df$id %in% selection.df$id & df$foldChange.gain == Inf),])
      }
      if(input$negInfGain){
        selection.df = rbind(selection.df, df[(!df$id %in% selection.df$id & df$foldChange.gain == -Inf),])
      }
    }
    assign("selection.df", selection.df, envir = .GlobalEnv)
    if(input$ChangeChoice == 'FPKM Difference'){
      output.df = selection.df[,c('gene_short_name', 'lossOrGain', 'qValueLoss', 'qValueGain', 'fpkmdiff.loss', 'fpkmdiff.gain')]
    } else{
      output.df = selection.df[,c('gene_short_name', 'lossOrGain', 'qValueLoss', 'qValueGain', 'foldChange.loss', 'foldChange.gain')]
    }
  })
  output$extraPoints = renderDataTable({
    if (is.null(input$file))
      return(NULL)
    newpoints = event_data("plotly_selected", source = "subset")
    if(is.null(newpoints)){
      if(input$ChangeChoice == "FPKM Difference"){
        newpoints.df = data.frame("fpkmdiff.loss" = numeric(), "fpkmdiff.gain"=numeric(),"gene_short_name"=character(), "qValueLoss" = numeric(), "qValueGain" = numeric(), "lossOrGain"=character(), "id" = integer())
      } else{
        newpoints.df = data.frame("foldChange.loss" = numeric(), "foldChange.gain"=numeric(),"gene_short_name"=character(), "qValueLoss" = numeric(), "qValueGain" = numeric(), "lossOrGain"=character(), "id" = integer())
      }
    } else{
      loss = subset(newpoints, curveNumber == 0)$pointNumber + 1
      gain = subset(newpoints, curveNumber == 1)$pointNumber + 1
      both = subset(newpoints, curveNumber == 2)$pointNumber + 1
      loss.df = df.loss[loss,]
      gain.df = df.gain[gain,]
      both.df = df.both[both,]
      newpoints.df = rbind(loss.df, gain.df, both.df)
    }
    assign("newpoints.df", newpoints.df, envir=.GlobalEnv)
    if(input$ChangeChoice == 'FPKM Difference'){
      output.df = newpoints.df[,c('gene_short_name', 'lossOrGain', 'qValueLoss', 'qValueGain', 'fpkmdiff.loss', 'fpkmdiff.gain')]
    } else{
      newpoints.df = newpoints.df[!(newpoints.df$foldChange.gain == Inf | newpoints.df$foldChange.gain == -Inf | newpoints.df$foldChange.loss == Inf | newpoints.df$foldChange.loss == -Inf),]
      output.df = newpoints.df[,c('gene_short_name', 'lossOrGain', 'qValueLoss', 'qValueGain', 'foldChange.loss', 'foldChange.gain')]
    }
  })
  output$downloadData <- downloadHandler(
    "selection",
    content = function(file) {
      write.table(selection.df[,c("gene_short_name","foldChange.loss", "foldChange.gain", "qValueLoss", "qValueGain")], file, row.names = F)
    },
    contentType = "text/txt"
  )
}