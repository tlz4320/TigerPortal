server = function(input, output, session) {
  observeEvent(input$doPredict, {
    output$button <- renderUI({
      
    })
    sgInfos <- session$userData$infos
    sel <- unlist(lapply(1 : nrow(sgInfos), function(i){
      input[[paste0("id_", i)]]
    }))
    sgInfos$id <- paste0("id_", 1 : nrow(sgInfos))
    sgInfos_sel <- sgInfos[sel,]
    for_indelphi <- lapply(1 : nrow(sgInfos_sel), function(i){
      
      strand <- sgInfos_sel$strand[i]
      start <- sgInfos_sel$start[i]
      if(strand == "+"){
        seq <- session$userData$seq
      } else {
        seq <- session$userData$reverseSeq
      }
      data.frame(id = sgInfos_sel$id[i],
                 seq = seq, pos = start - 4)
      
    })
    for_indelphi <- data.frame(do.call(rbind, for_indelphi))
    output_file <- tempfile(fileext = ".csv")
    td <- tempfile(pattern = "Result")
    dir.create(td, recursive = T)
    write.table(for_indelphi, file = output_file, col.names=T, row.names=F, quote = T, sep=",")
    system(paste("/home/bioinfo/micromamba/envs/bio/bin/python /home/bioinfo/code/inDelphi-model/script4R.py",
                 input$database, output_file, td))
    session$userData$tmpfiles <- c(td, output_file)
    session$userData$download <- paste0(td, ".tar.gz")
    readPredict(td, sgInfos_sel, session$userData$infos2, session$userData$cds_region)
    system(paste0("sh ~/code/earProject/gene_therapy/Model/tar.sh ", tempdir(), " ", basename(td)))
    output$button <- renderUI({
      tagList(
        actionButton("showRes", "Show Result"),
        downloadButton("downloadData", label = "Download")
      )
    })
    removeUI("#twarn")
    renderResult(output, td)
  })
  output$downloadData <- downloadHandler(
    filename <- function() {
      paste("output", "zip", sep=".")
    },
    
    content <- function(file) {
      print(session$userData$download)
      file.copy(session$userData$download, file)
    },
    contentType = "application/zip"
  )
  observeEvent(input$sel, {
    sgInfos <- session$userData$infos
    for(i in 1 : nrow(sgInfos)){
      updateCheckboxInput(session, paste0("id_", i), value = T)
    }
    output$selAll <- renderUI({
      actionButton("selNo", "Unselect All")
    })
  })
  observeEvent(input$selNo, {
    sgInfos <- session$userData$infos
    for(i in 1 : nrow(sgInfos)){
      updateCheckboxInput(session, paste0("id_", i), value = F)
    }
    output$selAll <- renderUI({
      actionButton("sel", "Select All")
    })
  })
  observeEvent(input$showRes, {
    updateNavbarPage(session=session,
                     inputId="pages",
                     selected="Result")
  })
  
  observeEvent(input$searchSg, {
    checked <- checkRS(input$seq_input)
    if(class(checked) == "logical"){
      if(!checked){
        shinyalert("Not a valid RS ID", "Not a valid RS ID", type="error")
        return()
      }
    }
    mutationCDS <- checked
    mutationCDS <- mutationCDS[order(str_length(mutationCDS$V7), decreasing = T),]
    inputseq <- mutationCDS$V7[1]
    cds_region <- c(mutationCDS$V8[1], mutationCDS$V9[1])
    session$userData$cds_region <- cds_region
    checked <- formatCheck(inputseq)
    if(class(checked) == "logical"){
      if(!checked){
        shinyalert("Not a valid input", "Not a valid input", type="error")
        return()
      }
    }
    indel_pos <- checked$pos
    seq <- checked$seq
    isIndel <- checked$isIndel
    indelNT <- checked$indelNT
    possiblePAM <- findPAM(seq, pam = input$ngg_seq)
    if(input$reverse){
      reverseSeq <- as.character(reverseComplement(DNAString(seq)))
      possiblePAM2 <- findPAM(reverseSeq, pam = input$ngg_seq)
    } else {
      #只是为了方便下面判断而已
      possiblePAM2 <- possiblePAM
    }
    if(nrow(possiblePAM) == 0 & nrow(possiblePAM2) == 0){
      shinyalert("No any valid gRNA", "No any valid gRNA", type="error")
      return()
    }
    if(nrow(possiblePAM) != 0){
      possiblePAM$dis <- abs(possiblePAM$start - 4 - indel_pos)
      coloredSeq_f <- colorSeq(seq, possiblePAM, indel_pos, isIndel,indelNT)
      coloredSeq_f$strand <- "+"
    } else {
      coloredSeq_f <- data.frame()
    }
    if(input$reverse & nrow(possiblePAM2) != 0){
      possiblePAM2$dis <- abs(possiblePAM2$start - 4 - (str_length(seq) - indel_pos + 1))
      coloredSeq_r <- colorSeq(reverseSeq, possiblePAM2, str_length(reverseSeq) - indel_pos + 1, isIndel,indelNT)
      coloredSeq_r$strand <- "-"
      session$userData$reverseSeq <- reverseSeq
    } else {
      coloredSeq_r <- data.frame()
    }
    
    if(nrow(coloredSeq_r) != 0 & nrow(coloredSeq_f) != 0){
      coloredSeq <- data.frame(rbind(coloredSeq_f, coloredSeq_r))
    } else if(nrow(coloredSeq_r) == 0){
      coloredSeq <- coloredSeq_f
    } else {
      coloredSeq <- coloredSeq_r
    }
    if(isIndel != 0){
      coloredSeq <- coloredSeq[order(coloredSeq$dis),]
      coloredSeq <- coloredSeq[coloredSeq$dis < 20,]
    }
    
    if(nrow(possiblePAM) == 0 & nrow(possiblePAM2) == 0){
      shinyalert("No any valid gRNA", "No any valid gRNA", type="error")
      return()
    }
    checkBoxs <- unlist(lapply(1 : nrow(coloredSeq), function(i){
      as.character(checkboxInput(paste0("id_", i),NULL, value = i <= 3, width = "10%"))
    }))
    df2 <- data.frame(Seq = coloredSeq$colored, 
                      Strand = coloredSeq$strand, 
                      choose=checkBoxs)
    output$tableOutput <- renderTable(df2, 
                                      sanitize.text.function = function(x) x, 
                                      bordered = T, 
                                      align = "c", 
                                      striped = TRUE)
    session$userData$infos <- coloredSeq
    session$userData$seq <- seq
    session$userData$infos2 <- checked
    output$button <- renderUI({
      actionButton("doPredict", "Do Predict")
    })
    output$selAll <- renderUI({
      actionButton("sel", "Select All")
    })
  })
  session$onSessionEnded(function() {
    unlink(session$userData$tmpfiles, recursive = T)
    unlink(session$userData$download)
  }) 
}