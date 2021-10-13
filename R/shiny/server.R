shinyServer(function(input, output, session){
  observeEvent(input$sendToJs, {
    if (input$how_input == 1) {
      if (input$axis == 1) {
        data_matrix <- array(0, dim = c(input$x_len, input$y_len, input$z_len))
        data_matrix <- aperm(data_matrix, perm=c(2, 3, 1))
        mask_dim <- dim(data_matrix)
        output$h_axis = renderText("Horizontal axis: y")
        output$v_axis = renderText("Vertical axis: z")
        session$sendCustomMessage("initData",
          message = list(dataMatrix = data_matrix, maskDim = mask_dim, ratio = c(input$y_ratio, input$z_ratio), alongAxis = "x"))
      } else if (input$axis == 2) {
        data_matrix <- array(0, dim = c(input$x_len, input$y_len, input$z_len))
        data_matrix <- aperm(data_matrix, perm=c(1, 3, 2))
        mask_dim <- dim(data_matrix)
        output$h_axis = renderText("Horizontal axis: x")
        output$v_axis = renderText("Vertical axis: z")
        session$sendCustomMessage("initData",
          message = list(dataMatrix = data_matrix, maskDim = mask_dim, ratio = c(input$x_ratio, input$z_ratio), alongAxis = "y"))
      } else {
        data_matrix <- array(0, dim = c(input$x_len, input$y_len, input$z_len))
        mask_dim <- dim(data_matrix)
        output$h_axis = renderText("Horizontal axis: z")
        output$v_axis = renderText("Vertical axis: y")
        session$sendCustomMessage("initData",
          message = list(dataMatrix = data_matrix, maskDim = mask_dim, ratio = c(input$x_ratio, input$y_ratio), alongAxis = "z"))
      }
    } else {
      inFile <- input$mask_file
      file <- inFile$datapath
      e <- new.env()
      name <- load(file, envir = e)
      data_matrix <- e[[name]]
      if (input$axis == 1) {
        data_matrix <- aperm(data_matrix, perm=c(2, 3, 1))
        mask_dim <- dim(data_matrix)
        output$h_axis = renderText("Horizontal axis: y")
        output$v_axis = renderText("Vertical axis: z")
        session$sendCustomMessage("existingData",
          message = list(dataMatrix = data_matrix, ratio = c(input$y_ratio, input$z_ratio), alongAxis = "x"))
      } else if (input$axis == 2) {
        data_matrix <- aperm(data_matrix, perm=c(1, 3, 2))
        mask_dim <- dim(data_matrix)
        output$h_axis = renderText("Horizontal axis: x")
        output$v_axis = renderText("Vertical axis: z")
        session$sendCustomMessage("existingData",
          message = list(dataMatrix = data_matrix, ratio = c(input$x_ratio, input$z_ratio), alongAxis = "x"))
      } else {
        mask_dim <- dim(data_matrix)
        output$h_axis = renderText("Horizontal axis: x")
        output$v_axis = renderText("Vertical axis: y")
        session$sendCustomMessage("existingData",
          message = list(dataMatrix = data_matrix, ratio = c(input$x_ratio, input$y_ratio), alongAxis = "x"))
      }
    }
  })

  observeEvent(input$hogefuga, {print(input$hogefuga)})

  output$downloadData = downloadHandler(
    filename = "mask.rda",
    content = function(file) {
      mask_matrix <- input$mask_matrix %>%
        array(dim = input$mask_dim) %>%
        aperm(perm = c(3,2,1))
      if(input$along_axis == "x"){
        mask_matrix <- aperm(mask_matrix, perm = c(3, 1, 2))
      } else if(input$along_axis == "y") {
        mask_matrix <- aperm(mask_matrix, perm = c(1, 3, 2))
      }
      save(mask_matrix, file = file)
    }
  )
})