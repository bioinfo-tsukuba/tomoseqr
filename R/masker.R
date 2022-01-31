#' Make mask
#' @importFrom shiny
#'    fluidPage
#'    tags
#'    titlePanel
#'    sidebarPanel
#'    h4
#'    inputPanel
#'    radioButtons
#'    fileInput
#'    conditionalPanel
#'    numericInput
#'    actionButton
#'    downloadButton
#'    mainPanel
#'    textOutput
#'    observeEvent
#'    renderText
#'    downloadHandler
#'    shinyApp
#'    addResourcePath
#'    removeResourcePath
#' @import jsonlite
#' @return NA
#' @examples
#' if (interactive()) {
#'     Masker()
#' }
#' @export
Masker <- function () {
    ui <- function () {
        addResourcePath(
            'www',
            system.file("shinyApp/www", package = "tomoseqr")
        )
        on.exit(removeResourcePath('www'))
        fluidPage(
            tags$head(
                tags$link(
                    rel = "stylesheet",
                    type = "text/css",
                    href = "www/custom_style.css"
                )
            ),
            titlePanel("Masker (tomoseqr mask maker)"),
            sidebarPanel(style = "overflow-y: scroll; max-height: 730px",
                h4("Input style"),
                inputPanel(
                    radioButtons(
                        "how_input",
                        label = "How input",
                        choices = list("New mask" = 1, "From file" = 2),
                        selected = 1
                    ),
                    fileInput("mask_file", label = "Mask file")
                ),
                conditionalPanel(condition = "input.how_input == 1",
                    h4("Number of intercept"),
                    inputPanel(
                        numericInput(
                            "x_len",
                            label = "x",
                            width = 100,
                            value = 10
                        ),
                        numericInput(
                            "y_len",
                            label = "y",
                            width = 100,
                            value = 10
                        ),
                        numericInput(
                            "z_len",
                            label = "z",
                            width = 100,
                            value = 10
                        )
                    )
                ),
                h4("Width of cells"),
                inputPanel(
                    numericInput(
                        "x_ratio",
                        label = "x",
                        width = 100,
                        value = 10
                    ),
                    numericInput(
                        "y_ratio",
                        label = "y",
                        width = 100,
                        value = 10
                    ),
                    numericInput(
                        "z_ratio",
                        label = "z",
                        width = 100,
                        value = 10
                    )
                ),
                h4("Along which axis?"),
                inputPanel(
                    radioButtons(
                        "axis",
                        choices = list("x" = 1, "y" = 2, "z" = 3),
                        selected = 1,
                        label = "axis"
                    )
                ),
                actionButton(inputId = "sendToJs", label = "make table"),
                downloadButton("downloadData", label = "Download as .rda")
            ),
            mainPanel(
                tags$p(
                    tags$div(class = "table_button",
                        tags$button("back", onclick = "back_click()"),
                        tags$button("next", onclick = "next_click()"),
                        tags$button("Pencil", onclick = "SetPencil()"),
                        tags$button("Eraser", onclick = "SetEraser()"),
                        tags$button("Copy", onclick = "CopyCanvas()"),
                        tags$button("Paste", onclick = "PasteCanvas()")
                    )
                ),
                textOutput("h_axis"),
                textOutput("v_axis"),
                tags$p(
                    tags$div(id = "position"),
                    tags$div(id="nowPage"),
                ),
                tags$canvas(id="sample"),
                tags$div(class="tmp_canvas",
                tags$canvas(id="tmpCanvas")
                ),
                tags$script(src="www/first_load.js")
            )
        )
    }

    server <- function (input, output, session) {
        observeEvent(input$sendToJs, {
            if (input$how_input == 1) {
                x_len <- input$x_len
                y_len <- input$y_len
                z_len <- input$z_len
                x_ratio <- input$x_ratio
                y_ratio <- input$y_ratio
                z_ratio <- input$z_ratio
                data_matrix <- array(0, dim = c(x_len, y_len, z_len))
                if (input$axis == 1) {
                    data_matrix <- aperm(data_matrix, perm=c(2, 3, 1))
                    mask_dim <- dim(data_matrix)
                    output$h_axis = renderText("Horizontal axis: y")
                    output$v_axis = renderText("Vertical axis: z")
                    session$sendCustomMessage(
                        "initData",
                        message = list(
                            dataMatrix = data_matrix,
                            maskDim = mask_dim,
                            ratio = c(y_ratio, z_ratio),
                            alongAxis = "x"
                        )
                    )
                } else if (input$axis == 2) {
                    data_matrix <- aperm(data_matrix, perm=c(1, 3, 2))
                    mask_dim <- dim(data_matrix)
                    output$h_axis = renderText("Horizontal axis: x")
                    output$v_axis = renderText("Vertical axis: z")
                    session$sendCustomMessage(
                        "initData",
                        message = list(
                            dataMatrix = data_matrix,
                            maskDim = mask_dim,
                            ratio = c(x_ratio, z_ratio),
                            alongAxis = "y"
                        )
                    )
                } else {
                    mask_dim <- dim(data_matrix)
                    output$h_axis = renderText("Horizontal axis: z")
                    output$v_axis = renderText("Vertical axis: y")
                    session$sendCustomMessage(
                        "initData",
                        message = list(
                            dataMatrix = data_matrix,
                            maskDim = mask_dim,
                            ratio = c(x_ratio, y_ratio),
                            alongAxis = "z"
                        )
                    )
                }
            } else {
                ## File-loading process
                inFile <- input$mask_file
                file <- inFile$datapath
                e <- new.env()
                name <- load(file, envir = e)
                data_matrix <- e[[name]]

                x_ratio <- input$x_ratio
                y_ratio <- input$y_ratio
                z_ratio <- input$z_ratio

                if (input$axis == 1) {
                    data_matrix <- aperm(data_matrix, perm = c(2, 3, 1))
                    mask_dim <- dim(data_matrix)
                    output$h_axis = renderText("Horizontal axis: y")
                    output$v_axis = renderText("Vertical axis: z")
                    session$sendCustomMessage(
                        "existingData",
                        message = list(
                            dataMatrix = data_matrix,
                            ratio = c(y_ratio, z_ratio),
                            alongAxis = "x"
                        )
                    )
                } else if (input$axis == 2) {
                    data_matrix <- aperm(data_matrix, perm = c(1, 3, 2))
                    mask_dim <- dim(data_matrix)
                    output$h_axis = renderText("Horizontal axis: x")
                    output$v_axis = renderText("Vertical axis: z")
                    session$sendCustomMessage(
                        "existingData",
                        message = list(
                            dataMatrix = data_matrix,
                            ratio = c(x_ratio, z_ratio),
                            alongAxis = "y"
                        )
                    )
                } else {
                    mask_dim <- dim(data_matrix)
                    output$h_axis = renderText("Horizontal axis: x")
                    output$v_axis = renderText("Vertical axis: y")
                    session$sendCustomMessage(
                        "existingData",
                        message = list(
                            dataMatrix = data_matrix,
                            ratio = c(x_ratio, y_ratio),
                            alongAxis = "z"
                        )
                    )
                }
            }
        })

        output$downloadData = downloadHandler(
            filename = "mask.rda",
            content = function(file) {
                mask_matrix <- input$mask_matrix %>%
                    array(dim = rev(input$mask_dim)) %>%
                    aperm(perm = c(3, 2, 1))
                if (input$along_axis == "x") {
                    mask_matrix <- aperm(mask_matrix, perm = c(3, 1, 2))
                } else if (input$along_axis == "y") {
                    mask_matrix <- aperm(mask_matrix, perm = c(1, 3, 2))
                }
                save(mask_matrix, file = file)
            }
        )
    }

    shinyApp(ui, server)
}