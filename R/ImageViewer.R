#' Output the reconstructed expression pattern as an image.
#' @param tomoObj tomoSeq object
#' @importFrom shiny
#'    reactive
#'    animationOptions
#'    plotOutput
#'    renderPlot
#'    sidebarLayout
#'    sliderInput
#'    tabPanel
#'    tabsetPanel
#'    textInput
#'    hr
#'    h3
#'    column
#'    selectInput
#'    fluidRow
#'    uiOutput
#'    checkboxInput
#'    renderUI
#' @importFrom grDevices
#'    png
#'    dev.off
#' @importFrom plotly
#'    plot_ly
#'    add_markers
#'    layout
#'    config
#'    plotlyOutput
#'    renderPlotly
#' @importFrom dplyr
#'    sym
#' @importFrom ggplot2
#'    geom_tile
#'    aes
#'    ggtitle
#'    ggsave
#'    ggplot
#'    scale_fill_gradientn
#'    ylim
#'    xlim
#'    theme_minimal
#'    theme
#'    element_rect
#'    element_text
#'    element_line
#'    element_blank
#' @return NA
#' @examples
#' if (interactive()) {
#'     data(tomoObj)
#'     imageViewer(tomoObj, "gene2")
#' }
#' @export
imageViewer <- function (tomoObj) {
    ui <- function () {
        fluidPage(
            titlePanel("Image viewer"),
            selectInput(
                "geneID",
                label = h4("Gene ID"),
                choices = names(tomoObj[["results"]])
            ),
            fluidRow(
                column(4,
                    textInput(
                        "xlab",
                        label = "x label",
                        value = "x"
                    )
                ),
                column(4,
                    textInput(
                        "ylab",
                        label = "y label",
                        value = "y"
                    )
                ),
                column(4,
                    textInput(
                        "zlab",
                        label = "z label",
                        value = "z"
                    )
                )
            ),
            h4("Section width (relative value)"),
            fluidRow(
                column(4,
                    numericInput(
                        "xWidth",
                        label = "x",
                        value = 1
                    )
                ),
                column(4,
                    numericInput(
                        "yWidth",
                        label = "y",
                        value = 1
                    )
                ),
                column(4,
                    numericInput(
                        "zWidth",
                        label = "z",
                        value = 1
                    )
                )
            ),
            tabsetPanel(
                type="tabs",

                tabPanel(
                    "2D view",
                    # hr(),
                    sidebarLayout(
                        sidebarPanel(
                            selectInput(
                                "angle",
                                label = h4("Select Angle"),
                                choices = c("x", "y", "z")
                            ),
                            uiOutput("sectionSlider"),
                            uiOutput("title"),
                            uiOutput("zlim_slider"),
                            downloadButton(
                                "XdownloadFig",
                                label = "Download figure"
                            ),
                            downloadButton(
                                "XgenerateGIF",
                                label = "Download GIF"
                            )
                        ),
                        mainPanel(
                            plotOutput("xplot", height = "800px", width="800px")
                        )
                    )
                ),
                tabPanel(
                    "3D view",
                    sidebarLayout(
                        sidebarPanel(
                            numericInput(
                                "threshold_3D",
                                label = h4("Threshold"),
                                value = 1
                            ),
                            numericInput(
                                "exp_size_3D",
                                label = h4("size of expression dots"),
                                value = 40
                            ),
                            numericInput(
                                "exp_opacity_3D",
                                label = h4("Opacity of expression dots"),
                                value = 80
                            ),
                            numericInput(
                                "mask_size_3D",
                                label = h4("size of mask dots"),
                                value = 2
                            ),
                            numericInput(
                                "mask_opacity_3D",
                                label = h4("Opacity of mask dots"),
                                value = 70
                            ),
                            textInput(
                                "mask_color_3D",
                                label = h4("color of mask dots"),
                                value = "#999999"
                            ),
                            checkboxInput(
                                "addMask_3D",
                                label = "Add Mask",
                                value = TRUE
                            )
                        ),
                        mainPanel(
                            plotlyOutput(
                                "choice_3D",
                                height = "800px",
                                width = "800px"
                            )
                        )
                    )
                )
            )
        )
    }
    server <- function (input, output) {
    maskResult <- tomoObj$mask
    axesOrder <- list(
        "x" = c("y", "z", "x"),
        "y" = c("x", "z", "y"),
        "z" = c("x", "y", "z")
    )
    maskAlong <- list()
    maskDim <- dim(maskResult)
    dimOrder <- list(
        "x" = c(maskDim[2], maskDim[3], maskDim[1]),
        "y" = c(maskDim[1], maskDim[3], maskDim[2]),
        "z" = maskDim
    )
    maskAlong[["x"]] <- aperm(maskResult, perm = c(2, 3, 1))
    maskAlong[["y"]] <- aperm(maskResult, perm = c(1, 3, 2))
    maskAlong[["z"]] <- aperm(maskResult, perm = c(1, 2, 3))

    recResult <- reactive({
        tomoObj[["results"]][[input$geneID]][["reconst"]]
    })
    recAlong <- reactive({
        recAlong <- list()
        recAlong[["x"]] <- aperm(recResult(), perm = c(2, 3, 1))
        recAlong[["y"]] <- aperm(recResult(), perm = c(1, 3, 2))
        recAlong[["z"]] <- aperm(recResult(), perm = c(1, 2, 3))
        return(recAlong)
    })
    
    output$sectionSlider <- renderUI({
        max <- dimOrder[[input$angle]][3]
        sliderInput(
            "xSec",
            "Section",
            min=1,
            max=max,
            step=1,
            value=1,
            animate = animationOptions(
                interval=400,
                loop=TRUE
            )
        )
    })
    output$title <- renderUI({
        textInput(
            "main",
            label="Title",
            value=input$geneID
        )
    })

    output$zlim_slider <- renderUI({
        sliderMax <- ceiling(max(recResult()))
        sliderInput(
            "x_zlim",
            "zlim",
            min=0,
            max=sliderMax,
            step=1,
            value=c(0, sliderMax),
        )
    })

        maskDf <- matrixToDataFrame(tomoObj[["mask"]])
        expDf <- reactive({
            (toDataFrame(tomoObj, input$geneID))[maskDf[, 4] == 1,]
        })

        labelOrder <- reactive({
            list(
            "x" = c(input$ylab, input$zlab),
            "y" = c(input$xlab, input$zlab),
            "z" = c(input$xlab, input$ylab)
            )
        })
        widthOrder <- reactive({
            list(
            "x" = c(input$yWidth, input$zWidth),
            "y" = c(input$xWidth, input$zWidth),
            "z" = c(input$xWidth, input$yWidth)
            )
        })
        base_plot <- reactive({
            makeBasePlot(
                expDf = expDf(),
                xAxis = axesOrder[[input$angle]][1],
                yAxis = axesOrder[[input$angle]][2],
                xMax = dimOrder[[input$angle]][1],
                yMax = dimOrder[[input$angle]][2],
                xAsp = dimOrder[[input$angle]][1] *
                    widthOrder()[[input$angle]][1],
                yAsp = dimOrder[[input$angle]][2] *
                    widthOrder()[[input$angle]][2],
                xlabel = labelOrder()[[input$angle]][1],
                ylabel = labelOrder()[[input$angle]][2],
                zlim = input$x_zlim
            )
        })
        plot2D <- reactive({
            if (is.null(input$xSec)) {
                message("Please wait...")
            } else {
                xAxis <- sym(axesOrder[[input$angle]][1])
                yAxis <- sym(axesOrder[[input$angle]][2])
                fill <- sym("value")
                return(
                base_plot() +
                    geom_tile(
                        data=expDf()[expDf()[[input$angle]] == input$xSec,],
                        aes(
                            x = !!xAxis,
                            y = !!yAxis,
                            fill = !!fill
                        )
                    ) +
                    ggtitle(str_c(input$main, " : ", input$xSec))
                )
            }
        })
        
        output$xplot <- renderPlot({
            if (is.null(input$xSec)) {
                message("Please wait...")
            } else {
                print(plot2D())
            }
        })

        output$XdownloadFig <- downloadHandler(
            filename = reactive({str_c(input$geneID, ".png")}),
            content = function (file) {
                ggsave(
                    file,
                    plot = plot2D(),
                    height=800,
                    width=800,
                    units = "px",
                    scale = 5
                )
            }
        )

        output$XgenerateGIF <- downloadHandler(
            filename = "xPlot.gif",
            content = function (file) {
                tomoObj[["mask"]] <- maskResult
                animate2d(
                    tomoObj,
                    input$geneID,
                    along = input$angle,
                    main=input$main,
                    xlab=labelOrder()[[input$angle]][1],
                    ylab=labelOrder()[[input$angle]][2],
                    file=file,
                    zlim=input$x_zlim,
                    interval=0.1,
                    aspectX = widthOrder()[[input$angle]][1],
                    aspectY = widthOrder()[[input$angle]][2]
                )
            }
        )

        tomoObj$mask <- plotsurface(tomoObj$mask)
        output$choice_3D <- renderPlotly(
            showExp(
                tomoObj,
                geneID = input$geneID,
                threshold = input$threshold_3D,
                xlab = input$xlab,
                ylab = input$ylab,
                zlab = input$zlab,
                addMask = input$addMask_3D,
                exp_size = input$exp_size_3D,
                exp_opacity = input$exp_opacity_3D,
                mask_size = input$mask_size_3D,
                mask_opacity = input$mask_opacity_3D,
                mask_color = input$mask_color_3D,
                aspX = input$xWidth,
                aspY = input$yWidth,
                aspZ = input$zWidth
            )
        )
    }
    shinyApp(ui, server)
}
