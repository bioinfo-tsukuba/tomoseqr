#' Output the reconstructed expression pattern as an image.
#' @param tomoObj tomoSeq object
#' @param geneID geneID as string
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
#' @importFrom grDevices
#'    png
#'    dev.off
#' @return NA
#' @examples
#' if (interactive()) {
#'     data(tomoObj)
#'     imageViewer(tomoObj, "gene2")
#' }
#' @export
imageViewer <- function (tomoObj, geneID) {
    checkParameters(tomoObj, geneID)
    recResult <- tomoObj[["results"]][[geneID]][["reconst"]]
    recAlongX <- aperm(recResult, perm = c(2, 3, 1))
    recAlongY <- aperm(recResult, perm = c(1, 3, 2))
    recAlongZ <- aperm(recResult, perm = c(1, 2, 3))

    maskResult <- tomoObj$mask
    maskAlongX <- aperm(maskResult, perm = c(2, 3, 1))
    maskAlongY <- aperm(maskResult, perm = c(1, 3, 2))
    maskAlongZ <- aperm(maskResult, perm = c(1, 2, 3))

    sliderMax <- ceiling(max(recResult) * 1.1)
    ui <- function () {
        fluidPage(
            titlePanel("Image viewer"),
            tabsetPanel(
                type="tabs",

                ## Along x
                tabPanel(
                    "Along x",
                    hr(),
                    sidebarLayout(
                        sidebarPanel(
                            radioButtons(
                                "x_type",
                                label = h3("Type"),
                                choices = list(
                                    "Expression" = "expression",
                                    "Mask" = "mask",
                                    "Unite" = "unite"
                                ), 
                                selected = "expression"
                            ),
                            sliderInput(
                                "xSec",
                                "Section",
                                min=1,
                                max=dim(recResult)[3],
                                step=1,
                                value=1,
                                animate = animationOptions(
                                    interval=100,
                                    loop=TRUE
                                )
                            ),
                            h4("Aspect ratio"),
                            inputPanel(
                                numericInput(
                                    "x_asp_x",
                                    label = "x",
                                    value = dim(recAlongX)[1],
                                    min=1,
                                    width = 100
                                ),
                                numericInput(
                                    "x_asp_y",
                                    value = dim(recAlongX)[2],
                                    label="y",
                                    min=1,
                                    width=100
                                )
                            ),
                            h4("Labels"),
                            inputPanel(
                                textInput(
                                    "x_main",
                                    label="Title",
                                    value=geneID
                                ),
                                textInput(
                                    "x_xlab",
                                    label="x label",
                                    value="y"
                                ),
                                textInput(
                                    "x_ylab",
                                    label="y label",
                                    value="z"
                                ),
                            ),
                            h4("zlim"),
                            inputPanel(
                                sliderInput(
                                    "x_zlim",
                                    "zlim",
                                    min=0,
                                    max=sliderMax,
                                    step=1,
                                    value=c(1, sliderMax),
                                )
                            ),
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
                            plotOutput("xplot")
                        )
                    )
                ),

                ## Along y
                tabPanel(
                    "Along y",
                    hr(),
                    sidebarLayout(
                        sidebarPanel(
                            radioButtons(
                                "y_type",
                                label = h3("Type"),
                                choices = list(
                                    "Expression" = "expression",
                                    "Mask" = "mask",
                                    "Unite" = "unite"
                                ), 
                                selected = "expression"
                            ),
                            sliderInput(
                                "ySec",
                                "Section",
                                min=1,
                                max=dim(recResult)[3],
                                step=1,
                                value=1,
                                animate = animationOptions(
                                    interval=100,
                                    loop=TRUE
                                )
                            ),
                            h4("Aspect ratio"),
                            inputPanel(
                                numericInput(
                                    "y_asp_x",
                                    label = "x",
                                    value = dim(recAlongY)[1],
                                    min=1,
                                    width = 100
                                ),
                                numericInput(
                                    "y_asp_y",
                                    value = dim(recAlongY)[2],
                                    label="y",
                                    min=1,
                                    width=100
                                )
                            ),
                            h4("Labels"),
                            inputPanel(
                                textInput(
                                    "y_main",
                                    label="Title",
                                    value=geneID
                                ),
                                textInput(
                                    "y_xlab",
                                    label="x label",
                                    value="x"
                                ),
                                textInput(
                                    "y_ylab",
                                    label="y label",
                                    value="z"
                                ),
                            ),
                            h4("zlim"),
                            inputPanel(
                                sliderInput(
                                    "y_zlim",
                                    "zlim",
                                    min=0,
                                    max=sliderMax,
                                    step=1,
                                    value=c(1, sliderMax),
                                )
                            ),
                            downloadButton(
                                "YdownloadFig",
                                label = "Download figure"
                            ),
                            downloadButton(
                                "YgenerateGIF",
                                label = "Download GIF"
                            )
                        ),
                        mainPanel(
                            plotOutput("yplot")
                        )
                    )
                ),

                ## Along z
                tabPanel(
                    "Along z",
                    hr(),
                    sidebarLayout(
                        sidebarPanel(
                            radioButtons(
                                "z_type",
                                label = h3("Type"),
                                choices = list(
                                    "Expression" = "expression",
                                    "Mask" = "mask",
                                    "Unite" = "unite"
                                ), 
                                selected = "expression"
                            ),
                            sliderInput(
                                "zSec",
                                "Section",
                                min=1,
                                max=dim(recResult)[3],
                                step=1,
                                value=1,
                                animate = animationOptions(
                                    interval=100,
                                    loop=TRUE
                                )
                            ),
                            h4("Aspect ratio"),
                            inputPanel(
                                numericInput(
                                    "z_asp_x",
                                    label = "x",
                                    value = dim(recAlongZ)[1],
                                    min=1,
                                    width = 100
                                ),
                                numericInput(
                                    "z_asp_y",
                                    value = dim(recAlongY)[2],
                                    label="y",
                                    min=1,
                                    width=100
                                )
                            ),
                            h4("Labels"),
                            inputPanel(
                                textInput(
                                    "z_main",
                                    label="Title",
                                    value=geneID
                                ),
                                textInput(
                                    "z_xlab",
                                    label="x label",
                                    value="x"
                                ),
                                textInput(
                                    "z_ylab",
                                    label="y label",
                                    value="y"
                                ),
                            ),
                            h4("zlim"),
                            inputPanel(
                                sliderInput(
                                    "z_zlim",
                                    "zlim",
                                    min=0,
                                    max=sliderMax,
                                    step=1,
                                    value=c(1, sliderMax),
                                )
                            ),
                            downloadButton(
                                "ZdownloadFig",
                                label = "Download figure"
                            ),
                            downloadButton(
                                "ZgenerateGIF",
                                label = "Download GIF"
                            )
                        ),
                        mainPanel(
                            plotOutput("zplot")
                        )
                    )
                )
            )
        )
    }
    server <- function (input, output) {

        ## Along x
        xPlotArray <- reactive({
            makePlotArray(
                reconstArray=recAlongX,
                maskArray=maskAlongX,
                zlim=input$x_zlim,
                type=input$x_type
            )
        })
        output$xplot <- renderPlot({
            basePlot(
                sourceArray=xPlotArray(),
                sectionNumber=input$xSec,
                main=input$x_main,
                xlab=input$x_xlab,
                ylab=input$x_ylab,
                aspectRatio=(input$x_asp_y / input$x_asp_x)
            )
        })
        output$XdownloadFig <- downloadHandler(
            filename = "xPlot.png",
            content = function (file) {
                png(file, height=1000, width=1000, res=144)
                basePlot(
                    sourceArray=xPlotArray(),
                    sectionNumber=input$xSec,
                    main=input$x_main,
                    xlab=input$x_xlab,
                    ylab=input$x_ylab,
                    aspectRatio=(input$x_asp_y / input$x_asp_x)
                )
                dev.off()
            }
        )
        output$XgenerateGIF <- downloadHandler(
            filename = "xPlot.gif",
            content = function (file) {
                shinyAnimate2d(
                    tomoObj,
                    geneID,
                    target=input$x_type,
                    xaxis=2,
                    yaxis=3,
                    main=input$x_main,
                    xlab=input$x_xlab,
                    ylab=input$x_ylab,
                    file=file,
                    zlim=input$x_zlim,
                    interval=0.1,
                    aspectRatio=c(input$x_asp_x, input$x_asp_y)
                )
            }
        )

        ## Along y
        yPlotArray <- reactive({
            makePlotArray(
                reconstArray=recAlongY,
                maskArray=maskAlongY,
                zlim=input$y_zlim,
                type=input$y_type
            )
        })
        output$yplot <- renderPlot({
            basePlot(
                sourceArray=yPlotArray(),
                sectionNumber=input$ySec,
                main=input$y_main,
                xlab=input$y_xlab,
                ylab=input$y_ylab,
                aspectRatio=(input$y_asp_y / input$y_asp_x)
            )
        })
        output$YdownloadFig <- downloadHandler(
            filename = "yPlot.png",
            content = function (file) {
                png(file, height=1000, width=1000, res=144)
                basePlot(
                    sourceArray=yPlotArray(),
                    sectionNumber=input$ySec,
                    main=input$y_main,
                    xlab=input$y_xlab,
                    ylab=input$y_ylab,
                    aspectRatio=(input$y_asp_y / input$y_asp_x)
                )
                dev.off()
            }
        )

        output$YgenerateGIF <- downloadHandler(
            filename = "yPlot.gif",
            content = function (file) {
                shinyAnimate2d(
                    tomoObj,
                    geneID,
                    target=input$y_type,
                    xaxis=1,
                    yaxis=3,
                    main=input$y_main,
                    xlab=input$y_xlab,
                    ylab=input$y_ylab,
                    file=file,
                    zlim=input$y_zlim,
                    interval=0.1,
                    aspectRatio=c(input$y_asp_x, input$y_asp_y)
                )
            }
        )

        ## Along z
        zPlotArray <- reactive({
            makePlotArray(
                reconstArray=recAlongZ,
                maskArray=maskAlongZ,
                zlim=input$z_zlim,
                type=input$z_type
            )
        })
        output$zplot <- renderPlot({
            basePlot(
                sourceArray=zPlotArray(),
                sectionNumber=input$zSec,
                main=input$z_main,
                xlab=input$z_xlab,
                ylab=input$z_ylab,
                aspectRatio=(input$z_asp_y / input$z_asp_x)
            )
        })
        output$ZdownloadFig <- downloadHandler(
            filename = "zPlot.png",
            content = function (file) {
                png(file, height=1000, width=1000, res=144)
                basePlot(
                    sourceArray=zPlotArray(),
                    sectionNumber=input$zSec,
                    main=input$z_main,
                    xlab=input$z_xlab,
                    ylab=input$z_ylab,
                    aspectRatio=(input$z_asp_y / input$z_asp_x)
                )
                dev.off()
            }
        )

        output$ZgenerateGIF <- downloadHandler(
            filename = "zPlot.gif",
            content = function (file) {
                shinyAnimate2d(
                    tomoObj,
                    geneID,
                    target=input$z_type,
                    xaxis=1,
                    yaxis=2,
                    main=input$z_main,
                    xlab=input$z_xlab,
                    ylab=input$z_ylab,
                    file=file,
                    zlim=input$z_zlim,
                    interval=0.1,
                    aspectRatio=c(input$z_asp_x, input$z_asp_y)
                )
            }
        )
    }
    shinyApp(ui, server)
}
