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
imageViewer <- function (tomoObj) {
    ui <- function () {
        fluidPage(
            titlePanel("Image viewer"),
            selectInput(
                "geneID",
                label = h3("Select gene ID"),
                choices = names(tomoObj[["results"]])
            ),
            tabsetPanel(
                type="tabs",

                tabPanel(
                    "2D view",
                    hr(),
                    selectInput(
                        "angle",
                        label = h3("Select Angle"),
                        choices = c("x", "y", "z")
                    ),
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
                            uiOutput("sectionSlider"),

                            h4("Aspect ratio"),
                            inputPanel(
                                uiOutput("asp_x"),
                                uiOutput("asp_y")
                            ),
                            h4("Labels"),
                            inputPanel(
                                uiOutput("title"),
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
                                uiOutput("zlim_slider")
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
                tabPanel(
                    "3D view",
                    sidebarLayout(
                        sidebarPanel(
                            numericInput("threshold_3D", label = h3("Threshold"), value = 1),
                            numericInput("exp_size_3D", label = h3("size of expression dots"), value = 80),
                            numericInput("exp_opacity_3D", label = h3("Opacity of expression dots"), value = 100),
                            numericInput("mask_size_3D", label = h3("size of mask dots"), value = 2),
                            numericInput("mask_opacity_3D", label = h3("Opacity of mask dots"), value = 70),
                            textInput("mask_color_3D", label = h3("color of mask dots"), value = "#999999"),
                            numericInput("aspX_3D", "Width of section: x", value = 1),
                            numericInput("aspY_3D", "Width of section: y", value = 1),
                            numericInput("aspZ_3D", "Width of section: z", value = 1),
                            textInput("xlab_3D", label = "x label", value = "x"),
                            textInput("ylab_3D", label = "y label", value = "y"),
                            textInput("zlab_3D", label = "z label", value = "z"),
                            checkboxInput("addMask_3D", label = "Add Mask", value = TRUE)
                        ),
                        mainPanel(
                            plotlyOutput("choice_3D", width = "1000px", height = "1000px"),
                        )
                    )
                )
            )
        )
    }
    server <- function (input, output) {
    maskResult <- tomoObj$mask
    maskAlong <- list()
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
        max <- dim(recAlong()[[input$angle]])[3]
        sliderInput(
            "xSec",
            "Section",
            min=1,
            max=max,
            step=1,
            value=1,
            animate = animationOptions(
                interval=100,
                loop=TRUE
            )
        )
    })
    output$title <- renderUI({
        textInput(
            "x_main",
            label="Title",
            value=input$geneID
        )
    })
    output$asp_x <- renderUI({
        numericInput(
            "x_asp_x",
            label = "x",
            value = dim(recAlong()[[input$angle]])[1],
            min=1,
            width = 100
        )
    })
    output$asp_y <- renderUI({
        numericInput(
            "x_asp_y",
            value = dim(recAlong()[[input$angle]])[2],
            label="y",
            min=1,
            width=100
        )
    })

    output$zlim_slider <- renderUI({
        sliderMax <- ceiling(max(recResult()) * 1.1)
        sliderInput(
            "x_zlim",
            "zlim",
            min=0,
            max=sliderMax,
            step=1,
            value=c(1, sliderMax),
        )
    })

        xPlotArray <- reactive({
            makePlotArray(
                reconstArray=recAlong()[[input$angle]],
                maskArray=maskAlong[[input$angle]],
                zlim=input$x_zlim,
                type=input$x_type
            )
        })
        output$xplot <- renderPlot({
            if (is.null(input$xSec)) {
                message("Please wait...")
            } else {
                basePlot(
                    sourceArray=xPlotArray(),
                    sectionNumber=input$xSec,
                    main=input$x_main,
                    xlab=input$x_xlab,
                    ylab=input$x_ylab,
                    aspectRatio=(input$x_asp_y / input$x_asp_x)
                )
            }
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
                animate2d(
                    tomoObj,
                    input$geneID,
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

        tomoObj$mask <- plotsurface(tomoObj$mask)
        output$choice_3D <- renderPlotly(
            showExp(
                tomoObj,
                geneID = input$geneID,
                threshold = input$threshold_3D,
                xlab = input$xlab_3D,
                ylab = input$ylab_3D,
                zlab = input$zlab_3D,
                addMask = input$addMask_3D,
                exp_size = input$exp_size_3D,
                exp_opacity = input$exp_opacity_3D,
                mask_size = input$mask_size_3D,
                mask_opacity = input$mask_opacity_3D,
                mask_color = input$mask_color_3D,
                aspX = input$aspX_3D,
                aspY = input$aspY_3D,
                aspZ = input$aspZ_3D
            )
        )
    }
    shinyApp(ui, server)
}
