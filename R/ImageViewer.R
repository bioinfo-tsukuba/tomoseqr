ImageViewer <- function (tomoObj, geneID) {
    recResult <- tomoObj$GetReconstructedResult(geneID)
    ui <- function () {
        fluidPage(
            titlePanel("Image viewer"),
            mainPanel(
                tabsetPanel(
                    type="tabs",
                    tabPanel(
                        "x",
                        plotOutput("xplot"),
                        sliderInput(
                            "xSec",
                            "Section",
                            min=1,
                            max=dim(recResult)[1],
                            step=1,
                            value=1,
                            # animate=TRUE,
                        animate = animationOptions(
                            interval=70,
                            loop=TRUE
                        )
                        ),
                        downloadButton("downloadFig", label = "Download figure")
                    ),
                    tabPanel(
                        "y",
                        sliderInput(
                            "ySec",
                            "Section",
                            min=1,
                            max=dim(recResult)[2],
                            step=1,
                            value=1
                        )
                    ),
                    tabPanel(
                        "z",
                        sliderInput(
                            "zSec",
                            "Section",
                            min=1,
                            max=dim(recResult)[3],
                            step=1,
                            value=1
                        )
                        )
                )
            )
        )
    }
    server <- function (input, output) {
        output$xplot <- renderPlot({
            filled.contour(recResult[,,input$xSec], zlim=c(min(recResult), max(recResult)))
            # plot(1, input$xSec)
        })
        output$downloadFig <- downloadHandler(
            filename = "xPlot.png",
            content = function (file) {
                png(file, height=1000, width=1000, res=144)
                filled.contour(recResult[,,input$xSec], zlim=c(min(recResult), max(recResult)))
                dev.off()
            }
        )
    }
    shinyApp(ui, server)
}