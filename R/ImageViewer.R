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
                        sliderInput(
                            "xSec",
                            "Section",
                            min=1,
                            max=dim(recResult)[1],
                            step=1,
                            value=1
                        )
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
    server <- function (input, output) {}
    shinyApp(ui, server)
}