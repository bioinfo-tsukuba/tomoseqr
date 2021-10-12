library(shiny)
library(jsonlite)
library(dplyr)
shinyUI(fluidPage(
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "custom_style.css")),
  titlePanel("Masker (tomoseqr mask maker)"),
  sidebarPanel(style = "overflow-y:scroll;max-height: 730px",
               h4("Input style"),
               inputPanel(
  radioButtons("how_input", label = "How input", choices = list("New mask" = 1, "From file" = 2), selected = 1),
  fileInput("mask_file", label = "Mask file")
  ),
    conditionalPanel(
      condition = "input.how_input == 1",
      h4("Number of intercept"),
      inputPanel(
        numericInput("x_len", label = "x", width = 100, value = 10),
        numericInput("y_len", label = "y", width = 100, value = 10),
        numericInput("z_len", label = "z", width = 100, value = 10)
      )),
    h4("Width of cells"),
    inputPanel(
      numericInput("x_ratio", label = "x", width = 100, value = 10),
      numericInput("y_ratio", label = "y", width = 100, value = 10),
      numericInput("z_ratio", label = "z", width = 100, value = 10)
    ),
  h4("Along which axis?"),
  inputPanel(radioButtons("axis", choices = list("x" = 1, "y" = 2, "z" = 3), selected = 1, label = "axis")),
  actionButton(inputId = "sendToJs", label = "make table"),
  downloadButton("downloadData", label = "Download as .rda")
  ),
  mainPanel(
    tags$p(
    tags$div(class="table_button",
    tags$button("back", onclick="back_click()"),
    tags$button("next", onclick="next_click()"),
    tags$button("Pencil", onclick="SetPencil()"),
    tags$button("Eraser", onclick="SetEraser()")
    )),
    textOutput("h_axis"),
    textOutput("v_axis"),
    tags$p(
    tags$div(id = "txtY"),
    tags$div(id="nowPage"),
    ),
  tags$canvas(id="sample"),
  tags$div(class="tmp_canvas",
  tags$canvas(id="tmpCanvas")
  ),
  tags$script(src = "first_load.js")
  )
    )
)