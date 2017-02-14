library(shinyjs)
fluidPage(
  useShinyjs(),
  titlePanel("Format SNP Group file from dbSNP ID"),
  tags$code("Author: Zhao Rui"),
  tags$code("Last update: 2017-02-14"),
  tags$link(rel = 'stylesheet', type = 'text/css', href = 'style.css'),
  tags$head(tags$script(src="jquery.highlight-5.js")),
  fluidRow(
    column(3, wellPanel(
      selectInput('format', 'Common SNP Mask:', choices = c("IUPAC", "Lower Case", "n"), selected = "Lower Case"),
      selectInput('dbsnp', 'Common dbSNP DB:', choices = c("147", "144", "141"), selected = "147"),
      sliderInput("FLANK_LEN", "FLANK LENGTH:",
                  min = 50, max = 300, value = 100, step= 10),
      textAreaInput("rstext", "dnSNP rs #", "", resize = "vertical"),
      actionButton("submit", "Submit"),
      actionButton("highlight", "Highlight common SNP"),
      downloadButton('downloadData', 'Download SNP Group File'),
      helpText("将已知dbSNP编号的SNP位点转为MassArray输入需要的SNP Group格式。
               dbSNP编号用换行分隔，需要加rs前缀。
               常见SNP（MAF>1%）的碱基默认标记为小写字母，亦可标记为简并碱基或n。")
    )),
    column(9,
      helpText("IUPAC Mapping Code:"),
      verbatimTextOutput("text"),
      tableOutput("tbl")
    )
  )
)
