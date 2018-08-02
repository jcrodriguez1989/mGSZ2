library('DT');
library('plotly');

resultsUi <- function() {
    tabPanel(title='mGSZ2 results', value='mGsz2Restab',
             br(),
             wellPanel(fluidRow(
                 column(2, actionButton(inputId='saveResbtn', label='Save results')),
                 column(4, fileInput(inputId='loadResFile', label=NULL, buttonLabel='Load results', accept='.RData'))
             )),
             wellPanel(fluidRow(DTOutput(outputId='mGsz2ResDtable'))),
             wellPanel(
                 fluidRow(actionButton(inputId='enrScorePlotbtn', label='Enrichment Plot')),
                 br(),
                 fluidRow(plotlyOutput(outputId='enrPlot'))
             )
    )
}
