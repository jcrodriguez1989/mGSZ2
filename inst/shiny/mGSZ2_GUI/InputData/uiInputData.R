library('DT');
library('plotly');
source('utils/geneSetsUtils.R');

inputDataUi <- function() {
    tabPanel(title='Input data', value='inputDatatab',
             br(),
             fluidRow(column(2, actionButton(inputId='runMgszbtn', label='RUN mGSZ2!!'))),
             br(),
             wellPanel(fluidRow(
                 column(2, numericInput('nPerm', label='Number of permutations', value=200, min=1, step=1)),
                 column(5, selectInput(inputId='gsetsSelection', label='Selected gene sets', multiple=!FALSE,
                                        choices=gsetsLib(), selected=gsetsLib())
                 ),
                 column(5, fileInput(inputId='gmtFile', label='Load gene set file (.gmt)',
                                     accept=c('.gmt'))
                 )
             )),
             wellPanel(
                 fluidRow(
                     column(7, fileInput(inputId='inputFile', label='Load expression matrix',
                                         accept=c('application/vnd.ms-excel', # xls
                                                  'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet', # xlsx
                                                  'text/csv', # csx
                                                  'text/tab-separated-values') # tsv
                                         )),
                     column(2,
                            checkboxInput(inputId='rnaSeqData', label='RNA-Seq data', value=FALSE),
                            checkboxInput(inputId='printboxplot', label='Show boxplot', value=!FALSE)),
                     column(3, checkboxInput(inputId='geneSymbol', label='Using gene symbol (if not, Entrez)'))
                 ),
                 fluidRow(
                     column(9, textInput(inputId='condsInput', label='Conditions', width='100%',
                                         placeholder='Each column\'s condition. Comma separated.')),
                     column(1, br(), actionButton(inputId='loadCondsbtn', label='Load')),
                     column(2, selectInput(inputId='selectedCtrst', label='Contrast', choices=list()))
                 )
             ),
             wellPanel(fluidRow(
                 column(10, DTOutput(outputId='inputMatrixDtable')),
                 # conditionalPanel('input.printboxplot', column(2, plotOutput(outputId='normboxplot')))
                 column(2, plotlyOutput(outputId='normboxplot'))
             ))
    )
}
