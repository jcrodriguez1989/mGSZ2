library('DT');
library('plotly');
source('utils/geneSetsUtils.R');

inputDataUi <- function() {
    tabPanel(title='Input data', value='inputDatatab',
             br(),
             fluidRow(column(2, actionButton(inputId='runMgszbtn', label='RUN mGSZ2!!'))),
             br(),
             wellPanel(fluidRow(
                 column(4, numericInput('nPerm', label='Number of permutations', value=200, min=1, step=1)),
                 column(6, selectInput(inputId='gsetsSelection', label='Selected gene sets', multiple=!FALSE,
                                        choices=gsetsLib(), selected=gsetsLib())
                            # todo: upload gset file
                 )
             )),
             wellPanel(
                 fluidRow(
                     column(7, fileInput(inputId='inputFile', label='Load expression matrix',
                                         accept=c('.xls', '.xlsx', '.csv', '.tsv'))),
                     column(2,
                            checkboxInput(inputId='rnaSeqData', label='RNA-Seq data', value=FALSE),
                            checkboxInput(inputId='printboxplot', label='Show boxplot', value=!FALSE)),
                     column(3, disabled(actionButton(inputId='symb2entrezbtn', label='Symbol2Entrez')))
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
