library('shiny');
library('shinyjs');
source('InputData/uiInputData.R');
source('Results/uiResults.R');

shinyUI(fluidPage(
    useShinyjs(),
    tabsetPanel(id='maintab', type='pills',
                inputDataUi(),
                resultsUi()
    )
))
