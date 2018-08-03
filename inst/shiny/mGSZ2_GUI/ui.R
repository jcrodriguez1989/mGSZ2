library('shiny');
library('shinyjs');
source('InputData/uiInputData.R');
source('Results/uiResults.R');

shinyUI(fluidPage(

    # Include IntroJS styling
    includeCSS("introjs.min.css"),

    # Include IntroJS library
    includeScript("intro.min.js"),

    # Include JavaScript code to make shiny communicate with introJS
    includeScript("app.js"),

    useShinyjs(),
    tabsetPanel(id='maintab', type='pills',
                inputDataUi(),
                resultsUi()
    )
))
