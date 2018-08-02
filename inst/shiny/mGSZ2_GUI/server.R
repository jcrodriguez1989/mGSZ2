library('shiny');
library('DT');
library('plotly');
source('InputData/serverInputData.R');
source('Results/serverResults.R');
source('utils/plots.R');
source('utils/runMGSZ2.R');

options(shiny.maxRequestSize=50*1024^2); # upload file top of 50MB

shinyServer(function(input, output, session) {
    ############## Input Data

    # Loads expression matrix
    exprMatrix <- reactive(exprMatrixServer(input));

    # Show expression matrix in data table
    output$inputMatrixDtable <- renderDT(exprMatrix(), selection='none');

    # plots subjects boxplot if desired
    output$normboxplot <- renderPlotly({
        if (input$printboxplot) {
            ggboxplot(exprMatrix(), classes());
        } else {
            plotly_empty();
        }
    });

    # parses the classes
    classes <- reactiveVal(NULL);
    observeEvent(input$loadCondsbtn, classes(parseConds(input, session)));

    # runs mGSZ if inputs are correct
    mGSZres <- reactiveVal(NULL);
    observeEvent(input$runMgszbtn, {
        if (is.na(input$nPerm) | input$nPerm < 1) {
            showNotification('Number of permutations should be > 0.', type='error');
            return();
        }
        if (is.null(input$gsetsSelection)) {
            showNotification('Select at least one gene set to test.', type='error');
            return();
        }
        mtrx <- exprMatrix();
        if (is.null(mtrx)) {
            showNotification('No expression matrix provided.', type='error');
            return();
        }
        conds <- classes();
        if (is.null(conds)) {
            showNotification('No conditions provided. One per column of the expression matrix, and at least two different.', type='error');
            return();
        }
        browser()
        showModal(modalDialog(title='Running mGSZ algorithm...', size='l', footer=NULL));
        mGSZres(runMGSZ2(mtrx, conds, input, session))
        removeModal();
    });

    ############## Results

    # Show mGSZ results in data table
    output$mGsz2ResDtable <- renderDT(renderMgszRes(mGSZres()), selection='single');

    # Do enrichment plot when button clicked
    observeEvent(input$enrScorePlotbtn, output$enrPlot <- renderPlotly(enrPlot(mGSZres(), input)));

    # Save mGSZ results
    output$saveResbtn <- downloadHandler(
        filename = function() paste0('mGSZ2res_', Sys.Date(), '.RData'),
        content = function(con) save(mGSZres(), file=con)
    )

    # Load mGSZ results
    observeEvent(input$loadResFile, mGSZres(loadRes(input)));

})
