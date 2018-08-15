library('mGSZ2');

runMGSZ2 <- function(exprMatrix, classes, ownGsets, input, session) {
    selGsets <- input$gsetsSelection;

    progress <- Progress$new(session);
    progress$set(message='Loading gene sets.', detail='');

    inSymbol <- input$geneSymbol;
    gsets <- loadGsets(selGsets, ownGsets, inSymbol); # load selected gene sets

    actCtrst <- input$selectedCtrst;
    actCtrst <- strsplit(actCtrst, ' -vs- ')[[1]];

    exprMatrix <- exprMatrix[, classes %in% actCtrst];
    classes <- classes[classes %in% actCtrst];

    progress$set(message='Starting mGSZ2 analysis.', detail='');

    rankFn <- 'MA';
    if (input$rnaSeqData)
        rankFn <- 'RNA';

    mGszRes <- mGSZ2(exprMatrix, gsets, classes, rankFn);
    # genesRank <- rankFn(exprMatrix, classes);
    progress$close();
    updateTabsetPanel(session, inputId='maintab', selected='mGsz2Restab');

    return(mGszRes);
}
