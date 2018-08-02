runMGSZ2 <- function(exprMatrix, classes, input, session) {
    selGsets <- input$gsetsSelection;

    progress <- Progress$new(session);
    progress$set(message='Loading gene sets.', detail='');
    gsets <- loadGsets(selGsets); # load selected gene sets

    actCtrst <- input$selectedCtrst;
    actCtrst <- strsplit(actCtrst, ' -vs- ')[[1]];

    exprMatrix <- exprMatrix[, classes %in% actCtrst];
    classes <- classes[classes %in% actCtrst];

    progress$set(message='Starting mGSZ2 analysis.', detail='');

    # todo: add rankFn
    rankFn <- mGSZ2:::mGszEbayes;
    if (input$rnaSeqData)
        rankFn <- 'RNA';

    mGszRes <- mGSZ2(exprMatrix, gsets, classes, rankFn);
    # genesRank <- rankFn(exprMatrix, classes);
    progress$close();
    updateTabsetPanel(session, inputId='maintab', selected='mGsz2Restab');

    return(mGszRes);
}
