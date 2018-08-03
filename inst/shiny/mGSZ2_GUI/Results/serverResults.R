renderMgszRes <- function(mGszRes) {
    if (is.null(mGszRes))
        return(mGszRes);

    gsets <- attr(mGszRes,'GSs');
    gsetNames <- attr(gsets, 'GeneSetNames');
    gsetSsz <- unlist(lapply(gsets, length));
    res <- mGszRes[, c('gene.sets', 'mGszScore', 'pvalue')];

    aux <- as.character(res$gene.sets);
    actNames <- gsetNames[aux];
    aux[!is.na(actNames)] <- actNames[!is.na(actNames)];
    res$gene.sets <- aux;

    res <- data.frame(res);
    colnames(res) <- c('Name', 'ES', 'P-value');
    rownames(res) <- mGszRes$gene.sets;
    res$FDR <- p.adjust(res[,'P-value'], method='fdr');
    res[,'#leading edge'] <- unlist(lapply(mGszRes$impGenes, function(x)
        length(strsplit(as.character(x), ', ')[[1]])));
    res[,'#genes'] <- gsetSsz[rownames(res)];

    return(res);
}

loadRes <- function(input) {
    infile <- input$loadResFile;
    if (is.null(infile)) {
        # User has not uploaded a file yet
        return(NULL);
    }
    mGSZres <- get(load(infile$datapath));

    if (is.null(mGSZres) || any(!c('GenesRanking', 'ESs', 'GSs') %in% names(attributes(mGSZres)))) {
        showNotification('Corrupted mGSZ results file.', type='error');
        mGSZres <- NULL;
    }
    return(mGSZres);
}

enrPlot <- function(mGSZres, input) {
    selectedRow <- isolate(input$mGsz2ResDtable_rows_selected);
    if (is.null(selectedRow))
        return(ggplot());
    res <- enrichmentPlot(mGSZres, rownames(mGSZres)[[selectedRow]]);
    res <- suppressWarnings(ggplotly(res));
    return(res);
}
