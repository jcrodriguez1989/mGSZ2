exprMatrixServer <- function(input) {
    infile <- input$inputFile;
    if (is.null(infile)) {
        # User has not uploaded a file yet
        return(NULL);
    }
    mtrx <- read.csv(infile$datapath);
    rownames(mtrx) <- mtrx[,1];
    mtrx <- mtrx[,-1];
    return(mtrx);
}

parseConds <- function(input, session) {
    condstext <- input$condsInput;
    conds <- trimws(strsplit(condstext, ',')[[1]]);
    if (length(conds) == 0 || (length(conds) == 1 && conds == '') || length(unique(conds)) < 2)
        conds <- NULL;

    contrasts <- list();
    if (!is.null(conds))
        contrasts <- apply(combn(unique(conds), 2),2,function(x) paste0(x, collapse=' -vs- '));
    updateSelectInput(session, inputId='selectedCtrst', choices=contrasts);

    return(conds);
}
