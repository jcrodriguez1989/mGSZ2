library('limma');

exprMatrixServer <- function(input) {
    infile <- input$inputFile;
    if (is.null(infile))
        return(NULL);
    datapath <- infile$datapath;
    type <- infile$type;

    if (type == 'application/vnd.ms-excel') {
        mtrx <- read_excel(datapath);
        mtrx <- data.frame(mtrx);
    } else if (type == 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet') {
        mtrx <- read_excel(datapath);
        mtrx <- data.frame(mtrx);
    } else if (type == 'text/csv') {
        mtrx <- read.csv(datapath);
    } else if (type == 'text/tab-separated-values') {
        mtrx <- read.csv(datapath, sep='\t');
    } else {
        showNotification('Invalid file extension. Allowed: xls, xlsx, csv, tsv.', type='error');
        mtrx <- NULL;
    }

    if (!is.null(mtrx)) {
        showNotification('File loaded. First column as gene names, first row as subject names.', type='message');
        ids <- mtrx[,1]; mtrx <- mtrx[,-1];
        if (length(ids) != length(unique(ids))) {
            showNotification('Duplicated genes detected, replacing by mean value.', type='warning');
            mtrx <- as.data.frame(avereps(mtrx, ID=ids));
        } else {
            rownames(mtrx) <- mtrx[,1];
        }
    }

    if (ncol(mtrx) < 2) {
        showNotification('Expression matrix must have at least two subjects (columns).', type='error');
        mtrx <- NULL;
    }

    if (nrow(mtrx) < 2) {
        showNotification('Expression matrix must have at least two genes (rows).', type='error');
        mtrx <- NULL;
    }

    if (!all(unlist(lapply(seq_len(ncol(mtrx)), function(i) is.numeric(mtrx[,i]))))) {
        showNotification('Expression matrix values must be numbers.', type='error');
        mtrx <- NULL;
    }

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

readGmt <- function(input, actGsets) {
    infile <- input$gmtFile;
    if (is.null(infile))
        return(NULL);
    datapath <- infile$datapath;
    type <- sub('.*\\.', '', datapath);

    if (type != 'gmt') {
        showNotification('Invalid file extension. Allowed: gmt.', type='error');
        return(actGsets);
    }

    gsetName <- sub('\\.[^\\.]*$', '', infile$name);
    if (gsetName %in% names(actGsets)) {
        showNotification('File already loaded.', type='error');
        return(actGsets);
    }

    newGset <- read_gmt(datapath);
    if (length(newGset) == 0) {
        showNotification('Could not read gmt file.', type='error');
        return(actGsets);
    }

    actGsets[[gsetName]] <- newGset;
    return(actGsets);
}
