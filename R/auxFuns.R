# exprData=expr.data; gSets=gene.sets; minGsetSize=min.cl.sz
filterInputs <- function(exprData, gSets, minGsetSize) {
    # keep only genes in exprData and gene sets
    exprDataGenes <- rownames(exprData);
    geneSetsGenes <- do.call(c, gSets);
    commonGenes <- intersect(exprDataGenes, geneSetsGenes);

    exprData <- exprData[commonGenes,];
    gSets <- lapply(gSets, intersect, commonGenes);

    # remove gene sets with less than minSz genes (detected in the experiment)
    gSets <- gSets[lapply(gSets, length) > minGsetSize];
    return(list(exprData=exprData, gSets=gSets));
}

#' @importFrom BiocParallel bpparam bplapply
#' @include rankFunctions.R
# classes=l
getRankings <- function(exprData, classes, nPerm, rankFn, rankInParallel) {
    if (is.character(rankFn) && rankFn == 'MA') {
        rankFn <- mGszEbayes;
    } else if (is.character(rankFn) && rankFn == 'RNA') {
        stop('Not developed yet');
        # rankFn <-
    } else if (!is.function(rankFn)) {
        stop('rankFn must be one of "MA", "RNA", or an R function');
    } ## todo: add the possibility of having a matrix with the perms

    # generate permutations
    perms <- replicate(nPerm, sample(1:ncol(exprData), replace=FALSE));
    permsConds <- do.call(cbind, lapply(seq_len(ncol(perms)), function(i) {
        as.character(classes[perms[,i]]);
    }))
    perms <- t(perms[,!duplicated(t(permsConds))]);
    rm(permsConds);

    nPerm <- nrow(perms);
    print(paste("Number of unique permutations:", nPerm));

    rankings <- matrix(0, nrow=nrow(exprData), ncol=nPerm+1);
    # gene scores for real data
    rankings[,1] <- rankFn(exprData, classes);

    # gene scores for permuted data
    # I am passing MIGSA's not exported functions to bplapply to avoid
    # SnowParam environment errors
    if (rankInParallel) {
        whichLapply <- bplapply;
        print(paste("Getting ranking at cores:", bpparam()$workers));
    } else {
        whichLapply <- lapply;
    }
    
    rankings[,-1] <- do.call(cbind,
        whichLapply(seq_len(nPerm), function(i) {
            actRank <- rankFn(exprData, classes[perms[i,]]);
            return(actRank);
        }));

    ## todo: maybe the best alternative would be delete the gene from
    # everywhere
    rankings[is.na(rankings)] <- mean(rankings, na.rm=TRUE);
    # genes which are NA are given the mean value of all genes,
    # so they dont contribute to enrichment score.

    return(rankings);
}

countHygeVarMean <- function (M, K) {
    N <- matrix(rep(seq_len(M), length(K)), ncol=length(K));
    K <- matrix(rep(K, M), byrow=TRUE, ncol=length(K));
    mean <- N * K/M;
    var <- N * (K * (M - K)/M^2) * ((M - N)/(M - 1));
    return(list(mean=mean, var=var));
}

countProbSum <- function (M, hyge.mean, hyge.var) {
    tulos <- matrix(0, nrow=length(hyge.mean));
    N_tab <- matrix(c(2:M), byrow=FALSE);
    tulos[1] <- hyge.mean[1];
    tulos[2:M] <- N_tab/(N_tab-1) * hyge.mean[2:M] - (hyge.mean[2:M]^2 +
                  hyge.var[2:M])/(N_tab-1);
    return(tulos);
}

sumVarMeanCalc <- function(ranking, preVar, normFactors) {
    divider <- seq_along(ranking);
    mean_table <- cumsum(ranking)/divider;
    mean_table_sq <- mean_table^2;

    # not sure if it is +2*preVar or just +preVar , in mGSZ it does +preVar
    # and again +preVar
    var_table <- cumsum(ranking^2)/divider - (mean_table)^2 + 2*preVar;

    z_means <- do.call(cbind, lapply(seq_len(ncol(normFactors$normMean)),
        function(j) {
            mean_table * (2 * normFactors$normMean[,j] - divider);
    }));
    colnames(z_means) <- colnames(normFactors$normMean);

    z_vars <- do.call(cbind, lapply(seq_len(ncol(normFactors$normVar)),
        function(j) {
            4 * (var_table * normFactors$probSum[,j] + mean_table_sq *
                normFactors$normVar[,j]);
    }));
    colnames(z_vars) <- colnames(normFactors$normVar);

    return(list(Z_means=z_means, Z_vars=z_vars));
}

zVarCalc <- function(z_vars, wgt2, varConstant) {
    medianPart <- matrix(matrixStats::rowMedians(t(z_vars)) *
                     wgt2, nrow=nrow(z_vars), ncol=ncol(z_vars), byrow=!FALSE);
    z_var <- (z_vars + medianPart + varConstant)^0.5;
    return(z_var);
}

getEnrichScore <- function(ranking, actGset, zM_d, zM_i, zV_d, zV_i) {
    startVal <- 5; # hard coded
    nMG <- !names(ranking) %in% actGset;
    ranking[nMG] <- -ranking[nMG];

    es_dec <- cumsum(ranking) - zM_d;
    es_inc <- cumsum(rev(ranking)) - zM_i;

    es_dec[1:startVal] <- 0; es_inc[1:startVal] <- 0;

    es_dec <- es_dec/zV_d;
    es_inc <- es_inc/zV_i;

    return(c(es_dec, es_inc));
}

#' @importFrom compiler cmpfun
getEnrichScore_c <- cmpfun(getEnrichScore);

#' @importFrom matrixStats colSds colVars
#' @importFrom ismev gum.fit
getPvalues <- function(realESs, permESs) {
    # some normalization
    mean.prof <- colMeans(permESs);
    std.prof <- matrixStats::colSds(permESs);
    std.prof[std.prof == 0] <- 0.1;
    realESs <- (realESs - mean.prof)/std.prof;
    permESs <- do.call(cbind, lapply(seq_len(ncol(permESs)), function(k) {
        (permESs[, k] - mean.prof[k])/std.prof[k];
    }))
    rm(mean.prof); rm(std.prof);

    col.ind <- matrixStats::colVars(permESs) > 0;
    permESs <- permESs[, col.ind];
    realESs <- realESs[col.ind];
    ev.p.val.class <- rep(0, length(realESs))

    pvals <- do.call(c, lapply(seq_along(realESs), function(k) {
        ev.param.class <- ismev::gum.fit(as.vector(permESs[, k]),
                                         show=FALSE)$mle
        logEVcdf(realESs[[k]], ev.param.class);
    }))
    return(pvals);
}

logEVcdf <- function (x, fitParams) {
    mu <- fitParams[[1]];
    sigma <- fitParams[[2]];
    x <- (mu - x)/sigma;
    out <- rep(0, length(x));
    if (!(is.vector(x))) {
        out <- matrix(out, nrow(x), ncol(x));
    }
    po1 <- which(x < 5);
    out[po1] <- -log(1 - exp(-exp(x[po1])));
    x <- x[-po1];
    out[-po1] <- -x + exp(x)/2 - exp(2 * x)/24 + exp(4 * x)/2880;
    out <- out/log(10);
    out <- 10^(-out);
    return(out);
}

