#' @rdname mGSZ2
#' @export
#'
#' @title
#' mGSZ2: Improved Modified Gene Set Z-score method for gene set analysis.
#'
#' @description
#' \code{mGSZ2} Gene set analysis based on Gene Set Z scoring function and
#' asymptotic p-value
#'
#' @param x Gene expression data matrix (rows as genes and columns as samples).
#' Raw counts for RNA-Seq data.
#' @param y Gene set data (dataframe/table/matrix/list)
#' @param l Vector of response values (example:1,2)
#' @param rankFn One of 'MA', 'RNA' if data comes from microarrays-compatible
#' or RNA-Seq technologies respectively. Or a function with inputs (x, l), and
#' returns a vector with the same length as the nrow(x).  Or a matrix of
#' rankings genes as rows (must have genes as rownames), first column for real
#' ranking, rest of columns for permuted rankings, in this case mGSZ2 will not
#' use x and l inputs (can be NA).
#' @param min.sz Minimum size of gene sets (number of genes in a gene set) to
#' be included in the analysis
#' @param pv Estimate of the variance associated with each observation
#' @param w1 Weight 1, parameter used to calculate the prior variance
#' obtained with class size var.constant. This penalizes especially small
#' classes and small subsets. Default is 0.2. Values around 0.1 - 0.5
#' are expected to be reasonable
#' @param w2 Weight 2, parameter used to calculate the prior variance
#' obtained with the same class size as that of the analyzed class. This
#' penalizes small subsets from the gene list. Default is 0.5. Values around
#' 0.3 and 0.5 are expected to be reasonable
#' @param vc Size of the reference class used with wgt1. Default is 10
#' @param p Number of permutations for p-value calculation
#' @param rankInParallel If FALSE, permutation gene rankings will be calculated
#' sequentially. Useful if rankFn is provided and is already parallelized.
#'
#' @details A function for Gene set analysis based on Gene Set Z-scoring
#' function and asymptotic p- value. It differs from GSZ (Toronen et al 2009)
#' in that it implements asymptotic p-values instead of empirical p-values.
#' Asymptotic p-values are based on fitting suitable distribution model to
#' the permutation data. Unlike empirical p-values, the resolution of
#' asymptotic p-values are independent of the number of permutations and
#' hence requires consideralbly fewer permutations. In addition to GSZ, this
#' function allows the users to carry out analysis with seven other scoring
#' functions and compare the results.
#'
#' @return Dataframe with gene sets (in decreasing order based on the
#' significance) reported by mGSZ method, scores, p-values, and list of genes
#' that contributed to the enrichment.
#'
#' @examples
#' \dontrun{
#' data(dummyData);
#' mGSZres <- mGSZ2(dummyData$x, dummyData$y, dummyData$l);
#' }
#' @author
#' Juan Cruz Rodriguez \email{jcrodriguez@@bdmg.com.ar}
#'
#' @importFrom BiocParallel bplapply
#' @importFrom futile.logger flog.info
#' @include auxFuns.R
#'

# rankFn='MA'; min.sz=5; pv=0; w1=0.2; w2=0.5; vc=10; p=200; rankInParallel=F;
mGSZ2 <- function(x, y, l, rankFn='MA', min.sz=5, pv=0, w1=0.2, w2=0.5, vc=10,
                 p=200, rankInParallel=!F) {
    if (length(unique(l)) != 2 && !is.matrix(rankFn))
        stop('l must have exactly two different categories.');
    if (!is.matrix(rankFn) && length(l) != ncol(x))
        stop('l must have the same length as columns of x.');
    if (is.matrix(rankFn)) {
        x <- rankFn;
        p <- ncol(rankFn)-1;
    }

    # filtering inputs, as required
    filteredInputs <- filterInputs(x, y, min.sz);
    exprData <- filteredInputs$exprData;
    gSets <- filteredInputs$gSets;
    attr(gSets, 'GeneSetNames') <- attr(y, 'GeneSetNames'); # if y has this attribute then keep it
    rm(filteredInputs); rm(x); rm(y);

    # get gene rankings (also for permuted data)
    nPerm <- p;
    preVar <- pv;
    wgt2 <- w2;
    varConstant <- vc;

    if (is.matrix(rankFn)) {
        flog.info("Rankings provided by user.");
        rankings <- exprData;
    } else {
        flog.info("Getting rankings.");
        rankings <- getRankings(exprData, l, nPerm, rankFn, rankInParallel);
    }

    if (any(is.na(rankings) | is.infinite(rankings))) {
        rankings[is.infinite(rankings)] <- NA;
        warning(paste('Ranking function returned', sum(is.na(rankings)),
            'non-finite values'));

        rankings <- t(apply(rankings, 1, function(x) {
            x[is.na(x)] <- mean(x, na.rm=!F);
            return(x);
        }));
    }
    rownames(rankings) <- rownames(exprData);
    rm(nPerm);

    # calculate normalizing values
    setSizes <- unlist(lapply(gSets, length));
    uniqSizes <- unique(setSizes);
    normValues <- countHygeVarMean(nrow(rankings), uniqSizes);
    normMean <- normValues$mean; colnames(normMean) <- uniqSizes;
    normVar <- normValues$var; colnames(normVar) <- uniqSizes;
    rm(normValues);

    probSum <- do.call(cbind, lapply(seq_along(uniqSizes), function(j) {
        countProbSum(nrow(rankings), normMean[, j], normVar[, j]);
    }));
    colnames(probSum) <- uniqSizes;
    normFactors <- list(normMean=normMean, normVar=normVar,
                        probSum=probSum);

    # get extra info for real data
    realRank <- rankings[,1];
    realRank <- sort(realRank, decreasing=!FALSE);
    sVMC_dec <- sumVarMeanCalc(realRank, preVar, normFactors);
    sVMC_inc <- sumVarMeanCalc(rev(realRank), preVar, normFactors);
    zVars_dec <- zVarCalc(sVMC_dec$Z_vars, wgt2, varConstant);
    zVars_inc <- zVarCalc(sVMC_inc$Z_vars, wgt2, varConstant);

    flog.info("Getting enrichment scores.");
    enrichScores <- do.call(c, bplapply(uniqSizes, function(actSize) {
        actGsets <- gSets[setSizes == actSize];
        actSize_c <- as.character(actSize);

        zM_d <- sVMC_dec$Z_means[, actSize_c];
        zM_i <- sVMC_inc$Z_means[, actSize_c];
        zV_d <- zVars_dec[, actSize_c];
        zV_i <- zVars_inc[, actSize_c];

        actESs <- lapply(actGsets, function(actGset) {
            diffScores <- getEnrichScore_c(realRank, actGset,
                                           zM_d, zM_i, zV_d, zV_i);
            # real data
            maxI <- which.max(abs(diffScores));
            rawES <- diffScores[[maxI]];
            nGenes <- length(realRank);

            if (maxI > nGenes) {
                impGenes <- names(diffScores)[(nGenes+1):maxI];
            } else {
                impGenes <- names(diffScores)[1:maxI];
            }
            impGenes <- intersect(impGenes, actGset);
            actES <- abs(rawES);

            return(list(actES=actES, rawES=rawES, impGenes=impGenes, ESs=diffScores));
        });
        names(actESs) <- names(actGsets);
        return(actESs);
    }))
    rm(sVMC_dec); rm(sVMC_inc); rm(zVars_dec); rm(zVars_inc);

    realESs <- do.call(c, lapply(enrichScores, function(x) x$actES));
    rawESs <- do.call(c, lapply(enrichScores, function(x) x$rawES));
    impGenes <- do.call(c, lapply(enrichScores, function(x)
                            paste(x$impGenes, collapse=', ')));
    ESmatrix <- do.call(rbind, lapply(enrichScores, function(x) x$ESs));
    stopifnot(all(c(names(realRank), rev(names(realRank))) == colnames(ESmatrix)));

    permESs <- do.call(cbind, bplapply(seq_len(ncol(rankings)-1)+1,
        function(i) {
            flog.info(paste0('Getting ESs for perm: ', i-1));
            ranking <- rankings[,i];
            ranking <- sort(ranking, decreasing=!FALSE);
            sVMC_dec <- sumVarMeanCalc(ranking, preVar, normFactors);
            sVMC_inc <- sumVarMeanCalc(rev(ranking), preVar, normFactors);

            zVars_dec <- zVarCalc(sVMC_dec$Z_vars, wgt2, varConstant);
            zVars_inc <- zVarCalc(sVMC_inc$Z_vars, wgt2, varConstant);

            # as norm factors are separated by gene set sizes, then in order to
            # consume less ram, lets separate by sizes.
            # In some cases it will not use the best of cores, but its for ram!
            enrichScores <- do.call(c, lapply(uniqSizes,
                function(actSize, sVMC_dec, sVMC_inc, zVars_dec, zVars_inc) {
                    actGsets <- gSets[setSizes == actSize];
                    actSize_c <- as.character(actSize);

                    act_z_means_dec <- sVMC_dec$Z_means[, actSize_c];
                    act_z_means_inc <- sVMC_inc$Z_means[, actSize_c];
                    act_zVars_dec <- zVars_dec[, actSize_c];
                    act_zVars_inc <- zVars_inc[, actSize_c];

                    actESs <- do.call(c, lapply(actGsets,
                         function(actGset, zM_d, zM_i, zV_d, zV_i) {
                             actES <- max(abs(getEnrichScore_c(ranking, actGset,
                                                   zM_d, zM_i, zV_d, zV_i)));
                             return(actES);
                         },
                         zM_d=act_z_means_dec, zM_i=act_z_means_inc,
                         zV_d=act_zVars_dec, zV_i=act_zVars_inc
                    ));
                    names(actESs) <- names(actGsets);
                    return(actESs);
            }, sVMC_dec=sVMC_dec, sVMC_inc=sVMC_inc,
            zVars_dec=zVars_dec, zVars_inc=zVars_inc))
            return(enrichScores);
        }
    ));

    stopifnot(all(names(realESs) == rownames(permESs)));

    permESs <- t(permESs);
    pvals <- getPvalues(realESs, permESs);
    names(pvals) <- names(realESs);

    result <- data.frame(gene.sets=names(pvals), pvalue=pvals,
                         mGszScore=rawESs, impGenes=impGenes);
    result <- result[order(pvals),];
    attr(result, 'GenesRanking') <- realRank;
    attr(result, 'ESs') <- ESmatrix;
    attr(result, 'GSs') <- gSets;

    return(result);
}
