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
#' @param x Gene expression data matrix (rows as genes and columns as samples)
#' @param y Gene set data (dataframe/table/matrix/list)
#' @param l Vector of response values (example:1,2)
#' @param rankFn One of 'MA', 'RNA' if data comes from microarrays-compatible
#' or RNA-Seq technologies respectively. Or a function with inputs
#' (x, l), and returns a vector with the same length as the nrow(x).
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
#' @include auxFuns.R
#'

# rankFn='MA'; min.sz=5; pv=0; w1=0.2; w2=0.5; vc=10; p=200
mGSZ2 <- function(x, y, l, rankFn='MA', min.sz=5, pv=0, w1=0.2, w2=0.5, vc=10,
                 p=200) {
    if (length(unique(l)) != 2)
        stop('l must have exactly two different categories.');
    if (length(l) != ncol(x))
        stop('l must have the same length as columns of x.');

    # filtering inputs, as required
    filteredInputs <- filterInputs(x, y, min.sz);
    exprData <- filteredInputs$exprData;
    gSets <- filteredInputs$gSets;
    rm(filteredInputs);

    # get gene rankings (also for permuted data)
    nPerm <- p;
    preVar <- pv;
    wgt2 <- w2;
    varConstant <- vc;

    rankings <- getRankings(exprData, l, nPerm, rankFn);
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

            return(list(actES=actES, rawES=rawES, impGenes=impGenes));
        });
        names(actESs) <- names(actGsets);
        return(actESs);
    }))
    rm(sVMC_dec); rm(sVMC_inc); rm(zVars_dec); rm(zVars_inc);

    realESs <- do.call(c, lapply(enrichScores, function(x) x$actES));
    rawESs <- do.call(c, lapply(enrichScores, function(x) x$rawES));
    impGenes <- do.call(c, lapply(enrichScores, function(x)
                            paste(x$impGenes, collapse=', ')));

    permESs <- do.call(cbind, bplapply(seq_len(ncol(rankings)-1)+1,
        function(i) {
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

    return(result);
}