#' @importFrom limma contrasts.fit lmFit treat
#' @importFrom stats model.matrix
maRank <- function(exprData, classes) {
    # should be the same, -1meanA + 1meanC
    # design <- model.matrix(~classes+0); contrast <- c(-1, 1);
    design <- model.matrix(~classes);
    contrast <- c(0, 1);
    fit <- lmFit(exprData, design);
    # the contrast is the classes difference (0 for intercept, 1 for diff)
    fit <- treat(contrasts.fit(fit, contrast));

    res <- fit$t;
    stopifnot(ncol(res) == 1);
    return(res);
}

## voom + limma
## counts a numeric matrix containing raw counts. Counts must be non-negative
## and NAs are not permitted
#' @importFrom limma contrasts.fit lmFit treat voom
#' @importFrom stats model.matrix
rnaRank <- function(exprData, classes) {
    #     should be the same, -1meanA + 1meanC
    #     design <- model.matrix(~classes+0); contrast <- c(-1, 1);
    ##     it is doing -1cond[[1]] + cond[[2]]
    design <- model.matrix(~classes);
    contrast <- c(0, 1);
    newExpr <- voom(exprMatrix, design);
    fit <- lmFit(newExpr, design);
    # the contrast is the classes difference (0 for intercept, 1 for diff)
    fit <- treat(contrasts.fit(fit, contrast));

    res <- fit$t;
    stopifnot(ncol(res) == 1);
    return(res);
}
