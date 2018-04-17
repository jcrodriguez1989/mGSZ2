#' @importFrom limma treat contrasts.fit lmFit
#' @importFrom stats model.matrix
mGszEbayes <- function(exprData, classes) {
#     should be the same, -1meanA + 1meanC
#     design <- model.matrix(~classes+0); contrast <- c(-1, 1); 
    design <- model.matrix(~classes);
    contrast <- c(0, 1);
    fit <- lmFit(exprData, design);
    # the contrast is the classes difference (0 for intercept, 1 for diff)
    fit <- treat(contrasts.fit(fit, contrast));

    res <- fit$t;
    return(res);
}
