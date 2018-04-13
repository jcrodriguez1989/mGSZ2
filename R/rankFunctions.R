#' @importFrom limma treat contrasts.fit lmFit
#' @importFrom stats model.matrix
mGszEbayes <- function(exprData, classes) {
    design <- model.matrix(~classes);
    fit <- lmFit(exprData, design);
    # the contrast is the classes difference (0 for intercept, 1 for diff)
    fit <- treat(contrasts.fit(fit, c(0,1)));

    res <- fit$t;
    return(res);
}
