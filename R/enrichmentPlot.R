#' @rdname enrichmentPlot
#' @export
#'
#' @title
#' Enrichment plot.
#'
#' @description
#' \code{enrichmentPlot} Enrichment plot of a gene set resulting from mGSZ2.
#'
#' @param mGSZres Result object obtained after applying \code{mGSZ2} function.
#' @param geneSetId character with one ID present in mGSZres (a rowname of it).
#'
#' @details It plots the running enrichment score of a gene set.
#'
#' @return A ggplot2 object.
#'
#' @examples
#' \dontrun{
#' data(dummyData);
#' mGSZres <- mGSZ2(dummyData$x, dummyData$y, dummyData$l);
#' enrichmentPlot(mGSZres, rownames(mGSZres)[[1]]);
#' }
#' @author
#' Juan Cruz Rodriguez \email{jcrodriguez@@bdmg.com.ar}
#'
#' @importFrom ggplot2 ggplot geom_line ylab labs aes scale_color_continuous theme
#' @importFrom ggplot2 element_blank geom_segment geom_hline ggtitle geom_point
#'

enrichmentPlot <- function(mGSZres, geneSetId) {
    gSet <- attr(mGSZres, 'GSs')[[geneSetId]];

    if (is.null(gSet))
        stop('geneSetId not in mGSZres.');

    EScore <- attr(mGSZres, 'ESs')[geneSetId,, drop=F];
    if (which.max(abs(EScore)) > ncol(EScore)/2) {
        EScore <- EScore[,((ncol(EScore)/2)+1):ncol(EScore), drop=F];
    } else {
        EScore <- EScore[,1:(ncol(EScore)/2), drop=F];
    }
    genesRank <- attr(mGSZres, 'GenesRanking');
    genesRank <- genesRank[colnames(EScore)];

    plotData <- data.frame(
        Gene=factor(names(genesRank), levels=names(genesRank)),
        ES=EScore[1,],
        GeneRank=genesRank,
        # GeneSet=names(genesRank) %in% gSet);
        InGSet=names(genesRank) %in% gSet);

    ESpoint <- which.max(abs(EScore));
    ES <- EScore[which.max(abs(EScore))];

    ggp <- ggplot(plotData);
    ggp <- ggp + geom_line(aes(x=Gene, y=ES, group=1)); # running ES points line
    ggp <- ggp + suppressWarnings(geom_point(aes(x=Gene, y=ES, color=GeneRank, size=InGSet))) +
        scale_color_continuous(low='red', high='blue');

    # ggp <- ggp + suppressWarnings(geom_point(aes(x=Gene, y=ES, alpha=InGSet))) +
    #     scale_color_continuous(low='red', high='blue');
    # ggp <- ggp + lapply(seq_len(nrow(plotData)-1), function(i) {
    #     actRow <- plotData[i,];
    #     nextRow <- plotData[i+1,];
    #     geom_segment(aes(x=actRow$Gene, xend=nextRow$Gene,
    #                      y=actRow$ES, yend=nextRow$ES,
    #                      color=actRow$GeneRank));
    # })

    # this does not work in ggplotly (colored line), however it looks better
    # ggp <- ggp + geom_line(aes(x=Gene, y=ES, color=GeneRank, group=1)) +
    #     scale_color_continuous(low='red', high='blue'); # running ES points line
    # ggp <- ggp + geom_point(data=plotData[plotData$InGSet,], aes(x=Gene, y=ES, color=GeneRank));

    ggp <- ggp + theme(axis.text.x=element_blank()); # remove gene names from axis
    ggp <- ggp + geom_segment(x=ESpoint, xend=ESpoint,
                              y=min(plotData$ES), yend=max(plotData$ES),
                              linetype='dashed', size=0.1); # real ES
    ggp <- ggp + geom_segment(y=ES, yend=ES,
                              x=ESpoint, xend=0,
                              linetype='dashed', size=0.1);
    ggp <- ggp + geom_hline(yintercept=0, alpha=.25); # zero x line
    ggp <- ggp + ggtitle(paste('Enrichment plot:', geneSetId));
    ggp <- ggp + ylab('Enrichment score (ES)');
    ggp <- ggp + labs(color='Gene rank');
    ggp <- ggp + labs(alpha='In gene set');
    ggp <- ggp + theme(legend.position='bottom');

    return(ggp);
}
