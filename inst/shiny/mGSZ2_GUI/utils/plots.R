library('ggplot2');
library('plotly');
library('reshape2');

ggboxplot <- function(matrix, conditions=NULL, colNames=c('Gene', 'Sample', 'Value')) {
    if (is.null(matrix))
        return(plotly_empty());
    meltedMatrix <- reshape2::melt(as.matrix(matrix));
    colnames(meltedMatrix) <- colNames;

    meltedMatrix[,1] <- as.character(meltedMatrix[,1]);
    meltedMatrix[,2] <- factor(as.character(meltedMatrix[,2]), levels=colnames(matrix)); # to keep samples order
    meltedMatrix[,3] <- as.numeric(as.character(meltedMatrix[,3]));

    if (!is.null(conditions) && length(conditions) == ncol(matrix))
        meltedMatrix$Conds <- as.character(rep(conditions, each=nrow(matrix)));

    p <- ggplot(data=meltedMatrix);

    if ('Conds' %in% colnames(meltedMatrix)) {
        p <- p + geom_boxplot(aes_string(x=colNames[[2]], y=colNames[[3]], color='Conds'));
    } else {
        p <- p + geom_boxplot(aes_string(x=colNames[[2]], y=colNames[[3]]));
    }

    p <- p + theme(axis.text.x=element_blank(), legend.position='top');
    p <- ggplotly(p);

    return(p);
}
