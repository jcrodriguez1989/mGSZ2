library('GO.db');
library('KEGG.db');
library('org.Hs.eg.db');

gsetsKeys <- function() {
    # add at the end!!
    gsetsKeys <- rbind(
        goKeys(),
        keggKeys()
    );
    return(gsetsKeys);
}

goKeys <- function() {
    gsetsKeys <- matrix(ncol=2, byrow=!FALSE, c(
        'GO Biological Process', 'BP',
        'GO Molecular Function', 'MF',
        'GO Cellular Component', 'CC'
    ));
    return(gsetsKeys);
}

keggKeys <- function() {
    gsetsKeys <- matrix(ncol=2, byrow=!FALSE, c(
        'KEGG', 'KEGG'
    ));
    return(gsetsKeys);
}

gsetsLib <- function(dbase=NULL) {
    aux <- gsetsKeys();
    res <- aux[,2]; names(res) <- aux[,1];
    if (!is.null(dbase)) {
        if (dbase == 'GO') {
            res <- res[1:3];
        } else if (dbase == 'KEGG') {
            res <- res[[4]];
        }
    }
    return(res);
}

loadGsets <- function(selGsets, ownGsets=list()) {
    goGsets <- loadGoGsets(selGsets); length(goGsets)
    keggGsets <- loadKeggGsets(selGsets); length(keggGsets)

    res <- c(goGsets, keggGsets); length(res)

    return(res);
}

loadGoGsets <- function(selGsets) {
    if (length(selGsets) == 0) return(list())

    go <- org.Hs.egGO2ALLEGS;
    goIds <- mappedkeys(go);

    ontologies <- Ontology(goIds);
    resGsets <- go[!is.na(ontologies) & ontologies %in% selGsets];
    resGsets <- as.list(resGsets);
    resGsets <- lapply(resGsets, unique);
    return(resGsets);
}

loadKeggGsets <- function(selGsets, specie='hsa') {
    if (!'KEGG' %in% selGsets) return(list())

    kegg <- KEGGPATHID2EXTID;
    keggIds <- mappedkeys(kegg);

    actKegg <- kegg[grep(specie, keggIds)];
    actKegg <- as.list(actKegg);
    actKegg <- lapply(actKegg, unique);

    return(actKegg);
}
