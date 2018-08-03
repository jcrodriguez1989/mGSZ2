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

loadGsets <- function(selGsets, ownGsets, inSymbol) {
    goGsets <- loadGoGsets(selGsets); length(goGsets)
    keggGsets <- loadKeggGsets(selGsets); length(keggGsets)
    ownGSs <- loadOwnGsets(selGsets, ownGsets); length(ownGSs)

    if (inSymbol) {
        # Entrez to Symbol in loaded gene sets.
        toSymb <- as.list(org.Hs.egSYMBOL);

        aux <- attr(goGsets, 'GeneSetNames');
        goGsets <- lapply(goGsets, function(actGset) unique(unlist(toSymb[actGset])));
        attr(goGsets, 'GeneSetNames') <- aux;

        aux <- attr(keggGsets, 'GeneSetNames')
        keggGsets <- lapply(keggGsets, function(actGset) unique(unlist(toSymb[actGset])));
        attr(keggGsets, 'GeneSetNames') <- aux;
    }

    res <- c(goGsets, keggGsets, ownGSs); length(res)
    res <- res[unique(names(res))];
    attr(res, 'GeneSetNames') <-
        c(attr(goGsets, 'GeneSetNames'),
          attr(keggGsets, 'GeneSetNames'),
          attr(ownGSs, 'GeneSetNames'));

    return(res);
}

loadGoGsets <- function(selGsets) {
    if (length(selGsets) == 0) return(list())

    go <- org.Hs.egGO2ALLEGS;
    goIds <- mappedkeys(go);

    ontologies <- Ontology(goIds);
    resGsets <- go[!is.na(ontologies) & ontologies %in% selGsets];
    resGsets <- as.list(resGsets);
    resGsets <- lapply(resGsets, unique); # remove duplicate genes
    goNames <- Term(names(resGsets));
    attr(resGsets, 'GeneSetNames') <- goNames;
    return(resGsets);
}

loadKeggGsets <- function(selGsets, specie='hsa') {
    if (!'KEGG' %in% selGsets) return(list())

    kegg <- KEGGPATHID2EXTID;
    keggIds <- mappedkeys(kegg);

    actKegg <- kegg[grep(specie, keggIds)];
    actKegg <- as.list(actKegg);
    actKegg <- lapply(actKegg, unique); # remove duplicate genes
    keggNames <- unlist(as.list(KEGGPATHID2NAME)[substring(names(actKegg), 4)]);
    attr(actKegg, 'GeneSetNames') <- keggNames;
    return(actKegg);
}

loadOwnGsets <- function(selGsets, ownGsets) {
    thisGsets <- intersect(selGsets, names(ownGsets));
    actGsets <- ownGsets[thisGsets];
    actGsets <- unname(actGsets);

    resGsets <- do.call(c, actGsets);

    attr(resGsets, 'GeneSetNames') <- do.call(c, lapply(actGsets, function(x) attr(x, 'GeneSetNames')));
    return(resGsets);
}

##Simple function to read in a .gmt file and return a list of pathways
read_gmt <- function(filename){
    if (!grepl("\\.gmt$",filename)[1]) {
        print("Pathway information must be a .gmt file");
        return(list());
    }
    geneSetDB <- readLines(filename);                               ##read in the gmt file as a vector of lines
    geneSetDB <- strsplit(geneSetDB,"\t");                          ##convert from vector of strings to a list
    names(geneSetDB) <- sapply(geneSetDB,"[",1);                    ##move the names column as the names of the list
    descr <- sapply(geneSetDB,"[",2);                               ##get gene sets description
    names(descr) <- names(geneSetDB);
    geneSetDB <- lapply(geneSetDB, "[",-1:-2);                      ##remove name and description columns
    geneSetDB <- lapply(geneSetDB, function(x){x[which(x!="")]});   ##remove empty strings
    descr <- descr[unique(names(geneSetDB))];                       ##remove duplicate gene sets
    geneSetDB <- geneSetDB[unique(names(geneSetDB))];               ##remove duplicate gene sets
    geneSetDB <- lapply(geneSetDB, unique);                         ##remove duplicate genes

    attr(geneSetDB, 'GeneSetNames') <- descr;
    return(geneSetDB)
}
