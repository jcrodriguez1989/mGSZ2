library('mGSZ');
library('mGSZ2');

context('mGSZ2 equals mGSZ');

test_that('mGSZ2 results equals mGSZ', {
    data('dummyData');

    set.seed(8818);
    mGSZres <- mGSZ(dummyData$x, dummyData$y, dummyData$l);
    mGSZres <- mGSZres$mGSZ;
    mGSZres <- mGSZres[, c("gene.sets", "gene.set.scores", "pvalue")];

    set.seed(8818);
    mGSZ2res <- mGSZ2(dummyData$x, dummyData$y, dummyData$l);
    mGSZ2res <- mGSZ2res[, c("gene.sets", "mGszScore", "pvalue")];

    expect(all.equal(mGSZres$gene.sets, mGSZ2res$gene.sets));
    expect(all.equal(mGSZres$gene.set.scores, abs(mGSZ2res$mGszScore)));
    expect(all.equal(mGSZres$pvalue, mGSZ2res$pvalue));
})
