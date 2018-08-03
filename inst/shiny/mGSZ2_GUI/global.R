# libraries
library(shiny)
library(jsonlite)

# create help data frame
steps <- data.frame(
    step=1:7,
    element=paste0('#helpStep', 1:7),
    intro=c(
        'Upload an expression matrix file. It can be a file with extension xls, xlsx, csv, or tsv. The first column will be used as gene identifier, the first row as subject identifier. If it has multiple samples for the same gene, then they will be averaged.',
        'Select whether the data comes from RNA-Seq technology or microarrays.',
        'Select whether the gene IDs of the expression matrix are in Symbol or Entrez. This will only be used if analyzing the GO or KEGG gene sets automatically loaded by this tool.',
        "Indicate to which category/group each subject in the expression matrix belongs to, separated by ','. There must be a number of categories equal to the number of columns in the expression matrix. There must be at least two different categories. Then click on the 'Load' button. And select the contrast you want to test.",
        'Upload as many gene set files (.gmt) as you want.',
        'Select at least one group of gene sets to be evaluated.',
        "And start the mGSZ2 analysis by clicking on the 'RUN mGSZ2!!' button. Be patient, this process can take up to three hours depending on the size of the expression matrix and the number of gene sets to be evaluated."
    ),
    position='auto'
);
