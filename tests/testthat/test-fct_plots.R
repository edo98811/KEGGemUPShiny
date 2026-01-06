library(SummarizedExperiment)

# Create a simple fake SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(counts = matrix(rnbinom(100, mu = 5, size = 2), nrow = 10, ncol = 10)),
  rowData = DataFrame(
    SYMBOL = c("BRCA1", "TP53", "EGFR", "MYC", "PTEN", "KRAS", "AKT1", "PIK3CA", "BRAF", "ALK"),
    ENTREZID = c("672", "7157", "1956", "4609", "5728", "3845", "207", "5290", "673", "238")
  )
)

# Set Ensembl IDs as rownames
rownames(se) <- c(
  "ENSG00000012048", "ENSG00000141510", "ENSG00000005471", "ENSG00000136997",
  "ENSG00000171862", "ENSG00000133703", "ENSG00000142208", "ENSG00000121879",
  "ENSG00000157764", "ENSG00000139083"
)