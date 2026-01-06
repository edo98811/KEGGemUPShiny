# Installing the package

You can install the KEGGemUPShiny package from GitHub using the remotes package. 

```R
install.packages("remotes")
remotes::install_github("yourusername/KEGGemUPShiny")
``` 

# Running the Application

To run the application, you will need to have a dde experiment object loaded in the R environement. You can then call the main function like this:

```
KEGGemUPShiny(dde_object)
```

This will launch the Shiny application where you can interactively explore KEGG pathway enrichment results. This package provides a Shiny application for visualizing KEGG pathway enrichment results from differential gene expression experiments. It is based on the KEGGemUP package and allows users to interactively explore and visualize the results of differential expression analyses on KEGG pathways.