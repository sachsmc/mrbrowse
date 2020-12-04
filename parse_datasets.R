## new datasets are added to "txtfiles", then run this before building the package.

library(data.table)

datapath <- "txtfiles"
datanames <- list.files(datapath, pattern = "*.txt")

qtlnames <- list.files(file.path(datapath, "pQTLs"))
ctnames <- list.files(file.path(datapath, "CVDI"))

dflist <- lapply(qtlnames, function(j) {

  dtmp <- read.table(file.path(datapath, "pQTLs", j), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  colnames(dtmp) <- tolower(colnames(dtmp))
  if("p" %in% colnames(dtmp) & !"pval" %in% colnames(dtmp)) {
    dtmp[["pval"]] <- dtmp[["p"]]
    dtmp[["p"]] <- NULL
  }
  dtmp <- subset(dtmp, exposure != "" & outcome != "")
  dtmp

})

qtldf <- rbindlist(dflist, use.names = TRUE, fill = TRUE)
save(qtldf, file = "qtldf.rda")


dflist2 <- lapply(ctnames, function(j) {

  dtmp <- read.table(file.path(datapath, "CVDI", j), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  colnames(dtmp) <- tolower(colnames(dtmp))
  if("p" %in% colnames(dtmp) & !"pval" %in% colnames(dtmp)) {
    dtmp[["pval"]] <- dtmp[["p"]]
    dtmp[["p"]] <- NULL
  }
  dtmp <- subset(dtmp, exposure != "" & outcome != "")
  dtmp

})

ctdf <- rbindlist(dflist2, use.names = TRUE, fill = TRUE)
save(ctdf, file = "ctdf.rda")


qtlivres <- data.table(read.table(file.path(datapath, "Table_CVDI_pQTLs_R_shiny_app.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE))
save(qtlivres, file = "qtlivres.rda")


metadata <- data.table(read.table(file.path(datapath, "table_disease_outcomes_reference.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE))
save(metadata, file = "metadata.rda")


