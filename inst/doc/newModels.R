## ----setup, include=FALSE-----------------------------------------------------
# suppressMessages(library(digitalDLSorteR))
knitr::opts_chunk$set(
  echo = TRUE, dpi = 80, fig.width = 8, fig.height = 4.5, fig.align = "center",
  eval = digitalDLSorteR:::.checkPythonDependencies(alert = "none")
)

## ----workflowfull_newModels, fig.cap = "Workflow to build new context-specific deconvolution models", echo = FALSE----
knitr::include_graphics("workflow_building_models.png")

## -----------------------------------------------------------------------------
## loading the packages
suppressMessages(library(digitalDLSorteR))
suppressMessages(library(SingleCellExperiment))

## set seed for reproducibility
set.seed(123)
sce <- SingleCellExperiment(
  matrix(
    stats::rpois(50000, lambda = 5), nrow = 500, ncol = 100, 
    dimnames = list(paste0("Gene", seq(500)), paste0("RHC", seq(100)))
  ),
  colData = data.frame(
    Cell_ID = paste0("RHC", seq(100)),
    Cell_Type = sample(
      x = paste0("CellType", seq(5)), size = 100, replace = TRUE
    )
  ),
  rowData = data.frame(
    Gene_ID = paste0("Gene", seq(500))
  )
)

## ----loadData-----------------------------------------------------------------
DDLSToy <- loadSCProfiles(
  single.cell = sce, 
  cell.ID.column = "Cell_ID",
  gene.ID.column = "Gene_ID",
  min.cells = 0,
  min.counts = 0,
  project = "ToyExample"
)
DDLSToy

## ----loadFromFile, eval=FALSE-------------------------------------------------
#  ## this code will not be run
#  toyFiles <- c("countsMatrix.tsv.gz",
#                "cellsMetadata.tsv.gz",
#                "genesMetadata.tsv.gz")
#  
#  DDLSToy <- loadSCProfiles(
#    single.cell = toyFiles,
#    cell.ID.column = "Cell_ID",
#    gene.ID.column = "external_gene_name",
#    min.cells = 0, min.counts = 0,
#    project = "ToyExampleBreast"
#  )

## ---- eval = FALSE------------------------------------------------------------
#  DDLSToy <- loadSCProfiles(
#    single.cell = toyFiles, cell.ID.column = "Cell_ID",
#    gene.ID.column = "external_gene_name",
#    min.cells = 0, min.counts = 0,
#    file.backend = "singlecell_data.h5",
#    project = "ToyExampleBreast"
#  )

## -----------------------------------------------------------------------------
DDLSToy <- estimateZinbwaveParams(
  object = DDLSToy,
  cell.ID.column = "Cell_ID",
  gene.ID.column = "Gene_ID",
  cell.type.column = "Cell_Type",
  subset.cells = 40,
  threads = 1,
  verbose = TRUE
)

## -----------------------------------------------------------------------------
DDLSToy

## -----------------------------------------------------------------------------
DDLSToy <- simSCProfiles(
  object = DDLSToy,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  n.cells = 10,
  suffix.names = "_Simul",
  verbose = TRUE
)

## -----------------------------------------------------------------------------
DDLSToy

## ---- eval = FALSE------------------------------------------------------------
#  DDLSToy <- simSCProfiles(
#    object = DDLSToy,
#    cell.ID.column = "Cell_ID",
#    cell.type.column = "Cell_Type",
#    n.cells = 10,
#    suffix.names = "_Simul",
#    file.backend = "simulated_singlecell_data.h5",
#    block.processing = TRUE,
#    block.size = 20, # number of single-cell profiles simulated per batch
#    verbose = TRUE
#  )

## -----------------------------------------------------------------------------
## for reproducibility
set.seed(123)

## prior knowledge for prob.design argument
probMatrix <- data.frame(
  Cell_Type = paste0("CellType", seq(5)),
  from = c(rep(1, 2), 1, rep(30, 2)),
  to = c(rep(15, 2), 50, rep(70, 2))
)

DDLSToy <- generateBulkCellMatrix(
  object = DDLSToy,
  cell.ID.column = "Cell_ID",
  cell.type.column = "Cell_Type",
  prob.design = probMatrix,
  num.bulk.samples = 250,
  n.cells = 100,
  verbose = TRUE
)

## -----------------------------------------------------------------------------
DDLSToy

## -----------------------------------------------------------------------------
head(getProbMatrix(DDLSToy, type.data = "train"))
tail(getProbMatrix(DDLSToy, type.data = "train"))

## ----showProbPlot_newModels---------------------------------------------------
lapply(
  1:6, function(x) {
    showProbPlot(
      DDLSToy, type.data = "train", set = x, type.plot = "boxplot"
    )
  }
)

## ----simBulkProfiles----------------------------------------------------------
DDLSToy <- simBulkProfiles(
  object = DDLSToy, type.data = "both", pseudobulk.function = "MeanCPM"
)

## -----------------------------------------------------------------------------
DDLSToy

## ---- eval = FALSE------------------------------------------------------------
#  DDLSToy <- simBulkProfiles(
#    object = DDLSToy,
#    type.data = "both",
#    file.backend = "pseudobulk_samples.h5",
#    block.processing = TRUE,
#    block.size = 1000,
#    threads = 2
#  )

## ---- warning = FALSE---------------------------------------------------------
DDLSToy <- trainDigitalDLSorterModel(object = DDLSToy, scaling = "standarize")

## -----------------------------------------------------------------------------
DDLSToy

## ---- eval = FALSE------------------------------------------------------------
#  DDLSToy <- trainDigitalDLSorterModel(object = DDLSToy, on.the.fly = TRUE)

## -----------------------------------------------------------------------------
DDLSToy <- calculateEvalMetrics(object = DDLSToy)

## ----distErr1_newModels-------------------------------------------------------
distErrorPlot(
  DDLSToy,
  error = "AbsErr",
  x.by = "CellType",
  color.by = "CellType", 
  error.labels = FALSE, 
  type = "boxplot",
  size.point = 1
)

## ----distErr2_newModels, warning=FALSE, fig.height=4--------------------------
distErrorPlot(
  DDLSToy,
  error = "AbsErr",
  facet.by = "CellType",
  color.by = "nCellTypes", 
  type = "violinplot",
  size.point = 1
)

## ----distErr3_newModels, fig.height=4-----------------------------------------
distErrorPlot(
  DDLSToy,
  error = "AbsErr",
  color.by = "CellType", 
  facet.by = "nCellTypes",
  type = "boxplot",
  size.point = 1
)

## ----barError_newModels-------------------------------------------------------
barErrorPlot(DDLSToy, error = "MAE", by = "CellType")

## ----corr1_newModels----------------------------------------------------------
corrExpPredPlot(
  DDLSToy,
  color.by = "CellType",
  size.point = 1,
  corr = "both"
)

## ----corr2_newModels----------------------------------------------------------
corrExpPredPlot(
  DDLSToy,
  color.by = "CellType",
  facet.by = "CellType",
  size.point = 1, 
  filter.sc = FALSE,
  corr = "both"
)

## ----corr3_newModels----------------------------------------------------------
corrExpPredPlot(
  DDLSToy,
  color.by = "CellType",
  facet.by = "nCellTypes",
  size.point = 1,
  corr = "both"
)

## ----bland1_newModels---------------------------------------------------------
blandAltmanLehPlot(
  DDLSToy, 
  color.by = "CellType",
  log.2 = FALSE,
  size.point = 1,
  filter.sc = TRUE,
  density = TRUE,
)

## ----bland2_newModels---------------------------------------------------------
blandAltmanLehPlot(
  DDLSToy, 
  color.by = "nCellTypes",
  facet.by = "nCellTypes",
  log.2 = FALSE,
  size.point = 1,
  filter.sc = TRUE,
  density = FALSE
)

## -----------------------------------------------------------------------------
countsBulk <- matrix(
  stats::rpois(100, lambda = sample(seq(4, 10), size = 100, replace = TRUE)), 
  nrow = 40, ncol = 15, 
  dimnames = list(paste0("Gene", seq(40)), paste0("Bulk", seq(15)))
)

## ----deconvoluteNewBulk_newModels---------------------------------------------
suppressMessages(library(SummarizedExperiment, quietly = TRUE))
seExample <- SummarizedExperiment(assay = list(counts = countsBulk))

DDLSToy <- loadDeconvData(
  object = DDLSToy,
  data = seExample, 
  name.data = "Simulated.example"
)

## ----resultsBarPlot_newModels, warning=FALSE----------------------------------
DDLSToy <- deconvDigitalDLSorterObj(
  object = DDLSToy, 
  name.data = "Simulated.example",
  normalize = TRUE,
  scaling = "standarize",
  verbose = FALSE
)
## plot results
barPlotCellTypes(
  DDLSToy, name.data = "Simulated.example", 
  rm.x.text = TRUE, color.line = "black"
)

## ----saveRDS, eval=FALSE------------------------------------------------------
#  ## this code will not be run
#  saveRDS(object = DDLSToy, file = "valid/path")

## ----saveHDF5Model, eval=FALSE------------------------------------------------
#  ## this code will not be run
#  saveTrainedModelAsH5(DDLSToy, file.path = "valid/path")
#  DDLSToy <- loadTrainedModelFromH5(DDLSToy)

