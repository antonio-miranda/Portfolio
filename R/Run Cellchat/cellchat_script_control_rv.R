print(getwd())

library(reticulate)

use_condaenv("personal")

library(CellChat)

library(patchwork)

options(stringsAsFactors = FALSE)

options(repr.plot.width = 10, repr.plot.height = 10, repr.plot.res = 300)

ad <- import("anndata", convert = FALSE)

ad_object <- ad$read_h5ad('/rds/general/user/aalmeid2/projects/cardiac_single_cell_biology/live/DCM_project/global_control_rv_oEC7.h5ad')

data=as.matrix(py_to_r(ad_object$X))

data.input <- t(data)

rownames(data.input) <- rownames(py_to_r(ad_object$var))
colnames(data.input) <- rownames(py_to_r(ad_object$obs))

meta.data <- py_to_r(ad_object$obs)
meta <- meta.data

sub <- c(AD="AD",CM="CM",EC="EC7.0",FB="FB",Lymphoid="LY",Mast="Mast",Mural="Mu",Myeloid="MY",N="NC")

meta$cell_type <- as.character(sub[meta$cell_type])

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cell_type")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 

future::plan("multiprocess", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

cellchat <- computeCommunProb(cellchat)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)

saveRDS(cellchat, file = "cellchat_control_lv_oEC7.rds")

cellchat <- computeCommunProbPathway(cellchat)

cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = "cellchat_control_lv_oEC7.rds")