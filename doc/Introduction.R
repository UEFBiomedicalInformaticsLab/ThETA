## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- load_rpack, cache.lazy = TRUE, echo = FALSE, message = FALSE, warning = FALSE----
library(dplyr)

## ---- load_theta, cache.lazy = TRUE, echo = TRUE, message = FALSE--------
library(ThETA)
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

## ---- t2d_genes, cache.lazy = TRUE, echo = TRUE--------------------------
T2DM_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:0001360")]]

## ---- t2d_tzs, cache.lazy = TRUE, echo = TRUE----------------------------
T2DM_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:0001360"),]

## ---- t2d_sco, cache.lazy = TRUE, echo = TRUE, results = 'hide'----------
T2DM_Tscores <- tissue.specific.scores(T2DM_genes$entrez[1:2], 
                                        ppi_network = ppi_strdb_700, 
                                        directed_network = FALSE, 
                                        tissue_expr_data = gtexv7_zscore,
                                        dis_relevant_tissues = T2DM_rel_tissue_scores, 
                                        W = centrality_score$borda.disc, 
                                        cutoff = 4, verbose = TRUE)

## ---- cache = TRUE, echo = FALSE, message = FALSE------------------------
knitr::kable(T2DM_Tscores[1:5,]) %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F) 

## ---- t2d_ord, cache.lazy = TRUE, echo = TRUE----------------------------
T2DM_top50 <- T2DM_Tscores[order(T2DM_Tscores$avg_tissue_score, 
                                  decreasing = TRUE)[1:50],]

## ---- t2d_plot1, fig.width=5, fig.height=3, cache.lazy = TRUE, echo = TRUE----
library(ggplot2)
library(reshape)
data_t2d50 <- reshape::melt(as.matrix(T2DM_top50), id = 0)
colnames(data_t2d50) <- c("EntrezID", "Tissue", "EfficacyScore")
ggplot(data_t2d50, aes(x = Tissue, y = EfficacyScore, fill = Tissue)) +
        geom_boxplot(alpha = 0.7) +
        ggtitle("Boxplot of the efficacy scores for tissues") +
        theme_bw() +
        theme(plot.title = element_text(size = 8, family = "Tahoma", face = "bold"),
                text = element_text(size = 7, family = "Tahoma"),
                axis.title = element_text(face="bold"),
                axis.text.x=element_text(size = 7)) 

## ---- t2d_plot2, fig.width=5, fig.height=3, cache.lazy = TRUE, echo = TRUE----
library(ggplot2)
library(reshape)
data_t2d5 <- reshape::melt(as.matrix(T2DM_top50[1:5,-ncol(T2DM_top50)]), id = 0)
colnames(data_t2d5) <- c("EntrezID", "Tissue", "EfficacyScore")
data_t2d5$EntrezID = factor(data_t2d5$EntrezID)
ggplot(data_t2d5, aes(x = EntrezID, y = EfficacyScore, fill = EntrezID)) +
        geom_bar(stat='identity', alpha = 0.7) +
        ggtitle("Barplot comparing the efficacy scores of genes acrsso significant tissues") +
        facet_wrap(~Tissue) +
        theme_bw() + 
        theme(plot.title = element_text(size = 8, family = "Tahoma", face = "bold"),
                    text = element_text(size = 7, family = "Tahoma"),
                    axis.title = element_text(face="bold"),
                    axis.text.x=element_text(size = 7)) 

## ---- load_gpl, cache.lazy = TRUE, echo = TRUE---------------------------
data(geo_gene_sets)

## ---- get_mods, cache.lazy = TRUE, echo = TRUE---------------------------
modulation_scores <- modulation.score(geneSets = geo_gene_sets)

## ---- modscores_tab, cache = TRUE, echo = FALSE--------------------------
knitr:::kable(modulation_scores[1:5, ]) %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F) 

## ---- load_conf, cache.lazy = TRUE, echo = TRUE--------------------------
enrichr_to_efo <- read.csv(system.file("conversion_enrichr_efo.csv", 
                                        package = "ThETA"), row.names = 1,
                            stringsAsFactors = F)
modulation_scores$disease.id <- enrichr_to_efo[modulation_scores$disease.id,'disease.id']

## ---- map_genes, cache.lazy = TRUE, echo = TRUE, message = FALSE, warning = FALSE----
library(org.Hs.eg.db)
modulation_scores$target.entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, modulation_scores$target.id,'ENTREZID','SYMBOL')
modulation_scores <- modulation_scores[modulation_scores$disease.id != '' &
                                       !is.na(modulation_scores$target.entrez),]

## ---- max_pert, cache.lazy = TRUE, echo = TRUE, message = FALSE----------
library(data.table)
modul_score <- data.table::as.data.table(modulation_scores)
modul_score <- as.data.frame(modul_score[, .SD[which.max(modscore)], 
                                           by=list(disease.id, target.entrez)])

## ---- t2d_mods, cache.lazy = TRUE, echo = TRUE---------------------------
T2DM_Mscores = data.frame(modul_score[modul_score$disease.id=='EFO:0001360', 
                                      c("target.entrez", "modscore")], row.names = 1)

## ---- t2d_mods_tab, cache = TRUE, echo = FALSE---------------------------
knitr:::kable(head(T2DM_Mscores)) %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F) %>%
  kableExtra::scroll_box(width = "200px", height = "300px")

## ---- int_scores, cache.lazy = TRUE, echo = TRUE-------------------------
common_t2d_genes <- intersect(rownames(T2DM_Mscores), rownames(T2DM_Tscores)) 
T2DM_Iscores <- data.frame("Mscore" = T2DM_Mscores[common_t2d_genes,], 
                            "TSEscore" = T2DM_Tscores[common_t2d_genes,],
                            row.names = common_t2d_genes)

## ---- int_scores_tab, cache = TRUE, echo = FALSE-------------------------
knitr:::kable(T2DM_Iscores[1:5, c(1:3, ncol(T2DM_Iscores))]) %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F) %>%
  kableExtra::scroll_box(width = "700px", height = "200px")

## ---- query_ot, cache.lazy = TRUE, echo = TRUE---------------------------
server <- 'https://platform-api.opentargets.io/v3/platform'
endpoint_prmtrs <- '/public/association/filter'
optional_prmtrs <- '?size=10000&disease=EFO_0001360&fields=disease.id&fields=target.gene_info.symbol&fields=association_score.overall&fields=disease.efo_info.label'
uri <- paste(server,endpoint_prmtrs,optional_prmtrs,sep='')

## ---- op_get, cache.lazy = TRUE, echo = TRUE, message = FALSE------------
if("httr" %in% rownames(installed.packages()) == FALSE) {install.packages("httr")}
if("jsonlite" %in% rownames(installed.packages()) == FALSE) {install.packages("jsonlite")}
library(httr)
library(jsonlite)

get_association_json <- httr::content(httr::GET(uri),'text')
get_association_usable <- jsonlite::fromJSON(get_association_json, flatten = TRUE)

OT_score <- get_association_usable$data[,c(2:3,1,4)]
OT_score$disease.id <- gsub('_',':',OT_score$disease.id)
colnames(OT_score)[c(1,4)] <- c('target.id', 'disease.name')

# remove duplicated gene symbols
OT_score = OT_score[-which(duplicated(OT_score$target.id)),]

## ---- op_map, cache.lazy = TRUE, echo = TRUE, message = FALSE, warning = FALSE----
library(org.Hs.eg.db)
OT_score$target.entrez <- AnnotationDbi::mapIds(org.Hs.eg.db,OT_score$target.id,'ENTREZID','SYMBOL')
OT_score <- OT_score[!is.na(OT_score$target.entrez),]

## ---- op_tab, cache = TRUE, echo = FALSE---------------------------------
knitr:::kable(head(OT_score)) %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F) %>%
  kableExtra::scroll_box(width = "700px", height = "150px")

## ---- add_op_scores, cache.lazy = TRUE, echo = TRUE----------------------
all_scores <- base::merge(OT_score, T2DM_Iscores, by.x = "target.entrez", by.y = "row.names", all = TRUE)

## ---- int_op_scores1, cache.lazy = TRUE, echo = TRUE---------------------
T2DM_allsc <- integrate.scores(all_scores, c("association_score.overall",
                                             "Mscore", 
                                             "TSEscore.avg_tissue_score"))
T2DM_allsc <- T2DM_allsc[order(T2DM_allsc$HS, decreasing = TRUE),]
rownames(T2DM_allsc) <- T2DM_allsc[,1]

# let's semplify the final table of the disease-gene association scores
tab_score <- T2DM_allsc[,c("target.id","association_score.overall", "Mscore", 
                           "TSEscore.avg_tissue_score", "HS","MAX")]
colnames(tab_score)[1:4] <- c("GeneTarget","OTScore","ModulationScore","TissueEfficacyScore")

## ---- op_int, cache = TRUE, echo = FALSE---------------------------------
knitr:::kable(head(tab_score)) %>%
  kableExtra::kable_styling(bootstrap_options = "striped", full_width = F) %>%
  kableExtra::scroll_box(width = "700px", height = "200px")

## ---- cache.lazy = TRUE, eval = FALSE------------------------------------
#  library(shiny)
#  library(visNetwork)
#  library(org.Hs.eg.db)
#  
#  visualize.graph(tissue_scores = T2DM_Tscores,
#                   disease_genes =T2DM_genes$entrez[1:5],
#                   ppi_network = ppi_strdb_700,
#                   tissue_expr_data = gtexv7_zscore,
#                   top_targets = rownames(T2DM_top50)[1:5],
#                   db='BP')

## ---- tsrwr_dat, cache.lazy = TRUE, echo = TRUE, warning=FALSE, message = FALSE----
tsrwr = build.tissue.specific.networks(tissue_scores = T2DM_Tscores, disease_genes = T2DM_genes$entrez,
                                       ppi_network = ppi_strdb_700, tissue_expr_data = gtexv7_zscore, 
                                       top_targets = rownames(T2DM_top50)[1:5], verbose = FALSE)


## ---- ora_dat, cache.lazy = TRUE, echo = TRUE, message = FALSE, warning = FALSE----
T2D_ora_data_shp = generate.ora.data(tsrwr$shp[[1]], databases = "KEGG")
T2D_ora_data_rwr = generate.ora.data(tsrwr$rwr, databases = "KEGG")

## ---- plot_rwr, cache.lazy = TRUE, echo = TRUE, message = FALSE, warning = FALSE----
T2D_ora_plot_rwr = generate.ora.plots(T2D_ora_data_rwr, set_plots = c("dotplot","cnetplot"), 
                                      showCategory = 5, font_size = 10)

## ---- mult_plot_rwr, fig.width=7, fig.height=9, cache.lazy = TRUE, echo = TRUE, warning = FALSE----
figure <- ggpubr::ggarrange(plotlist = T2D_ora_plot_rwr[1:2], nrow = 2, ncol = 1, 
                            common.legend = TRUE, legend = "bottom", labels=names(T2D_ora_plot_rwr)[1:2])
figure

## ---- novelty_plot, fig.width=6, fig.height=2, cache.lazy = TRUE, echo = TRUE, warning = FALSE, message = FALSE----
library(org.Hs.eg.db)
T2D_pmd_plot_top = novelty.plots(rownames(T2DM_allsc)[c(1:5)], orgdb = org.Hs.eg.db, font_size = 14, pubmed = c(2010,2018))
T2D_pmd_plot_top

