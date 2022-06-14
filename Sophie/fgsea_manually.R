library(ggplot2)
library(openxlsx)
library(dplyr)
library(org.Hs.eg.db)
library(fgsea)
library(tidyr)
library(tibble)
library(EnhancedVolcano)


# si besoin, changer data pour uniformiser les FC 
file <- "data/DE/SN-comparaison-1.xlsx"

wb_ct_1 <- loadWorkbook(file)
sheets_name <- names(wb_ct_1)
sheets_name
sheets_name <- sheets_name[-1]
sheets_name
# change FC and LogFC for NT vs...
ind_inv <- c(1,3,4,6)
new_wb = createWorkbook()
i = 1
for (s in sheets_name){
    dat = read.xlsx(wb_ct_1, sheet = s)
    dat$Fold.Change = as.numeric(dat$Fold.Change)
    dat[,'Log.(Fold.Change)'] = as.numeric(dat[,'Log.(Fold.Change)'])
    dat = na.omit(dat)
    if (i %in% ind_inv){
        dat$Fold.Change = 1/dat$Fold.Change
        dat[,'Log.(Fold.Change)'] = -1*dat[,'Log.(Fold.Change)']
        split_name = strsplit(s, "_", fixed = T)
        sub_split_name = strsplit(split_name[[1]][2], " ")
        new_sheet_name = paste(split_name[[1]][1], "_",
                                sub_split_name[[1]][3]," ",
                                sub_split_name[[1]][2]," ",
                                sub_split_name[[1]][1], sep = "")
        addWorksheet(new_wb, sheetName = new_sheet_name)
        # addWorksheet(new_wb, sheetName = s)
        writeData(new_wb, sheet = new_sheet_name, x = dat)
        # writeData(new_wb, sheet = s, x = dat)
    }
    else{
        addWorksheet(new_wb, sheetName = s)
        writeData(new_wb, sheet = s, x = dat)
    }
    i = i+1
}
split = strsplit( file, split = ".", fixed = T)

saveWorkbook(new_wb, file = paste(split[[1]][1],"_cor.xlsx", sep = ""),
            overwrite = T)




# Analyse GSEA et Volcano plot

# utiiser les xlsx "cor"
file <- "data/DE/SN-comparaison-2_cor.xlsx"

dir_name <- "SN_comp_2"
dir.create(dir_name)
wb_ct_1 <- loadWorkbook(file)
sheets_name <- names(wb_ct_1)
# sheets_name <- sheets_name[-1]
# sheets_name_new <- gsub(" ", "_", sheets_name)
sheets_name
# ind_inv <- c(1,3,4,6) #indices des fc à inverser
# recuperation des fc pour chaque comparaison
# vecteur indices des données à inverser : ind_inv = c()
list_comp <- list()
# hs <- org.Hs.eg.db
# i <- 1
for (s in sheets_name) {
    dat <- read.xlsx(wb_ct_1, sheet = s)
    l <- as.numeric(dat[, 5])
    names(l) <- dat[, 3]
    l <- na.omit(l)
    # if (i %in% ind_inv) {
    #     print(s)
    #     l <- -1 * l
    # }
    l <- sort(l, decreasing = T)
    # switch symbol to entrezid
    # e_id = select(hs, keys = names(l), columns = "ENTREZID", keytype = "SYMBOL")
    # e_id = e_id[!(duplicated(e_id[,1])),]
    # names(l) = e_id[,2]
    list_comp[[s]] <- l[1:1000]
    # i <- + 1
}
names(list_comp) <- gsub(" ", "_", names(list_comp))
# hs <- org.Hs.eg.db
pathways <- list.files("GSEA_sig/")
pathways



for (p in pathways) {
    name_pathways_list <- p
    print(name_pathways_list)
    pathways.hallmark <- gmtPathways(paste("GSEA_sig/", name_pathways_list, sep = ""))
    names(pathways.hallmark) = substring(names(pathways.hallmark), first = 1, last = 50)

    for (n in names(list_comp)) {
        if (!dir.exists(paste(dir_name, "/", n, sep = ""))) {
            dir.create(paste(dir_name, "/", n, sep = ""))
        }
        fgseaRes <- fgsea(pathways = pathways.hallmark, stats = list_comp[[n]])
        fgseaRes <- na.omit(fgseaRes)
        fgseaResTidy <- fgseaRes %>%
            as_tibble() %>%
            arrange(desc(NES))
        l = dim(fgseaResTidy)[1]
        g = ggplot(fgseaResTidy[c(1:20, (l - 20):l), ], aes(reorder(pathway, NES), NES)) +
        geom_col(aes(fill = padj < 0.05)) +
        coord_flip() +
        labs(x = "", y = "Normalized Enrichment Score") +
        ggtitle(paste(n, " ", p)) +
        theme(
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 7),
            plot.title.position = "plot",
            legend.position = "none"
            )

            png(paste(dir_name, "/", n, "/", name_pathways_list, ".png", sep = ""))
            print(g)
            dev.off()
        }

    }


# volcano

for (s in sheets_name) {
    dat <- read.xlsx(wb_ct_1, sheet = s)
    df = data.frame("LogFC" = as.numeric(dat[,'Log.(Fold.Change)']))
    df[,"pval"]= dat[,'p-value']
    rownames(df) = dat[,"Peak.Name"]
    df <- na.omit(df)
    v = EnhancedVolcano(toptable = df, lab = rownames(df),
                    x = "LogFC", y = "pval",
                    FCcutoff = c(-0.3,0.3), title=s)
    ns = gsub(" ", "_", s)
    png(paste(dir_name, "/", ns, "/volcano_", s ,".png", sep = ""))
    print(v)
    dev.off()
}
