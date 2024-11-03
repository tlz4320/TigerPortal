setwd("~/data/project/ear_project/gene_therapy_ll/batch2/")
load("~/data/project/ear_project/gene_therapy_ll/Result/sgCmp132.rda")
load("~/data/project/ear_project/gene_therapy_ll/Result/merged_indel_table_first_second.Rda")


###计算各个sgRNA相关性在ABC样本两两之间的结果
max_region <- c(-200, 30)
total_sg_cor <- lapply(sgCmp132$id2, function(y){
  sel_sg <- total_indel_table_first_second_rev[[y]]
  if(length(sel_sg) < 2){
    print(y)
  }
  cmb <- combn(length(sel_sg), 2)
  cor_name <- lapply(1 : ncol(cmb), function(x){
    paste(c("A", "B", "C")[unlist(cmb[,x])], collapse = "-")
  })
  cor_res <- lapply(1 : ncol(cmb), function(x){
    ids <- unlist(cmb[,x])
    sel_sg1 <- sel_sg[[ids[1]]]
    sel_sg2 <- sel_sg[[ids[2]]]
    max_del <- min(c(sel_sg1[,1], sel_sg2[,1]))
    max_del <- max(c(max_del, max_region[1]))
    max_in <- max(c(sel_sg1[,1], sel_sg2[,1]))
    max_in <- min(c(max_in, max_region[2]))
    region <- max_del : max_in
    tmp_region <- region[!region %in% sel_sg1[,1]]
    if(length(tmp_region) != 0){
      sel_sg1 <- data.frame(rbind(sel_sg1, data.frame(indel_size = tmp_region, fq = 0)))
    }
    tmp_region <- region[!region %in% sel_sg2[,1]]
    if(length(tmp_region) != 0){
      sel_sg2 <- data.frame(rbind(sel_sg2, data.frame(indel_size = tmp_region, fq = 0)))
    }
    sel_sg1[sel_sg1[,1] == 0, 2] <- 0
    sel_sg2[sel_sg2[,1] == 0, 2] <- 0
    common_indel <- intersect(sel_sg1[,1], sel_sg2[,1])
    sel_sg1 <- sel_sg1[sel_sg1[,1] %in% common_indel,]
    sel_sg2 <- sel_sg2[sel_sg2[,1] %in% common_indel,]
    sel_sg1 <- sel_sg1[order(sel_sg1[,1]),]
    sel_sg2 <- sel_sg2[order(sel_sg2[,1]),]
    # rm0 <- sel_sg1[,2] == 0 & sel_sg2[,2] == 0
    #cor.test(sel_sg1[!rm0,2], sel_sg2[!rm0,2])
    cor(sel_sg1[,2], sel_sg2[,2])
  })
  data.frame(name = unlist(cor_name), cor = unlist(cor_res), sg = y)
})
total_sg_cor <- data.frame(do.call(rbind, total_sg_cor))
total_sg_cor$sg <- factor(total_sg_cor$sg, 
                          levels = gtools::mixedsort(unique(total_sg_cor$sg)))

sel <- c("Sg_19_131", "Sg_13_89", "Sg_24_234", "Sg_24_233", "Sg_24_271", "Sg_24_272", "Sg_24_236", "Sg_24_235", 
         "Sg_10_67", "Sg_19_129")

total_sg_cor[total_sg_cor$sg %in% sel,]
