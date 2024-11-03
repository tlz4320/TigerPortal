library(seqinr)
totalRNA <- read.fasta(file = "~/database/annovar/humandb/hg38_refGeneMrna.fa", as.string = TRUE, seqtype = "DNA")
RNA_region <- read.table("~/database/annovar/humandb/hg38_refGene.txt", sep="\t", quote = '"')
RNA_region <- RNA_region[RNA_region$V8 > RNA_region$V7,]
RNA_region <- RNA_region[RNA_region$V14 == "cmpl",]
RNA_region$Len <- unlist(lapply(1 : nrow(RNA_region), function(i){
  sum(as.numeric(unlist(strsplit(RNA_region$V11[i], ","))) - as.numeric(unlist(strsplit(RNA_region$V10[i], ","))))
}))
RNA_region <- RNA_region[RNA_region$V2 %in% names(totalRNA),]
tmp <- unlist(lapply(1 : nrow(RNA_region), function(i){
  id <- RNA_region$V2[i]
  str_length(as.character(totalRNA[[id]])) == RNA_region$Len[i]
}))
RNA_region <- RNA_region[tmp,]
tmp <- mclapply(1 : nrow(RNA_region), function(i){
  id <- RNA_region$V2[i]
  seq <- as.character(totalRNA[[id]])
  strand <- RNA_region$V4[i]
  len <- RNA_region$Len[i]
  if(strand == "+"){
    starts <- as.numeric(unlist(strsplit(RNA_region$V10[i], ",")))
    ends <- as.numeric(unlist(strsplit(RNA_region$V11[i], ",")))
    
    cds_start <- RNA_region$V7[i]
    cds_end <- RNA_region$V8[i]
    rm_left <- 0
    rm_right <- 0
    for(step in 1 : length(starts)){
      if(ends[step] < cds_start){
        rm_left <- rm_left + ends[step] - starts[step]
      } else if(starts[step] <= cds_start){
        rm_left <- rm_left + cds_start - starts[step]
      }
      if(starts[step] > cds_end){
        rm_right <- rm_right + ends[step] - starts[step]
      } else if(ends[step] > cds_end){
        rm_right <- rm_right + ends[step] - cds_end
      }
    }
    if(((len - rm_left - rm_right) %% 3) != 0){
      #存在一些很奇怪的RNA不知道为什么数据库里面的位置是错的
      return(data.frame(left = "WTF", cds = "WTF", right="WTF"))
    }
    cds <- str_sub(seq, rm_left + 1, -rm_right - 1)
    left <- str_sub(seq, 1, rm_left)
    if(rm_right == 0){
      right <- ""
    } else {
      right <- str_sub(seq, -rm_right)
    }
    return(data.frame(left = left, cds = cds, right = right))
  } else {
    starts <- as.numeric(unlist(strsplit(RNA_region$V10[i], ",")))
    ends <- as.numeric(unlist(strsplit(RNA_region$V11[i], ",")))
    
    cds_start <- RNA_region$V7[i]
    cds_end <- RNA_region$V8[i]
    rm_left <- 0
    rm_right <- 0
    for(step in length(starts) : 1){
      if(starts[step] > cds_end){
        rm_left <- rm_left + ends[step] - starts[step]
      } else if(ends[step] > cds_end){
        rm_left <- rm_left + ends[step] - cds_end
      }
      if(ends[step] < cds_start){
        rm_right <- rm_right + ends[step] - starts[step]
      } else if(starts[step] <= cds_start){
        rm_right <- rm_right + cds_start - starts[step]
      }
    }
    if(((len - rm_left - rm_right) %% 3) != 0){
      #存在一些很奇怪的RNA不知道为什么数据库里面的位置是错的
      return(data.frame(left = "WTF", cds = "WTF", right="WTF"))
    }
    cds <- str_sub(seq, rm_left + 1, -rm_right - 1)
    left <- str_sub(seq, 1, rm_left)
    if(rm_right == 0){
      right <- ""
    } else {
      right <- str_sub(seq, -rm_right)
    }
    
    return(data.frame(left = left, cds = cds, right = right))
  }
  return(data.frame(left = "WTF", cds = "WTF", right="WTF"))
}, mc.cores = 10)
tmp <- data.frame(do.call(rbind, tmp))

RNA_region$Seq <- tmp$cds
RNA_region$Left <- tmp$left
RNA_region$Right <- tmp$right
RNA_region <- RNA_region[RNA_region$Seq != "WTF" & RNA_region$Seq != "WTF2",]
# save(RNA_region, file="~/data/project/ear_project/gene_therapy_ll/RNA_region.rda")
ToNX::write_tb(RNA_region, file="~/RNA_region.txt")


load("~/data/project/ear_project/gene_therapy_ll/RNA_region.rda")
dbsnp <- read.table("~/Indel1_mutation.txt")
colnames(dbsnp) <- c("Chr", "pos", "pos2", "Ref", "Alt", "ID")
dbsnp$Chr <- paste0("chr", dbsnp$Chr)
ref_len <- str_length(dbsnp$Ref)
dbsnp <- dbsnp[-which(is.na(ref_len)),]
dbsnp <- dbsnp[str_length(dbsnp$Ref) == 1,]


###注意一下RNA_region的坐标是0开始的  dbsnp是1开始的  验证过了
isInCDS <- mclapply(1 : nrow(dbsnp), function(i){
  chr <- dbsnp$Chr[i]
  start <- dbsnp$pos[i]
  end <- dbsnp$pos2[i]
  ref <- dbsnp$Ref[i]
  alt <- dbsnp$Alt[i]
  sel <- RNA_region[RNA_region$V3 == chr & 
                      RNA_region$V7 <= start & RNA_region$V8 > start,]
  if(nrow(sel) == 0){
    return("")
  }
  possible_res <- c()
  for(j in 1 : nrow(sel)){
    starts <- as.numeric(unlist(strsplit(sel$V10[j], ",")))
    ends <- as.numeric(unlist(strsplit(sel$V11[j], ",")))
    for(step in 1 : length(starts)){
      if(start > starts[step] & start <= ends[step]){
        possible_res <- c(possible_res, sel$V2[j])
        break
      }
        
    }
  }
  return(paste0(possible_res, collapse = ","))
}, mc.cores = 20)
table(unlist(isInCDS) == "")

save(dbsnp, file="~/data/project/ear_project/gene_therapy_ll/indel_snp.rda")

mutation_info <- read.table("~/Indel1_mutation_inCDS_anno.txt", sep="\t")

id <- paste(mutation_info$V6, mutation_info$V9)
mutation_info <- mutation_info[!duplicated(id),]
mutation_info$V7 <- str_to_upper(mutation_info$V7)
# save(mutation_info, file='~/data/project/ear_project/gene_therapy_ll/website_data.Rda')



largest_cds <- lapply(split(RNA_region, RNA_region$V13), function(x){
  x[which.max(str_length(x$Seq)),]
})
largest_cds <- data.frame(do.call(rbind, largest_cds))
# save(largest_cds, file="~/data/project/ear_project/gene_therapy_ll/largest_cds.rda")
