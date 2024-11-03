# sg <- "GCTACAGATTTAGGCCTGAA"
# edit_table <- read.table("~/code/inDelphi-model/tmp_test.csv", sep = ",", header = T)
fix_indel1_func <- function(edit_table, sg){
  leftnt2 <- str_sub(sg, -5, -4)
  rightnt2 <- str_sub(sg, -3, -2)
  leftnt <- str_sub(leftnt2, 2, 2)
  rightnt <- str_sub(rightnt2, 1, 1)
  #only fix when left and right nt is different
  if(leftnt != rightnt){
    del1_row <- which(edit_table$Category == "del" & edit_table$Length == 1)
    totalpct <- sum(edit_table$Predicted.frequency[del1_row])
    if(leftnt2 %in% names(del1_besides_2nt_stat)){
      if(rightnt2 %in% names(del1_besides_2nt_stat[[leftnt2]])){
        pctlist <- del1_besides_2nt_stat[[leftnt2]][[rightnt2]]
      } else {
        pctlist <- del1_besides_1nt_stat[[leftnt]][[rightnt]]
      }
    } else {
      pctlist <- del1_besides_1nt_stat[[leftnt]][[rightnt]]
    }
    pctlist <- pctlist / sum(pctlist) * totalpct
    for(i in del1_row){
      if(edit_table$Genotype.position[i] == 0){
        edit_table$Predicted.frequency[i] = pctlist[leftnt]
      } else {
        edit_table$Predicted.frequency[i] = pctlist[rightnt]
      }
    }
  }
  ins1_row <- which(edit_table$Category == "ins" & edit_table$Length == 1)
  totalpct <- sum(edit_table$Predicted.frequency[ins1_row])
  if(leftnt2 %in% names(ins1_besides_2nt_stat)){
    if(rightnt2 %in% names(ins1_besides_2nt_stat[[leftnt2]])){
      pctlist <- ins1_besides_2nt_stat[[leftnt2]][[rightnt2]]
    } else {
      pctlist <- ins1_besides_1nt_stat[[leftnt]][[rightnt]]
    }
  } else {
    pctlist <- ins1_besides_1nt_stat[[leftnt]][[rightnt]]
  }
  pctlist <- pctlist / sum(pctlist) * totalpct
  for(i in ins1_row){
    edit_table$Predicted.frequency[i] = pctlist[edit_table$Inserted.Bases[i]]
  }
  edit_table
}
