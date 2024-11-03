possiblePAM <- findPAM(seq, pam = pam)
if(nrow(possiblePAM) == 0){
  "No PAM Found"
}
