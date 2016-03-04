define_profiles <- function(){
  bg <- 1.5
  width <- 200
  peak.base <- 20
  hts <- c(bg, 3, 10, 60)
  profiles <- matrix(nrow=width, ncol=length(hts))
  for(i in 1:length(hts)){
    w <- peak.base/2
    w1 <- width/2 - peak.base/2
    w2 <- w1 + peak.base/2
    w3 <- w2 + peak.base/2
    w4 <- w3 + w1
    profiles[1:w1,i] <- bg
    profiles[(w3+1):w4, i] <- bg
    profiles[(w1+1):w2, i] <- bg + (1:w)*(hts[i]-bg)/w
    profiles[(w2+1):w3, i] <- hts[i]- (1:w)*(hts[i]-bg)/w
  }
  return(profiles)

}
