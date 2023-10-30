###Confirm of tRNA anti-codon location
for (i in 1:nrow(nc_trna)){
  string <- nc_trna[i,3]
  pattern <- nc_trna[i,2]
  tmp<-unlist(sapply(seq_len(nchar(string)), function(x) {
    if(substr(string, x, x + 3 - 1) == pattern){
      return (x)}
  }))
  tmp<-tmp[tmp>32&tmp<37]
  if(length(tmp)>1){
    print(i)
  }
}	
for (i in 1:nrow(nc_trna)){
  string <- nc_trna[i,3]
  pattern <- nc_trna[i,2]
  tmp<-unlist(sapply(seq_len(nchar(string)), function(x) {
    if(substr(string, x, x + 3 - 1) == pattern){
      return (x)}
  }))
  tmp<-tmp[tmp>32&tmp<37] 
  nc_trna$antic_pos[i]<-tmp[1]
}	
nc_trna[c(13,40,59,81),"antic_pos"]<-rep(34,4)

for (i in 1:nrow(mt_trna)){
  string <- mt_trna[i,3]
  pattern <- mt_trna[i,2]
  tmp<-unlist(sapply(seq_len(nchar(string)), function(x) {
    if(substr(string, x, x + 3 - 1) == pattern){
      return (x)}
  }))
  tmp<-tmp[tmp>20&tmp<40]
  mt_trna$antic_pos[i]<-tmp[1]
}	
mt_trna[c(2,4,11,15),"antic_pos"]<-c(rep(31,3),26)



