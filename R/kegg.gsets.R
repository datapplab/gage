kegg.gsets<-function(species="hsa", id.type="kegg", check.new=FALSE){
  id.type=tolower(id.type)
  egid=NA
  use.egid=FALSE
  if(species!="ko") {
    species.data=kegg.species.code(species, na.rm=T, code.only=F)
    species=species.data[1]
    use.egid=species.data[2]=="1"
    egid=species.data[4]
  } else if(id.type!="kegg"){
    id.type="kegg"
  }

  path.list <- keggLink("pathway",species)
  genes=gsub(paste(species,":",sep=""), "", names(path.list))

  if(id.type=="entrez" & !use.egid & !is.na(egid)){
    gid.map=keggConv("ncbi-geneid",species)
    kegg.ids=gsub(paste(species, ":", sep=""), "", names(gid.map))
    eg.ids=gsub("ncbi-geneid:", "", gid.map)
    id.idx=match(genes, kegg.ids)
    genes[!is.na(id.idx)]=eg.ids[id.idx[!is.na(id.idx)]]
  } else if (!(id.type %in% c("kegg", "entrez")) | (id.type=="entrez" & !use.egid)) {
    na.msg=sprintf("KEGG Gene ID instead of '%s' is used!", id.type)
    message("Note: ", na.msg)
  }

  paths=gsub(paste("path:",species,sep=""), "", path.list)
  kg.sets=split(genes, paths)

  if(!exists("khier")){
    if(check.new){
      si=try(load(url("https://pathview.uncc.edu/data/khier.rda")))
      if(class(si)[1]=="try-error") data(khier, package="gage")
    } else data(khier, package="gage")
  }

  kh.idx=match(names(kg.sets), substr(khier[,3], 1,5))
  kg.sets=kg.sets[!is.na(kh.idx)]
  kh.idx=kh.idx[!is.na(kh.idx)]
  ks.names=khier[kh.idx,3]

  aterms=unique(khier[,1])
  signals=khier[which(khier[,1] %in% aterms[2:5]),3]
  metabs=khier[which(khier[,1] == aterms[1]),3]
  diseases=khier[which(khier[,1] == aterms[6]),3]

  sig.idx=which(ks.names %in% signals)
  met.idx=which(ks.names %in% metabs)
  dise.idx=which(ks.names %in% diseases)
  sigmet.idx=c(sig.idx,met.idx)
  names(kg.sets)=paste(species, ks.names, sep="")
  res=list(kg.sets=kg.sets, sigmet.idx=sigmet.idx, sig.idx
    =sig.idx, met.idx=met.idx, dise.idx=dise.idx)
  return(res)
}
