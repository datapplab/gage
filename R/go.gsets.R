go.gsets <-  function (species = "human", pkg.name = NULL, id.type="eg", keep.evidence=FALSE)
  {
    if (is.null(pkg.name)){
      species=tolower(species)
      data(bods, package="gage")
      bods[,"species"]=tolower(bods[,"species"])
      idx=which(bods[,"species"] == species)
      if(length(idx)!=1) stop("bad species value")
      pkg.name = bods[idx, "package"]
      id.type=toupper(bods[idx, "id.type"])
    } else{
      species=tolower(species)
      id.type=toupper(id.type)
      warn.msg=paste("User specified annotation package, please make sure:\n",
        "-package is ready for loading\n-'species' and 'id.type' are consistent!", sep="")
      message(warn.msg)
    }

    pkg.on = require(pkg.name, character.only = TRUE)
    if (!pkg.on) {
      if (!requireNamespace("BiocManager", quietly=TRUE))
          install.packages("BiocManager")
      BiocManager::install(pkg.name, suppressUpdates =TRUE)
      pkg.on = require(pkg.name, character.only = TRUE)
      if (!pkg.on)
        stop(paste("Fail to install/load gene annotation package ",
                   pkg.name, "!", sep = ""))
    }

    pkg.name = gsub("[.]db", "", pkg.name)
    gid.msg=sprintf("Gene ID type for '%s' is: '%s'", species, id.type)
    message(gid.msg)
    bimap.go = eval(as.name(paste(pkg.name, "GO2ALL", id.type, "S", sep = "")))
    go= AnnotationDbi::as.list(bimap.go)
    #library(GO.db)
    #requireNamespace("GO.db", quiet=TRUE)
    go.names=sapply(AnnotationDbi::mget(names(go), GO.db::GOTERM), function(x) x@Term)
    names(go)=paste(names(go), go.names)
    if(keep.evidence) go.sets=go else go.sets=lapply(go, function(x) unique(x))
    go.sets.len=sapply(go.sets, length)

    idx=grep("cellular_component|biological_process|molecular_function",names(go.sets.len))
    gts=substr(names(go.sets.len),1,10)
    gts.mains=gts[idx]
    gotype=c('BP','CC','MF')
    go.subs=list()
    for(i in 1:3){
      #offsenv=eval(as.name(paste('GO.db::GO',gotype[i],'OFFSPRING',sep='')))
      offsenv=eval(parse(text=paste('GO.db::GO',gotype[i],'OFFSPRING',sep='')))
      branches=AnnotationDbi::get(gts.mains[i], env=offsenv)
      go.subs[[i]]=which(gts %in% branches)
    }
    names(go.subs)=gotype
    go.mains=idx

    res=list(go.sets=go.sets, go.subs=go.subs, go.mains=go.mains)
    return(res)
  }
