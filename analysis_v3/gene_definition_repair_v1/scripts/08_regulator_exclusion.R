# Reuse only old numerical helpers/read-only environment; redirect every output.
source(file.path(commandArgs(trailingOnly=TRUE)[1],"scripts/common.R"))
OUT<-normalizePath(commandArgs(trailingOnly=TRUE)[2],winslash="/",mustWork=TRUE)
put<-function(d,f){p<-file.path(OUT,f);dir.create(dirname(p),recursive=TRUE,showWarnings=FALSE);con<-if(grepl("gz$",f))gzfile(p,"wt")else file(p,"wt");write.table(d,con,sep="\t",row.names=FALSE,quote=FALSE,na="NA");close(con)}
suppressPackageStartupMessages(library(decoupleR))
prior<-readtab("results/regulators/REGULON_SNAPSHOT.tsv")
v<-read.delim(file.path(OUT,"inputs/REGULATOR_EXCLUSION_VOCABULARY.tsv"),check.names=FALSE)
excluded<-toupper(v$gene[tolower(as.character(v$exclude))=="true"])
allcov<-list();allact<-list();allres<-list()
for(ds in c("GSE24290","GSE148059")){
 a<-getdata(ds);rownames(a$y)<-toupper(rownames(a$y));id<-if(ds=="GSE24290")"C01"else"C02";dd<-getdesign(a$m,id)
 p<-prior;p$target<-toupper(p$target);p$tf<-toupper(p$tf);original<-tapply(p$target,p$tf,function(x)length(unique(x)))
 edge<-paste(p$tf,p$target,sep="::");valid<-tapply(p$mor,edge,function(x)length(unique(x))==1);p<-p[valid[edge],];p<-p[!duplicated(paste(p$tf,p$target)),]
 p<-p[!p$target%in%excluded & p$target%in%rownames(a$y),];measured<-table(p$tf)
 cv<-data.frame(dataset=ds,tf=names(original),original_targets=as.integer(original),measured_targets=as.integer(measured[names(original)]));cv$measured_targets[is.na(cv$measured_targets)]<-0;cv$coverage<-cv$measured_targets/cv$original_targets;cv$eligible<-cv$measured_targets>=10 & cv$coverage>=.2
 # The prior/TF universe is fixed. No TF ranking or significant-only selection.
 p<-p[p$tf%in%cv$tf[cv$eligible],];act<-as.data.frame(run_ulm(a$y,data.frame(source=p$tf,target=p$target,mor=p$mor),.source=source,.target=target,.mor=mor,minsize=10));names(act)[names(act)=="p_value"]<-"single_sample_target_fit_P"
 mat<-xtabs(score~source+condition,act);mat<-mat[,a$m$sample_id,drop=FALSE];z<-fitout(mat,dd$X,dd$c,id,trend=FALSE)$table;names(z)[1]<-"tf";z$dataset<-ds;z$prior<-"ABC_EXCLUDE_REPAIRED_PROGRAM"
 cv$prior<-z$prior[1];act$dataset<-ds;allcov[[ds]]<-cv;allact[[ds]]<-act;allres[[ds]]<-z;cat(ds,nrow(z),"completed\n")
}
put(do.call(rbind,allcov),"results/REGULATOR_REPAIRED_COVERAGE.tsv");put(do.call(rbind,allact),"results/REGULATOR_REPAIRED_ACTIVITY.tsv.gz");new<-do.call(rbind,allres);old<-readtab("results/regulators/TF_PROGRAM_EXCLUSION_SENSITIVITY.tsv");put(merge(new,old,by=c("dataset","tf"),all=TRUE,suffixes=c("_new","_old")),"results/REGULATOR_EXCLUSION_OLD_NEW.tsv")
capture.output(sessionInfo(),file=file.path(OUT,"logs/REGULATOR_SESSION_INFO.txt"))
