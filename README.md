# color
### FUNCTION TO COMPUTE CHROMATIC AND ACHROMATIC COLOUR VARIATION IN TETRA VISUAL SYSTEM  


tetra<-function(Ref1, Irr, Tetra, pcones=NULL, weber=NULL, doubleCone=NULL, wLong=0.05, idVars=NULL, wlRef1=NULL, wlIrr=1, irrIrr=2, wlTetra=1, conesTetra=2:5, wlDoubleCone=1, sensDoubleCone=2, wlInclude=TRUE)
 {
 require(plyr)
 require(reshape2)
 IrrArg<-deparse(substitute(Irr))
 irrName<-names(Irr)[names(Irr)== names(Irr[irrIrr])]
 TetraArg<-deparse(substitute(Tetra))
 if(is.null(idVars)){
 idVars<-names(Ref1)[!(gsub("[[:digit:]]", "", names(Ref1))%in% "wl")]
 }
 if(is.null(wlRef1)){
 wlRef1<-names(Ref1)[gsub("[[:digit:]]", "", names(Ref1))%in% "wl"]
 }
idRefl<-Ref1[idVars]
dfRefl<-Ref1[wlRef1]
ReflW1<-as.numeric(gsub("[[:punct:]]", "", gsub("[[:alpha:]]", "", names(dfRefl))))
names(Irr)[names(Irr)==names(Irr[wlIrr])]<-"wl"
names(Irr)[names(Irr)==names(Irr[irrIrr])]<-"irr"
names(Tetra)[names(Tetra)==names(Tetra[wlTetra])]<-"wl"
names(Tetra)[conesTetra]<-c("VS", "S", "M", "L")
commonWL<-Reduce(intersect, list(round(Irr$wl), round(Tetra$wl), round(ReflW1)))

if (!is.null(doubleCone)){
names(doubleCone)[names(doubleCone)==names(doubleCone[wlDoubleCone])]<-"wl"
names(doubleCone)[names(doubleCone)==names(doubleCone[sensDoubleCone])]<-"sens"
commonWL<-intersect(commonWL, round(doubleCone$wl))
doubleCone<-doubleCone[round(doubleCone$wl)%in% commonWL,
]
}
Irr<-Irr[round(Irr$wl)%in%commonWL, ]
Tetra<-Tetra[round(Tetra$wl)%in% commonWL, ]
dfRefl<-dfRefl[, round(ReflW1)%in% commonWL]
reso<-unique(diff(Irr$wl))
matrixRef1<-as.matrix(dfRefl)
if(!is.null(doubleCone)){
irrDC<-cbind(doubleCone$sens)*Irr$irr
dcCatch<-c((matrixRef1%*%irrDC)*reso)
}
irrcones<-as.matrix(Tetra[, c("VS","S","M","L")]*Irr$irr)
Tcatch<-(matrixRef1%*%irrcones)*reso
if(!is.null(pcones)|!is.null(weber)){
if(!is.null(pcones)){
W<-(wLong/sqrt(pcones))^2
}
if(!is.null(weber)){
W<-weber
}
H<-c(1/sqrt(W[1]+W[2]), sqrt((W[1]+W[2])/(W[1]*
W[2]+W[1]*W[3]+W[2]*W[3])), sqrt((W[1]* 
W[2]+W[1]*W[3]+W[2]*W[3])/(W[1]* W[2]* 
W[3]+W[1]*W[2]*W[4]+W[1]*W[3]*W[4]+  
W[2]*W[3]*W[4])))

h<-c(W[2]/(W[1]+W[2]), W[1]/(W[1]+W[2]), W[2]* 
W[3]/(W[1]* W[2]+W[2]*W[3]+W[1]*W[3]),W[1]* 
W[3]/(W[1]*W[2]+W[2]*W[3]+W[1]*W[3]),W[1]* 
W[2]/(W[1]*W[2]+W[2]*W[3]+W[1]*W[3]))
TcatchTot<-apply(Tcatch, 1, sum)

lcr<-log(Tcatch/matrix(rep(TcatchTot, 4), ncol=4))

x<-H[1]*(lcr[,1]-lcr[,2])
y<-H[2]*(lcr[,3]-(h[1]*lcr[,1]+h[2]*lcr[,
2]))
z<-H[3]*(lcr[,4]-(h[3]*lcr[,1]+h[4]*lcr[,2]+h[5]*lcr[,3]))
}
out<-data.frame(idRefl, wlResol=reso, Irr=paste(IrrArg, 
irrName, sep="$"), Tetra=TetraArg)
if(!is.null(pcones)){
out<-data.frame(out, conesProp=paste(round(pcones, 2), collapse=","))
}
if(!is.null(pcones) | !is.null(weber)){
out<-data.frame(out, WeberFrac=paste(round(sqrt(W), 4), collapse=","))
}
out<-data.frame(out, Tcatch)

if(!is.null(doubleCone)){
out<-data.frame(out, DC=dcCatch)
}

if(!is.null(pcones)|!is.null(weber)){
out<-data.frame(out, x=x, y=y, z=z )
}

if(wlInclude){
out<-cbind(out, dfRefl)
}
out
}

