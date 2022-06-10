stddiff.numeric<-function(data,gcol,vcol){
  data[,gcol]<-as.factor(data[,gcol])
  rst<-matrix(rep(0,9*length(vcol)),ncol=9)
  dimnames(rst)<-list(names(data)[vcol],
    c("mean.c","sd.c","mean.t","sd.t","missing.c","missing.t","stddiff","stddiff.l","stddiff.u"))
  for(i in 1:length(vcol)){
    data[,vcol[i]]<-as.numeric(data[,vcol[i]])
    na.c<-length(which(is.na(data[,vcol[i]][which(data[,gcol]==levels(data[,gcol])[1])])))
    na.t<-length(which(is.na(data[,vcol[i]][which(data[,gcol]==levels(data[,gcol])[2])])))
    temp<-na.omit(data[,c(gcol,vcol[i])])
    m<-aggregate(temp[,2],by=list(temp[,1]),FUN=mean)
    s<-aggregate(temp[,2],by=list(temp[,1]),FUN=sd)
    stddiff<-abs(m[2,2]-m[1,2])/sqrt((s[2,2]^2+s[1,2]^2)/2)
    n<-table(temp[,1])
    se<-sqrt(nrow(temp)/(n[1]*n[2])+stddiff^2/(2*nrow(temp)))
    stddiff.l<-stddiff-1.96*se
    stddiff.u<-stddiff+1.96*se
    rst[i,]<-c(m[1,2],s[1,2],m[2,2],s[2,2],na.c,na.t,stddiff,stddiff.l,stddiff.u)
  }
  rst<-round(rst,3)
  return(rst)
}

stddiff.binary<-function(data,gcol,vcol){
  for(i in 1:length(c(gcol,vcol))){
    data[,c(gcol,vcol)[i]]<-as.factor(data[,c(gcol,vcol)[i]])
  }
  rst<-matrix(rep(0,7*length(vcol)),ncol=7)
  dimnames(rst)<-list(names(data)[vcol],
    c("p.c","p.t","missing.c","missing.t","stddiff","stddiff.l","stddiff.u"))
  for(i in 1:length(vcol)){
    na.c<-length(which(is.na(data[,vcol[i]][which(data[,gcol]==levels(data[,gcol])[1])])))
    na.t<-length(which(is.na(data[,vcol[i]][which(data[,gcol]==levels(data[,gcol])[2])])))
    temp<-na.omit(data[,c(gcol,vcol[i])])
    temp[,2]<-as.numeric(temp[,2])
    p1<-mean(temp[,2][which(temp[,1]==levels(temp[,1])[2])])-1
    p2<-mean(temp[,2][which(temp[,1]==levels(temp[,1])[1])])-1
    stddiff<-abs(p1-p2)/sqrt((p1*(1-p1)+p2*(1-p2))/2)
    n<-table(temp[,1])
    se<-sqrt(nrow(temp)/(n[1]*n[2])+stddiff^2/(2*nrow(temp)))
    stddiff.l<-stddiff-1.96*se
    stddiff.u<-stddiff+1.96*se
    rst[i,]<-c(p2,p1,na.c,na.t,stddiff,stddiff.l,stddiff.u)
  }
  rst<-round(rst,3)
  return(rst)
}

stddiff.category<-function(data,gcol,vcol){
  for(i in 1:length(c(gcol,vcol))){data[,c(gcol,vcol)[i]]<-as.factor(data[,c(gcol,vcol)[i]])}
  nr<-NA
  for(i in 1:length(vcol)){nr[i]<-length(levels(data[,vcol[i]]))}
  rst<-matrix(rep(0,7*sum(nr)),ncol=7)
  rname<-NA
  for(i in 1:length(nr)){
    rname<-c(rname,paste(names(data)[vcol[i]],levels(data[,vcol[i]])))
  }
  dimnames(rst)<-list(rname[-1],
    c("p.c","p.t","missing.c","missing.t","stddiff","stddiff.l","stddiff.u"))
  for(i in 1:length(vcol)){
    na.c<-length(which(is.na(data[,vcol[i]][which(data[,gcol]==levels(data[,gcol])[1])])))
    na.t<-length(which(is.na(data[,vcol[i]][which(data[,gcol]==levels(data[,gcol])[2])])))
    temp<-na.omit(data[,c(gcol,vcol[i])])
    tbl<-table(temp[,2],temp[,1])
    prop<-prop.table(tbl,2)
    t<-prop[-1,2]#k=2,3,...,K
    c<-prop[-1,1]#k=2,3,...,K
    k<-nr[i]-1
    l<-k
    s<-matrix(rep(0,k*l),ncol=l)
    for(ii in 1:k){
      for(j in 1:l){
        if(ii==j){s[ii,j]<-0.5*(t[ii]*(1-t[ii])+c[ii]*(1-c[ii]))}
        if(ii!=j){s[ii,j]<--0.5*(t[ii]*t[j]+c[ii]*c[j])}
      }
    }
    e<-rep(1,k)
    e<-diag(e)
    s<-solve(s,e)
    tc1<-t-c
    tc2<-t-c
    stddiff<-sqrt(t(tc1) %*% s %*% tc2)
    n<-table(temp[,1])
    se<-sqrt(nrow(temp)/(n[1]*n[2])+stddiff^2/(2*nrow(temp)))
    stddiff.l<-stddiff-1.96*se
    stddiff.u<-stddiff+1.96*se
    if(i==1){
      rst[1:nr[i],]<-cbind(prop,c(na.c,rep(NA,k)),c(na.t,rep(NA,k)),c(stddiff,rep(NA,k)),c(stddiff.l,rep(NA,k)),c(stddiff.u,rep(NA,k)))
    }
    if(i>1){
      rst[(sum(nr[1:(i-1)])+1):(sum(nr[1:(i-1)])+nr[i]),]<-cbind(prop,c(na.c,rep(NA,k)),c(na.t,rep(NA,k)),c(stddiff,rep(NA,k)),c(stddiff.l,rep(NA,k)),c(stddiff.u,rep(NA,k)))
    }
  }
  rst<-round(rst,3)
  return(rst)
}

