"maboost.default" <-
  function(x,y,test.x,test.y=NULL,breg=c("entrop","l2"),C50tree=FALSE,iter=100, nu, bag.frac=0.5,random.feature=TRUE,random.cost=TRUE,
           maxmargin=FALSE,smooth=FALSE,smoothfactor=1,sparse=T,delta=10^(-10),verbose=FALSE,
           ...,na.action=na.rpart){
    cl<-match.call(expand.dots=TRUE)
    cl[[1]]<-as.name("maboost")
    multiclass=FALSE;
    
    
      type="discrete"
   
    if(missing(breg))
      breg="quad"
    breg=tolower(breg)
    if(breg=="expo"|breg=="exponential"){
      breg="entrop"
    }
    if(breg=="l"|breg=="quad"){
      breg="l2"
    }
    if(missing(nu)){
      if(breg=="l2")  nu=2;
      if(breg=="entrop")  nu=0.1;
    }
    if(missing(y) | missing(x)){
      stop("This procedure requires a response and a set of variables")
    }
    if( smooth ){
      if(breg=='entrop') {
        warning('for the moment smooth and sparse boosters are only available for l2. Bregman divergence was changed to l2');
        breg="l2";
      }
      
      if(missing(smoothfactor)|| smoothfactor>N || smoothfactor<1){
        stop("With smooth option, smoothfactor should be set to an integer between 1 and N so that w_i<1/smoothfactor")
      }
    
    }
    lev=levels(as.factor(y))
    K=length(lev)
    if(length(lev)>2){
      print('Multiclass boosting is selected'); flush.console();
      multiclass=TRUE;
 
    }
    if (C50tree)  na.action=na.pass
    
    
    x=as.data.frame(x)
    M=dim(x)[1];
    test.dat=FALSE
    if(!is.null(test.y)){
      tlev=NULL
      if(!missing(test.x)){
        tlev=levels(as.factor(as.vector(as.matrix(test.y))))
        test.x=as.data.frame(test.x)
        test.dat=TRUE
      }
      if(length(tlev)<1)
        warning("test.x must be the testing data and the response must have at least 2 levels")
    }else{
      test.x=NULL
    }

    extraArgs=list(...)
    if (length(extraArgs)) {
      if(!C50tree)    arg <- c('maxdepth','cp','minsplit','surrogatestyle','xval','usesurrogate','maxsurrogate','minbuckt','maxcompete')
      else  arg <-c('subset','bands', 'winnow','noGlobalPruning', 
                    'CF', 'minCases', 'fuzzyThreshold', 'sample', 
                    'seed', 'earlyStopping','label')
       
      indx <- match(names(extraArgs), arg, nomatch = 0)
      if (any(indx == 0)) 
        stop(paste("Error:  Argument", names(extraArgs)[indx == 0], 
                   "not matched"))
    }
    if(multiclass){
      ny=y;
      levels(ny) = 1 : ( length(levels(y)) )
      if(test.dat){
        test.y=as.factor(as.vector(as.matrix(test.y)))
        lev_inter=intersect(levels(test.y),levels(y))
        if(length(lev_inter)==K && length(lev_inter)==length(levels(test.y)) ){
        levels(test.y)<-levels(ny);
        }else{
          warning('testset has less or more classes than the training set')
         ind_tr=which(levels(y) %in% lev_inter) 
         ind_te=which(levels(test.y) %in% lev_inter)
         levels(test.y)[ind_te]=levels(ny)[ind_tr];
        }
      }
    }else{
      ny<-c(-1,1)[as.numeric(as.factor(y))]
      if(test.dat)
        test.y<-c(-1,1)[as.numeric(as.factor(as.vector(test.y)))]
    }
  
    
    ### Set Up Predictions for each boosting type
    if(!multiclass){
    method="class"
    if(type=="discrete"){
      predict.type<-function(f,dat){ 
        p <- predict(f,newdata=dat,type="class") 
        a=as.numeric(levels(p))[p] 
        a 
      }
    }
    if(type=="real"){
      predict.type<-function(f,dat){
        a=predict(f,newdata=dat,type="prob")[,2]
        f=(a-.5)*2
        f
      }
    }
    #d=eta*f where f is the last fit
    
    dcal=function(f,yf){
     d=yf*f;
      
    }
    if(breg=="entrop"){
      if(maxmargin){
      coefs=function(d,w,t){
       
        eta=max(sum(w*d),0);
        eta/t^.5;
      }
      }else{
        coefs=function(d,w,t){
          eta=max(sum(w*d),0);
          eta;
        }
      }
      wfun<-function(w,d,eta){
     
        w=w*exp(-eta*d)
        w=w/sum(w)
        w
      }
    }
    if(breg=="l2"){
      if(maxmargin){
      coefs=function(d,w,t){
        eta=max(sum(w*d)/M,0);
        eta/t^.5
      }
    }else{
      coefs=function(d,w,t){
        eta=max(sum(w*d)/M,0);
        eta
      }
    }
     wfun<-function(w,d,eta){
        w=projsplx(w-eta*d)
        w
      }
      if(smooth){
        wfun<-function(w,d,eta){
        w=projsplx_k(w-eta*d,smoothfactor)
        w
       }
      }
      if(sparse){
        wfun<-function(w,d,eta){
          
          l1_co=eta^2*M/2;
          w=pmin(pmax(w - eta*d - l1_co/2,0),1);
          print(c(length(which(w==0)),sum(w)))
          w/sum(w)
        }
      }
      if(sparse && smooth){
        wfun<-function(w,d,eta){
          l1_co=eta^2*M;
          y=pmin(pmax(w - eta*d - l1_co/2,0),1);
          w=projsplx_k(y,smoothfactor)
          w
        }
      }
      
    }
    
    
  }else{  #when it is multiclass
    method="class"
    #for the momenet multiclass has only been implemented for discrete outputs
    if(type=="discrete" || type=="real" ){
      predict.type<-function(f,dat){ 
        p <- predict(f,newdata=dat,type="class") 
        a =  as.numeric(levels(p))[p];
        a 
      }
    }
    dcal=function(f,yf){
      yf=as.numeric(levels(yf))[yf];
      ind_corr= which(f==yf);
      ind_notcorr=setdiff(1:M,ind_corr);
      col_zeros =f[ind_notcorr]-yf[ind_notcorr];
 
      D=matrix(0,M,K-1);
      D[ind_corr,]=-1;
      D [ cbind( ind_notcorr[col_zeros<0], f[ind_notcorr[col_zeros<0]]) ] = 1;
      D [ cbind( ind_notcorr[col_zeros>0], f[ind_notcorr[col_zeros>0]]-1) ] = 1;
      D
     
      
    }
    if(breg=="entrop"){
      if(maxmargin){
      coefs=function(D,W,t){
        eta=max(-sum(W*D),0)
        eta/t^.5
      }
      }else{
        coefs=function(D,W,t){
          eta=max(-sum(W*D),0)
          eta
        }
      }
      wfun<-function(W,D,eta){
        W=W*exp(eta*D)
        W=W/sum(W)
        W
        
      }
    }
    if(breg=="l2"){
      if(maxmargin){
      coefs=function(D,W,t){
        eta=max(-sum(W*D)/((K-1)*M),0)
        eta/t^.5;
      }
      }else{
          coefs=function(D,W,t){
          eta=max(-sum(W*D)/((K-1)*M),0)
          eta
      }
      }
      wfun<-function(W,D,eta){
        W=projsplx(W + eta * D)
        W 
      }
          if(smooth){stop("smooth boosting for multiclass setting has not yet been implemented")}
          if(sparse){
            
            wfun<-function(W,D,eta){
              l1_co=eta^2(K-1)*M;
              
              W=pmin(pmax(W + eta * D - l1_co/2,0),1);
              #W=projsplx(Y)
              W 
            }
          }
 
          
    }
  }#end else (multiclass else if)
    ### Set up coeficients

    lossObj=list(dcal=dcal,coefs=coefs,wfun=wfun,predict.type=predict.type,method=method,type=type,breg=breg,K=K,Ctree=C50tree)
    if(multiclass){
      result =maboost.machine.multi(x,ny,test.x,test.y,iter,maxmargin,smooth,smoothfactor,sparse,nu,bag.frac,random.feature,random.cost,lossObj,oldObj=NULL,na.action=na.action,...)
      
      g = factor ( apply(result$F[[1]],1,which.max),levels=1:K )
      levels(g)=lev;
      
    }else{
      
      result =maboost.machine.bin(x,ny,test.x,test.y,iter,maxmargin,smooth,smoothfactor,sparse,nu,bag.frac,random.feature,random.cost,lossObj,oldObj=NULL,na.action=na.action,...)
      g  =  factor(sign(result$F[[1]]),levels=c(-1,1));
      levels(g)=lev;
      
      
     
    }
   

    tab=table(as.factor(y),g)
    nm<-1:(dim(x)[2])
    if(is.data.frame(x)){
      nm=names(x)
    } 
    obj=structure(list(model=result,fit=g,call=cl,confusion=tab,iter=iter,
                       actual=as.factor(y),nu=nu,dim=dim(x),names=nm,bag.frac=bag.frac
                       ,na.action=na.action,
                       Ctree=C50tree,maxmargin=maxmargin
                       ,smooth=smooth,
                       smoothfactor=smoothfactor,
                       random.feature=random.feature,
                       random.cost=random.cost,sparse=sparse)
                  ,class="maboost")
    obj
  }
