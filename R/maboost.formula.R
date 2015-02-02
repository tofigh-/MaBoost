"maboost.formula" <-
  function(formula, data,...,subset,na.action=na.rpart){
    ## m = match.call(expand.dots = FALSE)
    ##m[[1]] = as.name("model.frame")
    ##m$...=NULL
    ##m =eval(m,parent.frame())
    if (!(require('rpart', character.only=T, quietly=T))) {
      install.packages('rpart')
      library('rpart', character.only=T)
    }
    if (!(require('C50', character.only=T, quietly=T))) {
      install.packages('C50')
      library('C50', character.only=T)
    }
    
    m <- match.call(expand.dots = FALSE)
    m$model <- m$method <- m$control <- NULL
    m$x <- m$y <- m$parms <- m$... <- NULL
    m$cost <- NULL
    m$na.action <- na.action
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())

    Terms = attr(m, "terms")
  
    y = as.factor(model.extract(m,"response"))
    if(class(y)!="factor") {
      warning('Class variable should be pass as factor. It is now automatically converted to factor');
      y=as.factor(y);
    }
    preds<-attr(attributes(m)$terms,"term.labels")
    x<-as.data.frame(m[,!is.na(match(names(m),preds))])

    res = maboost.default(x,y,...,na.action=na.action)
    res$terms = Terms
    cl = match.call()
    cl[[1]] = as.name("maboost")
    res$call = cl
    res
  }

