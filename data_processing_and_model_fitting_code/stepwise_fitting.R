
reg.stepwise.fit = function(K.grid, num.cores = 1, max.profs = 5, k.cv = 5, metric = "L2", jacc.cutoff = if (metric == "L2") NA) {
  set.seed(0)
  
  # load non-cluster data
  load("profiles500_train80_seed0.Rdata")
  
  # load clusters
  load("clusters/public500_naiveAA_kmeans.Rdata")
  load("clusters/public500_kmeans.Rdata")
  
  # list the profile names
  prof.names = c(paste("clusts", K.grid, sep="-"), paste("naiveAAclusts", K.grid, sep="-"), "naiveAA", "vsubgrp", "vgene", "neut")
  
  # setup CV folds
  cv.folds = split(1:nrow(dff2.full.train), sample(1:nrow(dff2.full.train)) %% k.cv)
  
  # make alpha assignments
  alphas.cols.assign = as.numeric(gsub("^p_([0-9]+)_a_([0-9]+)$", "\\1", colnames(dff2.full.train)))
  alpha.groups = length(unique(alphas.cols.assign))
  
  # setup cluster
  cl = makeCluster(num.cores)
  clusterExport(cl, c("dff2.subsamp.train", "dff2.full.train", "dff2.naiveAA.train", "dff2.vsubgrp.train",
                      "dff2.vgene.train", "dff2.neut.train", "dff2.nogap.train", "alpha.groups",
                      "alphas.cols.assign", "cv.folds", "dff.kmeans", "dff.naiveAA.kmeans"), envir=environment())
  clusterEvalQ(cl,{
    Rcpp::sourceCpp("stepwise_fitting.cpp")
    NULL})
  
  # start stepwise fitting
  curr.best.names = c()
  curr.best.results = c()
  
  for (i in 1:max.profs) {
    print(i)
    
    # update the potential profile list
    curr.names = lapply(prof.names, function(name) c(curr.best.names, name))
    
    # make this round's parameter grid
    param.grid = expand.grid("curr.names" = curr.names, "fold" = names(cv.folds))
    clusterExport(cl, c("i", "param.grid"), envir=environment())
    
    results = clusterApplyLB(cl, 1:nrow(param.grid), function(k) {
      # extract fold indices
      fold.inds = cv.folds[[param.grid[k,"fold"]]]
      
      # if clusters are needed, extract the right clusters
      name.inds = grep("[0-9]+", param.grid[k,"curr.names"][[1]])
      for (ind in name.inds) {
        clust.str = param.grid[k,"curr.names"][[1]][ind]
        clust.name = gsub("^([A-Za-z]*)clusts-([0-9]+)$", "\\1", clust.str)
        clust.K = gsub("^([A-Za-z]*)clusts-([0-9]+)$", "\\2", clust.str)
        
        centers = get(paste("dff.", clust.name, if (clust.name != "") ".", "kmeans", sep=""), envir=environment())[[clust.K]]$center
        assign(paste("dff2.", clust.str, ".train", sep=""),
               compute_kmeans_clust_profs(centers, dff2.subsamp.train), envir=environment())
        assign(paste("dff2.", clust.str, ".train.cv", sep=""),
               get(paste("dff2.", clust.str, ".train", sep=""), envir=environment())[-fold.inds,])
        assign(paste("dff2.", clust.str, ".test.cv", sep=""),
               get(paste("dff2.", clust.str, ".train", sep=""), envir=environment())[fold.inds,])
      }
      
      # split other profiles into CV-train and CV-test sets
      dff2.full.train.cv = dff2.full.train[-fold.inds,]
      dff2.subsamp.train.cv = dff2.subsamp.train[-fold.inds,]
      dff2.neut.train.cv = dff2.neut.train[-fold.inds,]
      dff2.naiveAA.train.cv = dff2.naiveAA.train[-fold.inds,]
      dff2.vsubgrp.train.cv = dff2.vsubgrp.train[-fold.inds,]
      dff2.vgene.train.cv = dff2.vgene.train[-fold.inds,]
      dff2.nogap.train.cv = dff2.nogap.train[-fold.inds,]
      
      dff2.full.test.cv = dff2.full.train[fold.inds,]
      dff2.subsamp.test.cv = dff2.subsamp.train[fold.inds,]
      dff2.neut.test.cv = dff2.neut.train[fold.inds,]
      dff2.naiveAA.test.cv = dff2.naiveAA.train[fold.inds,]
      dff2.vsubgrp.test.cv = dff2.vsubgrp.train[fold.inds,]
      dff2.vgene.test.cv = dff2.vgene.train[fold.inds,]
      dff2.nogap.test.cv = dff2.nogap.train[fold.inds,]
      
      # make the profile lists
      profs.train.cv = mget(paste("dff2.", param.grid[k,"curr.names"][[1]], ".train.cv", sep=""), envir=environment())
      profs.test.cv = mget(paste("dff2.", param.grid[k,"curr.names"][[1]], ".test.cv", sep=""), envir=environment())
      
      # optimize "alphas"
      optim.obj.alphas = optim(par=rep(0, i*alpha.groups), f, method="L-BFGS-B", lower=0, upper=1,
                               alpha_cols_assign=alphas.cols.assign, dff2_subsamp=dff2.subsamp.train.cv, dff2_full=dff2.full.train.cv,
                               dff2_profs=profs.train.cv, alpha_groups=alpha.groups, control=list(maxit=150))
      
      if (metric == "L2") {
        fold.metric = f_nogap(optim.obj.alphas$par, alphas.cols.assign, dff2.subsamp.test.cv, dff2.full.test.cv,
                              profs.test.cv, dff2.nogap.test.cv, alpha.groups)
      } else {
        fold.metric = f_nogap_jacc(optim.obj.alphas$par, alphas.cols.assign, dff2.subsamp.test.cv, dff2.full.test.cv,
                                   profs.test.cv, dff2.nogap.test.cv, alpha.groups, jacc.cutoff, 99999)
      }
      
      return(fold.metric)
    })
    
    # extract the model results
    cv.metric = unlist(results)
    metric.mat = cbind(param.grid, cv.metric)
    curr.names.cat = sapply(metric.mat[,"curr.names"], function(names) paste(names, collapse="_"))
    metric.results = c(by(metric.mat, curr.names.cat, function(mat) mean(mat[,3])))
    
    # cache the best model results
    best.ind = if (metric == "L2") which.min(metric.results) else which.max(metric.results)
    curr.best.names = strsplit(names(metric.results)[best.ind], "_")[[1]]
    curr.best.results = c(curr.best.results, metric.results[best.ind])
    
    # narrow down the profile search list
    prof.names = prof.names[prof.names != curr.best.names[length(curr.best.names)]]
  }
  
  # close the cluster
  stopCluster(cl)
  
  return(curr.best.results)
}






jacc.stepwise.fit = function(cutoff, B.grid, K.grid, num.cores = 1, max.profs = 5, k.cv = 5) {
  set.seed(0)
  
  # load non-cluster data
  load("profiles500_train80_seed0.Rdata")
  
  # load clusters
  load("clusters/public500_naiveAA_kmeans.Rdata")
  load("clusters/public500_kmeans.Rdata")
  
  # list the profile names
  prof.names = c(paste("clusts", K.grid, sep="-"), paste("naiveAAclusts", K.grid, sep="-"), "naiveAA", "vsubgrp", "vgene", "neut")
  
  # setup CV folds
  cv.folds = split(1:nrow(dff2.full.train), sample(1:nrow(dff2.full.train)) %% k.cv)
  
  # make alpha assignments
  alphas.cols.assign = as.numeric(gsub("^p_([0-9]+)_a_([0-9]+)$", "\\1", colnames(dff2.full.train)))
  alpha.groups = length(unique(alphas.cols.assign))
  
  # setup cluster
  cl = makeCluster(num.cores)
  clusterExport(cl, c("dff2.subsamp.train", "dff2.full.train", "dff2.naiveAA.train", "dff2.vsubgrp.train",
                      "dff2.vgene.train", "dff2.neut.train", "dff2.nogap.train", "alpha.groups",
                      "alphas.cols.assign", "cv.folds", "dff.kmeans", "dff.naiveAA.kmeans", "cutoff"), envir=environment())
  clusterEvalQ(cl,{
    Rcpp::sourceCpp("stepwise_fitting.cpp")
    NULL})
  
  # start stepwise fitting
  curr.best.names = c()
  curr.best.results = c()
  
  for (i in 1:max.profs) {
    print(i)
    
    # update the potential profile list
    curr.names = lapply(prof.names, function(name) c(curr.best.names, name))
    
    # make this round's parameter grid
    param.grid = expand.grid("curr.names" = curr.names, "B" = B.grid, "fold" = names(cv.folds))
    clusterExport(cl, c("i", "param.grid"), envir=environment())
    
    results = clusterApplyLB(cl, 1:nrow(param.grid), function(k) {
      # extract fold indices
      fold.inds = cv.folds[[param.grid[k,"fold"]]]
      
      # if clusters are needed, extract the right clusters
      name.inds = grep("[0-9]+", param.grid[k,"curr.names"][[1]])
      for (ind in name.inds) {
        clust.str = param.grid[k,"curr.names"][[1]][ind]
        clust.name = gsub("^([A-Za-z]*)clusts-([0-9]+)$", "\\1", clust.str)
        clust.K = gsub("^([A-Za-z]*)clusts-([0-9]+)$", "\\2", clust.str)
        
        centers = get(paste("dff.", clust.name, if (clust.name != "") ".", "kmeans", sep=""), envir=environment())[[clust.K]]$center
        assign(paste("dff2.", clust.str, ".train", sep=""),
               compute_kmeans_clust_profs(centers, dff2.subsamp.train), envir=environment())
        assign(paste("dff2.", clust.str, ".train.cv", sep=""),
               get(paste("dff2.", clust.str, ".train", sep=""), envir=environment())[-fold.inds,])
        assign(paste("dff2.", clust.str, ".test.cv", sep=""),
               get(paste("dff2.", clust.str, ".train", sep=""), envir=environment())[fold.inds,])
      }
      
      # split other profiles into CV-train and CV-test sets
      dff2.full.train.cv = dff2.full.train[-fold.inds,]
      dff2.subsamp.train.cv = dff2.subsamp.train[-fold.inds,]
      dff2.neut.train.cv = dff2.neut.train[-fold.inds,]
      dff2.naiveAA.train.cv = dff2.naiveAA.train[-fold.inds,]
      dff2.vsubgrp.train.cv = dff2.vsubgrp.train[-fold.inds,]
      dff2.vgene.train.cv = dff2.vgene.train[-fold.inds,]
      dff2.nogap.train.cv = dff2.nogap.train[-fold.inds,]
      
      dff2.full.test.cv = dff2.full.train[fold.inds,]
      dff2.subsamp.test.cv = dff2.subsamp.train[fold.inds,]
      dff2.neut.test.cv = dff2.neut.train[fold.inds,]
      dff2.naiveAA.test.cv = dff2.naiveAA.train[fold.inds,]
      dff2.vsubgrp.test.cv = dff2.vsubgrp.train[fold.inds,]
      dff2.vgene.test.cv = dff2.vgene.train[fold.inds,]
      dff2.nogap.test.cv = dff2.nogap.train[fold.inds,]
      
      # make the profile lists
      profs.train.cv = mget(paste("dff2.", param.grid[k,"curr.names"][[1]], ".train.cv", sep=""), envir=environment())
      profs.test.cv = mget(paste("dff2.", param.grid[k,"curr.names"][[1]], ".test.cv", sep=""), envir=environment())
      
      # optimize "alphas"
      optim.obj.alphas = optim(par=rep(0, i*alpha.groups), f_jacc, method="L-BFGS-B", lower=0, upper=1,
                               alpha_cols_assign=alphas.cols.assign, dff2_subsamp=dff2.subsamp.train.cv, dff2_full=dff2.full.train.cv,
                               dff2_profs=profs.train.cv, alpha_groups=alpha.groups,
                               cutoff=cutoff, B=param.grid[k,"B"], control=list(maxit=150))
      
      fold.jacc = f_nogap_jacc(optim.obj.alphas$par, alphas.cols.assign, dff2.subsamp.test.cv, dff2.full.test.cv,
                               profs.test.cv, dff2.nogap.test.cv, alpha.groups, cutoff, 99999)
      
      return(fold.jacc)
    })
    
    # extract the model results
    cv.jacc = unlist(results)
    jacc.mat = cbind(param.grid, cv.jacc)
    jacc.results.names = paste(sapply(jacc.mat[,"curr.names"], function(names) paste(names, collapse="_")), jacc.mat[,"B"], sep="_")
    jacc.results = by(jacc.mat, jacc.results.names, function(mat) mean(mat[,4]))
    jacc.results = c(jacc.results)
    
    # cache the best model results
    best.ind = which.min(jacc.results)
    curr.best.names = strsplit(names(jacc.results)[best.ind], "_")[[1]][-(i+1)]
    curr.best.results = c(curr.best.results, -jacc.results[best.ind])
    
    # narrow down the profile search list
    prof.names = prof.names[prof.names != curr.best.names[length(curr.best.names)]]
  }
  
  # close the cluster
  stopCluster(cl)
  
  return(curr.best.results)
}
