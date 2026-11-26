#!/usr/bin/env Rscript
# Additional secondary component analysis. Revision author: mjwang.
suppressPackageStartupMessages({library(mice); library(digest)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop('Pass project root')
root <- normalizePath(args[1]); version <- '20261023_001156'
paper <- file.path(root,'paper1_cmrf_pathology')
src <- file.path(paper,'results/versions/20260809_152651')
out <- file.path(paper,'results/versions',version)
if(file.exists(file.path(out,'source_data/steatosis_component_models.csv'))) stop('Results already exist')
checks <- list(); ck <- function(n,p){checks[[n]] <<- isTRUE(p); if(!isTRUE(p))stop(n)}
rd <- function(n)read.csv(file.path(src,n),stringsAsFactors=FALSE)
cc <- rd('canonical_adult_complete_component.csv'); eligible <- rd('eligible_adult_for_mi.csv')
imp <- readRDS(file.path(src,'multiple_imputation_object.rds'))
ck('CC N935 and MI N1085 m50', nrow(cc)==935 && nrow(eligible)==1085 && nrow(imp$data)==1085 && imp$m==50)
ck('Unique IDs',!anyDuplicated(cc$id) && !anyDuplicated(eligible$id))
ck('CC is subset eligible',all(cc$id %in% eligible$id))
ck('Three centers',identical(sort(unique(cc$center_code)),c('HZ','JL','RJ')))
ck('ZS mapped to RJ',all(cc$center_code[grepl('^ZS',cc$id)]=='RJ') && all(eligible$center_code[grepl('^ZS',eligible$id)]=='RJ'))
ck('No missing pathology',!anyNA(cc$steatosis) && !anyNA(imp$data$steatosis))
ck('No imputation of steatosis',imp$method['steatosis']=='' && all(imp$predictorMatrix['steatosis',]==0))
ck('Steatosis used as MI predictor',all(imp$predictorMatrix[c('wc','glu','hba1c','tg','hdl'),'steatosis']==1))
for(v in names(imp$data)) {
 a <- as.character(imp$data[[v]]); b <- as.character(eligible[[v]])
 if(is.numeric(imp$data[[v]])) ck(paste('MI input order',v),identical(is.na(imp$data[[v]]),is.na(eligible[[v]])) && max(abs(imp$data[[v]]-eligible[[v]]),na.rm=TRUE)<1e-10)
 else {
  # Reproduce the frozen MI factor coding (one nonbinary HCV value became NA).
  if(is.factor(imp$data[[v]])) b <- as.character(factor(eligible[[v]],levels=levels(imp$data[[v]])))
  ck(paste('MI input order',v),identical(a,b))
 }
}
# Extract only reviewed function definitions; do not execute the old pipeline.
expr <- parse(file.path(paper,'code/R/versions/20260809_152651/01_reanalysis.R'))
for(fn in c('derive_cmrf','pool_scalar_rubin')) {
 ix <- which(vapply(expr,function(e)is.call(e) && identical(e[[1]],as.name('<-')) && identical(e[[2]],as.name(fn)),logical(1)))
 ck(paste('Unique source function',fn),length(ix)==1L); eval(expr[[ix]],envir=.GlobalEnv)
}
components <- paste0('cmrf',1:5)
cc$sex_label <- factor(cc$sex,levels=c(0,1),labels=c('Male','Female'))
cc$center_code <- factor(cc$center_code,levels=c('HZ','JL','RJ'))
rederived <- derive_cmrf(cc)
ck('Canonical CMRF agrees with frozen adjudication',all(as.matrix(rederived[components])==as.matrix(cc[components])))
ck('CC adiposity negative 13',sum(cc$cmrf1==0)==13)
ck('CC steatosis distribution',identical(as.integer(table(factor(cc$steatosis,levels=0:3))),c(19L,271L,344L,301L)))
rhs <- c(components,'age','sex_label','center_code')
diags <- list(); estimates <- list(); fit_store <- list(); bfgs <- list()
fit_model <- function(d,outcome,set,imputation=0L) {
 warnings <- character()
 ctl <- if(outcome=='f_ge2')glm.control() else glm.control(epsilon=1e-10,maxit=100)
 fit <- withCallingHandlers(glm(reformulate(rhs,outcome),data=d,family=binomial(),control=ctl),warning=function(w){warnings <<- c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
 ck(paste('Convergence rank',set,outcome,imputation),fit$converged && fit$rank==10 && nobs(fit)==nrow(d) && all(is.finite(coef(fit))))
 X <- model.matrix(fit); mu <- fitted(fit); info <- crossprod(X,X*as.vector(mu*(1-mu)))
 diags[[length(diags)+1L]] <<- data.frame(dataset=set,imputation=imputation,outcome=outcome,N=nobs(fit),events=sum(d[[outcome]]),rank=fit$rank,df_residual=fit$df.residual,converged=fit$converged,iterations=fit$iter,min_probability=min(mu),max_probability=max(mu),information_condition=kappa(info,exact=TRUE),max_abs_score=max(abs(crossprod(X,d[[outcome]]-mu))),warnings=paste(warnings,collapse=' | '))
 ck(paste('No numerical separation warning',set,outcome,imputation),length(warnings)==0 && min(mu)>1e-8 && max(mu)<1-1e-8 && all(diag(vcov(fit))>0))
 fit
}
cc_rows <- list(); counts <- list(); center_counts <- list()
for(outcome in c('s_ge2','s_ge3')) {
 cc[[outcome]] <- as.integer(cc$steatosis>=if(outcome=='s_ge2')2 else 3)
 fit <- fit_model(cc,outcome,'CC')
 z <- summary(fit)$coefficients
 X <- model.matrix(fit); y <- cc[[outcome]]
 # Separate likelihood optimizer with analytic score; verify coefficients and Hessian SE.
 softplus <- function(x)pmax(x,0)+log1p(exp(-abs(x)))
 objective <- function(beta){eta<-as.vector(X%*%beta);sum(softplus(eta)-y*eta)}
 gradient <- function(beta)as.vector(crossprod(X,plogis(X%*%beta)-y))
 alt <- optim(rep(0,ncol(X)),objective,gradient,method='BFGS',control=list(maxit=20000,reltol=1e-13))
 mu2 <- as.vector(plogis(X%*%alt$par)); se2 <- sqrt(diag(solve(crossprod(X,X*(mu2*(1-mu2))))))
 ck(paste('Independent BFGS coefficients',outcome),alt$convergence==0 && max(abs(alt$par-coef(fit)))<2e-5)
 ck(paste('Independent Hessian standard errors',outcome),max(abs(se2-sqrt(diag(vcov(fit)))))<2e-5)
 bfgs[[outcome]] <- data.frame(outcome,coefficient=colnames(X),glm_beta=coef(fit),bfgs_beta=alt$par,glm_se=sqrt(diag(vcov(fit))),bfgs_se=se2)
 for(term in components) {
  q <- z[term,'Estimate']; se <- z[term,'Std. Error']
  cc_rows[[length(cc_rows)+1L]] <- data.frame(dataset='CC',outcome,term,N=nrow(cc),events=sum(y),m=NA_integer_,rank=fit$rank,estimate=q,std_error=se,df=Inf,fmi=NA_real_,OR=exp(q),ci_low=exp(q-1.96*se),ci_high=exp(q+1.96*se),p_value=z[term,'Pr(>|z|)'])
  for(value in 0:1) counts[[length(counts)+1L]] <- data.frame(dataset='CC',outcome,term,component_value=value,N=sum(cc[[term]]==value),events=sum(y[cc[[term]]==value]))
  for(center in levels(cc$center_code)) for(value in 0:1) {
   ix <- cc$center_code==center & cc[[term]]==value
   center_counts[[length(center_counts)+1L]] <- data.frame(center,outcome,term,component_value=value,N=sum(ix),events=sum(y[ix]))
  }
 }
}
# Reproduce frozen F>=2 component model as a linkage check (not a new analysis).
oldcc <- rd('component_models_complete_component.csv'); frozen <- oldcc[oldcc$outcome=='f_ge2',]
ref_fit <- fit_model(cc,'f_ge2','CC verification')
ck('Reproduces all five frozen CC fibrosis coefficients',max(abs(coef(ref_fit)[frozen$term]-frozen$estimate))<1e-8)
ck('Reproduces frozen CC fibrosis standard errors',max(abs(sqrt(diag(vcov(ref_fit)))[frozen$term]-frozen$std_error))<1e-8)
mi_cells <- list()
for(i in seq_len(imp$m)) {
 d <- mice::complete(imp,i)
 for(v in c('sex','dm','htn','alcohol','hcv')) d[[v]] <- as.integer(as.character(d[[v]]))
 d$sex_label <- factor(d$sex,levels=c(0,1),labels=c('Male','Female'))
 d$center_code <- factor(d$center_code,levels=c('HZ','JL','RJ'))
 d <- cbind(d,derive_cmrf(d))
 ck(paste('MI derived and pathology fixed',i),!anyNA(d[components]) && identical(d$steatosis,imp$data$steatosis) && all(d$cmrf_count>=1))
 for(outcome in c('s_ge2','s_ge3','f_ge2')) {
  d[[outcome]] <- if(outcome=='f_ge2')as.integer(d$fibrosis>=2) else as.integer(d$steatosis>=if(outcome=='s_ge2')2 else 3)
  fit <- fit_model(d,outcome,'MI',i)
  for(term in components) {
   z <- data.frame(dataset='MI',imputation=i,outcome,term,N=nobs(fit),events=sum(d[[outcome]]),rank=fit$rank,Q=unname(coef(fit)[term]),U=unname(vcov(fit)[term,term]))
   key <- paste(outcome,term,sep='|'); fit_store[[key]][[i]] <- z
   if(outcome!='f_ge2') estimates[[length(estimates)+1L]] <- z
   if(outcome!='f_ge2') for(value in 0:1)mi_cells[[length(mi_cells)+1L]] <- data.frame(imputation=i,outcome,term,component_value=value,N=sum(d[[term]]==value),events=sum(d[[outcome]][d[[term]]==value]))
  }
 }
}
mi_rows <- list(); pool_checks <- list()
oldmi <- read.csv(file.path(paper,'results/versions/20261003_165814/component_models_multiple_imputation.csv'))
for(outcome in c('s_ge2','s_ge3','f_ge2')) for(term in components) {
 z <- do.call(rbind,fit_store[[paste(outcome,term,sep='|')]])
 ck(paste('MI pooling 50 fixed N rank',outcome,term),nrow(z)==50 && length(unique(z$N))==1 && length(unique(z$rank))==1 && length(unique(z$events))==1)
 p <- pool_scalar_rubin(z$Q,z$U,n=z$N[1],k=z$rank[1])
 reference <- mice::pool.scalar(z$Q,z$U,n=z$N[1],k=z$rank[1])
 delta <- max(abs(c(p['estimate']-reference$qbar,p['std_error']-sqrt(reference$t),p['df']-reference$df,p['fmi']-reference$fmi)))
 ck(paste('Rubin vs mice',outcome,term),delta<1e-9)
 pool_checks[[length(pool_checks)+1L]] <- data.frame(outcome,term,max_difference=delta,k=z$rank[1],dfcom=z$N[1]-z$rank[1],df=p['df'])
 if(outcome=='f_ge2') {
  old <- oldmi[oldmi$outcome==outcome & oldmi$term==term,]
  ck(paste('Reproduce frozen corrected MI fibrosis',term),max(abs(p-unlist(old[names(p)])))<1e-7)
 } else {
  row <- data.frame(dataset='MI',outcome,term,N=z$N[1],events=z$events[1],m=50,rank=z$rank[1],estimate=p['estimate'],std_error=p['std_error'],df=p['df'],fmi=p['fmi'],OR=p['OR'],ci_low=p['ci_low'],ci_high=p['ci_high'],p_value=p['p_value'])
  rownames(row)<-NULL;mi_rows[[length(mi_rows)+1L]] <- row
 }
}
results <- rbind(do.call(rbind,cc_rows),do.call(rbind,mi_rows))
results$p_bh_10 <- NA_real_
for(set in c('CC','MI')) {
 ix <- results$dataset==set; results$p_bh_10[ix] <- p.adjust(results$p_value[ix],method='BH')
 ck(paste('BH10',set),sum(ix)==10 && all(results$p_bh_10[ix]>=results$p_value[ix]))
}
results$supported <- results$p_bh_10<.05
results$comparison <- 'Component present versus absent; conditional on other four domains'
results$adjustment <- 'All five CMRF components, age, sex and three-centre fixed effects; no additional BMI'
results$analysis_role <- 'Additional secondary component analysis specified 2026-10-23 before model fitting'
ck('20 finite OR intervals',nrow(results)==20 && all(is.finite(as.matrix(results[c('OR','ci_low','ci_high','p_value','p_bh_10')]))) && all(results$ci_low<results$OR & results$ci_high>results$OR))
ck('CC outcomes 645 and 301',all(results$events[results$dataset=='CC' & results$outcome=='s_ge2']==645) && all(results$events[results$dataset=='CC' & results$outcome=='s_ge3']==301))
tables <- list(steatosis_component_models=results,component_outcome_counts_CC=do.call(rbind,counts),component_outcome_counts_by_center_CC=do.call(rbind,center_counts),component_outcome_counts_MI=do.call(rbind,mi_cells),per_imputation_estimates=do.call(rbind,estimates))
for(name in names(tables))write.csv(tables[[name]],file.path(out,'source_data',paste0(name,'.csv')),row.names=FALSE)
write.csv(do.call(rbind,diags),file.path(out,'qa/model_diagnostics.csv'),row.names=FALSE)
write.csv(do.call(rbind,bfgs),file.path(out,'qa/independent_BFGS_checks.csv'),row.names=FALSE)
write.csv(do.call(rbind,pool_checks),file.path(out,'qa/pooling_checks.csv'),row.names=FALSE)
write.csv(data.frame(check=names(checks),pass=unlist(checks)),file.path(out,'qa/analysis_checks.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'qa/sessionInfo_analysis.txt'))
cat(length(checks),'analysis checks passed. New component associations:\n')
print(results[c('dataset','outcome','term','N','events','OR','ci_low','ci_high','p_value','p_bh_10')],row.names=FALSE)
