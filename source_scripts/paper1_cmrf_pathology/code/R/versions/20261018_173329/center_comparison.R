#!/usr/bin/env Rscript
# Descriptive overall comparisons for the existing 1085-person center table.
# Author: mjwang. Does not change eligibility, imputation or regression models.
script <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)))
paper <- normalizePath(file.path(dirname(script), '../../../..'))
version <- basename(dirname(script))
out <- file.path(paper, 'results/versions', version)
input <- file.path(paper, 'results/versions/20260809_152651/eligible_adult_for_mi.csv')
if (file.exists(file.path(out, 'center_comparisons_R.csv'))) stop('Output exists; use a new version')
d <- read.csv(input, check.names=FALSE)
d$analysis_center <- ifelse(d$source_prefix == 'ZS', 'RJ', d$source_prefix)
stopifnot(nrow(d)==1085, !anyDuplicated(d$id), all(d$analysis_center==d$center_code))
d$analysis_center <- factor(d$analysis_center, levels=c('RJ','HZ','JL'))
d$indeterminate <- as.integer(!d$complete_component)
stopifnot(all(table(d$analysis_center)==c(369,577,139)))
variables <- c('age','bmi','indeterminate','dm','f_ge2','f_ge3')
labels <- c('年龄','BMI','组分不确定','已确诊T2DM','F≥2','F≥3')
results <- list(); summaries <- list(); expected <- list()
for (i in seq_along(variables)) {
  v <- variables[i]; y <- d[[v]]
  stopifnot(!anyNA(y))
  if (v %in% c('age','bmi')) {
    fit <- oneway.test(y ~ analysis_center, data=d, var.equal=FALSE)
    method <- 'Welch one-way ANOVA'
    df1 <- unname(fit$parameter[1]); df2 <- unname(fit$parameter[2]); min_exp <- NA_real_
  } else {
    stopifnot(all(y %in% c(0,1)))
    tab <- table(d$analysis_center, factor(y,levels=c(0,1)))
    fit <- chisq.test(tab, correct=FALSE)
    min_exp <- min(fit$expected); stopifnot(min_exp>=5)
    method <- 'Pearson chi-square (3x2, no continuity correction)'
    df1 <- unname(fit$parameter); df2 <- NA_real_
    expected[[v]] <- data.frame(variable=v,center=rep(rownames(tab),2),
                               category=rep(colnames(tab),each=3),observed=as.vector(tab),
                               expected=as.vector(fit$expected))
  }
  results[[v]] <- data.frame(variable=v,label=labels[i],n=length(y),method=method,
                             statistic=unname(fit$statistic),df1=df1,df2=df2,
                             p_value=fit$p.value,min_expected=min_exp)
  for (cc in levels(d$analysis_center)) {
    vals <- y[d$analysis_center==cc]
    summaries[[paste(v,cc)]] <- data.frame(variable=v,center=cc,n=length(vals),
      mean=mean(vals),sd=sd(vals),positive_n=if(v %in% c('age','bmi')) NA_real_ else sum(vals),
      display=if(v %in% c('age','bmi')) sprintf('%.1f±%.1f',mean(vals),sd(vals)) else sprintf('%d (%.1f)',sum(vals),100*mean(vals)))
  }
}
res <- do.call(rbind,results)
res$p_holm_6 <- p.adjust(res$p_value,method='holm')
res$p_display <- ifelse(res$p_value<0.001,'<0.001',sprintf('%.3f',res$p_value))
write.csv(res,file.path(out,'center_comparisons_R.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(do.call(rbind,summaries),file.path(out,'center_descriptives_R.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(do.call(rbind,expected),file.path(out,'contingency_expected_R.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
print(res[,c('variable','statistic','p_value','p_holm_6','min_expected')],row.names=FALSE)
