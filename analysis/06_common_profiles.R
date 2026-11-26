# Complete CMRF patterns versus all other participants. Author: mjwang.
arg <- grep('^--file=', commandArgs(FALSE), value = TRUE)
script <- normalizePath(sub('^--file=', '', arg))
root <- normalizePath(file.path(dirname(script), '../../../../..'))
version <- basename(dirname(script))
out <- file.path(root, 'paper1_cmrf_pathology/results/versions', version)
dest <- file.path(out, 'source_data')
target <- file.path(dest, 'pattern_vs_rest_2x2.csv')
if (file.exists(target)) stop('Use a new version; existing results are read-only.')
d <- read.csv(file.path(root, 'paper1_cmrf_pathology/results/versions/20260809_152651/canonical_adult_complete_component.csv'))
stopifnot(nrow(d) == 935, sum(d$f_ge2) == 550, sum(d$f_ge3) == 246)
cm <- as.matrix(d[paste0('cmrf', 1:5)])
stopifnot(!anyNA(cm), all(cm %in% c(0, 1)), all(rowSums(cm) == d$cmrf_count))
d$analysis_center <- ifelse(d$source_prefix == 'ZS', 'RJ', d$center_code)
stopifnot(setequal(unique(d$analysis_center), c('HZ', 'JL', 'RJ')))
d$pattern <- apply(cm, 1, function(z) paste(c('A', 'G', 'H', 'T', 'L')[z == 1], collapse = '+'))
key <- aggregate(id ~ cmrf_count + pattern, d, length)
names(key)[3] <- 'n'
key <- key[order(key$cmrf_count, -key$n, key$pattern), ]
write.csv(key, file.path(dest, 'all_21_patterns.csv'), row.names = FALSE)
key <- key[key$n >= 30, ]
stopifnot(nrow(key) == 11, sum(key$n) == 839)
key$order <- seq_len(nrow(key))
wilson <- function(e, n) {
  z <- qnorm(.975); pr <- e/n; den <- 1 + z*z/n
  mid <- (pr + z*z/(2*n))/den
  half <- z * sqrt(pr*(1-pr)/n + z*z/(4*n*n))/den
  c(lo = mid-half, hi = mid+half)
}
rows <- list()
for (ep in c('f_ge2', 'f_ge3')) for (i in seq_len(nrow(key))) {
  k <- key[i, ]; member <- d$pattern == k$pattern; y <- d[[ep]]
  a <- sum(member & y == 1); b <- sum(member & y == 0)
  c <- sum(!member & y == 1); dd <- sum(!member & y == 0)
  tab <- matrix(c(a,b,c,dd), nrow=2, byrow=TRUE,
                dimnames=list(c('This pattern','All others'),c('Endpoint','No endpoint')))
  chi <- suppressWarnings(chisq.test(tab, correct=FALSE))
  use_fisher <- any(chi$expected < 5)
  test <- if(use_fisher) fisher.test(tab, alternative='two.sided') else chi
  or <- a*dd/(b*c)
  se <- sqrt(1/a + 1/b + 1/c + 1/dd)
  ci <- exp(log(or) + c(-1,1)*qnorm(.975)*se)
  w1 <- wilson(a,a+b); w0 <- wilson(c,c+dd)
  p <- test$p.value
  rows[[length(rows)+1]] <- data.frame(
    endpoint=if(ep=='f_ge2') 'F>=2' else 'F>=3', pattern=k$pattern,
    cmrf_count=k$cmrf_count, order=k$order, N_total=nrow(d),
    n_pattern=a+b, n_other=c+dd,
    a_pattern_endpoint=a, b_pattern_no_endpoint=b,
    c_other_endpoint=c, d_other_no_endpoint=dd,
    proportion_pattern=a/(a+b), lower_pattern=w1['lo'], upper_pattern=w1['hi'],
    proportion_other=c/(c+dd), lower_other=w0['lo'], upper_other=w0['hi'],
    difference_percentage_points=100*(a/(a+b)-c/(c+dd)),
    OR_unadjusted=or, OR_lower=ci[1], OR_upper=ci[2],
    OR_interval='Log-Wald 95% CI',
    chi_square=unname(chi$statistic), df=unname(chi$parameter),
    expected_a=chi$expected[1,1], expected_b=chi$expected[1,2],
    expected_c=chi$expected[2,1], expected_d=chi$expected[2,2],
    min_expected=min(chi$expected),
    test=if(use_fisher) 'Fisher exact, two-sided' else 'Pearson chi-square, no Yates correction',
    p_raw=p,
    significance=if(p<.001) '***' else if(p<.01) '**' else if(p<.05) '*' else 'ns',
    covariate_adjustment='None', multiplicity_adjustment='None',
    contrast='Exact CMRF pattern versus all other participants',
    interpretation='Exploratory unadjusted association', row.names=NULL)
}
res <- do.call(rbind, rows)
stopifnot(nrow(res)==22, all(res$n_pattern+res$n_other==935), all(res$min_expected>=5))
prior_path <- file.path(out,'inputs/prior_pattern_values.csv')
if(file.exists(prior_path)) {
prior <- read.csv(prior_path)
for(i in seq_len(nrow(res))) {
  ep <- if(res$endpoint[i]=='F>=2') 'F2' else 'F3'
  p <- prior[prior$pattern==res$pattern[i] & prior$outcome==ep,]
  stopifnot(nrow(p)==1, p$n==res$n_pattern[i], p$events==res$a_pattern_endpoint[i],
            abs(p$prop-res$proportion_pattern[i])<1e-12,
            abs(p$lo-res$lower_pattern[i])<1e-12, abs(p$hi-res$upper_pattern[i])<1e-12)
}
} else warning('Historical figure parity file is absent; the computed comparisons are unchanged.')
write.csv(res, target, row.names=FALSE, na='')
write.csv(key, file.path(dest,'plotted_patterns.csv'),row.names=FALSE)
write.csv(as.data.frame(table(d$analysis_center)),file.path(dest,'center_counts.csv'),row.names=FALSE)
qc <- data.frame(N=935,F_ge2=550,F_ge3=246,patterns_all=nrow(unique(d['pattern'])),
                 patterns_plotted=11,participants_in_plotted_patterns=839,
                 participants_in_other_rare_patterns=96,comparisons=22,
                 Pearson_tests=sum(grepl('Pearson',res$test)),Fisher_tests=sum(grepl('Fisher',res$test)),
                 minimum_expected=min(res$min_expected),excluded_participants=0,
                 original_plotted_values_preserved=TRUE)
write.csv(qc,file.path(dest,'analysis_qc.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(dest,'R_sessionInfo.txt'))
print(res[c('endpoint','pattern','n_pattern','a_pattern_endpoint','c_other_endpoint','proportion_pattern','proportion_other','p_raw','significance')],row.names=FALSE)
print(qc)
