# Aggregate-data figure revision and requested descriptive comparisons.
# Revision author: mjwang. Run with project root and phase (draw or pdf).
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
root <- normalizePath(args[1],mustWork=TRUE)
out <- file.path(root,'paper1_cmrf_pathology/results/versions/20261025_124331')
library(grid)
run <- function(command,args){
 z<-system2(command,shQuote(args),stdout=TRUE,stderr=TRUE)
 if(!is.null(attr(z,'status')) && attr(z,'status')!=0)stop(paste(z,collapse='\n'))
}
if(args[2]=='draw'){
 library(svglite);library(ragg);library(magick)
 d<-read.csv(file.path(out,'source_data/robustness_models.csv'))
 stopifnot(nrow(d)==30,all(d$ci_low<d$OR),all(d$ci_high>d$OR))
 labs<-c('Primary (complete-component)','Raw-variable MI (50 datasets)',
 'All ages (including <18)','Restricted to S≥1','Legacy CMRF count',
 'Uniform HDL-C <1.0','Non-adiposity count + BMI','Recorded alcohol adjustment',
 'Exclude non-zero HCV','Age natural spline','Leave out HZ','Leave out JL',
 'Leave out RJ','Indeterminate domains: negative','Indeterminate domains: positive')
 W<-183;H<-175;bounds<-list()
 textmm<-function(label,x,y,bold=FALSE,col='#202832',just='left',size=9){
  grid.text(label,unit(x,'mm'),unit(y,'mm'),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))
  pushViewport(viewport(gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain')))
  w<-convertWidth(stringWidth(label),'mm',valueOnly=TRUE);popViewport()
  lo<-x-if(just=='right')w else if(just=='centre')w/2 else 0
  bounds[[paste(label,x,y)]]<<-data.frame(label,x_low=lo,x_high=lo+w,y,size)
 }
 line<-function(x,y,col='#b7c1c9',lty=1,lwd=.7){grid.lines(unit(x,'mm'),unit(y,'mm'),gp=gpar(col=col,lty=lty,lwd=lwd))}
 pformat<-function(p)if(p<.001)formatC(p,format='e',digits=2) else formatC(p,format='fg',digits=3,flag='#')
 draw<-function(){
  grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
  textmm('Sensitivity of the CMRF count–fibrosis association',3,168,TRUE,size=11)
  textmm('Adjusted ORs with 95% CIs; P values are unadjusted.',3,161,col='#586572')
  textmm('Analysis scenario',3,149,TRUE)
  for(j in 1:2){
   x0<-c(69,130)[j];x1<-c(102,163)[j];xp<-c(119,180)[j]
   col<-c('#0072B5','#BC3C29')[j];lim<-c(.8,2.25)
   px<-function(v)x0+diff(c(x0,x1))*log(v/lim[1])/log(lim[2]/lim[1])
   textmm(paste0('F≥',j+1),mean(c(x0,x1)),151,TRUE,col,'centre',10)
   textmm('P',xp,149,TRUE,just='right')
   line(rep(px(1),2),c(27,145),'#6d7a84',2)
   z<-d[d$outcome==paste0('f_ge',j+1),];z<-z[order(z$scenario_order),]
   stopifnot(identical(as.integer(z$scenario_order),1:15),all(z$ci_low>=lim[1]),all(z$ci_high<=lim[2]))
   for(k in 1:15){
    y<-140-(k-1)*7.8
    if(j==1){
     if(k%%2==0)grid.rect(unit(W/2,'mm'),unit(y,'mm'),width=unit(179,'mm'),height=unit(7.8,'mm'),gp=gpar(fill='#f4f6f8',col=NA))
     textmm(labs[k],3,y,bold=k<=2)
    }
    line(c(px(z$ci_low[k]),px(z$ci_high[k])),rep(y,2),col,lwd=1)
    for(xx in px(c(z$ci_low[k],z$ci_high[k])))line(rep(xx,2),y+c(-.65,.65),col)
    grid.points(unit(px(z$OR[k]),'mm'),unit(y,'mm'),pch=if(j==1)16 else 17,size=unit(2.1,'mm'),gp=gpar(col=col))
    textmm(pformat(z$p_value[k]),xp,y,col=col,just='right')
   }
   line(c(x0,x1),c(23,23),'#202832')
   for(v in c(.8,1,1.5,2)){xx<-px(v);line(rep(xx,2),c(23,21.5),'#202832');textmm(as.character(v),xx,18,just='centre')}
  }
  textmm('Odds ratio per additional CMRF (log scale)',124,9,just='centre')
 }
 stem<-file.path(out,'figures/FigureS1')
 stopifnot(!file.exists(paste0(stem,'.pdf')))
 cairo_pdf(paste0(stem,'.pdf'),width=W/25.4,height=H/25.4,family='Arial');draw();invisible(dev.off())
 svglite(paste0(stem,'.svg'),width=W/25.4,height=H/25.4,bg='white');draw();invisible(dev.off())
 agg_png(paste0(stem,'.png'),width=W,height=H,units='mm',res=300,background='white');draw();invisible(dev.off())
 agg_tiff(paste0(stem,'.tiff'),width=W,height=H,units='mm',res=300,compression='lzw',background='white');draw();invisible(dev.off())
 b<-do.call(rbind,bounds);stopifnot(all(b$x_low>=0),all(b$x_high<=W),all(b$size>=9))
 write.csv(b,file.path(out,'qa/FigureS1_text_bounds.csv'),row.names=FALSE)
 d$display_p<-vapply(d$p_value,pformat,character(1));d$display_label<-labs[d$scenario_order]
 write.csv(d,file.path(out,'source_data/FigureS1_display_values.csv'),row.names=FALSE)
 image_write(image_resize(image_read(paste0(stem,'.png')),'1600x'),file.path(out,'qa/FigureS1_preview.png'))
 image_write(image_convert(image_read(paste0(stem,'.png')),colorspace='gray'),file.path(out,'qa/FigureS1_grayscale.png'))
 # Descriptive centre tests requested for the transposed Table S2.
 mix<-read.csv(file.path(out,'source_data/center_fibrosis_case_mix.csv'))
 stopifnot(identical(sort(unique(mix$center_code)),c('HZ','JL','RJ')),sum(mix$n)==935)
 mat<-xtabs(n~center_code+fibrosis_stage,mix)
 binary<-function(k)cbind(below=rowSums(mat[,as.numeric(colnames(mat))<k,drop=FALSE]),above=rowSums(mat[,as.numeric(colnames(mat))>=k,drop=FALSE]))
 tables<-list(GLOBAL_STAGE=mat,GLOBAL_F2=binary(2),GLOBAL_F3=binary(3))
 results<-do.call(rbind,lapply(names(tables),function(name){
  z<-chisq.test(tables[[name]],correct=FALSE)
  stopifnot(min(z$expected)>=5)
  data.frame(comparison=name,N=sum(tables[[name]]),statistic=unname(z$statistic),df=unname(z$parameter),p_value=z$p.value,min_expected=min(z$expected),p_display=formatC(z$p.value,format='e',digits=2),test='Pearson chi-square, two-sided, unadjusted')
 }))
 write.csv(results,file.path(out,'source_data/centre_comparison_tests.csv'),row.names=FALSE)
 write.csv(data.frame(check=c('30 source estimates','15 scenario labels','all confidence intervals shown','minimum9pt','text within canvas','three centres only','expected counts at least5'),pass=TRUE),file.path(out,'qa/R_checks.csv'),row.names=FALSE)
 capture.output(sessionInfo(),file=file.path(out,'qa/R_sessionInfo.txt'))
 print(results)
} else if(args[2]=='pdf'){
 model<-jsonlite::fromJSON(file.path(out,'review_items.json'),simplifyVector=FALSE)
 docpdf<-file.path(out,'qa/render_tables/Supplementary_Tables_S1-S4.pdf')
 stopifnot(file.exists(docpdf))
 # Four independently paginated editable tables must each occupy one page.
 run(Sys.which('pdfseparate'),c(docpdf,file.path(out,'pdf_parts/NewSupp-%d.pdf')))
 for(n in 2:4)file.copy(file.path(out,paste0('pdf_parts/NewSupp-',n,'.pdf')),file.path(out,paste0('pdf_parts/TableS',n,'.pdf')),overwrite=TRUE)
 cover<-file.path(out,'pdf_parts/Cover.pdf')
 pdf(cover,width=8.27,height=11.69,title='Paper 1 revised figures and tables',author='mjwang',family='Helvetica',useDingbats=FALSE)
 par(mar=c(0,0,0,0));plot.new();plot.window(xlim=c(0,1),ylim=c(0,1),xaxs='i',yaxs='i')
 text(.07,.94,'Paper 1',adj=c(0,1),cex=1.8,font=2,col='#183b58')
 text(.07,.88,'Revised figures and tables',adj=c(0,1),cex=1.45,font=2)
 text(.07,.83,'3 main figures / 1 main table / 5 supplementary figures / 4 supplementary tables',adj=0,cex=.72)
 text(.07,.79,'Eight items retain author confirmation; five revised items await review.',adj=0,cex=.78,col='#546875')
 text(.07,.735,'Item',adj=0,font=2);text(.91,.735,'PDF page',adj=1,font=2)
 index<-list()
 for(n in seq_along(model$items)){
  i<-model$items[[n]];y<-.70-(n-1)*.037
  label<-gsub('FigureS','Figure S',i$id,fixed=TRUE);label<-gsub('TableS','Table S',label,fixed=TRUE)
  if(i$id%in%c('Figure1','Figure2','Figure3','Table1'))label<-gsub('([0-9])',' \\1',label)
  text(.07,y,label,adj=0,cex=.79)
  text(.29,y,if(nchar(i$title)>61)paste0(substr(i$title,1,58),'...') else i$title,adj=0,cex=.60)
  text(.91,y,as.character(n+1),adj=1,cex=.80)
  index[[n]]<-data.frame(id=i$id,pdf_page=n+1)
 }
 text(.07,.16,'Figure 2 contains fibrosis only. Steatosis estimates remain in Figure S3b.',adj=0,cex=.75)
 text(.07,.125,'Original Table S12 is now compact Table S4; all 20 contrasts are retained.',adj=0,cex=.75)
 text(.07,.09,'Full numerical source tables remain available in the offline package.',adj=0,cex=.75)
 text(.07,.045,'Author review version / 25 October 2026 / revision author mjwang',adj=0,cex=.74,col='#546875')
 invisible(dev.off())
 inputs<-c(cover)
 for(i in model$items)inputs<-c(inputs,if(i$kind=='figure')file.path(out,i$pdf) else file.path(out,'pdf_parts',paste0(i$id,'.pdf')))
 run(Sys.which('pdfunite'),c(inputs,file.path(out,'Paper1_All_Figures_Tables.pdf')))
 run(Sys.which('pdfunite'),c(vapply(model$items[1:3],function(i)file.path(out,i$pdf),character(1)),file.path(out,'Main_Figures_1-3.pdf')))
 figs<-Filter(function(i)i$kind=='figure'&&i$group=='supplementary',model$items)
 run(Sys.which('pdfunite'),c(vapply(figs,function(i)file.path(out,i$pdf),character(1)),file.path(out,'Supplementary_Figures_S1-S5.pdf')))
 render<-Sys.getenv('CMRF_PDFTOPPM');stopifnot(nzchar(render),file.exists(render))
 dir.create(file.path(out,'qa/pdf_pages'),showWarnings=FALSE)
 run(render,c('-r','110','-png',file.path(out,'Paper1_All_Figures_Tables.pdf'),file.path(out,'qa/pdf_pages/page')))
 write.csv(do.call(rbind,index),file.path(out,'pdf_index.csv'),row.names=FALSE)
 cat('Assembled 14 pages including the index.\n')
} else stop('Unknown phase')
