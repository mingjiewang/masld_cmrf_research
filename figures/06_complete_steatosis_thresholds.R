#!/usr/bin/env Rscript
# Supplementary Figure S4; panel relabelled a; quantitative source only. Author: mjwang.
suppressPackageStartupMessages({library(grid);library(ggsci);library(svglite);library(ragg);library(magick)})
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=1)stop('Pass project root')
root<-normalizePath(args[1]);version<-'20261024_003114'
out<-file.path(root,'paper1_cmrf_pathology/results/versions',version)
d<-read.csv(file.path(out,'source_data/steatosis_component_models.csv'))
stopifnot(nrow(d)==20,all(d$ci_low>.08 & d$ci_high<3.2))
W<-183;H<-174;pal<-ggsci::pal_nejm('default')(8);cols<-pal[c(2,1)]
ink<-'#202832';grey<-'#586572';geometry<-list();drawn<-list()
txt<-function(label,x,y,size=9,bold=FALSE,col=ink,just='centre')grid.text(label,x=unit(x,'mm'),y=unit(y,'mm'),just=just,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.05))
ln<-function(x1,y1,x2,y2,col='#A5AFB7',lwd=.6,lty=1)grid.lines(x=unit(c(x1,x2),'mm'),y=unit(c(y1,y2),'mm'),gp=gpar(col=col,lwd=lwd,lty=lty))
pt<-function(x,y,col,pch,fill='white')grid.points(x=unit(x,'mm'),y=unit(y,'mm'),pch=pch,size=unit(2.2,'mm'),gp=gpar(col=col,fill=fill,lwd=.9))
measure<-function(z){pushViewport(viewport(gp=gpar(fontfamily='Arial',fontsize=9)));w<-convertWidth(stringWidth(z),'mm',TRUE);popViewport();w}
label<-function(key,s,x,y,maxwidth){txt(s,x,y);geometry[[key]]<<-data.frame(key,label=s,width_mm=measure(s),available_mm=maxwidth)}
pfmt<-function(p)if(p<.001)sprintf('%.2e',p) else sprintf('%.3f',p)
efmt<-function(z)sprintf('%.2f (%.2f–%.2f)',z$OR,z$ci_low,z$ci_high)
lim<-c(.08,3.2);xp<-function(v,x0,x1)x0+(x1-x0)*(log(v)-log(lim[1]))/diff(log(lim))
draw<-function(){
 grid.newpage();grid.rect(gp=gpar(fill='white',col=NA))
 txt('a',3,167,12,TRUE,just='left');txt('CMRF component associations with histological steatosis',11,167,10.5,TRUE,just='left')
 txt('Mutually adjusted for all five domains, age, sex and centre',11,160,col=grey,just='left')
 pt(5,153,ink,21);txt('CC: complete-component (n = 935)',9,153,just='left')
 pt(100,153,ink,24);txt('MI: imputed (n = 1,085; m = 50)',104,153,just='left')
 plotx<-list(c(40,66),c(114,140));ox<-c(82,156);pp<-c(103,176)
 txt('Component',3,129,9,TRUE,just='left');txt('Set',33,129,9,TRUE)
 for(j in 1:2){
  x0<-plotx[[j]][1];x1<-plotx[[j]][2]
  txt(c('Moderate–severe steatosis (S≥2)','Severe steatosis (S3)')[j],mean(c(x0,pp[j])),143,9.5,TRUE,col=cols[j])
  txt(c('CC: 645/935 · MI: 742/1,085','CC: 301/935 · MI: 343/1,085')[j],mean(c(x0,pp[j])),136,col=grey)
  txt('OR (95% CI)',ox[j],129,9,TRUE);txt('P (BH)',pp[j],129,9,TRUE)
  ln(xp(1,x0,x1),31,xp(1,x0,x1),122,grey,.75,2)
  for(i in 1:5){
   cy<-c(114,95,76,57,38)[i]
   if(j==1)txt(c('Adiposity','Glycaemia','Hypertension','High triglycerides','Low HDL-C')[i],3,cy,just='left')
   if(i<5)ln(x0-2,cy-9.5,pp[j]+6,cy-9.5,'#EDF0F2',.45)
   for(k in 1:2){
    z<-d[d$dataset==c('CC','MI')[k] & d$outcome==c('s_ge2','s_ge3')[j] & d$term==paste0('cmrf',i),];stopifnot(nrow(z)==1)
    y<-cy+if(k==1)3.1 else -3.1
    if(j==1)txt(c('CC','MI')[k],33,y,col=grey)
    a<-xp(z$ci_low,x0,x1);b<-xp(z$ci_high,x0,x1);x<-xp(z$OR,x0,x1)
    ln(a,y,b,y,cols[j],1);ln(a,y-.7,a,y+.7,cols[j],.7);ln(b,y-.7,b,y+.7,cols[j],.7)
    pt(x,y,cols[j],if(k==1)21 else 24,if(z$supported)cols[j] else 'white')
    key<-paste(j,i,k,sep='|');drawn[[key]]<<-data.frame(key,OR=z$OR,ci_low=z$ci_low,ci_high=z$ci_high,x_mm=x,low_mm=a,high_mm=b,plot_low=lim[1],plot_high=lim[2],supported=z$supported)
    label(paste0('OR',key),efmt(z),ox[j],y,26.5);label(paste0('P',key),pfmt(z$p_bh_10),pp[j],y,13.5)
   }
  }
  ln(x0,27,x1,27,ink)
  for(v in c(.1,.5,1,2)){x<-xp(v,x0,x1);ln(x,27,x,25.8,ink);txt(format(v,trim=TRUE),x,22.8)}
  txt('Adjusted OR (log scale)',mean(c(x0,pp[j])),16)
 }
 txt('Filled symbols: BH-adjusted P <0.05. Intervals: nominal 95% CI.',3,7,col=grey,just='left')
}
base<-file.path(out,'supplementary/FigureS4_steatosis_components')
if(file.exists(paste0(base,'.pdf')))stop('Figure already exists')
cairo_pdf(paste0(base,'.pdf'),width=W/25.4,height=H/25.4,family='Arial');draw();dev.off()
svglite(paste0(base,'.svg'),width=W/25.4,height=H/25.4,bg='white');draw();dev.off()
agg_png(paste0(base,'.png'),width=W,height=H,units='mm',res=300,background='white');draw();dev.off()
agg_tiff(paste0(base,'.tiff'),width=W,height=H,units='mm',res=300,compression='lzw',background='white');draw();dev.off()
g<-do.call(rbind,geometry);f<-do.call(rbind,drawn)
stopifnot(all(g$width_mm<g$available_mm),nrow(f)==20,all(f$ci_low>=f$plot_low & f$ci_high<=f$plot_high))
write.csv(g,file.path(out,'qa/figure_label_widths.csv'),row.names=FALSE)
write.csv(f,file.path(out,'qa/figure_geometry.csv'),row.names=FALSE)
im<-image_read(paste0(base,'.png'));image_write(image_convert(im,colorspace='gray'),file.path(out,'qa/FigureS4_grayscale.png'))
image_write(image_resize(im,'1500x'),file.path(out,'qa/FigureS4_preview.png'))
write.csv(data.frame(outcome=c('S≥2','S3'),hex=cols,shape_CC='circle',shape_MI='triangle'),file.path(out,'source_data/figure_palette.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'qa/sessionInfo_figure.txt'))
cat('All 20 model estimates rendered; label widths and CI limits checked.\n')
