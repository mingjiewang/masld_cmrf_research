# Figure 2 revision; revision author: mjwang. Frozen summary data only.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1)stop('Pass project root')
root<-normalizePath(args[1]);version<-'20261024_094342'
out<-file.path(root,'paper1_cmrf_pathology/results/versions',version)
library(grid);library(svglite);library(ragg);library(ggsci)
rd<-function(n)read.csv(file.path(out,'source_data',n))
st<-rd('fibrosis_distribution_by_cmrf.csv');obs<-rd('observed_fibrosis_wilson.csv')
pr<-rd('standardized_probabilities.csv');mods<-rd('main_models_complete_component.csv')
ink<-'#202832';grey<-'#586572';palette<-pal_nejm('default')(8)
stage_col<-palette[c(6,2,7,3,1)];stage_text<-c(ink,'white',ink,ink,'white')
outcome_col<-c(palette[2],palette[1]);outcome_text<-outcome_col;width_scale<-1;H<-246;y_offset<-0;panel_part<-'upper'
fy<-function(y)if(panel_part=='upper')115+(y-82.8)*125/131.2 else if(panel_part=='lower')y+31 else y
fh<-function(h)if(panel_part=='upper')h*125/131.2 else h
fmt<-function(p)sprintf('%.1f',100*p)
checks<-list()
ck<-function(name,pass){checks[[name]]<<-isTRUE(pass);if(!isTRUE(pass))stop(name)}
N<-sapply(1:5,function(k)sum(st$n[st$cmrf_count==k]))
ck('five denominators',identical(as.numeric(N),c(40,188,313,292,102)))
ck('total935',sum(N)==935)
ck('25 stage cells',nrow(st)==25)
for(k in 1:5){z<-st[st$cmrf_count==k,];ck(paste('stage percentages',k),max(abs(z$pct-100*z$n/N[k]))<1e-9)}
for(i in 1:nrow(obs)){
 z<-obs[i,];q<-qnorm(.975);p<-z$events/z$N
 den<-1+q^2/z$N;mid<-(p+q^2/(2*z$N))/den;half<-q*sqrt(p*(1-p)/z$N+q^2/(4*z$N^2))/den
 ck(paste('Wilson row',i),max(abs(c(p,mid-half,mid+half)-unlist(z[c('probability','ci_low','ci_high')])))<1e-12)
 ck(paste('nested stage events',i),z$events==sum(st$n[st$cmrf_count==z$count & st$fibrosis_stage>=if(z$outcome=='F≥2')2 else 3]))
}
for(i in 1:2){z<-mods[i,];actual<-c(exp(z$estimate),exp(z$estimate+c(-1,1)*1.96*z$std_error),2*pnorm(-abs(z$estimate/z$std_error)))
 ck(paste('OR CI P row',i),max(abs(actual-c(z$OR,z$ci_low,z$ci_high,z$p_value)))<1e-12)}
ck('bootstrap1000',all(pr$bootstrap_success==1000 & pr$bootstrap_replicates==1000))
ck('standardization935',all(pr$N==935))
ck('F3 raw dip retained',obs$probability[obs$outcome=='F≥3' & obs$count==4]<obs$probability[obs$outcome=='F≥3' & obs$count==3])
text_boxes<-list()
txt<-function(label,x,y,size=9,bold=FALSE,col=ink,just='centre',rot=0){
 grid.text(label,x=unit(x*width_scale,'mm'),y=unit(fy(y)-y_offset,'mm'),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col))
 pushViewport(viewport(gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain')))
 w<-convertWidth(stringWidth(label),'mm',TRUE);h<-convertHeight(stringHeight(label),'mm',TRUE);popViewport()
 if(rot==90){t<-w;w<-h;h<-t}
 xx<-x*width_scale;yy<-fy(y)-y_offset
 lowx<-xx-if(just=='left')0 else if(just=='right')w else w/2
 key<-paste(panel_part,label,x,y,size,just,rot)
 text_boxes[[key]]<<-data.frame(panel=panel_part,label=label,size_pt=size,x_low=lowx,x_high=lowx+w,y_low=yy-h/2,y_high=yy+h/2,canvas_width=183*width_scale)
}

ln<-function(x1,y1,x2,y2,col='#A5AFB7',lwd=.6){grid.lines(x=unit(c(x1,x2)*width_scale,'mm'),y=unit(fy(c(y1,y2))-y_offset,'mm'),gp=gpar(col=col,lwd=lwd))}
rc<-function(x,y,w,h,col,border=NA){grid.rect(x=unit(x*width_scale,'mm'),y=unit(fy(y)-y_offset,'mm'),width=unit(w*width_scale,'mm'),height=unit(fh(h),'mm'),gp=gpar(fill=col,col=border,lwd=.45))}
pt<-function(x,y,col,pch=16,fill=col){grid.points(x=unit(x*width_scale,'mm'),y=unit(fy(y)-y_offset,'mm'),pch=pch,size=unit(2.1,'mm'),gp=gpar(col=col,fill=fill,lwd=.9))}
poly<-function(x,y,col){grid.polygon(x=unit(x*width_scale,'mm'),y=unit(fy(y)-y_offset,'mm'),gp=gpar(fill=col,col=NA))}
measure<-function(s,bold=FALSE){pushViewport(viewport(gp=gpar(fontfamily='Arial',fontsize=9,fontface=if(bold)'bold' else 'plain')));w<-convertWidth(stringWidth(s),'mm',TRUE);popViewport();w}
geometry<-list();bar_label_boxes<-list()
pct_label<-function(label,x,y,col,bar_x){
 w<-measure(label,TRUE)+.9;h<-3.2;xx<-x*width_scale;yy<-fy(y)-y_offset
 left<-xx-w/2;right<-xx+w/2;barleft<-(bar_x-14.5)*width_scale;barright<-(bar_x+14.5)*width_scale
 ci_right<-(bar_x+1.5)*width_scale
 ck(paste('white label inside bar',label),left>=barleft && right<=barright && yy-h/2>=fy(112)-y_offset && yy+h/2<=fy(184)-y_offset)
 ck(paste('white label clears CI',label),left-ci_right>.15)
 grid.rect(x=unit(xx,'mm'),y=unit(yy,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),gp=gpar(fill='white',col=NA))
 txt(label,x,y,9,TRUE,col)
 key<-paste(label,x,y,y_offset)
 bar_label_boxes[[key]]<<-data.frame(label,x_left=left,x_right=right,y_bottom=yy-h/2,y_top=yy+h/2,bar_left=barleft,bar_right=barright,ci_right=ci_right,ci_clearance_mm=left-ci_right,y_offset=y_offset)
}

draw<-function(combo_only=FALSE){
 panel_part<<-'upper';grid.newpage();rc(91.5,110,183,220,'white')
 txt('a',3,214,12,TRUE,just='left');txt('Histological fibrosis and observed F≥2 / F≥3 proportions',13.5,214,10.5,TRUE,just='left')
 txt('n = 935 · Points: observed proportions · Whiskers: Wilson 95% CI of those proportions',13.5,207,col=grey,just='left')
 for(j in 1:5){x<-20+(j-1)*14;rc(x,199,3.2,3.2,stage_col[j]);txt(paste0('F',j-1),x+2.5,199,just='left')}
 for(j in 1:2){x<-113+(j-1)*31;ln(x-3,199,x+3,199,outcome_col[j],1.3);pt(x,199,outcome_col[j],if(j==1)21 else 24,'white');txt(paste0('F≥',j+1),x+5,199,just='left')}
 ax0<-14;ax1<-178;yb<-112;yt<-184;bh<-yt-yb;centres<-c(28.5,61.5,94.5,127.5,160.5);bw<-29
 pyobs<-function(p)yb+bh*p
 for(p in seq(0,1,.25)){y<-pyobs(p);ln(ax0,y,ax1,y,'#E1E5E8',.45);txt(paste0(100*p,'%'),ax0-2,y,just='right')}
 ln(ax0,yb,ax0,yt,ink);ln(ax0,yb,ax1,yb,ink);txt('Observed proportion',2,148,rot=90)
 stage_labels<-list()
 for(k in 1:5){
  z<-st[st$cmrf_count==k,];z<-z[order(z$fibrosis_stage,decreasing=TRUE),];lo<-0
  for(r in 1:5){j<-z$fibrosis_stage[r]+1;p<-z$pct[r]/100;mid<-yb+bh*(lo+p/2)
   if(p>0){rc(centres[k],mid,bw,bh*p,stage_col[j],'white');stage_labels[[length(stage_labels)+1]]<-list(label=fmt(p),x=centres[k]-5.8,y=mid,col=stage_text[j])}
   else txt(paste0('F4 = ',fmt(p),'%'),centres[k],190,col=stage_col[5])
   lo<-lo+p
  }
  txt(k,centres[k],107);txt(paste0('n = ',N[k]),centres[k],102,col=grey)
  # Mathematical identity between stacked boundaries and the observed endpoints.
  for(j in 1:2){zobs<-obs[obs$count==k & obs$outcome==paste0('F≥',j+1),]
   ck(paste('cumulative boundary',k,j),abs(sum(z$pct[z$fibrosis_stage>=j+1])/100-zobs$probability)<1e-12)
  }
 }
 percentage_labels<-list()
 # Draw uncertainty/curves first, then labels. White under-strokes separate marks from bar fills.
 for(j in 1:2){z<-obs[obs$outcome==paste0('F≥',j+1),];z<-z[order(z$count),];col<-outcome_col[j]
  grid.lines(x=unit(centres*width_scale,'mm'),y=unit(fy(pyobs(z$probability))-y_offset,'mm'),gp=gpar(col='white',lwd=3.8))
  grid.lines(x=unit(centres*width_scale,'mm'),y=unit(fy(pyobs(z$probability))-y_offset,'mm'),gp=gpar(col=col,lwd=1.6))
  for(k in 1:5){x<-centres[k];a<-pyobs(z$ci_low[k]);b<-pyobs(z$ci_high[k]);y<-pyobs(z$probability[k])
   for(stroke in c('white',col)){w<-if(stroke=='white')3.2 else 1.0
    ln(x,a,x,b,stroke,w);ln(x-1.5,a,x+1.5,a,stroke,w);ln(x-1.5,b,x+1.5,b,stroke,w)}
   pt(x,y,'white',if(j==1)16 else 17,'white');pt(x,y,col,if(j==1)21 else 24,'white')
   label<-paste0(fmt(z$probability[k]),'%')
   label_y<-pyobs(z$probability[k])-3
   percentage_labels[[paste(j,k)]]<-list(label=label,x=x+8,y=label_y,col=col,bar_x=x)
   geometry[[paste('observed',k,j,width_scale)]]<<-data.frame(variant=width_scale,label=label,width_mm=measure(label,TRUE),available_mm=29*width_scale)

  }
 }
 for(z in stage_labels)txt(z$label,z$x,z$y,col=z$col)
 for(z in percentage_labels)pct_label(z$label,z$x,z$y,z$col,z$bar_x)
 txt('CMRF count',91,96)
 txt('Events / n',3,92,9,TRUE,just='left')
 for(j in 1:2){y<-if(j==1)88.5 else 82.8;txt(paste0('F≥',j+1),13.5,y,9,TRUE,col=outcome_col[j],just='right')
  for(k in 1:5){z<-obs[obs$count==k & obs$outcome==paste0('F≥',j+1),];label<-sprintf('%d/%d',z$events,z$N);txt(label,centres[k],y,col=outcome_col[j])}
 }
 if(combo_only)return(invisible(NULL))
 panel_part<<-'lower'
 # b/c: frozen marginal-standardization estimates, linked at the 5 saved points.
 for(j in 1:2){left<-if(j==1)3 else 99;px0<-left+11;px1<-left+80;yb2<-23;yt2<-57
  col<-outcome_col[j];outcome<-paste0('f_ge',j+1);z<-pr[pr$outcome==outcome,];z<-z[order(z$cmrf_count),];m<-mods[mods$outcome==outcome,]
  txt(letters[j+1],left,74,12,TRUE,just='left');txt(paste0('Standardized F≥',j+1,' probability'),left+10.5,74,10.5,TRUE,just='left')
  txt('Age, sex and centre standardized',left+11,67,col=grey,just='left');txt(if(j==1)'Bootstrap 95% CI · Axis: 20–80%' else 'Bootstrap 95% CI · Axis: 0–60%',left+11,62.5,col=grey,just='left')
  xx<-seq(px0+2,px1-2,length.out=5);lims<-if(j==1)c(.2,.8) else c(0,.6)
  ck(paste('all CI in model axis',j),all(z$ci_low>=lims[1] & z$ci_high<=lims[2]))
  ck(paste('common60pp span',j),abs(diff(lims)-.6)<1e-12)
  py<-function(p)yb2+(yt2-yb2)*(p-lims[1])/diff(lims)
  for(p in seq(lims[1],lims[2],.1)){y<-py(p);ln(px0,y,px1,y,'#E1E5E8',.45);txt(paste0(round(100*p),'%'),px0-2,y,just='right')}
  poly(c(xx,rev(xx)),py(c(z$ci_low,rev(z$ci_high))),adjustcolor(col,alpha.f=.16))
  grid.lines(x=unit(xx*width_scale,'mm'),y=unit(fy(py(z$probability))-y_offset,'mm'),gp=gpar(col=col,lwd=1.5))
  for(k in 1:5){pt(xx[k],py(z$probability[k]),col,if(j==1)21 else 24,'white');txt(paste0(fmt(z$probability[k]),'%'),xx[k],py(z$ci_high[k])+3,col=outcome_text[j],just=if(k==1)'left' else if(k==5)'right' else 'centre')
   ln(xx[k],yb2,xx[k],yb2-1.1,ink);txt(k,xx[k],18.5)}
  ln(px0,yb2,px0,yt2,ink);ln(px0,yb2,px1,yb2,ink)
  txt('Standardized probability',left-1,38.5,rot=90);txt('CMRF count',mean(c(px0,px1)),13)

 }
 drawd()
}

# d contains the main complete-component estimates only; robustness remains S1.
drawd<-function(newpage=FALSE){
 panel_part<<-'forest'
 if(newpage){grid.newpage();rc(91.5,21,183,42,'white')}
 txt('d',3,36.5,12,TRUE,just='left')
 txt('Adjusted OR per additional CMRF',13.5,36.5,10.5,TRUE,just='left')
 txt('Main analysis · n = 935 · Adjusted for age, sex and centre',13.5,30.5,col=grey,just='left')
 txt('Outcome',3,24.5,9,TRUE,just='left');txt('Events / n',38,24.5,9,TRUE)
 txt('OR (log scale)',82,24.5,9,TRUE);txt('OR (95% CI)',139,24.5,9,TRUE);txt('P',170.5,24.5,9,TRUE)
 x0<-55;x1<-109;lim<-c(.8,2);px<-function(v)x0+(x1-x0)*log(v/lim[1])/log(lim[2]/lim[1])
 grid.lines(x=unit(rep(px(1)*width_scale,2),'mm'),y=unit(c(8.5,22),'mm'),gp=gpar(col=grey,lwd=.75,lty=2))
 for(j in 1:2){
  z<-mods[mods$outcome==paste0('f_ge',j+1),];y<-c(18,10.5)[j];col<-outcome_col[j]
  ck(paste('d complete-component denominator',j),nrow(z)==1 && z$N==935)
  ck(paste('d full confidence interval',j),z$ci_low>lim[1] && z$ci_high<lim[2])
  txt(paste0('F≥',j+1),3,y,9,TRUE,col,just='left');txt(sprintf('%d/%d',z$events,z$N),38,y)
  a<-px(z$ci_low);b<-px(z$ci_high);x<-px(z$OR)
  ln(a,y,b,y,col,1);ln(a,y-.7,a,y+.7,col,.7);ln(b,y-.7,b,y+.7,col,.7);pt(x,y,col,if(j==1)21 else 24,col)
  effect<-sprintf('%.2f (%.2f–%.2f)',z$OR,z$ci_low,z$ci_high)
  pvalue<-if(z$p_value<.001)sprintf('%.2e',z$p_value) else sprintf('%.5f',z$p_value)
  txt(effect,139,y);txt(pvalue,170.5,y)
  geometry[[paste('d OR',j)]]<<-data.frame(variant=width_scale,label=effect,width_mm=measure(effect),available_mm=33*width_scale)
  geometry[[paste('d P',j)]]<<-data.frame(variant=width_scale,label=pvalue,width_mm=measure(pvalue),available_mm=17*width_scale)
 }
 ln(x0,6.5,x1,6.5,ink)
 for(v in c(.8,1,1.25,1.5,2)){x<-px(v);ln(x,6.5,x,5.5,ink);txt(format(v,trim=TRUE),x,2.5)}
}
robust<-rd('robustness_models.csv');mi<-rd('main_CC_MI_models.csv')
ck('30 supplementary robustness estimates',nrow(robust)==30 && all(is.finite(robust$OR)) && all(robust$ci_low<robust$OR & robust$ci_high>robust$OR))
for(j in 1:2){z<-robust[robust$outcome==paste0('f_ge',j+1) & robust$analysis=='Adult raw-variable MI',];m<-mi[mi$outcome==paste0('f_ge',j+1) & mi$m==50 & !is.na(mi$m),]
 ck(paste('MI retained in supplementary',j),nrow(z)==1 && nrow(m)==1 && max(abs(unlist(z[c('OR','ci_low','ci_high','p_value')])-unlist(m[c('OR','ci_low','ci_high','p_value')])))<1e-12)
}
panel_part<-'lower';ck('b c same 34mm plot height',fh(57-23)==34)
if(file.exists(file.path(out,'figures/Figure2_fibrosis.pdf')))stop('Refusing to overwrite existing figure')
width_scale<-.85;y_offset<-0
probe<-tempfile(fileext='.pdf');cairo_pdf(probe,width=155.55/25.4,height=246/25.4,family='Arial');draw();dev.off()
bb0<-do.call(rbind,text_boxes)
if(any(bb0$x_low<0 | bb0$x_high>bb0$canvas_width))stop(paste('Text outside narrow width:',paste(bb0$label[bb0$x_low<0 | bb0$x_high>bb0$canvas_width],collapse='; ')))
if(any(vapply(geometry,function(g)g$width_mm>=g$available_mm,logical(1))))stop('Narrow table/point labels exceed available space')
for(scale in .85){
 y_offset<-0; width_scale<-scale;W<-183*scale;suffix<-'';stem<-file.path(out,'figures',paste0('Figure2_fibrosis',suffix))
 cairo_pdf(paste0(stem,'.pdf'),width=W/25.4,height=H/25.4,family='Arial');draw();dev.off()
 svglite(paste0(stem,'.svg'),width=W/25.4,height=H/25.4,bg='white');draw();dev.off()
 agg_png(paste0(stem,'.png'),width=W,height=H,units='mm',res=300,background='white');draw();dev.off()
 agg_tiff(paste0(stem,'.tiff'),width=W,height=H,units='mm',res=300,compression='lzw',background='white');draw();dev.off()
}
# Standalone a: same native font sizes and marks as the complete figure.
y_offset<-112;width_scale<-.85;stem<-file.path(out,'figures/Figure2a_combined')
cairo_pdf(paste0(stem,'.pdf'),width=155.55/25.4,height=132/25.4,family='Arial');draw(TRUE);dev.off()
svglite(paste0(stem,'.svg'),width=155.55/25.4,height=132/25.4,bg='white');draw(TRUE);dev.off()
agg_png(paste0(stem,'.png'),width=155.55,height=132,units='mm',res=300,background='white');draw(TRUE);dev.off()
agg_tiff(paste0(stem,'.tiff'),width=155.55,height=132,units='mm',res=300,compression='lzw',background='white');draw(TRUE);dev.off()
write.csv(do.call(rbind,bar_label_boxes),file.path(out,'qa/white_label_geometry.csv'),row.names=FALSE)
y_offset<-0
bb<-do.call(rbind,text_boxes);write.csv(bb,file.path(out,'qa/text_bounds.csv'),row.names=FALSE)
ck('all text within width',all(bb$x_low>=0 & bb$x_high<=bb$canvas_width))
ck('body minimum9pt',min(bb$size_pt)>=9)
g<-do.call(rbind,geometry);write.csv(g,file.path(out,'qa/label_widths.csv'),row.names=FALSE)
ck('label width fits',all(g$width_mm<g$available_mm))
write.csv(data.frame(check=names(checks),pass=unlist(checks)),file.path(out,'qa/numerical_checks.csv'),row.names=FALSE)
write.csv(data.frame(stage=paste0('F',0:4),hex=stage_col),file.path(out,'source_data/palette_mapping.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'qa/sessionInfo.txt'))
cat(length(checks),'checks passed; four-panel Figure2 and standalone a exported.\n')

library(magick)
im<-image_read(file.path(out,'figures/Figure2_fibrosis.png'))
image_write(image_resize(im,'1600x'),file.path(out,'preview/Figure2_preview.png'))
image_write(image_convert(im,colorspace='gray'),file.path(out,'qa/Figure2_grayscale.png'))

oldim<-image_read(file.path(root,'paper1_cmrf_pathology/results/versions/20261024_093625/figures/Figure2_fibrosis.png'))
info<-image_info(im);cut<-ceiling((246-111)*300/25.4);crop<-sprintf('%dx%d+0+%d',info$width,info$height-cut,cut)
comparison<-image_compare(image_crop(oldim,crop),image_crop(im,crop),metric='AE');ae<-as.numeric(attr(comparison,'distortion'))
stopifnot(ae==0);write.csv(data.frame(region='b c d unchanged',AE=ae),file.path(out,'qa/pixel_checks.csv'),row.names=FALSE)
image_write(image_resize(image_read(file.path(out,'figures/Figure2a_combined.png')),'1600x'),file.path(out,'preview/Figure2a_preview.png'))
