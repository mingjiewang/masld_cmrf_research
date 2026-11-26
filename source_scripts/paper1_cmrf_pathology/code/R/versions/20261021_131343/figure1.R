# Figure 1: cohort flow and centre profiles. Revision author: mjwang.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Pass project root')
root <- normalizePath(args[1]); version <- '20261021_131343'
paper <- file.path(root,'paper1_cmrf_pathology')
src <- file.path(paper,'results/versions/20260809_152651')
out <- file.path(paper,'results/versions',version)
library(grid); library(svglite); library(ragg); library(ggsci)
previous <- file.path(paper,'results/versions/20261021_125942')
frozen_files <- c('centre_counts.csv','fibrosis_stage_by_centre.csv','clinical_profiles.csv',
                  'clinical_overall_tests.csv','analysis_population_flow.csv','palette_mapping.csv')
for(name in frozen_files) {
 f <- file.path(previous,'source_data',name); dest <- file.path(out,'source_data',name)
 if(!file.exists(dest)) file.copy(f,dest)
 stopifnot(identical(readBin(f,'raw',n=file.info(f)$size),readBin(dest,'raw',n=file.info(dest)$size)))
}

read_summary <- function(name) read.csv(file.path(out,'source_data',name))
counts <- read_summary('centre_counts.csv')
st <- read_summary('fibrosis_stage_by_centre.csv')
descr <- read_summary('clinical_profiles.csv')
stats <- read_summary('clinical_overall_tests.csv')
centres <- c('HZ','JL','RJ'); names_c <- c(HZ='Hangzhou',JL='Jilin',RJ='Ruijin')
labels <- unique(descr$label)
width_scale <- 1.30; figure_width <- 183*width_scale; figure_height <- 205
palette <- ggsci::pal_nejm('default')(8)
ink <- '#202832'; navy <- palette[2]; magenta <- palette[1]; teal <- palette[4]
centre_col <- c(HZ=navy,JL=teal,RJ=magenta)
stage_col <- palette[c(6,2,7,3,1)]
stage_text <- c(ink,'white',ink,ink,'white')
write.csv(data.frame(stage=paste0('F',0:4),group=c('F0–F1','F0–F1',rep('F≥2',3)),
  palette='ggsci::pal_nejm(default)',index=c(6,2,7,3,1),hex=stage_col),file.path(out,'source_data/palette_mapping.csv'),row.names=FALSE)
txt <- function(label,x,y,size=9,bold=FALSE,col=ink,just='centre',lh=1.1) {
 grid.text(label,x=unit(x*width_scale,'mm'),y=unit(y,'mm'),just=just,
           gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=lh))
}
line <- function(x1,y1,x2,y2,col='#8C969F',lwd=.65,arrow=FALSE) {
 grid.lines(x=unit(c(x1,x2)*width_scale,'mm'),y=unit(c(y1,y2),'mm'),
            gp=gpar(col=col,lwd=lwd),arrow=if(arrow)grid::arrow(length=unit(1.5,'mm'),type='closed') else NULL)
}
rect <- function(x,y,w,h,fill='white',border=NA,lwd=.65) {
 grid.rect(x=unit(x*width_scale,'mm'),y=unit(y,'mm'),width=unit(w*width_scale,'mm'),height=unit(h,'mm'),
           gp=gpar(fill=fill,col=border,lwd=lwd))
}
box <- function(label,x,y,w,h,em=FALSE) {
 rect(x,y,w,h,if(em)'#EAF0F8' else 'white',if(em)navy else '#99A9B7',if(em).95 else .65)
 txt(label,x,y,size=9,bold=em,lh=1.12)
}
# Geometry uses physical millimetres so the pies stay circular in both layouts.
roundrect <- function(x,y,w,h,fill,border,lwd=.6) {
 grid.roundrect(x=unit(x*width_scale,'mm'),y=unit(y,'mm'),width=unit(w*width_scale,'mm'),
   height=unit(h,'mm'),r=unit(1.3,'mm'),gp=gpar(fill=fill,col=border,lwd=lwd))
}
stage_card <- function(step,heading,count,sub,y,kind='standard') {
 dark <- kind=='main'; accent <- if(kind=='eligible')teal else navy
 fill <- if(dark)navy else if(kind=='eligible')'#ECF5F0' else '#F5F8FB'
 roundrect(22.5,y,40,if(dark)28 else 26,fill,if(dark)navy else '#CCDAE4',.65)
 if(!dark)rect(3.0,y,1,22,accent)
 grid.circle(x=unit(7*width_scale,'mm'),y=unit(y+8,'mm'),r=unit(2.1,'mm'),
   gp=gpar(fill=if(dark)'white' else accent,col=NA))
 txt(sprintf('%02d',step),7,y+8,9,TRUE,col=if(dark)navy else 'white')
 txt(heading,11.7,y+8,9,TRUE,col=if(dark)'white' else ink,just='left')
 txt(paste0('n = ',count),7,y-1.0,17,TRUE,col=if(dark)'white' else accent,just='left')
 txt(sub,7,y-9,9,col=if(dark)'#E0EDF7' else '#536675',just='left')
}
branch_card <- function(heading,reason,count,y,retained=FALSE) {
 accent <- if(retained)teal else palette[3]
 roundrect(57.5,y,23,if(retained)28 else 26,if(retained)'#EFF6F1' else '#FFF6EE',
           if(retained)'#B7D6C3' else '#ECD3BE',.6)
 rect(46.5,y,1,if(retained)24 else 22,accent)
 txt(heading,57.5,y+8,9,TRUE,col=if(retained)teal else '#9A5916')
 txt(reason,57.5,y+1.4,9,col=ink,lh=1.0)
 txt(paste0('n = ',count),57.5,y-8.1,11,TRUE,col=if(retained)teal else ink)
}
fade <- function(col,f=.18) {z<-col2rgb(col)/255;rgb(1-f+f*z[1],1-f+f*z[2],1-f+f*z[3])}
pie_rows <- list()
centre_pie <- function(focus,cx,cy,radius=5.4) {
 share <- counts$analysed/sum(counts$analysed)
 bounds <- pi/2-c(0,cumsum(share))*2*pi
 for(j in 1:3) {
   focal <- centres[j]==focus; theta<-seq(bounds[j],bounds[j+1],length.out=181)
   mid<-mean(bounds[j:(j+1)]);offset<-if(focal).85 else 0
   x0<-cx*width_scale+offset*cos(mid);y0<-cy+offset*sin(mid)
   xp<-c(x0,x0+radius*cos(theta),x0);yp<-c(y0,y0+radius*sin(theta),y0)
   fill<-if(focal)centre_col[centres[j]] else fade(centre_col[centres[j]])
   grid.polygon(x=unit(xp,'mm'),y=unit(yp,'mm'),gp=gpar(fill=fill,col=if(focal)centre_col[focus] else 'white',lwd=if(focal).8 else .5))
   area<-abs(sum(xp[-length(xp)]*yp[-1]-xp[-1]*yp[-length(yp)]))/2
   pie_rows[[paste(focus,centres[j])]] <<- data.frame(pie_focus=focus,sector_centre=centres[j],
     n=counts$analysed[j],N=sum(counts$analysed),percentage=share[j]*100,
     angle_degrees=share[j]*360,radius_mm=radius,exploded_mm=offset,highlight=focal,
     polygon_area_mm2=area,theoretical_sector_area_mm2=pi*radius^2*share[j])
 }
 value<-counts$analysed[counts$centre==focus]/sum(counts$analysed)*100
 txt(sprintf('%.1f%%',value),cx,cy-8.5,9,TRUE,col=centre_col[focus])
}
callout_rows <- list()
draw <- function() {
 grid.newpage(); pushViewport(viewport(width=unit(figure_width,'mm'),height=unit(figure_height,'mm')))
 rect(91.5,102.5,183,205)
 # a: numbered screening stages with visually distinct analysis populations.
 txt('a',3,199,12,TRUE,just='left');txt('Patient selection',9,199,10.5,TRUE,just='left')
 ys<-c(181,147,113,79,41)
 stage_card(1,'Source records','1,100','Starting dataset',ys[1])
 stage_card(2,'Adult cohort','1,098','Age ≥18 years',ys[2])
 stage_card(3,'Recorded histology','1,087','S and F recorded',ys[3])
 stage_card(4,'Eligible cohort','1,085','≥1 definite CMRF · MI',ys[4],'eligible')
 stage_card(5,'Main analysis','935','Complete-component',ys[5],'main')
 yy<-c(164,130,96,60.5)
 for(k in 1:4) {
   top<-ys[k]-13;bottom<-ys[k+1]+if(k==4)14 else 13
   line(22.5,top,22.5,bottom,navy,.9,TRUE)
   grid.circle(x=unit(22.5*width_scale,'mm'),y=unit(yy[k],'mm'),r=unit(.65,'mm'),gp=gpar(fill=navy,col=NA))
   line(22.5,yy[k],45.5,yy[k],'#9AA8B3',.7,TRUE)
 }
 branch_card('Excluded','Age <18 years','2',yy[1])
 branch_card('Excluded','No S value*','11',yy[2])
 branch_card('Excluded','No definite\nCMRF','2',yy[3])
 branch_card('MI retained','Indeterminate\nCMRF domain','150',yy[4],TRUE)
 roundrect(22.5,19.5,40,9,'#EEF4F8',NA)
 txt(sprintf('%.1f%% of eligible adults',935/1085*100),22.5,19.5,9,TRUE,col=navy)
 txt('*S absent in 11; F recorded in 10.\nS, steatosis; F, fibrosis.\nMI, multiple imputation.',3,7.5,9,just='left')
 # b: actual counts and labelled fibrosis distribution.
 rx <- 74;txt('b',rx,199,12,TRUE,just='left')
 txt('Centre contribution and fibrosis',rx+6,199,10.5,TRUE,just='left')
 txt('Complete-component cohort, n = 935',rx+6,191,9,just='left',col='#586572')
 txt('F0–F1',104.4,186.5,9,TRUE,col=navy)
 txt('F≥2',141.4,186.5,9,TRUE,col=magenta)
 for(j in 1:5) {x <- 97+(j-1)*14.8;rect(x,182,3.2,3.2,stage_col[j]);txt(paste0('F',j-1),x+2.5,182,9,just='left')}
 txt('Centre',75,174,9,TRUE,just='left')
 txt('n / eligible',75,169.5,9,just='left',col='#586572')
 txt('Share',104,174,9,TRUE)
 txt('Stage (%)',139.5,174,9,TRUE)
 txt('F≥2',174,174,9,TRUE,col=magenta)
 bx <- 115; bw <- 49
 for(i in 1:3) {
   y <- c(163,145,127)[i];cc <- centres[i];a <- st[st$centre==cc,]; n <- counts[counts$centre==cc,]
   txt(names_c[cc],75,y+2,9,TRUE,just='left',col=centre_col[cc]);txt(sprintf('%d / %d',n$analysed,n$eligible),75,y-3.8,9,just='left',col='#586572')
   centre_pie(cc,104,y)
   pos <- bx
   for(j in 1:5) {
     w <- bw*a$pct[j]/100;rect(pos+w/2,y,w,8,stage_col[j],border='white',lwd=.35)
     if(a$pct[j]>=11) txt(sprintf('%.1f',a$pct[j]),pos+w/2,y,9,col=stage_text[j])
     else {
       anchor <- pos+w/2
       target <- if(cc=='JL') (if(width_scale==1) c(`2`=148,`3`=158,`4`=168) else
         c(`2`=150,`3`=157.8,`4`=165.6))[as.character(a$stage[j])] else
         if(a$stage[j]==0)116 else 164.5
       line(anchor,y+4.1,anchor,y+5.2,'#707B85',.6)
       line(anchor,y+5.2,target,y+7.0,'#707B85',.6)
       txt(sprintf('%.1f%%',a$pct[j]),target,y+8.7,9)
       callout_rows[[paste(cc,a$stage[j])]] <<- data.frame(centre=cc,stage=a$stage[j],
         n=a$n[j],N=a$N[j],pct=a$pct[j],label=sprintf('%.1f%%',a$pct[j]),
         anchor_x=anchor,label_x=unname(target),bar_y=y,label_y=y+8.7)
     }
     pos <- pos+w
   }
   cut <- bx+bw*sum(a$pct[1:2])/100
   line(cut,y-4.4,cut,y+4.4,ink,.85)
   txt(sprintf('%.1f%%',n$f_ge2_pct),174,y,9,TRUE,col=magenta)
 }
 line(bx,119,bx+bw,119,'#7D8893',.55)
 for(t in c(0,50,100)) {x<-bx+bw*t/100;line(x,119,x,117.8);txt(paste0(t,'%'),x,114.8,9)}
 # c: graphical table; numeric values remain primary and bars show actual percentages.
 txt('c',rx,104,12,TRUE,just='left')
 txt('Clinical profiles by centre',rx+6,104,10.5,TRUE,just='left')
 cx <- c(115,138,161)
 txt('Characteristic',75,95,9,TRUE,just='left')
 for(i in 1:3) {txt(centres[i],cx[i],96.5,9,TRUE,col=centre_col[i]);txt(paste0('n = ',counts$analysed[i]),cx[i],92.3,9)}
 txt('P*',177,95,9,TRUE);line(75,88,181,88,'#BCC7D0')
 rowy <- 82-(0:9)*7.7
 for(i in 1:10) {
   y <- rowy[i]; lab <- labels[i]
   if(lab=='Glycaemic abnormality') lab <- 'Glycaemic\nabnormality'
   if(lab=='High triglycerides') lab <- 'High\ntriglycerides'
   txt(lab,75,y,9,just='left',lh=1.0)
   for(k in 1:3) {
     a <- descr[descr$row==i & descr$centre==centres[k],]
     txt(a$display,cx[k],y+if(i>2)1.2 else 0,9)
     if(i>2) {
       rect(cx[k],y-2.5,18,1.0,'#E9EDF2')
       w <- 18*a$percentage/100
       rect(cx[k]-9+w/2,y-2.5,w,1.0,centre_col[k])
     }
   }
   txt(stats$p_display[i],177,y,9)
   line(75,y-4.0,181,y-4.0,'#E5E9EF',.4)
 }
 txt('Mean ± SD or n (%); short bars: 0–100%.',75,6.7,9,just='left',col='#586572')
 txt('*Overall P, Holm-adjusted across 10 comparisons.',75,2.2,9,just='left',col='#586572')
 popViewport()
}
variants <- data.frame(name=c('wide_master','layout_check_183mm'),scale=c(1.3,1.0),
  stem=c('Figure1_cohort_and_centre_profiles','Figure1_cohort_and_centre_profiles_183mm'))
for(vi in seq_len(nrow(variants))) {
 callout_rows <- list()
 width_scale <- variants$scale[vi];figure_width <- 183*width_scale
 base <- file.path(out,'figures',variants$stem[vi])
 cairo_pdf(paste0(base,'.pdf'),width=figure_width/25.4,height=figure_height/25.4,family='Arial');draw();dev.off()
 svglite(paste0(base,'.svg'),width=figure_width/25.4,height=figure_height/25.4);draw();dev.off()
 agg_png(paste0(base,'.png'),width=figure_width,height=figure_height,units='mm',res=300,background='white');draw();dev.off()
 agg_tiff(paste0(base,'.tiff'),width=figure_width,height=figure_height,units='mm',res=600,compression='lzw',background='white');draw();dev.off()
 callout_name <- if(vi==1)'small_segment_callouts.csv' else 'small_segment_callouts_183mm.csv'
 write.csv(do.call(rbind,callout_rows),file.path(out,'source_data',callout_name),row.names=FALSE)
}
variants$width_mm <- 183*variants$scale;variants$height_mm <- figure_height;variants$body_text_pt <- 9
write.csv(variants,file.path(out,'source_data/export_variants.csv'),row.names=FALSE)
labels_all <- st;labels_all$display_percent <- sprintf('%.1f',labels_all$pct)
labels_all$label_position <- ifelse(labels_all$pct>=11,'inside','external elbow leader')
write.csv(labels_all,file.path(out,'source_data/fibrosis_display_labels.csv'),row.names=FALSE)
write.csv(do.call(rbind,pie_rows),file.path(out,'source_data/centre_share_pies.csv'),row.names=FALSE)
write.csv(data.frame(eligible_n=1085,complete_component_n=935,complete_percentage=935/1085*100),file.path(out,'source_data/analysis_completeness.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,'sessionInfo.txt'))
cat('Figure 1: centre-share pies and numbered screening stages; frozen clinical values retained.\n')
print(stats[,c('variable','method','p_raw','p_holm_10')],row.names=FALSE)
