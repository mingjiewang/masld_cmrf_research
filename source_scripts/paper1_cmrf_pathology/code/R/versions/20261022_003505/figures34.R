# Figures3/4 redesign from frozen results. Revision author: mjwang.
args<-commandArgs(trailingOnly=TRUE);if(length(args)!=1)stop('Pass project root')
root<-normalizePath(args[1]);version<-'20261022_003505'
out<-file.path(root,'paper1_cmrf_pathology/results/versions',version)
library(grid);library(svglite);library(ragg);library(ggsci)
rd<-function(n)read.csv(file.path(out,'source_data',n),check.names=FALSE)
s<-rd('steatosis_distribution.csv');so<-rd('observed_sge2_by_cmrf.csv');sm<-rd('steatosis_secondary_models.csv')
comp<-rd('Figure4a_component_models.csv');gly<-rd('Figure4b_glycaemia_decomposition.csv')
checks<-list();ck<-function(n,p){checks[[n]]<<-isTRUE(p);if(!isTRUE(p))stop(n)}
N<-sapply(1:5,function(k)sum(s$n[s$CMRF==k]));ck('N935groups',all(N==c(40,188,313,292,102)))
ck('Sgrades20cells',nrow(s)==20);ck('S0is19',sum(s$n[s$S==0])==19)
for(i in 1:5){z<-so[so$cmrf_count==i,];p<-z$events/z$N;q<-qnorm(.975);den<-1+q^2/z$N;mid<-(p+q^2/(2*z$N))/den;hw<-q*sqrt(p*(1-p)/z$N+q^2/(4*z$N^2))/den
 ck(paste('Sge2 events',i),z$events==sum(s$n[s$CMRF==i & s$S>=2]) && z$N==N[i])
 ck(paste('Wilson',i),max(abs(c(p,mid-hw,mid+hw)-unlist(z[c('probability','ci_low','ci_high')])))<1e-12)
 ck(paste('cumulative S boundary',i),abs(sum(s$n[s$CMRF==i & s$S>=2])/N[i]-z$probability)<1e-12)
}
ck('Smodel endpoints',identical(as.character(sm$outcome),c('s_ge1','s_ge2','s_ge3')))
for(i in 1:3){z<-sm[i,];ck(paste('Smodel events',i),z$events==sum(s$n[s$S>=i]));ck(paste('Smodel ORCI',i),max(abs(c(exp(z$estimate),exp(z$estimate+c(-1,1)*1.96*z$std_error))-c(z$OR,z$ci_low,z$ci_high)))<1e-12)}
ck('BH3',max(abs(p.adjust(sm$p_value,'BH')-sm$p_fdr_3))<1e-12)
for(type in c('component','glycaemia')){
 d<-if(type=='component')comp else gly;adj<-if(type=='component')'p_fdr_all_10' else 'p_holm_all_4'
 for(set in c('Complete-component','Multiple imputation')){
  z<-d[d$dataset==set,];f<-if(type=='component')paste0('component_models_',if(set=='Complete-component')'complete_component.csv' else 'multiple_imputation.csv') else paste0('glycaemia_decomposition_',if(set=='Complete-component')'complete_component.csv' else 'multiple_imputation.csv')
  raw<-rd(f);key<-function(x)paste(x$outcome,x$term,x$level,sep='|');a<-raw[match(key(z),key(raw)),]
  ck(paste(type,set,'source rows'),all(!is.na(a$OR)) && nrow(a)==nrow(z))
  for(field in c('N','events','estimate','std_error','OR','ci_low','ci_high','p_value',adj))ck(paste(type,set,field),max(abs(z[[field]]-a[[field]]))<1e-12)
  ck(paste(type,set,'adjustedP'),max(abs(p.adjust(z$p_value,if(type=='component')'BH' else 'holm')-z[[adj]]))<1e-12)
  ck(paste(type,set,'flags'),all(z$supported==(z[[adj]]<.05)))
  for(i in 1:nrow(z)){u<-z[i,];crit<-if(set=='Complete-component')1.96 else qt(.975,u$df)
   ck(paste(type,set,'CI row',i),max(abs(exp(u$estimate+c(-1,1)*crit*u$std_error)-c(u$ci_low,u$ci_high)))<1e-10)
   ck(paste(type,set,'OR row',i),abs(exp(u$estimate)-u$OR)<1e-12)
  }
 }
}
ck('20components',nrow(comp)==20);ck('8glycaemia',nrow(gly)==8)
ck('componentCIwithin',all(comp$ci_low>=.06 & comp$ci_high<=5.2));ck('glycaemiaCIwithin',all(gly$ci_low>=.85 & gly$ci_high<=7.1))
pal<-pal_nejm('default')(8);ink<-'#202832';grey<-'#586572';blue<-pal[2];brick<-pal[1]
scol<-pal[c(6,2,7,3)];stext<-c(ink,'white',ink,ink)
W<-183;geometry<-list();drawn<-list();figure_name<-''
txt<-function(label,x,y,size=9,bold=FALSE,col=ink,just='centre',rot=0){grid.text(label,x=unit(x,'mm'),y=unit(y,'mm'),just=just,rot=rot,gp=gpar(fontfamily='Arial',fontsize=size,fontface=if(bold)'bold' else 'plain',col=col,lineheight=1.02))}
ln<-function(x1,y1,x2,y2,col='#A5AFB7',lwd=.6,lty=1){grid.lines(x=unit(c(x1,x2),'mm'),y=unit(c(y1,y2),'mm'),gp=gpar(col=col,lwd=lwd,lty=lty))}
rc<-function(x,y,w,h,fill='white',border=NA){grid.rect(x=unit(x,'mm'),y=unit(y,'mm'),width=unit(w,'mm'),height=unit(h,'mm'),gp=gpar(fill=fill,col=border,lwd=.45))}
pt<-function(x,y,col=ink,pch=21,fill='white'){grid.points(x=unit(x,'mm'),y=unit(y,'mm'),pch=pch,size=unit(2.2,'mm'),gp=gpar(col=col,fill=fill,lwd=.9))}
measure<-function(z){pushViewport(viewport(gp=gpar(fontfamily='Arial',fontsize=9)));w<-convertWidth(stringWidth(z),'mm',TRUE);popViewport();w}
fmt<-function(x)sprintf('%.1f',100*x)
pfmt<-function(p)if(p<.001)sprintf('%.2e',p) else formatC(p,format='fg',digits=3,flag='#')
ci_fmt<-function(x)if(abs(x-1)<.005 && x!=1)sprintf('%.3f',x) else sprintf('%.2f',x)
efmt<-function(z)sprintf('%.2f (%s–%s)',z$OR,ci_fmt(z$ci_low),ci_fmt(z$ci_high))
label_fit<-function(key,label,x,y,available){txt(label,x,y);geometry[[paste(figure_name,key)]]<<-data.frame(figure=figure_name,key=key,label=label,width_mm=measure(label),available_mm=available)}
logx<-function(p,x0,x1,lim)x0+(x1-x0)*(log(p)-log(lim[1]))/diff(log(lim))
forest<-function(z,y,x0,x1,lim,col,pch=21,supported=TRUE,key=''){
 a<-logx(z$ci_low,x0,x1,lim);b<-logx(z$ci_high,x0,x1,lim);x<-logx(z$OR,x0,x1,lim)
 ln(a,y,b,y,col,1);ln(a,y-.7,a,y+.7,col,.7);ln(b,y-.7,b,y+.7,col,.7);pt(x,y,col,pch,if(supported)col else 'white')
 drawn[[paste(figure_name,key)]]<<-data.frame(figure=figure_name,key=key,OR=z$OR,ci_low=z$ci_low,ci_high=z$ci_high,plot_low=lim[1],plot_high=lim[2],x=x,x_low=a,x_high=b,y=y)
}
draw3<-function(){
 figure_name<<-'Figure3';grid.newpage();rc(W/2,108,W,216)
 txt('a',3,210,12,TRUE,just='left');txt('Histological steatosis and observed S≥2 proportion',13.5,210,10.5,TRUE,just='left')
 txt('n = 935 · Points: observed S≥2 · Whiskers: Wilson 95% CI of that proportion',13.5,203,col=grey,just='left')
 for(j in 1:4){x<-20+(j-1)*17;rc(x,195,3.2,3.2,scol[j]);txt(paste0('S',j-1),x+2.5,195,just='left')}
 ln(119,195,125,195,blue,1.5);pt(122,195,blue);txt('S≥2',128,195,just='left')
 xx<-c(25,58,91,124,157);bw<-22;yb<-100;yt<-174;bh<-yt-yb;py<-function(p)yb+bh*p
 txt('S0 (%)',3,182,col=grey,just='left')
 for(p in seq(0,1,.25)){y<-py(p);ln(14,y,178,y,'#E1E5E8',.45);txt(paste0(100*p,'%'),12,y,just='right')}
 ln(14,yb,14,yt,ink);ln(14,yb,178,yb,ink);txt('Observed proportion',3,137,rot=90)
 labs<-list()
 for(k in 1:5){z<-s[s$CMRF==k,];z<-z[order(z$S,decreasing=TRUE),];lo<-0
  for(r in 1:4){j<-z$S[r]+1;p<-z$n[r]/N[k];mid<-py(lo+p/2)
   if(p>0)rc(xx[k],mid,bw,bh*p,scol[j],'white')
   if(j>1)labs[[length(labs)+1]]<-list(label=fmt(p),x=xx[k]-5.8,y=mid,col=stext[j])
   if(j==1)txt(fmt(p),xx[k],182,col=grey)
   lo<-lo+p
  }
  txt(k,xx[k],95);txt(paste0('n = ',N[k]),xx[k],90,col=grey)
 }
 grid.lines(x=unit(xx,'mm'),y=unit(py(so$probability),'mm'),gp=gpar(col='white',lwd=3.8))
 grid.lines(x=unit(xx,'mm'),y=unit(py(so$probability),'mm'),gp=gpar(col=blue,lwd=1.6))
 for(k in 1:5){z<-so[so$cmrf_count==k,];a<-py(z$ci_low);b<-py(z$ci_high);y<-py(z$probability)
  for(col in c('white',blue)){w<-if(col=='white')3.2 else 1;ln(xx[k],a,xx[k],b,col,w);ln(xx[k]-1.5,a,xx[k]+1.5,a,col,w);ln(xx[k]-1.5,b,xx[k]+1.5,b,col,w)}
  pt(xx[k],y,blue);txt(paste0(fmt(z$probability),'%'),xx[k]+16.5,y-3,9,TRUE,col=blue)
 }
 for(z in labs)txt(z$label,z$x,z$y,col=z$col)
 txt('CMRF count',91,84);txt('S≥2 n/N',3,78,9,TRUE,col=blue,just='left')
 for(k in 1:5){z<-so[so$cmrf_count==k,];txt(sprintf('%d/%d',z$events,z$N),xx[k],78,col=blue)}
 txt('b',3,65,12,TRUE,just='left');txt('Adjusted associations across steatosis thresholds',13.5,65,10.5,TRUE,just='left')
 txt('Age, sex and centre adjusted · Three formal secondary thresholds',13.5,58,col=grey,just='left')
 txt('Outcome',3,51,9,TRUE,just='left');txt('Events / n',40,51,9,TRUE);txt('OR (95% CI)',141,51,9,TRUE);txt('P (BH)',173,51,9,TRUE)
 lim<-c(.8,6);x0<-60;x1<-119;ln(logx(1,x0,x1,lim),17,logx(1,x0,x1,lim),47,grey,.75,2)
 for(i in 1:3){z<-sm[i,];y<-c(43,32,21)[i];txt(c('S≥1','S≥2','S3')[i],3,y,9,TRUE,just='left');txt(sprintf('%d/%d',z$events,z$N),40,y)
  forest(z,y,x0,x1,lim,blue,21,TRUE,paste0('S',i));label_fit(paste('OR',i),efmt(z),141,y,30);label_fit(paste('P',i),pfmt(z$p_fdr_3),173,y,18)
 }
 ln(x0,13,x1,13,ink);for(v in c(.8,1,2,4,6)){x<-logx(v,x0,x1,lim);ln(x,13,x,11.8,ink);txt(v,x,9)}
 txt('Adjusted OR per additional CMRF (log scale)',90,3)
}
draw4<-function(){
 figure_name<<-'Figure4';grid.newpage();rc(W/2,123,W,246)
 txt('a',3,240,12,TRUE,just='left');txt('CMRF component associations with fibrosis',11,240,10.5,TRUE,just='left')
 txt('Mutually adjusted for all five domains, age, sex and centre',11,233,col=grey,just='left')
 pt(14,226,ink,21);txt('CC: complete-component (n = 935)',18,226,just='left')
 pt(100,226,ink,24);txt('MI: imputed (n = 1,085; m = 50)',104,226,just='left')
 plotx<-list(c(40,66),c(114,140));ox<-c(82,156);pp<-c(103,176);col<-c(blue,brick)
 txt('Component',3,210,9,TRUE,just='left');txt('Set',33,210,9,TRUE)
 for(j in 1:2){txt(paste0('F≥',j+1),mean(c(plotx[[j]],ox[j],pp[j])),218,10.5,TRUE,col=col[j]);txt('OR (95% CI)',ox[j],210,9,TRUE);txt('P (BH)',pp[j],210,9,TRUE)}
 labels<-c('Adiposity','Glycaemia','Hypertension','High triglycerides','Low HDL-C');ys<-c(195,176,157,138,119);lim<-c(.06,5.2)
 for(j in 1:2){x0<-plotx[[j]][1];x1<-plotx[[j]][2];ln(logx(1,x0,x1,lim),112,logx(1,x0,x1,lim),203,grey,.75,2)
  for(i in 1:5){cy<-ys[i];if(j==1)txt(labels[i],3,cy,9,i==2,just='left')
   if(i<5)ln(x0-2,cy-9.5,pp[j]+6,cy-9.5,'#EDF0F2',.45)
   for(k in 1:2){set<-c('Complete-component','Multiple imputation')[k];z<-comp[comp$dataset==set & comp$outcome==paste0('f_ge',j+1) & comp$term==paste0('cmrf',i),];y<-cy+if(k==1)3.1 else -3.1
    if(j==1)txt(c('CC','MI')[k],33,y,col=grey)
    forest(z,y,x0,x1,lim,col[j],if(k==1)21 else 24,z$p_fdr_all_10<.05,paste('component',i,j,k));label_fit(paste('OR',i,j,k),efmt(z),ox[j],y,26.5);label_fit(paste('P',i,j,k),pfmt(z$p_fdr_all_10),pp[j],y,13.5)
   }
  }
  ln(x0,108,x1,108,ink);for(v in c(.1,.5,1,2,4)){x<-logx(v,x0,x1,lim);ln(x,108,x,106.8,ink);txt(v,x,103.7)}
  txt('Adjusted OR (log scale)',(x0+x1)/2+10,98)
 }
 txt('Filled symbols: adjusted P <0.05. Adiposity-negative: n = 13 in CC.',3,91,col=grey,just='left')
 txt('b',3,80,12,TRUE,just='left');txt('Glycaemic categories and fibrosis',11,80,10.5,TRUE,just='left')
 txt('Reference: no glycaemic abnormality · Adjusted for other four domains, age, sex and centre',11,73,col=grey,just='left')
 txt('Category',3,58,9,TRUE,just='left');txt('Set',33,58,9,TRUE)
 lim<-c(.85,7.1);levels<-c('Dysglycaemia without recorded T2DM','Recorded T2DM');ys<-c(43,24)
 for(j in 1:2){x0<-plotx[[j]][1];x1<-plotx[[j]][2]
  txt(paste0('F≥',j+1),mean(c(plotx[[j]],ox[j],pp[j])),65,10.5,TRUE,col=col[j]);txt('OR (95% CI)',ox[j],58,9,TRUE);txt('P (Holm)',pp[j],58,9,TRUE)
  ln(logx(1,x0,x1,lim),16,logx(1,x0,x1,lim),52,grey,.75,2)
  for(i in 1:2){cy<-ys[i];if(j==1)txt(if(i==1)'Dysglycaemia\nwithout recorded\nT2DM' else 'Recorded T2DM',3,cy,9,FALSE,just='left')
   if(i==1)ln(x0-2,33.5,pp[j]+6,33.5,'#EDF0F2',.45)
   for(k in 1:2){set<-c('Complete-component','Multiple imputation')[k];z<-gly[gly$dataset==set & gly$outcome==paste0('f_ge',j+1) & gly$level==levels[i],];y<-cy+if(k==1)3.1 else -3.1
    if(j==1)txt(c('CC','MI')[k],33,y,col=grey)
    forest(z,y,x0,x1,lim,col[j],if(k==1)21 else 24,z$p_holm_all_4<.05,paste('gly',i,j,k));label_fit(paste('glyOR',i,j,k),efmt(z),ox[j],y,26.5);label_fit(paste('glyP',i,j,k),pfmt(z$p_holm_all_4),pp[j],y,13.5)
   }
  }
  ln(x0,13,x1,13,ink);for(v in c(1,2,4,6)){x<-logx(v,x0,x1,lim);ln(x,13,x,11.8,ink);txt(v,x,9)}
 }
 txt('Adjusted OR vs no glycaemic abnormality (log scale)',105,3)
}
export<-function(draw,name,h){base<-file.path(out,'figures',name)
 cairo_pdf(paste0(base,'.pdf'),width=W/25.4,height=h/25.4,family='Arial');draw();dev.off()
 svglite(paste0(base,'.svg'),width=W/25.4,height=h/25.4,bg='white');draw();dev.off()
 agg_png(paste0(base,'.png'),width=W,height=h,units='mm',res=300,background='white');draw();dev.off()
 agg_tiff(paste0(base,'.tiff'),width=W,height=h,units='mm',res=300,compression='lzw',background='white');draw();dev.off()
}
export(draw3,'Figure3_steatosis',216);export(draw4,'Figure4_components_and_glycaemia',246)
g<-do.call(rbind,geometry);write.csv(g,file.path(out,'qa/label_widths.csv'),row.names=FALSE);ck('all table labels fit',all(g$width_mm<g$available_mm))
d<-do.call(rbind,drawn);write.csv(d,file.path(out,'qa/forest_geometry.csv'),row.names=FALSE)
ck('all31modelpoints',nrow(d)==31);ck('allCIinside',all(d$ci_low>=d$plot_low & d$ci_high<=d$plot_high))
write.csv(data.frame(check=names(checks),pass=unlist(checks)),file.path(out,'qa/numerical_checks.csv'),row.names=FALSE)
write.csv(data.frame(grade=paste0('S',0:3),hex=scol),file.path(out,'source_data/steatosis_palette.csv'),row.names=FALSE)
write.csv(data.frame(outcome=c('F≥2','F≥3','S≥2'),hex=c(blue,brick,blue)),file.path(out,'source_data/outcome_palette.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(out,'qa/sessionInfo.txt'));cat(length(checks),'checks passed; both figures exported.\n')
