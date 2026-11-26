# R-only figure for exact patterns versus all other participants. Author: mjwang.
arg <- grep('^--file=', commandArgs(FALSE), value=TRUE)
script <- normalizePath(sub('^--file=', '', arg))
root <- normalizePath(file.path(dirname(script), '../../../../..'))
version <- basename(dirname(script))
pass <- commandArgs(TRUE); pass <- if(length(pass)) pass[1] else 'v1'
out <- file.path(root,'paper1_cmrf_pathology/results/versions',version)
dest <- file.path(out,'figures',pass)
if(dir.exists(dest) && length(list.files(dest))) stop('Use a new figure pass.')
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(grid)})
r <- read.csv(file.path(out,'source_data/pattern_vs_rest_2x2.csv'))
key <- read.csv(file.path(out,'source_data/plotted_patterns.csv'))
key$y <- 12-key$order
r$y <- 12-r$order
ink <- '#23374A'; blue <- '#21669D'; orange <- '#C87925'; grey <- '#788794'
separators <- c(10.5,8.5,4.5,1.5)
base <- theme_classic(base_family='Arial',base_size=12)+
  theme(text=element_text(color=ink),axis.text=element_text(color=ink),
        plot.margin=margin(4,5,4,5),panel.grid=element_blank(),
        axis.ticks=element_line(color='#98A6B2',linewidth=.3),
        axis.line=element_line(color='#98A6B2',linewidth=.3))
mat <- do.call(rbind,lapply(seq_len(nrow(key)),function(i){
  letters <- c('A','G','H','T','L')
  data.frame(order=key$order[i],y=key$y[i],component=letters,
             x=seq(4.65,8.25,length.out=5),
             on=letters %in% strsplit(key$pattern[i],'+',fixed=TRUE)[[1]])
}))
conn <- do.call(rbind,lapply(split(mat,mat$order),function(z){
  z<-z[z$on,];data.frame(y=z$y[1],xmin=min(z$x),xmax=max(z$x))
}))
left <- ggplot()+
  geom_hline(yintercept=separators,color='#DFE5EB',linewidth=.35)+
  geom_segment(data=conn,aes(x=xmin,xend=xmax,y=y,yend=y),color=blue,linewidth=.6)+
  geom_point(data=mat,aes(x=x,y=y,color=on),size=2.5)+
  scale_color_manual(values=c('TRUE'=blue,'FALSE'='#DDE5EB'),guide='none')+
  geom_text(data=key,aes(x=0,y=y,label=pattern),family='Arial',hjust=0,size=3.8,color=ink)+
  geom_text(data=key,aes(x=3.62,y=y,label=n),family='Arial',hjust=.5,size=3.8,color=ink)+
  annotate('text',x=0,y=12.0,label='Pattern',hjust=0,family='Arial',fontface='bold',size=4.0,color=ink)+
  annotate('text',x=3.62,y=12.0,label='n',family='Arial',fontface='bold',size=4.0,color=ink)+
  annotate('text',x=seq(4.65,8.25,length.out=5),y=12.0,label=c('A','G','H','T','L'),family='Arial',fontface='bold',size=3.9,color=ink)+
  scale_x_continuous(limits=c(-.05,8.65),expand=c(0,0))+
  scale_y_continuous(limits=c(.2,12.55),expand=c(0,0))+
  labs(title='a  Complete CMRF combinations')+
  base+theme(axis.text=element_blank(),axis.title=element_blank(),axis.ticks=element_blank(),
             axis.line=element_blank(),plot.title=element_text(size=15,face='bold'))

endpoint_plot <- function(ep,color,title) {
  z <- r[r$endpoint==ep,]
  z$p_text <- ifelse(z$p_raw<.001,'<0.001',sprintf('%.3f',z$p_raw))
  z$p_text <- paste(z$p_text,z$significance)
  ggplot(z)+
    geom_hline(yintercept=separators,color='#DFE5EB',linewidth=.35)+
    geom_vline(xintercept=c(0,25,50,75,100),color='#EEF1F4',linewidth=.25)+
    geom_segment(aes(x=100*proportion_pattern,xend=100*proportion_other,y=y+.12,yend=y-.12),
                 color='#9AA7B3',linewidth=.6)+
    geom_errorbar(aes(xmin=100*lower_other,xmax=100*upper_other,y=y-.12),
                  orientation='y',width=0,color=grey,linewidth=.5)+
    geom_point(aes(x=100*proportion_other,y=y-.12),shape=21,fill='white',color=grey,size=2.3,stroke=.7)+
    geom_errorbar(aes(xmin=100*lower_pattern,xmax=100*upper_pattern,y=y+.12),
                  orientation='y',width=0,color=color,linewidth=.65)+
    geom_point(aes(x=100*proportion_pattern,y=y+.12),shape=16,color=color,size=2.3)+
    geom_text(aes(x=117,y=y,label=p_text,fontface=ifelse(p_raw<.05,'bold','plain')),
              family='Arial',size=3.9,hjust=.5,color=ink)+
    annotate('text',x=50,y=12,label='Prevalence (95% CI)',family='Arial',fontface='bold',size=4.0,color=ink)+
    annotate('text',x=117,y=12,label='Raw P',family='Arial',fontface='bold',size=4.0,color=ink)+
    scale_x_continuous(limits=c(-3,135),breaks=c(0,25,50,75,100),labels=function(x)paste0(x,'%'),expand=c(0,0))+
    scale_y_continuous(limits=c(.2,12.55),expand=c(0,0))+
    labs(title=title,x=NULL,y=NULL)+base+
    theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          axis.text.x=element_text(size=11),plot.title=element_text(size=16,face='bold',color=color))
}
p2 <- endpoint_plot('F>=2',blue,'b  Significant fibrosis: F≥2')
p3 <- endpoint_plot('F>=3',orange,'c  Advanced fibrosis: F≥3')
head <- ggplot()+
  annotate('text',x=0,y=.85,label='Each combination versus all other participants (N = 935)',
           hjust=0,family='Arial',fontface='bold',size=5.2,color=ink)+
  annotate('point',x=.012,y=.16,shape=16,size=2.5,color=ink)+
  annotate('text',x=.03,y=.16,label='This combination',hjust=0,family='Arial',size=4,color=ink)+
  annotate('point',x=.235,y=.16,shape=21,fill='white',color=grey,size=2.5,stroke=.8)+
  annotate('text',x=.255,y=.16,label='All other participants',hjust=0,family='Arial',size=4,color=ink)+
  annotate('text',x=.58,y=.16,label='Points and whiskers: observed prevalence and Wilson 95% CI',
           hjust=0,family='Arial',size=3.75,color=ink)+
  xlim(0,1)+ylim(0,1)+theme_void()+theme(plot.margin=margin(2,10,3,8))
foot <- ggplot()+
  annotate('text',x=0,y=.82,label='A, adiposity; G, glycaemia; H, hypertension; T, high triglycerides; L, low HDL-C.',
           hjust=0,family='Arial',size=3.55,color=ink)+
  annotate('text',x=0,y=.45,label='Pearson χ²: * P<0.05, ** P<0.01, *** P<0.001; ns P≥0.05. No covariate or multiplicity adjustment.',
           hjust=0,family='Arial',size=3.55,color=ink)+
  annotate('text',x=0,y=.08,label='11 common patterns (n≥30) cover 839/935 participants; each comparator includes every other pattern, including rare patterns.',
           hjust=0,family='Arial',size=3.55,color=ink)+
  xlim(0,1)+ylim(0,1)+theme_void()+theme(plot.margin=margin(3,8,3,8))
body <- (left|p2|p3)+plot_layout(widths=c(3.7,4.2,4.2))
fig <- head/body/foot+plot_layout(heights=c(.62,4.75,.64))
width <- 12.4; height <- 6.15
for(ext in c('svg','pdf','png')) {
  ggsave(file.path(dest,paste0('CMRF_pattern_vs_rest_chisq.',ext)),fig,width=width,height=height,
         device=switch(ext,svg=svglite::svglite,pdf=cairo_pdf,png=ragg::agg_png),dpi=300,bg='white')
}
ggsave(file.path(dest,'preview.png'),fig,width=width,height=height,device=ragg::agg_png,dpi=130,bg='white')
write.csv(r,file.path(dest,'plotted_values.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(dest,'R_figure_sessionInfo.txt'))
