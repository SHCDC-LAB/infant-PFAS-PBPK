# FAST175 post-hoc convergence and parameter-rank robustness (NO PBPK / NO ODE)
# Inputs: FAST175_population_metrics_by_draw.csv + FAST175_LHS_design.csv
# Optional: global_sensitivity_FAST175.xlsx is NOT required.
# Run from project root: source('R/10_FAST175_posthoc_convergence_noODE.R')
# Outputs: results/.../global_sensitivity_FAST175/full/posthoc_diagnostics/
options(stringsAsFactors=FALSE)
root <- normalizePath(Sys.getenv('PFAS_PBPK_ROOT', unset=getwd()),winslash='/',mustWork=FALSE)
mode <- Sys.getenv('PFAS_10_POSTHOC_MODE','full')
base <- Sys.getenv('PFAS_10_POSTHOC_DIR',unset=file.path(root,'results','PFAS_infant_PBPK_S1_renal_totalplasma_v2','global_sensitivity_FAST175',mode))
design_file <- file.path(base,'FAST175_LHS_design.csv')
endpoint_file <- file.path(base,'FAST175_population_metrics_by_draw.csv')
if(!file.exists(design_file) || !file.exists(endpoint_file)) stop('Missing already-computed CSV. Need:\n',design_file,'\n',endpoint_file,'\nNO simulation was started.')
out <- file.path(base,'posthoc_diagnostics');dir.create(out,showWarnings=FALSE,recursive=TRUE)
cat('[POSTHOC] CSV-only diagnostics. NO deSolve, NO PBPK ODE.\n')
design <- read.csv(design_file,check.names=FALSE)
eps <- read.csv(endpoint_file,check.names=FALSE)
params <- c('Milk_Intake','R_CM','R_BM','Tmc','Kt','fu','GFR','P_li','P_kd','P_re','Beta_Proxy')
metrics <- c('C_Plasma_12m','C_Plasma_6m','AUC_Plasma_0_12m','C_Liver_12m','C_Kidney_12m')
chem_order <- c('PFOA','PFNA','PFDA','PFUdA','PFHxS','PFOS')
need <- function(d,cols,label){z<-setdiff(cols,names(d));if(length(z))stop(label,' missing columns: ',paste(z,collapse=', '))}
need(design,c('Draw',params),'design');need(eps,c('Draw','Compound',metrics),'endpoints')
if(anyDuplicated(design$Draw) || anyDuplicated(eps[c('Compound','Draw')]))stop('Duplicate draw keys')
if(any(!is.finite(as.matrix(design[params]))))stop('Nonfinite design parameter')
chem <- intersect(chem_order,unique(eps$Compound))
if(!length(chem))stop('No matching PFAS')
STAGES <- c(300L,400L,500L,600L)
B <- as.integer(Sys.getenv('PFAS_10_POSTHOC_BOOT','200'))
if(is.na(B)||B<0||B>2000)stop('PFAS_10_POSTHOC_BOOT must be 0..2000')
SEED <- as.integer(Sys.getenv('PFAS_10_POSTHOC_SEED','20261106'))
set.seed(SEED)
active <- function(cc)if(cc %in% c('PFOA','PFHxS','PFOS'))setdiff(params,'Beta_Proxy')else params
# Partial rank correlation from inverse correlation matrix; identical definition to
# rank residualization (apart from roundoff). Faster for repeated resamples.
prcc <- function(x,y){
  rx<-apply(x,2,rank,ties.method='average');ry<-rank(y,ties.method='average')
  C<-stats::cor(cbind(rx,ry));if(any(!is.finite(C)))return(rep(NA_real_,ncol(x)))
  P<-tryCatch(solve(C),error=function(e)NULL)
  if(is.null(P))return(rep(NA_real_,ncol(x)))
  k<-ncol(x);p<--P[seq_len(k),k+1]/sqrt(diag(P)[seq_len(k)]*P[k+1,k+1]);as.numeric(p)
}
qfun<-function(y) unname(stats::quantile(y,c(.025,.5,.975),names=FALSE,type=7))
get_slice <- function(cc,n,p){
  di<-design[design$Draw<=n,c('Draw',p),drop=FALSE]
  ep<-eps[eps$Compound==cc & eps$Draw<=n,c('Draw',metrics),drop=FALSE]
  d<-merge(di,ep,by='Draw',sort=TRUE)
  if(nrow(d)!=n||!identical(d$Draw,seq_len(n)))stop('Missing cumulative draws for ',cc,' N=',n)
  if(any(!is.finite(as.matrix(d[c(p,metrics)]))))stop('Nonfinite observations for ',cc,' N=',n)
  d
}
stats_list<-list();quant_list<-list();compare_list<-list();rank_boot<-list();comp_summary<-list();top_list<-list();ix<-0L
for(cc in chem){
  p<-active(cc);available<-sort(unique(eps$Draw[eps$Compound==cc]));maxn<-max(available)
  stages<-STAGES[STAGES<=maxn]
  if(length(stages)<2)stop('Need >=2 stages for ',cc)
  if(!identical(available,seq_len(maxn)))stop('Non-consecutive draw IDs in ',cc)
  cat('[POSTHOC] ',cc,' | ',maxn,' existing draws\n',sep='')
  for(n in stages){
    d<-get_slice(cc,n,p)
    for(m in metrics){
      vals<-prcc(as.matrix(d[p]),d[[m]]);names(vals)<-p
      ord<-names(sort(abs(vals),decreasing=TRUE,na.last=TRUE))
      for(j in seq_along(p)){
        stats_list[[length(stats_list)+1L]]<-data.frame(Compound=cc,Metric=m,N_Draw=n,Parameter=p[j],PRCC=vals[j],Abs_PRCC=abs(vals[j]),Rank_Abs=match(p[j],ord))
      }
      q<-qfun(d[[m]])
      quant_list[[length(quant_list)+1L]]<-data.frame(Compound=cc,Metric=m,N_Draw=n,P2_5=q[1],Median=q[2],P97_5=q[3])
    }
  }
  # Bootstrap at FINAL N: empirical stability of rank under resampling existing draws.
  # These are Monte Carlo rank-stability diagnostics, not parameter CIs.
  n<-tail(stages,1);d<-get_slice(cc,n,p)
  if(B>0){
    for(m in metrics){
      boot<-matrix(NA_real_,nrow=B,ncol=length(p),dimnames=list(NULL,p))
      for(b in seq_len(B)){
        ind<-sample.int(n,n,replace=TRUE)
        boot[b,]<-prcc(as.matrix(d[ind,p,drop=FALSE]),d[[m]][ind])
      }
      final<-prcc(as.matrix(d[p]),d[[m]]);names(final)<-p
      for(j in seq_along(p)){
        topfreq<-mean(apply(abs(boot),1,function(v) p[j] %in% names(sort(v,decreasing=TRUE)[seq_len(min(5,length(p)))])),na.rm=TRUE)
        rank_boot[[length(rank_boot)+1L]]<-data.frame(Compound=cc,Metric=m,Parameter=p[j],N_Draw=n,Bootstrap_Replicates=B,PRCC=final[j],
                                                      Boot_PRCC_P2_5=unname(quantile(boot[,j],.025,na.rm=TRUE)),Boot_PRCC_P97_5=unname(quantile(boot[,j],.975,na.rm=TRUE)),
                                                      Top5_Inclusion_Frequency=topfreq,Sign_Agreement_Frequency=mean(sign(boot[,j])==sign(final[j]),na.rm=TRUE))
      }
    }
  }
}
prcc_table<-do.call(rbind,stats_list);quant_table<-do.call(rbind,quant_list)
for(cc in chem)for(m in metrics){
  d<-prcc_table[prcc_table$Compound==cc & prcc_table$Metric==m,,drop=FALSE]
  qq<-quant_table[quant_table$Compound==cc & quant_table$Metric==m,,drop=FALSE]
  stages<-sort(unique(d$N_Draw))
  for(i in 2:length(stages)){
    prev<-d[d$N_Draw==stages[i-1],,drop=FALSE];cur<-d[d$N_Draw==stages[i],,drop=FALSE]
    prev<-prev[match(cur$Parameter,prev$Parameter),,drop=FALSE]
    oldtop<-prev$Parameter[order(-prev$Abs_PRCC)][seq_len(min(5,nrow(prev)))]
    newtop<-cur$Parameter[order(-cur$Abs_PRCC)][seq_len(min(5,nrow(cur)))]
    uniontop<-union(oldtop,newtop);sel<-cur$Parameter %in% uniontop
    oldq<-as.numeric(unlist(qq[qq$N_Draw==stages[i-1],c('P2_5','Median','P97_5')],use.names=FALSE));newq<-as.numeric(unlist(qq[qq$N_Draw==stages[i],c('P2_5','Median','P97_5')],use.names=FALSE))
    delta<-abs(newq-oldq)/pmax(abs(newq),abs(oldq),1e-12)
    rankcorr<-suppressWarnings(cor(prev$Rank_Abs,cur$Rank_Abs,method='spearman'))
    gaps<-sort(cur$Abs_PRCC,decreasing=TRUE)
    k<-min(5,length(gaps)-1)
    compare_list[[length(compare_list)+1L]]<-data.frame(Compound=cc,Metric=m,Previous_N=stages[i-1],Current_N=stages[i],
                                                        Top5_Same=setequal(oldtop,newtop),Top5_Overlap=length(intersect(oldtop,newtop)),Top5_Boundary_Gap=gaps[k]-gaps[k+1],
                                                        Signs_Same=all(sign(prev$PRCC[sel])==sign(cur$PRCC[sel])),Max_Top_PRCC_Delta=max(abs(prev$PRCC[sel]-cur$PRCC[sel])),
                                                        AbsRank_Spearman=rankcorr,P2_5_RelChange=delta[1],Median_RelChange=delta[2],P97_5_RelChange=delta[3],
                                                        Max_Quantile_RelChange=max(delta),PRCC_Stable=setequal(oldtop,newtop)&&all(sign(prev$PRCC[sel])==sign(cur$PRCC[sel]))&&max(abs(prev$PRCC[sel]-cur$PRCC[sel]))<=.05,
                                                        Quantile_Stable=all(delta<=.05),Old_Top5=paste(oldtop,collapse=';'),New_Top5=paste(newtop,collapse=';'))
  }
  final<-d[d$N_Draw==max(stages),,drop=FALSE]
  final<-final[order(final$Rank_Abs),,drop=FALSE]
  top_list[[length(top_list)+1L]]<-head(final,6)
}
comparisons<-do.call(rbind,compare_list);top6<-do.call(rbind,top_list)
last<-do.call(rbind,lapply(split(comparisons,list(comparisons$Compound,comparisons$Metric),drop=TRUE),function(d)d[which.max(d$Current_N),,drop=FALSE]))
last$Endpoint_Both_Stable<-last$PRCC_Stable & last$Quantile_Stable
summary<-aggregate(cbind(PRCC_Stable,Quantile_Stable,Endpoint_Both_Stable)~Compound,data=last,FUN=sum)
summary$N_Endpoints<-5L
summary$Fully_Stable_Last_Comparison<-summary$Endpoint_Both_Stable==5L
summary$Note<-'Last-comparison only; not equivalent to two consecutive stable eligible comparisons'
boot_table<-if(length(rank_boot))do.call(rbind,rank_boot)else data.frame()
write.csv(prcc_table,file.path(out,'PRCC_by_stage.csv'),row.names=FALSE)
write.csv(quant_table,file.path(out,'Quantiles_by_stage.csv'),row.names=FALSE)
write.csv(comparisons,file.path(out,'Stage_pair_comparisons.csv'),row.names=FALSE)
write.csv(last,file.path(out,'Final_stage_endpoint_diagnostics.csv'),row.names=FALSE)
write.csv(summary,file.path(out,'Compound_summary.csv'),row.names=FALSE)
write.csv(top6,file.path(out,'Final_Top6_PRCC.csv'),row.names=FALSE)
if(nrow(boot_table))write.csv(boot_table,file.path(out,'Bootstrap_rank_robustness.csv'),row.names=FALSE)
metadata<-data.frame(Item=c('Method','Source','Stages','PRCC','Top5','Quantile criterion','Bootstrap','Interpretation','Important limitation'),
                     Value=c('Retrospective posthoc only, NO ODE','Existing FAST175 LHS and population metrics CSV',paste(STAGES,collapse=','),
                             'Partial Pearson correlation among ranks, inverse rank correlation matrix','Absolute PRCC rank; compare top5 set and sign',
                             'P2.5, median, P97.5 relative changes <=5% per adjacent cumulative stage',paste0(B,' ordinary draw bootstrap replicates, seed ',SEED),
                             'Bootstrap measures sampling stability of PRCC/ranks, not parameter probability/95% model interval',
                             'Cumulative stages share draws and 500 is an extra posthoc checkpoint; tests are descriptive, not independent'))
if(requireNamespace('writexl',quietly=TRUE)){
  wb<-list(README=metadata,Compound_summary=summary,Final_endpoint=last,Stage_comparisons=comparisons,
           PRCC_by_stage=prcc_table,Quantiles_by_stage=quant_table,Top6_PRCC=top6)
  if(nrow(boot_table))wb$Bootstrap_rank=boot_table
  writexl::write_xlsx(wb,file.path(out,'FAST175_posthoc_convergence.xlsx'))
}
# A small plot with base graphics, so no plotting package dependency.
png(file.path(out,'Final_12m_plasma_PRCC_top6.png'),width=1900,height=1250,res=180)
par(mfrow=c(2,3),mar=c(5,9,3,1))
for(cc in chem){
  z<-top6[top6$Compound==cc & top6$Metric=='C_Plasma_12m',,drop=FALSE]
  z<-z[order(z$PRCC),,drop=FALSE]
  if(!nrow(z)){plot.new();title(cc);next}
  barplot(z$PRCC,names.arg=z$Parameter,horiz=TRUE,las=1,cex.names=.8,
          xlim=c(-1,1),main=cc,col=ifelse(z$PRCC>=0,'#358A78','#BE8062'))
  abline(v=0,lty=2)
}
dev.off()
cat('[POSTHOC] DONE. Files: ',out,'\n',sep='')
print(summary)
