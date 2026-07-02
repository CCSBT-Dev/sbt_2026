
load("inst/extdata/grid_tmbfit.rda")
windows(height=3,width=9)
par(mfrow=c(1,4),cex.lab=1.5,cex.axis=1.2)
barplot(table(exp(grid_tmbfit$samples[,,'par_log_h'])),xlab='h')
barplot(table(exp(grid_tmbfit$samples[,,'par_log_m0'])),xlab='m0')
barplot(table(factor(exp(grid_tmbfit$samples[,,'par_log_m10']),
                     levels=c(.065,.085,.105))),xlab='m10')
barplot(table(exp(grid_tmbfit$samples[,,'par_log_psi'])),xlab='psi')

# get best fitting grid cell
load("inst/extdata/grid_check.rda")
load("inst/extdata/grid_list.rda")

grid_summary <- grid_check$grid_summary
index <- which(grid_summary$nll==min(grid_summary$nll))
report_best <- grid_list[[index]]$report()
save(report_best,HSPs,POPs, file="x:/work/ccsbt/ommp 2026/report_best.rda" )


# plots from sbt package:
plot_cpue(data = data, object = obj, nsim = 25)
plot_cpue_residuals(data = data, obj = obj, type = "OSA")

plot_aerial_survey(data = data, object = obj, nsim = 25)

plot_af(data = data, object = obj, fishery = "Australian")
plot_af(data = data, object = obj, fishery = "Indonesian")
plot_lf(data = data, object = obj, fishery = "LL1")
plot_lf(data = data, object = obj, fishery = "LL2")
#plot_lf(data = data, object = obj, fishery = "LL3")
#plot_lf(data = data, object = obj, fishery = "LL4")
plot_lf(data = data, object = obj, fishery = "CPUE")

p1 <- plot_selectivity(data = data, object = obj, years = 1969:2025, fisheries = "LL1")
p2 <- plot_selectivity(data = data, object = obj, years = 1969:2025, fisheries = "CPUE")
p1 + p2
plot_selectivity(data = data, object = obj, years = 1990:2025)
# plot_selectivity(data = data, object = obj, years = 1969:2022, fisheries = "LL3")
# plot_selectivity(data = data, object = obj, years = 1969:2022, fisheries = "LL4")

plot_biomass_spawning(data_list = list(data), object_list = list(obj))



# plot fit to GT data:
with(GTs, plot(RelYear, Nmatch, pch=16, ylim = c(0,90)))
#cv <- with(GTs, sqrt(1/Nmatch))
#ci.025 <- with(GTs,qbinom(.025,Nsam,Nmatch/Nsam))
#ci.975 <- with(GTs,qbinom(.975,Nsam,Nmatch/Nsam))
#segments(GTs$RelYear,ci.025,GTs$RelYear,ci.975,lty=2)
pred.match <- with(report_best, gt_nscan*gt_prob)
points(GTs$RelYear, pred.match, pch=4, col=4, lwd=2)
legend('topright',c("Observed","Predicted"), pch=c(16,4),col=c(1,4),pt.lwd=c(1,2),bty='n')

# plot fit to conventional data:

# aggregate over tagger groups
obs.recap <- apply(data$tag_recap_ctaa,c(1,3,4),sum)
pred.recap <- apply(report_best$tag_pred,c(1,3,4),sum)

par(mfrow=c(6,3),mai=c(.3,.3,.3,.2),omi=c(.5,.5,.5,.5))
for(k in 1:6) {
  for(i in 1:3) {
    x=barplot(obs.recap[k,i,(i+1):7],ylim=c(0,max(obs.recap[k,i,],pred.recap[k,i,])*1.2),names=(i+1):7)
    points(x, pred.recap[k,i,(i+1):7],pch=16,cex=1.3)
    if(k==1) mtext(paste("Rel age",i),side=3,outer=F)
  }
  if(i==3) mtext(c(1989:1994)[k],side=4, outer=F)
}
mtext("Recapture age",side=1,outer=T, line=0.8)
mtext("Number of recaptures",side=2, outer=T,line=0.8)

