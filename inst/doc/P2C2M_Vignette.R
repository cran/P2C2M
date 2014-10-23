### R code from vignette source 'P2C2M_Vignette.Rnw'

###################################################
### code chunk number 1: P2C2M_Vignette.Rnw:17-18
###################################################
  options(width=60, size="scriptsize")


###################################################
### code chunk number 2: P2C2M_Vignette.Rnw:42-50
###################################################
library(P2C2M)
data(viz_example_1)
inp = viz_example_1
alpha = 0.05
inData = qnts = df = titles = list()
df$lwr = df$upr = list()
titles$sorted = sprintf("gene%02d_sorted", c(1:10))
titles$unsorted = sprintf("gene%02d_unsorted", c(1:10))


###################################################
### code chunk number 3: P2C2M_Vignette.Rnw:60-66
###################################################
colnames(inp$sorted) = titles$sorted
inData$sorted = stack(as.data.frame(inp$sorted))
colnames(inData$sorted) = c("value", "gene")
qnts$sorted = apply(inp$sorted, 2, quantile, c(alpha, 1-alpha), na.rm=TRUE)
df$lwr$sorted = data.frame(lwrQntl=qnts$sorted[1,], gene=names(qnts$sorted[1,]))
df$upr$sorted = data.frame(uprQntl=qnts$sorted[2,], gene=names(qnts$sorted[2,]))


###################################################
### code chunk number 4: P2C2M_Vignette.Rnw:76-82
###################################################
colnames(inp$unsorted) = titles$unsorted
inData$unsorted = stack(as.data.frame(inp$unsorted))
colnames(inData$unsorted) = c("value", "gene")
qnts$unsorted = apply(inp$unsorted, 2, quantile, c(alpha, 1-alpha), na.rm=TRUE)
df$lwr$unsorted = data.frame(lwrQntl=qnts$unsorted[1,], gene=names(qnts$unsorted[1,]))
df$upr$unsorted = data.frame(uprQntl=qnts$unsorted[2,], gene=names(qnts$unsorted[2,]))


###################################################
### code chunk number 5: P2C2M_Vignette.Rnw:91-95
###################################################
inData = rbind(inData$sorted, inData$unsorted)
dfLwr = rbind(df$lwr$sorted, df$lwr$unsorted)
dfUpr = rbind(df$upr$sorted, df$upr$unsorted)
inData$gene = factor(inData$gene, levels = sort(c(titles$sorted, titles$unsorted)))


###################################################
### code chunk number 6: P2C2M_Vignette.Rnw:105-125
###################################################
library(ggplot2)
ggplot(data=inData, aes(x=value)) +
  geom_density() +
  facet_grid(gene~.) +
  labs(x="Difference values") + 
  ggtitle(expression(atop("Ranked vs. Unranked Distributions", 
                     atop(italic("Descriptive Statistic: RAY"), "")))) +

  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y=element_text(angle=0),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        strip.background=element_rect(fill="white")) +
  # Limits on the x-axis improve the visualization
  xlim(-500, 500) + 
  geom_vline(xintercept=0, linetype = "dashed") + 
  geom_vline(aes(xintercept=lwrQntl), dfLwr, color="grey") +
  geom_vline(aes(xintercept=uprQntl), dfUpr, color="grey")


###################################################
### code chunk number 7: P2C2M_Vignette.Rnw:131-151
###################################################
library(ggplot2)
ggplot(data=inData, aes(x=value)) +
  geom_density() +
  facet_grid(gene~.) +
  labs(x="Difference values") + 
  ggtitle(expression(atop("Ranked vs. Unranked Distributions", 
                     atop(italic("Descriptive Statistic: RAY"), "")))) +

  theme_bw() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        strip.text.y=element_text(angle=0),
        panel.grid.major.x=element_blank(),
        panel.grid.major.y=element_blank(),
        strip.background=element_rect(fill="white")) +
  # Limits on the x-axis helps the visualization
  xlim(-500, 500) + 
  geom_vline(xintercept=0, linetype = "dashed") + 
  geom_vline(aes(xintercept=lwrQntl), dfLwr, color="grey") +
  geom_vline(aes(xintercept=uprQntl), dfUpr, color="grey")


###################################################
### code chunk number 8: P2C2M_Vignette.Rnw:167-170
###################################################
library(P2C2M)
data(viz_example_2)
inp = viz_example_2


###################################################
### code chunk number 9: P2C2M_Vignette.Rnw:181-198
###################################################
myfunc = function(inData, simNum){
  handle = inData
  colnames(handle) = c("gtp", "ray", "ndc", "gsi")
  # Convert results into presence/absence matrix
  handle[!grepl("n.s.", handle)] = 1
  handle[grepl("n.s.", handle)] = 0
  # Stack the individual descriptive statistics
  handle = stack(data.frame(handle, stringsAsFactors=FALSE))
  colnames(handle)[1] = "value"
  colnames(handle)[2] = "stat"
  # Add gene identifiers (under the assumption that there are 10 genes)
  handle[,3] = rep(c(1:10), 4)
  colnames(handle)[3] = "gene"
  handle[,4] = simNum
  colnames(handle)[4] = "sim"
return(handle)
}


###################################################
### code chunk number 10: P2C2M_Vignette.Rnw:208-223
###################################################
highL = list()
sims = as.numeric(names(inp$High))
for (i in 1:length(inp$High)) {highL[[i]] = myfunc(inp$High[[i]], sims[i])}
High = do.call("rbind", highL)
High[,ncol(High)+1] = "High_Subst_Rate"
colnames(High)[ncol(High)] = "ratetype"

lowL = list()
sims = as.numeric(names(inp$Low))
for (i in 1:length(inp$Low)) {lowL[[i]] = myfunc(inp$Low[[i]], sims[i])}
Low = do.call("rbind", lowL)
Low[,ncol(Low)+1] = "Low_Subst_Rate"
colnames(Low)[ncol(Low)] = "ratetype"

inData = rbind(High, Low)


###################################################
### code chunk number 11: P2C2M_Vignette.Rnw:233-245
###################################################
library(ggplot2)
ggplot(data=inData, aes(x=sim,y=gene)) + 
  geom_point(aes(colour=value), size = 3) +
  scale_colour_manual(values = c(NA,'black')) + 
  facet_grid(stat~ratetype) + 
  ggtitle(expression(atop("Distribution of False Positives",
          atop(italic("Alpha=0.1") , "")))) +
  theme_bw() + 
  scale_x_discrete(breaks=c(1:5), labels=c(1:5)) +
  scale_y_discrete(breaks=c(10:1), labels=c(10:1)) +
  theme(strip.background = element_rect(fill="white")
  )


###################################################
### code chunk number 12: P2C2M_Vignette.Rnw:252-264
###################################################
library(ggplot2)
ggplot(data=inData, aes(x=sim,y=gene)) + 
  geom_point(aes(colour=value), size = 3) +
  scale_colour_manual(values = c(NA,'black')) + 
  facet_grid(stat~ratetype) + 
  ggtitle(expression(atop("Distribution of False Positives",
          atop(italic("Alpha=0.1") , "")))) +
  theme_bw() + 
  scale_x_discrete(breaks=c(1:5), labels=c(1:5)) +
  scale_y_discrete(breaks=c(10:1), labels=c(10:1)) +
  theme(strip.background = element_rect(fill="white")
  )


