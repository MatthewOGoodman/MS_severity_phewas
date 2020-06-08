Date.of.analysis = "2018-03-04"

# results = read.table('MSSS.lcl.logreg.60.no.excl.txt')
# results = read.csv('MSSS.lcl.logreg.120.no.excl.csv')
# results = read.csv('PheWAS output_July 25th - no notes.csv')
# results = read.csv('Meta.analysis.Fixed.effects.MSSS.lcl.logreg.MA120.no.excl.csv')
# results = read.csv('2018-03-03 primary BC.replicated.significant.w.meta.analysis.(Table2).csv', header = TRUE)
result.name = paste(Date.of.analysis, 'primary.MA.complete.list.(TableS1v2).csv', sep = '.')
MA.samp.size = 3439; MA.thresh = .03 * 3439
samp.size = MA.samp.size; thresh = MA.thresh
cohort.name = 'MA Cohort' #Meta Analysis

result.name = paste(Date.of.analysis, 'primary.BC.complete.list.(TableS2v2).csv', sep = '.')
BC.samp.size = 4867; BC.thresh = .03 * 4867
samp.size = BC.samp.size; thresh = BC.thresh
cohort.name = 'BC Cohort' #Meta Analysis

results = read.csv( result.name, colClasses = c('phecode' = 'character'))
colnames(results)[colnames(results) == 'phecode'] = 'phenotype'

head(results)
# rename columns
results$prev = results[,grepl("prev", colnames(results))] / 100
results$eff = results$OR
results$bonferroni = results$p<.05/sum(!is.na(results$p)) & !is.na(results$p)

# create color map
map = data.frame(
  groups = c('INFECTIOUS', 'NEOPLASTIC', 'ENDOCRINE&METABOLIC', 'HEMATOPOIETIC',
    'PSYCHIATRIC', 'NEUROLOGIC', 'CARDIOVASCULAR', 'PULMONARY', 'DIGESTIVE',
    'GENITOURINARY', 'DERMATOLOGIC', 'MUSCULOSKELETAL', 'SYMPTOMS&SIGNS', 'INJURIES'),
  colors = c('blue', 'cyan4', 'brown', 'darkorange1', 'magenta1', 'blue4', 'red', 'coral4',
    'chartreuse4', 'black', 'firebrick', 'darkolivegreen', 'purple', 'gray50'), 
  numbers = c(1:14))

  breakpoints = c(0,140,240,280,
                290,320,390,460,520,
                580, 680, 710, 780, 800,1000) #ends at 1000

groups = cut(as.numeric(as.character(results$phenotype)), breaks = breakpoints, right = F)
levels(groups) = map$colors
results$colors = groups

color = colors()[grep('orange',colors())]

savepdf<-function(x,  w = 7.5, h = 3.5, file, display=TRUE) {
  if (display){
    x;
    dev.copy2pdf(file=file, width = w, height = h, out.type = "pdf")
  }
  else {
    pdf(file=file, width = w, height = h, out.type = "pdf")
    x;
    dev.off()
  }
}

###############################################################################\

pdf(file = paste('OR vs Prev', cohort.name, '.pdf'), width = 7, height = 5)
xlim = c(0.02,max(1.25,results$prev)); ylim = c(0.8,max(1.8,results$eff))   #MA cohort

pch.str = 1:length(results$bonferroni)
sig.symb = 16; insig.symb = 1
pch.str[results$bonferroni] = sig.symb
pch.str[!results$bonferroni] = insig.symb

plot(results$prev, results$eff, xlim = xlim, ylim = ylim, log = 'x', 
     col = as.character(results$colors),
     pch = pch.str, cex = 1.0,
     main = cohort.name, xlab = 'Prevalence', ylab = 'OR')

label.codes = results$phewas_code
label.codes[!results$bonferroni] = NA
text(results$prev, results$eff, xlim = xlim, ylim = ylim,
     labels = label.codes, cex = .75, adj = c(-.1,-.4))
legend('topright', inset=c(-0.0,-.0), legend = map$groups, fill=as.character(map$colors), cex = 0.65)
legend('bottomright', legend = c('significant', 'not significant'), pch = c(sig.symb,insig.symb), cex = .85, bty = 'n')
abline(h=1.0, lty = 2)
abline(v=thresh/samp.size, lty = 1, col = 'blue')

##################################################################
getwd()
p = recordPlot()
dev.off()
p

# savepdf(p, w = 7, h = 5, file = paste('OR vs Prev', cohort.name, '.pdf'))

write.table(
  data.frame(results$phewas_code, results$phewas_description)[results$bonferroni & !is.na(results$bonferroni),],
  file = paste('OR vs Prev', cohort.name, 'key.txt')
)


################
# To be able to color a table of bonferroni results easily
plot(as.numeric(results$phewas_code[results$bonferroni & !is.na(results$bonferroni)]), 
     col = as.character(results$colors[results$bonferroni & !is.na(results$bonferroni)]),
     pch = c(19), cex =1)
text(results$phewas_code[results$bonferroni & !is.na(results$bonferroni)],
     labels = results$phewas_code[results$bonferroni & !is.na(results$bonferroni)], cex = .3, adj = c(.3,-3))

################
# The starred version of the bonferroni OR vs prev
# pch.str = 1:length(results$bonferroni)
# pch.str[results$bonferroni] = 8
# pch.str[!results$bonferroni] = 20
# 
# plot(results$prev, results$eff, xlim = xlim, ylim = ylim, log = 'x', 
#      col = as.character(results$colors),
#      pch = pch.str, cex = 1.1,
#      main = cohort.name, xlab = 'Prevalence', ylab = 'OR')
# par(new = T)
# plot(results$prev, results$eff, xlim = xlim, ylim = ylim, log = 'x', 
#      col = as.character(results$colors),
#      pch =20, cex = 1.5,
#      main = cohort.name, xlab = 'Prevalence', ylab = 'OR')
# label.codes = results$phewas_code
# label.codes[!results$bonferroni] = NA
# text(results$prev, results$eff, xlim = xlim, ylim = ylim, 
#      labels = label.codes, cex = .65, adj = c(-.4,-.4))
# legend('topright', inset=c(-0.0,-.0), legend = map$groups, fill=as.character(map$colors), cex = 0.75)
# legend('bottomright', legend = c('significant', 'not significant'), pch = c(8,20), cex = 0.75, bty = 'n')
# legend('bottomright', legend = c('significant', 'not significant'), pch = c(20,20), cex = 0.75, bty = 'n')
# abline(h=1.0, lty = 2)
# abline(v=120/4118, lty = 1, col = 'blue')

