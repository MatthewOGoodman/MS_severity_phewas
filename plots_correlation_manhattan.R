# results 6.0

# results = addPhewasDescription(phewas.results.allMSSS.25)
# results = addPhewasDescription(phewas.results.allMSSS.60.no.excl)
# results = addPhewasDescription(pmlr.results)

Date.of.analysis = "2018-03-03"
###################################################
# Load results from past Analysis
results = read.csv( paste(Date.of.analysis, 'MSSS.lcl.min.rec.25.no.excl.csv', sep = '.'),
                    header = TRUE , colClasses = c('phecode' = "character", 'description = "character'))
)
dir(); load('msss.phewas.phenotypes')
phenotypes = phewas.phenotypes.msss

MA.all = read.csv( 
  paste(Date.of.analysis, 'primary.MA.complete.list.(TableS1v2).csv', sep = '.'), 
  header = TRUE , colClasses = c('phecode' = "character", 'description' = "character"))
head(MA.all)
MA.rep = read.csv( 
  paste(Date.of.analysis, 'primary.BC.replicated.significant.w.meta.analysis.(Table2).csv', sep = '.'), 
  header = TRUE , colClasses = c('phecode' = "character", 'description.MA' = "character"))
head(MA.rep)



##################################
# Define Results
##################################
slct.res = MA.rep
selected.codes = as.character(slct.res$phecode)
pheno.select = phenotypes[ , colnames(phenotypes) %in% selected.codes]

length(names(pheno.select)) 
length(selected.codes)


cor.select = cor(pheno.select, use = 'pairwise.complete.obs', method = 'pearson')


#######################################################
# Get Phenotype Groupings
map = data.frame(
  groups = c('INFECTIOUS', 'NEOPLASTIC', 'ENDOCRINE&METABOLIC', 'HEMATOPOIETIC',
             'PSYCHIATRIC', 'NEUROLOGIC', 'CARDIOVASCULAR', 'PULMONARY', 'DIGESTIVE',
             'GENITOURINARY', 'DERMATOLOGIC', 'MUSCULOSKELETAL', 'SYMPTOMS&SIGNS', 'INJURIES'),
  colors = c('blue', 'cyan4', 'brown', 'darkorange', 'deeppink', 'navy', 'red', 'saddlebrown',
             'green3', 'black', 'red4', 'darkolivegreen', 'darkviolet', 'dimgrey'), 
  numbers = c(1:14))

breakpoints = c(0,140,240,280,
                290,320,390,460,520,
                580, 680, 710, 780, 800,1000) #ends at 1000

groups = cut(as.numeric(selected.codes), breaks = breakpoints, right = F)
levels(groups) = map$colors
slct.res$colors = groups

##############################################################
hc = hclust(as.dist(1-abs(cor.select)), method  = 'complete')
length(names(pheno.select))


pdf(file = 'p.vals.hclust.pdf', width = 7, height = 5)

plot.new()
par(mar=c(4,4,4,2))
par(fig=c(0.05,1.0,0.1,0.75), new=TRUE)

mybar = barplot(
  -log(slct.res$p.MA[hc$order]), 
  xaxs="i",  xlim = c(0,20),
  yaxs="i", ylim = c(0,175), #c(0,129),
  ylab = NULL, xlab = NULL,
  col = paste(slct.res$colors[hc$order]))

par(fig=c(0.05,1.0,0.1,0.75), new=TRUE)
plot(mybar, -log(slct.res$p.MA[hc$order])+3,
     xaxs="i",  xlim = c(0,20),
     yaxs="i", ylim = c(0,175), #c(0,129),
     ylab = "-log(p)", xlab = NA, axes = FALSE,
     col = paste(slct.res$colors[hc$order]),
     pch = (24*(slct.res$OR.MA>1)+(25*(slct.res$OR.MA<=1)))[hc$order])



axis(1, at=mybar,labels=names(pheno.select[hc$order]), las = 2)
legend('topright', legend = c('pos. assoc.', 'neg. assoc.'), pch = c(24,25), cex = 1, bty = 'n')

mtext("phenotypes", side=1, line=4)


par(mar=c(0,4,2,2))
par(fig=c(0.05,1.0,0.6,1.0), new=TRUE)
plot(hc, hang = -1, ylim = c(1,0), ylab = "absolute corr", 
     xlab=NA, sub="", main = "", axes=FALSE, xlim = c(0,sum(slct.res$bonferroni)+1),xaxs = "r" )
axis(2, at=seq(0,1,0.25),labels=seq(1,0,-.25))

dev.off()




###################################################################################################
# Do not include exclusion-induced 100% correlations
duplicate.pairs = which(cor.select>.999, arr.ind=TRUE)
unique = unique(rownames(duplicate.pairs[!duplicated(duplicate.pairs[,'col']),]))
duplicates = unique(rownames(duplicate.pairs[duplicated(duplicate.pairs[,'col']),]))
cor.select.reduced = cor(pheno.select[,!(colnames(pheno.select) %in% duplicates)],
  use = 'pairwise.complete.obs')

cor.select.1 = cor.select.reduced
cor.grps = apply(cor.select.1, 2, sum)
names(cor.grps) = paste('PW', names(cor.grps), sep = '' )



###################################################################################################
#Correlation Matrix Plots

# install.packages('corrplot')
library(corrplot)

pdf(file = 'corrplot.pdf', width = 7, height = 7)

# corrplot(cor.select, order="original", method="color", 
#          tl.pos="lt", type="full", diag = F,    
#          tl.col="black", tl.cex=0.6, tl.srt=90
# #          addCoef.col="black", addCoefasPercent = T
# #          p.mat = rcorrelations$P, sig.level=.05/(m*(m-1)/2), insig = "blank"
#          )

corrplot(cor.select, 
         order="original", method="circle", type="upper", 
         tl.col=paste(groups), tl.cex=.8, tl.offset =0.5, tl.srt = 90,
         # addCoef.col="black", addCoefasPercent = T, diag = F,    
         mar = c(1,2,1,1)
         )

mtext('plot of phenotype correlations', side=1, line=3)

dev.off()



###################################################################################################
# install.packages('PheWAS')
suppressWarnings(library(PheWAS))
# phenotypeManhattan
# phewasManhattan

# result.name = paste(Date.of.analysis, 'primary.MA.complete.list.(TableS1v2).csv', sep = '.')
# result.name = paste(Date.of.analysis, 'primary.BC.complete.list.(TableS2v2).csv', sep = '.')

results = read.csv( result.name, colClasses = c('phecode' = 'character'))
colnames(results)[colnames(results) == 'phecode'] = 'phenotype'
#############
head(results)
head(results[order(results$phenotype),])
head(results[(results$description) == 'Pneumonia',])
max(-log10(results$p))


pdf(file = paste('manhattan', result.name, 'pdf', sep = '.'), width = 10, height = 7)
# plot.new()
# par(mar=c(2,4,2,2))
# par(fig=c(0.15,0.95,0.05,0.95), new=TRUE)
phewasManhattan(results,
                # y.axis.interval = 5, max.y = 55, 
                y.axis.interval = 10, max.y = 80, max.x = 111,
                size.x.labels = 11,
                size.y.labels = 11,
                annotate.phenotype = F, 
                # suggestive.line = .05,
                # OR.direction = T,
                # annotate.level = .001,
                # annotate.angle = 0, annotate.size = 5, annotate.level, 
                # annotate.phenotype = APD, annotate.angle = 0, annotate.size = 5, annotate.level, 
                # annotate.phenotype = F, annotate.snp.w.phenotype = F, annotate.snp = F,
                # annotate.snp.angle = 0, annotate.list, annotate.only.largest = T,
                # lc.labels = F, x.group.labels = T, x.phenotype.labels = F, sizes = F,
                # direction = F, point.size = 3, use.color = T, color.palette, 
                # title = paste0("Phenotype Plot ", date()), x.axis.label = "Phenotypes", 
                # y.axis.label = "Values", y.axis.interval = 5
                title = "Phenotype Associations with MSSS"#,
                # ICD9 count threshold = 3 (inclusive), 
                #covariates = (SEX, RACE, SMK, AGE@onset, PROGRESSIVE@onset,
                #disease.duration, short.followup, ICD9\'s from PCP or HxPhys)"
)
dev.off()

head(results)
APD = data.frame(as.character(results$phenotype),as.character(results$phenotype))
colnames(APD) = c('phenotype','description')
dim(results)
dim(APD)
colnames(APD)


pdf(file = paste('manhattan.rep.annotated', result.name, 'pdf', sep = '.'), width = 10, height = 7)

phewasManhattan(results,
                # y.axis.interval = 5, max.y = 55, 
                y.axis.interval = 10, max.y = 80, max.x = 111,
                size.x.labels = 11,
                size.y.labels = 11,
                # annotate.phenotype = F, 
                suggestive.line = .05,
                OR.direction = T, point.size = 1,
                annotate.level = 10^-100,
                annotate.list = as.list(colnames(pheno.select)),
                # annotate.phenotype = APD, annotate.angle = 0, annotate.size = 5, 
                # annotate.angle = 0, annotate.size = 5, annotate.only.largest = T,
                annotate.phenotype = APD, annotate.snp.w.phenotype = F, annotate.snp = F,
                annotate.size = 2.5, annotate.angle = 30, # annotate.only.largest = T,
                # annotate.snp.angle = 0, 
                # lc.labels = F, x.group.labels = T, x.phenotype.labels = F, sizes = F,
                # direction = F, point.size = 3, use.color = T, color.palette, 
                # title = paste0("Phenotype Plot ", date()), x.axis.label = "Phenotypes", 
                # y.axis.label = "Values", y.axis.interval = 5
                title = "Phenotype Associations with MSSS"#,
                # ICD9 count threshold = 3 (inclusive), 
                #covariates = (SEX, RACE, SMK, AGE@onset, PROGRESSIVE@onset,
                #disease.duration, short.followup, ICD9\'s from PCP or HxPhys)"
)

dev.off()


pdf(file = paste('manhattan.all.annotated', result.name, 'pdf', sep = '.'), width = 10, height = 7)

phewasManhattan(results,
                # y.axis.interval = 5, max.y = 55, 
                y.axis.interval = 10, max.y = 80, max.x = 111,
                size.x.labels = 11,
                size.y.labels = 11,
                # annotate.phenotype = F, 
                suggestive.line = .05,
                OR.direction = T, point.size = 1,
                # annotate.level = 10^-100,
                # annotate.list = as.list(colnames(pheno.select)),
                # annotate.phenotype = APD, annotate.angle = 0, annotate.size = 5, 
                # annotate.angle = 0, annotate.size = 5, annotate.only.largest = T,
                annotate.phenotype = APD, annotate.snp.w.phenotype = F, annotate.snp = F,
                annotate.size = 2.5, annotate.angle = 30, annotate.only.largest = T,
                # annotate.snp.angle = 0, 
                # lc.labels = F, x.group.labels = T, x.phenotype.labels = F, sizes = F,
                # direction = F, point.size = 3, use.color = T, color.palette, 
                # title = paste0("Phenotype Plot ", date()), x.axis.label = "Phenotypes", 
                # y.axis.label = "Values", y.axis.interval = 5
                title = "Phenotype Associations with MSSS"#,
                # ICD9 count threshold = 3 (inclusive), 
                #covariates = (SEX, RACE, SMK, AGE@onset, PROGRESSIVE@onset,
                #disease.duration, short.followup, ICD9\'s from PCP or HxPhys)"
)

dev.off()

