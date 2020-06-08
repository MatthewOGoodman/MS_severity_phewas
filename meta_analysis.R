# install.packages('metafor')
library('metafor')
Date.of.analysis = "2018-03-03"


setwd('C:\\Users\\MattyOrlando\\Desktop\\Dropbox\\De Jager Lab\\2014-2. Zongqi Xia\\Figures and Results 6.0')
resultsBC = read.csv('PheWAS output_July 25th.csv')
col.classes = sapply(resultsBC,class)
col.classes["phewas_code"] = "character"
resultsBC.all = read.csv('PheWAS output_July 25th.csv', colClasses = col.classes)
head(resultsBC.all)
min(resultsBC.all$n_cases)
resultsBC = resultsBC.all[resultsBC.all$n_cases>=.03*4867 & 
                            resultsBC.all$n_controls>=.03*4867 & 
                            !is.na(resultsBC.all$n_cases) &
                            !is.na(resultsBC.all$n_controls), ]
resultsBC$phecode = resultsBC$phewas_code
resultsBC$description = resultsBC$phewas_description
resultsBC$OR_lo = exp(resultsBC$beta + qnorm(0.025,mean = 0, sd = resultsBC$SE))
resultsBC$OR_hi = exp(resultsBC$beta + qnorm(0.975,mean = 0, sd = resultsBC$SE))
resultsBC$prev.pct = round(resultsBC$n_cases/resultsBC$n_total*100,1)


setwd('C:\\Users\\MattyOrlando\\Desktop\\Dropbox\\De Jager Lab\\2014-2. Zongqi Xia\\Figures and Results 6.0')

resultsMA.all = read.csv(  paste(Date.of.analysis,'MSSS.lcl.min.rec.25.no.excl.csv', sep = '.') )
col.classes = sapply(resultsMA.all,class)
col.classes[1:2] = 'character'
resultsMA.all = read.csv( paste(Date.of.analysis,'MSSS.lcl.min.rec.25.no.excl.csv', sep = '.'),
                          colClasses = col.classes)
head(resultsMA.all)
resultsMA = resultsMA.all[resultsMA.all$n_cases>=.03* 3439 & 
                            resultsMA.all$n_controls>=.03* 3439 & 
                            !is.na(resultsMA.all$n_cases) &
                            !is.na(resultsMA.all$n_controls), ]
resultsMA$OR_lo = exp(resultsMA$beta + qnorm(0.025,mean = 0, sd = resultsMA$SE))
resultsMA$OR_hi = exp(resultsMA$beta + qnorm(0.975,mean = 0, sd = resultsMA$SE))
resultsMA$prev.pct = round(resultsMA$n_cases/resultsMA$n_total*100,1)


both.phewas.codes = intersect(resultsMA$phecode, resultsBC$phecode)
resultsMA$phecode[!(resultsMA$phecode %in% resultsBC$phecode)]
resultsBC$phecode[!(resultsBC$phecode %in% resultsMA$phecode)]
resultsMA$prev.pct = round(resultsMA$n_cases/resultsMA$n_total*100,1)

names(resultsMA)
names(resultsBC)
code = both.phewas.codes[1]
r = length(both.phewas.codes)

meta.result = data.frame(logOR = NA, OR = NA, SE = NA, zval = NA, p = NA)
for (i in 1:r){
  code = both.phewas.codes[i]
  data = rbind(resultsMA[resultsMA$phecode == code,c('phecode','beta','SE')],
               resultsBC[resultsBC$phecode == code,c('phecode','beta','SE')]
               )
  # model = rma(yi = log(OR), vi = SE^2*(1/log(OR))^2, method = 'FE', data = data)
  model = rma(yi = beta, sei = SE, method = 'FE', data = data)
  names(model)
  meta.result[i,] = c(model$b,exp(model$b),model$se,model$zval,model$pval)
  # rma(yi = log(OR), vi = SE, method = 'FE', data = data)
}

meta.result$OR_lo = exp(meta.result$logOR + qnorm(0.025,mean = 0, sd = meta.result$SE))
meta.result$OR_hi = exp(meta.result$logOR + qnorm(0.975,mean = 0, sd = meta.result$SE))

colnames(meta.result) = paste(colnames(meta.result), 'meta', sep = '.')
meta.result = data.frame(phecode = both.phewas.codes,meta.result)

sapply(as.data.frame(meta.result),class)

MAinBC.TF = resultsMA$phecode %in% resultsBC$phecode
BCinMA.TF = resultsBC$phecode %in% resultsMA$phecode

meta.df = merge(meta.result, 
               merge(resultsMA[MAinBC.TF, c('phecode', 'description', 'beta', 'OR', 'OR_lo', 'OR_hi', 'SE', 'p', 
                                            'n_total', 'n_cases', 'n_controls', 'prev.pct') ],
                     resultsBC[BCinMA.TF, c('phecode', 'description', 'beta', 'OR', 'OR_lo', 'OR_hi', 'SE', 'p', 
                                            'n_total', 'n_cases', 'n_controls', 'prev.pct')], 
                     by = 'phecode', suffixes = c('.MA','.BC')),
               by = 'phecode'
               )
meta.df$prev.pct.meta = round( (meta.df$n_cases.MA + meta.df$n_cases.BC)/ (meta.df$n_total.MA + meta.df$n_total.BC)*100,1)

output.columns =c('phecode', 'description.MA', 
                  paste(c('prev.pct', 'OR', 'OR_lo', 'OR_hi', 'p'), 'meta', sep = '.' ),
                  paste(c('prev.pct', 'OR', 'OR_lo', 'OR_hi', 'p'), 'MA', sep = '.' ),
                  paste(c('prev.pct', 'OR', 'OR_lo', 'OR_hi', 'p'), 'BC', sep = '.' ))

(meta.df[order(meta.df$p.meta), output.columns ])
cor(meta.df[,grep('OR\\.', colnames(meta.df))])
cor(meta.df[,grep('logOR|beta', colnames(meta.df))])

sapply(meta.df,class)

#Bonferroni significant:
meta.df$BONF.indep.meta = meta.df$p.meta < .05/nrow(meta.df)

meta.df[ meta.df$BONF.indep.meta, c( 'description.MA', 'prev.pct.meta', 'OR.meta', 'OR_lo.meta', 'OR_hi.meta', 'p.meta')]
sum(meta.df$BONF.indep.meta)
dim(meta.df)

write.csv(meta.df[,c(output.columns, 'BONF.indep.meta')], 
          paste(Sys.Date(),'Meta.analysis.Fixed.effects.Complete.list.(TableS3).csv', sep = '.'), row.names = FALSE)


meta.df[((meta.df$OR.meta > 1) != (meta.df$OR.MA > 1)) & meta.df$BONF.indep.meta,'description.MA']
meta.df[((meta.df$OR.meta > 1) != (meta.df$OR.MA > 1)) & meta.df$BONF.indep.meta,c('description.MA', 'OR.meta', 'OR.MA')]

