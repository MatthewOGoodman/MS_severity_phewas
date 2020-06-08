analysis = 'primary'
# analysis = 'secondary'
Date.of.analysis = "2018-03-03"

if(analysis == 'primary'){
  BC.raw = read.csv('PheWAS output_July 25th.csv',
                    header = TRUE , colClasses = c("phewas_code" = "character", "phewas_description" = "character"))
  BC.raw$phecode = BC.raw$phewas_code
  BC.raw$description = BC.raw$phewas_description
  BC.case.threshold = .03*4867 ## 4,867: Note that in the paper the number reported is 4,876
  
  MA.raw = read.csv(  paste(Date.of.analysis,'MSSS.lcl.min.rec.25.no.excl.csv', sep = '.'),
                      header = TRUE , colClasses = c("phecode" = "character", "description" = "character"))
  MA.case.threshold = .03* 3439
  
  # BC.raw.primary = BC.raw
  # MA.raw.primary = MA.raw

  } else if(analysis == 'secondary'){
  BC.raw = read.csv('phewas secondary analysis - output.csv',
                    header = TRUE , colClasses = c("phenotype" = "character", "description" = "character"))
  BC.raw$phecode = BC.raw$phenotype
  BC.case.threshold =  .03* 3439 #using same threshold as above for MA, can't use 3% because lower sample size.
  
  MA.raw = read.csv(  paste(Date.of.analysis,'edss.6.at.7.min.rec.15.no.excl.csv', sep = '.'),
                      header = TRUE , colClasses = c("phecode" = "character", "description" = "character"))
  MA.case.threshold =  .03* 3439 #using same threshold as above for MA, can't use BC threshold since too stringent
  MA.raw = MA.raw[,!(colnames(MA.raw) == 'X')]
  

}

names(BC.raw); names(MA.raw)
nrow(BC.raw); nrow(MA.raw)

# Subset MA:
MA.raw = MA.raw[!is.na(MA.raw$p),] 
nrow(MA.raw)
MA = MA.raw[MA.raw$n_cases >  MA.case.threshold, ]
max(MA$n_total)
min(MA$n_cases)
dim(MA) #[1] 109  22 (analyzed phecodes)
MA$prev.pct = round(MA$n_cases/MA$n_total*100,1)
MA$OR_lo = exp(MA$beta + qnorm(0.025,mean = 0, sd = MA$SE))
MA$OR_hi = exp(MA$beta + qnorm(0.975,mean = 0, sd = MA$SE))
MA$BONF.indep = MA$p < 0.05 / nrow(MA)
write.csv(MA[, c("phecode", "description", 
               "prev.pct", "OR", "OR_lo", "OR_hi", "p", "BONF.indep")],
          file = paste(Sys.Date(),analysis,'MA.complete.list.(TableS1v2).csv', sep = '.'),
          row.names = FALSE)

# Subset BC:
BC.raw = BC.raw[!is.na(BC.raw$p), ]
nrow(BC.raw) #601
BC = BC.raw[BC.raw$n_cases > BC.case.threshold,  ]
max(BC$n_total)
min(BC$n_cases)
dim(BC) #[1] 365  25 (analyzed phecodes)
BC$prev.pct =  round(BC$n_cases/BC$n_total*100,1)
BC$OR_lo = exp(BC$beta + qnorm(0.025,mean = 0, sd = BC$SE))
BC$OR_hi = exp(BC$beta + qnorm(0.975,mean = 0, sd = BC$SE))
BC$BONF.indep = BC$p < 0.05 / nrow(BC)

write.csv(BC[, c("phecode", "description", 
               "prev.pct", "OR", "OR_lo", "OR_hi", "p", "BONF.indep")],
          file = paste(Sys.Date(),analysis,'BC.complete.list.(TableS2v2).csv', sep = '.'),
          row.names = FALSE)

# Find intersection:
overlap = intersect(MA$phecode, BC$phecode)

# MA$description[!MA$phecode %in% overlap]
MA.not.analyzed = MA[!MA$phecode %in% overlap,]
nrow(MA.not.analyzed) #[1] 31
write.csv(MA.not.analyzed, file = paste(Sys.Date(),analysis,'MA.not.analyzed.in.replication.csv', sep = '.'))

# BC$description[!BC$phecode %in% overlap]
BC.not.analyzed = BC[!BC$phecode %in% overlap,]
nrow(BC.not.analyzed) #[1] 287
write.csv(BC.not.analyzed, file = paste(Sys.Date(),analysis,'BC.not.analyzed.in.replication.csv', sep = '.'))


# Subset to intersection of analyzed codes:
MA.rep = MA[MA$phecode %in% overlap, ]
BC.rep = BC[BC$phecode %in% overlap, ]
dim(MA.rep);dim(BC.rep)



#add bonferroni for replication:
#MA$BONF.rep:
MA.rep$BONF.rep = MA.rep$p < .05/nrow(MA.rep)
nrow(MA.rep)
MA.raw[MA.raw$p < 0.05/78 & MA.raw$p > 0.05/381 ,]

#MA significant results
MA.significant = MA.rep[MA.rep$BONF.rep,]
N.MA.rslts = nrow(MA.significant)

#BC$BONF.rep:
min(BC.rep$n_cases)
BC.rep$BONF.rep = BC.rep$p < .05/nrow(BC.rep)
all(BC.rep$BONF.rep  == BC.rep$bonferroni)
#BC significant results
BC.significant = BC.rep[BC.rep$BONF.rep ,]
N.BC.rslts = nrow(BC.significant)

MA.BC.rep = merge(MA.rep,BC.rep, by = 'phecode', suffixes = c('.MA','.BC'))
MA.BC.meta.rep = 
  merge(MA.BC.rep[,c("phecode", "description.MA", 
                                    "prev.pct.MA", "OR.MA", "OR_lo.MA", "OR_hi.MA", "p.MA",
                                    "prev.pct.BC", "OR.BC", "OR_lo.BC", "OR_hi.BC", "p.BC")],
                       meta.df[,c("phecode", "prev.pct.meta", "OR.meta", "OR_lo.meta", "OR_hi.meta", "p.meta")], by = "phecode" )

write.csv(MA.BC.meta.rep,
  file = paste(Sys.Date(),analysis,'MA.BC.all.analyzed.in.replication.(Table.S1,S2,S3).csv', sep = '.'),
  row.names = FALSE) #all replication results in overlap with BONF.rep


# The next two steps are only relevant in cases where the MA and BC analyzed in replication are not identical
# i.e. if we don't first reduce to the intersection of the analyzed codes.

#BC codes that are available for replication
BC.for.rep = BC.rep[BC.rep$phecode %in% MA.significant$phecode,]
N.BC.for.rep = nrow(BC.for.rep)
N.BC.for.rep

#Combined data frame of significant MA results with BC results for comparison
MA.sig.BC.for.rep = merge(MA.significant,BC.for.rep, by.x = 'phecode', by.y = 'phecode', suffixes = c('.MA','.BC'))
nrow(MA.sig.BC.for.rep)

#Logical Test: significant replication p-value?
BC.rep.sig = MA.sig.BC.for.rep$p.BC < 0.05/nrow(BC.rep) 

#Logical Test: consistent replication direction?
BC.rep.dir = (MA.sig.BC.for.rep$beta.MA > 0) == (MA.sig.BC.for.rep$beta.BC > 0)

sum(BC.rep.sig)
sum(BC.rep.dir)
sum(BC.rep.sig & BC.rep.dir)


BC.rep.dir.sig = MA.sig.BC.for.rep[BC.rep.sig & BC.rep.dir,]
BC.rep.dir.sig[,c("phecode", "description.MA", "prev.pct.MA", "OR.MA", "p.MA", "prev.pct.BC", "OR.BC", "p.BC")]
nrow(BC.rep.dir.sig)
write.csv(BC.rep.dir.sig[,
                         c("phecode", "description.MA", 
                           "prev.pct.MA", "OR.MA", "OR_lo.MA", "OR_hi.MA", "p.MA",
                           "prev.pct.BC", "OR.BC", "OR_lo.BC", "OR_hi.BC", "p.BC")],
          file = paste(Sys.Date(),analysis,'BC.replicated.significant.csv', sep = '.'),
          row.names = FALSE)


if(analysis == "primary"){
  Table2 = merge( meta.df[,output.columns],  BC.rep.dir.sig[, "phecode", drop = FALSE
                                                            # c("phecode", "OR.MA", "OR.BC") # to check with meta.df
                                                            ],by = 'phecode')
  write.csv(Table2, file = paste(Sys.Date(),analysis,'BC.replicated.significant.w.meta.analysis.(Table2).csv', sep = '.'),
            row.names = FALSE)
}


BC.not.sig = MA.sig.BC.for.rep[!BC.rep.sig & BC.rep.dir,]
nrow(BC.not.sig)
write.csv(BC.not.sig, file = paste(Sys.Date(),analysis,'BC.not.replicated(not.significant).csv', sep = '.'),
          row.names = FALSE)

BC.not.consistent = MA.sig.BC.for.rep[!BC.rep.dir,]
nrow(BC.not.consistent)
write.csv(BC.not.consistent, file = paste(Sys.Date(),analysis,'BC.not.replicated(inconsistent.effect).csv', sep = '.'),
          row.names = FALSE)


subset.col = c('phecode','n_cases','OR','p')
merge(MA.significant[MA.significant$phecode %in% BC.rep.dir.sig$phecode,       c('description', subset.col)],
      BC            [BC$phecode             %in% BC.rep.dir.sig$phecode,       subset.col],            
      by.x = 'phecode', by.y = 'phecode', suffixes = c('.MA','.BC'))

plot(MA.sig.BC.for.rep$beta.MA,MA.sig.BC.for.rep$beta.BC)
abline(0,1)
cor(MA.sig.BC.for.rep$beta.MA, MA.sig.BC.for.rep$beta.BC)

