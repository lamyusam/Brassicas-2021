library("mice")
library("missForest")
library("tidyverse")

setwd("/home/benjamin/Documents/Brassicas_repo")
#source("Analysis/Phenotypic_analysis/morphology_modeltesting_september.R")

#create object for mice imputation (removing root shoot ratio since this is not a direct measurement)
#phenodata.gen2.clean.brass.test = select(phenodata.gen2.clean.brass, measure.vars.brass) %>% select(-c("Log_root_shoot_ratio"))

#first select predictors- we don't want to impute between variables unless they're correlated (in this case by at least 0.5)
imp = mice(phenodata.gen2.clean.brass, pred=quickpred(phenodata.gen2.clean.brass, mincor=.5), print=F, seed=123)
imp$pred
#some traits don't have a high enough correlation with any other trait and so won't be imputed: 
#Height & Leaf Length 2326, Num Leaves 1820, and SLA & SLDW 1820 (and of course the non-pheno variables)
rowSums(imp$pred)
#For the remaining variables, predictive mean matching will be used, which is appropriate since they're all numeric:
imp$meth
#Check whether imputed values are plausible (blue=original values, red=imputed values)
stripplot(imp, as.formula(paste0(paste(measure.vars.brass,collapse = "+"),"~.imp")),  pch=20, cex=2)
#yes, imputed values look like real values (not surprising given that we're using PMM!)














#narrow down to non-biomass data and remove 1820 data since these are mostly mutually exclusive with non-1820 data
impute = phenodata.gen2.clean.raph %>% 
  select(measure.vars) %>%
  select(-c("Aboveground_dw","Root_dw","Root_to_shoot_ratio","Days_germ")) %>%  # also remove days germ, because the distribution is so massively unimodal
  select(-grep(names(phenodata.gen2.clean.raph), pattern = "1820",value = T)) %>%
  mutate_all(scale) 
  

#display missing values: 
#left side shows number of rows with that combination missing
#bottom shows total number missing for each variable
md.pattern(impute, plot = F)
#also check number of complete cases manually
impute %>% complete.cases() %>% table()
#remove incomplete cases to get testing dataset
impute.test = impute[complete.cases(impute),]
#introduce some NAs
set.seed(123)
impute.test.miss = prodNA(impute.test, noNA = 0.1)

#create mice output the simple way
init = mice(as.matrix(impute.test.miss), m=5, method = 'pmm')
meth = init$method
predM = init$predictorMatrix
imputed = mice(as.matrix(impute.test.miss), method=meth, predictorMatrix=predM, m=5)
imputed = complete(imputed)

for(i in 1:length(names(impute.test))) {
  
  actual <- impute.test[is.na(impute.test.miss[,i]),i]
  predicted <- imputed[is.na(impute.test.miss[,i]),i]
  
  foo = data.frame(actual = actual, predicted = predicted, var = names(impute.test)[i])
  
  if(i==1){
    imp.plotdata = foo
  } else {
    imp.plotdata = rbind(imp.plotdata, foo)
    }
}

miceplot = ggplot(imp.plotdata, aes(x = actual, y = predicted)) +
  geom_point() +
  facet_wrap(~var, scales = "free")
miceplot

ggsave(miceplot, 
       filename = "imputation_scatters_mice.png",
       device = "png", path = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Images/",
       width =  40, height = 25, units = "cm")

# #create mice output a more complicated way
# # Build regression model
# model_fit <- with(data = imputed, exp = lm(Days_germ ~ Height_2326 + Leaf_length_2326)) 
# 
# # combining results of all 5 models using pool() function
# pooled_output <- pool(model_fit)
# summary(pooled_output)

