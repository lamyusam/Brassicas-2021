library("mice")
library("missForest")
library("tidyverse")
library("sinkr")
library("pcaMethods")

#narrow down to non-biomass data and remove 1820 data since these are mostly mutually exclusive with non-1820 data
impute = phenodata.gen2.clean.raph %>% 
  select(measure.vars) %>%
  select(-c("Aboveground_dw","Root_dw","Root_to_shoot_ratio","Days_germ")) %>%  # also remove days germ, because the distribution is so massively unimodal
  select(-grep(names(phenodata.gen2.clean.raph), pattern = "1820",value = T))

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
init = mice(impute.test.miss, m=5, method = 'pmm')
meth = init$method
predM = init$predictorMatrix
imputed = mice(impute.test.miss, method=meth, predictorMatrix=predM, m=5)
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





library(sinkr)

imputed.dineof = dineof(Xo = as.matrix(impute.test.miss), delta.rms = 1e-05)

for(i in 1:length(names(impute.test))) {
  
  actual <- impute.test[is.na(impute.test.miss[,i]),i]
  predicted <- imputed.dineof$Xa[is.na(impute.test.miss[,i]),i]
  
  foo = data.frame(actual = actual, predicted = predicted, var = names(impute.test)[i])
  
  if(i==1){
    imp.plotdata = foo
  } else {
    imp.plotdata = rbind(imp.plotdata, foo)
  }
}

dineofplot = ggplot(imp.plotdata, aes(x = actual, y = predicted)) +
  geom_point() +
  facet_wrap(~var, scales = "free")
dineofplot

ggsave(dineofplot, 
       filename = "imputation_scatters_dineof.png",
       device = "png", path = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Images/",
       width =  40, height = 25, units = "cm")

#imputed.dineof = dineof(Xo = as.matrix(impute), delta.rms = 1e-05)


imputed.methods.pca = pca(impute.test.miss, nPcs=2, method="ppca")
imputed.methods = completeObs(imputed.methods.pca)

for(i in 1:length(names(impute.test))) {
  
  actual <- impute.test[is.na(impute.test.miss[,i]),i]
  predicted <- imputed.methods[is.na(impute.test.miss[,i]),i]
  
  foo = data.frame(actual = actual, predicted = predicted, var = names(impute.test)[i])
  
  if(i==1){
    imp.plotdata = foo
  } else {
    imp.plotdata = rbind(imp.plotdata, foo)
  }
}

methodsplot = ggplot(imp.plotdata, aes(x = actual, y = predicted)) +
  geom_point() +
  facet_wrap(~var, scales = "free")
methodsplot

ggsave(methodsplot, 
       filename = "imputation_scatters_pcamethods.png",
       device = "png", path = "/home/benjamin/Dropbox/Soton Postdoc/Meeting notes/Images/",
       width =  40, height = 25, units = "cm")

md.pattern(impute, plot = F)
table(complete.cases(impute.test))
table(complete.cases(select(impute.test, -c("Sla_2326","Leaf_length_2326"))))
