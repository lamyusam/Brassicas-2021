library("mice")
library("missForest")

impute.test = phenodata.gen2.clean.raph %>% 
  select(measure.vars) %>%
  select(-c("Aboveground_dw","Root_dw","Root_to_shoot_ratio")) %>% 
  select(-grep(names(phenodata.gen2.clean.raph), pattern = "1820",value = T))

impute.test %>% complete.cases() %>% table()

impute.test = impute.test[complete.cases(impute.test),]

impute.test.miss = prodNA(impute.test, noNA = 0.2)

init = mice(impute.test.miss, maxit=0) 
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

ggplot(imp.plotdata, aes(x = actual, y = predicted)) +
  geom_point() +
  facet_wrap(~var, scales = "free")
