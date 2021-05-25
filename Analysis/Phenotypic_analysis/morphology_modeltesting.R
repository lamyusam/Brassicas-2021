library(rcompanion)

phenodata.gen2.clean.brass.trans = phenodata.gen2.clean.brass %>%
  mutate(Root_dw = transformTukey(Root_dw))

response_var = measure.vars[8]
response_var
ggqqplot(phenodata.gen2.clean.brass.trans$Root_dw)

foo = lmer(scale(get(response_var))~Environment + Wild_Dom + Environment:Wild_Dom + (1|Wild_Dom/Population), 
           data=phenodata.gen2.clean.brass.trans, 
           control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
check_model(foo)




goo = glmer(get(response_var)~Environment + Wild_Dom + Environment:Wild_Dom + (1|Wild_Dom/Population),
      data=phenodata.gen2.clean.brass,family="Gamma",glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))

summary(goo)
