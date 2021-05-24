library(rcompanion)

phenodata.gen2.clean.brass.trans = phenodata.gen2.clean.brass %>%
  mutate(Days_germ = transformTukey(Days_germ))

response_var = measure.vars[1]
response_var
foo = lmer(scale(get(response_var))~Environment + Wild_Dom + Environment:Wild_Dom + (1|Wild_Dom/Population), 
           data=phenodata.gen2.clean.brass.trans, 
           control = lmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))
check_model(foo)




glmer(Aboveground_dw~Environment + Wild_Dom + Environment:Wild_Dom + (1|Wild_Dom/Population), data=brp_stand_off,family="Gamma",glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 10000)))