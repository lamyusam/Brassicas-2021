library(rcompanion)

###attempting to model fit between num leaves and abg dry weight for raphanus

#pearson correlation is significant
cor.test(phenodata.gen2.clean.raph$Leaf_length_1820,phenodata.gen2.clean.raph$Aboveground_dw,)

#lm is significant but fit is bad
imp.lm = glm(Aboveground_dw~Leaf_length_1820, family = gaussian(), data=phenodata.gen2.clean.raph)
summary(imp.lm)
performance::check_model(imp.lm)

#lm with log link
imp.lm.log = glm(Aboveground_dw~Leaf_length_1820, family = gaussian(link="log"), data=phenodata.gen2.clean.raph)
summary(imp.lm.log)
performance::check_model(imp.lm.log)

#glm with gamma
imp.glm.gamma = glm(Aboveground_dw~Leaf_length_1820, data=phenodata.gen2.clean.raph, family = "Gamma")
summary(imp.glm.gamma)
performance::check_model(imp.glm.gamma)

#glm with quasi
imp.glm.quasi = glm(Aboveground_dw~Leaf_length_1820, data=phenodata.gen2.clean.raph, family = "quasi", )
summary(imp.glm.quasi)
performance::check_model(imp.glm.quasi)

compareGLM(imp.lm, imp.lm.log, imp.glm, imp.glm.quasi)
