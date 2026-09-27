# Species distribution modelling with biomod2.

library(biomod2)
sp_name <- "Abaeis_nicippe"
myRespName <- sp_name

myRespXY <- pixel_level_df %>%
  filter(species == sp_name) %>%
  dplyr::select(x, y)

myResp <- rep(1, nrow(myRespXY))

myExpl <- c(env_cont, landuse_ras)

myExpl$landuse <- as.factor(myExpl$landuse)

myBiomodData <- BIOMOD_FormatingData(
  resp.var = myResp,
  expl.var = myExpl,
  resp.xy = myRespXY,
  resp.name = myRespName,
  PA.nb.rep = 1,
  PA.nb.absences = 1000,
  PA.strategy = 'random'
)

myBiomodData
plot(myBiomodData)
