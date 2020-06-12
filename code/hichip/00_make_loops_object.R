library(diffloop)

# Import / process loops data
loopsRaw <- loopsMake.mango("../../data/hichipper", mergegap = 0) # takes about 1 minute; need to gunzip the raw loop data
loopsRaw.mango <- mangoCorrection(loopsRaw) # takes ~3 minutes
saveRDS(loopsRaw.mango, "../../output/11JUNE2020_hichip_loopsRaw.mango.rds")

