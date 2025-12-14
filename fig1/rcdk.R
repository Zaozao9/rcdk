## Derive chemistry ----
# extract smiles and parse to correct format
smilesOnly <- as.vector(compReady$smiles)
names(smilesOnly) <- compReady$SID
smilesParsed <- rcdk::parse.smiles(smilesOnly, omit.nulls = T)
# extract parameters (all; need to uncomment descriptors above if subset)
redundant <- c(2, 7, 8, 11, 15, 17, 18, 20, 21, 24, 29, 33:38, 41, 43:45)
descriptors <- rcdk::get.desc.names()[-redundant]
rawChem <- rcdk::eval.desc(smilesParsed, descriptors)
# format output
compOut <- compReady %>%
  left_join(., rownames_to_column(rawChem, "SID")) %>%
  select(SID, smiles, Fsp3:nAcid)
# make rownames
rownames(compOut) <- compOut$SID