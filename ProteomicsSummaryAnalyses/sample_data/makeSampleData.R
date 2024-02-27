
library (data.table)


protQuant <- fread ("/Users/ben/Library/CloudStorage/Box-Box/Project_Coronavirus/00_data/proteomics/Level2/UK6_PH_results-mss-ProteinLevelData_v2.txt")
protQuant <- protQuant[grep("VIC|Mock|DeltaA|Omicron1", GROUP)]

results <- fread ("/Users/ben/Library/CloudStorage/Box-Box/Project_Coronavirus/00_data/proteomics/Level2/UK6_PH_results_v2.txt")
allLabels <- results$Label |> unique()
goodLabels <- grep("(VIC|Mock|DeltaA|Omicron1)_(10h|24h)-Mock_\\2", allLabels, value = TRUE)
results <- results[Label %in% goodLabels]


bigProteins <- protQuant[LogIntensities > 17, .N, by = Protein][N > 15, Protein] 

fwrite (protQuant[Protein %in% bigProteins], "Phospho_ProteinLevelData.csv")
fwrite (results[Protein %in% bigProteins], "Phospho_GroupComparisonResult.csv")
