# DEG_Pattern
Mfuzz clustering for DEGs and thier patterns 
#command 
Rscript DEG_PatternsV2_Genotype_Corder_allGeneLineCOlorCodeMfuzzMember1.r -d DesignMatrix_fish2.txt -c contrast2.txt -i CountsFolder -p Data3label -a 0.05 -l 0.263 -m any

-d DesignMatrix_fish2.txt (deisgnmatrix) 
-c contrast2.txt (compare contrasts)
-i CountsFolder (path to individual counts) 
-p Data3label (label for output files) 
-a 0.05 (padj cutoff)
-l 0.263 (log2FoldChange cutoff) 
-m any (run mode)
