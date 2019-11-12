
require(stringr)
#Format fichier chromHMM.bed 
#1 : bin / 2:chromosome / 3:start / 4:stop / 5:name / 6:score / 7:strand / 8: thickStart / 9: thickEnd / 10:itemRgb

tab1 <- read.table("/home/theophile/Documents/fichierDL/wgEncodeBroadHmmHuvecHMM.bed")
colnames(tab1) <- c("chromosome" ,"start", "stop","name",  "score", "strand", "thickStart","thickEnd","itemRgb")
tab2 <- tab1[1:4] 
tab3 <- tab2[str_detect(tab2$name,"_Strong_Enhancer"),]
write.table(tab3, "/home/theophile/Documents/huvec/Huvec_strong_enhancer.bed", row.names = FALSE)

tab3 <- tab2[str_detect(tab2$name,"_Weak_Enhancer"),]
write.table(tab3, "/home/theophile/Documents/huvec/Huvec_weak_enhancer.bed",row.names = FALSE )

tab3 <- tab2[str_detect(tab2$name,"_Active_Promoter"),]
write.table(tab3, "/home/theophile/Documents/huvec/Huvec_active_promoter.bed",row.names = FALSE)

tab3 <- tab2[str_detect(tab2$name,"_Weak_Promoter"),]
write.table(tab3, "/home/theophile/Documents/huvec/Huvec_weak_promoter.bed",row.names = FALSE)

tab3 <- tab2[str_detect(tab2$name,"_Poised_Promoter"),]
write.table(tab3, "/home/theophile/Documents/huvec/Huvec_inactiv_promoter.bed",row.names = FALSE)	