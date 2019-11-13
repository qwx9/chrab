H3K4me3 <- read.table("/home/theophile/Documents/Projet/chrab/H3K4me3.bed")
H3K4me32<- H3K4me3[1:5]
colnames(H3K4me32) <- c("chromosome" ,"start", "stop","name",  "score")
write.table(H3K4me32, "/home/theophile/Documents/Projet/chrab/H3K4me3_count.bed",row.names = FALSE)

H3K27ac <- read.table("/home/theophile/Documents/Projet/chrab/H3K27ac.bed")
H3K27ac2 <- H3K27ac[1:5]
colnames(H3K4me32) <- c("chromosome" ,"start", "stop","name",  "score")
write.table(H3K27ac2, "/home/theophile/Documents/Projet/chrab/H3K27ac_count.bed",row.names = FALSE)

DHS <- read.table("/home/theophile/Documents/Projet/chrab/DHS.bed")
DHS2 <- DHS[1:5]
colnames(DHS2) <- c("chromosome" ,"start", "stop","name",  "score")
write.table(DHS2, "/home/theophile/Documents/Projet/chrab/DHS_count.bed",row.names = FALSE)