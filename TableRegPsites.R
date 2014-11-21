# Raw data:
tab <- read.csv("SILAC1018TotV5_FAAposition.csv", as.is = T)

# Packages:
library("ggplot2")
library("dplyr")

# Function:
FilledTable <- function(List, data.frame){
        mat <- matrix("", nrow=length(List), ncol=max(sapply(List, length)))
        for(x in 1:length(List))
                if(length(List[[x]]))
                {mat[x,1:length(List[[x]])] <- List[[x]]}
        cbind(data.frame, mat)
}
replace.vector <-function(x, tochange=unique(x), toreplace=1:length(unique(x))){
        for (y in 1:length(tochange))  x=replace(x, list=(x==tochange[y]), toreplace[y])
        x
}

##-------
# Retrieve sample information:

sample <- read.csv("Samples.csv")
tab <- read.csv("SILAC1018TotV5_FAAposition.csv")
tab <- tab[,c(1:ncol(tab)-1)] # remove the last column that contains only "NA".

m <- match(tab$SpectrumFile,sample$SpectrumFile)
df1 <- data.frame(tab, "Experiment"=sample$Experiment[m], "SampleName"=sample$SampleName[m], "Fraction"=sample$Fraction[m], "Cell"=sample$LabelledCell[m])
tab <- df1
tab <- tab[tab$qValue<=.01,]
tab <- tab[!is.na(tab$pRSSiteProbabilities),]
tab <- tab[tab$pRSSiteProbabilities!="",]
val1 <- length(unique(tab$Sequence))/length(unique(df1$Sequence))*100
Ambiguous <- tab[grepl(";", tab$ProteinGroupAccessions),]
UAmbSeq <- length(unique(Ambiguous$Sequence))
# write.csv(Ambiguous, "Ambiguous.csv", row.names=F) # if you want to have a 
# check at the ambiguous peptides (not attributable to a unique protein)
tab$ProteinGroupAccessions <- as.character(tab$ProteinGroupAccessions)
tab$ProteinGroupAccessions[(grepl("CLIP1_HUMAN", as.character(tab$ProteinDescriptions)))&(grepl("CLIP2_HUMAN", as.character(tab$ProteinDescriptions)))] <- "P30622"
tab$ProteinGroupAccessions[(grepl("TMCC2_HUMAN", as.character(tab$ProteinDescriptions)))&(grepl("TMCC3_HUMAN", as.character(tab$ProteinDescriptions)))] <- "O75069"
tab$FirstAminoAcidPosition <- as.character(tab$FirstAminoAcidPosition)
tab$FirstAminoAcidPosition[(grepl("CLIP1_HUMAN", as.character(tab$ProteinDescriptions)))&(grepl("CLIP2_HUMAN", as.character(tab$ProteinDescriptions)))] <- "346"
tab$FirstAminoAcidPosition[(grepl("TMCC2_HUMAN", as.character(tab$ProteinDescriptions)))&(grepl("TMCC3_HUMAN", as.character(tab$ProteinDescriptions)))] <- "436"
# Manual assessment of the "Ambiguous" table draw our attention of these two phosphopeptides. There log2(H/M) values may indicate phosphoregulation. Thus, I label them differently before I remove the ambiguous peptides from the analysis. These ones will be kept. This way, the localisation of the phosphosite will be based only on CLIP1 and TMCC2 sequences, respectively.
tab <- tab[!grepl(";", tab$ProteinGroupAccessions),]

##------
# Phosphosites ID:

prob <- as.character(tab$pRSSiteProbabilities)
P <- NULL
fp <- as.character(tab$FirstAminoAcidPosition)
for (i in 1:length(prob)){
        x <- strsplit(prob[i], split='; ', fix=T)[[1]]
        mat <- matrix(0, ncol=3,nrow=length(x))
        for (j in 1:length(x)){
                y <- strsplit(x[j], split=':', fix=T)[[1]]
                mat[j,3] <- y[2]
                y <- strsplit(y[1], split='(', fix=T)[[1]]
                y[2] <- strsplit(y[2], split=')', fix=T)[[1]][1]
                mat[j, 1:2] <- y
                mat[j,2] <- as.numeric(mat[j,2])+as.numeric(fp[i])-1
        }
        mat <- mat[as.numeric(mat[,3])>=60,,drop=F]
        mat <- mat[sort.list(as.numeric(mat[,3]),decreasing=T), ]
        mat <- as.vector(t(mat))
        P <- c(P, list(mat))
}
Uniprot <- as.character(tab$ProteinDescriptions)
Uniprot <- sapply(Uniprot, function(x) strsplit(x, split='[', fixed=T)[[1]][2])
Uniprot <- sapply(Uniprot, function(x) strsplit(x, split='_', fixed=T)[[1]][1])
names(Uniprot) <- NULL
Uniprot[grepl("CLIP2", tab$ProteinDescriptions)] <- "CLIP1/2"
Uniprot[(grepl("TMCC2", tab$ProteinDescriptions))&(grepl("TMCC3", tab$ProteinDescriptions))] <- "TMCC2/3"
# This way, the ambiguous proteins are considered as one. 
tab <- FilledTable(P, data.frame(tab, nSite=sapply(P,function(x) length(x)/3), Uniprot=Uniprot))
tab <- tab[(tab$nSite!=0), ] # Removes the peptides with ambiguous phosphosite position (pRS probability <= 60%).
ID <- vector(length=nrow(tab))
tab$Uniprot <- as.character(tab$Uniprot)
tab$"1" <- as.character(tab$"1")
tab$"2" <- as.numeric(as.character(tab$"2"))
tab$"4" <- as.character(tab$"4")
tab$"5" <- as.numeric(as.character(tab$"5"))
tab$"7" <- as.character(tab$"7")
tab$"8" <- as.numeric(as.character(tab$"8"))
for (i in 1:nrow(tab)){
        if (tab$nSite[i]%in%2) {ID[i] <- c(paste0(tab$Uniprot[i],'-',tab$"1"[i],tab$"2"[i],'-',tab$"4"[i],tab$"5"[i]))}
        if (tab$nSite[i]%in%3) {ID[i] <- c(paste0(tab$Uniprot[i],'-',tab$"1"[i],tab$"2"[i],'-',tab$"4"[i],tab$"5"[i],'-',tab$"7"[i],tab$"8"[i]))}
        else {ID[i] <- c(paste0(tab$Uniprot[i],'-',tab$"1"[i],tab$"2"[i]))}
}
tab <- data.frame(tab, ID=ID, stringsAsFactors = FALSE)

##-----
# Then I perform the same quantification as I have done with the total analysis:

# Split HUVEC and MDA:

tab <- tab[!is.na(tab$log2HM),]
tMDA <- tab[tab$Cell=="MDA-MB-231",]
tHUV <- tab[tab$Cell=="HUVEC",] 

##-------
# Filter out the PSMs from proteins that are not in the MDA or HUV regulated list:

MDAreg <- read.csv("/Users/Yoko/Documents/CRUKMI/PaperHUVMDA/PaperStats/MDARegProt.csv", as.is = T)
MDAreg <- MDAreg[,1]
HUVreg <- read.csv("/Users/Yoko/Documents/CRUKMI/PaperHUVMDA/PaperStats/HUVRegProt.csv", as.is = T)
HUVreg <- HUVreg[,1]
MDAreg <- sapply(1:length(MDAreg), function(x) strsplit(as.character(MDAreg)[x], "_", fixed = T)[[1]][1])
HUVreg <- sapply(1:length(HUVreg), function(x) strsplit(as.character(HUVreg)[x], "_", fixed = T)[[1]][1])
tMDA <- tMDA[tMDA$Uniprot %in% MDAreg,]
tHUV <- tHUV[tHUV$Uniprot %in% HUVreg,]

###================
# Checking the spectra: I have already some files with checked spectra. I need 
# to make sure that I keep only the "checked" regulation fold value for all 
# targets I talk about, not just the ones that are in our "high confidence"
# list.

spMDA <- read.csv("MDAsignV5spF.csv", as.is = T)
spMDAFLNA <- read.csv("MDAsignV5spFLMNA.csv", as.is = T)
spMDA <- select(spMDA, X, Sequence, ProteinDescriptions, SpectrumFile, log2HM, normalisedratioHM, NormalisationFactor, UniqueID, ID, Experiment, IDQuan)
spMDAFLNA <- select(spMDAFLNA, X, Sequence, ProteinDescriptions, SpectrumFile, log2HM, normalisedratioHM, NormalisationFactor, UniqueID, ID, Experiment, IDQuan)
spMDA <- rbind(spMDA, spMDAFLNA)
Uniprot <- sapply(1:nrow(spMDA), function(x) strsplit(spMDA$ID[x], "-")[[1]][1])
spMDA <- data.frame(spMDA, "Uniprot" = Uniprot, stringsAsFactors = F)
# proteins I have to check the spectra of for MDAs:
MDAprot <- setdiff(unique(tMDA$Uniprot), unique(spMDA$Uniprot))

## Same for the HUVEC:
spHUV <- read.csv("HUVsignV5spF.csv", as.is = T)
spHUV <- select(spHUV, X, Sequence, ProteinDescriptions, SpectrumFile, log2HM, normalisedratioHM, NormalisationFactor, UniqueID, ID, Experiment, IDQuan)
Uniprot <- sapply(1:nrow(spHUV), function(x) strsplit(spHUV$ID[x], "-")[[1]][1])
spHUV <- data.frame(spHUV, "Uniprot" = Uniprot, stringsAsFactors = F)
# proteins I have to check the spectra of for HUVs:
HUVprot <- setdiff(unique(tHUV$Uniprot), unique(spHUV$Uniprot))

###----------------
# Export the list of spectra to check for each cell line:

tHUVtoCheck2 <- tHUV[tHUV$Uniprot %in% HUVprot,]
write.csv(tHUVtoCheck2, "Tables/HUVtoCheck.csv", row.names = F)

tMDAtoCheck2 <- tMDA[tMDA$Uniprot %in% MDAprot,]
write.csv(tMDAtoCheck2, "Tables/MDAtoCheck.csv", row.names = F)


##=======
# HUVEC:

tab <- tHUV
spchecked <- read.csv("Tables/HUVChecked.csv", as.is = T)
spchecked <- unique(spchecked$UniqueID[spchecked$X == "g"])
spchecked1 <- read.csv("HUVsignV5spF.csv", as.is = T)
spchecked1 <- unique(spchecked1$UniqueID[spchecked1$X == "g"])
spchecked <- c(spchecked, spchecked1)
tab <- tab[tab$UniqueID %in% spchecked,]
tab <- tab[!is.na(tab$log2HM),]
IDQuan <- paste0(tab$QuanResultID,"-",tab$SpectrumFile) # creates a unique quantification ID per phosphosite. 
tab <- data.frame(tab,IDQuan=IDQuan,stringsAsFactors=F)

QID <- lapply(unique(tab$IDQuan), function(y) tab[tab$IDQuan==y, 'X3']) 
names(QID) <- unique(tab$IDQuan) # list with the pRS probability of the first phosphosite for each QuanResultID
SID <- lapply(unique(tab$IDQuan), function(y) tab[tab$IDQuan==y, 'SpectrumFile']) 
names(SID) <- unique(tab$IDQuan) # list with the sample files for each QuanResultID
UID <- lapply(unique(tab$IDQuan), function(y) tab[tab$IDQuan==y, 'UniqueID']) 
names(UID) <- unique(tab$IDQuan) # list with the unique row ID for each QuanResultID


MaxpRS <- lapply(1:length(QID), function(x) aggregate(as.numeric(as.character(QID[[x]])), list(as.character(SID[[x]])), max)[,2]) 
names(MaxpRS) <- unique(tab$IDQuan) # determines which quantification result corresponds to the highest pRS probability for a phosphosite in a given sample.

v <- c(rep(0,length(MaxpRS)))
for (i in 1:length(MaxpRS)) {
        for (j in 1:length(UID[[i]])) {
                if (MaxpRS[[i]]==as.numeric(as.character(QID[[i]][j]))) {
                        v[i] <- as.character(UID[[i]][j])
                }
        }
}

tab1 <- tab[tab$UniqueID%in%v,] ## This table contain one quan per peptide + the peptide of this quanvalue with the max pRS

#### Then, I perform quantification (average of all median per sample) on the phosphosites 

tab1$normalisedratioHM <- as.numeric(as.character(tab1$normalisedratioHM))
tab1$Experiment <- as.character(tab1$Experiment)
dat <- lapply(unique(tab1$ID), function(y) tab1[tab1$ID==y, 'normalisedratioHM'])
dat2 <- lapply(unique(tab1$ID), function(y) tab1[tab1$ID==y, 'Experiment'])
names(dat) <- unique(tab1$ID)



dat <- lapply(1:length(dat), function(x) aggregate(dat[[x]], list(as.character(dat2[[x]])), median)[,2]) # Keep the median of log2(H/M) per biological replicate in order to perform the t-test.
names(dat) <- unique(tab1$ID)

# s <- lapply(dat, length)
# dat <- dat[s>=2] # Here I keep even the psites quantified in only one experiment
# p5 <- sapply(1:length(dat), function(y) t.test(dat[[y]], alternative='two.sided')$p.value)
fc <- sapply(1:length(dat), function(y) mean(dat[[y]]))
Uniprot <- tab1$Uniprot[match(names(dat), tab1$ID)]


# Create matrix with median value per experiment for each phosphosite in the list:
mat <- matrix(ncol = (length(unique(tab1$Experiment)) + 1), nrow = length(dat))
mat[,1] <- names(dat)
colnames(mat) <- c("ID", unique(as.character(tab1$Experiment)))
for (i in 1:nrow(mat)) {
        dat2[[i]] <- unique(dat2[[i]])
        for (j in 1:length(dat[[i]])) {
                for (k in 2:(length(unique(tab1$Experiment)) + 1)) {
                        if (colnames(mat)[k] == dat2[[i]][j]) mat[i,k] <- dat[[i]][j]
                }
        }
}

tabHUV <- data.frame(mat, "FC" = fc, "Uniprot" = Uniprot)
write.csv(tabHUV, "Tables/TopHitsHUVECAverages.csv", row.names = F)

# Now, I want to have a figure with the values: I plot "normalisedlog2HM" from tab1 (filtered table with only checked spectra):

g <- ggplot(data = tab1, aes(x = ID, y = normalisedratioHM, shape = Experiment, col = Fraction)) + geom_hline(yintercept=0,) + geom_hline(yintercept=-1,linetype="dashed",col="grey") + geom_hline(yintercept=1,linetype="dashed",col="grey")  
g + geom_point() + theme_bw()  + theme(axis.text.x=element_text(angle=-90), axis.title.x=element_blank()) + labs(title="Regulated phosphosite in the HUVEC", y="Fold Change (Log2)")  + scale_colour_manual(values = c("membrane"="springgreen4", "cytoplasm"="darkorange2"))

