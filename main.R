###packages/library#####
packages <- c("pheatmap", "naniar", "FactoMineR", "factoextra", 
              "corpcor", "fields", "ggplot2", "lattice", 
              "vegan", "rgl", "cluster", "ggfortify",
              "coRanking", "MASS", "clValid", "e1071", "outliers")
lapply(packages, library, character.only = T)
###讀資料###########################
musk <- read.csv("data/musk.csv",header = T,stringsAsFactors = F)
muskdata <- musk[2:nrow(musk),]
rownames(muskdata) <- muskdata[, 1]
muskdata <- muskdata[,2:170]
str(muskdata)

###資料處理###########################
dim(muskdata)
#continuous
#遺失值
na.col <- apply(muskdata, 2, function(x) sum(is.na(x)) > 0)
which(na.col == T) #沒有遺失值
miss_var_summary(muskdata)
#outliers
md.scal <- decostand(muskdata[,3:168], "standardize")
box.musk.scal <- apply(md.scal, 1, scale)
box.musk <- rev(order(apply(muskdata[,3:168], 2, median)))
boxplot(md.scal, main="Musk") #outliers-標準化
boxplot(muskdata[,3:168][box.musk], main="Musk") #outliers-原始
outlier(muskdata[,3:168])
#標準化
boxplot(muskdata[,c(18,165)], main="Musk",col=c("orange","grey")) 
boxplot(md.scal[,c(16,163)], main="Musk",col=c("orange","grey")) 


###EDA
a <- unique(muskdata[,1])#確實為102個分子
b <- table(muskdata[, 1])
range(muskdata[,3:168])
summary(muskdata[,3:168])

#categorical
barplot(table(muskdata[, 1]), main = "conformation of molecule", col="white",
        horiz = T, ylab = "molecule name", xlab = "conformation",cex.names = 0.2,las=1)

b <- as.data.frame(b)
colnames(b) <- c("molecule", "conformation")
#前10小的分子長條圖
head(b[order(b[,2]),], 10)
ggplot(head(b[order(b[,2]),], 10),aes(molecule,conformation))+geom_bar(stat="identity",fill="steelblue")
#前10大的分子長條圖
tail(b[order(b[,2]),], 10)
ggplot(tail(b[order(b[,2]),], 10),aes(reorder(molecule, -conformation),conformation))+geom_bar(stat="identity",fill="steelblue")
#https://zhuanlan.zhihu.com/p/27298708

#class
table(muskdata[, 169])
barplot(table(muskdata[, 169]), main = "Class of musk", col = c("orange","blue"))
legend("topright", legend = c("non-musk", "musk"), fill =  c("orange","blue"), cex = 1.2)


##看變數間的correlation
data_subset <- muskdata[, 3:168]
is.positive.definite(cov(data_subset))
corr <- cor(data_subset)
#畫corr的圖
cor.col <- two.colors(start="blue", middle="white",end="red") 
range.col <- floor((1+range(corr))*127+1)
corr.image <- image(t(corr)[, ncol(corr):1], main = "Correlation plot",
                                        col = cor.col[range.col[1]: range.col[2]], xaxt = "n", yaxt = "n")
#顏色條
rb.color.spectrum <- function(color, label){
       plot(0, 0, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", type="n")
       rasterImage(t(color), -0.8, 0, 0.8, 0.3, interpolate=FALSE)
       text(x=c(-0.8, 0, 0.8), y=c(-0.2, -0.2, -0.2), label= label,
                     cex = 1.8)}
col.table <- rb.color.spectrum(cor.col, c("-1", "0", "1"))

###3 Data Set###########################
###train data
splitdf <- function(df, train.ratio, seed=NULL) {
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(df)
  id <- sample(index, trunc(length(index)*train.ratio))
  train <- df[id, ]
  test <- df[-id, ]
  list(trainset=train,testset=test)
}

splits <- splitdf(muskdata, 0.2, 12345)
lapply(splits, dim)
DataSetA <- splits$trainset

splits_n_equal_p <- splitdf(DataSetA, 0.126, 12345)
DataSetB <- splits_n_equal_p$trainset

splits_n_small_p <- splitdf(DataSetA, 0.042, 12345)
DataSetC <- splits_n_small_p$trainset

table(DataSetA[, 169])
barplot(table(DataSetA[, 169]), main = "Class of musk A", col = c("orange","blue"))
legend("topright", legend = c("non-musk", "musk"), fill = c("orange","blue"), cex = 1.2)
table(DataSetB[, 169])
barplot(table(DataSetB[, 169]), main = "Class of musk B", col = c("orange","blue"))
legend("topright", legend = c("non-musk", "musk"), fill = c("orange","blue"), cex = 1.2)
table(DataSetC[, 169])
barplot(table(DataSetC[, 169]), main = "Class of musk C", col = c("orange","blue"))
legend("topright", legend = c("non-musk", "musk"), fill = c("orange","blue"), cex = 1.2)
###Dimension Reduction###########################
###DataSetA###########################
###PCA
dataA_subset <- as.matrix(DataSetA[,3:168])
mypca <- princomp(dataA_subset, cor=TRUE, scores=TRUE)
# 2D plot for first two components
pca.dim1a <- mypca$scores[,1]
pca.dim2a<- mypca$scores[,2]
pca.dim3a <- mypca$scores[,3]
mask_class <- unique(DataSetA[,169])
plot(pca.dim1a, pca.dim2a, main="PCA for musk A", xlab="1st PCA
     Component", ylab="2nd PCA Component",col=c(mask_class)+2, pch=c(mask_class)+2)
legend("topright", legend = c("non-musk", "musk"), pch=c(mask_class)+2, col = c(mask_class)+2, cex = 1)
summary(mypca)
eig.val <- get_eigenvalue(mypca)
fviz_eig(mypca, addlabels = TRUE, ylim = c(0, 50))

###LCMC
pca.lcmc.a <- coranking(as.matrix(dist(dataA_subset)), mypca$scores[,1:2], input = "data")
imageplot(pca.lcmc.a, main="co-ranking matrix of Musk A")
lcmc.a.pca <- LCMC(pca.lcmc.a, K = 5:10)
lcmc.a.pca

###ISOMAP
dataA_subset <- DataSetA[,3:168]
a.dist <- dist(dataA_subset)
iso.a.col <- DataSetA[,169]
open3d()
plot3d(pca.dim1a, pca.dim2a, pca.dim3a, col=iso.a.col+1, size=3,
       xlab="", ylab="", zlab="", axes = T)
sra.isomap <-  isomap(a.dist, ndim=2, k=10)
plot(sra.isomap, col=iso.a.col+1, main="ISOMAP of Musk Data A")

###LCMC
iso.lcmc.a <- coranking(as.matrix(a.dist), as.matrix(sra.isomap$points), input = "data")
imageplot(iso.lcmc.a, main="co-ranking matrix of Musk A")
lcmc.a <- LCMC(iso.lcmc.a, K = 5:10)
lcmc.a
###DataSetB###########################
###PCA
dataB_subset <- as.matrix(DataSetB[,3:168])
mypca <- princomp(dataB_subset, cor=TRUE, scores=TRUE)
# 2D plot for first two components
pca.dim1b <- mypca$scores[,1]
pca.dim2b <- mypca$scores[,2]
pca.dim3b <- mypca$scores[,3]
mask_class <- unique(DataSetB[,169])
plot(pca.dim1b, pca.dim2b, main="PCA for musk B", xlab="1st PCA
     Component", ylab="2nd PCA Component",col=c(mask_class)+2, pch=c(mask_class)+2)
legend("topright", legend = c("non-musk", "musk"), pch=c(mask_class)+2, col = c(mask_class)+2, cex = 1)
biplot(mypca)
summary(mypca)
eig.val <- get_eigenvalue(mypca)
fviz_eig(mypca, addlabels = TRUE, ylim = c(0, 50))

#lcmc
pca.lcmc.b <- coranking(as.matrix(dist(dataB_subset)), mypca$scores[,1:2], input = "data")
imageplot(pca.lcmc.b, main="co-ranking matrix of Musk B")
lcmc.b <- LCMC(pca.lcmc.b, K = 5:10)
lcmc.b

###ISOMAP
b.dist <- dist(dataB_subset)
iso.b.col <- DataSetB[,169]
open3d()
plot3d(pca.dim1b, pca.dim2b, pca.dim3b, col=iso.b.col+1, size=3,
       xlab="", ylab="", zlab="", axes = T)
srb.isomap <-  isomap(b.dist, ndim=2, k=35)
plot(srb.isomap, col=iso.b.col+1, main="ISOMAP of Musk Data B")

###LCMC
iso.lcmc.b <- coranking(as.matrix(b.dist), as.matrix(srb.isomap$points), input = "data")
imageplot(iso.lcmc.b, main="co-ranking matrix of Musk B")
lcmc.b <- LCMC(iso.lcmc.b, K = 5:10)
lcmc.b
###DataSetC###########################
###PCA
dataC_subset <- as.matrix(DataSetC[,3:168])
mypca <- prcomp(dataC_subset)
# 2D plot for first two components
pca.dim1c <- mypca$x[,1]
pca.dim2c <- mypca$x[,2]
pca.dim3c <- mypca$x[,3]
mask_class <- unique(DataSetC[,169])
plot(pca.dim1c, pca.dim2c, main="PCA for musk C", xlab="1st PCA
     Component", ylab="2nd PCA Component",col=c(mask_class)+2, pch=c(mask_class)+2)
legend("bottomright", legend = c("non-musk", "musk"), pch=c(mask_class)+2, col = c(mask_class)+2, cex = 1)
biplot(mypca)
summary(mypca)
eig.val <- get_eigenvalue(mypca)
fviz_eig(mypca, addlabels = TRUE, ylim = c(0, 50))

# n<<p prcomp
# https://c3h3notes.wordpress.com/2010/11/16/r%E7%9A%84%E5%85%A7%E5%BB%BApca%E6%8C%87%E4%BB%A4-princomp/
# https://cran.r-project.org/web/packages/ggfortify/vignettes/plot_pca.html

#lcmc
pca.lcmc.c <- coranking(as.matrix(dist(dataC_subset)), t(mypca$rotation)[,1:2], input = "data")
imageplot(pca.lcmc.c, main="co-ranking matrix of Musk C")
lcmc.c <- LCMC(pca.lcmc.c, K = 5:10)
lcmc.c
###ISOMAP
c.dist <- dist(dataC_subset)
iso.c.col <- DataSetC[,169]
open3d()
plot3d(pca.dim1c, pca.dim2c, pca.dim3c, col=iso.c.col+1, size=3,
       xlab="", ylab="", zlab="", axes = T)
src.isomap <-  isomap(c.dist, ndim=2, k=15)
plot(src.isomap, col=iso.c.col+1, main="ISOMAP of Musk Data C")

###LCMC
iso.lcmc.c <- coranking(as.matrix(c.dist), as.matrix(src.isomap$points), input = "data")
imageplot(iso.lcmc.c, main="co-ranking matrix of Musk C")
lcmc.c <- LCMC(iso.lcmc.c, K = 5:10)
lcmc.c
###cluster###########################
###DataSetA###########################
##kmeans
#k=2
dataA.kmeans2 <- kmeans(scale(dataA_subset), 2, 10)
plot(pca.dim1a, pca.dim2a, main="PCA for Musk Data with K-means Result",
     xlab="PCA-1", ylab="PCA-2", col=dataA.kmeans2$cluster)
#k=3
dataA.kmeans3 <- kmeans(scale(dataA_subset), 3, 10)
plot(pca.dim1a, pca.dim2a, main="PCA for Musk Data A with K-means Result",
     xlab="PCA-1", ylab="PCA-2", col=dataA.kmeans3$cluster)


##PCA2維
pca_dimemsion2 <- cbind(pca.dim1a, pca.dim2a)
#pca分3群熱圖
dataA.pca2kmeans2 <- kmeans(pca_dimemsion2, 3, 10)
pca2kmeans2 <- as.data.frame(dataA.pca2kmeans2$cluster)
names(pca2kmeans2) <- c("cluster")
pca2kmeans2$class <- DataSetA[,169]
pheatmap(pca_dimemsion2, show_rownames = FALSE,annotation_row = pca2kmeans2, cluster_cols = FALSE, cluster_rows = FALSE,main="Data A with K-means Result")

#階層式分群
DataSetA.ar <- as.data.frame(DataSetA[,169])
rownames(DataSetA.ar) <- rownames(DataSetA)
names(DataSetA.ar) <- c("class")
pheatmap(pca_dimemsion2, show_rownames = FALSE, annotation_row = DataSetA.ar, cluster_cols = FALSE, clustering_method = "average", main = "Hierarchical of Dataset A")
pheatmap(pca_dimemsion2, show_rownames = FALSE, annotation_row = DataSetA.ar, cluster_cols = FALSE, clustering_method = "complete",main = "Hierarchical of Dataset A")

###DataSetB###########################
##kmeans
#k=2
dataB.kmeans2 <- kmeans(scale(dataB_subset), 2, 10)
plot(pca.dim1b, pca.dim2b, main="PCA for Musk Data B with K-means Result",
     xlab="PCA-1", ylab="PCA-2", col=dataB.kmeans2$cluster)
#k=3
dataB.kmeans3 <- kmeans(scale(dataB_subset), 3, 10)
plot(pca.dim1b, pca.dim2b, main="PCA for Musk Data B with K-means Result",
     xlab="PCA-1", ylab="PCA-2", col=dataB.kmeans3$cluster)


#PCA2維分3群熱圖
pca_dimemsion2b <- cbind(pca.dim1b, pca.dim2b)
dataB.pca2kmeans3 <- kmeans(pca_dimemsion2b, 3, 10)
pca2kmeans3 <- as.data.frame(dataB.pca2kmeans3$cluster)
names(pca2kmeans3) <- c("cluster")
pca2kmeans3$class <- DataSetB[,169]
pheatmap(pca_dimemsion2b, show_rownames = FALSE, annotation_row = pca2kmeans3, cluster_cols = FALSE, cluster_rows = FALSE,main="Data B with K-means")

#階層式分群
DataSetB.ar <- as.data.frame(DataSetB[,169])
rownames(DataSetB.ar) <- rownames(DataSetB)
names(DataSetB.ar) <- c("class")
pheatmap(pca_dimemsion2b, show_rownames = FALSE, annotation_row = DataSetB.ar, cluster_cols = FALSE, clustering_method = "average", main = "Hierarchical of Dataset B")
pheatmap(pca_dimemsion2b, show_rownames = FALSE, annotation_row = DataSetB.ar, cluster_cols = FALSE, clustering_method = "complete",main = "Hierarchical of Dataset B")

###DataSetC###########################
##kmeans
#k=2
dataC.kmeans2 <- kmeans(scale(dataC_subset), 2, 10)
plot(pca.dim1, pca.dim2, main="PCA for Musk Data C with K-means Result",
     xlab="PCA-1", ylab="PCA-2", col=dataC.kmeans2$cluster)
#k=3
dataC.kmeans3 <- kmeans(scale(dataC_subset), 3, 10)
plot(pca.dim1, pca.dim2, main="PCA for Musk Data C with K-means Result",
     xlab="PCA-1", ylab="PCA-2", col=dataC.kmeans3$cluster)

##PCA2維
pca_dimemsion2c <- cbind(pca.dim1c, pca.dim2c)
#pca分3群熱圖
dataC.pca2kmeans3 <- kmeans(pca_dimemsion2c, 3, 10)
pca2kmeans3 <- as.data.frame(dataC.pca2kmeans3$cluster)
names(pca2kmeans3) <- c("cluster")
pca2kmeans3$class <- DataSetC[,169]
pheatmap(pca_dimemsion2c, show_rownames = FALSE, annotation_row = pca2kmeans3, cluster_cols = FALSE, cluster_rows = FALSE,main="Data C with K-means")

#階層式分群
DataSetC.ar <- as.data.frame(DataSetC[,169])
rownames(DataSetC.ar) <- rownames(DataSetC)
names(DataSetC.ar) <- c("class")
pheatmap(pca_dimemsion2c, show_rownames = FALSE, annotation_row = DataSetC.ar, cluster_cols = FALSE, clustering_method = "average", main = "Hierarchical of Dataset C")
pheatmap(pca_dimemsion2c, show_rownames = FALSE, annotation_row = DataSetC.ar, cluster_cols = FALSE, clustering_method = "complete",main = "Hierarchical of Dataset C")

###cluster validation####
###DatasetA####
internA <- clValid(dataA_subset, 2:6, clMethods=c("hierarchical", "kmeans", "pam"),
                   validation="internal")
y
summary(internA)
par(mfrow=c(1, 3))
plot(internA)

stabA <- clValid(dataA_subset, 2:6, clMethods=c("hierarchical", "kmeans"),
                   validation="stability")
y
summary(stabA)
par(mfrow=c(1, 4))
plot(stabA)
###DatasetB####
internB <- clValid(dataB_subset, 2:6, clMethods=c("hierarchical", "kmeans"),
                  validation="internal")
y
summary(internB)
par(mfrow=c(1, 3))
plot(internB)

stabB <- clValid(dataB_subset, 2:6, clMethods=c("hierarchical", "kmeans"), validation="stability")
y
summary(stabB)
par(mfrow=c(1, 4))
plot(stabB)
###DatasetC####
internC <- clValid(dataC_subset, 2:6, clMethods=c("hierarchical", "kmeans"),
                   validation="internal")
y
summary(internC)
par(mfrow=c(1, 3))
plot(internC)

stabC <- clValid(dataC_subset, 2:6, clMethods=c("hierarchical", "kmeans"),
                 validation="stability")
y
summary(stabC)
par(mfrow=c(1, 4))
plot(stabC)

####classification####
#datasetA
modelA <- svm(dataA_subset, DataSetA[,169], type="C-classification")
print(modelA)
summary(modelA)
# test with train data and report accuracy
pred <- predict(model, dataA_subset)
table(pred, DataSetA[,169])
sum(diag(table(pred, DataSetA[,169])))/length(DataSetA[,169])

#datasetB
modelB <- svm(dataB_subset, DataSetB[,169], type="C-classification")
print(modelB)
summary(modelB)
# test with train data and report accuracy
pred <- predict(modelB, dataB_subset)
table(pred, DataSetB[,169])
sum(diag(table(pred, DataSetB[,169])))/length(DataSetB[,169])

#datasetC
modelC <- svm(dataC_subset, DataSetC[,169], type="C-classification")
print(modelC)
summary(modelC)
# test with train data and report accuracy
pred <- predict(modelC, dataC_subset)
table(pred, DataSetC[,169])
sum(diag(table(pred, DataSetC[,169])))/length(DataSetC[,169])

###Testing data for classification####
splits <- splitdf(muskdata, 0.7, 123)
lapply(splits, dim)
musk_Testset <- splits$testset
musk_Testsubset <- musk_Testset[,3:168]

table(musk_Testset[, 169])
barplot(table(musk_Testset[, 169]), main = "Class of musk_Testset", col = c("orange","blue"))
#modelA
pred <- predict(modelA, musk_Testsubset)
table(pred, musk_Testset[,169])
sum(diag(table(pred, musk_Testset[,169])))/length(musk_Testset[,169])
#modelB
pred <- predict(modelB, musk_Testsubset)
table(pred, musk_Testset[,169])
sum(diag(table(pred, musk_Testset[,169])))/length(musk_Testset[,169])
#modelC
pred <- predict(modelC, musk_Testsubset)
table(pred, musk_Testset[,169])
sum(diag(table(pred, musk_Testset[,169])))/length(musk_Testset[,169])




