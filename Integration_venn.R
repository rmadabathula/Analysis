library(readxl)
library(tibble)
library(BiocManager)
BiocManager::install("preprocessCore")
library(preprocessCore)

#########################Common&exclusive metabolites#############################

sigMetabolitesSummary <- read.csv("VennMetabolites.csv", header = TRUE,  row.names = 1)


########################## Gene Expression Matrix ##################################

Gene_expression <- read_excel("Venn_genes.xlsx")

Gene_expression <- column_to_rownames(Gene_expression, var = names(Gene_expression)[1])

Gene_expression <- Gene_expression[ apply(Gene_expression, 1, min) >= 1, ]



######################## omics ##########################################

BiocManager::install("mixOmics")
library(mixOmics)

TGene_expression <- t(Gene_expression)

TsigMetaboliteSummary <- t(sigMetabolitesSummary)

## use the gene expression data as the X matrix ##

x <- TGene_expression

## use the metabolite data as the Y matrix ##

y <- TsigMetaboliteSummary

pls.result <- pls(x, y)   ### run the partial least square method ###

plotIndiv(pls.result)  ### plot the samples ###

plotVar(pls.result)   ### plot the variables ###

spls.result <- spls(x, y, keepX = c(10, 20), keepY = c(3, 2))   ### run the method ###

plotIndiv(spls.result) ### plot the samples ### 

plotVar(spls.result)    ### plot the variables ###

spls.grafting <- spls(X = x, Y = y, ncomp = 5, mode = 'regression', near.zero.var = TRUE)

perf.pls.grafting <- perf(spls.grafting,
                          validation = 'Mfold', folds = 7,
                          nrepeat = 50, progressBar = TRUE)

plot(perf.pls.grafting, criterion = 'Q2.total')

# set range of test values for number of variables to use from X dataframe
list.keepX <- c(seq(20, 50, 5))

# set range of test values for number of variables to use from Y dataframe
list.keepY <- c(3:10)

tune.spls.grafting <- tune.spls(x, y, ncomp = 2,
                                test.keepX = list.keepX,
                                test.keepY = list.keepY,
                                nrepeat = 1, folds = 10, # use 10 folds
                                mode = 'regression', measure = 'cor')

plot(tune.spls.grafting)         # use the correlation measure for tuning

tune.spls.grafting$choice.keepX
tune.spls.grafting$choice.keepY

optimal.keepX <- tune.spls.grafting$choice.keepX

# extract optimal number of variables for Y datafram
optimal.keepY <- tune.spls.grafting$choice.keepY

optimal.ncomp <-  length(optimal.keepX)

final.spls.grafting <- spls(x, y, ncomp = optimal.ncomp,
                            keepX = optimal.keepX,
                            keepY = optimal.keepY,
                            mode = "regression")
metadata <- read.csv(file = "ColDataWMWF.csv", header = TRUE, row.names = 1)


Xvariateplot <- plotIndiv(final.spls.grafting, ind.names = FALSE,
                          rep.space = "X-variate", # plot in X-variate subspace
                          group = metadata$Infestation_Type, # colour by time group
                          pch = as.factor(metadata$Infestation),
                          col.per.group = color.mixo(1:10),
                          legend = TRUE, legend.title = 'Type', legend.title.pch = 'Infestation')

Yvariateplot <- plotIndiv(final.spls.grafting, ind.names = FALSE,
                          rep.space = "Y-variate", # plot in X-variate subspace
                          group = metadata$Infestation_Type, # colour by time group
                          pch = as.factor(metadata$Infestation),
                          col.per.group = color.mixo(1:10),
                          legend = TRUE, legend.title = 'Type', legend.title.pch = 'Infestation')
XYvariateplot <- plotIndiv(final.spls.grafting, ind.names = FALSE,
                           rep.space = "XY-variate", # plot in X-variate subspace
                           group = metadata$Infestation_Type, # colour by time group
                           pch = as.factor(metadata$Infestation),
                           col.per.group = color.mixo(1:10),
                           legend = TRUE, legend.title = 'Type', legend.title.pch = 'Infestation')

col.tox <- color.mixo(as.numeric(as.factor(metadata$Infestation_Type))) # create set of colours

plotVar(final.spls.grafting, cex = c(3,4), var.names = c(FALSE, TRUE))

x11()
plotArrow(final.spls.grafting, ind.names = FALSE,
          group = metadata$Infestation_Type, # colour by time group
          col.per.group = color.mixo(1:10),
          legend.title = 'Type')
X11()
cim(final.spls.grafting, comp = 1:2, xlab = "metabolites", ylab = "genes")
final.spls.grafting$Y
color.edge <- color.GreenRed(50)
x11()
network(final.spls.grafting, comp = 1:2,
        cutoff = 0.5, # only show connections with a correlation above 0.7
        shape.node = c("rectangle", "circle"), size.node = 0.001,
        color.node = c("cyan", "pink"), interactive = TRUE,
        color.edge = color.edge, # save as a png to the current working directory
        name.save = 'sPLS Liver Toxicity Case Study Network Plot', save = "pdf")
x11()
imgCor(x, y, sideColors = c("purple", "green"))

#set grid search values for each regularisation parameter
cor_X <- cor(x)
cor_Y <- cor(y)

# Visualize correlation matrices
heatmap(cor_X, symm = TRUE)
heatmap(cor_Y, symm = TRUE)

det(cor_X)  # If determinant is near 0, X is nearly singular
det(cor_Y)
install.packages("caret")
library(caret)

highCorr_X <- findCorrelation(cor_X, cutoff = 0.9)  # Remove variables with correlation > 0.9
highCorr_Y <- findCorrelation(cor_Y, cutoff = 0.9)

X_filtered <- x[, -highCorr_X]
Y_filtered <- y[, -highCorr_Y]
X_scaled <- scale(X_filtered)
Y_scaled <- scale(Y_filtered)

pca_X <- prcomp(X_scaled, center = TRUE, scale. = TRUE)
pca_Y <- prcomp(Y_scaled, center = TRUE, scale. = TRUE)

# Choose number of PCs that explain 95% variance
explained_var_X <- cumsum(pca_X$sdev^2) / sum(pca_X$sdev^2)
explained_var_Y <- cumsum(pca_Y$sdev^2) / sum(pca_Y$sdev^2)

num_pcs_X <- which(explained_var_X >= 0.95)[1]
num_pcs_Y <- which(explained_var_Y >= 0.95)[1]

X_pca <- pca_X$x[, 1:num_pcs_X]
Y_pca <- pca_Y$x[, 1:num_pcs_Y]


  cv.tune.rcc.grafting <- tune.rcc(x, y,
                                 grid1 = seq(0.001, 1, length.out = 10),
                                 grid2 = seq(0.001, 1, length.out = 10),
                                 validation = "Mfold", folds = 5)

cv.tune.rcc.grafting


opt.l1 <- cv.tune.rcc.grafting$opt.lambda1 # extract the optimal lambda values
opt.l2 <- cv.tune.rcc.grafting$opt.lambda2

# formed optimised CV rCCA

CV.rcc.grafting <- rcc(X_scaled, Y_scaled, method = "ridge",
                       lambda1 = opt.l1, lambda2 = opt.l2)

shrink.rcc.grafting <- rcc(X_scaled,Y_scaled, method = 'shrinkage')

# examine the optimal lambda values after shrinkage

shrink.rcc.grafting$lambda


plot(CV.rcc.grafting, type = "barplot", main = "Cross Validation")

# barplot of shrinkage method rCCA canonical correlations
plot(shrink.rcc.grafting, type = "barplot", main = "Shrinkage")

plotIndiv(CV.rcc.grafting, comp = 1:2,
          ind.names = metadata$Type,
          group = metadata$Infestation_Type, rep.space = "XY-variate",
          legend = TRUE, title = 'Watermelon Grafting, rCCA CV XY-space')

plotIndiv(shrink.rcc.grafting, comp = 1:2,
          ind.names = metadata$Type,
          group = metadata$Infestation_Type, rep.space = "XY-variate",
          legend = TRUE, title = '(a) Grafting, rCCA CV XY-space')

plotArrow(CV.rcc.grafting, group = metadata$Infestation_Type,
          col.per.group = color.mixo(1:10),
          title = '(a) Nutrimouse, CV method')


plotArrow(shrink.rcc.grafting, group = metadata$Infestation_Type,
          col.per.group = color.mixo(1:10),
          title = '(a) Nutrimouse, CV method')

plotVar(CV.rcc.grafting, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = '(a) Nutrimouse, rCCA CV comp 1 - 2')

plotVar(shrink.rcc.grafting, var.names = c(TRUE, TRUE),
        cex = c(4, 4), cutoff = 0.5,
        title = '(b) Nutrimouse, rCCA shrinkage comp 1 - 2')
x11()
network(CV.rcc.grafting, comp = 1:2, interactive = TRUE, size.node = 0.1,
        lwd.edge = 2,
        cutoff = 0.5)
x11()
cim(shrink.rcc.grafting, comp = 1:2, xlab = "genes", ylab = "metabolome", save = "pdf" )


data = list(mRNA = x,
            Metabolomics = y)

ZZ = metadata$Infestation_Type
design = matrix(0.1, ncol = length(data), nrow = length(data),
                dimnames = list(names(data), names(data)))
list.keepX = list(mRNA = c(6,14), Metabolomics = c(5,18))

sgccda.res = block.splsda(X = data, Y = ZZ, ncomp = 2,
                          keepX = list.keepX, design = design)
x11()
circosPlot(sgccda.res, cutoff = 0.4, line = TRUE,
           color.blocks= c('darkorchid', 'orange1'),
           color.cor = c("lightcoral","blue"), size.labels = 1, size.legend = 1, save="pdf" )
sgccda.res$design



































































































