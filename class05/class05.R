#' ---
#' title: "Class 5: Data Visualization and Graphs in R"
#' author: "Christopher Cheng"
#' date: "January 23rd, 2020"
#' ---

# Class 5
# Data Visualization and graphs in R

#Need to import/read input data file first
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)

# A basic plot of the weight vs age
plot(weight$Age, weight$Weight, xlab = "Age (months)", ylab = "Weight (kg)", main = "Baby weight with age", col = "blue", type = "o", pch = 15, ylim = c(2,10))

#Solution 
plot(weight$Age, weight$Weight, typ="o", 
     pch=15, cex=1.5, lwd=2, ylim=c(2,10), 
     xlab="Age (months)", ylab="Weight (kg)", 
     main="Baby weight with age", col="blue") 


# A silly example of 'pch' plot character and 'cex' size
plot(1:5, cex = 1:5, pch = 1:5)

# Reading in feature_counts.txt 
mouse <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE, sep = "\t")

par(mar = c(5,11,2,1))
barplot(mouse$Count, horiz = TRUE, names.arg = mouse$Feature,col = "lightblue", main = "Number of features in the mouse GRCm38 genome", las =1)


par(mar = c(5,4,2,2))
plot(1:10)

#Solution
par(mar=c(3.1, 11.1, 4.1, 2))
barplot(mouse$Count, names.arg=mouse$Feature, 
        horiz=TRUE, ylab="", 
        main="Number of features in the mouse GRCm38 genome", 
        las=1, xlim=c(0,80000))

#Section 2C: Histograms
x <- c(rnorm(10000), rnorm(10000)+4)
hist(x, breaks = 80)

# Section 3  
# Using Color in plots

sex <- read.delim("bimm143_05_rstats/male_female_counts.txt", header = TRUE, sep = '\t')

barplot(sex$Count, names.arg = sex$Sample, col = rainbow(nrow(sex)), las = 2, ylab = "Counts")

#Solution
mf <- read.delim("bimm143_05_rstats/male_female_counts.txt")

barplot(mf$Count, names.arg=mf$Sample, col=rainbow(nrow(mf)), 
        las=2, ylab="Counts")

# another plot of the same t hing with different colors
barplot(sex$Count, names.arg = sex$Sample, col = c("red", "blue", "darkgreen"), las = 2, ylab = "Counts", main= "Bert")

#Solution Part 2
barplot(mf$Count, names.arg=mf$Sample, col=c("blue2","red2"), 
        las=2, ylab="Counts")

#Section 3B: Coloring by Value
genes <- read.delim("bimm143_05_rstats/up_down_expression.txt")
nrow(genes)

table(genes$State)

plot(genes$Condition1, genes$Condition2, col=genes$State, 
     xlab="Expression condition 1", ylab="Expression condition 2")

palette()
levels(genes$State)

#Solution
palette(c("blue","gray","red"))
plot(genes$Condition1, genes$Condition2, col=genes$State, xlab="Expression condition 1", ylab="Expression condition 2")

#Section 3C: Dynamic use of color

# Lets plot expresion vs gene regulation
meth <- read.delim("bimm143_05_rstats/expression_methylation.txt")

nrow(meth)

plot(meth$gene.meth, meth$expression)

dcols <- densCols(meth$gene.meth, meth$expression)

# Plot changing the plot character ('pch') to a solid circle
plot(meth$gene.meth, meth$expression, col = dcols, pch = 20)

# Find the indices of genes with above 0 expresion
inds <- meth$expression > 0

# Plot just these genes
plot(meth$gene.meth[inds], meth$expression[inds])

## Make a desnisty color vector for these genes and plot
dcols <- densCols(meth$gene.meth[inds], meth$expression[inds])

plot(meth$gene.meth[inds], meth$expression[inds], col = dcols, pch = 20)

dcols.custom <- densCols(meth$gene.meth[inds], meth$expression[inds],
                         colramp = colorRampPalette(c("blue2",
                                                      "green2",
                                                      "red2",
                                                      "yellow")) )

plot(meth$gene.meth[inds], meth$expression[inds], 
     col = dcols.custom, pch = 20)





















