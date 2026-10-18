lib <- file.path(getwd(), "rlibs")

if (!requireNamespace("glmnet", lib.loc = lib, quietly = TRUE) ||
    packageVersion("glmnet", lib.loc = lib) < "5.1") {
  install.packages("glmnet", lib = lib, repos="https://cran.r-project.org/")
}
library(glmnet, lib.loc = lib)

set.seed(123)

#Rscript to generate single gene graph for IKBKG data
#Usage: Rscript plot_ikbkg_cn.r infile

MER67B_COORDINATES = "refs/beds/ikbkg_mer67b_coordinates.bed"
HOM_REGION_BED = "refs/beds/ikbkg_ikbkgp1_COMPLETE_HOMOLOGOUS_REGIONS.bed"

args = commandArgs(trailingOnly=TRUE)
infile = args[1]
id = gsub("\\.csv$|Fx-", "", basename(infile))
gene = "IKBKG"

# only interested in a small portion of the genome 
ROI_START = 154541239
ROI_END = 154605971

hom_region_df <- read.table(HOM_REGION_BED, sep = "\t", header = FALSE, nrows = 1)

# remember, start coords of bed files are 0-based
# getting only first hom region
HOM_START = hom_region_df[1,2] + 1
HOM_END = hom_region_df[1,3]

mer67b_df <- read.table(MER67B_COORDINATES, sep = "\t", header = FALSE, nrows = 2)
 
MER67B_1_Start = mer67b_df[1,2] + 1
MER67B_1_End = mer67b_df[1,3]

MER67B_2_Start = mer67b_df[2,2] + 1
MER67B_2_End = mer67b_df[2,3]
 
dummy1 <- read.csv(infile)
 
dummy2 <- dummy1[dummy1$Chr == "chrX" &
                 dummy1$Position >= ROI_START &
                 dummy1$Position <= ROI_END, ]
 
png(file = paste("Rplot","_",id,"_",gene,".png", sep=""), width = 7, height = 7, units = "in", res = 300)
par(mar=c(6,5,5,5)+0.1, oma = c(1,1,1,1))
plot(dummy2$Position,dummy2$copynumber, type = "n", main=paste(id,gene),ylab="Copy Number", xlab=paste("chrX (hg38)"), ylim=c(0,6), cex = 0.3)
abline(h=1.0,col="red",lwd=1)
abline(h=2.0,col="green",lwd=1)
abline(h=3.0,col="red",lwd=1)
 
# added indicator for where IKBKG hom region should be
abline(v=HOM_START,col="purple",lwd=1)
abline(v=HOM_END,col="purple",lwd=1)
 
# Mer67bs
abline(v=MER67B_1_Start,col="blue",lwd=1)
abline(v=MER67B_1_End,col="blue",lwd=1)
abline(v=MER67B_2_Start,col="blue",lwd=1)
abline(v=MER67B_2_End,col="blue",lwd=1)
 
# Compute mean CN for the following regions: 
# ROI_Start - HOM_START
# HOM_START - Mer67B 1
# MER67B_1 - MER67B_2
# MER67B_2 - Hom_1 End
# Hom_1 End - ROI_END

region = 1
for (vec in list(c(ROI_START, HOM_START - 1),
                 c(HOM_START, MER67B_1_Start - 1),
                 c(MER67B_1_End + 1, MER67B_2_Start - 1), # DELETION REGION
                 c(MER67B_2_End + 1, HOM_END),
                 c(HOM_END + 1, ROI_END))) {
    df_subset <- dummy2[dummy2$Position >= vec[1] & dummy2$Position <= vec[2], ]
 
    mean = mean(df_subset$copynumber, na.rm = TRUE)
    ci  <- t.test(df_subset$copynumber)$conf.int

    # ignore the c(HOM_START, MER67B_1_Start - 1) region
    if (region == 2) {
        region = region + 1
        next
    }

    region_text = paste("Region", as.character(region), sep=" ")
    if (region == 3) {
        region_text = "Deletion Region"
    }

    text(x = (vec[1] + vec[2]) / 2, y = 3.85,
            labels = region_text, cex = 0.6, col = "black", font = 2)

    segments(x0 = vec[1], y0 = mean,
             x1 = vec[2], y1 = mean,
             col = "black", lwd = 2)

    text(x = (vec[1] + vec[2]) / 2, y = 3.5,
            labels = paste("Mean CN:", as.character(round(mean, 2))), cex = 0.6, col = "black", font = 2)

    segments(x0 = vec[1], y0 = ci[2],
             x1 = vec[2], y1 = ci[2],
             col = "black", lwd = 1,  lty = "dashed")

    segments(x0 = vec[1], y0 = ci[1],
             x1 = vec[2], y1 = ci[1],
             col = "black", lwd = 1,  lty = "dashed")

    region = region + 1
}
 
# Welch two-sample t-tests: deletion region (region 3) vs region 1 and region 4
del_cn  <- dummy2[dummy2$Position >= MER67B_1_End + 1 & dummy2$Position <= MER67B_2_Start - 1, "copynumber"]
reg1_cn <- dummy2[dummy2$Position >= ROI_START    & dummy2$Position <= HOM_START - 1,      "copynumber"]
reg4_cn <- dummy2[dummy2$Position >= MER67B_2_End + 1 & dummy2$Position <= HOM_END,        "copynumber"]
reg5_cn <- dummy2[dummy2$Position >= HOM_END + 1      & dummy2$Position <= ROI_END,            "copynumber"]

# Need to fix the below, giving us very low P-values despite whatever method gets applied to it
# This is because the t-test treats every individual base pair as independent, (adjacent base pairs are not independent)
# The test therefore thinks that it has > 10,000 degrees of freedom
# In other words, you are feeding it thousands of correlated measurements and it treats all of them as independent
# Hence, take a fixed window (1000bp) and get the mean of that, and pass that to the test

to_windows <- function(cn_vec, w = 1000) {
    starts <- seq(1, length(cn_vec), by = w)
    result <- numeric(length(starts))
    for (i in seq_along(starts)) {
        window <- cn_vec[starts[i] : min(starts[i] + w - 1, length(cn_vec))]
        result[i] <- mean(window, na.rm = TRUE)
    }
    return(result)
}

del_w  <- to_windows(del_cn)
reg1_w <- to_windows(reg1_cn)
reg4_w <- to_windows(reg4_cn)
reg5_w <- to_windows(reg5_cn)

t_vs_reg1 <- t.test(del_w, reg1_w, var.equal = FALSE)
t_vs_reg4 <- t.test(del_w, reg4_w, var.equal = FALSE)

p_adjusted <- p.adjust(c(t_vs_reg1$p.value, t_vs_reg4$p.value), method = "bonferroni")

cat("\n--- Welch two-sample t-test: Deletion region vs Region 1 (unique upstream) ---\n")
print(p_adjusted[1])
cat("\n--- Welch two-sample t-test: Deletion region vs Region 4 (MER67B_2 to HOM_END) ---\n")
print(p_adjusted[2])

par(new=TRUE)

# FINAL STEP, GLMNET PREDICTION:

get_ratios <- function(del_w, reg_cn, reg_name) {
    ratios <- del_w / mean(reg_cn, na.rm = TRUE)
    names(ratios) <- paste0("del_", 1:length(del_w), "_", reg_name)
    return(ratios)
}

sample_ratios <- c(get_ratios(del_w, reg1_cn, "r1"),
                   get_ratios(del_w, reg4_cn, "r4"),
                   get_ratios(del_w, reg5_cn, "r5"))

sample_df <- as.data.frame(t(sample_ratios))
rownames(sample_df) <- id

model <- readRDS("refs/models/cvfit_category.rds")
probability <- predict(model, newx = sample_df, s = "lambda.min", type="response")

print(probability)

mtext(bquote(bold("Probablity of sample having deletion:") ~ .(probability[1])),
      outer = TRUE,
      side = 1,
      line = -1.6,
      adj = 1,
      cex = 0.8)

# adj = horizontal position
# line = vertical position
# side = which edge of page

del_vs_r1_FC <- log2(mean(del_cn, na.rm = TRUE)/mean(reg1_cn, na.rm = TRUE))
del_vs_r4_FC <- log2(mean(del_cn, na.rm = TRUE)/mean(reg4_cn, na.rm = TRUE))

# in the bquote function, most of the expression is treated literally but anything inside .(xxx) is evaluated as normal R code.
# ~ means insert a space character in between the 2 pieces
mtext(bquote(bold("p-value(Del vs R1):") ~ .(p_adjusted[1]) * ", " *  bold("log2FC(Del vs R1):") ~ .(del_vs_r1_FC)), # .(fmt_p(t_vs_reg1$p.value))),
      outer = TRUE,
      side = 1,
      line = -0.7,
      adj = 1,
      cex = 0.8)

mtext(bquote(bold("p-value(Del vs R4):") ~ .(p_adjusted[2]) * ", " *  bold("log2FC(Del vs R4):") ~ .(del_vs_r4_FC)), # .(fmt_p(t_vs_reg4$p.value))),
      outer = TRUE,
      side = 1,
      line = 0.2,
      adj = 1,
      cex = 0.8)

# adding exons to plot

exons_df <- read.table("refs/beds/ikbkg_exons.bed", sep = "\t", header = FALSE)

# Exons (bed start is 0-based, so add 1)
exons_df <- read.table("refs/beds/ikbkg_exons.bed", sep = "\t", header = FALSE)

# added an offset so line doesnt stick out from exons 1-1 respectively
ikbkg_start = exons_df[1, 2] + 100
ikbkg_end = exons_df[10, 3] - 100

for (i in 1:nrow(exons_df)) {
    exon_start = exons_df[i, 2] + 1
    exon_end = exons_df[i, 3]
    rect(xleft = exon_start, ybottom = 0.35, xright = exon_end, ytop = 0.65,
         col = "darkorange", border = "darkorange")

    # even exons labelled above the box, odd exons below
    if (i %% 2 == 0) {
        label_y = 0.2
    } else {
        label_y = 0.8
    }
    if (i == 1) {
        text(x = (exon_start + exon_end) / 2, y = label_y,
            labels = paste("Exon", i, sep = " "), cex = 0.5, col = "black")
    } else {
        text(x = (exon_start + exon_end) / 2, y = label_y,
            labels = i, cex = 0.5, col = "black")
    }
}

segments(x0 = ikbkg_start, y0 = 0.5,
         x1 = ikbkg_end, y1 = 0.5,
         col = "darkorange", lwd = 2)

text(x = (ikbkg_start + ikbkg_end) / 2, y = 0.1,
            labels = "IKBKG", cex = 0.5, col = "black")

dev.off()

