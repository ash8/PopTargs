
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

## This is a function to plot the DAF / TAF 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
testTAB <- read.table("Paper_miRPopGen/Population_work/PopParts_2/HiLoExpMiR_TAF_DAF_byPop_real-ALL.txt", skip = 1)
tafs <- c(testTAB[, 10], testTAB[, 12], testTAB[, 14], testTAB[, 16], testTAB[, 18])
dafs <- c(testTAB[, 11], testTAB[, 13], testTAB[, 15], testTAB[, 17], testTAB[, 19])
PlotAF(tafs, dafs, "segf", "warg", "TAF")

## Function to plot the DAF / TAF 
PlotAF <- function(targFreq_1, targFreq_2, group1, group2, AF, title = T, KS = F, WX = F, Interact = T, freqMin = 0.01, freqMax = 0.99 ){

# Plot parameters:
# targFreq_1 -- a list of allele frequencies
# targFreq_2 -- a list of allele frequencies (to comapre with the above)
# group1 -- name the group in targFreq_1
# group2 -- name the group in targFreq_2
# AF -- make this either 'TAF' for target allele frequency or 'DAF' for derived allele freq
# Title = T -- default is to plot: paste(group1, " vs. ", group2), otherwise specify the string
# KS = F  -- make it 'T' to plot the KS test results
# WX = F  -- make it 'T' to plot the Wilcox test results
# freqMin =0.01  -- allele frequency cutoff
# freqMax = 0.99  -- allele frequency cutoff

    if (AF == 'TAF'){
            titleAF <- "Target Allele Frequency"
        }else if(AF == 'DAF'){
            titleAF <- "Derived Allele Frequency"
        }else{
            return("Not acceptable AF parameter, choose 'TAF' or 'DAF' ")
    }

    # take the selected cutoff in frequency (default is 0.01 - 0.99)
    targFreq_1_cut <- targFreq_1[c((targFreq_1 <= freqMax) & (targFreq_1 >= freqMin))]
    targFreq_2_cut <- targFreq_2[c((targFreq_2 <= freqMax) & (targFreq_2 >= freqMin))] 

    # Get the densities of each using hist() and then averaging the counts
    hist_TargFreq_1 <- hist(targFreq_1_cut, breaks = seq(0, 1, by = 0.05), plot = F)
    hist_p1 <- hist_TargFreq_1$counts / sum(hist_TargFreq_1$counts)

    hist_TargFreq_2 <- hist(targFreq_2_cut, breaks = seq(0, 1, by = 0.05), plot = F)
    hist_p2 <- hist_TargFreq_2$counts / sum(hist_TargFreq_2$counts)

    maxFreq = max(hist_p1, hist_p2) + 0.05
    
    if( title == T){
        mainName <- paste(group1, " vs. ", group2)
    }else{
        mainName <- title
    }
    
    OvInter <- length(targFreq_1_cut)
    UnInter <- length(targFreq_2_cut)

    # Do the plot 
    barplot(rbind(hist_p1, hist_p2), 
        beside = T,
        ylim = c(0, maxFreq), 
        xlab = titleAF, 
        ylab = "Density",
        main = mainName)
            
    # Draw the x axis:
    breaks <- seq(0, 1, by = 0.1)
    at <- axTicks(side=1, axp=c(par("xaxp")[1:2], length(breaks)-1))
    labels <- seq(min(breaks), max(breaks), length.out=1+par("xaxp")[3])
    labels <- round(labels, digits=1)
    axis(side=1, at=at, labels = breaks)
    
    # Add the legend
    legend(10, maxFreq -.01,
            c(group1, group2),
            pch=15,
            col=c("gray30", "gray95"),
            bty = "n")
     legend(10, maxFreq -.01,
            c(group1, group2),
            pch=22,
            col=c(1, 1),
            bty = "n")         

    # Add the number of interactions
    if (Interact == T ){
    mtext(paste(group1, "Interactions:", OvInter, "    ", group2, " Interactions:", UnInter), 
        line = .5, 
        cex = .6)
    }
    
    # Add the KS values if requested        
    if(KS == T){
    a <- targFreq_1_cut
    b <- targFreq_2_cut
    ksResult_pos <- ks.test(a, b, alternative = "greater")
    ksResult_neg <- ks.test(a, b, alternative = "less")
    ks_D_pos <- ksResult_pos[1]
    ks_p_pos <- ksResult_pos[2]
    ks_D_neg <- ksResult_neg[1]
    ks_p_neg <- ksResult_neg[2]
    
    mtext(paste("KS Test:", " D+ = ", round(-log10(as.numeric(ks_D_pos)), digits = 2), 
        " p-val = ", round(as.numeric(ks_p_pos), digits = 3)), 
        line = -.25, 
        cex = .6)
    mtext(paste("               ", "D- = ", round(-log10(as.numeric(ks_D_neg)), digits = 2), 
        " p-val = ", round(as.numeric(ks_p_neg), digits = 3)), 
        line = -.8, 
        cex = .6)
    }
    
    # Add the wilcox test values if requested        
    if(WX == T){
    a <- targFreq_1_cut
    b <- targFreq_2_cut
    wxResult_pos <- wilcox.test(a, b, alternative = "greater")
    wxResult_neg <- wilcox.test(a, b, alternative = "less")
    wx_D_pos <- wxResult_pos[1]
    wx_p_pos <- wxResult_pos$p.value
    wx_D_neg <- wxResult_neg[1]
    wx_p_neg <- wxResult_neg$p.value
    
    mtext(paste("Wilcox Test:", " D+ = ", round(-log10(as.numeric(wx_D_pos)), digits = 2), 
        " p-val = ", round(as.numeric(wx_p_pos), digits = 3)), 
        line = -.25, 
        cex = .6)
    mtext(paste("               ", "D- = ", round(-log10(as.numeric(wx_D_neg)), digits = 2), 
        " p-val = ", round(as.numeric(wx_p_neg), digits = 3)), 
        line = -.8, 
        cex = .6)
    }        
    
}





