interf = data.frame(read.csv("aa_school/biophysics_project/MutationalAnalysis_IntefaceOfDistance_8.0.tsv", header = T))
seq = data.frame(read.csv("aa_school/biophysics_project/SequentialIntroductionThroughDistance.tsv", header = T))

library(ggplot2)

paste( seq$Residue[97:112], seq$ResID[97:112])

View(seq)
View(interf)

sum(seq$Contributing_Interaction_Energy)

neg = as.factor(seq$Contributing_Interaction_Energy < 0)

par(mar = c(8, 8, 4, 8) )

max(seq$Contributing_Interaction_Energy)

plot(seq$Contributing_Interaction_Energy~seq$Threshold_Distance, 
     type = "h",
     ylim = c(-10, 5),
     xlim = c(0, 70),
     ylab = "Contribution to Interaction Energy\n(absolute value)", 
     main = "Energetic Contribution to Interaction Energy\nas a Function of Distance", 
     xlab = "Distance in Angstrom to Nearest\nNeighbour of Opposite Chain", 
     frame.plot=F
     )


attach(interf)
sum(Non.ALA_Electrostatics + Non.ALA_Vdw + Non.ALA_AB.Solvation - Non.ALA_Chain_Solvation)

sum(seq$Contributing_Interaction_Energy)
length(interf$ResID)

pairs(data.frame(seq$ΔG_of_ALA_mutation, 
        seq$Threshold_Distance, 
        seq$Non.ALA_Chain_Solvation, 
        seq$Non.ALA_Electrostatics, 
        seq$Non.ALA_Vdw, 
        seq$Non.ALA_AB_Solvation),
        col = as.factor(seq$Chain)
        )

plot(seq$ΔG_of_ALA_mutation ~ seq$,
     ylab = "Contribution to Interaction Energy\n(absolute value)", 
     main = "Energetic Contribution to Interaction Energy\nas a Function of Distance", 
     xlab = "Distance in Angstrom to Nearest\nNeighbour of Opposite Chain", 
     frame.plot=T)

plot(sort(seq$ΔG_of_ALA_mutation),
          type = "h",
          ylab = "Contribution to Interaction Energy\n(absolute value)", 
          main = "Energetic Contribution to Interaction Energy\nas a Function of Distance", 
          xlab = "Distance in Angstrom to Nearest\nNeighbour of Opposite Chain", 
          frame.plot=F
)


library(tidyr)
library(ggplot2)
library(GGally)

names(interf)

ploting = data.frame(
                     "Distance" = seq$Threshold_Distance[interf$ResID],
                     "Chain Solvation" = interf$Non.ALA_Chain_Solvation, 
                     "Electrostatic Attraction" = interf$Non.ALA_Electrostatics, 
                     "Van der Waals" = interf$Non.ALA_Vdw, 
                     "Complex Solvation" = interf$Non.ALA_AB.Solvation,
                     "ABS(ΔG of Ala-mutation)" = abs(interf$ΔG_of_ALA_mutation)
                     )

ggpairs(ploting, 
        aes(color = interf$Chain, alpha = 0.5),
        lower = list(continuous = "points"),
        )
ggpairs(ploting, 
        lower = list(continuous = "points"),
        )

par(mar = c(8, 8, 4, 8) )
barplot(interf$ΔG_of_ALA_mutation, 
        main = "Energetic Effect of Mutation",
        ylab = "ΔG of mutation",
        ylim = c(-15, 15),
        xlab = "Residues",
        #names.arg = paste(interf$Residue, interf$ResID, interf$Chain),
        las = 2, 
        cex.names = 0.7,
        cex.main = 2,
        cex.lab = 1.5
        )

plot(abs(interf$ΔG_of_ALA_mutation), 
        type = "h",
        main = "Energetic Effect of Mutation",
        ylab = "Absolute ΔG of mutation",
        ylim = c(0, 16),
        xlab = "Interface Residues",
        xlim = c(0, 120),
        names.arg = paste(interf$Residue, interf$ResID, interf$Chain),
        las = 2, 
        cex.names = 0.6,
        space = 0,
        frame.plot = F
        )
install.packages("plotrix")
library("plotrix")

extrema = c(seq(0, 15, by = 1), seq(length(interf$ΔG_of_ALA_mutation) - 15, length(interf$ΔG_of_ALA_mutation), by = 1))


barplot(interf$ΔG_of_ALA_mutation[extrema], 
        main = "Energetic Effect of Mutation",
        ylab = "ΔG of mutation",
        ylim = c(-15, 15),
        #xlab = "Residues",
        names.arg = paste(interf$Residue[extrema], interf$ResID[extrema], interf$Chain[extrema]),
        las = 2, 
        cex.names = 1,
        cex.main = 2,
        cex.lab = 1.5
)
