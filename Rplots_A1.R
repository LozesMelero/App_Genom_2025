setwd("~/Desktop/quantinemo_mac/Exercises/Exercice_A/A1")
main='~/Desktop/quantinemo_mac/Exercises/Exercice_A/A1/'
library(plotrix)
N=c(50,1000)
m=c('0.0','0.5')
par(mfrow=c(2,2))
for (i in c(1:length(N))){
        for (j in c(1:length(m))){
                setwd(paste0(main,paste0('N',N[i],'_mig',m[j])))    
                means=read.table('quantinemo_mean.txt',h=TRUE)
                vars=read.table('quantinemo_var.txt',h=TRUE)
                plotCI(x=means$generation,
                       y=means$n.adlt.rs,
                       uiw=0.5*(qnorm(0.975)*sqrt(vars$n.adlt.rs))/sqrt(unique(means$alive.rpl)),
                       sfrac=0.001,
                       pch=21,
                       pt.bg=par("bg"),
                       ylim=c(0,20),
                       xlab="Generations",
                       ylab="Mean number of alleles per locus",
                       main=paste0('N',N[i],'_mig',m[j]))
                abline(h=means[means$generation==2000,]$n.adlt.rs,lty=2,col='red')
        }
}
# Gen 0 = generation 10 car il commence a gen 10
# pas de temps = 10 gen
# pq on est pas à 1? y aura tjrs une petite mutation (ici = 0.00005) equilibre mutation/dérive
# attention un taux de migration  de 0.5 c'est le max entre 2 pops (simulation), 3 pops = 0.33 etc...


for (i in c(1:length(N))){
  for (j in c(1:length(m))){
                setwd(paste0(main,paste0('N',N[i],'_mig',m[j])))    
                means=read.table('quantinemo_mean.txt',h=TRUE)
                vars=read.table('quantinemo_var.txt',h=TRUE)
                plotCI(x=means$generation,
                       y=means$n.adlt.fst.wc_p1.2,
                       uiw=0.5*(qnorm(0.975)*sqrt(vars$n.adlt.fst.wc_p1.2))/sqrt(unique(means$alive.rpl)),
                       sfrac=0.001,
                       pch=21,
                       pt.bg=par("bg"),
                       ylim=c(0,1),
                       xlab="Generations",
                       ylab="Pairwise FST",
                       main=paste0('N',N[i],'_mig',m[j]))
                abline(h=means[means$generation==2000,]$n.adlt.fst.wc_p1.2,lty=2,col='red')
        }
}

# pq pas = 1? pcq en proba on peut pas avoir deux populations absolus ou y'a rien de pareil (hasard, mutation)
# ATTENTION, une mesure de Fst n'est pas une mesure de flux de gènes
# dans la nature on ne connait pas la nature de la valeur du Fst (2 grosses pop ou bien flux de gènes )
# et on connait pas la taille de pop (Fst, on fait l'hypothèse que ce sont les mm tailles de pop)

# Fst de 0 = population sont indiférenciable quasiment (fréq all sont quasi les mêmes entre 2 pops)
