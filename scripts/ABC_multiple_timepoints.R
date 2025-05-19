library(abc)
library(spatstat)
library(weights)

setwd("//wsl.localhost/Ubuntu/home/vivak/nfds/")
demog = "severe_contraction"

#stats without filtering
t_abc <- head(read.table(file.path("ABC", paste("readyForABC_", demog, ".txt", sep=''))
                         , h=T, sep='\t'), 3000)
t_par <- head(read.table("params/fp_concat.txt",
                         sep='\t', col.names = c('param','Tb','Feq')), 3000)
#get time and make stats:
t_par <- cbind(t_par[,c(2:3)])
t_sim_stats <- cbind(t_abc[,c(2:181)])

t_emp_stats <- read.table(file.path("empirical_sims", paste(demog, "_ABC.txt", sep=''))
                          , h=T)

cv_ridge <- cv4abc(t_par, t_sim_stats, nval=100, tols=c(0.05, 0.08, 0.1), statistic="median", method="ridge", transf="none")
summary(cv_ridge)
plot(cv_ridge)

cv_nnet <- cv4abc(t_par, t_sim_stats,nval=100, tols=c(0.05, 0.08, 0.1), statistic="median", method="neuralnet", transf="none")
summary(cv_nnet)
plot(cv_nnet)

rej <- abc(target=t_emp_stats, param=t_par, sumstat=t_sim_stats, tol=0.08, method="rejection")
summary(rej)
df <- rej$unadj.values
write.table(df, "//wsl.localhost/Ubuntu/home/vivak/nfds/posteriors/eq_0.1_rej.txt", 
            sep='\t', quote=FALSE, row.names=FALSE)
hist(rej, breaks=10)

nnet <- abc(target=t_emp_stats, param=t_par, sumstat=t_sim_stats, tol=0.05, 
            method="neuralnet")
summary(nnet)


df <- nnet$adj.values


par(mfrow=c(2,2))
m_adj_values <- nnet$adj.values
m_weights <- nnet$weights

#Inference performed multiple times and then taking a mean:
num_reps <- 50
Tb <- c() 
Feq <- c()

i = 1
while(i<=num_reps){
  print (i)
  nnet <- abc(target=t_emp_stats, param=t_par, sumstat=t_sim_stats, tol=0.05, method="neuralnet")
  if (i==1){
    m_adj_values <- nnet$adj.values
    m_weights <- nnet$weights
  }
  else{
    m_adj_values <- rbind(m_adj_values, nnet$adj.values)
    m_weights <- c(m_weights, nnet$weights)
  }
  #calculating errors:
  est1 <- weighted.median(nnet$adj.values[,1], nnet$weights, na.rm = TRUE)
  est2 <- weighted.median(nnet$adj.values[,2], nnet$weights, na.rm = TRUE)
  
  if (est1 < 0){Tb <- c(Tb, 0.0)}
  else{Tb <- c(Tb, est1)}
  if (est2 < 0){Feq <- c(Feq, 0.0)}
  else{Feq <- c(Feq, est2)}
  i = i + 1
}

Tb_pos <- wtd.hist(x=m_adj_values[,1], weight=m_weights, plot=F)
Feq_pos <- wtd.hist(x=m_adj_values[,2], weight=m_weights, plot=F)

write.table(m_adj_values, file.path("//wsl.localhost/Ubuntu/home/vivak/nfds/posteriors/nnet", 
                                    paste(demog, ".txt", sep='')), 
            sep='\t', quote=FALSE, row.names=FALSE)

df <- data.frame(head(Tb_pos$breaks,-1), Tb_pos$density)
colnames(df) <- c("Tb", "density")
write.table(df, "//wsl.localhost/Ubuntu/home/vivak/nfds/posteriors/nnet/Tb_severe_contraction.txt", 
            sep='\t', quote=FALSE, row.names=FALSE)

df <- data.frame(head(Feq_pos$breaks,-1), Feq_pos$density)
colnames(df) <- c("Feq", "density")
write.table(df, "//wsl.localhost/Ubuntu/home/vivak/nfds/posteriors/nnet/Feq_severe_contraction.txt", 
            sep='\t', quote=FALSE, row.names=FALSE)