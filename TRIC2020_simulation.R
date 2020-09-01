set.seed(47)
library('MASS') 
library('rethinking')
library('ggplot2')
library('reshape2')

## fitting function 
fitMAP <- function(data){
  m <- map(
    alist(
      Y ~ dnorm(mu ,sigma),
      mu <- a + b * X,
      a ~ dnorm(0,2),
      b ~ dnorm(1,2),
      sigma ~ dunif(0,10)
    ), 
    data=data
  )
  post <- extract.samples(m)
  beta <- coef(m)['b']
  upper <- HPDI(post[,2])['0.89|']
  lower <- HPDI(post[,2])['|0.89']
  return(list(beta, upper, lower))
}

## simulatetion parameters
n <- 110
n_samples <- 40
rhos <- rvals <- seq(0,.99,0.01)

## data storage 
coefs  <- array(dim= c(n, ncol=length(rhos)))
ups <-coefs
lows  <- coefs

## simulate
for (i in seq(1,n)) {
    for (k in seq(1,length(rhos))) {
      rho <- rhos[k]
      df <- data.frame(mvrnorm(n = n_samples, mu = c(0,0),Sigma=matrix(c(1, rho, rho, 1), nrow=2), empirical=TRUE))
      names(df) <- c('Y','X')
      #df$X <- scale(df$X)
      estimates <- try(fitMAP(df))
      estimates
      if (class(estimates) == "try-error"){
        coefs[i,k] <- NA
        ups[i,k] <- NA
        lows[i,k] <- NA
      } else {
        coefs[i,k] <- estimates[[1]]
        ups[i,k] <- estimates[[2]]
        lows[i,k] <- estimates[[3]]
    }
    }}

write.csv(coefs, 'coefs.csv')
write.csv(ups, 'ups.csv')
write.csv(lows, 'lows.csv')

## remove NA's and average sims

Calc_average <- function(data, finalr=100){
  data <- data[complete.cases(data),]
  data <- data[1:finalr,]
  m_data <- apply(data, 2, mean)
  return(m_data)
}

m_beta <- Calc_average(coefs)
m_lows <- Calc_average(lows)
m_ups <- Calc_average(ups)

minval <- rhos[which.min(abs(m_lows))]

mat <- data.frame(cbind(rhos, m_lows, m_ups))

## plot results 
ggplot(mat, aes(rhos)) + 
  geom_ribbon(aes(ymin = m_lows, ymax = m_ups, alpha =.85), fill = "grey70") + 
  scale_x_continuous(breaks=seq(0, 1.2, 0.25))  + 
  scale_y_continuous(breaks = seq(-0.4, 1, 0.2)) +
  theme(axis.line = element_line(colour = "black"), text = element_text(size=16), legend.position = "none") +
  xlab('Pearson Correlation r') + ylab('89 % Credibility Interval of Beta') +
  geom_path(aes(rhos,m_lows)) +  geom_path(aes(rhos,m_ups)) +
  geom_segment(aes(x = 0, y = 0, xend = minval, yend = 0), colour = 'blue', linetype=2) +
  geom_segment(aes(x = minval, y = 0, xend = minval, yend = -.3),colour = 'blue', linetype=2) +
  geom_segment(aes(x = minval+0.2, y = -0.08, xend = minval+0.02 , yend = -0.01), 
               lineend = 'round', linejoin = 'round',
               size = 0.5, arrow = arrow(length = unit(0.1, "inches"))) +
  annotate('text' , x = 0.46, y = -0.1, label = paste('Pearson r =', round(minval,3)),size = 5, hjust = 0 )  +
  coord_cartesian(xlim = c(0.0, 1), ylim = c(-0.3,1), expand = FALSE  ) 

