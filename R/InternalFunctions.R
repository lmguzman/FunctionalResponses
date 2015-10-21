
## Ben Bolker
## April 19, 2012
## Rogers random predator equation: extensions and estimation by numerical integration
## http://ms.mcmaster.ca/~bolker/misc/rogers2.pdf


##       attack
# This is a function that models the attack rate. Input is b, N and q.
# b is what is typically known as the attack rate per capita. N is the density.
#In this case the attack rate is not linear but is dependent on N.
#Parameter q determines the type of the functional response.
#When q = 0, Nq = 1 and the attack rate is the same regardless of prey density,
#so the equation pro- duces a Type II functional response (Hammill et al., 2010b).
#When q > 0, a Type III curve is produced as the attack rate is multiplied by Nq and so increases with prey density (Hammill et al., 2010b).


attack <- function(b,q,N) {
  b*N^q
}

####   grafun
#This function models a type II or III functional response depending on what q is
# THis function does not consider that depletion may occur

gradfun <- function(t,y,parms) {
  with(as.list(c(y,parms)),
       { A  <- attack(b,q,N)
       grad <- -A*N/(1+A*h*N)
       list(grad,NULL)
       })
}

####    fun2
#Using numerical integration, we obtain the values of prey eaten when depletion occurs
#that is, when prey are not being replaced throughout the experiment

fun2 <- function(b,q,h,T, P,N0) {
  library(deSolve)
  L <- lsoda(c(N=N0),
             times=c(0,T),
             func=gradfun,
             parms=c(b=b,q=q,h=h,P=P))
  N.final <- L[2,2]
  N0-N.final
}





#### ONE TREATMENT ######

## Likelyhood function for one treatment scenario
## Initial=dd$offered Initial is a vector with offered inidividuals

NLL.oneT = function(b, q, h, T = 4, P = 1) {
  b <- b
  q <- q
  h <- h
  prop.exp <- numeric(length=length(Initial))
  prop.exp <- sapply(Initial, fun2,
                     b=b, q=q, h=h, P=P, T=T)/Initial
  - sum(dbinom(Killed, prob = prop.exp, size = Initial,
               log = TRUE))
}





##### FULL MODEL TWO TREATMENTS ####

## Likelyhood function for two treatment scenario. Induced treatment == 0
# ind is the induced vector, Initial is a vector with offered inidivuals
# where all variables can change
NLL.real.full <- function(b.int,b.ind,
                          q.int, q.ind,
                          h.int, h.ind,
                          T = 4, P = 1,
                          ind) {

  ## NOTE ### here the parameters are not the same as before!!!
  ## IMPORTANT -- the parameters are the actual values, not the usual for R
  ## This makes constraining the values easier
  b <- ifelse(ind==0, b.int, b.ind)
  q <- ifelse(ind==0, q.int, q.ind)
  h <- ifelse(ind==0, h.int, h.ind)

  print(c(b.int, b.ind, h.int, h.ind, q.int, q.ind))

  prop.exp <- numeric(length=length(Initial))
  prop.exp[ind==0] <- sapply(Initial[ind==0], fun2,
                             b=b.int, q=q.int, h=h.int, P=P, T=T)/Initial[ind==0]
  prop.exp[ind==1] <- sapply(Initial[ind==1], fun2,
                             b=b.ind, q=q.ind, h=h.ind, P=P, T=T)/Initial[ind==1]

  print(sum(prop.exp))

  - sum(dbinom(Killed, prob = prop.exp, size = Initial,
               log = TRUE))
}





###### B IS CONSTANT BETWEEN TREATMENTS#######

## Likelyhood function for two treatment scenario. Induced treatment == 0
# ind is the induced vector, Initial is a vector with offered inidivuals
# ### b is constant between treatments

NLL.real.b <- function(b,
                       q.int, q.ind,
                       h.int, h.ind,
                       T = 4, P = 1,
                       ind) {

  ## NOTE ### here the parameters are not the same as before!!!
  ## IMPORTANT -- the parameters are the actual values, not the usual for R
  ## This makes constraining the values easier
  b <- b
  q <- ifelse(ind==0, q.int, q.ind)
  h <- ifelse(ind==0, h.int, h.ind)

  print(c(b, h.int, h.ind, q.int, q.ind))

  prop.exp <- numeric(length=length(Initial))
  prop.exp[ind==0] <- sapply(Initial[ind==0], fun2,
                             b=b, q=q.int, h=h.int, P=P, T=T)/Initial[ind==0]
  prop.exp[ind==1] <- sapply(Initial[ind==1], fun2,
                             b=b, q=q.ind, h=h.ind, P=P, T=T)/Initial[ind==1]

  print(sum(prop.exp))

  - sum(dbinom(Killed, prob = prop.exp, size = Initial,
               log = TRUE))
}





###### H IS CONSTANT BETWEEN TREATMENTS#######

## Likelyhood function for two treatment scenario. Induced treatment == 0
# ind is the induced vector, Initial is a vector with offered inidivuals
# ### h is constant between treatments

NLL.real.h <- function(b.int,b.ind,
                       q.int, q.ind,
                       h,
                       T = 4, P = 1,
                       ind) {

  ## NOTE ### here the parameters are not the same as before!!!
  ## IMPORTANT -- the parameters are the actual values, not the usual for R
  ## This makes constraining the values easier
  b <- ifelse(ind==0, b.int, b.ind)
  q <- ifelse(ind==0, q.int, q.ind)
  h <- h

  prop.exp <- numeric(length=length(Initial))
  prop.exp[ind==0] <- sapply(Initial[ind==0], fun2,
                             b=b.int, q=q.int, h=h, P=P, T=T)/Initial[ind==0]
  prop.exp[ind==1] <- sapply(Initial[ind==1], fun2,
                             b=b.ind, q=q.ind, h=h, P=P, T=T)/Initial[ind==1]


  - sum(dbinom(Killed, prob = prop.exp, size = Initial,
               log = TRUE))
}





###### Q IS CONSTANT BETWEEN TREATMENTS#######

## Likelyhood function for two treatment scenario. Induced treatment == 0
# ind is the induced vector, Initial is a vector with offered inidivuals
# ### q is constant between treatments

NLL.real.q <- function(b.int,b.ind,
                       q,
                       h.int, h.ind,
                       T = 4, P = 1,
                       ind) {

  ## NOTE ### here the parameters are not the same as before!!!
  ## IMPORTANT -- the parameters are the actual values, not the usual for R
  ## This makes constraining the values easier
  b <- ifelse(ind==0, b.int, b.ind)
  q <- q
  h <- ifelse(ind==0, h.int, h.ind)

  prop.exp <- numeric(length=length(Initial))
  prop.exp[ind==0] <- sapply(Initial[ind==0], fun2,
                             b=b.int, q=q, h=h.int, P=P, T=T)/Initial[ind==0]
  prop.exp[ind==1] <- sapply(Initial[ind==1], fun2,
                             b=b.ind, q=q, h=h.ind, P=P, T=T)/Initial[ind==1]


  - sum(dbinom(Killed, prob = prop.exp, size = Initial,
               log = TRUE))
}





###### B and H IS CONSTANT BETWEEN TREATMENTS#######

## Likelyhood function for two treatment scenario. Induced treatment == 0
# ind is the induced vector, Initial is a vector with offered inidivuals
# ### b and h is constant between treatments


NLL.real.bh <- function(b,
                        q.int, q.ind,
                        h,
                        T = 4, P = 1,
                        ind) {

  ## NOTE ### here the parameters are not the same as before!!!
  ## IMPORTANT -- the parameters are the actual values, not the usual for R
  ## This makes constraining the values easier
  b <- b
  q <- ifelse(ind==0, q.int, q.ind)
  h <- h

  prop.exp <- numeric(length=length(Initial))
  prop.exp[ind==0] <- sapply(Initial[ind==0], fun2,
                             b=b, q=q.int, h=h, P=P, T=T)/Initial[ind==0]
  prop.exp[ind==1] <- sapply(Initial[ind==1], fun2,
                             b=b, q=q.ind, h=h, P=P, T=T)/Initial[ind==1]


  - sum(dbinom(Killed, prob = prop.exp, size = Initial,
               log = TRUE))
}



###### B and Q IS CONSTANT BETWEEN TREATMENTS#######

## Likelyhood function for two treatment scenario. Induced treatment == 0
# ind is the induced vector, Initial is a vector with offered inidivuals
# ### b and q is constant between treatments

NLL.real.bq <- function(b,
                        q,
                        h.int, h.ind,
                        T = 4, P = 1,
                        ind) {

  ## NOTE ### here the parameters are not the same as before!!!
  ## IMPORTANT -- the parameters are the actual values, not the usual for R
  ## This makes constraining the values easier
  b <- b
  q <- q
  h <- ifelse(ind==0, h.int, h.ind)

  prop.exp <- numeric(length=length(Initial))
  prop.exp[ind==0] <- sapply(Initial[ind==0], fun2,
                             b=b, q=q, h=h.int, P=P, T=T)/Initial[ind==0]
  prop.exp[ind==1] <- sapply(Initial[ind==1], fun2,
                             b=b, q=q, h=h.ind, P=P, T=T)/Initial[ind==1]


  - sum(dbinom(Killed, prob = prop.exp, size = Initial,
               log = TRUE))
}


###### H and Q IS CONSTANT BETWEEN TREATMENTS#######

## Likelyhood function for two treatment scenario. Induced treatment == 0
# ind is the induced vector, Initial is a vector with offered inidivuals
# ### H and q is constant between treatments


NLL.real.hq <- function(b.int,b.ind,
                        q,
                        h,
                        T = 4, P = 1,
                        ind) {

  ## NOTE ### here the parameters are not the same as before!!!
  ## IMPORTANT -- the parameters are the actual values, not the usual for R
  ## This makes constraining the values easier
  b <- ifelse(ind==0, b.int, b.ind)
  q <- q
  h <- h

  prop.exp <- numeric(length=length(Initial))
  prop.exp[ind==0] <- sapply(Initial[ind==0], fun2,
                             b=b.int, q=q, h=h, P=P, T=T)/Initial[ind==0]
  prop.exp[ind==1] <- sapply(Initial[ind==1], fun2,
                             b=b.ind, q=q, h=h, P=P, T=T)/Initial[ind==1]


  - sum(dbinom(Killed, prob = prop.exp, size = Initial,
               log = TRUE))
}



### ALL CONTANT###
## Likelyhood function for two treatment scenario. Induced treatment == 0
# ind is the induced vector, Initial is a vector with offered inidivuals
# ### ALL constant between treatments

NLL.real.bhq <- function(b,
                         q,
                         h,
                         T = 4, P = 1,
                         ind) {

  ## NOTE ### here the parameters are not the same as before!!!
  ## IMPORTANT -- the parameters are the actual values, not the usual for R
  ## This makes constraining the values easier
  b <- b
  q <- q
  h <- h

  prop.exp <- numeric(length=length(Initial))
  prop.exp[ind==0] <- sapply(Initial[ind==0], fun2,
                             b=b, q=q, h=h, P=P, T=T)/Initial[ind==0]
  prop.exp[ind==1] <- sapply(Initial[ind==1], fun2,
                             b=b, q=q, h=h, P=P, T=T)/Initial[ind==1]


  - sum(dbinom(Killed, prob = prop.exp, size = Initial,
               log = TRUE))
}

