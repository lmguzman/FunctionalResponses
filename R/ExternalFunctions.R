
source('InternalFunctions.R')
detach()

 OneTreatmentModel <- function(offered, eaten){

  #packages required
  library(bbmle)

  #Making sure data has at least 2 columns
  #if(nrow(offered) != nrow(eaten)) stop ("Prey offered and eaten dont have the same lenght")

  #We are organizing the data
  dd.ml = list(T = 4, P = 1, Initial=offered,
               Killed=eaten)

  # we set up some plausible starting values with upper and lower limits
  start.vals = list(b = 0.004, q = 1.6, h = 0.5)
  lower.lims = c(b = 1e-10, q = -3, h = 1e-10)
  upper.lims = c(b = 3, q = 3, h = 3)

  #Using maximum likelyhood and NLL.oneT
  full.mod <- mle2(NLL.oneT, start=start.vals, data = dd.ml,
                   method="L-BFGS-B",
                   lower=lower.lims,
                   upper=upper.lims)
  return(summary(full.mod))
}








 FullModel <- function(offered, eaten, treatment){
   #packages required
   library(bbmle)

   #Making sure data has at least 2 columns
   #if(nrow(offered) != nrow(eaten)) stop ("Prey offered and eaten dont have the same lenght")

   #We are organizing the data
   dd.ml = list(T = 4, P = 1, Initial=offered,
                Killed=eaten, ind=ifelse(treatment=="no", 0, 1))
   #Start values and limits
   start.vals = list(b.int=0.14, b.ind=0.004,
                     q.int=0.01, q.ind=1.6,
                     h.int=0.357, h.ind=1.68)
   lower.lims=c(b.int=1e-10, b.ind=1e-10,
                q.int=-2, q.ind=-2,
                h.int=1e-10, h.ind=1e-10)
   upper.lims=c(b.int=1, b.ind=1,
                q.int=3, q.ind=3,
                h.int=3, h.ind=3)

   full.mod <- mle2(NLL.real.full, start=start.vals, data=dd.ml,
                    method="L-BFGS-B",
                    lower=lower.lims,
                    upper=upper.lims)
   Full.Mod<-list('Coefficients' = attributes(summary(full.mod))$coef, 'AIC' = AIC(full.mod))

   return(Full.Mod)
 }




 BConstant <- function(offered, eaten, treatment){
   #packages required
   library(bbmle)

   #Making sure data has at least 2 columns
   #if(nrow(offered) != nrow(eaten)) stop ("Prey offered and eaten dont have the same lenght")

   #We are organizing the data
   dd.ml = list(T = 4, P = 1, Initial=offered,
                Killed=eaten, ind=ifelse(treatment=="no", 0, 1))
   start.vals = list(b=0.14,
                     q.int=0.01, q.ind=1.6,
                     h.int=0.357, h.ind=1.68)
   lower.lims=c(b=1e-10,
                q.int=-2, q.ind=-2,
                h.int=1e-10, h.ind=1e-10)
   upper.lims=c(b=1,
                q.int=3, q.ind=3,
                h.int=3, h.ind=3)
   b.fixed.mod <- mle2(NLL.real.b, start=start.vals, data=dd.ml,
                       method="L-BFGS-B",
                       lower=lower.lims,
                       upper=upper.lims)
   B.Fixed.Mod<-list('Coefficients' = attributes(summary(b.fixed.mod))$coef, 'AIC' = AIC(b.fixed.mod))

   return(B.Fixed.Mod)
 }


 HConstant <- function(offered, eaten, treatment){
   #packages required
   library(bbmle)

   #Making sure data has at least 2 columns
   #if(nrow(offered) != nrow(eaten)) stop ("Prey offered and eaten dont have the same lenght")

   #We are organizing the data
   dd.ml = list(T = 4, P = 1, Initial=offered,
                Killed=eaten, ind=ifelse(treatment=="no", 0, 1))

   start.vals = list(b.int=0.1, b.ind=0.04,
                     q.int=0.01, q.ind=1.6,
                     h=1.3)
   lower.lims=c(b.int=1e-10, b.ind=1e-10,
                q.int=-2, q.ind=-2,
                h=1e-10)
   upper.lims=c(b.int=1, b.ind=1,
                q.int=3, q.ind=3,
                h=3)
   h.fixed.mod <- mle2(NLL.real.h, start=start.vals, data=dd.ml,
                       method="L-BFGS-B",
                       lower=lower.lims,
                       upper=upper.lims)
   H.Fixed.Mod<-list('Coefficients' = attributes(summary(h.fixed.mod))$coef, 'AIC' = AIC(h.fixed.mod))

   return(H.Fixed.Mod)
 }



 QConstant <- function(offered, eaten, treatment){
   #packages required
   library(bbmle)

   #Making sure data has at least 2 columns
   #if(nrow(offered) != nrow(eaten)) stop ("Prey offered and eaten dont have the same lenght")

   #We are organizing the data
   dd.ml = list(T = 4, P = 1, Initial=offered,
                Killed=eaten, ind=ifelse(treatment=="no", 0, 1))
   start.vals = list(b.int=0.14, b.ind=0.004,
                     q=1.6,
                     h.int=0.357, h.ind=1.68)
   lower.lims=c(b.int=1e-10, b.ind=1e-10,
                q=-2,
                h.int=1e-10, h.ind=1e-10)
   upper.lims=c(b.int=1, b.ind=1,
                q=3,
                h.int=3, h.ind=3)
   q.fixed.mod <- mle2(NLL.real.q, start=start.vals, data=dd.ml,
                       method="L-BFGS-B",
                       lower=lower.lims,
                       upper=upper.lims)
   Q.Fixed.Mod<-list('Coefficients' = attributes(summary(q.fixed.mod))$coef, 'AIC' = AIC(q.fixed.mod))

   return(Q.Fixed.Mod)
 }




BHConstant <- function(offered, eaten, treatment){

  #packages required
  library(bbmle)

  #Making sure data has at least 2 columns
  #if(nrow(offered) != nrow(eaten)) stop ("Prey offered and eaten dont have the same lenght")

  #We are organizing the data
  dd.ml = list(T = 4, P = 1, Initial=offered,
               Killed=eaten, ind=ifelse(treatment=="no", 0, 1))

  start.vals = list(b=0.1,
                    q.int=0.01, q.ind=0.6,
                    h=0.35)
  lower.lims=c(b=1e-10,
               q.int=-2, q.ind=-2,
               h=1e-10)
  upper.lims=c(b=1,
               q.int=3, q.ind=3,
               h=3)
  bh.fixed.mod <- mle2(NLL.real.bh, start=start.vals, data=dd.ml,
                       method="L-BFGS-B",
                       lower=lower.lims,
                       upper=upper.lims)
  BH.Fixed.Mod<-list('Coefficients' = attributes(summary(bh.fixed.mod))$coef, 'AIC' = AIC(bh.fixed.mod))

  return(BH.Fixed.Mod)
}



BQConstant <- function(offered, eaten, treatment){

  #packages required
  library(bbmle)

  #Making sure data has at least 2 columns
  #if(nrow(offered) != nrow(eaten)) stop ("Prey offered and eaten dont have the same lenght")

  #We are organizing the data
  dd.ml = list(T = 4, P = 1, Initial=offered,
               Killed=eaten, ind=ifelse(treatment=="no", 0, 1))

  start.vals = list(b=0.14,
                    q=1.6,
                    h.int=0.357, h.ind=1.68)
  lower.lims=c(b=1e-10,
               q=-2,
               h.int=1e-10, h.ind=1e-10)
  upper.lims=c(b=1,
               q=3,
               h.int=3, h.ind=3)
  bq.fixed.mod <- mle2(NLL.real.bq, start=start.vals, data=dd.ml,
                       method="L-BFGS-B",
                       lower=lower.lims,
                       upper=upper.lims)
  BQ.Fixed.Mod<-list('Coefficients' = attributes(summary(bq.fixed.mod))$coef, 'AIC' = AIC(bq.fixed.mod))

  return(BQ.Fixed.Mod)
}


HQConstant <- function(offered, eaten, treatment){
  #packages required
  library(bbmle)

  #Making sure data has at least 2 columns
  #if(nrow(offered) != nrow(eaten)) stop ("Prey offered and eaten dont have the same lenght")

  #We are organizing the data
  dd.ml = list(T = 4, P = 1, Initial=offered,
               Killed=eaten, ind=ifelse(treatment=="no", 0, 1))

  start.vals = list(b.int=0.14, b.ind=0.004,
                    q=1.6,
                    h=1.01)
  lower.lims=c(b.int=1e-10, b.ind=1e-10,
               q=-2,
               h=1e-10)
  upper.lims=c(b.int=1, b.ind=1,
               q=3,
               h=3)
  hq.fixed.mod <- mle2(NLL.real.hq, start=start.vals, data=dd.ml,
                       method="L-BFGS-B",
                       lower=lower.lims,
                       upper=upper.lims)
  HQ.Fixed.Mod<-list('Coefficients' = attributes(summary(hq.fixed.mod))$coef, 'AIC' = AIC(hq.fixed.mod))

  return(HQ.Fixed.Mod)
}



BQHConstant <- function(offered, eaten, treatment){
  #packages required
  library(bbmle)

  #We are organizing the data
  dd.ml = list(T = 4, P = 1, Initial=offered,
               Killed=eaten, ind=ifelse(treatment=="no", 0, 1))
  start.vals = list(b=0.14,
                    q=1.6,
                    h=1.01)
  lower.lims=c(b=1e-10,
               q=-2,
               h=1e-10)
  upper.lims=c(b=1,
               q=3,
               h=3)
  bqh.fixed.mod <- mle2(NLL.real.bhq, start=start.vals, data=dd.ml,
                        method="L-BFGS-B",
                        lower=lower.lims,
                        upper=upper.lims)
  BQH.Fixed.Mod<-list('Coefficients' = attributes(summary(bqh.fixed.mod))$coef, 'AIC' = AIC(bqh.fixed.mod))

  return(BQH.Fixed.Mod)

}





  ## Treament has to be a vector of y and n

TwoTreatmentModel<- function(offered, eaten, treatment, type){
  offered <- offered
  eatern <- eaten
  treatment <- treatment

if(type == "FULL") {
  FullModel(offered, eaten, treatment)
} else if(type == "B") {
    BConstant(offered, eaten, treatment)
}else if(type == "H"){
    HConstant(offered, eaten, treatment)
} else if(type == "Q"){
    QConstant(offered, eaten, treatment)
} else if(type == "BH"){
    BHConstant(offered, eaten, treatment)
} else if(type == "BQ"){
    BQConstant(offered, eaten, treatment)
} else if(type == "HQ"){
    HQConstant(offered, eaten, treatment)
} else if(type == "BQH"){
    BQHConstant(offered, eaten, treatment)
} else if(type == "ALL"){
    FM <- FullModel(offered, eaten, treatment)
    print('full mode done!')

    B <- BConstant(offered, eaten, treatment)
    print('B mode done!')

    H <- HConstant(offered, eaten, treatment)
    print('H mode done!')

    Q <- QConstant(offered, eaten, treatment)
    print('Q mode done!')

    BH <- BHConstant(offered, eaten, treatment)
    print('B and H mode done!')

    BQ <- BQConstant(offered, eaten, treatment)
    print('B and Q mode done!')

    HQ <- HQConstant(offered, eaten, treatment)
    print('H and Q mode done!')

    BQH <- BQHConstant(offered, eaten, treatment)
    print('All done!')

    results <- list('Full Model' = FM, 'B Fixed' = B, 'H Fixed' = H, 'Q Fixed' = Q,
                    'BH Fixed' = BH, 'BQ Fixed' = BQ, 'HQ Fixed' = HQ, 'BQH Fixed' = BQH)
    return(results)
  }else {stop("You need to input a type")}
}







