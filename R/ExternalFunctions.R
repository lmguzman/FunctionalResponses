

OneTreatmentModel <- function(offered, eaten){

  #packages required
  require(bbmle)

  #Making sure data has at least 2 columns
  #if(nrow(offered) != nrow(eaten)) stop ("Prey offered and eaten dont have the same lenght")

  #We are organizing the data
  dd.ml = list(T = 4, P = 1, Initial=offered,
               Killed=eaten)

  # we set up some plausible starting values with upper and lower limits
  start.vals = list(b = 0.004, q = 1.6, h = 0.5)
  lower.lims = c(b = 1e-10, q = -2, h = 1e-10)
  upper.lims = c(b = 1, q = 3, h = 3)

  #Using maximum likelyhood and NLL.oneT
  full.mod <- mle2(NLL.oneT, start=start.vals, data = dd.ml,
                   method="L-BFGS-B",
                   lower=lower.lims,
                   upper=upper.lims)
  return(summary(full.mod))
}



simple1<- dd$offered
simple2<-dd$ne.int.0


OneTreatmentModel(simple1,simple2)


