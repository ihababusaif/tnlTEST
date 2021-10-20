
#' @import partitions
#' @import mathjaxr
tnl<-function(n,l){
  if (any(n < (2*l+1)))
    stop(paste("n must be > 2l", "\n",""))
  x<-NULL
  for(j in 0:n) {x<-rbind(x,t(partitions::compositions(j,n)))}
  ss<-0; prob<-NULL
   for(v in 1:n) {
    ss<-ss+1
    nn<-nrow(x)
    count<-0
    for (kk in 1:nn) {
      zz2<-NULL
      m<-x[kk,]
      for(i in 1:l) {if(sum(m[1:(i+l)])>=i) zz2[i]<-1 else zz2[i]<-0}
      for(i in (l+1):(n-l)) {
        if(sum(m[1:(i-l)])<i&sum(m[1:(i+l)])>=i) zz2[i]<-1 else zz2[i]<-0}
      for(i in (n-l+1):n) {if(sum(m[1:(i-l)])<i) zz2[i]<-1 else zz2[i]<-0}
      if(sum(zz2)==v) count<-count+1
    }
    prob[ss]<-count
  }
  prob<-prob/choose(2*n,n)
  res<-NULL
  for (t in 1:n) {
    res[t]<-sum(prob[1:t])
  }
  result<-list(method="exact",pmf=prob,cdf=res)
  return(result)
}

#' @import plyr
#' @import stats
tnl.sim<-function(n,l,trial=100000){
  if (any(n < (2*l+1)))
    stop(paste("n must be > 2l", "\n",""))
  stest<- function(a,b,l) {
    a<-sort(a)
    b<-sort(b)
    t<-0
    for(i in 1:l)         {if (b[i]<a[i+l]) t<-t+1}
    for(i in (l+1):(n-l)) {if (b[i]>=a[i-l]&b[i]<a[i+l]) t<-t+1}
    for(i in (n-l+1):n)   {if (b[i]>=a[i-l]) t<-t+1}
    return(t)
  }
  x<-y<-NULL
  statistic<-NULL
  for(i in 1:trial){
    x<-rnorm(n);  y<-rnorm(n)
    statistic[i]<-stest(x,y,l)
  }
  statistict<-(plyr::count(statistic)/trial)$freq
  statistict<-c(rep(0,(l-1)),statistict)
  res<-NULL
  for (t in 1:n) {
    res[t]<-sum(statistict[1:t])
  }
  result<-list(method="Monte Carlo simulation",pmf=statistict,cdf=res)
  return(result)
}


#' Non-parametric tests for the two-sample problem based on order statistics and power comparisons
#' @export
#' @rdname tnl.test
#' @param x the first (non-empty) numeric vector of data values.
#' @param y the second (non-empty) numeric vector of data values.
#' @param l class parameter of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}.
#' @param exact the method that will be used. "NULL" or a logical indicating whether an exact should be computed. See ‘Details’ for the meaning of NULL.
#' @description  \code{\link{tnl.test}} performs a nonparametric test for two sample test on vectors of data.
#' @return
#' \subsection{\code{\link{tnl.test}} returns a list with the following components}{
#'    \describe{
#'      \item{\code{statistic:}}{the value of the test statistic.}
#'      \item{\code{p.value:}}{the p-value of the test.}
#'    }
#' }
#' @details A non-parametric two-sample test is performed for testing null hypothesis \loadmathjax \mjeqn{H_0: F = G}{ascii} against the alternative hypothesis \loadmathjax \mjeqn{H_1:F \not= G}{ascii}. The assumptions of the \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}  test are that both samples should come from a continuous distribution and the samples should have the same sample size.\cr
#'      Missing values are silently omitted from x and y.\cr
#'      Exact and simulated p-values are available for the \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} test. If exact ="NULL" (the default) the p-value is computed based on exact distribution when the sample size is less than 11. Otherwise, p-value is computed based on a Monte Carlo simulation. If exact ="TRUE", an exact p-value is computed. If exact="FALSE", a Monte Carlo simulation is performed to compute the p-value. It is recommended to calculate the p-value by a Monte Carlo simulation (use exact="FALSE"), as it takes too long to calculate the exact p-value when the sample size is greater than 10. \cr
#'      The probability mass function (pmf), cumulative density function (cdf) and quantile function of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} are also available in this package, and the above-mentioned conditions about exact ="NULL", exact ="TRUE" and exact="FALSE" is also valid for these functions.\cr
#'      Exact distribution of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} test is also computed under Lehman alternative.\cr
#'      Random number generator of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} test statistic are provided under null hypothesis in the library.
#' @references Karakaya K. et al. (2021). *A Class of Non-parametric Tests for the Two-Sample Problem based on Order Statistics and Power Comparisons*. Submitted paper.\cr
#'             Aliev F. et al. (2021). *A Nonparametric Test for the Two-Sample Problem based on Order Statistics*. Submitted paper.
#' @examples
#' require(stats)
#' x=rnorm(7,2,0.5)
#' y=rnorm(7,0,1)
#' tnl.test(x,y,l=2)
#' #$statistic
#' #[1] 2
#' #
#' #$p.value
#' #[1] 0.02447552

tnl.test<-function(x,y,l,exact="NULL"){
  if( any(is.na(x)) ) warning('Since the data should not contain missing values, we exclude the missing values from the data')
  x=x[!is.na(x)]
  y=y[!is.na(y)]
  if (length(x)!=length(y))
    stop(paste("The length of x and y must be equal", "\n",""))
  n<-length(x)
  stest<- function(a,b,l) {
    a<-sort(a)
    b<-sort(b)
    t<-0
    for(i in 1:l)         {if (b[i]<a[i+l]) t<-t+1}
    for(i in (l+1):(n-l)) {if (b[i]>=a[i-l]&b[i]<a[i+l]) t<-t+1}
    for(i in (n-l+1):n)   {if (b[i]>=a[i-l]) t<-t+1}
    return(t)
  }
  stat<-stest(x,y,l)
  if(exact=="TRUE") {p.value<-tnl(n,l)$cdf[stat]}
  if(exact=="FALSE") {p.value<-tnl.sim(n,l)$cdf[stat]}
  if (exact=="NULL" & n<=10){ p.value<-tnl(n,l)$cdf[stat]}
  if (exact=="NULL" & n>10)  {p.value<-tnl.sim(n,l)$cdf[stat]}
  result<-list(statistic=stat,p.value=p.value)
  return(result)
}


#' Distribution function of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} against the specified quantiles
#' @export
#' @rdname tnl.test
#' @param k,q vector of quantiles.
#' @param n sample size.
#' @param l class parameter of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}
#' @param exact the method that will be used. "NULL" or a logical indicating whether an exact should be computed. See ‘Details’ for the meaning of NULL.
#' @param trial number of trials for simulation.
#' @description  \code{\link{ptnl}} gives the distribution function of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} against the specified quantiles.
#' @return
#' \subsection{\code{\link{ptnl}} returns a list with the following components}{
#'    \describe{
#'      \item{\code{method}:}{The method that was used (exact or simulation). See ‘Details’.}
#'      \item{\code{cdf}:}{distribution function of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} against the specified quantiles.}
#'    }
#' }
#' @examples
#' ##############
#' ptnl(q=2,n=6,l=2,trial = 100000)
#' #$method
#' #[1] "exact"
#' #
#' #$cdf
#' #[1] 0.03030303
#'
ptnl<-function(q,n,l,exact="NULL",trial = 100000){
  if (any(q < l))
    stop(paste("q must be >= l", "\n",""))
  if (any(n < (2*l+1)))
    stop(paste("n must be > 2l", "\n",""))
  if(exact=="TRUE") {ptnl<-tnl(n,l)$cdf[q];method<-"exact"}
  if(exact=="FALSE") {ptnl<-tnl.sim(n,l,trial)$cdf[q];method<-"Monte Carlo simulation"}
  if (exact=="NULL" & n<=10){ ptnl<-tnl(n,l)$cdf[q];method<-"exact"}
  if (exact=="NULL" & n>10)  {ptnl<-tnl.sim(n,l,trial)$cdf[q];method<-"Monte Carlo simulation"}
  result<-list(method=method,cdf=ptnl)
  return(result)
}

#' Density of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} against the specified quantiles
#' @export
#' @rdname tnl.test
#' @param k,q vector of quantiles.
#' @param n sample size.
#' @param l class parameter of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}
#' @param trial number of trials for simulation.
#' @param exact the method that will be used. "NULL" or a logical indicating whether an exact should be computed. See ‘Details’ for the meaning of NULL.
#' @description  \code{\link{dtnl}} gives the density of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} against the specified quantiles.
#' @return
#' \subsection{\code{\link{dtnl}} returns a list with the following components}{
#'    \describe{
#'      \item{\code{method}:}{The method that was used (exact or simulation). See ‘Details’.}
#'      \item{\code{pmf}:}{density of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} against the specified quantiles.}
#'    }
#' }
#' @examples
#' ##############
#' dtnl(k=3,n=7,l=2,trial = 100000)
#' #$method
#' #[1] "exact"
#' #
#' #$pmf
#' #[1] 0.05710956
dtnl<-function(k,n,l,exact="NULL",trial=100000){
  if (any(k < l))
    stop(paste("k must be >= l", "\n",""))
  if (any(n < (2*l+1)))
    stop(paste("n must be > 2l", "\n",""))

  if(exact=="TRUE") {dtnl<-tnl(n,l)$pmf[k];method<-"exact"}
  if(exact=="FALSE") {dtnl<-tnl.sim(n,l,trial)$pmf[k];method<-"Monte Carlo simulation"}
  if (exact=="NULL" & n<=10){ dtnl<-tnl(n,l)$pmf[k];method<-"exact"}
  if (exact=="NULL" & n>10)  {dtnl<-tnl.sim(n,l,trial)$pmf[k];method<-"Monte Carlo simulation"}

  result<-list(method=method,pmf=dtnl)
  return(result)
}


#' quantile function of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}
#' @export
#' @rdname tnl.test
#' @param p vector of probabilities.
#' @param n sample size.
#' @param l class parameter of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}.
#' @param exact the method that will be used. "NULL" or a logical indicating whether an exact should be computed. See ‘Details’ for the meaning of NULL.
#' @param trial number of trials for simulation.
#' @description  \code{\link{qtnl}} gives the quantile function of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} against the specified probabilities.
#' @return
#' \subsection{\code{\link{qtnl}} returns a list with the following components}{
#'    \describe{
#'      \item{\code{method}:}{The method that was used (exact or simulation). See ‘Details’.}
#'      \item{\code{quantile}:}{quantile function  against the specified probabilities.}
#'    }
#' }
#' @examples
#' \dontrun{
#'  qtnl(.3,4,1,exact="FALSE")
#' #$method
#' #[1] "Monte Carlo simulation"
#' #
#' #$quantile
#' #[1] 2
#' }
qtnl<-function(p,n,l,exact="NULL",trial=100000){
  if (any(n < (2*l+1)))
    stop(paste("n must be > 2l", "\n",""))
  if (any(p < 0) | any(p > 1))
    stop(paste("p must be between 0 and 1", "\n", ""))
  if(exact=="TRUE") {cdf=tnl(n,l)$cdf;method<-"exact"}
  if(exact=="FALSE") {cdf=tnl.sim(n,l,trial)$cdf;method<-"Monte Carlo simulation"}
  if (exact=="NULL" & n<=10){ cdf=tnl(n,l)$cdf;method<-"exact"}
  if (exact=="NULL" & n>10)  {cdf=tnl.sim(n,l,trial)$cdf;method<-"Monte Carlo simulation"}
  for (i in 1:n) {
    if(cdf[i]>p) {break}
  }
q=i
result<-list(method=method,quantile=q)
return(result)
  }


#' Random generation for the \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}
#' @export
#' @rdname tnl.test
#' @param N number of observations. If length(N) > 1, the length is taken to be the number required.
#' @param n sample size.
#' @param l class parameter of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}
#' @description \code{\link{rtnl}} generates random values from \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}.
#' @return \code{\link{rtnl}} return *N* of the generated random values.
#' @examples
#' ##############
#' rtnl(N=15,n=7,l=2)
#' #[1] 6 7 6 6 3 5 7 7 7 7 7 6 5 7 4
rtnl<-function(N,n,l){
  stest<- function(a,b,l) {
    a<-sort(a)
    b<-sort(b)
    t<-0
    for(i in 1:l)         {if (b[i]<a[i+l]) t<-t+1}
    for(i in (l+1):(n-l)) {if (b[i]>=a[i-l]&b[i]<a[i+l]) t<-t+1}
    for(i in (n-l+1):n)   {if (b[i]>=a[i-l]) t<-t+1}
    return(t)
  }
  x<-y<-NULL
  statistic<-NULL
  for(i in 1:N){
    x<-rnorm(n);  y<-rnorm(n)
    statistic[i]<-stest(x,y,l)
  }
  return(statistic)
}

#' The density of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} under Lehmann alternatives.
#' @export
#' @rdname tnl.test
#' @param l class parameter of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}
#' @param n sample size.
#' @param gamma parameter of Lehmann alternative
#' @import partitions
#' @description \code{\link{dtnl.lehmann}}  gives the density of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} under Lehmann alternatives.
#' @return \code{\link{dtnl.lehmann}} return vector of the density under Lehmann alternatives against the specified gamma.
#' @examples
#' ##############
#' dtnl.lehmann(l=2,n=6,gamma=0.8)
#' #[1] 0.00000000 0.03757203 0.08230829 0.14229514 0.23276690 0.50505764
dtnl.lehmann<-function(l,n,gamma){
  if (any(n < (2*l+1)))
    stop(paste("n must be > 2l", "\n",""))
  x<-NULL
  for(j in 0:n) {x<-rbind(x,t(partitions::compositions(j,n)))}
  ss<-0; prob<-lehmann<-NULL
  for(v in 1:n) {
    ss<-ss+1
    nn<-nrow(x)
    count<-0
    for (kk in 1:nn) {
      zz2<-NULL
      m<-x[kk,]
      for(i in 1:l) {if(sum(m[1:(i+l)])>=i) zz2[i]<-1 else zz2[i]<-0}
      for(i in (l+1):(n-l)) {
        if(sum(m[1:(i-l)])<i&sum(m[1:(i+l)])>=i) zz2[i]<-1 else zz2[i]<-0}
      for(i in (n-l+1):n) {if(sum(m[1:(i-l)])<i) zz2[i]<-1 else zz2[i]<-0}
      if(sum(zz2)==v) {
        plam <- 1
        for (jj in 1:(n-1)){
          plam <- plam*gamma(sum(m[1:jj]) + jj*gamma)/gamma(sum(m[1:(jj+1)]) + jj*gamma + 1)
        }
        plam <- plam*gamma(sum(m) + n*gamma)/gamma(n + n*gamma + 1)
        const <- (factorial(n)*factorial(n)*(gamma^n))/factorial(m[1])
        lehm<- plam* const
        count<-count+lehm
      }    }
    lehmann[v]<-count
  }
  return(lehmann)
}

#' The distribution function of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} under Lehmann alternatives.
#' @export
#' @rdname tnl.test
#' @param l class parameter of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii}
#' @param n sample size.
#' @param gamma parameter of Lehmann alternative
#' @import partitions
#' @description \code{\link{ptnl.lehmann}}  gives the  distribution function of \loadmathjax \mjeqn{T_n^{(\ell)}}{ascii} under Lehmann alternatives.
#' @return \code{\link{ptnl.lehmann}} return vector of the distribution under Lehmann alternatives against the specified gamma.
#' @examples
#' ##############
#' ptnl.lehmann(l=2,5,gamma=1.2)
#' #[1] 0.0000000 0.0444164 0.1529147 0.3694389 1.0000000
ptnl.lehmann<-function(l,n,gamma){
  if (any(n < (2*l+1)))
    stop(paste("n must be > 2l", "\n",""))
  x<-NULL
  for(j in 0:n) {x<-rbind(x,t(partitions::compositions(j,n)))}
  ss<-0; prob<-lehmann<-NULL
  for(v in 1:n) {
    ss<-ss+1
    nn<-nrow(x)
    count<-0
    for (kk in 1:nn) {
      zz2<-NULL
      m<-x[kk,]
      for(i in 1:l) {if(sum(m[1:(i+l)])>=i) zz2[i]<-1 else zz2[i]<-0}
      for(i in (l+1):(n-l)) {
        if(sum(m[1:(i-l)])<i&sum(m[1:(i+l)])>=i) zz2[i]<-1 else zz2[i]<-0}
      for(i in (n-l+1):n) {if(sum(m[1:(i-l)])<i) zz2[i]<-1 else zz2[i]<-0}
      if(sum(zz2)==v) {
        plam <- 1
        for (jj in 1:(n-1)){
          plam <- plam*gamma(sum(m[1:jj]) + jj*gamma)/gamma(sum(m[1:(jj+1)]) + jj*gamma + 1)
        }
        plam <- plam*gamma(sum(m) + n*gamma)/gamma(n + n*gamma + 1)
        const <- (factorial(n)*factorial(n)*(gamma^n))/factorial(m[1])
        lehm<- plam* const
        count<-count+lehm
      }    }
    lehmann[v]<-count
    res<-NULL
    for (t in 1:n) {
      res[t]<-sum(lehmann[1:t])
    }  }
  return(res)
}
