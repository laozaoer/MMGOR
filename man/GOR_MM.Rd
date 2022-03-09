\name{GOR_MM}
\alias{GOR_MM}
\title{Estimating parameters of the generalized odds rate model with clustered current status data.}
\description{
GOR_MM is used to analyze clustered current status data by fitting the generalized odds rate model. This is the main function of this package. It is assumed that a clustered dataset has \eqn{n} units and there are \eqn{n_i} subjects within the \eqn{i}th unit.  Let \eqn{N=\sum_in_i} be the total number of subject.  The response and explanatory variables are measured for every subjects. Response consists of \eqn{\Delta}, the censoring indicator, and \eqn{C}, the inspection time.   It is assumed that there is also a set of unit level covariates.  
}
\usage{
GOR_MM(Delta,X,Z,n,ni,r,C,knotsnum,order,cluster.ind=TRUE,
pen.ind=FALSE,lambda = 0,itermax=500,tol=1e-7,quadnum=30)

}
\arguments{
  

  \item{Delta}{
        A list of \eqn{n} components, the \eqn{i}th component of the list is an \eqn{n_i\times 1} vector containing  indicators for censoring. The values 1 and 0 indicate left censored and right censored, respectively.}
        
  \item{X}{
  A list of \eqn{n} components, the \eqn{i}th component is an \eqn{n_i\times p} matrix of  subject-level covariates. All columns of each matrix is assumed to be binary or numeric. For including a nominal categorical variable, the corresponding dummy variables must be included in this matrix. If there is no subject-level covariate (that means \eqn{p=0}), then set \eqn{X=}NULL.

  }
  \item{Z}{
  An \eqn{n\times q} matrix of unit-level (cluster-level) covariates. All columns of this matrix is assumed to be binary or numeric. For including a nominal categorical variable, the corresponding dummy variables must be included in this matrix. If there is no unit-level covariate (that means \eqn{q=0}), then set \eqn{Z=}NULL.
}
  \item{n}{An integer that denotes the number of units (clusters). }
  \item{ni}{
  An \eqn{n\times 1} vector of integers. The sum of the elements of this vector equals to \eqn{N}. Each entry of this vector must be greater than one. If all \eqn{n} entries are one, then it is no longer a clustered data, and cluster.ind is set to FALSE.}
  \item{r}{
  A scalar that specifies the  generalized odds rate model. The value of \eqn{r} must be non-negative. In this function \eqn{r=0} denotes the GOR model with \eqn{r\rightarrow 0}. }
 \item{C}{
 A list of \eqn{n+1} components, the \eqn{i}th component of the first \eqn{n} components is a \eqn{n_i\times 1} vector containing the inspection time of all \eqn{n_i} subjects within the \eqn{i}th unit. The inspection time must be positive for all subjects. The \eqn{(n+1)}th component of \eqn{C} is the lower and upper bounds of inspection times. The bounds should not be reached.
 
 }
    \item{knotsnum}{
    An integer indicating the number of equidistant interior knots for the integrated B-spline approximation of the nonparametric component of the GOR model. 
   }
   \item{order}{
   An integer indicating the  degree of integrated B-splines. }
   \item{cluster.ind}{
   Logical. If TRUE, the parameters will be estimated along with clustering effect. If it is FALSE, then all subjects are treated as independent and there is no clustering effect.  The default value is TRUE.
   }
 \item{pen.ind}{
 Logical. If FALSE, there is no penalty term. The default value is FALSE. If TRUE, users should specify a  positive value for the penalty term lambda.
 }
 \item{lambda}{
 A scalar containing the penalty parameter. The default value is 0.
 }
  \item{itermax}{
  The maximum number of iterations. The default value is 500.
  }
   \item{tol}{
   This denotes the sum of the absolute relative differences of the estimates at two
successive iterations. It is used to check the convergence of the parameter estimates.
The default value is 1e-7.

   }

\item{quadnum}{
The number of Gauss-Hermite quadrature nodes used in numerical integration. The default value is 30.
}
}

\value{
Function GOR_MM returns a list containing the following components:


\item{result}{The first column contains the estimate of the regression parameters. The second column contains the standard errors of these estimators.}
\item{AIC}{
The AIC value.
}


}

\references{
Wang, T., He, K., Ma, W., Bandyopadhyay, D., and Sinha, S. (2022). Minorize-Maximize Algorithm  for the Generalized Odds Rate Model for Clustered Current Status Data. \emph{Canadian Journal of Statistics}, in press.
}


\examples{

library(MMGOR)
set.seed(1)
H=function(t) log(1+t)+t^(3/2) # specify the H function

#This function generates a cohort with 300 subjects, 
#where each subject has up to 8 observations. 
#Both the subject-level covariates and unit-level covariates
# are generated from Uniform(-1,1) distribution. 
#Users can use their own H function.

#Scenario 1
#3 subject-level covariates, 2 unit-level covariates (p=3,q=2)
data=SIM_DATA(r=0,beta=c(-1,1,0.5),gamma=c(-1,1),theta=1,n=300,H=H,up.quantile=0.85,
par_size=c(exp(1.7),1,8))
#estimate the parameters
result=GOR_MM(data$Delta,data$X,data$Z,data$n,data$ni,r=data$r,data$C,knotsnum = 2,order=2)

#Scenario 2
#only 2 unit-level covariates (p=0,q=2)
data=SIM_DATA(r=0,beta=c(0),gamma=c(-1,1),theta=1,n=300,H=H,up.quantile=0.85,
par_size=c(exp(1.7),1,8))
#estimate the parameters
result=GOR_MM(data$Delta,X=NULL,data$Z,data$n,data$ni,r=data$r,data$C,knotsnum = 2,order=2)

#Scenario 3
#only 3 subject-level covariates (p=3,q=0)
data=SIM_DATA(r=0,beta=c(-1,1,0.5),gamma=c(0),theta=1,n=300,H=H,up.quantile=0.85,
par_size=c(exp(1.7),1,8))
result=GOR_MM(data$Delta,data$X,Z=NULL,data$n,data$ni,r=data$r,data$C,knotsnum = 2,order=2)

#Scenario 4
#3 subject-level covariates, 2 unit-level covariates (p=3,q=2)
#no clustering effect
data=SIM_DATA(r=0,beta=c(-1,1,0.5),gamma=c(-1,1),theta=0,n=300,H=H,up.quantile=0.85,
par_size=c(exp(1.7),1,8))
result=GOR_MM(data$Delta,data$X,data$Z,data$n,data$ni,r=data$r,data$C,cluster.ind = FALSE,
knotsnum = 2,order=2)

# Scenario 5
# Data generation without the function SIM_DATA
n=100
ni= rpois(n, 4)
ni[ni<1]<-1
ni[ni>10]<-10
p=2;
q=1;
gamma=rep(-0.5, q)
beta=rep(0.5, p)
theta=0.8
xcov<-list()
zcov<-c()
mytime<-list()
mydelta<-list()
myinspection<-list()
zcov<-matrix(runif(n*q),nrow = n,ncol = q)

for( j in 1:n){  
  
  xcov[[j]]<- matrix( runif((p*ni[j]), -1, 1), ncol=p )  
  
  myrate= exp(as.numeric(zcov[j]%*%gamma)+
                as.numeric(xcov[[j]]%*%beta)-3+theta*rnorm(1))
  mytime[[j]]<-40+ rexp(ni[j],myrate)
  myinspection[[j]]<- runif(ni[j], 40, 90)
  mydelta[[j]]<- as.numeric(mytime[[j]]<myinspection[[j]] )
}
myinspection[[n+1]]=c(40,90)
# End of data generation

### Analysis of the data using the GOR_MM function 
result=GOR_MM(mydelta,xcov,zcov,n,ni,r=0,
              myinspection,knotsnum = 2,order=2)



}