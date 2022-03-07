\name{GOR_MM}
\alias{GOR_MM}
\title{Estimating parameters of the generalized odds rate model with clustered current status data.}
\description{
GOR_MM is used to analyze clustered current status data by fitting the generalized odds rate model. This is the main function of this package. It is assumed that a clustered dataset has \eqn{n} units and there are \eqn{n_i}subjects within the \eqn{i}th unit.  Let \eqn{N=\sum_in_i} be the total number of subject.  The response and explanatory variables are measured for every subjects. Response consists of \eqn{\Delta}, the censoring indicator, and \eqn{C}, the inspection time.   It is assumed that there is also a set of unit level covariates.  
}
\usage{
GOR_MM(Delta,X,Z,n,ni,r,C,knotsnum,order,cluster.ind=TRUE,
pen.ind=FALSE,lambda = 0,itermax=500,tol=1e-7,quadnum=30)

}
\arguments{
  

  \item{Delta}{
        A list of \eqn{n} components, the \eqn{i}th component of the list is a \eqn{n_i\times 1} vector containing  indicators for censoring. The values 0 and 1 indicate left censored and right censored, respectively.}
        
  \item{X}{
  A list of \eqn{n} components, the \eqn{i}th component is a \eqn{n_i\times p} matrix of  subject-level covariates. All columns of each matrix is assumed to be binary or numeric. For including a nominal categorical variable, the corresponding dummy variables must be included in this matrix. If there is no subject-level covariate (that means \eqn{p=0}), then set \eqn{X=}NULL.

  }
  \item{Z}{
  An \eqn{n\times q} matrix of unit-level (cluster-level) covariates. All columns of this matrix is assumed to be binary or numeric. For including a nominal categorical variable, the corresponding dummy variables must be included in this matrix. If there is no unit-level covariate (that means \eqn{q=0}), then set \eqn{Z=}NULL.
}
  \item{n}{An integer that denotes the number of units (clusters). }
  \item{ni}{
  An \eqn{n\times 1} vector of integers. The sum of the elements of this vector equals to \eqn{N}. }
  \item{r}{
  A scalar that specifies the  generalized odds rate model.}
 \item{C}{
 A list of \eqn{n} components, the \eqn{i}th component is a \eqn{n_i\times 1} vector containing the inspection time of all \eqn{n_i} subjects within the \eqn{i}th units. The inspection time must be positive for all subjects. 
 
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
 Logical. If FALSE, there is no penalty term. The default value is FALSE. If TRUE, users should specify a a positive value for the penalty term lambda.
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

library(gaussquad)
library(extraDistr)
library(fda)
library(splines2)
library(nleqslv)

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
data=data_for_est(r=0,beta=c(-1,1,0.5),gamma=c(-1,1),theta=1,n=300,H=H)
#estimate the parameters
result=GOR_MM(data$Delta,data$X,data$Z,data$n,data$ni,r=data$r,data$C,knotsnum = 2,order=2)

#Scenario 2
#only 2 unit-level covariates (p=0,q=2)
data=data_for_est(r=0,beta=c(0),gamma=c(-1,1),theta=1,n=300,H=H)
#estimate the parameters
result=GOR_MM(data$Delta,X=NULL,data$Z,data$n,data$ni,r=data$r,data$C,knotsnum = 2,order=2)

#Scenario 3
#only 3 subject-level covariates (p=3,q=0)
data=data_for_est(r=0,beta=c(-1,1,0.5),gamma=c(0),theta=1,n=300,H=H)
result=GOR_MM(data$Delta,data$X,Z=NULL,data$n,data$ni,r=data$r,data$C,knotsnum = 2,order=2)

#Scenario 4
#3 subject-level covariates, 2 unit-level covariates (p=3,q=2)
#no clustering effect
data=data_for_est(r=0,beta=c(-1,1,0.5),gamma=c(-1,1),theta=0,n=300,H=H)
result=GOR_MM(data$Delta,data$X,data$Z,data$n,data$ni,r=data$r,data$C,cluster.ind = FALSE,
knotsnum = 2,order=2)

}