\name{GOR_MM}
\alias{GOR_MM}
\title{Estimating parameters of the generalized odds rate model with clustered current status data.}
\description{
GOR_MM is used to analyze clustered current status data for the generalized odds rate model.
}
\usage{
GOR_MM(Delta,X,Z,n,ni,r,C,knotsnum,order,cluster.ind=TRUE,
pen.ind=FALSE,lambda = 0,itermax=500,tol=1e-7,quadnum=30)

}
\arguments{
  

  \item{Delta}{
        A list, each component of the list is a vector containing  indicators of the censoring for each observation with 1 or 0 indicating left or right censoring, respectively.}
  \item{X}{
  A list, each component is a matrix containing the within-subject level covariates. If there is no within-subject level covariates, let X=NULL.
  }
  \item{Z}{
  A matrix containing the subject level covariates. If there is no subject level covariate, let Z=NULL.
}
  \item{n}{An integer containing the numbers of total subjects.}
  \item{ni}{
  A vector containing numbers of observations within each subject.}
  \item{r}{
  A scalar which specifies the the generalized odds rate model.}
 \item{C}{
 A list, each component is a vector containing the inspection time for each subject.
 }
    \item{knotsnum}{
    An integer indicating the number of equidistant interior knots for the integrated B-spline approximation of the nonparametric component of the GOR model. 
   }
   \item{order}{
   An integer indicating the  degree of integrated B-splines. }
   \item{cluster.ind=TRUE}{
   Logical. If TRUE, the parameters will be estimated along with clustering effect. The default value is TRUE.
   }
 \item{pen.ind}{
 Logical. If FALSE, there is no penalty term. The default value is FALSE. If TRUE, users should specify a lambda value.
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
set.seed(1)
H=function(t) log(1+t)+t^(3/2) # specify the H function

#This function generates a cohort with 300 subjects, 
#where each subject has up to 8 observations. 
#Both the covariates of subject level and within-subject
#level are generated from uniform(-1,1) distribution. 
#Users can use their own H function by changing the specific H function form.
data=data_for_est(r=0,beta=c(-1),gamma=c(-1),theta=1,n=300,H=H)

#estimate the parameters
result=GOR_MM(data$Delta,data$X,data$Z,data$n,data$ni,r=data$r,data$C,knotsnum = 2,order=2)

}