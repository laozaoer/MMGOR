\name{SIM_DATA}
\alias{SIM_DATA}
\title{
Generating datasets for the main function GOR_MM
}
\description{
SIM_DATA is used to simulate clustered or non-clustered current status data. This function makes use of the generalized odds rate (GOR) model to generate the time-to-event. 
}
\usage{
SIM_DATA(r, beta, gamma, theta, n, H, par_insp,par_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{r}{
A scalar that specifies the  generalized odds rate model.
}
  \item{beta}{
A \eqn{p\times 1} vector representing the regression parameter for subject-level covariates.
}
  \item{gamma}{
A \eqn{q\times 1 } vector representing the regression parameter for unit-level covariates.
}
  \item{theta}{
A positive scalar representing the standard error of the clustering effect. When it is zero, the subjects are all independent (no clustering). 
}
  \item{n}{
An integer that denotes the number of units (clusters).
}
  \item{H}{
A function representing the \eqn{H} function included the generalized odds rate model. \eqn{H} should be a non-negative and non-decreasing function with \eqn{H(0)=0}.
}
  \item{up.quantile}{
In this function, the inspection time is generated from a uniform distribution on \eqn{(0,b)} with \eqn{b} being the up.quantile\% quantile of the simulated survival times. Users can change this value to adjust the percentage of right censoring.
}
 \item{par_size}{
In this function, the size of each unit is generated from a truncated Poisson distribution with the mean \eqn{\lambda}, and truncated above \eqn{\alpha} and below \eqn{\beta}. Users should input the specific par_size=c(\eqn{\lambda,\alpha,\beta}).
 }
}

\value{
Function SIM_DATA returns a list containing the following components:

  \item{X}{
  A list of \eqn{n} components, the \eqn{i}th component is an \eqn{n_i\times p} matrix of  subject-level covariates. All covariates are generated from Uniform\eqn{(-1,1)}.
  }
  
    \item{Z}{
  An \eqn{n\times q} matrix of unit-level (cluster-level) covariates. All covariates are generated from Uniform\eqn{(-1,1)}.
}

  \item{n}{An integer that denotes the number of units (clusters). }
  \item{ni}{
  An \eqn{n\times 1} vector of integers. Each elements of this vector represents the size of each unit. }
  \item{Delta}{
        A list of \eqn{n} components, the \eqn{i}th component of the list is an \eqn{n_i\times 1} vector containing  indicators for censoring. The values 1 and 0 indicate left censored and right censored, respectively.}
 \item{C}{
 A list of \eqn{n+1} components, the \eqn{i}th component of the first \eqn{n} components is a \eqn{n_i\times 1} vector containing the inspection time of all \eqn{n_i} subjects within the \eqn{i}th unit. The inspection time must be positive for all subjects. The \eqn{(n+1)}th component of \eqn{C} is the lower and upper bounds of inspection times.
 
 }

}
