\name{Sim_data}
\alias{Sim_data}
\title{
Generating datasets for the main function GOR_MM
}
\description{
Sim_data is used to simulate the applicable datasets for the main function GOR_MM. Users can check the outputs of this function to adjust the forms of  their own datasets.
}
\usage{
Sim_data(r, beta, gamma, theta, n, H, C_inte)
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
A \eqn{q\times 1 }vector representing the regression parameter for unit-level covariates.
}
  \item{theta}{
A positive scalar representing the standard error of the clustering effect.
}
  \item{n}{
An integer that denotes the number of units (clusters).
}
  \item{H}{
A function representing the \eqn{H} function included the generalized odds rate model. \eqn{H} should be a non-negative and non-decreasing function with \eqn{H(0)=0}.
}
  \item{C_inte}{
A two-dimensional vector representing the interval of the uniform distribution for generating the inspection times.
}
}

\value{
Function Sim_data returns a list containing the following components:

  \item{X}{
  A list of \eqn{n} components, the \eqn{i}th component is an \eqn{n_i\times p} matrix of  subject-level covariates. All covariates are generated from Uniform\eqn{(-1,1)}.
  }
  
    \item{Z}{
  An \eqn{n\times q} matrix of unit-level (cluster-level) covariates. ll covariates are generated from Uniform\eqn{(-1,1)}.
}

  \item{n}{An integer that denotes the number of units (clusters). }
  \item{ni}{
  An \eqn{n\times 1} vector of integers. Each elements of this vector represents the size of each unit. }
  \item{Delta}{
        A list of \eqn{n} components, the \eqn{i}th component of the list is an \eqn{n_i\times 1} vector containing  indicators for censoring. The values 0 and 1 indicate left censored and right censored, respectively.}
 \item{C}{
 A list of \eqn{n} components, the \eqn{i}th component is a \eqn{n_i\times 1} vector containing the inspection time of all \eqn{n_i} subjects within the \eqn{i}th unit. The inspection time must be positive for all subjects. 
 
 }

}
