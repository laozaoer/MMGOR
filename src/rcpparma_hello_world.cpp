// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#define ARMA_DONT_PRINT_ERRORS
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
arma::mat rcpparma_hello_world() {
    arma::mat m1 = arma::eye<arma::mat>(3, 3);
    arma::mat m2 = arma::eye<arma::mat>(3, 3);
    
    return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector, 
// returning a matrix
//
// [[Rcpp::export]]
arma::mat rcpparma_outerproduct(const arma::colvec & x) {
    arma::mat m = x * x.t();
    return m;
}

// and the inner product returns a scalar
//
// [[Rcpp::export]]
double rcpparma_innerproduct(const arma::colvec & x) {
    double v = arma::as_scalar(x.t() * x);
    return v;
}


// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcpparma_bothproducts(const arma::colvec & x) {
    arma::mat op = x * x.t();
    double    ip = arma::as_scalar(x.t() * x);
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}


// [[Rcpp::export]]
arma::vec converNaN(const arma::vec&x){
    arma::uvec indices = find(x!=x);
    // std::cout << x.is_finite() << std::endl;
    arma::vec result=x;
    result.elem(indices)=arma::zeros<arma::vec>(indices.n_elem);
    // arma::vec r1=x.subvec(result);
    return result;
}




// [[Rcpp::export]]
double likelihoodi(const double&b,const arma::vec&parameters,
                   const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                   const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
    int totaldim=parameters.n_elem;
    double result=0;double S;arma::mat midresult;midresult.zeros(1,1);
    arma::vec covariate(betadim+gammadim);
    if(betadim==0){
        for(int j=0;j<ni;j++){
            
            covariate.subvec(betadim,betadim+gammadim-1)=Z;
            if(r!=0){
                S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                                *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b)),-1/r);
            }
            else{
                S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                               *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b));
            }
            
            if(S>0.99999999999){
                S=0.99999999999;
            }
            
            if((S)<std::pow(10,-30)){
                
                S=std::pow(10,-30);
            }
            result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
        } 
    }else if(gammadim==0){
        for(int j=0;j<ni;j++){
            covariate.subvec(0,betadim-1)=trans(X.row(j));
            if(r!=0){
                S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                                *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b)),-1/r);
            }
            else{
                S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                               *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b));
            }
            
            if(S>0.99999999999){
                S=0.99999999999;
            }
            
            if((S)<std::pow(10,-30)){
                
                S=std::pow(10,-30);
            }
            result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
        }
    }else{
        for(int j=0;j<ni;j++){
            covariate.subvec(0,betadim-1)=trans(X.row(j));
            covariate.subvec(betadim,betadim+gammadim-1)=Z;
            if(r!=0){
                S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                                *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b)),-1/r);
            }
            else{
                S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                               *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b));
            }
            
            if(S>0.99999999999){
                S=0.99999999999;
            }
            
            if((S)<std::pow(10,-30)){
                
                S=std::pow(10,-30);
            }
            result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
        }
    }
    
    
    result=std::exp(result);
    // result=result*std::exp(b*b);
    return result;
}


// [[Rcpp::export]]
double likelihoodinocluster(const arma::vec&parameters,
                            const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                            const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
    int totaldim=parameters.n_elem;
    double result=0;double S;arma::mat midresult;midresult.zeros(1,1);
    arma::vec covariate(betadim+gammadim);
    for(int j=0;j<ni;j++){
        covariate.subvec(0,betadim-1)=trans(X.row(j));
        covariate.subvec(betadim,betadim+gammadim-1)=Z;
        if(r!=0){
            S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim,totaldim-1)))*blC.col(j))
                            *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate))),-1/r);
        }
        else{
            S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim,totaldim-1)))*blC.col(j))
                           *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)));
        }
        
        if(S>0.99999999999){
            S=0.99999999999;
        }
        
        if((S)<std::pow(10,-30)){
            
            S=std::pow(10,-30);
        }
        result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
    }
    
    result=std::exp(result);
    // result=result*std::exp(b*b);
    return result;
}



// [[Rcpp::export]]
arma::mat weightfunction(const arma::vec&parameters,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                         const arma::field<arma::mat>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                         const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim){
    int order=rules.n_rows;
    arma::vec normvec=arma::zeros(order);
    for(int i=0;i<order;i++){
        normvec(i)=R::dnorm(rules(i,0),0,1,0);
    }
    arma::vec weightvec=rules.col(1)%exp(pow(rules.col(0),2))%normvec;
    arma::mat result=arma::zeros(order,n);
    arma::vec functionvalue(order);
    for(int i=0;i<n;i++){
        for(int k=0;k<order;k++){
            functionvalue(k)=likelihoodi(rules(k,0),parameters,Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
            
        }
        result.col(i)=weightvec%functionvalue/sum(weightvec%functionvalue);
    }
    return result;
}



// [[Rcpp::export]]
arma::field<arma::mat> DerivCal(const double&b,const arma::vec&lastpar,
                                const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                                const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
    arma::field<arma::mat> result(4);
    if(betadim==0){
        if(r==0){
            arma::vec lastpsi=lastpar.subvec(betadim+gammadim+1,lastpar.n_elem-1);
            arma::vec lastgamma=lastpar.subvec(betadim,betadim+gammadim-1);
            double lasttheta=std::exp(lastpar(betadim+gammadim));
            arma::vec onesvec=arma::ones(ni);
            arma::vec exppsi=exp(lastpsi);
            arma::vec H1=(trans(blC)*exppsi);
            arma::vec lastpara=H1%exp(onesvec*sum(trans(lastgamma)*Z)+(lasttheta)*b);
            arma::vec A1=exp(-lastpara)%pow((1-exp(-lastpara)),-1);
            arma::vec A2=0.5*(A1+pow(A1,2));
            arma::vec p1,p2,p3;
            p1=Delta%((A1+2*A2%lastpara)%lastpara-onesvec);
            p2=-2*Delta%A2%pow(lastpara,2);
            p3=-(1-Delta)%lastpara;
            arma::vec ParFirstshare,ParSecondshare;
            
            ParFirstshare=p1+p2+p3+Delta;
            ParSecondshare=4*p2+2*p3-2*Delta;
            arma::mat covariate=join_rows((onesvec)*trans(Z),onesvec*b*lasttheta);
            arma::mat logHderiv=trans(blC)%repmat(trans(exppsi),ni,1)%repmat(pow(H1,-1),1,lastpsi.n_elem);
            arma::mat ParFirstDeriv=trans(trans(ParFirstshare)*covariate);
            arma::mat PsiFirstDeriv=trans(trans(ParFirstshare)*logHderiv);
            arma::mat ParSecondDeriv=arma::zeros(betadim+gammadim+1,betadim+gammadim+1);
            arma::mat PsiSecondDeriv=arma::zeros(lastpsi.n_elem,lastpsi.n_elem);
            arma::mat D;
            arma::mat mat1,mat2;
            
            for(int i=0;i<ni;i++){
                ParSecondDeriv=ParSecondDeriv+ParSecondshare(i)*trans(covariate.row(i))*covariate.row(i);
                D=diagmat(blC.col(i));
                mat1=D*diagmat(exppsi)/H1(i);
                mat2=D*exppsi*trans(exppsi)*D/std::pow(H1(i),2);
                PsiSecondDeriv=PsiSecondDeriv+ParFirstshare(i)*(mat1-mat2)+ParSecondshare(i)*mat2;
            }
            ParSecondDeriv(betadim+gammadim,betadim+gammadim)=ParSecondDeriv(betadim+gammadim,betadim+gammadim)+ParFirstDeriv(betadim+gammadim);
            
            result(0)=ParFirstDeriv;
            result(1)=ParSecondDeriv;
            result(2)=PsiFirstDeriv;
            result(3)=PsiSecondDeriv;
        }
        else{
            arma::vec lastpsi=lastpar.subvec(betadim+gammadim+1,lastpar.n_elem-1);
            arma::vec lastgamma=lastpar.subvec(betadim,betadim+gammadim-1);
            double lasttheta=std::exp(lastpar(betadim+gammadim));
            arma::vec onesvec=arma::ones(ni);
            double Const=1/r;
            arma::vec exppsi=exp(lastpsi);
            arma::vec H1=(trans(blC)*exppsi);
            arma::vec A1U,A2U2;
            
            
            
            arma::vec lastpara=H1%exp(onesvec*sum(trans(lastgamma)*Z)+(lasttheta)*b);
            
            arma::vec term1=pow((1+r*lastpara),1/r);
            arma::vec term2=r + 1/lastpara;
            A1U=1/(term1-1)%(1/term2);
            A2U2=1/(term1%term2%term2)%(1+r*(1-1/term1))/pow(1-1/term1,2)/2;
            arma::vec p1,p2,p3;
            p1=Delta%(A1U-onesvec*Const);
            p2=-2*Delta%A2U2;
            p3=-(1-Delta)%pow((pow(lastpara,-1)+r),-1);
            arma::vec ParFirstshare,ParSecondshare;
            
            ParFirstshare=p1+p3+Const*Delta;
            ParSecondshare=4*p2+2*p3-2*Const*Delta;
            arma::mat covariate=join_rows((onesvec)*trans(Z),onesvec*b*lasttheta);
            arma::mat Hminus=repmat(lastpsi,1,lastpsi.n_elem);
            for(int k=0;k<lastpsi.n_elem;k++){
                Hminus.col(k)=Hminus.col(k)-lastpsi(k);
            }
            
            arma::mat blCHminusinverse=pow((trans(blC)*exp(Hminus)),-1);
            arma::mat logHderiv=trans(blC)%blCHminusinverse;
            
            arma::mat ParFirstDeriv=trans(trans(ParFirstshare)*covariate);
            arma::mat PsiFirstDeriv=trans(trans(ParFirstshare)*logHderiv);
            arma::mat ParSecondDeriv=arma::zeros(betadim+gammadim+1,betadim+gammadim+1);
            arma::mat PsiSecondDeriv=arma::zeros(lastpsi.n_elem,lastpsi.n_elem);
            arma::mat D;
            arma::mat mat1,mat2;
            
            for(int i=0;i<ni;i++){
                ParSecondDeriv=ParSecondDeriv+ParSecondshare(i)*trans(covariate.row(i))*covariate.row(i);
                D=diagmat(blC.col(i));
                mat1=D*diagmat(blCHminusinverse.row(i));
                // mat1=D*diagmat(exppsi)/H1(i);
                // mat2=D*exppsi*trans(exppsi)*D/std::pow(H1(i),2);
                mat2=D*trans(blCHminusinverse.row(i))*blCHminusinverse.row(i)*D;
                PsiSecondDeriv=PsiSecondDeriv+ParFirstshare(i)*(mat1-mat2)+ParSecondshare(i)*mat2;
                
            }
            ParSecondDeriv(betadim+gammadim,betadim+gammadim)=ParSecondDeriv(betadim+gammadim,betadim+gammadim)+ParFirstDeriv(betadim+gammadim);
            
            result(0)=ParFirstDeriv;
            result(1)=ParSecondDeriv;
            result(2)=PsiFirstDeriv;
            result(3)=PsiSecondDeriv;
            
            
        }
    }else if(gammadim==0){
        if(r==0){
            arma::vec lastpsi=lastpar.subvec(betadim+gammadim+1,lastpar.n_elem-1);
            arma::vec lastbeta=lastpar.subvec(0,betadim-1);
            double lasttheta=std::exp(lastpar(betadim+gammadim));
            arma::vec onesvec=arma::ones(ni);
            arma::vec exppsi=exp(lastpsi);
            arma::vec H1=(trans(blC)*exppsi);
            arma::vec lastpara=H1%exp(X*lastbeta+(lasttheta)*b);
            arma::vec A1=exp(-lastpara)%pow((1-exp(-lastpara)),-1);
            arma::vec A2=0.5*(A1+pow(A1,2));
            arma::vec p1,p2,p3;
            p1=Delta%((A1+2*A2%lastpara)%lastpara-onesvec);
            p2=-2*Delta%A2%pow(lastpara,2);
            p3=-(1-Delta)%lastpara;
            arma::vec ParFirstshare,ParSecondshare;
            
            ParFirstshare=p1+p2+p3+Delta;
            ParSecondshare=4*p2+2*p3-2*Delta;
            arma::mat covariate=join_rows(X,onesvec*b*lasttheta);
            arma::mat logHderiv=trans(blC)%repmat(trans(exppsi),ni,1)%repmat(pow(H1,-1),1,lastpsi.n_elem);
            arma::mat ParFirstDeriv=trans(trans(ParFirstshare)*covariate);
            arma::mat PsiFirstDeriv=trans(trans(ParFirstshare)*logHderiv);
            arma::mat ParSecondDeriv=arma::zeros(betadim+gammadim+1,betadim+gammadim+1);
            arma::mat PsiSecondDeriv=arma::zeros(lastpsi.n_elem,lastpsi.n_elem);
            arma::mat D;
            arma::mat mat1,mat2;
            
            for(int i=0;i<ni;i++){
                ParSecondDeriv=ParSecondDeriv+ParSecondshare(i)*trans(covariate.row(i))*covariate.row(i);
                D=diagmat(blC.col(i));
                mat1=D*diagmat(exppsi)/H1(i);
                mat2=D*exppsi*trans(exppsi)*D/std::pow(H1(i),2);
                PsiSecondDeriv=PsiSecondDeriv+ParFirstshare(i)*(mat1-mat2)+ParSecondshare(i)*mat2;
            }
            ParSecondDeriv(betadim+gammadim,betadim+gammadim)=ParSecondDeriv(betadim+gammadim,betadim+gammadim)+ParFirstDeriv(betadim+gammadim);
            
            result(0)=ParFirstDeriv;
            result(1)=ParSecondDeriv;
            result(2)=PsiFirstDeriv;
            result(3)=PsiSecondDeriv;
        }
        else{
            arma::vec lastpsi=lastpar.subvec(betadim+gammadim+1,lastpar.n_elem-1);
            arma::vec lastbeta=lastpar.subvec(0,betadim-1);
            double lasttheta=std::exp(lastpar(betadim+gammadim));
            arma::vec onesvec=arma::ones(ni);
            double Const=1/r;
            arma::vec exppsi=exp(lastpsi);
            arma::vec H1=(trans(blC)*exppsi);
            arma::vec A1U,A2U2;
            
            
            
            arma::vec lastpara=H1%exp(X*lastbeta+(lasttheta)*b);
            
            arma::vec term1=pow((1+r*lastpara),1/r);
            arma::vec term2=r + 1/lastpara;
            A1U=1/(term1-1)%(1/term2);
            A2U2=1/(term1%term2%term2)%(1+r*(1-1/term1))/pow(1-1/term1,2)/2;
            arma::vec p1,p2,p3;
            p1=Delta%(A1U-onesvec*Const);
            p2=-2*Delta%A2U2;
            p3=-(1-Delta)%pow((pow(lastpara,-1)+r),-1);
            arma::vec ParFirstshare,ParSecondshare;
            
            ParFirstshare=p1+p3+Const*Delta;
            ParSecondshare=4*p2+2*p3-2*Const*Delta;
            arma::mat covariate=join_rows(X,onesvec*b*lasttheta);
            arma::mat Hminus=repmat(lastpsi,1,lastpsi.n_elem);
            for(int k=0;k<lastpsi.n_elem;k++){
                Hminus.col(k)=Hminus.col(k)-lastpsi(k);
            }
            
            // arma::mat logHderiv=trans(blC)%repmat(trans(exppsi),ni,1)%repmat(pow(H1,-1),1,lastpsi.n_elem);
            arma::mat blCHminusinverse=pow((trans(blC)*exp(Hminus)),-1);
            arma::mat logHderiv=trans(blC)%blCHminusinverse;
            
            arma::mat ParFirstDeriv=trans(trans(ParFirstshare)*covariate);
            arma::mat PsiFirstDeriv=trans(trans(ParFirstshare)*logHderiv);
            arma::mat ParSecondDeriv=arma::zeros(betadim+gammadim+1,betadim+gammadim+1);
            arma::mat PsiSecondDeriv=arma::zeros(lastpsi.n_elem,lastpsi.n_elem);
            arma::mat D;
            arma::mat mat1,mat2;
            
            for(int i=0;i<ni;i++){
                ParSecondDeriv=ParSecondDeriv+ParSecondshare(i)*trans(covariate.row(i))*covariate.row(i);
                D=diagmat(blC.col(i));
                mat1=D*diagmat(blCHminusinverse.row(i));
                // mat1=D*diagmat(exppsi)/H1(i);
                // mat2=D*exppsi*trans(exppsi)*D/std::pow(H1(i),2);
                mat2=D*trans(blCHminusinverse.row(i))*blCHminusinverse.row(i)*D;
                PsiSecondDeriv=PsiSecondDeriv+ParFirstshare(i)*(mat1-mat2)+ParSecondshare(i)*mat2;
                
            }
            ParSecondDeriv(betadim+gammadim,betadim+gammadim)=ParSecondDeriv(betadim+gammadim,betadim+gammadim)+ParFirstDeriv(betadim+gammadim);
            
            result(0)=ParFirstDeriv;
            result(1)=ParSecondDeriv;
            result(2)=PsiFirstDeriv;
            result(3)=PsiSecondDeriv;
            
            
        }
    }else{
        if(r==0){
            arma::vec lastpsi=lastpar.subvec(betadim+gammadim+1,lastpar.n_elem-1);
            arma::vec lastbeta=lastpar.subvec(0,betadim-1);
            arma::vec lastgamma=lastpar.subvec(betadim,betadim+gammadim-1);
            double lasttheta=std::exp(lastpar(betadim+gammadim));
            arma::vec onesvec=arma::ones(ni);
            arma::vec exppsi=exp(lastpsi);
            arma::vec H1=(trans(blC)*exppsi);
            arma::vec lastpara=H1%exp(X*lastbeta+sum(trans(lastgamma)*Z)+(lasttheta)*b);
            arma::vec A1=exp(-lastpara)%pow((1-exp(-lastpara)),-1);
            arma::vec A2=0.5*(A1+pow(A1,2));
            arma::vec p1,p2,p3;
            p1=Delta%((A1+2*A2%lastpara)%lastpara-onesvec);
            p2=-2*Delta%A2%pow(lastpara,2);
            p3=-(1-Delta)%lastpara;
            arma::vec ParFirstshare,ParSecondshare;
            
            ParFirstshare=p1+p2+p3+Delta;
            ParSecondshare=4*p2+2*p3-2*Delta;
            arma::mat covariate=join_rows(join_rows(X,(onesvec)*trans(Z)),onesvec*b*lasttheta);
            arma::mat logHderiv=trans(blC)%repmat(trans(exppsi),ni,1)%repmat(pow(H1,-1),1,lastpsi.n_elem);
            arma::mat ParFirstDeriv=trans(trans(ParFirstshare)*covariate);
            arma::mat PsiFirstDeriv=trans(trans(ParFirstshare)*logHderiv);
            arma::mat ParSecondDeriv=arma::zeros(betadim+gammadim+1,betadim+gammadim+1);
            arma::mat PsiSecondDeriv=arma::zeros(lastpsi.n_elem,lastpsi.n_elem);
            arma::mat D;
            arma::mat mat1,mat2;
            
            for(int i=0;i<ni;i++){
                ParSecondDeriv=ParSecondDeriv+ParSecondshare(i)*trans(covariate.row(i))*covariate.row(i);
                D=diagmat(blC.col(i));
                mat1=D*diagmat(exppsi)/H1(i);
                mat2=D*exppsi*trans(exppsi)*D/std::pow(H1(i),2);
                PsiSecondDeriv=PsiSecondDeriv+ParFirstshare(i)*(mat1-mat2)+ParSecondshare(i)*mat2;
            }
            ParSecondDeriv(betadim+gammadim,betadim+gammadim)=ParSecondDeriv(betadim+gammadim,betadim+gammadim)+ParFirstDeriv(betadim+gammadim);
            
            result(0)=ParFirstDeriv;
            result(1)=ParSecondDeriv;
            result(2)=PsiFirstDeriv;
            result(3)=PsiSecondDeriv;
        }
        else{
            arma::vec lastpsi=lastpar.subvec(betadim+gammadim+1,lastpar.n_elem-1);
            arma::vec lastbeta=lastpar.subvec(0,betadim-1);
            arma::vec lastgamma=lastpar.subvec(betadim,betadim+gammadim-1);
            double lasttheta=std::exp(lastpar(betadim+gammadim));
            arma::vec onesvec=arma::ones(ni);
            double Const=1/r;
            arma::vec exppsi=exp(lastpsi);
            arma::vec H1=(trans(blC)*exppsi);
            arma::vec A1U,A2U2;
            
            
            
            arma::vec lastpara=H1%exp(X*lastbeta+sum(trans(lastgamma)*Z)+(lasttheta)*b);
            
            arma::vec term1=pow((1+r*lastpara),1/r);
            arma::vec term2=r + 1/lastpara;
            A1U=1/(term1-1)%(1/term2);
            A2U2=1/(term1%term2%term2)%(1+r*(1-1/term1))/pow(1-1/term1,2)/2;
            arma::vec p1,p2,p3;
            p1=Delta%(A1U-onesvec*Const);
            p2=-2*Delta%A2U2;
            p3=-(1-Delta)%pow((pow(lastpara,-1)+r),-1);
            arma::vec ParFirstshare,ParSecondshare;
            
            ParFirstshare=p1+p3+Const*Delta;
            ParSecondshare=4*p2+2*p3-2*Const*Delta;
            arma::mat covariate=join_rows(join_rows(X,(onesvec)*trans(Z)),onesvec*b*lasttheta);
            arma::mat Hminus=repmat(lastpsi,1,lastpsi.n_elem);
            for(int k=0;k<lastpsi.n_elem;k++){
                Hminus.col(k)=Hminus.col(k)-lastpsi(k);
            }
            
            // arma::mat logHderiv=trans(blC)%repmat(trans(exppsi),ni,1)%repmat(pow(H1,-1),1,lastpsi.n_elem);
            arma::mat blCHminusinverse=pow((trans(blC)*exp(Hminus)),-1);
            arma::mat logHderiv=trans(blC)%blCHminusinverse;
            
            arma::mat ParFirstDeriv=trans(trans(ParFirstshare)*covariate);
            arma::mat PsiFirstDeriv=trans(trans(ParFirstshare)*logHderiv);
            arma::mat ParSecondDeriv=arma::zeros(betadim+gammadim+1,betadim+gammadim+1);
            arma::mat PsiSecondDeriv=arma::zeros(lastpsi.n_elem,lastpsi.n_elem);
            arma::mat D;
            arma::mat mat1,mat2;
            
            for(int i=0;i<ni;i++){
                ParSecondDeriv=ParSecondDeriv+ParSecondshare(i)*trans(covariate.row(i))*covariate.row(i);
                D=diagmat(blC.col(i));
                mat1=D*diagmat(blCHminusinverse.row(i));
                // mat1=D*diagmat(exppsi)/H1(i);
                // mat2=D*exppsi*trans(exppsi)*D/std::pow(H1(i),2);
                mat2=D*trans(blCHminusinverse.row(i))*blCHminusinverse.row(i)*D;
                PsiSecondDeriv=PsiSecondDeriv+ParFirstshare(i)*(mat1-mat2)+ParSecondshare(i)*mat2;
                
            }
            ParSecondDeriv(betadim+gammadim,betadim+gammadim)=ParSecondDeriv(betadim+gammadim,betadim+gammadim)+ParFirstDeriv(betadim+gammadim);
            
            result(0)=ParFirstDeriv;
            result(1)=ParSecondDeriv;
            result(2)=PsiFirstDeriv;
            result(3)=PsiSecondDeriv;
            
            
        }
    }
    
    
    
    
    return result; 
}



// [[Rcpp::export]]
arma::field<arma::mat> DerivCalnocluster(const arma::vec&lastpar,
                                         const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                                         const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
    arma::field<arma::mat> result(4);
    if(r==0){
        arma::vec lastpsi=lastpar.subvec(betadim+gammadim,lastpar.n_elem-1);
        arma::vec lastbeta=lastpar.subvec(0,betadim-1);
        arma::vec lastgamma=lastpar.subvec(betadim,betadim+gammadim-1);
        arma::vec onesvec=arma::ones(ni);
        arma::vec exppsi=exp(lastpsi);
        arma::vec H1=(trans(blC)*exppsi);
        arma::vec lastpara=H1%exp(X*lastbeta+sum(trans(lastgamma)*Z));
        arma::vec A1=exp(-lastpara)%pow((1-exp(-lastpara)),-1);
        arma::vec A2=0.5*(A1+pow(A1,2));
        arma::vec p1,p2,p3;
        p1=Delta%((A1+2*A2%lastpara)%lastpara-onesvec);
        p2=-2*Delta%A2%pow(lastpara,2);
        p3=-(1-Delta)%lastpara;
        arma::vec ParFirstshare,ParSecondshare;
        
        ParFirstshare=p1+p2+p3+Delta;
        ParSecondshare=4*p2+2*p3-2*Delta;
        arma::mat covariate=join_rows(X,(onesvec)*trans(Z));
        arma::mat logHderiv=trans(blC)%repmat(trans(exppsi),ni,1)%repmat(pow(H1,-1),1,lastpsi.n_elem);
        arma::mat ParFirstDeriv=trans(trans(ParFirstshare)*covariate);
        arma::mat PsiFirstDeriv=trans(trans(ParFirstshare)*logHderiv);
        arma::mat ParSecondDeriv=arma::zeros(betadim+gammadim,betadim+gammadim);
        arma::mat PsiSecondDeriv=arma::zeros(lastpsi.n_elem,lastpsi.n_elem);
        arma::mat D;
        arma::mat mat1,mat2;
        
        for(int i=0;i<ni;i++){
            ParSecondDeriv=ParSecondDeriv+ParSecondshare(i)*trans(covariate.row(i))*covariate.row(i);
            D=diagmat(blC.col(i));
            mat1=D*diagmat(exppsi)/H1(i);
            mat2=D*exppsi*trans(exppsi)*D/std::pow(H1(i),2);
            PsiSecondDeriv=PsiSecondDeriv+ParFirstshare(i)*(mat1-mat2)+ParSecondshare(i)*mat2;
        }
        
        result(0)=ParFirstDeriv;
        result(1)=ParSecondDeriv;
        result(2)=PsiFirstDeriv;
        result(3)=PsiSecondDeriv;
    }
    else{
        arma::vec lastpsi=lastpar.subvec(betadim+gammadim,lastpar.n_elem-1);
        arma::vec lastbeta=lastpar.subvec(0,betadim-1);
        arma::vec lastgamma=lastpar.subvec(betadim,betadim+gammadim-1);
        arma::vec onesvec=arma::ones(ni);
        double Const=1/r;
        arma::vec exppsi=exp(lastpsi);
        arma::vec H1=(trans(blC)*exppsi);
        arma::vec A1U,A2U2;
        
        
        
        arma::vec lastpara=H1%exp(X*lastbeta+sum(trans(lastgamma)*Z));
        
        arma::vec term1=pow((1+r*lastpara),1/r);
        arma::vec term2=r + 1/lastpara;
        A1U=1/(term1-1)%(1/term2);
        A2U2=1/(term1%term2%term2)%(1+r*(1-1/term1))/pow(1-1/term1,2)/2;
        arma::vec p1,p2,p3;
        p1=Delta%(A1U-onesvec*Const);
        p2=-2*Delta%A2U2;
        p3=-(1-Delta)%pow((pow(lastpara,-1)+r),-1);
        arma::vec ParFirstshare,ParSecondshare;
        
        ParFirstshare=p1+p3+Const*Delta;
        ParSecondshare=4*p2+2*p3-2*Const*Delta;
        arma::mat covariate=join_rows(X,(onesvec)*trans(Z));
        arma::mat Hminus=repmat(lastpsi,1,lastpsi.n_elem);
        for(int k=0;k<lastpsi.n_elem;k++){
            Hminus.col(k)=Hminus.col(k)-lastpsi(k);
        }
        
        // arma::mat logHderiv=trans(blC)%repmat(trans(exppsi),ni,1)%repmat(pow(H1,-1),1,lastpsi.n_elem);
        arma::mat blCHminusinverse=pow((trans(blC)*exp(Hminus)),-1);
        arma::mat logHderiv=trans(blC)%blCHminusinverse;
        
        arma::mat ParFirstDeriv=trans(trans(ParFirstshare)*covariate);
        arma::mat PsiFirstDeriv=trans(trans(ParFirstshare)*logHderiv);
        arma::mat ParSecondDeriv=arma::zeros(betadim+gammadim,betadim+gammadim);
        arma::mat PsiSecondDeriv=arma::zeros(lastpsi.n_elem,lastpsi.n_elem);
        arma::mat D;
        arma::mat mat1,mat2;
        
        for(int i=0;i<ni;i++){
            ParSecondDeriv=ParSecondDeriv+ParSecondshare(i)*trans(covariate.row(i))*covariate.row(i);
            D=diagmat(blC.col(i));
            mat1=D*diagmat(blCHminusinverse.row(i));
            // mat1=D*diagmat(exppsi)/H1(i);
            // mat2=D*exppsi*trans(exppsi)*D/std::pow(H1(i),2);
            mat2=D*trans(blCHminusinverse.row(i))*blCHminusinverse.row(i)*D;
            PsiSecondDeriv=PsiSecondDeriv+ParFirstshare(i)*(mat1-mat2)+ParSecondshare(i)*mat2;
            
        }
        
        result(0)=ParFirstDeriv;
        result(1)=ParSecondDeriv;
        result(2)=PsiFirstDeriv;
        result(3)=PsiSecondDeriv;
        
        
    }
    return result; 
}



// [[Rcpp::export]]
arma::vec UpdateOnce(const arma::vec&lastpar,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                     const arma::field<arma::mat>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                     const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                     const bool&penind,const double&lambda,const arma::mat&R){
    int order=rules.n_rows;
    arma::mat weightmat=weightfunction(lastpar,rules,Delta,X,Z,n,ni,r,blC,betadim,gammadim);
    arma::mat ParFirstDeriv,PsiFirstDeriv,ParSecond,PsiSecond;
    ParFirstDeriv=arma::zeros(betadim+gammadim+1,1);
    PsiFirstDeriv=arma::zeros(lastpar.n_elem-betadim-gammadim-1,1);
    ParSecond=arma::zeros(betadim+gammadim+1,betadim+gammadim+1);
    PsiSecond=arma::zeros(lastpar.n_elem-betadim-gammadim-1,lastpar.n_elem-betadim-gammadim-1);
    arma::mat onesvec=arma::ones(n);
    arma::field<arma::mat> Derivresult;
    arma::mat penfirst;
    arma::mat pensecond;
    arma::mat exppsi;
    for(int i=0;i<n;i++){
        if(ni(i)==1){
            Derivresult=DerivCal(0,lastpar,Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
            ParFirstDeriv=ParFirstDeriv+Derivresult(0);
            PsiFirstDeriv=PsiFirstDeriv+Derivresult(2);
            ParSecond=ParSecond+Derivresult(1);
            PsiSecond=PsiSecond+Derivresult(3);
        }else{
            for(int k=0;k<order;k++){
                Derivresult=DerivCal(rules(k,0),lastpar,Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
                ParFirstDeriv=ParFirstDeriv+weightmat(k,i)*Derivresult(0);
                PsiFirstDeriv=PsiFirstDeriv+weightmat(k,i)*Derivresult(2);
                ParSecond=ParSecond+weightmat(k,i)*Derivresult(1);
                PsiSecond=PsiSecond+weightmat(k,i)*Derivresult(3);
                // std::cout<<i<<k<<ParFirstDeriv;
            }
        }

    }
    if(penind){
        exppsi=exp(lastpar.subvec(betadim+gammadim+1,lastpar.n_elem-1));
        penfirst=2*diagmat(exppsi)*(R*exppsi);
        arma::mat extramat;
        pensecond=2*(R%(exppsi*trans(exppsi))+diagmat(R*exppsi*trans(exppsi)));
        PsiSecond=PsiSecond-lambda*pensecond;
        PsiFirstDeriv=PsiFirstDeriv-lambda*penfirst;
    }
    arma::mat Parresult=lastpar.subvec(0,betadim+gammadim)-solve(ParSecond,ParFirstDeriv);
    
    arma::mat Psiresult=lastpar.subvec(betadim+gammadim+1,lastpar.n_elem-1)-solve(PsiSecond,PsiFirstDeriv);
    arma::mat result=join_cols(Parresult,Psiresult);
    
    return result;
}

// [[Rcpp::export]]
arma::vec UpdateOncenocluster(const arma::vec&lastpar,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                              const arma::field<arma::mat>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                              const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                              const bool&penind,const double&lambda,const arma::mat&R){
    arma::mat ParFirstDeriv,PsiFirstDeriv,ParSecond,PsiSecond;
    ParFirstDeriv=arma::zeros(betadim+gammadim,1);
    PsiFirstDeriv=arma::zeros(lastpar.n_elem-betadim-gammadim,1);
    ParSecond=arma::zeros(betadim+gammadim,betadim+gammadim);
    PsiSecond=arma::zeros(lastpar.n_elem-betadim-gammadim,lastpar.n_elem-betadim-gammadim);
    arma::mat onesvec=arma::ones(n);
    arma::field<arma::mat> Derivresult;
    arma::mat penfirst;
    arma::mat pensecond;
    arma::mat exppsi;
    
    for(int i=0;i<n;i++){
        Derivresult=DerivCalnocluster(lastpar,Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
        ParFirstDeriv=ParFirstDeriv+Derivresult(0);
        PsiFirstDeriv=PsiFirstDeriv+Derivresult(2);
        ParSecond=ParSecond+Derivresult(1);
        PsiSecond=PsiSecond+Derivresult(3);
        
        
    }
    if(penind){
        exppsi=exp(lastpar.subvec(betadim+gammadim,lastpar.n_elem-1));
        penfirst=2*diagmat(exppsi)*(R*exppsi);
        arma::mat extramat;
        pensecond=2*(R%(exppsi*trans(exppsi))+diagmat(R*exppsi*trans(exppsi)));
        PsiSecond=PsiSecond-lambda*pensecond;
        PsiFirstDeriv=PsiFirstDeriv-lambda*penfirst;
    }
    arma::mat Parresult=lastpar.subvec(0,betadim+gammadim-1)-solve(ParSecond,ParFirstDeriv);
    
    arma::mat Psiresult=lastpar.subvec(betadim+gammadim,lastpar.n_elem-1)-solve(PsiSecond,PsiFirstDeriv);
    arma::mat result=join_cols(Parresult,Psiresult);
    
    return result;
}

// [[Rcpp::export]]
arma::vec MainFunc(const arma::vec&lastpar,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                   const arma::field<arma::mat>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                   const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,const int&itermax,const double&criterion,
                   const bool&penind,const double&lambda,const arma::mat&R){
    arma::mat result,lastresult;
    lastresult=result=lastpar;
    double difference;
    int itertime=1;
    do{
        result=UpdateOnce(lastresult,rules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,penind,lambda,R);
        difference=accu(abs(result-lastresult)%abs(pow(lastresult,-1)));
        lastresult=result;
        itertime++;
    } while (difference>criterion&&itertime<itermax);
    
    
    return result;
}


// [[Rcpp::export]]
arma::vec MainFuncnocluster(const arma::vec&lastpar,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                            const arma::field<arma::mat>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                            const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,const int&itermax,const double&criterion,
                            const bool&penind,const double&lambda,const arma::mat&R){
    arma::mat result,lastresult;
    lastresult=result=lastpar;
    double difference;
    int itertime=1;
    do{
        result=UpdateOncenocluster(lastresult,rules,Delta,X,Z,n,ni,r,blC,betadim,gammadim,penind,lambda,R);
        difference=accu(abs(result-lastresult)%abs(pow(lastresult,-1)));
        lastresult=result;
        itertime++;
    } while (difference>criterion&&itertime<itermax);
    
    
    return result;
}


// [[Rcpp::export]]
double likelihoodfunc1current(const double&b,const arma::vec&parameters,
                              const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                              const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
    int totaldim=parameters.n_elem;
    double result=0;double S;arma::mat midresult;midresult.zeros(1,1);
    double Sderiv;
    arma::vec covariate(betadim+gammadim);
    double uncurerate;
    if(betadim==0){
        for(int j=0;j<ni;j++){
            covariate.subvec(betadim,betadim+gammadim-1)=Z;
            
            // uncurerate=1/(1+std::exp(-parameters(0)-sum(parameters.subvec(1,zetadim-1)%covariate)+std::exp(parameters(zetadim))*b));
            
            
            if(r==0){
                S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                               *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b));
            }
            else{
                S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                                *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b)),-1/r);
            }
            
            
            if(S>0.99999999999){
                S=0.99999999999;
            }
            if(Sderiv<std::pow(10,-30)){
                Sderiv=std::pow(10,-30);
            }
            if(uncurerate<std::pow(10,-30)){
                uncurerate=std::pow(10,-30);
            }
            if((S)<std::pow(10,-30)){
                
                S=std::pow(10,-30);
            }
            result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
        }
    }else if(gammadim==0){
        for(int j=0;j<ni;j++){
            covariate.subvec(0,betadim-1)=trans(X.row(j));
            
            // uncurerate=1/(1+std::exp(-parameters(0)-sum(parameters.subvec(1,zetadim-1)%covariate)+std::exp(parameters(zetadim))*b));
            
            
            if(r==0){
                S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                               *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b));
            }
            else{
                S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                                *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b)),-1/r);
            }
            
            
            if(S>0.99999999999){
                S=0.99999999999;
            }
            if(Sderiv<std::pow(10,-30)){
                Sderiv=std::pow(10,-30);
            }
            if(uncurerate<std::pow(10,-30)){
                uncurerate=std::pow(10,-30);
            }
            if((S)<std::pow(10,-30)){
                
                S=std::pow(10,-30);
            }
            result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
        }
    }else{
        for(int j=0;j<ni;j++){
            covariate.subvec(0,betadim-1)=trans(X.row(j));
            covariate.subvec(betadim,betadim+gammadim-1)=Z;
            
            // uncurerate=1/(1+std::exp(-parameters(0)-sum(parameters.subvec(1,zetadim-1)%covariate)+std::exp(parameters(zetadim))*b));
            
            
            if(r==0){
                S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                               *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b));
            }
            else{
                S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim+1,totaldim-1)))*blC.col(j))
                                *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)+std::exp(parameters(betadim+gammadim))*b)),-1/r);
            }
            
            
            if(S>0.99999999999){
                S=0.99999999999;
            }
            if(Sderiv<std::pow(10,-30)){
                Sderiv=std::pow(10,-30);
            }
            if(uncurerate<std::pow(10,-30)){
                uncurerate=std::pow(10,-30);
            }
            if((S)<std::pow(10,-30)){
                
                S=std::pow(10,-30);
            }
            result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
        }  
    }
    
    
    
    result=result+R::dnorm(b,0,1,true);
    result=std::exp(result);
    result=result*std::exp(b*b);
    return result;
}

// This function is used to caculate the likelihood value using Hermit quadrature.


// [[Rcpp::export]]
double likelihoodfunc1currentnocluster(const double&b,const arma::vec&parameters,
                                       const arma::vec&Delta,const arma::mat&X,const arma::vec&Z,
                                       const int&ni,const double&r,const arma::mat&blC,const int&betadim,const int&gammadim){
    int totaldim=parameters.n_elem;
    double result=0;double S;arma::mat midresult;midresult.zeros(1,1);
    double Sderiv;
    arma::vec covariate(betadim+gammadim);
    double uncurerate;
    
    for(int j=0;j<ni;j++){
        covariate.subvec(0,betadim-1)=trans(X.row(j));
        covariate.subvec(betadim,betadim+gammadim-1)=Z;
        
        // uncurerate=1/(1+std::exp(-parameters(0)-sum(parameters.subvec(1,zetadim-1)%covariate)+std::exp(parameters(zetadim))*b));
        
        
        if(r==0){
            S=std::exp(-sum(trans(exp(parameters.subvec(betadim+gammadim,totaldim-1)))*blC.col(j))
                           *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate)));
        }
        else{
            S=std::pow((1+r*sum(trans(exp(parameters.subvec(betadim+gammadim,totaldim-1)))*blC.col(j))
                            *std::exp(sum(parameters.subvec(0,betadim+gammadim-1)%covariate))),-1/r);
        }
        
        
        if(S>0.99999999999){
            S=0.99999999999;
        }
        if(Sderiv<std::pow(10,-30)){
            Sderiv=std::pow(10,-30);
        }
        if(uncurerate<std::pow(10,-30)){
            uncurerate=std::pow(10,-30);
        }
        if((S)<std::pow(10,-30)){
            
            S=std::pow(10,-30);
        }
        result=result+Delta(j)*(log(1-S))+(1-Delta(j))*std::log(S);
    }  
    
    
    
    
    result=std::exp(result);
    return result;
}





// [[Rcpp::export]]
double testquadrature1current(const arma::vec&parameters,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                              const arma::field<arma::mat>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                              const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                              const bool&penind,const double&lambda,const arma::mat&R){
    int totaldim=parameters.n_elem;
    int order=rules.n_rows;double result=0;
    arma::vec weightvec;weightvec=rules.col(1);double term1;
    arma::vec functionvalue(order);
    arma::mat penalty;
    arma::mat exppsi;
    arma::vec estpsi=parameters.subvec(betadim+gammadim+1,totaldim-1);
    for(int i=0;i<n;i++){
        
        if(ni(i)==1){
            term1=likelihoodi(0,parameters,
                                         Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
            if(term1<std::pow(10,-30)){
                term1=std::pow(10,-30);}
            result=result+std::log(term1);
        }else{
            for(int k=0;k<order;k++){
                functionvalue(k)=likelihoodfunc1current(rules(k,0),parameters,
                              Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
            }
            
            term1=sum(functionvalue%weightvec);
            if(term1<std::pow(10,-30)){
                term1=std::pow(10,-30);}
            result=result+std::log(term1);
        }
        
        
    }
    if(penind){
        exppsi=exp(parameters.subvec(betadim+gammadim+1,totaldim-1));
        penalty=trans(exppsi)*R*exppsi;
        result=result-lambda*penalty(0,0);
        
    }
    return -result;
}


// [[Rcpp::export]]
double testquadrature1currentnocluster(const arma::vec&parameters,const arma::mat&rules,const arma::field<arma::vec>&Delta,
                                       const arma::field<arma::mat>&X,const arma::mat&Z,const int&n,const arma::vec&ni,
                                       const double&r,const arma::field<arma::mat>&blC,const int&betadim,const int&gammadim,
                                       const bool&penind,const double&lambda,const arma::mat&R){
    int totaldim=parameters.n_elem;
    int order=rules.n_rows;double result=0;
    arma::vec weightvec;weightvec=rules.col(1);double term1;
    arma::vec functionvalue(order);
    arma::mat penalty;
    arma::mat exppsi;
    arma::vec estpsi=parameters.subvec(betadim+gammadim,totaldim-1);
    for(int i=0;i<n;i++){
        
        for(int k=0;k<order;k++){
            functionvalue(k)=likelihoodfunc1currentnocluster(rules(k,0),parameters,
                          Delta(i),X(i),trans(Z.row(i)),ni(i),r,blC(i),betadim,gammadim);
        }
        
        term1=sum(functionvalue%weightvec);
        if(term1<std::pow(10,-30)){
            term1=std::pow(10,-30);}
        result=result+std::log(term1);
        
    }
    if(penind){
        exppsi=exp(parameters.subvec(betadim+gammadim,totaldim-1));
        penalty=trans(exppsi)*R*exppsi;
        result=result-lambda*penalty(0,0);
        
    }
    return -result;
}

