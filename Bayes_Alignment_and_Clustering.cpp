// [[Rcpp::depends(RcppGSL)]]
 //#include <RcppGSL.h>
# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_vector.h>
#include <math.h>        
#include <RcppArmadilloExtensions/sample.h>
#include <vector>

/*------------------------------------------------------------------------------------------------------
This is for Bayesian clustering and registration of functional observation.
The function BayesRegCluster runs only in the burn-in peroid and is for choosing C and intial values MCMC.
The function MCMC is the regular MCMC whose initial values are picked by BayesRegCluster.
Last modified 01/09/2015
--------------------------------------------------------------------------------------------------------*/
using namespace Rcpp;

//Truncated normal N(y;mu,s)I(a<y<b); where b can be R_PosInf
// Rejection algorithm with a truncated expoential proposal for N(0,1)I(a<x<b) when a is very large: |a| < |b|
double rtexp(double a, double b){
  int stop = false;
  double twoasp = 2*std::pow(a,2);
  double expab = std::exp(-a*(b-a)) - 1;
  double z, e;
  while(!stop){
    R_CheckUserInterrupt();
    z = std::log(1 + unif_rand()*expab);
    e = -std::log(unif_rand());
    stop = (twoasp*e > std::pow(z,2));
  }
  return (a - z/a);
}

 
double trun_rnorm(const double mu, const double s, double a, double b){
  double xmin = -2.00443204036;                 // Left bound
  double xmax =  3.48672170399;                 // Right bound
  int stop = false;
  double r;
  //if( mu+ELARGE<0 ) return(a);
  //scalling
  if(mu!=0 || s!=1){
    a=(a-mu)/s;
    b=(b-mu)/s;
  }
  // Check if a < b
  if(a>=b){
    Rprintf( "*** B must be greater than A ! ***" ); return(NA_REAL);
  }
  else if(std::abs(a)>std::abs(b)) r = -trun_rnorm(0, 1, -b, -a);
  // If a in the right tail (a > xmax), use rejection algorithm with a truncated exponential proposal  
  else if(a>xmax) r = rtexp(a,b);
  // If a in the left tail (a < xmin), use rejection algorithm with a Gaussian proposal
  else if(a<xmin){
    while(!stop){
      R_CheckUserInterrupt();
      r = norm_rand();
      stop = (r>=a) && (r<=b);
    }
  }  
  // In other cases (xmin < a < xmax)
  else{
    double CDFa = Rf_pnorm5(a, 0, 1.0, true, false);
    double CDFb = Rf_pnorm5(b, 0, 1.0, true, false);
    double u = unif_rand();
    double CDFxi = CDFa + u*(CDFb - CDFa);
    r = Rf_qnorm5(CDFxi, 0, 1, true, false);
  }
  // Scaling
  if(mu!=0 || s!=1)
  r = r*s + mu;
  return r;
}

//[[Rcpp::export]]
double linearapprox(double x, arma::mat vec){
   if (x!=1){
        double N = vec.n_cols;
        double k = floor(x*(N-1));
        double  x0 = k/(N-1); double y0 = vec(0,k);
        double x1 = (k+1)/(N-1); double y1 = vec(0,k+1);
        double y = y0 + (y1-y0)*(x-x0)/(x1-x0);
        return(y);
   }else{return(1);}
}

//evaluate the B-spline basis function on T=[0,1]
//mix of GSL and Rcpp vector??
// [[Rcpp::export]]
 arma::mat evalbasis(double t, const size_t degree, const size_t nbasis){
  const size_t nbreak = nbasis + 2 - degree;
  size_t j;
  gsl_bspline_workspace *bw; //initialize the B-spline working space
  bw = gsl_bspline_alloc(degree, nbreak); //specify the B spline basis
  gsl_vector *B;
  B = gsl_vector_alloc(nbasis); //store basis function value at t in this vector
  gsl_bspline_knots_uniform(0.0,1.0,bw); // use uniform breakpoints on [0,1]
  gsl_bspline_eval (t, B, bw);
  arma::mat temp(1,nbasis);
  for ( j=0;j<nbasis;++j){
    double Bj = gsl_vector_get(B,j);
    temp(0,j) = Bj;
  }
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  return(temp);
}





// [[Rcpp::export]]
IntegerVector oneMultinomC(NumericVector probs) {
    int k = probs.size();
    IntegerVector ans(k);
    rmultinom(1, probs.begin(), k, ans.begin());
    return(ans);
}


//  The following function identify the 1 postion in a vector 
// [[Rcpp::export]]
int position1(IntegerVector vec){
  int temp = -1;
  int pos=0;
  while(temp<0){
   pos= pos+1;
   if (vec(pos-1)==1) temp=1;
  }
  return(pos); 
}
    
// [[Rcpp::export]]
double psum (arma::mat temp){
    double sum = 0;
    int n = temp.n_cols-1;
    for (int i=0;i<n;i++){
      sum = sum + temp(0,i);
    }
    return sum;
}

 


// [[Rcpp::export]]
arma::mat add0(arma::mat myvector){
    int k =myvector.n_cols;
    arma::mat temp(1,k+1);
    temp(0,0) = 0;
   for (int i=1; i<k+1;i++){
     temp(0,i) = myvector(0,i-1);
   }
   return(temp);
}

// [[Rcpp::export]]
double normsq(arma::mat vec){
  int k = vec.n_cols;
  double temp=0;
  for (int i=0;i<k;i++){
    temp=temp + vec(0,i)*vec(0,i);
  }
  return(temp);
}

 

// [[Rcpp::export]]
arma::mat cumsumC(arma::mat x){
    // initialize an accumulator variable
    double acc = 0;
        
    // initialize the result vector
    int l=x.n_cols;
    arma::mat res(1,l);
        
    for(int i = 0; i < l; i++){
         acc += x(0,i);
         res(0,i) = acc;
    }
    return res;
}
// [[Rcpp::export]]
arma::mat rdiri(NumericVector a){
   int n=a.size();
   arma::mat result(1,n);
   NumericVector vec(n);
   double tempsum=0;
   for (int i=0;i<n;i++){
     RNGScope scope; // ensure RNG gets set/reset
     double temp = R::rgamma(a(i),1);
     vec(i)=temp;
    tempsum=tempsum+temp;
   }
   for (int j=0;j<n;j++){
     
     result(0,j)=vec(j)/tempsum;
  
   }
   
 return(result);  
  
}
// [[Rcpp::export]]
NumericVector probparam(arma::mat member, double eta,int C){
   int n=member.n_cols;
   NumericVector result(C);
   for (int j=0;j<C;j++){
     int count=0;
     for (int i=0;i<n;i++){
         if(member(0,i)==j+1){count=count+1;}
      }
    result(j) = count+eta;
   }
   return(result);
}

 // [[Rcpp::export]]
IntegerVector csample ( IntegerVector x,
                           int size,
                           bool replace,
                           NumericVector prob = NumericVector::create()
                           ) {
  RNGScope scope;
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}


//The following function deletes the ith element of a NumericVector
// [[Rcpp::export]]
IntegerVector delete1(IntegerVector vec,int i){
     int n=vec.size();
     IntegerVector temp(n-1);
      for (int k=0;k<i-1;k++){
         temp(k)=vec(k);        
      }
     for (int j=i;j<n;j++){
       temp(j-1)=vec(j);
     }
     return(temp);
}

//This is equivalent to sort function in base R
// [[Rcpp::export]]
IntegerVector stl_sort(IntegerVector x) {
   IntegerVector y = clone(x);
   std::sort(y.begin(), y.end());
   return y;
}

NumericVector stl_sortnum(NumericVector x) {
   NumericVector y = clone(x);
   std::sort(y.begin(), y.end());
   return y;
}


// [[Rcpp::export]]
IntegerVector add1element(IntegerVector vec,int a){
     int n=vec.size();
     IntegerVector temp(n+1);
     for (int i=0;i<n;i++ ){temp(i)=vec(i);}
     temp(n)=a;
     return(temp);
}

// [[Rcpp::export]]
IntegerVector whichC(IntegerVector vec,int a){
    int n = vec.size();
    std::vector<int> temp;
    for (int i=0;i<n;i++){
      
      if(vec(i)==a){
        
         temp.push_back(i+1);
        }
    }
    return(wrap(temp));    
}

IntegerVector arma2vec(arma::mat vec){
  int n = vec.n_cols;
  IntegerVector temp(n);
  for (int i=0;i<n;i++){
    temp(i) = vec(0,i);
  }
  return(temp);
}


NumericVector arma2vecnum(arma::mat vec){
  int n = vec.n_cols;
  NumericVector temp(n);
  for (int i=0;i<n;i++){
    temp(i) = vec(0,i);
  }
  return(temp);
}



// [[Rcpp::export]]
int whichC2(IntegerVector vec,int a){
    int n = vec.size();
    int count = 0;
    
    for (int i=0;i<n;i++){
      count=count+1;
      if(vec(i)==a){ break;}
    }
    return(count);    
}


//the following function delete entire vec small from vec big
// [[Rcpp::export]]
IntegerVector deletevec(IntegerVector big,IntegerVector small){
    int n = small.size();
    IntegerVector temp=clone(big);
    //int count = 0;
    for (int i=0;i<n;i++){
      //count = count+1;
      int index = whichC2(temp,small(i));
      temp = delete1(temp,index);
    }
    return(temp);   
}

//the following function combine 2 vectors
// [[Rcpp::export]]
IntegerVector combinevec(IntegerVector big,IntegerVector small){
    int n = small.size();
   //int N = big.size();
    IntegerVector temp=clone(big);
    for (int i=0;i<n;i++){
       temp = add1element(temp,small(i));
    }
    return(temp);   
    
}
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::mat mu, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}
// [[Rcpp::export]]
double truncnormal(double a,double b,double mean,double var){
    double temp = b+1;
    while(temp>b||temp<a){ 
      RNGScope scope; // ensure RNG gets set/reset
      temp = R::rnorm(mean,sqrt(var));
      }
    return temp;
}

// [[Rcpp::export]]
NumericVector vecsubsetter(NumericVector vec,IntegerVector index){
   int n = index.size();
   NumericVector temp(n);
   for (int i=0;i<n;i++){
     temp[i] = vec[index(i)-1];
   }
   return(temp);
}



// [[Rcpp::export]]
IntegerVector tempfunc(int totdraws,int wait){
     std::vector<int> temp;
    for (int v=0;v< floor(totdraws/wait);v++){
      
     for (int a=5;a<wait;a++ ){
         temp.push_back(v*10+a);
        }
    }
    return(wrap(temp));    
}


// [[Rcpp::export]]
List BayesRegclust(const arma::mat dat, int totdraw,int C,int M,double kappa,double theta, double alpha, double phi, double warpingsd, 
                    double sigmaasq, double eta, double varbeta,int nbasis, int degree, NumericVector clocktime,arma::mat SigmaB,int wait,double percent){//func
  //J is number of obs
  const int J = dat.n_rows;
  //K is number of measurements
  const int K = dat.n_cols;
  
  //define the matrices store posterior samples
  arma::mat shiftsample(totdraw+1,J);
  shiftsample.fill(0.0);
 
  for (int w=0;w<J;w++){
     RNGScope scope; // ensure RNG gets set/reset
    shiftsample(0,w) = R::runif(-1*phi,phi);
  }
  
  //tausample is a vector
  arma::mat  tausample(totdraw+1,1);
  tausample.fill(1.0);
  tausample(0,0) = 1.0;
  
  //matrix for probabilities
  arma::mat probsample(totdraw+1,C);
  probsample.fill(0.0);
  for (int ii=0;ii<C;ii++){
     probsample(0,ii) = 1.0/C;
  }
  //matrix for membership
  arma::mat membersample(totdraw+1,J);
  membersample.fill(0.0);
  NumericVector initialprob = rep(1.0/C,C);
  for (int g=0;g<J; g++){
      membersample(0,g) = position1(oneMultinomC(initialprob));

  } 
  //matrix for factor a
  arma::mat factorasample(totdraw+1,J);
  factorasample.fill(1.0);
  
  //matrix for warping
  arma::mat warpingsample((totdraw+1)*J,M);
  warpingsample.fill(1.0/M);
  
   //matrix for spline coefficients
  arma::mat betasample((totdraw+1)*C,nbasis);
  betasample.fill(1.0);
  
 


 
 Rcpp::List betanulllist(C);
    for (int h=0;h<C; h++){
   
    arma::mat temp22(1,nbasis);
    temp22.fill(0.0);
     
     for (int s=0;s<nbasis;s++){
        RNGScope scope; // ensure RNG gets set/reset
     temp22(0,s)=R::rnorm(0,1);
   }
    betanulllist[h] = temp22;
 }
 
 
 
 //start main loop here
 for (int i=2; i<totdraw+2;i++){//i
   
   // member.sample[i,] <- member.sample[i-1,]
   //for (int m1=0; m1<J)
   
   //-----------------------------------------------------------
   //sample warping function

 for (int j=1; j<J+1;j++){//f
      double factoratemp = factorasample(i-2,j-1);
      arma::mat warpingveccurrent = warpingsample.row((i-2)*J+j-1);
      
    for (int k=1;k<M;k++){//r
         double warpingproposed = R::rnorm(warpingveccurrent(0,k-1),warpingsd); //may need to define double warpingproposed
   
    
    if (warpingproposed>0){//e
      arma::mat warpingvecproposed = warpingveccurrent;
      warpingvecproposed(0,k-1) = warpingproposed;
      warpingvecproposed(0,M-1) = 1-psum(warpingvecproposed);
    
    if (warpingvecproposed(0,M-1)>0){//b
         int q = membersample(i-2,j-1);
         arma::mat eval1(1,K);
         eval1.fill(1);
         arma::mat eval2(1,K);
         eval2.fill(1);

         for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t1 = evalbasis(linearapprox(temp,cumsumC(add0(warpingvecproposed))),degree, nbasis)*betasample.row((i-2)*C+q-1).t();
          double tt1 = t1(0,0);
          arma::mat t2 = evalbasis(linearapprox(temp,cumsumC(add0(warpingveccurrent))),degree, nbasis)*betasample.row((i-2)*C+q-1).t();
          double tt2 = t2(0,0); 
          eval1(0,f-1) = dat(j-1,f-1)-factoratemp*tt1-shiftsample(i-2,j-1);
          eval2(0,f-1) = dat(j-1,f-1)-factoratemp*tt2-shiftsample(i-2,j-1);

    }//g
    double logacceptratio = -0.5*tausample(i-2,0)*normsq(eval1)+ (alpha-1)*(log(warpingproposed)+log(warpingvecproposed(0,M-1)))+0.5*tausample(i-2,0)*normsq(eval2)-
                           (alpha-1)*(log(warpingveccurrent(0,k-1))+log(warpingveccurrent(0,M-1)));
                            RNGScope scope; // ensure RNG gets set/reset
        double rand = R::runif(0,1);
    if (logacceptratio>log(rand)){ //a
       warpingveccurrent = warpingvecproposed;
    }//a
}//b
   

}//e

}//r
warpingsample.row((i-1)*J+j-1) = warpingveccurrent;
 }//f
 
 
 //------------------------------------------------
 //sample cluster member EXACTLY
  for (int j=1; j<J+1;j++){ 
    double factoratemp = factorasample(i-2,j-1);
    NumericVector tempprobs(C);
    double sumtemp = 0;
    for (int jj=1; jj<C+1;jj++){
       arma::mat eval3(1,K);
       eval3.fill(1);
       for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t3 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*betasample.row((i-2)*C+jj-1).t();
          double tt3 = t3(0,0);
          eval3(0,f-1) = dat(j-1,f-1)-factoratemp*tt3-shiftsample(i-2,j-1);
    }//g
      double q = probsample(i-2,jj-1)*exp(-0.5*tausample(i-2,0)*normsq(eval3));
      tempprobs(jj-1)=q;
      sumtemp = sumtemp +q;
    }
    
    for (int jj=1; jj<C+1;jj++){
       tempprobs(jj-1) =  tempprobs(jj-1)/sumtemp;
        
    }
  membersample(i-1,j-1) = position1(oneMultinomC(tempprobs));
  //membersample.row(i-1)=membersample.row(i-2);
 
  }

//-----------------------------------------------
//sample probability
  probsample.row(i-1)=rdiri(probparam(membersample.row(i-1),eta,C));
 //probsample.row(i-1) = probsample.row(0);
 //---------------------------------------------

//make switch

Rcpp::List clustermemberlist(C);
Rcpp::List switchlist(C);
if (i%wait==0){
   
  IntegerVector Cluster(C);
  for (int mm=0;mm<C;mm++){Cluster(mm)=mm+1;}
  IntegerVector ordering = csample(Cluster,C,false);
  ordering = add1element(ordering,ordering(0));
  
  for (int m=1;m<C+1;m++){
    
    clustermemberlist[m-1]=whichC(arma2vec(membersample.row(i-1)),m);
    NumericVector aa(2);
    aa(0)=1;
    NumericVector tt = clustermemberlist[m-1];
     if (tt.size()!=0){
    aa(1)= tt.size()*percent;
    switchlist[m-1]=csample( clustermemberlist[m-1],max(aa),false);
    clustermemberlist[m-1] = deletevec(clustermemberlist[m-1], switchlist[m-1]);
     }
  }
  
  for (int n=1;n<C+1;n++){
     int from = ordering(n-1);
     int to = ordering(n);
    clustermemberlist[to-1] = combinevec(clustermemberlist[to-1],switchlist[from-1]);
  }

}else{ for (int m=1;m<C+1;m++){
    clustermemberlist[m-1]=whichC(arma2vec(membersample.row(i-1)),m);
    
  }
  }


//----------------------------------------------------------------
 //sample betas 
 for (int m=1;m<C+1;m++){
   IntegerVector tt = clustermemberlist[m-1];
    if (tt.size()==0){
      betasample.row((i-1)*C+m-1) = betasample.row((i-2)*C+m-1);
    }else{
     
     arma::mat D(nbasis,1);
     D.fill(0);
     arma::mat A(nbasis,nbasis);
     A.fill(0);
   for (int h=1;h<tt.size()+1;h++){
       arma::mat Beval(K,nbasis);
       Beval.fill(0);  
       for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          Beval.row(f-1) = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+tt(h-1)-1)))),degree, nbasis);
           
    }//g
       A = A + factorasample(i-2,tt(h-1)-1)*factorasample(i-2,tt(h-1)-1)*Beval.t()*Beval;
       arma::mat mattemp(1,K); mattemp.fill(1);
       arma::mat difftemp = dat.row(tt(h-1)-1)-shiftsample(i-2,tt(h-1)-1)*mattemp(1);
       D = D + factorasample(i-2,tt(h-1)-1)*Beval.t()*difftemp.t();
       
   }
    double  precisiontemp = tausample(i-2,0);
   arma::mat betanulltemp = betanulllist[m-1];
    A = precisiontemp*A + SigmaB.i();
    D = precisiontemp*D + SigmaB.i()*betanulltemp.t();
    arma::mat postmean = A.i()*D;  arma::mat postvar = A.i();
    betasample.row((i-1)*C+m-1) = mvrnormArma(1,postmean,postvar);
     }

 }

 
 

//sample shift by truncated normal
for (int j=1;j<J+1;j++){
       double factoratemp = factorasample(i-2,j-1);
       int k = membersample(i-1,j-1);
       arma::mat eval5(1,K);
       eval5.fill(1);
      for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t5 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*betasample.row((i-1)*C+k-1).t();
          double tt5 = t5(0,0);       
          eval5(0,f-1) = dat(j-1,f-1)-factoratemp*tt5;

    }//g  
      double tempmean = mean(arma2vecnum(eval5));//using rcpp sugar mean
      //double Ktemp=K;
       double tempvar = (1.0/tausample(i-2,0))/K;
      shiftsample(i-1,j-1) = trun_rnorm(tempmean,sqrt(tempvar),-phi,phi);
 
   }



//sample tau by gamma
 double tmp = 0;

for (int j=1;j<J+1;j++){
   double factoratemp = factorasample(i-2,j-1);
    
        arma::mat eval6(1,K);
        eval6.fill(1.0);
        int k = membersample(i-1,j-1);

        for (int f = 1; f<K+1;f++){//g
          double temp2 = clocktime(f-1); 
          arma::mat tempbeta=betasample.row((i-1)*C+k-1);
          arma::mat t6 = evalbasis(linearapprox(temp2,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*tempbeta.t();
          double tt6 = t6(0,0);
          eval6(0,f-1) = dat(j-1,f-1)-factoratemp*tt6-shiftsample(i-1,j-1);
 
    }//g
    
       tmp = tmp + normsq(eval6);
      
  }

   
      tausample(i-1,0) = R::rgamma(K/2.0* J+kappa,1.0/(0.5*tmp + theta));  
      //tausample(i-1,0) = 1.0;

//----------------------------------------
//sample shrinking/streching factor

for (int j=1;j<J+1;j++){
      int k = membersample(i-1,j-1);
      NumericVector tempvec1(K);
       NumericVector tempvec2(K);
       
      for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat tempbeta=betasample.row((i-1)*C+k-1);
          arma::mat t6 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*tempbeta.t();
          double tt6 = t6(0,0);
          tempvec1(f-1) = tt6*(dat(j-1,f-1)-shiftsample(i-1,j-1));
          tempvec2(f-1) = tt6*tt6;
    }//g
    double temp1 = sum(tempvec1); 
    double temp2 = sum(tempvec2);
    double tempmean = (tausample(i-1,0)*temp1 + 1.0/sigmaasq)/(1.0/sigmaasq+tausample(i-1,0)*temp2);
    double tempvar = 1.0/(1.0/sigmaasq+tausample(i-1,0)*temp2);
    factorasample(i-1,j-1) = R::rnorm(tempmean,sqrt(tempvar));
    //factorasample(i-1,j-1) = 1.0;
}


}//i
Rcpp::List results =
      Rcpp::List::create(Rcpp::Named("warpingsample")=warpingsample,
                         Rcpp::Named("membersample")= membersample,
                         Rcpp::Named("betasample")=betasample,
                         Rcpp::Named("probsample")=probsample,
                         Rcpp::Named("factorasample")=factorasample,
                         Rcpp::Named("tausample")=tausample,
                         Rcpp::Named("shiftsample")=shiftsample);

 return(results);
}//func





//BayesRegCluster serves as a function to pick the best C and for burn-in period only!
// [[Rcpp::export]]
List BayesRegCluster(const arma::mat dat, int totdraw,int Z,int M,double kappa,double theta, double alpha, double phi, 
                    double sigmaasq, double eta,  int nbasis, int degree, NumericVector clocktime,arma::mat SigmaB,int wait,double percent){//func
  //J is number of obs
  const int J = dat.n_rows;
  //K is number of measurements
  const int K = dat.n_cols;
  
  //define the matrices store posterior samples
  arma::mat shiftsample(totdraw+1,J);
  shiftsample.fill(0.0);
 
  for (int w=0;w<J;w++){
     RNGScope scope; // ensure RNG gets set/reset
    shiftsample(0,w) = R::runif(-1.0*phi,phi);
  }
  
  //tausample is a vector
  arma::mat  tausample(totdraw+1,1);
  tausample.fill(1.0);
  tausample(0,0) = 1.0;
  
  //matrix for probabilities
  arma::mat probsample(totdraw+1,Z);
  probsample.fill(0.0);
  for (int ii=0;ii<Z;ii++){
     probsample(0,ii) = 1.0/Z;
  }
  //matrix for membership
  arma::mat membersample(totdraw+1,J);
  membersample.fill(0.0);
  NumericVector initialprob = rep(1.0/Z,Z);
  for (int g=0;g<J; g++){
      membersample(0,g) = position1(oneMultinomC(initialprob));

  } 
  //matrix for factor a
  arma::mat factorasample(totdraw+1,J);
  factorasample.fill(1.0);
  
  //matrix for warping
  arma::mat warpingsample((totdraw+1)*J,M);
  warpingsample.fill(1.0/M);
  
   //matrix for spline coefficients
  arma::mat betasample((totdraw+1)*Z,nbasis);
  betasample.fill(1.0);
  
 
//likelihoodsample is a row vector
  arma::mat  likelihoodsample(1,totdraw+1);
  likelihoodsample.fill(0.0);
  
  //clusternumbersample is a row vector
  arma::mat  clusternumbersample(1,totdraw+1);
  clusternumbersample.fill(Z);
 
 Rcpp::List betanulllist(Z);
    for (int h=0;h<Z; h++){
   
    arma::mat temp22(1,nbasis);
    temp22.fill(0.0);
     
     for (int s=0;s<nbasis;s++){
        RNGScope scope; // ensure RNG gets set/reset
     temp22(0,s)=R::rnorm(0,1);
   }
    betanulllist[h] = temp22;
 }
 
 
//*********************************
IntegerVector tempset = tempfunc(totdraw,wait);

//*********************************
 int ff = 0;
 
 //start main loop here
 for (int i=2; i<totdraw+2;i++){//i
   
   // member.sample[i,] <- member.sample[i-1,]
   //for (int m1=0; m1<J)
   
   //-----------------------------------------------------------
   //sample warping function
 int C=clusternumbersample(0,i-2);
 for (int j=1; j<J+1;j++){//f
      double factoratemp = factorasample(i-2,j-1);
      arma::mat warpingveccurrent = warpingsample.row((i-2)*J+j-1);
      
    for (int k=1;k<M;k++){//r
         double warpingproposed = R::rnorm(warpingveccurrent(0,k-1),0.07); //may need to define double warpingproposed
   
    
    if (warpingproposed>0){//e
      arma::mat warpingvecproposed = warpingveccurrent;
      warpingvecproposed(0,k-1) = warpingproposed;
      warpingvecproposed(0,M-1) = 1-psum(warpingvecproposed);
    
    if (warpingvecproposed(0,M-1)>0){//b
         int q = membersample(i-2,j-1);
         arma::mat eval1(1,K);
         eval1.fill(1);
         arma::mat eval2(1,K);
         eval2.fill(1);

         for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t1 = evalbasis(linearapprox(temp,cumsumC(add0(warpingvecproposed))),degree, nbasis)*betasample.row((i-2)*Z+q-1).t();
          double tt1 = t1(0,0);
          arma::mat t2 = evalbasis(linearapprox(temp,cumsumC(add0(warpingveccurrent))),degree, nbasis)*betasample.row((i-2)*Z+q-1).t();
          double tt2 = t2(0,0); 
          eval1(0,f-1) = dat(j-1,f-1)-factoratemp*tt1-shiftsample(i-2,j-1);
          eval2(0,f-1) = dat(j-1,f-1)-factoratemp*tt2-shiftsample(i-2,j-1);

    }//g
    double logacceptratio = -0.5*tausample(i-2,0)*normsq(eval1)+ (alpha-1)*(log(warpingproposed)+log(warpingvecproposed(0,M-1)))+0.5*tausample(i-2,0)*normsq(eval2)-
                           (alpha-1)*(log(warpingveccurrent(0,k-1))+log(warpingveccurrent(0,M-1)));
                            RNGScope scope; // ensure RNG gets set/reset
        double rand = R::runif(0,1);
    if (logacceptratio>log(rand)){ //a
       warpingveccurrent = warpingvecproposed;
    }//a
}//b
   

}//e

}//r
warpingsample.row((i-1)*J+j-1) = warpingveccurrent;
 }//f
 

 //------------------------------------------------
 //sample cluster member EXACTLY
  for (int j=1; j<J+1;j++){ 
    double factoratemp = factorasample(i-2,j-1);
    NumericVector tempprobs(C);
    double sumtemp = 0;
    for (int jj=1; jj<C+1;jj++){
       arma::mat eval3(1,K);
       eval3.fill(1);
       for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t3 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*betasample.row((i-2)*Z+jj-1).t();
          double tt3 = t3(0,0);
          eval3(0,f-1) = dat(j-1,f-1)-factoratemp*tt3-shiftsample(i-2,j-1);
    }//g
      double q = probsample(i-2,jj-1)*exp(-0.5*tausample(i-2,0)*normsq(eval3));
      tempprobs(jj-1)=q;
      sumtemp = sumtemp +q;
    }
    
    for (int jj=1; jj<C+1;jj++){
       tempprobs(jj-1) =  tempprobs(jj-1)/sumtemp;
        
    }
  membersample(i-1,j-1) = position1(oneMultinomC(tempprobs)); 
  }
  
//find the nonempty clusters
ff = ff + 1;
IntegerVector nonzerocluster = stl_sort(unique(arma2vec(membersample.row(i-1))));
int Ztemp = nonzerocluster.size();
clusternumbersample(0,i-1) = Ztemp;

 
int sizeC = stl_sort(unique(arma2vec(clusternumbersample.submat(0,0,0,i-1)))).size();

//d=# of cluster(s) dropped at current iteration
int d = clusternumbersample(0,i-1)-clusternumbersample(0,i-2);
//whichC returns a IntegerVector of positions
if ((d==-1)&(sizeC>2)){
   ff = 0;
  //X.submat( first_row, first_col, last_row, last_col )
  //allnumber = recorded # of clusters up to current iteration
   IntegerVector allnumber =  stl_sort(unique(arma2vec(clusternumbersample.submat(0,0,0,i-1))));
   int previousnumber = allnumber[2];
  //set1 is the iterations with number of clusters equals to last C (# of clusters before drop)
  IntegerVector set1 = whichC(arma2vec(clusternumbersample),clusternumbersample(0,i-2));
  //set2 is the iterations with number of clusters equals to current (C-2)
  IntegerVector set2 = whichC(arma2vec(clusternumbersample),previousnumber);
  
  //only compare the avg. likelihood when the chain stay in both session long enough (>30)
  if((set1.size()>30)&(set2.size()>30)){//a
     int startiter1;
     int startiter2;
     if((1+sum(diff(set1)))==set1.size()){
        startiter1 = 1;
     }else{
       IntegerVector tempvec(set1.size());
       for (int k=0;k<set1.size();k++){tempvec[k]=k+1;}
       //only take the consecutive block; startiter1 is the starting pt of last consecutive block.
        startiter1 = deletevec(tempvec,whichC(diff(set1),1))[deletevec(tempvec,whichC(diff(set1),1)).size()-2]+1;  
     }
     int enditer1 = set1.size();
     IntegerVector temset1 = set1[(set1>(set1[startiter1-1]-1))&(set1<(set1[enditer1-1]+1))];
     //do the same to the member associated with previousnumber.
     if((1+sum(diff(set2)))==set2.size()){
         startiter2 = 1;
     }else{
     IntegerVector tempvec(set2.size());
       for (int k=0;k<set2.size();k++){tempvec[k]=k+1;}
       //only take the consecutive block
        startiter2 = deletevec(tempvec,whichC(diff(set2),1))[deletevec(tempvec,whichC(diff(set2),1)).size()-2]+1;  
     }
     int enditer2 = set2.size();
     
     IntegerVector temset2 = set2[(set2>(set2[startiter2-1]-1))&(set2<(set2[enditer2-1]+1))];
     //must have more than 30 iterations in each block to compare the likelihood   
     //tempset ignores 5 iterations right after the switch 
     if((temset1.size()>30)&(temset2.size()>30)){
       IntegerVector pick1 = intersect(temset1,tempset);
       double meannew = mean(vecsubsetter(arma2vecnum(likelihoodsample),pick1)); 
       IntegerVector pick2 = intersect(temset2,tempset);
       double meanold = mean(vecsubsetter(arma2vecnum(likelihoodsample),pick2)); 
       if(meanold > meannew){
         IntegerVector tempvec(set2.size());
         for (int k=0;k<set2.size();k++){tempvec[k]=k+1;}
         //tempindex is the iteration, where the most recent block of C=previousnumber happened. 
         int tempindex = deletevec(tempvec,whichC(diff(set2),1))[deletevec(tempvec,whichC(diff(set2),1)).size()-2]+1; 
         int oldindex = whichC(arma2vec(clusternumbersample),previousnumber)(tempindex);
         membersample.row(i-1) =membersample.row(oldindex-1);
         clusternumbersample(0,i-1) = previousnumber;
       }
     }
     
  }//a
}

if((d<-1)&(i>2)&(sizeC>2)){
  ff = 0;
  IntegerVector allnumber =  stl_sort(unique(arma2vec(clusternumbersample.submat(0,0,0,i-1))));
   int previousnumber = allnumber[2];
  IntegerVector set1 = whichC(arma2vec(clusternumbersample),clusternumbersample(0,i-2));
  IntegerVector set2 = whichC(arma2vec(clusternumbersample),previousnumber);
  int startiter1;
  int startiter2;
   if((1+sum(diff(set1)))==set1.size()){
        startiter1 = 1;
     }else{
       IntegerVector tempvec(set1.size());
       for (int k=0;k<set1.size();k++){tempvec[k]=k+1;}
       //only take the consecutive block
         startiter1 = deletevec(tempvec,whichC(diff(set1),1))[deletevec(tempvec,whichC(diff(set1),1)).size()-2]+1;  
     }
     int enditer1 = set1.size();
     IntegerVector temset1 = set1[(set1>(set1[startiter1-1]-1))&(set1<(set1[enditer1-1]+1))];
     //do the same to the member associated with previousnumber.
     if((1+sum(diff(set2)))==set2.size()){
         startiter2 = 1;
     }else{
     IntegerVector tempvec(set2.size());
       for (int k=0;k<set2.size();k++){tempvec[k]=k+1;}
       //only take the consecutive block
        startiter2 = deletevec(tempvec,whichC(diff(set2),1))[deletevec(tempvec,whichC(diff(set2),1)).size()-2]+1;  
     }
     int enditer2 = set2.size();
     
     IntegerVector temset2 = set2[(set2>(set2[startiter2-1]-1))&(set2<(set2[enditer2-1]+1))];
     //must have more than 30 iterations in each block to compare the likelihood    
     if((temset1.size()>30)&(temset2.size()>30)){  
         IntegerVector tempvec(set1.size());
         for (int k=0;k<set1.size();k++){tempvec[k]=k+1;}
         //tempindex is the iteration, where the most recent block of C=previousnumber happened. 
         int tempindex = deletevec(tempvec,whichC(diff(set2),1))[deletevec(tempvec,whichC(diff(set2),1)).size()-2]+1; 
         int oldindex = whichC(arma2vec(clusternumbersample),clusternumbersample(0,i-2))(tempindex);
         membersample.row(i-1) =membersample.row(oldindex-1);
         clusternumbersample(0,i-1) = clusternumbersample(0,i-2);      
     }
}

IntegerVector allnumber =  stl_sort(unique(arma2vec(clusternumbersample.submat(0,0,0,i-1))));
int previousnumber = allnumber[1];
IntegerVector settemp = whichC(arma2vec(clusternumbersample),previousnumber);
  //X.submat( first_row, first_col, last_row, last_col )
if (ff>300 & settemp.size()>30& ff%50==0 ){//a
     IntegerVector allnumber =  stl_sort(unique(arma2vec(clusternumbersample.submat(0,0,0,i-1))));
   int previousnumber = allnumber[1];
  IntegerVector set1 = whichC(arma2vec(clusternumbersample),clusternumbersample(0,i-1));
  IntegerVector set2 = whichC(arma2vec(clusternumbersample),previousnumber);

     int startiter1;
     int startiter2;
     if((1+sum(diff(set1)))==set1.size()){
        startiter1 = 1;
     }else{
       IntegerVector tempvec(set1.size());
       for (int k=0;k<set1.size();k++){tempvec[k]=k+1;}
       //only take the consecutive block
        startiter1 = deletevec(tempvec,whichC(diff(set1),1))[deletevec(tempvec,whichC(diff(set1),1)).size()-2]+1;  
     }
     int enditer1 = set1.size();
    IntegerVector temset1 = set1[(set1>(set1[startiter1-1]-1))&(set1<(set1[enditer1-1]))];
       
     //do the same to the member associated with previousnumber.
     if((1+sum(diff(set2)))==set2.size()){
         startiter2 = 1;
     }else{
     IntegerVector tempvec(set2.size());
       for (int k=0;k<set2.size();k++){tempvec[k]=k+1;}
       
       //only take the consecutive block
        startiter2 = deletevec(tempvec,whichC(diff(set2),1))[deletevec(tempvec,whichC(diff(set2),1)).size()-2]+1;  
     }
     int enditer2 = set2.size();
     
     IntegerVector temset2 = set2[(set2>(set2[startiter2-1]-1))&(set2<(set2[enditer2-1]+1))];
      
       IntegerVector pick1 = intersect(temset1,tempset);
       //vecsubsetter takes the subset indexed by the 2nd arg.
       double meannew = mean(vecsubsetter(arma2vecnum(likelihoodsample),pick1)); 
        IntegerVector pick2 = intersect(temset2,tempset);
       double meanold = mean(vecsubsetter(arma2vecnum(likelihoodsample),pick2)); 
     if(meanold > meannew){
         IntegerVector tempvec(set2.size());
         for (int k=0;k<set2.size();k++){tempvec[k]=k+1;}
         //tempindex is the iterations, where the most recent block of C=previousnumber happened. 
         int tempindex = deletevec(tempvec,whichC(diff(set2),1))[deletevec(tempvec,whichC(diff(set2),1)).size()-2]+1; 
         int oldindex = whichC(arma2vec(clusternumbersample),previousnumber)(tempindex);
         membersample.row(i-1) =membersample.row(oldindex-1);
         clusternumbersample(0,i-1) = previousnumber;
       }  
                                  
}//a

//relabel according to nonzerocluster
//Ztemp is the current number of non-empty clusters after adjustment. 
Ztemp = clusternumbersample(0,i-1);
//stl_sort return a integer vector
IntegerVector tempvec = arma2vec(membersample.row(i-1));

for (int b=1;b<Ztemp+1;b++){ 
   IntegerVector tempsort = stl_sort(unique(tempvec));
   IntegerVector tempindex = whichC(tempvec,tempsort[b-1]);
   int n = tempindex.size();
   for (int k=0;k<n;k++){
     membersample(i-1,tempindex(k)-1) = b;
   }
}

//now the non-empty clusters change through iterations
  C = Ztemp;
  
//-----------------------------------------------
//sample probability
 arma::mat probtemp = rdiri(probparam(membersample.row(i-1),eta,C));
   for (int k=0;k<C;k++){
    probsample(i-1,k)= probtemp(0,k);
  } //probsample.row(i-1) = probsample.row(0);
 //---------------------------------------------

//make switch

Rcpp::List clustermemberlist(C);
Rcpp::List switchlist(C);
if (i%wait==0){
   
  IntegerVector Cluster(C);
  for (int mm=0;mm<C;mm++){Cluster(mm)=mm+1;}
  IntegerVector ordering = csample(Cluster,C,false);
  ordering = add1element(ordering,ordering(0));
  
  for (int m=1;m<C+1;m++){
    
    clustermemberlist[m-1]=whichC(arma2vec(membersample.row(i-1)),m);
    NumericVector aa(2);
    aa(0)=1;
    NumericVector tt = clustermemberlist[m-1];
     if (tt.size()!=0){
    aa(1)= tt.size()*percent;
    switchlist[m-1]=csample( clustermemberlist[m-1],max(aa),false);
    clustermemberlist[m-1] = deletevec(clustermemberlist[m-1], switchlist[m-1]);
     }
  }
  
  for (int n=1;n<C+1;n++){
     int from = ordering(n-1);
     int to = ordering(n);
    clustermemberlist[to-1] = combinevec(clustermemberlist[to-1],switchlist[from-1]);
  }

}else{ for (int m=1;m<C+1;m++){
    clustermemberlist[m-1]=whichC(arma2vec(membersample.row(i-1)),m);
    
  }
  }


//----------------------------------------------------------------
 //sample betas 
 for (int m=1;m<C+1;m++){
   IntegerVector tt = clustermemberlist[m-1];
 
     arma::mat D(nbasis,1);
     D.fill(0);
     arma::mat A(nbasis,nbasis);
     A.fill(0);
   for (int h=1;h<tt.size()+1;h++){
       arma::mat Beval(K,nbasis);
       Beval.fill(0);  
       for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          Beval.row(f-1) = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+tt(h-1)-1)))),degree, nbasis);
           
    }//g
       A = A + factorasample(i-2,tt(h-1)-1)*factorasample(i-2,tt(h-1)-1)*Beval.t()*Beval;
       arma::mat mattemp(1,K); mattemp.fill(1);
       arma::mat difftemp = dat.row(tt(h-1)-1)-shiftsample(i-2,tt(h-1)-1)*mattemp(1);
       D = D + factorasample(i-2,tt(h-1)-1)*Beval.t()*difftemp.t();
       
   }
    double  precisiontemp = tausample(i-2,0);
   arma::mat betanulltemp = betanulllist[m-1];
    A = precisiontemp*A + SigmaB.i();
    D = precisiontemp*D + SigmaB.i()*betanulltemp.t();
    arma::mat postmean = A.i()*D;  arma::mat postvar = A.i();
    betasample.row((i-1)*Z+m-1) = mvrnormArma(1,postmean,postvar);
     

 }

 
 

//sample shift by truncated normal
for (int j=1;j<J+1;j++){
       double factoratemp = factorasample(i-2,j-1);
       int k = membersample(i-1,j-1);
       arma::mat eval5(1,K);
       eval5.fill(1);
      for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t5 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*betasample.row((i-1)*Z+k-1).t();
          double tt5 = t5(0,0);       
          eval5(0,f-1) = dat(j-1,f-1)-factoratemp*tt5-shiftsample(i-2,j-1);

    }//g  
      double tempmean = mean(arma2vecnum(eval5));//using rcpp sugar mean
      //double Ktemp=K;
       double tempvar = (1.0/tausample(i-2,0))/K;
      shiftsample(i-1,j-1) = trun_rnorm(tempmean,sqrt(tempvar),-phi,phi);
      //shiftsample(i-1,j-1) = 0;
   }



//sample tau by gamma
 double tmp = 0;

for (int j=1;j<J+1;j++){
   double factoratemp = factorasample(i-2,j-1);
    
        arma::mat eval6(1,K);
        eval6.fill(1.0);
        int k = membersample(i-1,j-1);

        for (int f = 1; f<K+1;f++){//g
          double temp2 = clocktime(f-1); 
          arma::mat tempbeta=betasample.row((i-1)*Z+k-1);
          arma::mat t6 = evalbasis(linearapprox(temp2,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*tempbeta.t();
          double tt6 = t6(0,0);
          eval6(0,f-1) = dat(j-1,f-1)-factoratemp*tt6-shiftsample(i-1,j-1);
 
    }//g
    
       tmp = tmp + normsq(eval6);
      
  }

   
      tausample(i-1,0) = R::rgamma(K/2.0* J+kappa,1.0/(0.5*tmp + theta));  
      //tausample(i-1,0) = 1.0;

//----------------------------------------
//sample shrinking/streching factor

for (int j=1;j<J+1;j++){
      int k = membersample(i-1,j-1);
      NumericVector tempvec1(K);
       NumericVector tempvec2(K);
       
      for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat tempbeta=betasample.row((i-1)*Z+k-1);
          arma::mat t6 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*tempbeta.t();
          double tt6 = t6(0,0);
          tempvec1(f-1) = tt6*(dat(j-1,f-1)-shiftsample(i-1,j-1));
          tempvec2(f-1) = tt6*tt6;
    }//g
    double temp1 = sum(tempvec1); 
    double temp2 = sum(tempvec2);
    double tempmean = (tausample(i-1,0)*temp1 + 1.0/sigmaasq)/(1.0/sigmaasq+tausample(i-1,0)*temp2);
    double tempvar = 1.0/(1.0/sigmaasq+tausample(i-1,0)*temp2);
    factorasample(i-1,j-1) = R::rnorm(tempmean,sqrt(tempvar));
    //factorasample(i-1,j-1) = 1.0;
}

//calculate the likelihood clustermemberlist
double logL = 0.0;
for (int l=1;l<Ztemp+1;l++){
  IntegerVector tt = clustermemberlist[l-1];
   for (int h=1;h<tt.size()+1;h++){
     arma::mat eval1(1,K);
     eval1.fill(1);
     for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t1 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+tt(h-1)-1)))),degree, nbasis)*betasample.row((i-1)*Z+l-1).t();
          double tt1 = t1(0,0);
    
          eval1(0,f-1) = dat(tt(h-1)-1,f-1)-factorasample(i-1,tt(h-1)-1)*tt1-shiftsample(i-1,tt(h-1)-1);
    }//g
     logL = logL + normsq(eval1);
   }
  
}
  logL = (-0.5)*tausample(i-1,0)*logL+ K*J/2.0*log(tausample(i-1,0));
  likelihoodsample(0,i-1) = logL;


//Rcpp::Rcout << i << std::endl;
//Rcpp::Rcout << clusternumbersample(0,i-1) << std::endl;


}//i
Rcpp::List results =
      Rcpp::List::create(Rcpp::Named("warpingsample")=warpingsample,
                         Rcpp::Named("membersample")= membersample,
                         Rcpp::Named("betasample")=betasample,
                         Rcpp::Named("probsample")=probsample,
                         Rcpp::Named("factorasample")=factorasample,
                         Rcpp::Named("tausample")=tausample,
                         Rcpp::Named("shiftsample")=shiftsample,
                         Rcpp::Named("clusternumbersample")=clusternumbersample,
                         Rcpp::Named("likelihoodsample")=likelihoodsample);

 return(results);
}//func

//the function MCMC give the standard MCMC samplong w/o switching membership or changing C.
// [[Rcpp::export]]
List MCMC(const arma::mat dat, int totdraw,int C,int M,double kappa,double theta, double alpha, double phi, 
                    double sigmaasq, double eta,int nbasis, int degree, NumericVector clocktime,arma::mat SigmaB,
                    arma::mat factoraest,arma::mat warpingest,arma::mat shiftest,arma::mat memberest,double tauest,
                    arma::mat probest,arma::mat betaest){//func
  //J is number of obs
  const int J = dat.n_rows;
  //K is number of measurements
  const int K = dat.n_cols;
  
  //define the matrices store posterior samples
  arma::mat shiftsample(totdraw+1,J);
  shiftsample.row(0) = shiftest;
  
  //tausample is a vector
  arma::mat  tausample(totdraw+1,1);
  tausample.fill(1.0);
  tausample(0,0) = tauest;
  
  //matrix for probabilities
  arma::mat probsample(totdraw+1,C);
  probsample.row(0)=probest;
  
  //matrix for membership
  arma::mat membersample(totdraw+1,J);
  membersample.row(0)=memberest;
  
  //matrix for factor a
  arma::mat factorasample(totdraw+1,J);
  factorasample.row(0) = factoraest;
  
    //X.submat( first_row, first_col, last_row, last_col )
  //matrix for warping
  arma::mat warpingsample((totdraw+1)*J,M);
  warpingsample.submat(0,0,J-1,M-1)=warpingest;
  
   //matrix for spline coefficients
  arma::mat betasample((totdraw+1)*C,nbasis);
  betasample.submat(0,0,C-1,nbasis-1)=betaest;
  
 

 Rcpp::List betanulllist(C);
    for (int h=0;h<C; h++){
   
    arma::mat temp22(1,nbasis);
    temp22.fill(0.0);
     
     for (int s=0;s<nbasis;s++){
        RNGScope scope; // ensure RNG gets set/reset
     temp22(0,s)=R::rnorm(0,1);
   }
    betanulllist[h] = temp22;
 }
 

 
 //start main loop here
 for (int i=2; i<totdraw+2;i++){//i
   
   // member.sample[i,] <- member.sample[i-1,]
   //for (int m1=0; m1<J)
   
   //-----------------------------------------------------------
   //sample warping function

 for (int j=1; j<J+1;j++){//f
      double factoratemp = factorasample(i-2,j-1);
      arma::mat warpingveccurrent = warpingsample.row((i-2)*J+j-1);
      
    for (int k=1;k<M;k++){//r
         double warpingproposed = R::rnorm(warpingveccurrent(0,k-1),0.07); //may need to define double warpingproposed
   
    
    if (warpingproposed>0){//e
      arma::mat warpingvecproposed = warpingveccurrent;
      warpingvecproposed(0,k-1) = warpingproposed;
      warpingvecproposed(0,M-1) = 1-psum(warpingvecproposed);
    
    if (warpingvecproposed(0,M-1)>0){//b
         int q = membersample(i-2,j-1);
         arma::mat eval1(1,K);
         eval1.fill(1);
         arma::mat eval2(1,K);
         eval2.fill(1);

         for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t1 = evalbasis(linearapprox(temp,cumsumC(add0(warpingvecproposed))),degree, nbasis)*betasample.row((i-2)*C+q-1).t();
          double tt1 = t1(0,0);
          arma::mat t2 = evalbasis(linearapprox(temp,cumsumC(add0(warpingveccurrent))),degree, nbasis)*betasample.row((i-2)*C+q-1).t();
          double tt2 = t2(0,0); 
          eval1(0,f-1) = dat(j-1,f-1)-factoratemp*tt1-shiftsample(i-2,j-1);
          eval2(0,f-1) = dat(j-1,f-1)-factoratemp*tt2-shiftsample(i-2,j-1);

    }//g
    double logacceptratio = -0.5*tausample(i-2,0)*normsq(eval1)+ (alpha-1)*(log(warpingproposed)+log(warpingvecproposed(0,M-1)))+0.5*tausample(i-2,0)*normsq(eval2)-
                           (alpha-1)*(log(warpingveccurrent(0,k-1))+log(warpingveccurrent(0,M-1)));
                            RNGScope scope; // ensure RNG gets set/reset
        double rand = R::runif(0,1);
    if (logacceptratio>log(rand)){ //a
       warpingveccurrent = warpingvecproposed;
    }//a
}//b
   
}//e

}//r
warpingsample.row((i-1)*J+j-1) = warpingveccurrent;
 }//f
 
 
 //------------------------------------------------
 //sample cluster member EXACTLY
  for (int j=1; j<J+1;j++){ 
    double factoratemp = factorasample(i-2,j-1);
    NumericVector tempprobs(C);
    double sumtemp = 0;
    for (int jj=1; jj<C+1;jj++){
       arma::mat eval3(1,K);
       eval3.fill(1);
       for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t3 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*betasample.row((i-2)*C+jj-1).t();
          double tt3 = t3(0,0);
          eval3(0,f-1) = dat(j-1,f-1)-factoratemp*tt3-shiftsample(i-2,j-1);
    }//g
      double q = probsample(i-2,jj-1)*exp(-0.5*tausample(i-2,0)*normsq(eval3));
      tempprobs(jj-1)=q;
      sumtemp = sumtemp +q;
    }
    
    for (int jj=1; jj<C+1;jj++){
       tempprobs(jj-1) =  tempprobs(jj-1)/sumtemp;
        
    }
  membersample(i-1,j-1) = position1(oneMultinomC(tempprobs));
  //membersample.row(i-1)=membersample.row(i-2);
 
  }

//-----------------------------------------------
//sample probability
  probsample.row(i-1)=rdiri(probparam(membersample.row(i-1),eta,C));
 //probsample.row(i-1) = probsample.row(0);
 //---------------------------------------------

Rcpp::List clustermemberlist(C);
//find cluster member
 for (int m=1;m<C+1;m++){
    clustermemberlist[m-1]=whichC(arma2vec(membersample.row(i-1)),m);
    
  }
 


//----------------------------------------------------------------
 //sample betas 
 for (int m=1;m<C+1;m++){
   IntegerVector tt = clustermemberlist[m-1];
    if (tt.size()==0){
      betasample.row((i-1)*C+m-1) = betasample.row((i-2)*C+m-1);
    }else{
     
     arma::mat D(nbasis,1);
     D.fill(0);
     arma::mat A(nbasis,nbasis);
     A.fill(0);
   for (int h=1;h<tt.size()+1;h++){
       arma::mat Beval(K,nbasis);
       Beval.fill(0);  
       for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          Beval.row(f-1) = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+tt(h-1)-1)))),degree, nbasis);
           
    }//g
       A = A + factorasample(i-2,tt(h-1)-1)*factorasample(i-2,tt(h-1)-1)*Beval.t()*Beval;
       arma::mat mattemp(1,K); mattemp.fill(1);
       arma::mat difftemp = dat.row(tt(h-1)-1)-shiftsample(i-2,tt(h-1)-1)*mattemp(1);
       D = D + factorasample(i-2,tt(h-1)-1)*Beval.t()*difftemp.t();
       
   }
    double  precisiontemp = tausample(i-2,0);
   arma::mat betanulltemp = betanulllist[m-1];
    A = precisiontemp*A + SigmaB.i();
    D = precisiontemp*D + SigmaB.i()*betanulltemp.t();
    arma::mat postmean = A.i()*D;  arma::mat postvar = A.i();
    betasample.row((i-1)*C+m-1) = mvrnormArma(1,postmean,postvar);
     }

 }

 
//sample shift by truncated normal
for (int j=1;j<J+1;j++){
       double factoratemp = factorasample(i-2,j-1);
       int k = membersample(i-1,j-1);
       arma::mat eval5(1,K);
       eval5.fill(1);
      for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat t5 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*betasample.row((i-1)*C+k-1).t();
          double tt5 = t5(0,0);       
          eval5(0,f-1) = dat(j-1,f-1)-factoratemp*tt5-shiftsample(i-2,j-1);

    }//g  
      double tempmean = mean(arma2vecnum(eval5));//using rcpp sugar mean
      //double Ktemp=K;
       double tempvar = (1.0/tausample(i-2,0))/K;
      shiftsample(i-1,j-1) = trun_rnorm(tempmean,sqrt(tempvar),-phi,phi);
      //shiftsample(i-1,j-1) = 0;
   }



//sample tau by gamma
 double tmp = 0;

for (int j=1;j<J+1;j++){
   double factoratemp = factorasample(i-2,j-1);
    
        arma::mat eval6(1,K);
        eval6.fill(1.0);
        int k = membersample(i-1,j-1);

        for (int f = 1; f<K+1;f++){//g
          double temp2 = clocktime(f-1); 
          arma::mat tempbeta=betasample.row((i-1)*C+k-1);
          arma::mat t6 = evalbasis(linearapprox(temp2,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*tempbeta.t();
          double tt6 = t6(0,0);
          eval6(0,f-1) = dat(j-1,f-1)-factoratemp*tt6-shiftsample(i-1,j-1);
 
    }//g
    
       tmp = tmp + normsq(eval6);
      
  }

   
      tausample(i-1,0) = R::rgamma(K/2.0* J+kappa,1.0/(0.5*tmp + theta));  
      //tausample(i-1,0) = 1.0;

//----------------------------------------
//sample shrinking/streching factor

for (int j=1;j<J+1;j++){
      int k = membersample(i-1,j-1);
      NumericVector tempvec1(K);
       NumericVector tempvec2(K);
       
      for (int f = 1; f<K+1;f++){//g
          double temp = clocktime(f-1); 
          arma::mat tempbeta=betasample.row((i-1)*C+k-1);
          arma::mat t6 = evalbasis(linearapprox(temp,cumsumC(add0(warpingsample.row((i-1)*J+j-1)))),degree, nbasis)*tempbeta.t();
          double tt6 = t6(0,0);
          tempvec1(f-1) = tt6*(dat(j-1,f-1)-shiftsample(i-1,j-1));
          tempvec2(f-1) = tt6*tt6;
    }//g
    double temp1 = sum(tempvec1); 
    double temp2 = sum(tempvec2);
    double tempmean = (tausample(i-1,0)*temp1 + 1.0/sigmaasq)/(1.0/sigmaasq+tausample(i-1,0)*temp2);
    double tempvar = 1.0/(1.0/sigmaasq+tausample(i-1,0)*temp2);
    factorasample(i-1,j-1) = R::rnorm(tempmean,sqrt(tempvar));
    //factorasample(i-1,j-1) = 1.0;
}


}//i
Rcpp::List results =
      Rcpp::List::create(Rcpp::Named("warpingsample")=warpingsample,
                         Rcpp::Named("membersample")= membersample,
                         Rcpp::Named("betasample")=betasample,
                         Rcpp::Named("probsample")=probsample,
                         Rcpp::Named("factorasample")=factorasample,
                         Rcpp::Named("tausample")=tausample,
                         Rcpp::Named("shiftsample")=shiftsample);

 return(results);
}//func

 

