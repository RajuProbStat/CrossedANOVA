#### Simultaneous ordering in a two-way crossed ANOVA for real data sets
library(Iso)
#####
#' Function to evaluate the observed, critical and P-values of the test statistics R_1
#'
#' @param error Numeric, the level of significance/type-I error
#' @param a Numeric, the number of levels of the factor A
#' @param b Numeric, the number of levels of the factor B
#' @param cellsize matrix having a rows and b columns, and (i,j)-th element is the number of observations in the (i,j)-th cell
#' @param cellmean matrix having a rows and b columns, and (i,j)-th element is the sample mean of the observations in the (i,j)-th cell
#' @param cellvariance matrix having a rows and b columns, and (i,j)-th element is the sample variance of the observations in the (i,j)-th cell
#'
#' @return numeric vector consisting of test statistic value, critical value and P-value of the test R_1. Also, the decision of the test R_1.
#' @export
#'
#' @examples
#' R_1(0.05,3,2,rbind(c(35,35),c(56,56),c(51,51)), rbind(c(52.7429,99.9143),c(58.8750,126.7679),c(95.8039,261.5490)), rbind(c(67.2555,274.0807),c(222.8023,955.0542),c(815.3608,8597.5325)))

#Evaluate the observed test statistic value, critical value and P-value for the test statistic R_1 from a given data
R_1<-function(error,a,b,cellsize,cellmean,cellvariance){
  a<-a
  b<-b
  n<-cellsize
  sam_mean<-cellmean
  sam_var<-cellvariance
  error<-error
  fun_LRT<-function(a,b,cellsize,cellmean,cellvariance){
    a<-a
    b<-b
    n<-cellsize
    sam_mean<-cellmean
    sam_var<-cellvariance
    mu_hat<-sum(n*sam_mean)/sum(n)
    alpha_hat<-apply(sam_mean,1,mean)-mu_hat
    beta_hat<-apply(sam_mean,2,mean)-mu_hat
    nu_hat<-sam_mean-matrix(c(alpha_hat),nrow=a,ncol=b,byrow = FALSE)-matrix(c(beta_hat),nrow=a,ncol=b,byrow=TRUE)-mu_hat*matrix(1,a,b)
    ## MLEs of the parameters under Omega 0R
    beta_00A<-apply(sam_mean,2,mean)
    var_00A<-sam_var
    repeat
    {
      beta_10A<-apply(((n*sam_mean)/(var_00A)),2,sum)/apply(n/var_00A,2,sum)
      var_10A<-sam_var+(sam_mean-matrix(c(beta_10A),nrow=a,ncol=b,byrow=TRUE))*(sam_mean-matrix(c(beta_10A),nrow=a,ncol=b,byrow=TRUE))

      if(max(c(max(abs(beta_10A-beta_00A)),max(abs(var_10A-var_00A))))<0.000001)
      {
        break
      }
      beta_00A<-beta_10A
      var_00A<-var_10A
    }
    ## MLEs of the parameters under Omega R
    nu_0A<-nu_hat
    var_0A<-sam_var
    repeat
    {
      nu_1A<-array(NA,dim=c(a,b))
      for(j in 1:b)
      {
        nu_1A[,j]<-pava(sam_mean[,j],(n/var_0A)[,j], decreasing=FALSE, long.out=FALSE, stepfun=FALSE)
      }
      var_1A<-sam_var+(sam_mean-nu_1A)*(sam_mean-nu_1A)

      if(max(c(max(abs(nu_1A-nu_0A)),max(abs(var_1A-var_0A))))<0.000001)
      {
        break
      }
      var_0A<-var_1A
      nu_0A<-nu_1A
    }
    ## MLEs of the parameters under Omega 0C
    alpha_00B<-apply(sam_mean,1,mean)
    var_00B<-sam_var
    repeat
    {
      alpha_10B<-apply(((n*sam_mean)/(var_00B)),1,sum)/apply(n/var_00B,1,sum)
      var_10B<-sam_var+(sam_mean-matrix(c(alpha_10B),nrow=a,ncol=b,byrow=FALSE))*(sam_mean-matrix(c(alpha_10B),nrow=a,ncol=b,byrow=FALSE))
      if(max(c(max(abs(alpha_10B-alpha_00B)),max(abs(var_10B-var_00B))))<0.000001)
      {
        break
      }
      alpha_00B<-alpha_10B
      var_00B<-var_10B
    }
    ## MLEs of the parameters under Omega C
    nu_0B<-nu_hat
    var_0B<-sam_var
    repeat
    {
      nu_1B<-array(NA,dim=c(a,b))
      for(i in 1:a)
      {
        nu_1B[i,]<-pava(sam_mean[i,],(n/var_0B)[i,], decreasing=FALSE, long.out=FALSE, stepfun=FALSE)
      }
      var_1B<-sam_var+(sam_mean-nu_1B)*(sam_mean-nu_1B)

      if(max(c(max(abs(nu_1B-nu_0B)),max(abs(var_1B-var_0B))))<0.000001)
      {
        break
      }
      var_0B<-var_1B
      nu_0B<-nu_1B
    }
    LRT_value<-max(prod((var_1A/var_10A)^(n/2)),prod((var_1B/var_10B)^(n/2)))
    return(c(LRT_value,t(sam_var)))
  }
  values<-fun_LRT(a,b,cellsize,cellmean,cellvariance)
  value_LRT<-values[1]
  ##Bootstrap step
  B<-10000 #number of bootstrap samples
  value_LRT_boot<-rep(NA,B)
  for(r in 1:B)
  {
    set.seed(17*r)
    sam_var_boot<-matrix(values[2:(a*b+1)],nrow=a,ncol=b,byrow=TRUE)
    sam_mean<-array(NA,dim=c(a,b))
    sam_var<-array(NA,dim=c(a,b))
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        data<-rnorm(n[i,j],0,sqrt(sam_var_boot[i,j]))
        sam_mean[i,j]<-mean(data)
        sam_var[i,j]<-var(data)
      }
    }
    mu_hat<-sum(n*sam_mean)/sum(n)
    alpha_hat<-apply(sam_mean,1,mean)-mu_hat
    beta_hat<-apply(sam_mean,2,mean)-mu_hat
    nu_hat<-sam_mean-matrix(c(alpha_hat),nrow=a,ncol=b,byrow = FALSE)-matrix(c(beta_hat),nrow=a,ncol=b,byrow=TRUE)-mu_hat*matrix(1,a,b)
    ## MLEs of the parameters under Omega 0R for bootstrap samples
    beta_00A<-apply(sam_mean,2,mean)
    var_00A<-sam_var
    repeat
    {
      beta_10A<-apply(((n*sam_mean)/(var_00A)),2,sum)/apply(n/var_00A,2,sum)
      var_10A<-sam_var+(sam_mean-matrix(c(beta_10A),nrow=a,ncol=b,byrow=TRUE))*(sam_mean-matrix(c(beta_10A),nrow=a,ncol=b,byrow=TRUE))

      if(max(c(max(abs(beta_10A-beta_00A)),max(abs(var_10A-var_00A))))<0.000001)
      {
        break
      }
      beta_00A<-beta_10A
      var_00A<-var_10A
    }
    ## MLEs of the parameters under Omega R for bootstrap samples
    nu_0A<-nu_hat
    var_0A<-sam_var
    repeat
    {
      nu_1A<-array(NA,dim=c(a,b))
      for(j in 1:b)
      {
        nu_1A[,j]<-pava(sam_mean[,j],(n/var_0A)[,j], decreasing=FALSE, long.out=FALSE, stepfun=FALSE)
      }
      var_1A<-sam_var+(sam_mean-nu_1A)*(sam_mean-nu_1A)

      if(max(c(max(abs(nu_1A-nu_0A)),max(abs(var_1A-var_0A))))<0.000001)
      {
        break
      }
      var_0A<-var_1A
      nu_0A<-nu_1A
    }
    ## MLEs of the parameters under Omega 0C for bootstrap samples
    alpha_00B<-apply(sam_mean,1,mean)
    var_00B<-sam_var
    repeat
    {
      alpha_10B<-apply(((n*sam_mean)/(var_00B)),1,sum)/apply(n/var_00B,1,sum)
      var_10B<-sam_var+(sam_mean-matrix(c(alpha_10B),nrow=a,ncol=b,byrow=FALSE))*(sam_mean-matrix(c(alpha_10B),nrow=a,ncol=b,byrow=FALSE))
      if(max(c(max(abs(alpha_10B-alpha_00B)),max(abs(var_10B-var_00B))))<0.000001)
      {
        break
      }
      alpha_00B<-alpha_10B
      var_00B<-var_10B
    }
    ## MLEs of the parameters under Omega C for bootstrap samples
    nu_0B<-nu_hat
    var_0B<-sam_var
    repeat
    {
      nu_1B<-array(NA,dim=c(a,b))
      for(i in 1:a)
      {
        nu_1B[i,]<-pava(sam_mean[i,],(n/var_0B)[i,], decreasing=FALSE, long.out=FALSE, stepfun=FALSE)
      }
      var_1B<-sam_var+(sam_mean-nu_1B)*(sam_mean-nu_1B)

      if(max(c(max(abs(nu_1B-nu_0B)),max(abs(var_1B-var_0B))))<0.000001)
      {
        break
      }
      var_0B<-var_1B
      nu_0B<-nu_1B
    }
    ####
    value_LRT_boot[r]<-max(prod((var_1A/var_10A)^(n/2)),prod((var_1B/var_10B)^(n/2)))
  }
  critical_value<-quantile(value_LRT_boot,error,names=F)
  cat("Test statistic value, critical value and P-value for the LRT statistic R_1 are respectively \n ",value_LRT, critical_value, mean(value_LRT>value_LRT_boot))
  if(value_LRT<critical_value)
  {
    cat("\n R_1 rejects the null hypothesis\n")
  }
  else
  {
    cat("\n R_1 do not reject the null hypothesis\n")
  }
  cat("*********************************************\n")
}

####
#Evaluate the observed test statistic values, critical values and P values for the test statistics W_2 and W_3 from a given data
#' Function to evaluate the observed, critical and P-values of the test statistics R_2 and R_3
#'
#' @param error Numeric, the level of significance/type-I error
#' @param a Numeric, the number of levels of the factor A
#' @param b Numeric, the number of levels of the factor B
#' @param cellsize matrix having a rows and b columns, and (i,j)-th element is the number of observations in the (i,j)-th cell
#' @param cellmean matrix having a rows and b columns, and (i,j)-th element is the sample mean of the observations in the (i,j)-th cell
#' @param cellvariance matrix having a rows and b columns, and (i,j)-th element is the sample variance of the observations in the (i,j)-th cell
#'
#'
#' @return numeric vector consisting of test statistic values, critical values and P-values of the test R_2 and R_3. Also, the decisions of the tests R_2 and R_3.
#' @export
#'
#' @examples
#' R_2_R_3(0.05,3,2, rbind(c(35,35),c(56,56),c(51,51)), rbind(c(52.7429,99.9143),c(58.8750,126.7679),c(95.8039,261.5490)), rbind(c(67.2555,274.0807),c(222.8023,955.0542),c(815.3608,8597.5325)))

R_2_R_3<-function(error,a,b,cellsize,cellmean,cellvariance){
  a<-a
  b<-b
  n<-cellsize
  sam_mean<-cellmean
  sam_var<-cellvariance
  error<-error
  fun_simul<-function(a,b,cellsize,cellmean,cellvariance){
    a<-a
    b<-b
    n<-cellsize
    sam_mean<-cellmean
    sam_var<-cellvariance
    T1<-array(NA,c((a-1),b))
    T2<-array(NA,c(a,(b-1)))
    for(i in 1:(a-1))
    {
      for(j in 1:b)
      {
        T1[i,j]<-(sam_mean[i+1,j]-sam_mean[i,j])/sqrt((sam_var[i+1,j]/n[i+1,j]+sam_var[i,j]/n[i,j]))
      }
    }
    for(i in 1:a)
    {
      for(j in 1:(b-1))
      {
        T2[i,j]<-(sam_mean[i,j+1]-sam_mean[i,j])/sqrt((sam_var[i,j+1]/n[i,j+1]+sam_var[i,j]/n[i,j]))
      }
    }
    value_1<-min(c(min(c(apply(T1,2,max))),min(c(apply(T2,1,max)))))
    value_2<-min(c(min(c(apply(T1,2,min))),min(c(apply(T2,1,min)))))
    return(c(value_1,value_2,t(sam_var)))
  }
  values<-fun_simul(a,b,cellsize,cellmean,cellvariance)
  value_simul_test1<-values[1]
  value_simul_test2<-values[2]
  ##Bootstrap step
  B<-10000 #number of bootstrap samples
  value_simul_test1_boot<-rep(NA,B)
  value_simul_test2_boot<-rep(NA,B)
  for(r in 1:B)
  {
    set.seed(17*r)
    sam_var_boot<-matrix(values[3:(a*b+2)],nrow=a,ncol=b,byrow=TRUE)
    sam_mean<-array(NA,dim=c(a,b))
    sam_var<-array(NA,dim=c(a,b))
    for(i in 1:a)
    {
      for(j in 1:b)
      {
        data<-rnorm(n[i,j],0,sqrt(sam_var_boot[i,j]))
        sam_mean[i,j]<-mean(data)
        sam_var[i,j]<-var(data)
      }
    }
    T1<-array(NA,c((a-1),b))
    T2<-array(NA,c(a,(b-1)))
    for(i in 1:(a-1))
    {
      for(j in 1:b)
      {
        T1[i,j]<-(sam_mean[i+1,j]-sam_mean[i,j])/sqrt((sam_var[i+1,j]/n[i+1,j]+sam_var[i,j]/n[i,j]))
      }
    }
    for(i in 1:a)
    {
      for(j in 1:(b-1))
      {
        T2[i,j]<-(sam_mean[i,j+1]-sam_mean[i,j])/sqrt((sam_var[i,j+1]/n[i,j+1]+sam_var[i,j]/n[i,j]))
      }
    }
    value_simul_test1_boot[r]<-min(c(min(c(apply(T1,2,max))),min(c(apply(T2,1,max)))))
    value_simul_test2_boot[r]<-min(c(min(c(apply(T1,2,min))),min(c(apply(T2,1,min)))))
  }
  critical_value_test1<-quantile(value_simul_test1_boot,1-error,names=F)
  critical_value_test2<-quantile(value_simul_test2_boot,1-error,names=F)

  cat("Test statistic value, critical value and P-value for the test statistic R_2 are respectively \n", value_simul_test1, critical_value_test1, mean(value_simul_test1<value_simul_test1_boot))
  if(value_simul_test1>critical_value_test1)
  {
    cat("\n R_2 rejects the null hypothesis")
  }
  else
  {
    cat("\n R_2 do not reject the null hypothesis")
  }
  cat("\n Test statistic value, critical value and P-value for the test statistic R_3 are respectively \n", value_simul_test2, critical_value_test2, mean(value_simul_test2<value_simul_test2_boot))
  if(value_simul_test2>critical_value_test2)
  {
    cat("\n R_3 reject the null hypothesis")
  }
  else
  {
    cat("\n R_3 do not reject the null hypothesis")
  }
  cat("\n*********************************************")
}



