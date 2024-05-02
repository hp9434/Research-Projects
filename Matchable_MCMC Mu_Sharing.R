# install.packages("truncnorm")
# install.packages("invgamma")
# install.packages("poisbinom", dependencies=T)
# install.packages("openxlsx",  dependencies=T)
# install.packages("xlsx", dependencies = T)

library(invgamma)
library(truncnorm)
library(poisbinom)

# Given parameters
MC=1
xu = 92.062*MC; xl = 92.024*MC; tx = 92; xspec = c(xl, tx, xu); mux = tx; sigx = 0.054*MC; 
yu = 92.052*MC; yl = 91.988*MC; ty = 92; yspec = c(yl, ty, yu); 
nc = 5; 
ST = 10000; 
MT=108; 
Dem = 10000
muy0 = 92*MC;
sigy0 = 0.01; 
a=3; b=0.015; #informative parameters (for true value)  

# M-H algorithm step values to be adjusted
step_sigy = 0.0005*MC; 
#step_muy = 0.0005*MC;

# Cylinder-liner: bin limits (given)
limx = c(92.024, 92.037, 92.050, 92.063, 92.076, 92.088)*MC;
# Cylinder-liner: area (in probability) for each bin
total_px = pnorm(xu, mux, sigx) - pnorm(xl, mux, sigx); 
px = array(0, nc); 
for (k in 1:nc) { 
  px[k] = (pnorm(limx[k+1], mux, sigx) - pnorm(limx[k], mux, sigx))/total_px };

  #Make sure that sum(px)=1
if(sum(px)>1){px<-(px-(sum(px)-1)/nc)}
sum(px)

#Pistons: bin limits (given)
limy = c(91.988, 92.001, 92.014, 92.027, 92.040, 92.052)*MC

# Generate data for Pistons Supplier
# For storing true mu and sig for all periods
mt_data_t_muy = array(0,MT);
mt_data_t_sigy = array(0,MT);
mt_data_t_pi = vector("list", length = MT);
mt_data_t_py = vector("list", length = MT);

# For storing posterior parameter values and densities FOR ALL PERIODS
mt_data_sigy = vector("list", length = MT);
mt_data_PB =  vector("list", length = MT);
mt_data_all_pi = vector("list", length = MT);
mt_data_post_dens = vector("list", length = MT);
mt_data_acpt = vector("list", length = MT);
mt_data_fcst_pi = vector("list", length = MT);

set.seed(1);
for (mt in 1:MT){
  ## Generating true process mean and variance  
  repeat {t_sigy = rinvgamma(1, a, b)*MC;
    if ((yu-yl)/(6*t_sigy)>= 1.33 && (yu-yl)/(6*t_sigy)<=2.00 ) {break} }; #check that the process is capable (i.e., cp>=1.33) and feasible (i.e., cp<=2.00)
  mt_data_t_sigy[mt] = t_sigy;
  if (mt == 1) { t_muy = ty; 
  } else{ t_muy = rtruncnorm(1, yl, yu, ty, t_sigy)} # muy
  mt_data_t_muy[mt] = t_muy;
  # Pistons: area (in probability) for each bin
  total_py = pnorm(yu, t_muy, t_sigy) - pnorm(yl, t_muy, t_sigy); 
  py = array(0, nc); 
  for (k in 1:nc) { 
      py[k] = (pnorm(limy[k+1], t_muy, t_sigy) - pnorm(limy[k], t_muy, t_sigy))/total_py};
  t_pi = rep(0,nc)
  for (j in 1:nc){
      t_pi[j]= min(py[j], px[j])
      }
  mt_data_t_py[[mt]] = py
  mt_data_t_pi[[mt]] = round(t_pi,3);
  }
#Check for capability and feasibility of the resulting process
cp.val=(yu-yl)/(6*mt_data_t_sigy)
cpk.val=array(0,MT)
for (j in 1:MT){
    cpk.val[j]<-min((mt_data_t_muy[j]-yl)/(3*mt_data_t_sigy[j]),(yu-mt_data_t_muy[j])/(3*mt_data_t_sigy[j]))
    }
sum(cpk.val>1.3) # check how many cpks are >1.3

#Check pi sum for each period
probsum=array(0,MT)
  for (mt in 1:MT) {
     probsum[mt]=sum(mt_data_t_pi[[mt]])
     }
sum(probsum<1)

# M-H Analysis 
set.seed(1)
for (mt in 1:MT){
  # For storing posterior parameter values and densities FOR EACH PERIOD
  data_sigy = array(0, ST); 
  data_pi = matrix(0, nc, ST);
  data_PB = array(0,ST);
  data_post_dens = array(0, ST); 
  data_acpt = array(0, ST); 
  data_while_check = array(0, ST); 
  
    if (mt==1){
     sigy = sigy0; # start with initial sigy value
     muy = muy0; # start with initial muy value
     a_0=2; b_0=0.01
     a=a_0; b=b_0; # non-informative prior
     } else {
             m_n=mean(mt_data_t_muy[1:(mt-1)]);
             #muy=(muy0 + (mt-1)*m_n)/(mt-1);
             a=a_0+(mt-1)/2;
             b=b_0 +(1/2)*sum((mt_data_t_muy[1:(mt-1)]-m_n)^2)+((mt-1)*(m_n-muy0)^2)/(2*(mt-1));
             sigy=(b/(a+1));
             }
     muy = mt_data_t_muy[mt]; # True muy value
     acpt_cnt = 0; # to check the acceptance rate
  
  # Prior and posterior densities for current (or "old") paramaters
  p_sigy = dinvgamma(sigy, a, b);
  # Pistons: area (in probability) for each bin
  total_pry = pnorm(yu, muy, sigy) - pnorm(yl, muy, sigy); 
  pry = array(0, nc); 
  for (k in 1:nc) { 
      pry[k] = (pnorm(limy[k+1], muy, sigy) - pnorm(limy[k], muy, sigy))/total_pry};
  t_api=array(0,nc)
  for (k in 1:nc){
      t_api[k]= min(pry[k], px[k])
      }
  p_PB = dpoisbinom(sum(t_api!=0), t_api);
  post_dens = p_PB*p_sigy;

  # M-H algorithm
  for (tt in 1:ST){

    # generate new sigy (i.e., sigy_new): should not be negative
    repeat {sigy_new = sigy + rnorm(1,0,step_sigy);
    if (sigy_new>=0) {break} };
    acpt = 0; # for checking whether the new parameter is accepted

    # Prior and posterior densities for "new" parameters
    p_sigy_new = dinvgamma(sigy_new, a, b);
        
    # Pistons: area (in probability) for each bin
    total_pry_new = pnorm(yu, muy, sigy_new) - pnorm(yl, muy, sigy_new); 
    pry_new = array(0, nc); 
    for (k in 1:nc) { 
      pry_new[k] = (pnorm(limy[k+1], muy, sigy_new) - pnorm(limy[k], muy, sigy_new))/total_pry_new};
    
    t_api_new=array(0,nc)
    for (k in 1:nc){
      t_api_new[k]= min(pry_new[k], px[k])
      }
    p_PB_new = dpoisbinom(sum(t_api_new!=0), t_api_new);
    post_dens_new = p_PB_new*p_sigy_new;
    alpha = min(1, post_dens_new/post_dens);

    # Update parameters / probabilities when accepted
    if (runif(1,0,1) <= alpha) {
      sigy = sigy_new; p_sigy = p_sigy_new; t_api = t_api_new;
      p_PB = p_PB_new;
      post_dens = post_dens_new;
      acpt = 1;
      acpt_cnt = acpt_cnt+1;
    };

    data_sigy[tt] = sigy;
    data_pi[, tt] = t_api;
    data_PB[tt] = p_PB
    data_post_dens[tt] = post_dens;
    data_acpt[tt] = acpt;
  }
  mt_data_sigy[[mt]] = data_sigy;
  mt_data_all_pi[[mt]] = data_pi;
  mt_data_PB[[mt]] = data_PB;
  mt_data_post_dens[[mt]] = data_post_dens;
  mt_data_acpt[[mt]] = data_acpt;
  mt_data_fcst_pi[[mt]] = rowMedians(data_pi);
}


#Find parameter predictions for each period
mt_fcst_sigy = vector("list", length = MT)
for (mt in 1:MT){
    mt_fcst_sigy[[mt]]=median(mt_data_sigy[[mt]])
    }

###Reactive approach
#Evaluate posterior piston proportion distribution per period
mt_fcst_pi_R = vector("list", length = MT)
for (mt in 1:MT){
    #Initial values for py
    total_py_new = pnorm(yu, mt_data_t_muy[mt], mt_fcst_sigy[[mt]]) - pnorm(yl, mt_data_t_muy[mt], mt_fcst_sigy[[mt]]); 
    py_new = array(0, nc); 
    for (k in 1:nc) { 
      py_new[k] = (pnorm(limy[k+1], mt_data_t_muy[mt], mt_fcst_sigy[[mt]]) - pnorm(limy[k], mt_data_t_muy[mt], mt_fcst_sigy[[mt]]))/total_py_new};
      #Make sure that sum(py_new)<=1
      if(sum(py_new)>1){py_new<-(py_new-(sum(py_new)-1)/nc)};
      #Evaluate predicted proportion matched in each bin     
      t_p_new1 = rep(0,nc)
      for (j in 1:nc){
          t_p_new1[j]= min(py_new[j], px[j])
      }
      mt_fcst_pi_R[[mt]] = round(t_p_new1,3)
      }

    #Pistons shortage and surplus amounts in each period for reactive forecasts
      Pi_ShortR=array(0,c(MT,nc))
      Pi_SrplsR=array(0,c(MT,nc))
      for (mt in 1:MT) {
              Pi_ShortR[mt,]= (mt_fcst_pi_R[[mt]]-mt_data_t_py[[mt]])
              Pi_SrplsR[mt,]= (mt_data_t_py[[mt]]-mt_fcst_pi_R[[mt]])
              for (k in 1:nc) {
                  Pi_ShortR[mt,k]=max(Pi_ShortR[mt,k],0)
                  Pi_SrplsR[mt,k]=max(Pi_SrplsR[mt,k],0)}              
           }
           

#Export raw reactive approach surplus/shortage results to Excel
require(xlsx)
require(openxlsx)
list_of_datasets <- list("Pistons_Short" = as.data.frame(Pi_ShortR), "Pistons_Surplus" = as.data.frame(Pi_SrplsR), "Actual_Del_Pistons" =as.data.frame(mt_data_t_py))
write.xlsx(list_of_datasets, file="Reactive_True_Mu_Shared.xlsx", append=T)


####Proactive approach####
#Objective function to minimize
        func <- function(mux_new){
        t_p=rep(0,nc)
        total_p = pnorm(xu, mux_new, sigx) - pnorm(xl, mux_new, sigx); 
        px_new = array(0, nc); 
        for (k in 1:nc) { 
            px_new[k] = (pnorm(limx[k+1], mux_new, sigx) - pnorm(limx[k], mux_new, sigx))/total_p};
        #Make sure that sum(px_new)=1
            if(sum(px_new)>1){px_new<-(px_new-(sum(px_new)-1)/nc)};
        for (k in 1:nc){
            t_p[k]= min(py_new[k], px_new[k])
            }
            return(sum(t_p))
            }

#Evaluate mismatch degree per period
mt_fcst_pi_P = vector("list", length = MT)
for (mt in 1:MT){
    #Initial values for py
    total_py_new = pnorm(yu, mt_data_t_muy[mt], mt_fcst_sigy[[mt]]) - pnorm(yl, mt_data_t_muy[mt], mt_fcst_sigy[[mt]]); 
    py_new = array(0, nc); 
    for (k in 1:nc) { 
      py_new[k] = (pnorm(limy[k+1], mt_data_t_muy[mt], mt_fcst_sigy[[mt]]) - pnorm(limy[k], mt_data_t_muy[mt], mt_fcst_sigy[[mt]]))/total_py_new};
    #Make sure that sum(py_new)<=1
    if(sum(py_new)>1){py_new<-(py_new-(sum(py_new)-1)/nc)};
 
    #Perform optimization and set current mux to optimal mux
     result.1<-optimize(f=func, interval=c(xl,xu), maximum=T)
     mux_curr=result.1$maximum
     
    #Current px based on optimized mux
     total_curr_x = pnorm(xu, mux_curr, sigx) - pnorm(xl, mux_curr, sigx); 
       px_curr = array(0, nc); 
        for (k in 1:nc) { 
            px_curr[k] = (pnorm(limx[k+1], mux_curr, sigx) - pnorm(limx[k], mux_curr, sigx))/total_curr_x};
      #Make sure that sum(px_new)=1
         if(sum(px_curr)>1){px_curr<-(px_curr-(sum(px_curr)-1)/nc)}
         
    #Matchable degree based on proactive forecasts per period 
    a_py=mt_data_t_py[[mt]]
    t_p_new=rep(0,nc)
    for (k in 1:nc){
            t_p_new[k]= min(px_curr[k], a_py[k])
    }
    mt_fcst_pi_P[[mt]]=round(t_p_new,3)
        }
 
 #Pistons shortage and surplus amounts in each period for proactive forecasts
      Pi_ShortP=array(0,c(MT,nc))
      Pi_SrplsP=array(0,c(MT,nc))
      for (mt in 1:MT) {
              Pi_ShortP[mt,]= (mt_fcst_pi_P[[mt]]-mt_data_t_py[[mt]])
              Pi_SrplsP[mt,]= (mt_data_t_py[[mt]]-mt_fcst_pi_P[[mt]])
              for (k in 1:nc) {
                  Pi_ShortP[mt,k]=max(Pi_ShortP[mt,k],0)
                  Pi_SrplsP[mt,k]=max(Pi_SrplsP[mt,k],0)}              
           }


#Export raw proactive approach surplus/shortage results to Excel
list_of_datasets <- list("Pistons_Short" = as.data.frame(round(Pi_ShortP,5)), "Pistons_Surplus" = as.data.frame(round(Pi_SrplsP,5)))
write.xlsx(list_of_datasets, file="Proactive_True_Mu_Shared.xlsx", append=T)

#Export raw forecasts results to Excel
list_of_datasets2 <- list("Reactive_Forecasts" = as.data.frame(mt_fcst_pi_R), "Proactive_Forecasts" = as.data.frame(mt_fcst_pi_P))
write.xlsx(list_of_datasets2, file="Forecasts_True_Mu_Shared.xlsx", append=T)