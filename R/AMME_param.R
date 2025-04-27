#'Function to calculate AMME
#'This function calls the appropriate algorithm for estimating the MEMS. The dependent variable is assumed to be a dyad or dyad equivalent (e.g., dyad-time-period, dyad-group)
#' not intended for external use
#'

AMME_param<-function(micro_model=micro_model,
                       macro_model=macro_model,
                       micro_process=micro_process,
                       mediator=mediator,
                       macro_function=macro_function,
                       link_id=link_id,
                       object_type=object_type,
                       controls=controls,
                       control_functions=control_functions,
                       interval=interval,
                       nsim=nsim,
                       algorithm=algorithm,
                       silent=silent,
                       full_output=full_output,
                       SAOM_data=SAOM_data,
                       SAOM_var=SAOM_var,
                      time_interval=time_interval,
                      covar_list=covar_list,
                       edgelist=edgelist,
                       net_logit_y=net_logit_y,
                       net_logit_x=net_logit_x,
                       group_id=group_id,
                       node_numbers=node_numbers,
                       sensitivity_ev=sensitivity_ev){

  if(class(macro_model)[1]%in%c("plm","lnam")){
    link_fun<-list(link="identity")
  }else{
    link_fun<-family(macro_model)
  }
  if(!link_fun$link%in%c("logit","probit","identity")){
    stop("netmediate currently only supports logit, probit, and linear GLMs for parametric AMME estimation.")
  }


  if(class(link_id)[1]!="list"){
    ID_vec<-link_id
    link_id<-list()
    for(i in 1:(length(controls)+1)){
      link_id[[i]]<-ID_vec
    }
  }

  if(length(link_id)<length(controls)){

    stop("Dimensions of link_id do not match number of controls specified.")

  }

  ####step 1: get MEMS output
  message("Getting output for micro-macro relationship.")
  MEMS_output<-MEMS(model=micro_model,
                    micro_process=micro_process,
                    macro_function=macro_function,
                    link_id=link_id,
                    object_type=object_type,
                    controls=controls,
                    control_functions=control_functions,
                    interval=interval,
                    nsim=nsim,
                    algorithm=algorithm,
                    silent=silent,
                    full_output=TRUE,
                    SAOM_data=SAOM_data,
                    SAOM_var = SAOM_var,
                    time_interval=time_interval,
                    covar_list=covar_list,
                    edgelist=edgelist,
                    net_logit_y=net_logit_y,
                    net_logit_x=net_logit_x,
                    group_id=group_id,
                    node_numbers=node_numbers)



  ##step 2 set up parameters to estimate AMME from macro model
  if(class(macro_model)[1]%in%c("glmerMod","lmerMod")){
    coef<-lme4::fixef(macro_model)

  }else{
    coef<-coef(macro_model)
  }

  if(class(macro_model)[1]=="lnam"){

    cov_mat<-macro_model$acvm[-c(nrow(macro_model$acvm)),-c(ncol(macro_model$acvm))]
    rownames(cov_mat)<-colnames(cov_mat)<-names(coef)

  }else{

    cov_mat<-stats::vcov(macro_model)

  }

  if(!mediator%in%names(coef)){
    stop("Mediator does not match model output. Check spelling and respecify.")
  }
  if(!any(is.null(controls)) & !any(controls%in%names(coef))){
    message("One more controls do not match model output. Check spelling and respecify.")
  }

  if(any(is.na(coef))|
     any(is.infinite(coef))){
    stop("Infinite or missing values in parameter estimates. Algorithm cannot continue.")

  }

  if(any(is.na(cov_mat))|
     any(is.infinite(cov_mat))){
    stop("Infinite or missing values in covariance matrix estimates. Algorithm cannot continue.")

  }
  interval<-sort(interval) #order from lowest to highest


  theta<-MASS::mvrnorm(n=nsim,
                       mu=coef,
                       Sigma=cov_mat,
                       empirical = TRUE)

  if(class(macro_model)[1]=="lnam"){

        #no intercept with LNAM, so DV can be ignored
      model_mat<-data.frame(macro_model$x)
      if(!is.null(macro_model$W1)){
        model_mat$rho1.1<-macro_model$W1[,,]%*%macro_model$y
      }
      if(!is.null(macro_model$W2)){
        model_mat$rho2.1<-macro_model$W2[,,]%*%macro_model$residuals
      }
      model_mat<-as.matrix(model_mat)

  }else{

      model_mat<-stats::model.matrix(macro_model)

  }

  if(class(macro_model)[1]%in%c("glmerMod","lmerMod")){
    model_mat<-model_mat[,!colnames(model_mat)%in%names(lme4::ranef(macro_model))]
  }



  ##step 3: start simulation

  ##now need to loop over parameter vector and predict
    #so, need to be able to match MEMS output to parameter iteration
  message("Beginning AMME estimation.")
  #get change stat values
  X_vals<-length(interval)
  M_vals<-length(interval)
  AMME_vec<-vector(length=nsim)# output data
  prop_vec<-vector(length=nsim) #proportion explained data
  k_mat <- list()  # Store k_vec results for RR computation
  crosswalk<-list() #only used for nested data
  for(i in 1:length(link_id)){
    crosswalk[[i]]<-data.frame(link_id=link_id[[i]])
  }

  for(j in 1:nsim){
    sim_model_mat<-model_mat

    #get predictions at each X by M interval
    for(i in 1:X_vals){

      k_vec<-matrix(NA,nrow=nrow(sim_model_mat),ncol=M_vals)

      for(m in 1:M_vals){ #loop over m intervals

        #use merges if observations don't align (i.e., nested data)
      if(nrow(sim_model_mat)!=nrow(MEMS_output$out_dat_main[[j]])){
          sim_dat<-data.frame(link_id[[1]]==rownames(MEMS_output$out_dat_main[[j]]),
                                 M=MEMS_output$out_dat_main[[j]][,m])
          sim_dat<-suppressMessages(plyr::join(sim_dat,crosswalk[[1]]))
          sim_model_mat[,mediator]<-sim_dat$M

      }else{
        sim_model_mat<-sim_model_mat[order(link_id[[1]]),]
        sim_model_mat[,mediator]<-MEMS_output$out_dat_main[[j]][,m]
      }

      if(!is.null(controls)){
        for(k in 1:length(MEMS_output$out_dat_controls)){
          control<-names(MEMS_output$out_dat_controls)[k]

          if(nrow(sim_model_mat)!=nrow(MEMS_output$out_dat_controls[[k]][[j]])){
            sim_dat<-data.frame(link_id[[k+1]]==rownames(MEMS_output$out_dat_main[[j]]),
                                M=MEMS_output$out_dat_main[[j]][,m])
            sim_dat<-suppressMessages(plyr::join(sim_dat,crosswalk[[k+1]]))
            sim_model_mat[,control]<-sim_dat$M
          }else{
            sim_model_mat<-sim_model_mat[order(link_id[[k+1]]),]
            sim_model_mat[,control]<-MEMS_output$out_dat_controls[[k]][[j]][,i] #for control k, simulation j, interval i
          }

        }
      }

      #get prediction
      lp<-as.matrix(sim_model_mat)%*%theta[j,]
      if(link_fun$link=="logit"){
        lp <- 1/(1 + exp(-lp))
      }
      if(link_fun$link=="probit"){
        lp<-VGAM::probitlink(lp,inverse=TRUE)
      }
      k_vec[,m]<-lp


    }#close M_vals loop

    #calculate differences in columns
    change_mat<-as.matrix(k_vec[,-c(1)])
    prop_mat<-change_mat
    k_mat[[j]] <- k_vec  # Store k_vec for later RR calculation

    for(k in 1:ncol(change_mat)){
      upper<-k+1
      lower<-k
      change_mat[,k]<-k_vec[,upper]-k_vec[,lower]
      prop_mat[,k]<-change_mat[,k]/k_vec[,upper]

    }

    #store all change stats
    if(i ==1){
      x_vec<-change_mat
    }else{
      x_vec<-cbind(x_vec,change_mat)
    }


  } #close x_vals loop

    #compute mean and store result
    AMME_vec[j]<-mean(x_vec,na.rm=T)
    prop_mat<-prop_mat[!is.infinite(prop_mat),]
    prop_vec[j]<-mean(prop_mat,na.rm=T)

    if(silent==FALSE){
      message("Simulation ", j, " of ",nsim," complete.")
    }
  }#close j loop



  ##package output
  summary_dat<-matrix(NA,nrow=2,ncol=5)
  rownames(summary_dat)<-c("AMME","Prop. Change in M")
  colnames(summary_dat)<-c("Estimate","Std. Dev.","lower 95% CI","Upper 95% CI","MC p-val")

  summary_dat[1,1]<-mean(AMME_vec,na.rm=TRUE)
  summary_dat[1,2]<-sd(AMME_vec,na.rm=TRUE)

  if(summary_dat[1,1]<0){
    summary_dat[1,5]<-length(AMME_vec[which(AMME_vec>=0)])/nsim
  }else{
    summary_dat[1,5]<-length(AMME_vec[which(AMME_vec<=0)])/nsim

  }
  summary_dat[1,3]<-quantile(AMME_vec,.025,na.rm=TRUE)
  summary_dat[1,4]<-quantile(AMME_vec,.975,na.rm=TRUE)
  summary_dat[2,1]<-mean(prop_vec,na.rm=TRUE)


  if(sensitivity_ev==TRUE){
    nsim <- length(k_mat)
    RR_vec <- numeric(nsim)
    n_boot <- 1000

    for (j in 1:nsim) {
      k_vec <- k_mat[[j]]
      if (is.null(k_vec)) next

      M_T1 <- mean(k_vec[, 2], na.rm = TRUE)  # E[M | T=1]
      M_T0 <- mean(k_vec[, 1], na.rm = TRUE)  # E[M | T=0]

      RR_vec[j] <- ifelse(M_T0 > 0, M_T1 / M_T0, NA)
    }

    RR_M_given_T <- mean(RR_vec, na.rm = TRUE)
    RR_sd <- sd(RR_vec, na.rm = TRUE)
    RR_CI <- quantile(RR_vec, c(0.025, 0.975), na.rm = TRUE)
    N1 <- sum(!is.na(RR_vec))

    boot_samples <- replicate(n_boot, {
      sample_RR <- sample(RR_vec, N1, replace = TRUE)
      mean(sample_RR, na.rm = TRUE)
    })
      CI_lower_boot <- quantile(boot_samples, 0.025, na.rm = TRUE)
      CI_upper_boot <- quantile(boot_samples, 0.975, na.rm = TRUE)

    # Compute E-value
      E_fun<-function(x){
        if (x < 1) {
          RR_inv <- 1 / x
          ev <- RR_inv + sqrt(RR_inv * (RR_inv - 1))
        } else {
          ev <- x + sqrt(x * (x - 1))
        }
        return(ev)
      }

      E_value<-E_fun(RR_M_given_T)
      CI_lower_boot<-E_fun(CI_lower_boot)
      CI_upper_boot<-E_fun(CI_lower_boot)


    message(sprintf("95%% CI (Bootstrap): [%.3f, %.3f]", CI_lower_boot, CI_upper_boot))
    message(sprintf("E-value: %.3f", E_value))

    # Return results
    ev_data <- list(
      RR_M_given_T = RR_M_given_T,
      CI_bootstrap = c(CI_lower_boot, CI_upper_boot),
      E_value = E_value
    )

  }


  if(full_output==FALSE){
    return(summary_dat)
  }else{
    out_dat<-list(summary_dat=summary_dat,
                  AMME_obs=AMME_vec,
                  prop_explained_obs=prop_vec,
                  ev_data=ev_data)
    return(out_dat)
  }


}#close function
