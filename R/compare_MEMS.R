#
# This function estimates the  difference in MEMS between two models
#
#
# partial_model is the model object without intervening variables
# full_model is the model with intervening variables.
#micro_process is the micro process of interest, provided as a character string. For SAOM, it should be as per the effectName parameter in your SAOM effects object. For example, effects_obj$effectName
#macro_function is the function to be calculated on the network object. NOTE: currently only supports statistics calculated on network and igraph objects.
#object type is the type of object used in the function to calculate macro statistics. Currently only supports igraph and network objects. If left NULL, the object is assumed to be a network object.
#the interval provides the values over which to calculate the MEMS. It defaults to 0 and 1
#nsim is the number of Monte Carlo samples
# silent tells R whether to provide updates on the algorithm's progress
#full_output tells R whether to return simulated distribution in addition to summary statistic


compare_MEMS <- function(partial_model,
                         full_model,
                         micro_process,
                         micro_process2=NULL,
                         macro_function,
                         macro_function2=NULL,
                         object_type=NULL,
                         interval=c(0,1),
                         nsim=500,
                         algorithm="parametric",
                         silent=FALSE,
                         full_output=FALSE,
                         sensitivity_ev = TRUE,
                         SAOM_data=NULL,
                         SAOM_var=NULL,
                         time_interval=NULL,
                         covar_list=NULL,
                         edgelist=NULL,
                         net_logit_y=NULL,
                         net_logit_x=NULL,
                         group_id=NULL,
                         node_numbers=NULL,
                         mediator=NULL,
                         link_id=NULL,
                         controls=NULL,
                         control_functions=NULL){

  diff_models<-FALSE
  if(!class(partial_model)[1]%in%class(full_model)){
    message("Distinct model types provided. Returning sensitivity test. Difference in MEMS should be interpreted as change in MEMS due to change in functional form specification.")
    if(is.null(micro_process2)){
      micro_process2<-micro_process
      message("micro_process2 not provided despite distinct model types. micro_process being recycled.")
    }
    diff_models<-TRUE
  }

  if(length(interval)>2){
    stop("compare_MEMS currently only allows interval specifications of length 2. Respecify and re-estimate.")
  }

  if(silent==FALSE){
    message("Computing MEMS for partial model.")
  }
  p_MEMS<-MEMS(partial_model,
               micro_process=micro_process,
               macro_function=macro_function,
               object_type=object_type,
               interval=interval,
               nsim=nsim,
               algorithm=algorithm,
               silent=silent,
               full_output=TRUE,
               sensitivity_ev = sensitivity_ev,
               SAOM_data=SAOM_data,
               SAOM_var=SAOM_var,
               time_interval=time_interval,
               covar_list=covar_list,
               edgelist=edgelist,
               net_logit_y=net_logit_y,
               net_logit_x=net_logit_x,
               group_id=group_id,
               node_numbers=node_numbers)

  if(silent==FALSE){
    message("Computing MEMS for full model.")
  }

  if(diff_models==TRUE){
    micro_process<-micro_process2
  }

  if(!is.null(macro_function2)){
    macro_function<-macro_function2
  }

  f_MEMS<-MEMS(full_model,
               micro_process=micro_process,
               macro_function=macro_function,
               object_type=object_type,
               interval=interval,
               nsim=nsim,
               algorithm=algorithm,
               silent=silent,
               full_output=TRUE,
               sensitivity_ev = sensitivity_ev,
               SAOM_data=SAOM_data,
               SAOM_var=SAOM_var,
               time_interval=time_interval,
               covar_list=covar_list,
               edgelist=edgelist,
               net_logit_y=net_logit_y,
               net_logit_x=net_logit_x,
               group_id=group_id,
               node_numbers=node_numbers)

  diff_MEMS_dat <- matrix(NA, nrow=3, ncol=7)
  colnames(diff_MEMS_dat) <- c(colnames(f_MEMS$summary_dat),
                               "Prop. of partial MEMS explained",
                               "E-value")

  r_names <- c("Partial MEMS", "Full MEMS", "Diff. in MEMS")
  if(diff_models==TRUE){
    r_names <- c(paste(class(partial_model)[1]),
                 paste(class(full_model)[1]),
                 "Diff. in MEMS")
  } else if(!is.null(micro_process2)) {
    r_names <- c(paste(micro_process), paste(micro_process2), "Diff. in MEMS")
  }
  rownames(diff_MEMS_dat) <- r_names

  diff_MEMS_dat[1, 1:5] <- p_MEMS$summary_dat[1, ]
  diff_MEMS_dat[1, 6] <- NA
  diff_MEMS_dat[1, 7] <- if (!is.null(p_MEMS$ev_data)) p_MEMS$ev_data$E_value else NA

  diff_MEMS_dat[2, 1:5] <- f_MEMS$summary_dat[1, ]
  diff_MEMS_dat[2, 6] <- NA
  diff_MEMS_dat[2, 7] <- if (!is.null(f_MEMS$ev_data)) f_MEMS$ev_data$E_value else NA

  diff_MEMS <- p_MEMS$mems_samples - f_MEMS$mems_samples
  diff_MEMS_dat[3, 1] <- mean(diff_MEMS)
  diff_MEMS_dat[3, 2] <- sd(diff_MEMS)
  diff_MEMS_dat[3, 3] <- quantile(diff_MEMS, 0.025, na.rm = TRUE)
  diff_MEMS_dat[3, 4] <- quantile(diff_MEMS, 0.975, na.rm = TRUE)
  diff_MEMS_dat[3, 5] <- if (diff_MEMS_dat[3, 1] < 0) {
    length(diff_MEMS[diff_MEMS > 0]) / nsim
  } else {
    length(diff_MEMS[diff_MEMS < 0]) / nsim
  }
  diff_MEMS_dat[3, 6] <- mean(diff_MEMS) / p_MEMS$summary_dat[1, 1]

  if (sensitivity_ev) {
    m_part <- mean(f_MEMS$output_data[,2], na.rm = TRUE)
    m_full <- mean(p_MEMS$output_data[,2], na.rm = TRUE)
    RR_raw <- m_part / m_full

    if (is.finite(RR_raw) && RR_raw > 0) {
      RR <- ifelse(RR_raw < 1, 1 / RR_raw, RR_raw)
      e_val_delta <- RR + sqrt(RR * (RR - 1))
      if (e_val_delta < 1) e_val_delta <- NA
    } else {
      e_val_delta <- NA
    }

    if(!class(partial_model)[1]%in%class(full_model)){
      e_val_delta<-NA
    }
    diff_MEMS_dat[3, 7] <- e_val_delta
  } else {
    diff_MEMS_dat[3, 7] <- NA
  }

  if (full_output == FALSE) {
    return(diff_MEMS_dat)
  } else {
    diff_results <- list(
      diff_summary = diff_MEMS_dat,
      diff_samples = diff_MEMS
    )
    all_results <- list(
      diff_MEMS_results = diff_results,
      p_MEMS_results = p_MEMS,
      f_MEMS_results = f_MEMS
    )
    return(all_results)
  }
}
