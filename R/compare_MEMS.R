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



compare_MEMS<-function(partial_model,
                       full_model,
                       micro_process,
                       macro_function,
                       object_type=NULL,
                       interval=c(0,1),
                       nsim=500,
                       algorithm="parametric",
                       silent=FALSE,
                       full_output=FALSE,
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


  if(!class(partial_model)[1]%in%class(full_model)){
    message("Distinct model types provided. Returning sensitivity test. Difference in MEMS should be interpreted as change in MEMS due to change in functional form specification.")
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

  f_MEMS<-MEMS(full_model,
               micro_process=micro_process,
               macro_function=macro_function,
               object_type=object_type,
               interval=interval,
               nsim=nsim,
               algorithm=algorithm,
               silent=silent,
               full_output=TRUE,
               SAOM_data=SAOM_data,
               SAOM_var=SAOM_var,
               time_interval=time_interval,
               covar_list=covar_list,
               edgelist=edgelist,
               net_logit_y=net_logit_y,
               net_logit_x=net_logit_x,
               group_id=group_id,
               node_numbers=node_numbers)


  diff_MEMS_dat<-matrix(NA,nrow=3,ncol=6)
  colnames(diff_MEMS_dat)<-c(colnames(f_MEMS$summary_dat),"Prop. of partial MEMS explained")
  rownames(diff_MEMS_dat)<-c("Partial MEMS","Full MEMS","Diff. in MEMS")
  diff_MEMS_dat[1,]<-c(p_MEMS$summary_dat[1,],NA)
  diff_MEMS_dat[2,]<-c(f_MEMS$summary_dat[1,],NA)

  diff_MEMS<-p_MEMS$mems_samples-f_MEMS$mems_samples
  diff_MEMS_dat[3,1]<-mean(diff_MEMS)
  diff_MEMS_dat[3,2]<-sd(diff_MEMS)
  diff_MEMS_dat[3,3]<-quantile(diff_MEMS,.025,na.rm=TRUE)
  diff_MEMS_dat[3,4]<-quantile(diff_MEMS,.975,na.rm=TRUE)
  if(diff_MEMS_dat[3,1]<0){
    diff_MEMS_dat[3,5]<-length(diff_MEMS[which(diff_MEMS>0)])/(nsim)
  }else{
    diff_MEMS_dat[3,5]<-length(diff_MEMS[which(diff_MEMS<0)])/(nsim)

  }
  diff_MEMS_dat[3,6]<-mean(diff_MEMS)/p_MEMS$summary_dat[1,1]

 if(full_output==FALSE){
   return(diff_MEMS_dat)
 }else{

   diff_results<-list(diff_summary=diff_MEMS_dat,diff_samples=diff_MEMS)

   all_results<-list(diff_MEMS_results=diff_results,
                     p_MEMS_results=p_MEMS,
                     f_MEMS_results=f_MEMS)
   return(all_results)

 }


}
