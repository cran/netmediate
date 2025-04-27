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
#model_comparison tells R whether to test the robustness of direct, total, and indirect MEMS estimates to choice of model
#partial_model2 is an optional parameter used only when model_comparison is TRUE. It should be the alternative model to compare partial MEMS estimates
#full_model 2 is an optional parameter used only when model_comparison is TRUE.  It should be a alternative model  to compare the full MEMS estimates
#micro_process2 is an optional parameter used only when model_comparison is set to TRUE. It should be the character string denoting the micro process in full_model2 and partial_model2

mediate_MEMS <- function(partial_model,
                         full_model,
                         micro_process,
                         macro_function,
                         model_comparison = FALSE,
                         partial_model2 = NULL,
                         full_model2 = NULL,
                         micro_process2 = NULL,
                         macro_function2 = NULL,
                         object_type = NULL,
                         interval = c(0, 1),
                         nsim = 500,
                         algorithm = "parametric",
                         silent = FALSE,
                         full_output = FALSE,
                         sensitivity_ev = TRUE,  # <-- NEW argument added
                         SAOM_data = NULL,
                         SAOM_var = NULL,
                         time_interval = NULL,
                         covar_list = NULL,
                         edgelist = NULL,
                         net_logit_y = NULL,
                         net_logit_x = NULL,
                         group_id = NULL,
                         node_numbers = NULL,
                         mediator = NULL,
                         link_id = NULL,
                         controls = NULL,
                         control_functions = NULL) {

  if (model_comparison == FALSE) {
    result <- compare_MEMS(partial_model = partial_model,
                           full_model = full_model,
                           micro_process = micro_process,
                           micro_process2 = micro_process2,
                           macro_function = macro_function,
                           macro_function2 = macro_function2,
                           object_type = object_type,
                           interval = interval,
                           nsim = nsim,
                           algorithm = algorithm,
                           silent = silent,
                           full_output = full_output,
                           sensitivity_ev = sensitivity_ev,  # <-- Pass new argument
                           SAOM_data = SAOM_data,
                           SAOM_var = SAOM_var,
                           time_interval = time_interval,
                           covar_list = covar_list,
                           edgelist = edgelist,
                           net_logit_y = net_logit_y,
                           net_logit_x = net_logit_x,
                           group_id = group_id,
                           node_numbers = node_numbers,
                           mediator = mediator,
                           link_id = link_id,
                           controls = controls,
                           control_functions = control_functions)
    return(result)
  }

  return_samples <- FALSE
  if (full_output == TRUE) {
    return_samples <- TRUE
  }

  if (isTRUE(model_comparison)) {
    if (is.null(partial_model2) | is.null(full_model2) | is.null(micro_process2)) {
      stop("partial_model2, full_model2, and micro_process2 must be provided when model_comparison is set to TRUE")
    }
  }

  message("Estimating mediation results for model set 1")
  model_set1 <- compare_MEMS(partial_model = partial_model,
                             full_model = full_model,
                             micro_process = micro_process,
                             macro_function = macro_function,
                             macro_function2 = NULL,
                             object_type = object_type,
                             interval = interval,
                             nsim = nsim,
                             algorithm = algorithm,
                             silent = silent,
                             full_output = TRUE,
                             sensitivity_ev = sensitivity_ev,  # <-- Pass new argument
                             SAOM_data = SAOM_data,
                             SAOM_var = SAOM_var,
                             time_interval = time_interval,
                             covar_list = covar_list,
                             edgelist = edgelist,
                             net_logit_y = net_logit_y,
                             net_logit_x = net_logit_x,
                             group_id = group_id,
                             node_numbers = node_numbers,
                             mediator = mediator,
                             link_id = link_id,
                             controls = controls,
                             control_functions = control_functions)

  message("Estimating mediation results for model set 2")
  if (!is.null(macro_function2)) {
    macro_function <- macro_function2
  }
  model_set2 <- compare_MEMS(partial_model = partial_model2,
                             full_model = full_model2,
                             micro_process = micro_process2,
                             macro_function = macro_function,
                             macro_function2 = NULL,
                             object_type = object_type,
                             interval = interval,
                             nsim = nsim,
                             algorithm = algorithm,
                             silent = silent,
                             full_output = TRUE,
                             sensitivity_ev = sensitivity_ev,  # <-- Pass new argument
                             SAOM_data = SAOM_data,
                             SAOM_var = SAOM_var,
                             time_interval = time_interval,
                             covar_list = covar_list,
                             edgelist = edgelist,
                             net_logit_y = net_logit_y,
                             net_logit_x = net_logit_x,
                             group_id = group_id,
                             node_numbers = node_numbers,
                             mediator = mediator,
                             link_id = link_id,
                             controls = controls,
                             control_functions = control_functions)

  sens_test_dat <- matrix(NA, nrow = 3, ncol = 5)
  colnames(sens_test_dat) <- c("Change in estimate between models", "Std. Dev.", "Low 95% CI", "Upp. 95% CI", "MC p-val")
  rownames(sens_test_dat) <- c("Partial MEMS", "Full MEMS", "Indirect MEMS")

  # Estimates
  sens_test_dat[, 1] <- model_set1$diff_MEMS_results$diff_summary[, 1] - model_set2$diff_MEMS_results$diff_summary[, 1]
  sens_test_dat[, 2] <- c(sd(model_set1$p_MEMS_results$mems_samples - model_set2$p_MEMS_results$mems_samples, na.rm = TRUE),
                          sd(model_set1$f_MEMS_results$mems_samples - model_set2$f_MEMS_results$mems_samples, na.rm = TRUE),
                          sd(model_set1$diff_MEMS_results$diff_samples - model_set2$diff_MEMS_results$diff_samples, na.rm = TRUE))

  sens_test_dat[, 3] <- c(quantile(model_set1$p_MEMS_results$mems_samples - model_set2$p_MEMS_results$mems_samples, .025, na.rm = TRUE),
                          quantile(model_set1$f_MEMS_results$mems_samples - model_set2$f_MEMS_results$mems_samples, .025, na.rm = TRUE),
                          quantile(model_set1$diff_MEMS_results$diff_samples - model_set2$diff_MEMS_results$diff_samples, .025, na.rm = TRUE))

  sens_test_dat[, 4] <- c(quantile(model_set1$p_MEMS_results$mems_samples - model_set2$p_MEMS_results$mems_samples, .975, na.rm = TRUE),
                          quantile(model_set1$f_MEMS_results$mems_samples - model_set2$f_MEMS_results$mems_samples, .975, na.rm = TRUE),
                          quantile(model_set1$diff_MEMS_results$diff_samples - model_set2$diff_MEMS_results$diff_samples, .975, na.rm = TRUE))

  d1 <- model_set1$p_MEMS_results$mems_samples - model_set2$p_MEMS_results$mems_samples
  sens_test_dat[1, 5] <- if (sens_test_dat[1, 1] < 0) {
    length(d1[d1 > 0]) / nsim
  } else {
    length(d1[d1 < 0]) / nsim
  }

  d2 <- model_set1$f_MEMS_results$mems_samples - model_set2$f_MEMS_results$mems_samples
  sens_test_dat[2, 5] <- if (sens_test_dat[2, 1] < 0) {
    length(d2[d2 > 0]) / nsim
  } else {
    length(d2[d2 < 0]) / nsim
  }

  d3 <- model_set1$diff_MEMS_results$diff_samples - model_set2$diff_MEMS_results$diff_samples
  sens_test_dat[3, 5] <- if (sens_test_dat[3, 1] < 0) {
    length(d3[d3 > 0]) / nsim
  } else {
    length(d3[d3 < 0]) / nsim
  }

  results <- list(model_set1 = model_set1$diff_MEMS_results$diff_summary,
                  model_set2 = model_set2$diff_MEMS_results$diff_summary,
                  model_comparison = sens_test_dat)

  if (return_samples == TRUE) {
    sample_dat <- list(model_set1 = model_set1,
                       model_set2 = model_set2,
                       sens_test = list(diff_MEMS_results = d3,
                                        p_MEMS_results = d1,
                                        f_MEMS_results = d2))
    results <- list(summary_results = results,
                    sample_results = sample_dat)
  }

  return(results)
}
