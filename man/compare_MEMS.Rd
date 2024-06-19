\name{compare_MEMS}
\alias{compare_MEMS}
\title{
Function to compare micro effect on macro structure (MEMS) estimates between models.}
\description{
\code{compare_MEMS} implements parametric and nonparametric routines to compare \code{MEMS} estimate between models. When compared between nested models, \code{compare_MEMS} results can be interpreted as the portion of a MEMS explained by a mediating or confounding variable. When compared between models with distinct functional forms and the same specification, \code{compare_MEMS} results can be interpreted as the sensitivity of \code{MEMS} results to decision about model functional form.

The difference in MEMS is the change in MEMS after one or more micro-processes are included into a model or, in the case of sensitivity tests, when the functional form is changed. Let \eqn{MEMS_p} represent the MEMS obtained from a model that omits one or more intervening variables and \eqn{MEMS_f} be the MMES obtained from a model that includes the intervening variable(s). The change in MEMS is given \deqn{\Delta MEMS=MEMS_p-MEMS_f}. \eqn{MEMS_p} and \eqn{MEMS_f} may also be have the same specification but use distinct functional forms or other modeling decisions in the case of sensitivity tests. Tuning parameters can be assigned to toggle the strength of \eqn{\theta} in model-implied estimates of \eqn{MEMS}. \code{MEMS} currently accepts \code{glm}, \code{glmer}, \code{ergm}, \code{btergm}, \code{sienaFit}, \code{rem.dyad}, and \code{netlogit} objects and implements both parametric and nonparametric estimation. Pooled estimation for multiple network models is also implemented for \code{ergm} and \code{sienaFit} objects.

}
\usage{
compare_MEMS(partial_model,
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
      control_functions=NULL)}

\arguments{
  \item{partial_model}{the micro-model excluding one or more intervening or confounding variables of interest. May also be a fully specified model with a distinct functional form in the case of sensitivity tests. Currently accepts \code{glm}, \code{glmer}, \code{ergm}, \code{btergm}, \code{sienaFit}, \code{rem.dyad}, and \code{netlogit} objects. Pooled estimation for multiple network models is also implemented for \code{ergm} and \code{sienaFit} objects. To implement pooled estimation, \code{model} should be provided as a list of \code{ergm} or \code{sienaFit} objects.
}
  \item{full_model}{the micro-model including one or more intervening or confounding variables of interest. May also be a fully specified model with a distinct functional form in the case of sensitivity tests. Currently accepts \code{glm}, \code{glmer}, \code{ergm}, \code{btergm}, \code{sienaFit}, \code{rem.dyad}, and \code{netlogit} objects. Pooled estimation for multiple network models is also implemented for \code{ergm} and \code{sienaFit} objects. To implement pooled estimation, \code{model} should be provided as a list of \code{ergm} or \code{sienaFit} objects.
}
\item{micro_process}{a character string containing the name of the micro process of interest. The character string should exactly match coefficient names in \code{model} output.}
\item{macro_function}{a \code{function} that calculates the macro statistic of interest. Currently accepts user defined functions as well as functions inherent in the \code{\link[igraph]{igraph}} and \code{statnet} packages for \R.}
\item{object_type}{A character string that tells netmediate the type of object to apply the
\code{macro_function} to. Currently accepts \code{\link[igraph]{igraph}} and \code{\link[network]{network}} objects. If left \code{NULL}, \code{\link[network]{network}} objects are assumed. Can be over-ridden to use other object types with a user-function by defining a function that accepts either a \code{\link[network]{network}} or \code{\link[igraph]{igraph}} object and returns a numeric value or vector of numeric values (see examples).}
\item{interval}{The value of tuning parameters to assign to \eqn{\theta}. Should be provided as a vector of numeric values with 2 entries.}
\item{nsim}{The number of simulations or bootstrap samples to use during estimation.}
\item{algorithm}{The estimation algorithm to be used. Currently accepts \code{"parametric"} and \code{"nonparametric"}. If \code{"parametric"}, estimation is obtained with Monte Carlo sampling. If \code{"nonparametric"}, estimation uses bootstrap resampling.}
\item{silent}{logical parameter. Whether to provide updates on the progress of the simulation or not.}
\item{full_output}{logical parameter. If set to \code{TRUE}, \code{compare_MEMS} will return all sampled statistics and complete results for \eqn{MEMS_p} and \eqn{MEMS_f}.}
\item{SAOM_data}{required when the model is a \code{\link[RSiena]{sienaFit}} object; ignored otherwise. If a \code{\link[RSiena]{sienaFit}} object is provided, \code{SAOM_data} should be the \code{siena} object that contains the data for SAOM estimation. If using pooled estimation on multiple \code{\link[RSiena]{sienaFit}} objects (i.e., providing a list of \code{\link[RSiena]{sienaFit}} objects), then \code{SAOM_data} should be provided as an ordered list with each entry containing the \code{\link[RSiena]{siena}} object corresponding to list of \code{\link[RSiena]{sienaFit}} objects.}
\item{SAOM_var}{optional parameter when the model is a \code{\link[RSiena]{sienaFit}} object. \code{SAOM_var} is a list of of the \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}} objects used to assign time varying node and dyad covariates when calling \code{\link[RSiena]{sienaDataCreate}}. If provided, \code{netmediate} assigns the varying node covariates and dyad covariates to each simulated network. This parameter is required when \code{macro_function} computes a statistic that varies as a function of time varying node or dyad covariates (i.e., network segregation, assorativity). Time invariant characteristics (\code{\link[RSiena]{coCovar}} and \code{\link[RSiena]{coDyadCovar}}) are handled internally by \code{MEMS} and should not be provided. When providing a list of \code{\link[RSiena]{sienaFit}} objects for pooled estimation, \code{SAOM_var} should be provided as a list of lists, where each entry in the list contains a list of \code{\link[RSiena]{varCovar}} and \code{\link[RSiena]{varDyadCovar}} objects associated with corresponding \code{\link[RSiena]{sienaFit}} object. }
\item{time_interval}{an optional parameter to be used with \code{\link[relevent]{rem.dyad}} objects. May be provided as a numeric vector or the character string \code{"aggregate"}. If a numeric vector is provided unique network snapshots at each interval. For example, \code{time_interval=c(0,2,3)} would induce two networks, one for the 0 - 2 time period and one for the 2 - 3 time period. If specified as \code{"aggregate"}, the MEMS is calculated by creating an aggregated cross-sectional representation of the entire event sequence. If left \code{NULL}, defaults to |\code{"aggregate"}.}
\item{covar_list}{an optional list of sender/receiver covariates used in \code{\link[relevent]{rem.dyad}} estimation. Only required for \code{\link[relevent]{rem.dyad}} objects when covariates are included. The list format should correspond to the format required by \code{\link[relevent]{rem.dyad}}}.
\item{edgelist}{an optional three column edgelist providing the sender, receiver, and time of event occurrence when using rem.\code{\link[relevent]{rem.dyad}}. Only required when \code{time_interval} is set to \code{NULL} or \code{"aggregate"}. Ignored for other types of models.}
\item{net_logit_y}{the dependent variable for \code{\link[sna]{netlogit}} objects. Should be provided as a vector. Only required when model is a \code{\link[sna]{netlogit}} object.}
\item{net_logit_x}{the matrix of independent variables for \code{\link[sna]{netlogit}} type objects. Only required when model is a \code{\link[sna]{netlogit}} object.}
\item{group_id}{optional vector of group identifiers to use when estimating a \code{glm} or \code{\link[lme4]{glmer}} on grouped data (i.e., multiple time periods, multiple networks). When specified,  \code{MEMS} will induce unique networks for each grouping factor. If left unspecified, all groups/time periods are pooled. If using \code{\link[lme4]{glmer}}, the grouping factor does not have to be provided as part of the model or used as a random effect.}
\item{node_numbers}{a numeric vector containing the number of nodes in each group_id when using \code{glm} or \code{\link[lme4]{glmer}}. If estimating MEMS aggregated over all networks (i.e., \code{group_id=NULL}), this shoud be the total number of nodes in all networks. Required when using \code{glm} or \code{\link[lme4]{glmer}}, ignored otherwise.}
\item{mediator}{a character string detailing the mediator of interest. Intended for internal use with the \code{\link{AMME}} function; not intended for end users.}
\item{link_id}{a vector or list of vectors corresponding to unique identifiers. Intended for internal use with the \code{\link{AMME}} function; not intended for end users. }
\item{controls}{a vector of character strings listing the controls to be calculated when using \code{\link{AMME}}. Intended for internal use with the \code{\link{AMME}} function; not intended for end users. }
\item{control_functions}{a list of functions to calculate the macro control variables provided in controls. Intended for internal use with the \code{\link{AMME}} function; not intended for end users. }

}
\details{
Compares MEMS estimates between two models. If one or more confounding or intervening variables are excluded or included between models, the change in MEMS can be interpreted as the portion of the MEMS explained by one or more confounding or intervening variable. If two models are provided with the same specification but a distinct functional form, the change in MEMS is a sensitivty test of how much the MEMS estimate changes because of a model decision. This can be useful, for example, when comparing TERGM and SAOM estimates as each models make distinct assumptions about sources of network change and the temporal ordering of tie changes.

\code{compare_MEMS} functionality inherits directly from the \code{MEMS} command. See the \code{MEMS} page for more details.



}
\value{
If \code{full_output=FALSE}, then a table is returned with the change MEMS, its standard error, confidence interval, and p-value, and the same results for the partial and complete MEMS.

If \code{full_output=TRUE}, then a list is returned with the following three elements.
\item{diff_MEMS_results}{is the table of summary output containing the MEMS, its standard error, confidence interval, and p-value, and a list of the simulated values of the change in MEMS.}
\item{p_MEMS_results}{contains the summary statistics for the partial MEMS along with all simulated statistics.}
\item{f_MEMS_results}{contains the summary statistics for the full MEMS along with all simulated statistics.}
}
\references{
Duxbury, Scott W. 2024. "Micro Effects on Macro Structure in Social Networks." Sociological Methodology.

Wertsching, Jenna, and Scott W. Duxbury. Working paper. "Comparing Micro Effects on Macro Structure between Nested Models."
}
\author{
Duxbury, Scott W. Associate Professor, University of North Carolina--Chapel Hill, Department of Sociology.
}

\seealso{
\code{\link{AMME}}
\code{\link{MEMS}}
\code{\link[ergMargins]{ergm.mma}}
\code{\link[mediation]{mediate}}
}
\examples{
\dontshow{
require(ergm)
require(network)
require(sna)


set.seed(21093)
a1<-network::as.network(matrix(c(rbinom(10, 1,.3),
            rbinom(10, 1,.3),
             rbinom(10, 1,.3),
              rbinom(10, 1,.3),
               rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3),
            rbinom(10, 1,.3)),
          nrow=10,ncol=10))

network::set.vertex.attribute(a1,"var.1",rbinom(10,1,.3))
a<-ergm(a1~edges+nodeifactor("var.1")+nodeofactor("var.1"))

compare_MEMS(partial_model=a,
            full_model=a,
      micro_process="nodeifactor.var.1.1",
      macro_function=gtrans,
      nsim=20,
      silent=TRUE)

}

\donttest{

##############
# Not run
###############
library(statnet)
library(igraph)
data("faux.mesa.high")

  #how much of the effect of racial homophily on transitivity
    #is explained by triadic closure effects?

model<-ergm(faux.mesa.high~edges+nodecov("Grade")+nodefactor("Race")+
               nodefactor("Sex")+nodematch("Race")+nodematch("Sex")+absdiff("Grade"))

model2<-ergm(faux.mesa.high~edges+nodecov("Grade")+nodefactor("Race")+
               nodefactor("Sex")+nodematch("Race")+nodematch("Sex")+absdiff("Grade")+
               gwesp(.5,fixed=TRUE))


compare_MEMS(partial_model=model,
              full_model=model2,
              micro_process="nodematch.Race",
             macro_function=transitivity,
             object_type = "igraph",
             silent=FALSE,
             algorithm="parametric")


}

}
\keyword{ ~mediation }
\keyword{ ~networks }
\keyword{ ~macro }
\keyword{ ~micro }
