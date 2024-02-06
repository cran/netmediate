#'Moran_dv is a function used internally to calculate the amount of network autocorrelation in a co-evolution SAOM using the endogenous behavioral dependent variable and the network dependent variable
#
#'@network is the network object

Moran_dv<-function(network){

    #note to Scott--need to find thesims object i parent frame
 #return(print(names(parent.frame())))

  if(class(network)[1]=="igraph"){
    network<-intergraph::asNetwork(network)
  }
  thesims<-parent.frame()$sim_nets
  j<-parent.frame()$j

  return(sna::nacf(net=network,
                   y=thesims$sims[[j]][[1]][[2]][[1]],
                   type='moran',lag.max=1)[2])
}
