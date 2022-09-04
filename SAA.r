#'@importFrom cluster pam
#'@export
SAA<-function (strata, errors, suggestions = NULL,
                         Temp=NA,initialStrata, decrement_constant=0.95, end_time =140,
                         showSettings = FALSE, jsize=100,length_of_markov_chain =50,
                         verbose = FALSE, dominio=NULL,minnumstrat=NULL,kmax_percent=0.025,ProbNewStratum=NA,
                         strcens=FALSE,writeFiles=FALSE, showPlot=TRUE, minTemp = 0.0005, realAllocation=TRUE)
{



  stringMin <- rep(1, nrow(strata))
  if(is.na(initialStrata)){
    stringMax <- (rep(ceiling(sqrt(nrow(strata)) *1), nrow(strata)))}else{  stringMax <- (rep(ceiling(initialStrata), nrow(strata)))}


  vars = length(stringMin)

  if (verbose)
    cat("Testing the sanity of parameters...\n")
  if (length(stringMin) != length(stringMax)) {
    stop("The vectors stringMin and stringMax must be of equal length.")
  }

  if (jsize < 1) {
    stop("The number of Markov Chains must be at least 1.")
  }


  if (showSettings) {
    if (verbose)
      cat("The start conditions:\n")
    result = list(stringMin = stringMin, stringMax = stringMax,
                  suggestions = suggestions, solSize = 1, jsize = jsize)
    class(result) = "rbga"
    cat(summary(result))
  }else {
    if (verbose)
      cat("Not showing GA settings...\n")
  }

  start_time <- proc.time()

  if (vars > 0) {
    nvar <- length(grep("CV", names(errors)))
    strata_means<-as.matrix(strata[,3:(nvar+2)])

    ndom <- nrow(errors)
    vars = length(stringMin)
    solution<-NULL

    if (!is.null(suggestions)){
      solution<-suggestions$suggestion
    }else {
      solution=sample.int(stringMax,vars,replace=TRUE,prob=NULL)
      #    solution<-cluster::pam(strata[,3:(nvar+2)],ceiling(sqrt(nrow(strata))),pamonce = 3)$cluster
    }

    kmax=ceiling(nrow(strata)*kmax_percent)




    if(is.na(length_of_markov_chain)){L_size<-2*(nrow(strata))}else{L_size<-length_of_markov_chain}

    soluz <- NULL
    v <- NULL
    dimens <- NULL
    censiti <- 0
    solution<-floor(solution)
    if( "matrix"%in%class(solution)){solution<-solution[1,]}
    strcor <- SimAnn::aggrStrata(strata, nvar, solution, censiti,
                                 dominio=dominio)
    dimsamp <- nrow(strcor)
    if (strcens == TRUE)
      strcor <- rbind(strcor, cens)
    dimens <- nrow(strcor)
    # alfa<-c(rep(1/nvar, nvar))
    res<-SimAnn::bethel_alfa(strcor, errors,realAllocation = TRUE)
    soluz <- res[[1]]
    alfa<- res[[2]]

    tot <- sum(soluz)

    best_tot<-tot
    best_sol<-solution
    iter<-0
    x_axis <- iter     # x axis
    y_axis <- tot  # y axis
    deltaList<-0
    j<-1
    k<-0

    L_chain<-tot
    Prev_tot<-tot
    newTemp<-Temp

    current_time<-proc.time() - start_time



    while(j < (jsize+1) & current_time[3] < end_time & newTemp > minTemp){

      #
      #       if (j ==1){(k=kmax)}else{k=k-1}
      #       if(k==0){k=1}
      #       if ( min(L_chain) >= Prev_tot){k=k+1}else{k=k-1}
      #       if(k >kmax){k=kmax}
      Prev_tot<-min(L_chain)
      #  if (j ==1){(k=kmax)}else{k=k-1}

      if (j ==1){(k<-kmax)}else{k=1}
      if (runif(1)<=1/(length_of_markov_chain)){
        groups<-unique(solution)
        newgroup<-max(groups)+1
        for (i in 1:length(solution)){
          if (runif(1)<=ProbNewStratum){
            solution[i]<-newgroup
          }}
      }

      for(i in 1:L_size){


        iter<-iter+1

        newsolution<-solution

        groups<-unique(newsolution)
        grps<-groups[sample(length(groups),2,replace=FALSE)]
        orig_group<-grps[1]

        orig_group_records<- strata_means[which(newsolution %in% orig_group),]

        if (is.null(nrow(orig_group_records))){
          num_rec_orig<-length(orig_group_records)
        } else{
          num_rec_orig<-nrow(orig_group_records)
        }

        orig_loc <- sample(num_rec_orig,k,replace=TRUE)
        #  cat("k is now: ",k,"\n")
        # cat("Sample is now: ",tot,"\n")


        replace_group<-grps[2]
        newsolution[which(newsolution %in% orig_group)[orig_loc]]<-replace_group

        str<-NULL
        str <- cbind(strata, newsolution)
        strataReplace <- str[which(str$newsolution==replace_group),]
        strataOrig <-str[which(str$newsolution==orig_group),]
        strataDelta<-rbind(strataReplace,strataOrig)


        strcorDelta <- SimAnn::aggrStrata(strataDelta, nvar=nvar,strataDelta$newsolution, censiti=censiti,
                                          dominio=dominio)
        newstrcor<-NULL
        newstrcor<-strcor[-which(strcor$STRATO %in% c(orig_group,replace_group)),]
        newstrcor<-rbind(newstrcor,strcorDelta)
        newstrcor<-newstrcor[order(newstrcor$STRATO),]
        dimsamp <- nrow(newstrcor)
        if (strcens == TRUE)
          strcor <- rbind(newstrcor, cens)
        dimens <- nrow(newstrcor)

        # strcor <- aggrStrata(strata, nvar, newsolution, censiti, dominio=dominio)

        # dimsamp <- nrow(strcor)
        # if (strcens == TRUE)
        #   strcor <- rbind(strcor, cens)
        # dimens <- nrow(strcor)
        # alfa<-res[[2]]
        # res<-second_bethel_alfa(newstrcor, errors,alfa,realAllocation = realAllocation)
        res<-SimAnn::bethel_alfa2(newstrcor, errors,alfa=res[[2]],realAllocation = TRUE)
        soluz <- res[[1]]


        newtot <- sum(soluz)

        delta<- newtot - tot



        if (delta <=0) {
          solution<-newsolution
          tot<-newtot
          strcor<-newstrcor
        }else{ if ((exp((-delta/(newTemp))) > runif(1))){
          deltaList <- c(deltaList, delta)

          solution<-newsolution
          tot<-newtot
          strcor<-newstrcor
        }else{
          solution<-solution
          tot<-tot
          strcor<-strcor
        }
        }




        x_axis <- c(x_axis, iter)
        y_axis <- c(y_axis, tot)
        L_chain<-c(L_chain,tot)
        if(tot==min(y_axis)){
          best_tot<-tot
          best_sol<-solution
        }


        current_time<-proc.time() - start_time

        i<-i+1

        if (j ==1){(k=ceiling(k*0.99))}else{k=1}
      }



      j<-j+1
      #plot(y_axis,type="l",xlab="Number of Solutions",ylab="Sample Size")

      newTemp<-decrement_constant*newTemp
    }

  }


  #dev.off()
  if (verbose)
    cat(" done.\n")

  result = list(maxiterations= jsize, j_reached=j-1, solutions_generated =iter,
                solution = best_sol,
                best = best_tot,samplesize=y_axis, deltaList=deltaList, Final_temperature=newTemp)


  return(result)
}
