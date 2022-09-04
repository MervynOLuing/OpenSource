parallelSAA<-function(strata,errors,
                                sugg,
                                Temp=0.01,initialStrata=NA, decrement_constant=0.99, end_time =Inf,
                                showSettings = FALSE, jsize=3,length_of_markov_chain =100,
                                verbose = FALSE, dominio=NULL,minnumstrat=NULL,kmax_percent=0.025,ProbNewStratum=NA,
                                strcens=FALSE,writeFiles=FALSE, showPlot=TRUE, minTemp = 0.0005, realAllocation=TRUE){

  if (writeFiles == TRUE) {
    dire <- getwd()
    direnew <- paste(dire, "/output", sep = "")
    if (dir.exists(direnew))
      unlink(direnew, recursive = TRUE)
    if (!dir.exists(direnew))
      dir.create(direnew)
  }

  #thanks to: https://towardsdatascience.com/getting-started-with-parallel-programming-in-r-d5f801d43745
  dom<-unique(strata$DOM1)
  ndom<-length(dom)
  cores<-(parallel::detectCores())-1
  if (ndom < cores) {cores <-ndom}
  cl <- parallel::makeCluster(cores)
  # Activate cluster for foreach library
  doParallel::registerDoParallel(cl)

  #ptm <- proc.time()
  r <- foreach::foreach(i = 1:ndom,
                        .combine = rbind
                        #,.packages = c("SimAnn")
  ) %dopar% {

    SAA(strata[which(strata$DOM1==dom[i]),],errors[i,],
                           #suggestion = sugg[which(sugg$domainvalue==dom[i]),],
                           sugg[which(sugg$domainvalue==dom[i]),],
                           Temp,initialStrata[i], decrement_constant, end_time,
                           showSettings = FALSE, jsize,length_of_markov_chain,
                           verbose = FALSE, dominio=dom[i],minnumstrat,kmax_percent,ProbNewStratum,
                           strcens=FALSE,writeFiles=FALSE, showPlot=TRUE, minTemp, realAllocation=TRUE)


  }
  #  time_foreach[3]
  # Stop cluster to free up resources
  parallel::stopCluster(cl)
  #  proc.time() - ptm
  if(showPlot==TRUE){

    if(ndom==1){
      plot(unlist(r$samplesize),type="l",xlab="Number of Solutions",ylab="Sample Size")
      title(paste("Domain #", dom[1], " - Sample cost",

                  round(min(unlist(r$samplesize)), 2)),

            col.main = "red")

    }else{
      for (i in 1:ndom){
        plot(unlist(r[,6][i]),type="l",xlab="Number of Solutions",ylab="Sample Size")

        title(paste("Domain #", dom[i], " - Sample cost",

                    round(min(unlist(r[,6][i])), 2)),

              col.main = "red")
      }
    }

  }

  if (writeFiles == TRUE) {

    if(ndom==1){

      stmt <- paste("png(filename = file.path(direnew, 'plotdom", 1, ".png'),height=5, width=7, units='in', res=144)", sep = "")

      eval(parse(text = stmt))



      plot(unlist(r$samplesize),type="l",xlab="Number of Solutions",ylab="Sample Size")

      title(paste("Domain #", i, " - Sample cost",

                  round(min(unlist(r$samplesize)), 2)),

            col.main = "red")

      if (writeFiles == TRUE)  dev.off()

      v <- as.vector(unlist(r$solution))

      if (class(v) == "matrix") v <- as.vector(v[1, ])

      stmt <- paste("write.table(v, file.path(direnew, 'solution", dom[i],

                    ".txt'),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)",

                    sep = "")

      eval(parse(text = stmt))

      s <-r$best


      stmt <- paste("write.table(s, file.path(direnew, 'Sample Size", dom[i],

                    ".txt'),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)",

                    sep = "")

      eval(parse(text = stmt))


    }else{

      for (i in 1:ndom){


        stmt <- paste("png(filename = file.path(direnew, 'plotdom", i, ".png'),height=5, width=7, units='in', res=144)", sep = "")

        eval(parse(text = stmt))



        plot(unlist(r[,6][i]),type="l",xlab="Number of Solutions",ylab="Sample Size")

        title(paste("Domain #", i, " - Sample cost",

                    round(min(unlist(r[,6][i])), 2)),

              col.main = "red")

        if (writeFiles == TRUE)  dev.off()

        v <- as.vector(unlist(r[,4][i]))

        if (class(v) == "matrix") v <- as.vector(v[1, ])

        stmt <- paste("write.table(v, file.path(direnew, 'solution", dom[i],

                      ".txt'),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)",

                      sep = "")

        eval(parse(text = stmt))

        s <-unlist(r[,5])[i]


        stmt <- paste("write.table(s, file.path(direnew, 'Sample Size", dom[i],

                      ".txt'),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)",

                      sep = "")

        eval(parse(text = stmt))

      }
    }


  }





  # solution<-NULL
  # for (i in 1:ndom){
  #  solution<-c(solution,as.vector(unlist(r[,4][1])))

  #  }

  if (ndom==1){
    result = list(maxiterations= as.vector(unlist(r$maxiterations)), j_reached=as.vector(unlist(r$j_reached)), solutions_generated =as.vector(unlist(r$solutions_generated)),
                  solution = r$solution,
                  best = as.vector(unlist(r$best)),samplesize=r$samplesize, deltaList=r$deltaList, Final_temperature=r$Final_temperature)
  }else{

    result = list(maxiterations= as.vector(unlist(r[,1])), j_reached=as.vector(unlist(r[,2])), solutions_generated =as.vector(unlist(r[,3])),
                  solution = r[,4],
                  best = as.vector(unlist(r[,5])),samplesize=r[,6], deltaList=r[,7], Final_temperature=r[,8])
  }

  return(result)
}

