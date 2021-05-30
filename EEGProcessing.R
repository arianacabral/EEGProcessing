# -------------------------------------------------------------------------
# Desenvolvido por: Adriano de Oliveira Andrade (adriano@ufu.br)
#                   Ariana Moura Cabral (arianacabral57@gmail.com)
# Descricao: Analise e processamento de sinais de EEG 
# Data: 29 de Maio de 2021
# -------------------------------------------------------------------------

EEG.Processing <- function(filename){
  
  # INTAN Toolbox
  source("read_Intan_RHD2000_file.R")
  
  # -------------------------------------------------------------------------
  # Packages
  
  # packages necessarias
  packages <- c("dygraphs","ggplot2", "stringr", "dplyr", "psd", "Rlibeemd", "tictoc", "openxlsx")
  
  # instalando as packages ainda nao instaladas
  if (any(packages %in% rownames(installed.packages()) == FALSE)){install.packages(packages[!packages %in% rownames(installed.packages())])}
  
  # lendo todas as packages
  invisible(lapply(packages, library, character.only = TRUE))
  invisible(rm(packages))
  
  # -------------------------------------------------------------------------
  SplitOverlap <- function(df, Nw, fs) {
    
    N <- (Nw/1000)*fs # Window with the number of samples 
    
    overlap = ceiling(N/2) #overlap of window 
    
    # splitting and extracting the features
    
    spliptOver <- function(y, seg.length, overlap) {
      
      vec <- 1:length(y)
      
      starts = seq(1, length(vec), by=seg.length-overlap)
      
      # handling exceptions
      #remove <- starts[starts == length(vec)]
      remove <- starts[length(vec) - starts < 100]
      starts %in% remove
      starts <- starts [! starts %in% remove]
      
      ends   = starts + seg.length - 1
      
      # handling exceptions
      ends[ends >= length(vec)] = length(vec) #garantindo que o indice final nao ultrapasse o tamanho do vetor
      
      
      v.mav <- vector(length = length(starts),mode="integer")
      v.cv <- vector(length = length(starts),mode="integer")
      v.zcr <- vector(length = length(starts),mode="integer")
      v.var <- vector(length = length(starts),mode="integer")
      v.mob <- vector(length = length(starts),mode="integer")
      v.comp <- vector(length = length(starts),mode="integer")
      
      index <- data.frame(starts,ends)
      
      
      for(i in 1:length(v.comp)){
        
        v.mav[i] <- sum(abs(y[starts[i]:ends[i]]))/length(y[starts[i]:ends[i]])
        
        v.cv[i] <- sd(y[starts[i]:ends[i]])/mean(y[starts[i]:ends[i]])
        
        v.zcr[i] <- zcr(y[starts[i]:ends[i]], wl = NULL, f = fs)
        
        v.var[i] <- var(y[starts[i]:ends[i]])
        
        v.mob[i] <- sqrt((var(diff(y[starts[i]:ends[i]])))/var(y[starts[i]:ends[i]]))
        
        v.comp[i] <- sqrt(var(diff(diff(y[starts[i]:ends[i]])))/ var(diff(y[starts[i]:ends[i]])))/v.mob[i]
        
      }
      
      
      ftr <- data.frame(v.mav,v.cv,v.zcr,v.var,v.mob,v.comp)
      ftr <- apply(ftr, 2, median)%>%t()
      
      return(ftr)
    }
    
    result <- spliptOver(df,N,overlap)
    
    return(result)
  }
  
  psf <- function(vec,fs){
    
    sss <- pspectrum(vec, verbose = FALSE, 
                     niter=10, AR=TRUE, x.frqsamp=fs, plot=FALSE) ##library(psd)
    
    return(data.frame("freq" = sss$freq, "spec" = sss$spec))
    
  }
  
  # -------------------------------------------------------------------------
  # Abrir o arquivo de EEG 
  
  X <- OpenIntanFile(filename = filename) 
  
  # -------------------------------------------------------------------------
  # Selecionar o canal desejado 
  
  df <- data.frame()
  
  for(c in 1:16){
    
    tic(".....")
    
    n_chan <- c # selecione o canal desejado (ha 16 canais de EEG, voce deve
    #                                         inserir um valor de 1 a 16)
    
    cat("\n Processando canal: ", n_chan,"\n\n")
    
    # Selecionando o canal 
    chan <- X[,c(1,n_chan+1)]
    chan[,1] <- chan[,1] - chan[1,1]
    
    cat("Reamostrando o sinal para 200 Hz \n\n")
    
    # Interpolando o sinal
    t <- seq(0,chan[length(chan[,1]),1],1/200) 
    chan.i <- spline(x = chan[,1], y = chan[,2], xout = t)
    chan.i <- data.frame(time = t, chan = chan.i$y)
    
    cat("Decompondo o sinal \n\n")
    
    # https://arxiv.org/pdf/1707.00487.pdf
    # https://www.ncl.ucar.edu/Document/Functions/Built-in/ceemdan.shtml
    
    imf <- ceemdan(chan.i[,2], num_imfs = 12, ensemble_size = 50)
    
    imf <- imf%>%as.data.frame()
    
    cat("Estimando as caracteristicas \n\n")
    
    # Calculo das caracteristicas
    
    featureExtract <- sapply(imf, function(i){
      SplitOverlap(df = i, Nw = 1000, fs = 200)
    })
    
    featureExtract <- featureExtract%>%as.data.frame()
    row.names(featureExtract) <- c("mav","cv","zcr","var","mob","comp")
    featureExtract <- featureExtract%>%t()%>%as.data.frame()
    featureExtract$chan <- n_chan
    featureExtract$filename <- basename(filename)
    featureExtract$imf <- ifelse(str_sub(row.names(featureExtract),1,3)=="IMF",
                                 str_c(str_sub(row.names(featureExtract),1,3),
                                       str_sub(row.names(featureExtract),5,-1)),
                                 row.names(featureExtract))
    
    cat("Caracteristicas estimadas para:", ifelse(str_sub(row.names(featureExtract),1,3)=="IMF",
                                                  str_c(str_sub(row.names(featureExtract),1,3),
                                                        str_sub(row.names(featureExtract),5,-1)),
                                                  row.names(featureExtract)), "\n\n")
    
    df <- rbind(df,featureExtract)
    
    toc()
    
    cat("\n")
    
    if(c == 16){
      write.xlsx(df,str_c(str_sub(basename(filename),1,-5),".xlsx"), colnames = TRUE, rownames = TRUE)
      cat("\n Arquivo salvo em:",getwd())
    }
  }
}