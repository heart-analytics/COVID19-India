#rm(list=ls())
library(tidyverse)
library(readxl)
library(readr)
library(magrittr)
library(glmnet)
theme_set(theme_minimal())

options(digits=12)
# Use a preset seed so test values are reproducable. 
#test.rng <- list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=c(979,1479,1542))
#old.seed <- setRNG(test.rng)

###########################################
# Set working directory
setwd("C:/Users/ssouyris/Box/COVID19/CodeAndData/India/")
###########################################

################################################
####### SIR + Mobility Regression Models #######
################################################

df_All <- read_csv("DistrictWise Data/dfIndia_to_models_districts_with_Mobility_for_paper.csv") 


df_All %>% ggplot(aes(x = ))


df_All %>% select(District_ID) %>% unique() %>% dim()

retail_and_recreation_percent_change_from_baseline 
grocery_and_pharmacy_percent_change_from_baseline 
parks_percent_change_from_baseline
transit_stations_percent_change_from_baseline 
workplaces_percent_change_from_baseline 

colnames(df_All)
colsMobility = seq(12,16)
View(df_All[ , colsMobility])

df_All$Mob = rowMeans(df_All[ , colsMobility], na.rm=TRUE)
df_All %>% ggplot(aes(Date, y = Mob, colour = District_ID)) +
  geom_line(size = 1)

##########################

mob_all_plot_quantile <- df_All %>% 
  mutate_at(c("id","block_id", "beta_null","beta_ext"), as.numeric) %>% 
  mutate(date = ymd(date)) %>% 
  group_by(date) %>% 
  do(data.frame(t(quantile(.$beta_ext, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))))) %>% 
  ggplot(aes(x = date, group = 1)) + theme_bw() +
  geom_ribbon(aes(ymin = X5., ymax = X95., fill = "05%-95%"), alpha = .25) + 
  geom_ribbon(aes(ymin = X25., ymax = X75., fill = "25%-75%"), alpha = .25) +
  geom_line(aes(y = X50.)) +
  scale_fill_manual(name = "", values = c("25%-75%" = "red", "05%-95%" = "blue")) + 
  labs(x ="Date", y = "Beta External") + theme(text = element_text(face="bold", size=sir_text_size),
                                               axis.text.x = element_text(face="bold", size = sir_text_size), axis.text.y = element_text(face="bold", size = sir_text_size)) +
  theme(legend.position="none")


##########################

df_All %>% write_csv("DistrictWise Data/dfIndia_to_models_districts_with_Mobility_for_paper.csv")

df_distance <- read_csv("DistrictWise Data/df_distance.csv")
df_Airports_edges <- read_csv("DistrictWise Data/District_Airport_Matched_Edges.csv") %>% select(id_j, id_i)

df_Airports_edges %>% dim()
df_Airports_edges %<>% drop_na() 
df_Airports_edges %>% dim()

df_All_locs <- df_All %>% select(id, State,District, P) %>% unique()

Dates <- unique(df_All$Date)

df_All %>% ggplot(aes(x = Date, y = dI, colour = District_ID)) +
                    geom_line()




locIDs <- unique(df_All$id)
locs_names <- df_All_locs %>% unite("id_State_District", id:District, remove = FALSE) %>% pull(id_State_District)
population <- df_All_locs %>% select(P) %>% pull(P)

nlocIDs <- locIDs %>% length()
nDates <- Dates %>% length()

df_lat_lon <- df_All %>% select(id, lat,lon) %>% unique()


#################################
#################################
df <- df_All
df$I_P <- df$I * 1 / df$P 
df$S_P <- df$S / df$P
df$dI_P <- df$dI / df$P

### t to t+1
for (i in locIDs){
  df_i <- df %>% filter(id == i)
  cS_P <- df_i %>% pull(S_P)
  cI_P <- df_i %>% pull(I_P)
  for(t in seq(nDates,2)){
    cS_P[t] <- cS_P[t-1]
    cI_P[t] <- cI_P[t-1]
  }
  df$S_P %<>% replace(df$id == i, cS_P)
  df$I_P %<>% replace(df$id == i, cI_P)
}

df %<>% select(id, Date, dI_P, S_P, I_P)
#df %>% group_by(id) %>% ggplot(aes(x = Date, y = S_P, colour = id)) + geom_line()

#######################################
############# Edges ###################
#######################################

d_cutoff <- 0.2
df_distance_edges <- df_distance %>% 
            filter(d <= d_cutoff) %>% 
              select(!d) 
              # %>% 
              # rbind(df_Airports_edges) %>%
              # unique()

df_distance_edges %>% group_by(id_i) %>% 
              summarise(n = n()) %>% 
              summary()


df_n_Airports_edges <- df_Airports_edges %>% group_by(id_i) %>% 
                      summarise(n = n()) %>% 
                      arrange(desc(n))

sum(df_n_Airports_edges$n)
  
# df_n_Airports_edges %>% filter(n == max(n))
# df_All_locs %>% filter(id == 2)

df_edges <- df_distance_edges %>% rbind(df_Airports_edges) %>% unique() %>% arrange(id_i)



df_n_edges <- df_edges %>% group_by(id_i) %>% 
              summarise(n = n()) 

df_n_edges %<>% rename(i = id_i) %>% arrange(desc(n))
df_n_edges

df_distance_edges %>% dim()
df_Airports_edges %>% dim()
df_edges %>% dim()








# 
# %>% 
# summary()

#######################################################
### Creation of matrices Iji for posterior reports ####
#######################################################
## Time windows
tw = 14
taus = Dates[seq(from = 13, to = nDates, by = tw)]
#taus = Dates[seq(from = tw+13, to = nDates, by = tw)]
ntaus <- length(taus)
nlocIDs_original <- nlocIDs
#alphas = c(0.5,0.9,0.99)
alphas = c(1)


Iji <- list()
for(idx in 1:nlocIDs){
  #idx <- 1
  i <- locIDs[idx]
  df_Lasso <- df
  df_i <- df_Lasso %>% filter(id == i)
  df_Lasso %<>% select(id, Date, I_P)
  df_xi <- pivot_wider(df_Lasso, names_from = c(id), values_from = c(I_P)) 
  
  ## select columns with edges
  edges_ji <-  df_edges %>% filter(id_i == i)  %>% arrange(id_j) %>% pull(id_j)
  df_xi <- df_xi[, c(1, (edges_ji +1)) ] ## first column is Date
  n_edges_ji <- length(edges_ji)
  
  ### Create df_Xi starting from taus[1]
  df_Xi_tau <- df_xi
  df_Xi_tau[(df_Xi_tau$Date > taus[1]) , seq(2, (n_edges_ji+1)) ] <- 0
  ## add names to columns 
  suffix = paste0("_",taus[1])
  colnames(df_Xi_tau) <- c("Date", paste0(edges_ji,suffix))
  df_Xi <- df_Xi_tau
  
  ## Complete  df_Xi
  for(t in 2:ntaus){
    df_Xi_tau <- df_xi
    df_Xi_tau[(df_Xi_tau$Date <= taus[t-1]) | (df_Xi_tau$Date > taus[t]) , seq(2, (n_edges_ji+1)) ] <- 0
    ## add names to columns 
    suffix = paste0("_",taus[t])
    colnames(df_Xi_tau) <- c("Date", paste0(edges_ji,suffix))
    df_Xi %<>% left_join(df_Xi_tau, by = c("Date"))
  }
  #View(df_Xi)
  ### Add row constraints 0 = Bb(t) - Bb(t+1)  
  ### Error_newrow = (0 - (Bb(t) - Bb(t+1)))
  ### Error_newrow = Bb(t+1) - Bb(t)
  
  ncols <- dim(df_Xi)[2]
  
  Iji[[idx]] <- as.matrix(df_Xi[, seq(2,ncols)])
}


#################################################################################################################
######  Regression   ###############
#################################################################################################################

if(1){
  print("### Regressions ###")
  list_YX <- list() ## object to save regression matrices 

  list_X <- list() ## object to save regression matrices 
  list_Y <- list() ## object to save regression matrices 
  
  cvfit <- list()
  Betas <- list()
  Predictions <- list()
  df_Betas <- list()
  df_Predictions <- list()
  model_reports <- list() 
  df_model_reports <- list()

  ######### Regressions #################
  for(idx in 1:nlocIDs){
    paste("Prepare Matrix list_YX:", idx,"****************************** ") %>% print()
    
    #idx <- 1
    i <- locIDs[idx]
  
    df_Lasso <- df
    df_i <- df_Lasso %>% filter(id == i)
    cS_P = df_i %>% pull(S_P)
    for (j in locIDs){
      df_Lasso$S_P %<>% replace(df_Lasso$id == j, cS_P)
    }
    df_Lasso$S_P.I_P <- df_Lasso$S_P * df_Lasso$I_P 
    df_Lasso %<>% select(id, Date, S_P.I_P)
    df_xi <- pivot_wider(df_Lasso, names_from = c(id), values_from = c(S_P.I_P)) 
  
    ## select columns with edges
    edges_ji <-  df_edges %>% filter(id_i == i)  %>% arrange(id_j) %>% pull(id_j)
    df_xi <- df_xi[, c(1, (edges_ji +1)) ] ## first column is Date
    n_edges_ji <- length(edges_ji)
    
    ### Create df_Xi starting from taus[1]
    df_Xi_tau <- df_xi
    df_Xi_tau[(df_Xi_tau$Date > taus[1]) , seq(2, (n_edges_ji+1)) ] <- 0
    ## add names to columns 
    suffix = paste0("_",taus[1])
    colnames(df_Xi_tau) <- c("Date", paste0(edges_ji,suffix))
    df_Xi <- df_Xi_tau
    
    ## Complete  df_Xi
    for(t in 2:ntaus){
      df_Xi_tau <- df_xi
      df_Xi_tau[(df_Xi_tau$Date <= taus[t-1]) | (df_Xi_tau$Date > taus[t]) , seq(2, (n_edges_ji+1)) ] <- 0
      ## add names to columns 
      suffix = paste0("_",taus[t])
      colnames(df_Xi_tau) <- c("Date", paste0(edges_ji,suffix))
      df_Xi %<>% left_join(df_Xi_tau, by = c("Date"))
    }
  
    ### Add row constraints 0 = Bb(t) - Bb(t+1)  
    ### Error_newrow = (0 - (Bb(t) - Bb(t+1)))
    ### Error_newrow = Bb(t+1) - Bb(t)
    
    list_Xi_links <- list()
    for (j in edges_ji){
      df_Xi_links <- df_Xi %>% filter(Date %in% taus)
      aux <- seq(2, dim(df_Xi_links)[2]) 
      df_Xi_links[,aux] <- 0
      for(t in 1:(ntaus-1)){
        suffix = paste0("_",taus[t])
        j_t <- paste0(j, suffix)
        suffix = paste0("_",taus[t+1])
        j_tplus <- paste0(j, suffix)
        # print(paste(j_t, j_tplus))
        
        df_Xi_links[ (df_Xi_links$Date == taus[t]), j_t] <- 1
        df_Xi_links[ (df_Xi_links$Date == taus[t]), j_tplus] <- -1 
      }
      df_Xi_links$Date <- paste0(j,"_",df_Xi_links$Date)
      list_Xi_links[[j]] <- df_Xi_links
    }
    
    jidx <- 1
    df_Xi_links <- list_Xi_links[[edges_ji[jidx]]]
    for (jidx in 2:length(edges_ji)){
      df_Xi_links %<>% rbind(list_Xi_links[[edges_ji[jidx]]])
    }
    
    ## Row bind Bind X covid and X links
    df_Xi %>% dim()
    df_Xi_links %>% dim()
    df_Xi$Date %<>% as.character()
    df_Xi %<>% rbind(df_Xi_links)
    
    X <- as.matrix(df_Xi[,seq(2,dim(df_Xi)[2])])
    Y <- df_i %>% pull(dI_P) 
    y_links <- rep(0, dim(df_Xi_links)[1])
    Y <- c(Y, y_links)
    
    list_X[[idx]] <- X
    list_Y[[idx]] <- Y
  
    df_Xi$Y <- Y
    
    df_Xi %<>% select(Date,Y,colnames(df_Xi)[!(colnames(df_Xi) %in% c("Date","Y"))])
    list_YX[[idx]] <-  df_Xi 
  }
  
  # lambda2 <- 0.05 ## good candidate
  # lambda2 <- 0.04 ## good candidate
  # lambda2 <- 0.035 ## good candidate
  # lambda2 <- 0.033 ## good candidate
  lambda2 <- 0.032 ## good candidate
  # lambda2 <- 0.03 ## good candidate
  # lambda2 <- 0.01 ## too small
  #lambda2 <- 0.1 ## too small
  lambda2p2 <- lambda2^2 
  # nlocIDs <- nlocIDs_original
  #nlocIDs <- 16 ## to run only for the first 10 districts
  ###### List to save Regressions 

  #whatlambda = "0"
  whatlambda = "min"
  #whatlambda = "min_most_div2"
  #whatlambda = "most_regularized"
  #plot(cvfit[[a]][[1]], xvar = "lambda", label = TRUE)
  
  

  
  
  ###########
  ###########
  if(1){
    
    for (a in 1:length(alphas)){
      cvfit[[a]] <- list()
    }

    for(idx in 1:nlocIDs){
      Y <- list_Y[[idx]]
      X <- list_X[[idx]]
      r1 <- (nDates + 1)
      r2 <- dim(X)[1]
      X[r1:r2, ] = X[r1:r2, ] * lambda2p2
      for (a in 1:length(alphas)){
        paste("Reg:", "a",a, "idx",idx, "****************************** ") %>% print()
        cvfit[[a]][[idx]] <- cv.glmnet(X, Y, lower.limits = 0.0, intercept = FALSE,alpha = alphas[a]) #, # penalty.factor = c(0,1,1,1)
      }
    }

    ## Save Betas and Predictions 
    
    paste("B. Get coefficients  Betas[[a]]") %>% print()
    for (a in 1:length(alphas)){
      Betas[[a]] <- list()
      Predictions[[a]] <- list()
      for(idx in 1:nlocIDs){
        i <- locIDs[idx]
        # paste("B. tw",tw,"a", a,  "idx", idx) %>% print()
        f <- cvfit[[a]][[idx]]
        if (sum(is.na(f)) == 0){
          lambdamin <- f$lambda.min
          lambda1se <- f$lambda.1se
          lambdahalf <- (lambdamin + lambda1se)/2
          #print(paste("lambdamin:", lambdamin))
          #print(paste("lambda1se:", lambda1se))
          #print(paste("lambdahalf:", lambdahalf))
          if(whatlambda == "0"){
            lambdaAux <- 0
          }else if(whatlambda == "min"){
            lambdaAux <- lambdamin
          }else if(whatlambda == "min_most_div2"){
            lambdaAux <- lambdahalf
          }else if(whatlambda == "most_regularized"){
            lambdaAux <- lambda1se
          }
          # print(paste("lambda:", lambdaAux))
          
          coefs = as.vector(coef(f, s = lambdaAux)) 
          coefs <- coefs[2:length(coefs)] # first coefficient is the intercept that is not used
          
          ## Save coefficients in data frame
          edges_ji <-  df_edges %>% filter(id_i == i)  %>% arrange(id_j) %>% pull(id_j)
          n_edges_ji <- length(edges_ji)
          
          ci <- c() 
          cDistricti <- c()
          cPi <- c()
          ctau <- taus
          cj <- c()
          cDistrictj <- c()
          cPj <- c()
          k <- 0 
          for(t in 1:ntaus){
            for(j in edges_ji){
              k <- k + 1 
              ctau[k] <- taus[t]
              ci[k] <- i
              cj[k] <- j
              cDistricti[k] <- locs_names[i]
              cDistrictj[k] <- locs_names[j]
              cPi[k] <- population[i]
              cPj[k] <- population[j]
            }
          }
          Betas[[a]][[idx]] <- data.frame("tau" = ctau, "i" = ci, "j" = cj, "District_i" = cDistricti,
                                            "District_j" = cDistrictj, "pi"= cPi, "pj"= cPj, "B" = coefs)
          #######################################################
          ##################### Save predictions ##################
          df_Xi <- list_YX[[idx]]
          y <- df_Xi %>% pull(Y)
          X <- as.matrix(df_Xi %>% select(colnames(df_Xi)[!(colnames(df_Xi) %in% c("Date","Y"))]))
          prediction <- predict(f, X, s = lambdaAux)
          df_Pred <- df_Xi %>% select(Date, Y)
          df_Pred$i <- idx
          df_Pred$Location <- locs_names[idx]
          df_Pred$Prediction <- c(prediction)
          df_Pred[is.na(df_Pred)] <- 0
          df_Pred %<>% filter(Date <= taus[ntaus])
          df_Pred$Date %<>% as.Date()
          df_Pred %<>% rename(Actual = Y)
          df_Pred$Error <- df_Pred$Actual - df_Pred$Prediction
          df_Pred %<>% pivot_longer(c(Actual, Prediction,Error), names_to = "Series")
          Predictions[[a]][[idx]] <- df_Pred
          
        }else{
          paste("f == NA, idx", idx) %>% print()
        }
      }
    }
  
    for (a in 1:length(alphas)){
      params_str <- paste0("_tw_",tw,"_alpha_",alphas[a],"_lambda1_",whatlambda,"_lambda2_",lambda2)
      
      df_Betas[[a]] <- Betas[[a]][[1]]
      for(idx in 2:nlocIDs){
        df_Betas[[a]] %<>% rbind(Betas[[a]][[idx]])
      }
      df_Betas[[a]] %<>% arrange(i, tau, j)
      
      file_name = paste0("Betas/Betas",params_str,".csv") 
      df_Betas[[a]] %>% write_csv(file_name)
      df_Betas[[a]] %<>% as_tibble()
      
      
      df_Predictions[[a]] <- Predictions[[a]][[1]]
      for(idx in 2:nlocIDs){
        df_Predictions[[a]] %<>% rbind(Predictions[[a]][[idx]])
      }
      df_Predictions[[a]] %<>% arrange(i, Date)
      
      df_Predictions[[a]]$wave <- "None"
      df_Predictions[[a]][df_Predictions[[a]]$Date <= "2020-12-31", ]$wave <- "First"
      df_Predictions[[a]][df_Predictions[[a]]$Date >= "2021-03-01", ]$wave <- "Second"
    
      file_name = paste0("Betas/Predictions",params_str,".csv") 
      
      df_Predictions[[a]] %>% write_csv(file_name)
      df_Predictions[[a]] %<>% as_tibble()
      
      
      # df_Predictions[[a]] %>% pivot_wider(names_from  = Series, values_from = value) %>% 
      #                           mutate(Error_Pct = Error/Actual) %>% 
      #                          ggplot(aes(x =Date, y = Error_Pct)) + geom_point() + ylim(c(-10,3))
      #   
        
      
      if(1){
        print("###C. Plot max_edges Betas and Predictions  ###")
        
        max_edges = 3
        df_Betas_i <- df_Betas[[a]] %>% left_join(df_n_edges, by = c("i")) %>% 
                                 filter( n <= max_edges) %>%
                                 arrange(desc(n), desc(pi), tau, desc(pj)) 
        
        df_Betas_i_aux <- df_Betas_i %>% select(i,pi) %>% unique()
        df_Betas_i_aux %<>% mutate(aux_order = seq(1,dim(df_Betas_i_aux)[1])) %>% select(i, aux_order)
        df_Betas_i %<>% left_join(df_Betas_i_aux, by = c("i"))
          
        df_Betas_i %<>% filter(aux_order <= 16) %>% arrange(i, tau, j)  %>% 
                                    mutate(District_i=as.factor(District_i))
        
        file_name = paste0("Betas/Betas_edges_",max_edges,params_str) 
        df_Betas_i %>% ggplot(aes(x = tau, y = B, colour = District_j)) + 
                      geom_line(size=1, show.legend = FALSE) + labs(title=file_name) +
                      facet_wrap(~District_i, ncol = 4, scales = "free_y") +
                      theme(text = element_text(size = 40))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 49, height = 30)
        
        ### Predictions 
        df_Predictions_i <- df_Predictions[[a]] %>% 
                          left_join(df_Betas_i_aux, by = c("i")) %>%
                            filter(aux_order <= 16)
        
        file_name = paste0("Betas/Preds_edges_",max_edges,params_str) 
        df_Predictions_i %>% filter(Series != "Error")%>% ggplot(aes(x = Date, y= value, colour =  Series)) + 
                              geom_line() +
                            geom_line(size=2) + labs(title=file_name) +
                              facet_wrap(~Location, ncol = 4, scales = "free_y") +
                              theme(text = element_text(size = 40))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 49, height = 30)
        
        file_name = paste0("Betas/Errors_edges_",max_edges,params_str) 
        df_Predictions_i %>% filter(Series == "Error") %>% ggplot(aes(x = Date, y= value, colour =  Series)) + 
          geom_point(size=1, show.legend = FALSE) +
          labs(title=file_name, y = "Residual") +
          facet_wrap(~Location, ncol = 4, scales = "free_y") +
          theme(text = element_text(size = 40))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 49, height = 30)
        
        file_name = paste0("Betas/Hist_errors_waves_edges_",max_edges,params_str) 
        df_Predictions_i %>% filter(wave != "None") %>% 
          filter(Series == "Error")  %>% 
          ggplot(aes(x= value, colour = wave, fill = wave)) + 
          geom_density(size=1, alpha = 0.1) +
          labs(title=file_name, x = "Residual") +
          facet_wrap(~Location, ncol = 4, scales = "free_y") +
          theme(text = element_text(size = 40))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 49, height = 30)
        
      }
      if(1){
        print("### D. Plot all districts Betas and Predictions  ###")
        
        df_Betas_i <- df_Betas[[a]] %>% 
                      arrange(desc(pi), tau, desc(pj)) 
        
        df_Betas_i %<>% filter(i <= 16) %>% arrange(i, tau, j)  %>% 
          mutate(District_i=as.factor(District_i))
        
        file_name = paste0("Betas/Betas_biggest_districts", params_str) 
        
        df_Betas_i %>% ggplot(aes(x = tau, y = B, colour = District_j)) + 
          geom_line(size=2, show.legend = FALSE) + labs(title=paste0(file_name)) +
          facet_wrap(~District_i, ncol = 4, scales = "free_y") +
          theme(text = element_text(size = 40))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 49, height = 30)
        
        ### Predictions 
        df_Predictions_i <- df_Predictions[[a]] %>% 
                          filter(i <= 16)
    
        file_name = paste0("Betas/Preds_biggest_districts", params_str) 
        df_Predictions_i %>% filter(Series != "Error")%>% ggplot(aes(x = Date, y= value, colour =  Series)) + 
          geom_line() +
          geom_line(size=1) + labs(title=file_name) +
          facet_wrap(~Location, ncol = 4, scales = "free_y") +
          theme(text = element_text(size = 40))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 49, height = 30)
        
        file_name = paste0("Betas/Errors_biggest_districts", params_str) 
        df_Predictions_i %>% filter(Series == "Error")%>% ggplot(aes(x = Date, y= value, colour =  Series)) + 
          geom_point(size=1, show.legend = FALSE) +
          labs(title=file_name, y = "Error") +
          facet_wrap(~Location, ncol = 4, scales = "free_y") +
          theme(text = element_text(size = 40))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 49, height = 30)
        
        file_name = paste0("Betas/Errors_all_districts", params_str) 
        df_Predictions_i %>% filter(Series == "Error")%>% ggplot(aes(x = Date, y= value, colour =  Series)) + 
          geom_point(size=1, show.legend = FALSE) +
          #geom_ribbon(aes(ymin=min(df_Predictions_i$value), max(df_Predictions_i$value)), linetype=2, alpha=0.1) +
          labs(title=file_name, y = "Residual") +
          theme(text = element_text(size = 20))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 15, height = 15)
        
        file_name = paste0("Betas/Hist_errors_waves_biggest_districts",params_str) 
        df_Predictions_i %>% filter(wave != "None") %>% 
          filter(Series == "Error")  %>% 
          ggplot(aes(x= value, colour = wave, fill = wave)) + 
          geom_density(size=1, alpha = 0.1) +
          labs(title=file_name, x = "Residual") +
          facet_wrap(~Location, ncol = 4, scales = "free_y") +
          theme(text = element_text(size = 40))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 49, height = 30)

        file_name = paste0("Betas/Hist_errors_waves_alldistricts",params_str) 
        df_Predictions_i %>% filter(wave != "None") %>% 
          filter(Series == "Error")  %>% 
          ggplot(aes(x= value, colour = wave, fill = wave)) + 
          geom_density(size=1, alpha = 0.1) +
          labs(title = file_name, x = "Residual") + # xlim() + 
          theme(text = element_text(size = 40))
        file_name = paste0(file_name,".pdf") 
        ggsave(file_name, width = 15, height = 15)
      }
      
      #####
      if(1){
        # df_Betas[[a]]
        print("### E. Create Reports ###")
        
        ############################
        #### Create reports ########
        ############################
        tau_to_Date <- Dates
        idxtaus <- 1
        for(t in 1:nDates){
          if(tau_to_Date[t] <  taus[idxtaus]){
            tau_to_Date[t] = taus[idxtaus]
          }else{
            tau_to_Date[t] = taus[idxtaus]
            idxtaus <- idxtaus + 1
          }
        }
        model_reports[[a]] <- list()
        for(idx in 1:nlocIDs){
        
          #idx <- 5
          i <- locIDs[idx]
          df_i <- df %>% filter(id == i)
          Y <- df_i %>% pull(dI_P) 
          cS_P <- df_i %>% pull(S_P)
          Iactual <- Y / cS_P 
      
          edges_ji <-  df_edges %>% filter(id_i == i)  %>% arrange(id_j) %>% pull(id_j)
          Iaux <- Iji[[idx]]
          betas <- df_Betas[[a]] %>% filter(i == idx) %>% pull(B)
          
          ncolsIaux <- dim(Iaux)[2]
        
          ledges <-  length(edges_ji)
          i_in_edges_ji <- which(edges_ji %in% i)
          columns_index_i_in_Iaux <- rep(0,length(taus))
          columns_index_i_in_Iaux[1] <- i_in_edges_ji
          for(tau in 2:length(taus)){
            columns_index_i_in_Iaux[tau] <- columns_index_i_in_Iaux[tau-1]  +  ledges
          }
          
          indexi <- (seq(1, ncolsIaux) %in% columns_index_i_in_Iaux)
          indexj <- !(seq(1, ncolsIaux) %in% columns_index_i_in_Iaux)
          
          Icolnames_i <- colnames(Iaux)[indexi]
          Icolnames_j <- colnames(Iaux)[indexj]
          
          Ii <- Iaux[, indexi]
          betasi <- betas[indexi]
          Iinternal <- c(Ii %*% betasi)
          
          Ij <- Iaux[, indexj]
          betasj <- betas[indexj]
          Iexternal <- c(Ij %*% betasj)
          Residual <- Iactual - Iinternal - Iexternal
          Prediction <- df_Predictions[[a]] %>% filter((i == idx) & (Series == "Prediction")) %>% pull(value)
        
          model_reports[[a]][[idx]] <- data.frame(i = rep(i, nDates), District = rep(locs_names[i], nDates), 
                                    tau_to_Date,Date = Dates, Iactual, Iinternal, Iexternal, Residual, Prediction, S = cS_P, Prediction_S = (Prediction/cS_P)) %>% as_tibble()

        }
        
        df_model_reports[[a]] <- model_reports[[a]][[1]] 
        for(idx in nlocIDs){
          df_model_reports[[a]] %<>% rbind(model_reports[[a]][[idx]]) 
        }
        df_model_reports[[a]]$Date %<>% as.Date()
        df_model_reports[[a]]$wave <- "None"
        df_model_reports[[a]][df_model_reports[[a]]$Date <= "2020-12-31", ]$wave <- "First"
        df_model_reports[[a]][df_model_reports[[a]]$Date >= "2021-03-01", ]$wave <- "Second"
        file_name = paste0("Betas_Models_Comparation/Disaggregate_edges_model_eq8",params_str,".csv") 
        df_model_reports[[a]] %>% write_csv(file_name)
        
        # df_to_plot <- df_model_reports[[a]] %>% group_by(Date) %>% summarise(Iinternal = sum(Iinternal), Iexternal = sum(Iexternal), 
        #                                                        Iactual = sum(Iactual), Prediction_S = sum(Prediction_S)) %>%
        #                       pivot_longer(cols=c(Iinternal, Iexternal, Iactual, Prediction_S), names_to = "I")
        # 
        # boxplot(Iinternal ~ tau_to_Date, data = df_model_reports[[a]], ylim = c(0, 0.0015))
        # boxplot(Iexternal ~ tau_to_Date, data = df_model_reports[[a]], ylim = c(0, 0.0015))
        # boxplot(Residual * S ~ tau_to_Date, data = df_model_reports[[a]], ylim = c(0, 0.0015))
        # 
        # 
        # df_to_plot <- df_model_reports[[a]] %>% group_by(Date) %>% summarise(Iinternal = sum(Iinternal), Iexternal = sum(Iexternal)) %>%
        #     pivot_longer(cols=c(Iinternal, Iexternal ), names_to = "I")
        # 
        # View(df_model_reports[[a]])
        # 
        # df_to_plot %>% ggplot(aes(Date, value, colour = I)) + geom_line()
        # 
        # df_Betas[[a]]$type <- "internal"
        # 
        # 
        # df_Betas[[a]]   $type <- "internal"
        # 
        # df_BetasInternal <-  df_Betas[[a]] %>% filter(i == j) 
        # df_BetasExternal <-  df_Betas[[a]] %>% filter(i != j)
        # 
        # df_BetasExternal$type <- "external"
        #  
        # df_to_plot <- df_BetasInternal %>% rbind(df_BetasExternal)
        # 
        # df_to_plot$tau
        # 
        # df_to_plot %>% ggplot(aes(y= B, colour = tau)) + 
        #               geom_boxplot() + 
        #               facet_wrap(~type)
        #               
        # boxplot(B ~ tau, data = df_BetasInternal)
        # boxplot(B ~ tau, data = df_BetasExternal, ylim = c(0, 0.025))
      }
    }
  }
}
