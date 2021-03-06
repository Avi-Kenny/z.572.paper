# Title: "572 Paper Reproduction"
# Author: Avi Kenny
# Date: 2020-06-19



#################.
##### Setup #####
#################.

# Set working directory
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  setwd("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Course materials/BIOST 572/Final project/z.572.paper/R")
} else {
  setwd("z.572.paper/R")
}

# Load packages
{
  library(sp)
  library(rgeos)
  library(sqldf)
  library(spdep)
  library(parallel)
  library(rgdal)
  library(simba)
  library(ggplot2)
  library(dplyr)
  library(magrittr)
  library(permute)
  library(maps)
  library(maptools)
  library(mvtnorm)
  library(rjags) # !!!!! On hutch load module: JAGS/4.2.0-foss-2016b
  library(MASS)
  library(pbapply)
  library(ngspatial)
  library(Matrix)
  library(INLA)

  # INLA:::inla.dynload.workaround()
}

# Load functions
{
  source("adj_from_gis_us.R")
  source("adj_from_graph.R")
  source("generate_dataset.R")
  source("generate_graph_grid.R")
  source("generate_graph_path.R")
  source("generate_matern_cov.R")
  source("generate_mtx.R")
}

# Set code blocks to run
{
  run_testing <- FALSE
  run_fig1 <- FALSE
  run_fig2 <- FALSE
  run_fig345 <- FALSE
  run_fig6 <- FALSE
  run_fig7 <- FALSE
  run_fig8 <- TRUE
  run_graphs <- FALSE
  run_ima <- FALSE
}



###########################################################.
##### Compile simulation results into a single object #####
###########################################################.

if (FALSE) {

  # Merge *.simba files
  sims <- list.files(
    path = "../simba.out/simba.out",
    pattern = "*.simba",
    full.names = TRUE,
    recursive = FALSE
  )
  print(length(sims))
  sim <- NULL
  for (s in sims) {
    s <- readRDS(s)
    if (is.null(sim)) { sim <- s } else { sim <- merge(sim, s) }
  }
  saveRDS(sim, file="../simba.out/sim_fig8_16200_6.16.simba")

}



#############################.
##### Generate figure 1 #####
#############################.

if (run_fig1) {

  rhos <- seq(from=0.1, to=0.9, by=0.1)
  rho_out <- c()
  adj_mtx_path <- generate_graph_path() %>% adj_from_graph()
  adj_mtx_grid <- generate_graph_grid() %>% adj_from_graph()
  adj_mtx_US <- adj_from_gis_us()

  for (adj_mtx in list(adj_mtx_path, adj_mtx_grid, adj_mtx_US)) {

    for (rho in rhos) {

      # Generate average neighbor-pair correlation for CAR
      mtx_car <- generate_mtx(model="CAR", adj_mtx=adj_mtx, rho=rho)
      cov_car <- solve(mtx_car$Q)
      cor <- 0
      count <- 0
      for (i in 1:length(mtx_car$neighbors)) {

        nbs <- mtx_car$neighbors[[i]]
        if (nbs[1]!=0) {
          for (nb in 1:length(nbs)) {
            j <- nbs[nb]
            cor <- cor + ( cov_car[i,j] /
                             (sqrt(cov_car[i,i]) * sqrt(cov_car[j,j])) )
            count <- count+1
          }
        }

      }
      anpc_car <- cor/count
      rho_out <- c(rho_out, anpc_car)

      # Generate average neighbor-pair correlation for DAGAR
      mtx_dagar <- generate_mtx(model="DAGAR", adj_mtx=adj_mtx, rho=rho)
      cov_dagar <- solve(mtx_dagar$Q)
      cor <- 0
      count <- 0
      for (i in 2:length(mtx_dagar$neighbors)) {

        nbs <- mtx_dagar$neighbors[[i]]
        if (length(nbs)!=0) {
          for (nb in 1:length(nbs)) {
            j <- nbs[nb]
            cor <- cor + ( cov_dagar[i,j] /
                             (sqrt(cov_dagar[i,i]) * sqrt(cov_dagar[j,j])) )
            count <- count+1
          }
        }

      }
      anpc_dagar <- cor/count
      rho_out <- c(rho_out, anpc_dagar)

    }

  }

  # Plot results
  # Export 900w x 400h
  ggplot(
    data = data.frame(
      x = rep(rep(rhos, each=2), 3),
      y = rho_out,
      grp = rep(c("CAR","DAGAR"), 27),
      type = rep(c(
        "(a) Path graph of length 100",
        "(b) 10 x 10 grid",
        "(c) 48 contiguous US states"
      ), each=18)
    ),
    aes(x=x, y=y, group=grp, color=as.factor(grp), shape=as.factor(grp),
        linetype=as.factor(grp))
  ) +
    geom_segment(x=0, y=0, xend=1, yend=1, color="grey") +
    geom_line() +
    geom_point(size=3) +
    scale_color_manual(values=c("firebrick", "darkolivegreen4")) +
    xlim(0.1,0.9) +
    ylim(0,1) +
    facet_wrap(~type, ncol=3) +
    scale_shape_manual(values=c(16,17)) +
    scale_linetype_manual(values=c(1,2)) +
    labs(x="rho", y="Average neighbor correlation", color="Model",
         shape="Model", linetype="Model") +
    theme(legend.position = "bottom")

}



#############################.
##### Generate figure 2 #####
#############################.

if (run_fig2) {

  rd_path <- function(p) {
    sqrt(
      (4*p^8 + 2*p^4) / ( (3+6*p^2+p^4)^2 + (18*p^2*(1+p^2)^2) + (2*p^4) )
    )
  }

  s <- function(p) {
    1 / ( 1 + (1-1)*p^2 ) +
    2 / ( 1 + (2-1)*p^2 ) +
    3 / ( 1 + (3-1)*p^2 ) +
    4 / ( 1 + (4-1)*p^2 )
  }

  rd_grid <- function(p) {
    sqrt(
      (
        p^4 * ( s(p)/5 - 2/(1+p^2) )^2 +
        2 * ( 1/3 - s(p)/30 - p^2/(1+p^2) )^2 +
        12 * ( 1/6 - s(p)/60 )^2
      ) /
      (
        ( 1 + p^2 + p^2*(s(p)/5) )^2 +
        4 * p^2 +
        20 * ( 1/6 - s(p)/60 )^2
      )
    )
  }

  # Generate data
  n <- 101
  rhos <- seq(from=0, to=1, length.out=n)

  # Plot results
  # Export 800w x 450h
  ggplot(
    data = data.frame(
      x = rep(rhos, 2),
      y = c(rd_path(rhos), rd_grid(rhos)),
      grp = rep(c("(a) Path graph","(b) Grid graph"), each=n)
    ),
    aes(x=x, y=y, group=grp)
  ) +
    geom_line() +
    geom_text(
      data = data.frame(
        x = rep(seq(from=0, to=1, by=0.25), 2),
        y = c(
          round(rd_path(c(0, 0.25, 0.5, 0.75, 1)),2),
          round(rd_grid(c(0, 0.25, 0.5, 0.75, 1)),2)
        ),
        grp = rep(c("(a) Path graph", "(b) Grid graph"), each=5)
      ),
      aes(x=x, y=y, label=y),
      size = 2
    ) +
    xlim(0,1) +
    facet_wrap(~grp, ncol=2) +
    labs(x="rho", y="Relative difference")

}



#######################################.
##### Testing: temp C and L lists #####
#######################################.

if (run_testing) {

  C = list(
    tau_w = 0.25,
    adj_mtx_path = generate_graph_path() %>% adj_from_graph(),
    adj_mtx_grid = generate_graph_grid() %>% adj_from_graph(),
    adj_mtx_US = adj_from_gis_us()
  )

  L = list(
    type = "path",
    rho = 0.9,
    cov_type = "Matern",
    model = "iCAR",
    p_order = "none" # "NE", "ZZ"
  )

}



#############################.
##### Simulation script #####
#############################.

if (run_fig345 || run_fig6 || run_fig7 || run_fig8) {

  script_mainsim <- function(L,C) {

    # Set variables depending on `type`
    switch(
      L$type,
      "path" = {
        k <- 100
        adj_mtx <- C$adj_mtx_path
      },
      "grid" = {
        k <- 100
        adj_mtx <- C$adj_mtx_grid
      },
      "US" = {
        k <- 48
        adj_mtx <- C$adj_mtx_US
      }
    )

    # Generate data
    data <- generate_dataset(
      type = L$type,
      tau_w = C$tau_w,
      cov_type = L$cov_type,
      phi = -1*log(L$rho),
      rho = L$rho
    )

    # Create JAGS code and objects

    if (L$model == "iCAR") {

      # # Run INLA model
      # data_inla <- data.frame(
      #   id=c(1:k), y=data$y, w=data$w, x1=data$x1, x2=data$x2
      # )
      # # hyper <- list(prec=list(param=c(2,1)))
      # inla_model <- inla(
      #   y ~ x1 + x2 -1 + f(id, model="besag", graph=adj_mtx),
      #   # y ~ x1 + x2 -1 + f(id, model="besag", graph=adj_mtx, hyper=hyper),
      #   data = data_inla
      # )
      #
      # # Get summary stats
      # w_hat <- inla_model$summary.random$id[["0.5quant"]]
      # rho_hat <- NA
      # ci_l_rho <- NA
      # ci_u_rho <- NA
      # beta1_hat <- inla_model$summary.fixed["x1","mean"]
      # beta1_se <- inla_model$summary.fixed["x1","sd"]
      # beta2_hat <- inla_model$summary.fixed["x2","mean"]
      # beta2_se <- inla_model$summary.fixed["x2","sd"]
      # sigma2_e_hat <- 1/inla_model$summary.hyperpar[1,"0.5quant"]
      # sigma2_e_se <- (
      #   sigma2_e_hat - (1/inla_model$summary.hyperpar[1,"0.975quant"])
      # ) / 1.96
      # mse <- mean((data$w-w_hat)^2)

      mtx <- generate_mtx(L$model, adj_mtx)
      D <- mtx$D
      n_i <- mtx$n_i

      # JAGS model code
      jags_code <- quote("
          model {
            for (i in 1:k) {
              y[i] ~ dnorm(x1[i]*beta1 + x2[i]*beta2 + w[i], tau_e)
            }

            w ~ dmnorm(mu0, Q)
            Q <- tau_w * (D-(0.99*A))

            beta1 ~ dnorm(0, 0.0001)
            beta2 ~ dnorm(0, 0.0001)
            tau_w ~ dgamma(2, 1)
            tau_e ~ dgamma(2, 0.1)
          }
        ")

    }

    if (L$model == "Scaled iCAR") {

      mtx <- generate_mtx(L$model, adj_mtx)
      D <- mtx$D
      n_i <- mtx$n_i

      sigma2ref <- exp(mean(log(diag(ginv(D - adj_mtx)))))
      sigma2ref <- 1/sigma2ref

      # JAGS model code
      jags_code <- quote("
          model {
            for (i in 1:k) {
              y[i] ~ dnorm(x1[i]*beta1 + x2[i]*beta2 + w[i], tau_e)
            }

            w ~ dmnorm(mu0, Q)
            Q <- tau_w * (D-(0.99*A))

            beta1 ~ dnorm(0, 0.0001)
            beta2 ~ dnorm(0, 0.0001)
            tau_w ~ dgamma(2, sigma2ref)
            tau_e ~ dgamma(2, 0.1)
          }
        ")

    }

    if (L$model == "CAR") {

      mtx <- generate_mtx(L$model, adj_mtx)
      D <- mtx$D
      n_i <- mtx$n_i

      # JAGS model code
      jags_code <- quote("
          model {
            for (i in 1:k) {
              y[i] ~ dnorm(x1[i]*beta1 + x2[i]*beta2 + w[i], tau_e)
            }

            w ~ dmnorm(mu0, Q)
            Q <- tau_w * (D - (rho*A))

            beta1 ~ dnorm(0, 0.0001)
            beta2 ~ dnorm(0, 0.0001)
            tau_w ~ dgamma(2, 1)
            tau_e ~ dgamma(2, 0.1)
            rho ~ dunif(0, 1)
          }
        ")

    }

    if (L$model=="DAGAR") {

      if (L$p_order!="none") {

        # Create order vectors (see `State_order.xlsx`, from QGIS centroids)
        order_sw <- c(21,2,17,1,10,43,40,22,25,8,27,31,24,15,29,14,48,38,45,
                      37,28,18,20,12,16,3,47,41,7,42,33,23,34,13,4,39,44,30,
                      19,26,9,6,46,36,5,35,32,11)
        order_se <- c(18,39,25,45,37,4,5,1,12,46,26,24,32,33,21,22,6,11,3,27,
                      35,20,28,44,36,43,10,8,34,15,9,40,23,30,47,16,2,7,38,
                      19,29,42,14,13,48,17,31,41)

        p_o <- case_when(
          L$p_order=="SW" ~ order_sw,
          L$p_order=="SE" ~ order_se,
          L$p_order=="NW" ~ rev(order_se),
          L$p_order=="NE" ~ rev(order_sw),
          TRUE ~ as.numeric(c(1:48))
        )

        mtx_p <- as(as.integer(p_o), "pMatrix")

        data$y <- Matrix::t(mtx_p) %*% data$y
        data$w <- Matrix::t(mtx_p) %*% data$w
        data$x1 <- Matrix::t(mtx_p) %*% data$x1
        data$x2 <- Matrix::t(mtx_p) %*% data$x2
        adj_mtx <- Matrix::t(mtx_p) %*% adj_mtx %*% mtx_p

      }

      mtx <- generate_mtx(model=L$model, adj_mtx=adj_mtx)
      n_i <- mtx$n_i

      # JAGS model code
      jags_code <- quote("
          model {
            for (i in 1:k) {
              y[i] ~ dnorm(x1[i]*beta1 + x2[i]*beta2 + w[i], tau_e)
            }

            w ~ dmnorm(mu0, Q)
            Q <- tau_w * (t(L) %*% FF %*% L)
            L <- I - B

            for (i in 1:k) {
              for (j in 1:k) {

                FF[i,j] <- ifelse(i==j, tau[i], 0)
                B[i,j] <- ifelse(
                  j<i && A[i,j]==1,
                  rho / (1+(n_i[i]-1)*(rho^2)),
                  0
                )

              }
            }

            tau <- (1+(n_i-1)*(rho^2)) / (1-rho^2)

            beta1 ~ dnorm(0, 0.0001)
            beta2 ~ dnorm(0, 0.0001)
            tau_w ~ dgamma(2, 1)
            tau_e ~ dgamma(2, 0.1)
            rho ~ dunif(0, 1)
          }
        ")

    }

    if (L$model=="DAGAR_OF") {

      # This is an approximation to the DAGAR_OF done by sampling permutations

      p1 <- c(k:1)
      p2 <- sample(1:k)
      p3 <- sample(1:k)
      mtx_p1 <- as(as.integer(p1), "pMatrix")
      mtx_p2 <- as(as.integer(p2), "pMatrix")
      mtx_p3 <- as(as.integer(p3), "pMatrix")
      mtx_p1 <- as.matrix(mtx_p1) * as.matrix(mtx_p1)
      mtx_p2 <- as.matrix(mtx_p2) * as.matrix(mtx_p2)
      mtx_p3 <- as.matrix(mtx_p3) * as.matrix(mtx_p3)
      adj_mtx_p1 <- t(mtx_p1) %*% adj_mtx %*% mtx_p1
      adj_mtx_p2 <- t(mtx_p2) %*% adj_mtx %*% mtx_p2
      adj_mtx_p3 <- t(mtx_p3) %*% adj_mtx %*% mtx_p3
      mtx_nbrs0 <- adj_mtx * lower.tri(adj_mtx)
      mtx_nbrs1 <- adj_mtx_p1 * lower.tri(adj_mtx_p1)
      mtx_nbrs2 <- adj_mtx_p2 * lower.tri(adj_mtx_p2)
      mtx_nbrs3 <- adj_mtx_p3 * lower.tri(adj_mtx_p3)
      n_i <- as.numeric(mtx_nbrs0 %*% rep(1,k))
      n_i1 <- as.numeric(mtx_nbrs1 %*% rep(1,k))
      n_i2 <- as.numeric(mtx_nbrs2 %*% rep(1,k))
      n_i3 <- as.numeric(mtx_nbrs3 %*% rep(1,k))

      # JAGS model code
      jags_code <- quote("
          model {
            for (i in 1:k) {
              y[i] ~ dnorm(x1[i]*beta1 + x2[i]*beta2 + w[i], tau_e)
            }

            w ~ dmnorm(mu0, Q)

            Q <- (Q0+Q1+Q2+Q3)/4

            Q0 <- tau_w * (t(L0) %*% FF0 %*% L0)
            Q1 <- tau_w * (mtx_p1 %*% t(L1) %*% FF1 %*% L1 %*% t(mtx_p1))
            Q2 <- tau_w * (mtx_p2 %*% t(L2) %*% FF2 %*% L2 %*% t(mtx_p2))
            Q3 <- tau_w * (mtx_p3 %*% t(L3) %*% FF3 %*% L3 %*% t(mtx_p3))

            L0 <- I - B0
            L1 <- I - B1
            L2 <- I - B2
            L3 <- I - B3

            for (i in 1:k) {
              for (j in 1:k) {

                FF0[i,j] <- ifelse(i==j, tau0[i], 0)
                FF1[i,j] <- ifelse(i==j, tau1[i], 0)
                FF2[i,j] <- ifelse(i==j, tau2[i], 0)
                FF3[i,j] <- ifelse(i==j, tau3[i], 0)

                B0[i,j] <- ifelse(
                  j<i && A[i,j]==1, rho / (1+(n_i[i]-1)*(rho^2)), 0
                )
                B1[i,j] <- ifelse(
                  j<i && adj_mtx_p1[i,j]==1, rho / (1+(n_i1[i]-1)*(rho^2)), 0
                )
                B2[i,j] <- ifelse(
                  j<i && adj_mtx_p2[i,j]==1, rho / (1+(n_i2[i]-1)*(rho^2)), 0
                )
                B3[i,j] <- ifelse(
                  j<i && adj_mtx_p3[i,j]==1, rho / (1+(n_i3[i]-1)*(rho^2)), 0
                )

              }
            }

            tau0 <- (1+(n_i-1)*(rho^2)) / (1-rho^2)
            tau1 <- (1+(n_i1-1)*(rho^2)) / (1-rho^2)
            tau2 <- (1+(n_i2-1)*(rho^2)) / (1-rho^2)
            tau3 <- (1+(n_i3-1)*(rho^2)) / (1-rho^2)

            beta1 ~ dnorm(0, 0.0001)
            beta2 ~ dnorm(0, 0.0001)
            tau_w ~ dgamma(2, 1)
            tau_e ~ dgamma(2, 0.1)
            rho ~ dunif(0, 1)
          }
        ")

    }

    objs <- c("n_i", "n_i1", "n_i2", "n_i3", "mtx_p1", "mtx_p2", "mtx_p3",
              "adj_mtx_p1", "adj_mtx_p2", "adj_mtx_p3", "sigma2ref")
    for (obj in objs) {
      if (!exists(obj)) {
        assign(obj, NULL)
      }
    }

    # Run SGLMM
    if (L$model == "SGLMM") {

      attr <- case_when(
        L$type == "path" ~ 40,
        L$type == "grid" ~ 40,
        L$type == "US" ~ 10
      )

      sm <- sparse.sglmm(
        y ~ x1 + x2 -1,
        family = "gaussian",
        data = data.frame(y=data$y, x1=data$x1, x2=data$x2),
        A = adj_mtx,
        method = "RSR",
        attractive = attr,
        minit = 5000,
        maxit = 20000
        # tune = list(sigma.s = 0.02),
        # hyper = list(
        #   a.h = 2,
        #   b.h = 10,
        #   sigma.b = 100
        # )
      )

      w_hat <- data$y - (
        # sm$coefficients[[1]] +
          (data$x1 * sm$coefficients["x1"]) +
          (data$x2 * sm$coefficients["x2"])
      ) - sm$residuals

      beta1_hat <- mean(sm$beta.sample[,1])
      beta1_se <- sd(sm$beta.sample[,1])
      beta2_hat <- mean(sm$beta.sample[,2])
      beta2_se <- sd(sm$beta.sample[,2])
      sigma2_e_hat <- mean(1/sm$tau.h.sample)
      sigma2_e_se <- sd(1/sm$tau.h.sample)

      mse <- mean((data$w-w_hat)^2)
      rho_hat <- NA
      ci_l_rho <- NA
      ci_u_rho <- NA

    }

    # Run others
    if (L$model != "SGLMM") {

      # Run JAGS model
      jm <- jags.model(
        file = textConnection(jags_code),
        data = list(
          k=k, x1=data$x1, x2=data$x2, y=data$y, D=D, I=diag(rep(1,k)),
          A=adj_mtx, n_i=n_i, sigma2ref=sigma2ref, mu0=rep(0,k), n_i1=n_i1,
          n_i2=n_i2, n_i3=n_i3, mtx_p1=mtx_p1, mtx_p2=mtx_p2, mtx_p3=mtx_p3,
          adj_mtx_p1=adj_mtx_p1, adj_mtx_p2=adj_mtx_p2, adj_mtx_p3=adj_mtx_p3
        ),
        n.chains = 1,
        n.adapt = 1000
      )
      if (L$model %in% c("CAR", "DAGAR", "DAGAR_OF")) {
        variable.names <- c("w", "rho", "beta1", "beta2", "tau_e")
      } else {
        variable.names <- c("w", "beta1", "beta2", "tau_e")
      }
      output <- coda.samples(
        model = jm,
        variable.names = variable.names,
        n.iter = 5000,
        thin = 1
      )

      # Calculate w_hat and rho_hat
      s <- summary(output)
      if (L$model %in% c("CAR", "DAGAR", "DAGAR_OF")) {

        logit_rho <- (function(x){log(x/(1-x))})(output[[1]][,"rho"])
        ci_l_logit_rho <- mean(logit_rho) - 1.96*sd(logit_rho)
        ci_u_logit_rho <- mean(logit_rho) + 1.96*sd(logit_rho)
        ci_l_rho <- 1/(1+exp(-ci_l_logit_rho))
        ci_u_rho <- 1/(1+exp(-ci_u_logit_rho))

        w_hat <- as.numeric(s$quantiles[,"50%"])[5:(k+4)]
        rho_hat <- s$statistics["rho","Mean"]

      } else {
        w_hat <- as.numeric(s$quantiles[,"50%"])[4:(k+3)]
        rho_hat <- NA
        ci_l_rho <- NA
        ci_u_rho <- NA
      }

      # Get other summary stats
      beta1_hat <- s$statistics["beta1","Mean"]
      beta1_se <- s$statistics["beta1","SD"]
      beta2_hat <- s$statistics["beta2","Mean"]
      beta2_se <- s$statistics["beta2","SD"]

      # Get sigma2_e summary stats
      tau_e_samples <- output[[1]][,"tau_e"]
      sigma2_e_hat <- mean(1/tau_e_samples)
      sigma2_e_se <- sd(1/tau_e_samples)

      # Calculate MSE
      mse <- mean((data$w-w_hat)^2)

    }

    return (list(
      "mse" = mse,
      "rho_hat" = rho_hat,
      "ci_l_rho" = ci_l_rho,
      "ci_u_rho" = ci_u_rho,
      "beta1_hat" = beta1_hat,
      "beta1_se" = beta1_se,
      "beta2_hat" = beta2_hat,
      "beta2_se" = beta2_se,
      "sigma2_e_hat" = sigma2_e_hat,
      "sigma2_e_se" = sigma2_e_se
    ))

  }

}



############################.
##### Simulation setup #####
############################.

if (run_fig345 || run_fig6 || run_fig7 || run_fig8) {

  set.seed(.tid)

  sim <- new_sim()

  sim %<>% set_config(
    num_sim = 100, # !!!!!
    parallel = "outer",
    packages = c("z.572.paper", "sp", "rgeos", "spdep", "parallel", "rgdal",
                 "simba", "ggplot2", "dplyr", "magrittr", "permute", "maps",
                 "maptools", "mvtnorm", "rjags", "MASS", "ngspatial", "Matrix",
                 "INLA")
  )

  sim %<>% add_constants(
    tau_w = 0.25,
    adj_mtx_path = generate_graph_path() %>% adj_from_graph(),
    adj_mtx_grid = generate_graph_grid() %>% adj_from_graph(),
    adj_mtx_US = adj_from_gis_us()
  )

  sim %<>% add_creator(generate_dataset)
  sim %<>% add_script(script_mainsim)

}



################################.
##### Generate figures 3-5 #####
################################.

if (run_fig345) {

  sim %<>% set_levels(
    type = c("path", "grid", "US"),
    rho = seq(from=0.1, to=0.9, by=0.1),
    model = c("iCAR", "CAR", "DAGAR", "DAGAR_OF", "SGLMM", "Scaled iCAR"),
    cov_type = "Matern",
    p_order = "none"
  )

  # Run simulation and save output
  sim %<>% run("script_mainsim", sim_uids=.tid)
  if (!file.exists("../simba.out")) { dir.create("../simba.out") }
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Create graphs
  if (run_graphs) {

    sim <- readRDS("../simba.out/sim_fig345_16200_5.31.simba")

    summ <- sim %>% summary() %>% arrange(model, type, rho)
    summ$type %<>% recode(
      "path" = "(a) Path",
      "grid" = "(b) Grid",
      "US" = "(c) USA"
    )

    summ2 <- summary(
      sim_obj = sim,
      coverage = list(
        list(name="cov_beta1", truth=1, estimate="beta1_hat", se="beta1_se"),
        list(name="cov_beta2", truth=5, estimate="beta2_hat", se="beta2_se"),
        list(name="cov_rho", truth="rho", estimate="rho_hat",
             lower="ci_l_rho", upper="ci_u_rho"),
        list(name="cov_sigma2", truth=0.4,
             estimate="sigma2_e_hat", se="sigma2_e_se")
      )
    )
    summ2$type %<>% recode(
      "path" = "(a) Path",
      "grid" = "(b) Grid",
      "US" = "(c) USA"
    )

    summ2_stacked <- sqldf("
    SELECT type, 'beta1' AS `type2`, rho, model, cov_beta1 AS `cov` FROM summ2
    UNION SELECT type, 'beta2', rho, model, cov_beta2 FROM summ2
    UNION SELECT type, 'sigma2', rho, model, cov_sigma2 FROM summ2
    UNION SELECT type, 'rho', rho, model, cov_rho FROM summ2
  ")
    summ2_stacked$type2 <- factor(summ2_stacked$type2, labels=c(
      bquote(beta[1]), bquote(beta[2]), bquote(rho), bquote(sigma^2)
    ))
    summ2_stacked$type2 <- factor(
      summ2_stacked$type2,
      levels=c("beta[1]","beta[2]","sigma^2","rho")
    )

    # Plot figure 3
    # Export 900w x 300h
    ggplot(
      data = summ,
      aes(x=rho, y=mean_mse, color=model, shape=model, linetype=model)) +
      geom_line() +
      geom_point(size=1) +
      xlim(0.1,0.9) +
      facet_wrap(~type, ncol=3) +
      labs(x="rho", y="MSE", color="Model", shape="Model", linetype="Model") +
      scale_color_manual(values=c("orangered3", "chartreuse3", "royalblue1",
                                  "plum4", "gray32", "lightsalmon3")) +
      scale_shape_manual(values=c(16,17,15,3,7,8)) +
      scale_linetype_manual(values=c(1,2,2,5,3,4)) +
      theme(legend.position = "right")


    # # !!!!! !!!!! !!!!! !!!!! !!!!!
    # # New prior
    # ggplot(
    #   data = summ %>% filter(model!="SGLMM" & model!="Scaled iCAR"),
    #   aes(x=rho, y=mean_mse, color=model, shape=model, linetype=model)) +
    #   geom_line() +
    #   geom_point(size=1) +
    #   xlim(0.1,0.9) +
    #   ylim(0,4) +
    #   facet_wrap(~type, ncol=3) +
    #   labs(x="rho", y="MSE", color="Model", shape="Model", linetype="Model") +
    #   scale_color_manual(values=c("orangered3", "chartreuse3", "royalblue1",
    #                               "plum4")) +
    #   scale_shape_manual(values=c(16,17,15,3)) +
    #   scale_linetype_manual(values=c(1,2,2,5)) +
    #   theme(legend.position = "right")
    # # !!!!! !!!!! !!!!! !!!!! !!!!!


    # Plot figure 4
    # Export 900w x 300h
    ggplot(
      data = summ %>%
        filter(model %in% c("CAR", "DAGAR", "DAGAR_OF")),
      aes(x=rho, y=mean_rho_hat)
    ) +
      geom_line(aes(color=model)) +
      geom_point(aes(color=model), size=1) +
      geom_ribbon(
        aes(ymin=mean_ci_l_rho, ymax=mean_ci_u_rho, fill=model),
        alpha = 0.2
      ) +
      geom_segment(x=0, y=0, xend=1, yend=1, linetype=5) +
      xlim(0.1,0.9) +
      ylim(0,1) +
      facet_wrap(~type, ncol=3) +
      labs(x="rho", y="Estimate", color="Model", fill="Model") +
      theme(legend.position = "right")

    # Plot figure 5
    # Export 700w x 700h
    ggplot(
      data = summ2_stacked,
      aes(x=rho, y=cov, color=model, shape=model, linetype=model)
    ) +
      geom_line() +
      geom_point(size=1) +
      xlim(0.1,0.9) +
      facet_grid(
        cols = vars(type),
        rows = vars(type2),
        labeller = labeller(.rows=label_parsed)
      ) +
      labs(x="rho",y="Coverage",color="Model",shape="Model",linetype="Model") +
      scale_color_manual(values=c("orangered3", "chartreuse3", "royalblue1",
                                  "plum4", "gray32", "lightsalmon3")) +
      scale_shape_manual(values=c(16,17,15,3,7,8)) +
      scale_linetype_manual(values=c(1,2,2,5,3,4)) +
      theme(legend.position = "right")

  }

}



#############################.
##### Generate figure 6 #####
#############################.

if (run_fig6) {

  sim %<>% set_levels(
    type = "US",
    rho = seq(from=0.1, to=0.9, by=0.1),
    model = "DAGAR",
    cov_type = "Matern",
    p_order = c("NW", "NE", "SW", "SE")
  )

  # Run simulation and save output
  sim %<>% run("script_mainsim", sim_uids=.tid)
  if (!file.exists("../simba.out")) { dir.create("../simba.out") }
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  # Create graphs
  if (run_graphs) {

    sim <- readRDS("../simba.out/sim_fig6_3600.simba")

    summ <- sim %>% summary(
      mean = list(all=TRUE, na.rm=TRUE)
    ) %>% arrange(model, type, rho)

    summ$p_order %<>% recode(
      "NW" = "Northwest",
      "NE" = "Northeast",
      "SW" = "Southwest",
      "SE" = "Southeast"
    )

    # Plot figure 6.1
    # Export 500w x 350h
    ggplot(
      data = summ %>% mutate(header="(a) USA: MSE"),
      aes(x=rho, y=mean_mse, color=p_order, shape=p_order, linetype=p_order)
    ) +
      geom_line() + geom_point(size=1) +
      facet_wrap(~header) +
      xlim(0.1,0.9) +
      labs(x="rho", y="MSE", color="Order", shape="Order", linetype="Order") +
      scale_color_manual(values=c("orangered3", "chartreuse3", "royalblue1",
                                  "plum4")) +
      scale_shape_manual(values=c(16,17,15,3)) +
      scale_linetype_manual(values=c(1,2,2,5)) +
      theme(legend.position = "right")

    # Plot figure 6.2
    # Export 500w x 350h
    ggplot(
      data = summ %>%
        mutate(header="(b) USA: Estimate and confidence bands of rho"),
        # mutate(
        #   header = "(b) USA: Estimate and confidence bands of rho",
        #   mean_ci_l_rho = pmax(mean_rho_hat-1.96*mean_rho_se,0),
        #   mean_ci_u_rho = pmin(mean_rho_hat+1.96*mean_rho_se,1)
        # ),
      aes(x=rho, y=mean_rho_hat)
    ) +
      geom_line(aes(color=p_order, linetype=p_order)) +
      geom_point(aes(color=p_order, shape=p_order), size=1) +
      geom_ribbon(
        aes(ymin=mean_ci_l_rho, ymax=mean_ci_u_rho, fill=p_order),
        alpha = 0.2
      ) +
      geom_segment(x=0, y=0, xend=1, yend=1, linetype=5) +
      facet_wrap(~header) +
      xlim(0.1,0.9) +
      ylim(0,1) +
      labs(x="rho", y="Estimate", color="Order", shape="Order",
           linetype="Order", fill="Order") +
      scale_color_manual(values=c("orangered3", "chartreuse3", "royalblue1",
                                  "plum4")) +
      scale_shape_manual(values=c(16,17,15,3)) +
      scale_linetype_manual(values=c(1,2,2,5)) +
      theme(legend.position = "right")

  }

}



#############################.
##### Generate figure 7 #####
#############################.

if (run_fig7) {

  sim %<>% set_levels(
    type = c("path", "grid", "US"),
    rho = seq(from=0.1, to=0.9, by=0.1),
    model = c("iCAR", "CAR", "DAGAR", "DAGAR_OF", "SGLMM", "Scaled iCAR"),
    cov_type = "DAGAR",
    p_order = "none"
  )

  # Run simulation and save output
  sim %<>% run("script_mainsim", sim_uids=.tid)
  if (!file.exists("../simba.out")) { dir.create("../simba.out") }
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  if (run_graphs) {

    sim <- readRDS("../simba.out/sim_fig7_9720.simba")

    summ <- sim %>% summary() %>% arrange(model, type, rho)
    summ$type %<>% recode(
      "path" = "(a) Path",
      "grid" = "(b) Grid",
      "US" = "(c) USA"
    )

    # Plot figure 7
    # Export 900w x 300h
    ggplot(
      data = summ,
      aes(x=rho, y=mean_mse, color=model, shape=model, linetype=model)
    ) +
      geom_line() +
      geom_point(size=1) +
      xlim(0.1,0.9) +
      facet_wrap(~type, ncol=3) +
      labs(x="rho", y="MSE", color="Model", shape="Model", linetype="Model") +
      scale_color_manual(values=c("orangered3", "chartreuse3", "royalblue1",
                                  "plum4", "gray32", "lightsalmon3")) +
      scale_shape_manual(values=c(16,17,15,3,7,8)) +
      scale_linetype_manual(values=c(1,2,2,5,3,4)) +
      theme(legend.position = "right")

  }

}



#############################.
##### Generate figure 8 #####
#############################.

if (run_fig8) {

  sim %<>% set_levels(
    type = c("path", "grid", "US"),
    rho = seq(from=0.1, to=0.9, by=0.1),
    model = c("iCAR", "CAR", "DAGAR", "DAGAR_OF", "SGLMM", "Scaled iCAR"),
    cov_type = "CAR",
    p_order = "none"
  )

  # Run simulation and save output
  sim %<>% run("script_mainsim", sim_uids=.tid)
  if (!file.exists("../simba.out")) { dir.create("../simba.out") }
  saveRDS(sim, file=paste0("../simba.out/sim_",.tid,".simba"))

  if (run_graphs) {

    sim <- readRDS("../simba.out/sim_fig8_9720_6.01.simba")

    summ <- sim %>% summary() %>% arrange(model, type, rho)
    summ$type %<>% recode(
      "path" = "(a) Path",
      "grid" = "(b) Grid",
      "US" = "(c) USA"
    )

    # Plot figure 8
    # Export 900w x 300h
    ggplot(
      data = summ,
      aes(x=rho, y=mean_mse, color=model, shape=model, linetype=model)
    ) +
      geom_line() +
      geom_point(size=1) +
      xlim(0.1,0.9) +
      ylim(0,4) +
      facet_wrap(~type, ncol=3) +
      labs(x="rho", y="MSE", color="Model", shape="Model", linetype="Model") +
      scale_color_manual(values=c("orangered3", "chartreuse3", "royalblue1",
                                  "plum4", "gray32", "lightsalmon3")) +
      scale_shape_manual(values=c(16,17,15,3,7,8)) +
      scale_linetype_manual(values=c(1,2,2,5,3,4)) +
      theme(legend.position = "right")

  }

}



#####################################.
##### Infant mortality analysis #####
#####################################.

if (run_ima) {

  # Sparse sGLMM: load/transform data
  data(infant)
  data(A)
  infant %<>% mutate(low_weight = low_weight/births)
  infant$id <- c(1:3071)

  # Sparse sGLMM: run model
  set.seed(12)
  fit = sparse.sglmm(
    deaths ~ low_weight + black + hispanic + gini + affluence +
             stability + offset(log(births)),
    data = infant,
    family = poisson,
    A = A,
    method = "RSR",
    tune = list(sigma.s = 0.02),
    verbose = TRUE,
    minit = 5000,
    maxit = 20000
  )
  summary(fit)
  print(fit$tau.s.est)
  print(quantile(fit$tau.s.sample, c(0.025,0.975)))

  # iCAR: run model using INLA
  set.seed(12)
  tau_prior <- list(prec=list(param=c(2, 1)))
  model <- inla(
    deaths ~ low_weight + black + hispanic + gini + affluence +
             stability + f(id, model="besag", graph=A, hyper=tau_prior),
    data = infant,
    family = "poisson",
    E = births,
    control.compute = list(dic=TRUE)
  )
  # print(summary(model))
  print(model$summary.fixed)
  print(model$summary.hyperpar)
  print(model$dic$deviance.mean)

  # # TESTING: declare model using generic INLA (CAR)
  # inla.rgeneric.CAR.model <- function(
  #   cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
  #           "log.prior", "quit"),
  #   theta = NULL) {
  #
  #   envir = parent.env(environment())
  #
  #   interpret.theta <- function() {
  #     return(list(
  #       tau_w = exp(theta[1L]),
  #       rho = 1 / (1 + exp(-theta[2L]))
  #     ))
  #   }
  #
  #   graph <- function(){
  #     require(Matrix)
  #     return(Diagonal(nrow(W), x = 1) + W)
  #   }
  #
  #   Q <- function() {
  #     require(Matrix)
  #     param <- interpret.theta()
  #     return(param$tau_w * (Diagonal(nrow(W), x = 1) - param$rho * W) )
  #   }
  #
  #   mu = function() { return(numeric(0)) }
  #
  #   log.norm.const <- function() { return(numeric(0)) }
  #
  #   log.prior <- function() {
  #
  #     param = interpret.theta()
  #     res <- dgamma(param$tau_w, 1, 5e-05, log = TRUE) + log(param$tau_w) +
  #            log(1) + log(param$rho) + log(1 - param$rho)
  #
  #     return(res)
  #
  #   }
  #
  #   initial <- function() {
  #     return(c(0,0))
  #   }
  #
  #   quit <- function() { return(invisible()) }
  #
  #   res <- do.call(match.arg(cmd), args = list())
  #
  #   return(res)
  #
  # }
  #
  # # TESTING: run model using generic INLA (CAR)
  # set.seed(12)
  # W <- as(A, "sparseMatrix")
  # e.values <- eigen(W)$values
  # rho.min <- min(e.values)
  # rho.max <- max(e.values)
  # W <- W / rho.max
  # CAR.model <- inla.rgeneric.define(inla.rgeneric.CAR.model, W=W)
  # model <- inla(
  #   deaths ~ low_weight + black + hispanic + gini + affluence +
  #            stability + f(id, model=CAR.model),
  #   data = infant,
  #   family = "poisson",
  #   E = births,
  #   control.compute = list(dic=TRUE)
  # )
  # print(summary(model))
  # print(model$dic$deviance.mean)

  # DAGAR: declare model using INLA
  inla.rgeneric.DAGAR.model <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const",
            "log.prior", "quit"),
    theta = NULL) {

    envir = parent.env(environment())

    interpret.theta <- function() {
      return(list(
        tau_w = exp(theta[1L]),
        rho = 1 / (1 + exp(-theta[2L]))
      ))
    }

    graph <- function() {

      require(Matrix)

      # Get number of neighbors
      dim <- length(W[1,])
      mtx_nbrs <- W * lower.tri(W)
      n_i <- as.numeric(mtx_nbrs %*% rep(1,dim))

      # Create matrices
      b_i <- 0.5 / (1+((n_i-1)*0.25))
      B <- mtx_nbrs * matrix(rep(b_i,dim), nrow=dim)
      L <- diag(rep(1,dim)) - B
      G <- t(L) %*% L

      return(as(G, "sparseMatrix"))

    }

    Q <- function() {

      require(Matrix)
      param <- interpret.theta()

      # Get number of neighbors
      dim <- length(W[1,])
      mtx_nbrs <- W * lower.tri(W)
      n_i <- as.numeric(mtx_nbrs %*% rep(1,dim))

      # Create matrices
      tau <- (1+(n_i-1)*(param$rho^2)) / (1-(param$rho^2))
      FF <- diag(tau)
      b_i <- param$rho / (1+((n_i-1)*(param$rho^2)))
      B <- mtx_nbrs * matrix(rep(b_i,dim), nrow=dim)
      L <- diag(rep(1,dim)) - B
      Q <- param$tau_w * ( t(L) %*% FF %*% L )

      return(as(Q, "sparseMatrix"))

    }

    mu = function() { return(numeric(0)) }

    log.norm.const <- function() { return(numeric(0)) }

    log.prior <- function() {

      param = interpret.theta()
      # res <- dgamma(param$tau_w, 2, 0.01, log=TRUE) + log(param$tau_w) +
      #        log(1) + log(param$rho) + log(1-param$rho)
      res <- dgamma(param$tau_w, 2, 1, log=TRUE) + log(param$tau_w) +
             log(1) + log(param$rho) + log(1-param$rho)

      return(res)

    }

    initial <- function() {
      return(c(0,0))
    }

    quit <- function() { return(invisible()) }

    res <- do.call(match.arg(cmd), args = list())

    return(res)

  }

  # DAGAR: run model using INLA
  set.seed(13)
  W <- as(A, "sparseMatrix")
  DAGAR.model <- inla.rgeneric.define(inla.rgeneric.DAGAR.model, W=W)
  model <- inla(
    deaths ~ low_weight + black + hispanic + gini + affluence +
             stability + f(id, model=DAGAR.model),
    data = infant,
    family = "poisson",
    E = births,
    control.compute = list(dic=TRUE)
  )
  summary(model)
  print(model$summary.fixed)
  print(model$summary.hyperpar)
  print(model$dic$deviance.mean)

}



#########################.
##### Testing: INLA #####
#########################.

if (FALSE) {

  set.seed(12)

  # Generate data
  data <- generate_dataset(
    type = "grid",
    tau_w = 0.25,
    phi = -1*log(0.7)
  )
  df <- data.frame(
    id = c(1:100),
    y = data$y,
    x1 = data$x1,
    x2 = data$x2
  )
  adj_mtx_grid <- generate_graph_grid() %>% adj_from_graph()

  # Run INLA CAR model
  model <- inla(
    y ~ x1 + x2 + f(id, model="besagproper", graph=adj_mtx_grid),
    data = df
  )
  summary(model)

  # Create rgeneric() ICAR model
  inla.rgeneric.CAR.model <- function(
    cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL) {

    envir = parent.env(environment())

    interpret.theta <- function() {
      return(list(
        tau_w = exp(theta[1L]),
        rho = 1 / (1 + exp(-theta[2L]))
      ))
    }

    graph <- function(){

      # !!!!! modify
      require(Matrix)
      return(Diagonal(nrow(W), x = 1) + W)

    }

    Q <- function() {

      # !!!!! modify
      require(Matrix)
      param <- interpret.theta()
      return(param$tau_w * (Diagonal(nrow(W), x = 1) - param$rho * W) )

    }

    mu = function() { return(numeric(0)) }

    log.norm.const <- function() { return(numeric(0)) }

    log.prior <- function() {

      param = interpret.theta()
      res <- dgamma(param$tau_w, 1, 5e-05, log = TRUE) + log(param$tau_w) +
             log(1) + log(param$rho) + log(1 - param$rho)

      return(res)

    }

    initial <- function() {
      return(c(0,0))
    }

    quit <- function() { return(invisible()) }

    res <- do.call(match.arg(cmd), args = list())

    return(res)

  }
  W <- as(adj_mtx_grid, "sparseMatrix")
  e.values <- eigen(W)$values
  rho.min <- min(e.values)
  rho.max <- max(e.values)
  W <- W / rho.max

  CAR.model <- inla.rgeneric.define(inla.rgeneric.CAR.model, W=W)

  # Run CAR model
  set.seed(12)
  model <- inla(
    y ~ x1 + x2 + f(id, model=CAR.model),
    data = df
  )
  summary(model)

}



####################################.
##### Testing: adj mtx for PPT #####
####################################.

if (FALSE) {

  # Adj matrix for presentation
  A <- rbind(
    c(0,1,1,0,0,0,0,0),
    c(1,0,1,1,0,0,0,0),
    c(1,1,0,1,0,0,0,1),
    c(0,1,1,0,1,0,0,1),
    c(0,0,0,1,0,1,0,1),
    c(0,0,0,0,1,0,1,1),
    c(0,0,0,0,0,1,0,1),
    c(0,0,1,1,1,1,1,0)
  )
  mtx_dagar <- generate_mtx(model="DAGAR", adj_mtx=A, rho=0.8)
  cov_dagar <- solve(mtx_dagar$Q)
  mtx_car <- generate_mtx(model="CAR", adj_mtx=A, rho=0.8)
  cov_car <- solve(mtx_car$Q)

  round(cov_dagar/mean(diag(cov_dagar)),2)
  round(cov_car/mean(diag(cov_car)),2)

}



###################.
##### Archive #####
###################.

if (FALSE) {

  # iCAR: Generate D matrix
  set.seed(12)
  mtx <- generate_mtx("iCAR", A)
  D <- mtx$D
  # D2 <- matrix(as.numeric(D), ncol=3071) # !!!!!

  # iCAR: JAGS model code
  jags_code <- quote("
    model {
      for (i in 1:k) {
        y[i] ~ dpois(
          log(births[i]) + alpha + beta1*low_weight[i] + beta2*black[i] +
          beta3*hispanic[i] + beta4*gini[i] + beta5*affluence[i] +
          beta6*stability[i] + w[i]
        )
      }

      w ~ dmnorm(mu0, Q)
      Q <- tau_w * ((D-A)+(0.001*I))

      alpha ~ dnorm(0, 0.000001)
      beta1 ~ dnorm(0, 0.000001)
      beta2 ~ dnorm(0, 0.000001)
      beta3 ~ dnorm(0, 0.000001)
      beta4 ~ dnorm(0, 0.000001)
      beta5 ~ dnorm(0, 0.000001)
      beta6 ~ dnorm(0, 0.000001)

      tau_w ~ dgamma(2, 1)
    }
  ")

  # iCAR: Generate D matrix
  # Run JAGS model
  system.time({
    jm <- jags.model(
      file = textConnection(jags_code),
      data = list(
        y = infant$deaths,
        births = infant$births,
        low_weight = infant$low_weight,
        black = infant$black,
        hispanic = infant$hispanic,
        gini = infant$gini,
        affluence = infant$affluence,
        stability = infant$stability,
        A = A,
        k = 3071,
        I = diag(rep(1,3071)),
        mu0 = rep(0,3071),
        D = D
      ),
      n.chains = 1,
      n.adapt = 10
      # n.adapt = 1000
    )
    output <- coda.samples(
      model = jm,
      variable.names = c("alpha", "beta1", "beta2", "beta3", "beta4", "beta5",
                         "beta6", "tau_w"),
      n.iter = 10,
      # n.iter = 1000,
      thin = 1
    )
    summary(output)
  })

}
