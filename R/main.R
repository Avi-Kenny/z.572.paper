# Title: "572 Paper Reproduction"
# Author: Avi Kenny
# Date: 2020-05-20



#################.
##### Setup #####
#################.

# Load packages
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

# Load functions
source(adj_from_gis_us.R)
source(adj_from_graph.R)
source(generate_dataset.R)
source(generate_graph_grid.R)
source(generate_graph_path.R)
source(generate_matern_cov.R)
source(generate_mtx.R)

# Set code blocks to run
run_fig1 <- FALSE
run_fig2 <- FALSE
run_fig35 <- TRUE
run_fig68 <- FALSE
run_ima <- FALSE



#########################.
##### Testing: INLA #####
#########################.

if (FALSE) {

  # INLA:::inla.dynload.workaround()

  # beta <- c(2,3)
  # x <- rnorm(1000, mean=6,sd=2)
  # y <- rnorm(1000, mean=(beta[1]+beta[2]*x),sd=1)
  # model_1 <- inla(
  #   y~x,
  #   family = c("gaussian"),
  #   data = data.frame(x=x, y=y),
  #   control.predictor = list(link=1),
  #   verbose = TRUE
  # )
  # summary(model_1)

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

      # Generate average neighbor-pair correlation for DAGAR
      mtx_dagar <- generate_mtx(model="DAGAR", adj_mtx=adj_mtx, rho=rho)
      cov_dagar <- solve(mtx_dagar$Q)

      cor <- 0
      count <- 0
      for (i in 2:length(mtx_dagar$neighbors)) {

        nbs <- mtx_dagar$neighbors[[i]]
        if (nbs[1]!=0) {
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

    }

  }

  # Plot results
  ggplot(
    data = data.frame(
      x = rep(rep(rhos, each=2), 3),
      y = rho_out,
      grp = rep(c("DAGAR","CAR"), 27),
      type = rep(c(
        "(a) Path graph of length 100",
        "(b) 10 x 10 grid",
        "(c) 48 contiguous US states"
      ), each=18)
    ),
    aes(x=x, y=y, group=grp,
        color=as.factor(grp))) +
    geom_segment(x=0, y=0, xend=1, yend=1, color="grey") +
    geom_line() + geom_point(size=3) +
    scale_color_manual(values=c("firebrick", "darkolivegreen4")) +
    xlim(0,1) +
    ylim(0,1) +
    facet_wrap(~type, ncol=3) +
    labs(
      title = paste0("Figure 1. Average neighbor pair correlations as a function",
                     " of rho for proper CAR and DAGAR model"),
      x = "rho",
      y = "Average neighbor correlation",
      color = "Model"
    ) +
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
  ggplot(
    data = data.frame(
      x = rep(rhos, 2),
      y = c( rd_path(rhos), rd_grid(rhos) ),
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
    # ylim(0,1) +
    facet_wrap(~grp, ncol=2) +
    labs(
      title = paste0("Figure 2. Asymptotic relative difference between the DAGAR",
                     "model and the order free DAGAR model"),
      x = "rho",
      y = "Relative difference"
    )

}



################################.
##### Generate figures 3-5 #####
################################.

if (run_fig35) {

  # # !!!!! Testing ("resume" after `# Generate data`)
  # type <- "path"
  # rho <- 0.5
  # model <- "SGLMM"
  # tau_w = 0.25
  # adj_mtx_path <- generate_graph_path() %>% adj_from_graph()
  # adj_mtx_grid <- generate_graph_grid() %>% adj_from_graph()
  # adj_mtx_US <- adj_from_gis_us()
  # adj_mtx <- adj_mtx_path
  # k <- 100
  # data <- generate_dataset(
  #   type = type,
  #   tau_w = tau_w,
  #   phi = -1*log(rho)
  # )

  set.seed(12)

  start_time <- Sys.time()

  sim <- new_sim()

  sim %<>% set_config(
    num_sim = 20, # !!!!!
    parallel = "outer",
    packages = c("z.572.paper", "sp", "rgeos", "spdep", "parallel", "rgdal",
                 "simba", "ggplot2", "dplyr", "magrittr", "permute", "maps",
                 "maptools", "mvtnorm", "rjags", "MASS", "ngspatial", "Matrix")
  )

  sim %<>% add_constants(
    tau_w = 0.25,
    adj_mtx_path = generate_graph_path() %>% adj_from_graph(),
    adj_mtx_grid = generate_graph_grid() %>% adj_from_graph(),
    adj_mtx_US = adj_from_gis_us()
  )

  sim %<>% set_levels(
    type = c("path", "grid", "US"),
    rho = c(0.1, 0.3, 0.5, 0.7, 0.9),
    # rho = seq(from=0.1, to=0.9, by=0.1),
    model = c("iCAR", "CAR", "DAGAR", "DAGAR_OF", "SGLMM", "Scaled iCAR")
    # type = "US",
    # rho = c(0.2,0.8),
    # model = c("iCAR", "CAR", "DAGAR", "DAGAR_OF", "SGLMM", "Scaled iCAR")
  )

  sim %<>% add_creator(generate_dataset)

  sim %<>% add_script(
    "script_mainsim",
    function(L,C) {

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
        phi = -1*log(L$rho)
      )

      # Create JAGS code and objects

      if (L$model == "iCAR") {

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
            Q <- tau_w * (D-(0.999999*A))

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

        # JAGS model code
        jags_code <- quote("
          model {
            for (i in 1:k) {
              y[i] ~ dnorm(x1[i]*beta1 + x2[i]*beta2 + w[i], tau_e)
            }

            w ~ dmnorm(mu0, Q)
            Q <- tau_w * (D-(0.999999*A))

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

        mtx <- generate_mtx(L$model, adj_mtx)
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
        adj_mtx_p1 <- mtx_p1 %*% adj_mtx %*% t(mtx_p1)
        adj_mtx_p2 <- mtx_p2 %*% adj_mtx %*% t(mtx_p2)
        adj_mtx_p3 <- mtx_p3 %*% adj_mtx %*% t(mtx_p3)
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
            Q1 <- tau_w * (t(mtx_p1) %*% t(L1) %*% FF1 %*% L1 %*% mtx_p1)
            Q2 <- tau_w * (t(mtx_p2) %*% t(L2) %*% FF2 %*% L2 %*% mtx_p2)
            Q3 <- tau_w * (t(mtx_p3) %*% t(L3) %*% FF3 %*% L3 %*% mtx_p3)

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
          L$type == "US" ~ 15
        )

        sm <- sparse.sglmm(
          y ~ x1 + x2,
          data = data.frame(y=data$y, x1=data$x1, x2=data$x2),
          A = adj_mtx,
          method = "RSR",
          attractive = attr,
          minit = 4000,
          maxit = 8000
        )

        w_hat <- data$y - (
          sm$coefficients[[1]] +
          (data$x1 * sm$coefficients["x1"]) +
          (data$x2 * sm$coefficients["x2"])
        ) - sm$residuals

        beta1_hat <- mean(sm$beta.sample[,2])
        beta1_se <- sd(sm$beta.sample[,2])
        beta2_hat <- mean(sm$beta.sample[,3])
        beta2_se <- sd(sm$beta.sample[,3])
        sigma2_e_hat <- mean(sm$residuals^2) # !!!!! Placeholder
        sigma2_e_se <- sd(sm$residuals^2) # !!!!! Placeholder

        mse <- mean((data$w-w_hat)^2)
        rho_hat <- NA
        rho_se <- NA

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
          n.iter = 1000,
          thin = 1
        )

        # Calculate w_hat
        s <- summary(output)
        if (L$model %in% c("CAR", "DAGAR", "DAGAR_OF")) {

          w_hat <- as.numeric(s$quantiles[,"50%"])[5:(k+4)]
          rho_hat <- s$statistics["rho","Mean"]
          rho_se <- s$statistics["rho","SD"]

        } else {
          w_hat <- as.numeric(s$quantiles[,"50%"])[4:(k+3)]
          rho_hat <- NA
          rho_se <- NA
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
        "rho_se" = rho_se,
        "beta1_hat" = beta1_hat,
        "beta1_se" = beta1_se,
        "beta2_hat" = beta2_hat,
        "beta2_se" = beta2_se,
        "sigma2_e_hat" = sigma2_e_hat,
        "sigma2_e_se" = sigma2_e_se
      ))

    }
  )

  sim %<>% run("script_mainsim")

  saveRDS(sim, file=paste("sim",Sys.time()))
  # sim <- readRDS("../sim 2020-05-16 18_02_37")

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
      list(name="cov_rho", truth="rho", estimate="rho_hat", se="rho_se"),
      list(name="cov_sigma2", truth=(1/2.5),
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

  # Plot figure 3
  ggplot(
    data = summ,
    aes(x=rho, y=mean_mse, color=model)) +
    geom_line() + geom_point(size=1) +
    xlim(0,1) +
    # ylim(0,4) +
    facet_wrap(~type, ncol=3) +
    labs(
      title = paste0("Figure 3. MSE as a function of the true rho for the ",
                     "exponential GP simulation data"),
      x = "rho",
      y = "MSE",
      color = "Model"
    ) +
    theme(legend.position = "right")

  # Plot figure 4
  ggplot(
    data = summ %>%
      filter(model %in% c("CAR", "DAGAR", "DAGAR_OF")) %>%
      mutate(
        # Think about calculating SEs of expit(rho)
        ci_l = pmax(mean_rho_hat-1.96*mean_rho_se,0),
        ci_u = pmin(mean_rho_hat+1.96*mean_rho_se,1)
      ),
    aes(x=rho, y=mean_rho_hat)) +
    geom_line(aes(color=model)) +
    geom_point(aes(color=model), size=1) +
    geom_ribbon(
      aes(ymin=ci_l, ymax=ci_u, fill=model),
      alpha = 0.2
    ) +
    geom_segment(x=0, y=0, xend=1, yend=1, linetype=5) +
    xlim(0,1) +
    ylim(0,1) +
    facet_wrap(~type, ncol=3) +
    labs(
      title = paste0("Figure 4. Estimate and confidence bands of rho as a ",
                     "function of the true rho"),
      x = "rho",
      y = "Estimate",
      color = "Model",
      fill = "Model"
    ) +
    theme(legend.position = "right")

  # Plot figure 5
  ggplot(
    data = summ2_stacked,
    aes(x=rho, y=cov, color=model)) +
    geom_line() + geom_point(size=1) +
    xlim(0,1) +
    # ylim(0,4) +
    facet_grid(
      cols = vars(type),
      rows = vars(type2),
      labeller = labeller(.rows=label_parsed)
    ) +
    labs(
      title = paste0("Figure 5. Coverage probabilities of the parameters as a",
                     " function of the true rho"),
      x = "rho",
      y = "Coverage",
      color = "Model"
    ) +
    theme(legend.position = "right")

  print(round(difftime(Sys.time(), start_time),2))

}



################################.
##### Generate figures 6-8 #####
################################.

if (run_fig68) {

  # Create order vectors (see `State_order.xlsx`, from QGIS centroids)
  order_sw <- c(21,2,17,1,10,43,40,22,25,8,27,31,24,15,29,14,48,38,45,
                37,28,18,20,12,16,3,47,41,7,42,33,23,34,13,4,39,44,30,
                19,26,9,6,46,36,5,35,32,11)
  order_se <- c(18,39,25,45,37,4,5,1,12,46,26,24,32,33,21,22,6,11,3,27,
                35,20,28,44,36,43,10,8,34,15,9,40,23,30,47,16,2,7,38,
                19,29,42,14,13,48,17,31,41)
  order_nw <- rev(order_se)
  order_ne <- rev(order_sw)

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
  # set.seed(12)
  # fit = sparse.sglmm(
  #   deaths ~ low_weight + black + hispanic + gini + affluence +
  #            stability + offset(log(births)),
  #   data = infant,
  #   family = poisson,
  #   A = A,
  #   method = "RSR",
  #   tune = list(sigma.s = 0.02),
  #   verbose = TRUE,
  #   minit = 4000,
  #   maxit = 8000
  # )
  # summary(fit)

  # iCAR: run model using INLA
  set.seed(12)
  model <- inla(
    deaths ~ low_weight + black + hispanic + gini + affluence +
             # stability + f(id, model="besagproper", graph=A),
             stability + f(id, model="besag", graph=A),
    data = infant,
    family = "poisson",
    E = births,
    control.compute = list(dic=TRUE)
    # verbose = TRUE
  )
  print(summary(model))
  print(model$dic$deviance.mean)

  # TESTING: declare model using generic INLA (CAR)
  inla.rgeneric.CAR.model <- function(
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

    graph <- function(){
      require(Matrix)
      return(Diagonal(nrow(W), x = 1) + W)
    }

    Q <- function() {
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

  # TESTING: run model using generic INLA (CAR)
  set.seed(12)
  W <- as(A, "sparseMatrix")
  e.values <- eigen(W)$values
  rho.min <- min(e.values)
  rho.max <- max(e.values)
  W <- W / rho.max
  CAR.model <- inla.rgeneric.define(inla.rgeneric.CAR.model, W=W)
  model <- inla(
    deaths ~ low_weight + black + hispanic + gini + affluence +
             stability + f(id, model=CAR.model),
    data = infant,
    family = "poisson",
    E = births,
    control.compute = list(dic=TRUE)
  )
  print(summary(model))
  print(model$dic$deviance.mean)
  # started 10:01 AM

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

    graph <- function(){

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
      res <- dgamma(param$tau_w, 2, 1, log = TRUE) + log(param$tau_w) +
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

  # DAGAR: run model using INLA
  set.seed(12)
  W <- as(A, "sparseMatrix")
  # e.values <- eigen(W)$values
  # rho.min <- min(e.values)
  # rho.max <- max(e.values)
  # W <- W / rho.max
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
