# jonashaslbeck@protonmail.com; May 16th, 2023

# --------------------------------------------------------
# ------- What is happening here? ------------------------
# --------------------------------------------------------

# These are auxiliary functions for Simulation 2

# --------------------------------------------------------
# ------- Data Generation --------------------------------
# --------------------------------------------------------

# VARY:
# 1) Number of total subjects (50, 175, 300)
# 2) Percentage in first group (0.5, 0.2)
# 3) Sizes of differences: c(0.05, 0.15, 0.25)
# 4) number of variables: c(6, 12)

SimData <- function(p,
                    Nsubj,
                    prop,
                    Delta,
                    Nt = 100,
                    seed = 1,
                    pbar = FALSE,
                    var_burnin = 500) {

  # ----- Fixed Parameters -----
  # Number of time points per subject
  # Nt <- 100
  # Random effects variances (for: phi, theta, random intercepts [gamma])
  sigma <- 0.05 # keep constant for now
  # sigma_btw <- 0.10

  # Reproducibility
  set.seed(seed)

  # ----- Define Group 1 -----
  phi_diag <- 0.4
  phi_offdiag <- 0.6/p

  ## phi-matrix
  phi1 <- matrix(NA, p, p)
  # diagonal blocks: within-valence
  phi1[1:(p/2), 1:(p/2)] <- phi_offdiag
  phi1[(p/2+1):p, (p/2+1):p] <- phi_offdiag
  # off-diagonal blocks: between valence
  phi1[1:(p/2), (p/2+1):p] <- -phi_offdiag
  phi1[(p/2+1):p, 1:(p/2)] <- -phi_offdiag
  diag(phi1) <- phi_diag

  ## Random effect variances of phi
  phi1_RE_sd <- matrix(sigma, p, p)

  ## residual partial correlations
  theta1 <- matrix(NA, p, p)
  theta_offdiag <- 0.6/p
  # diagonal blocks: within-valence
  theta1[1:(p/2), 1:(p/2)] <- phi_offdiag
  theta1[(p/2+1):p, (p/2+1):p] <- phi_offdiag
  # off-diagonal blocks: between valence
  theta1[1:(p/2), (p/2+1):p] <- -phi_offdiag
  theta1[(p/2+1):p, 1:(p/2)] <- -phi_offdiag
  diag(theta1) <- 1

  ## Random effects of theta
  theta1_RE_sd <- matrix(sigma, p, p)

  ## Partial correlations across intercepts
  gamma1 <- matrix(0, p, p) # all zero
  gamma_offdiag <- 0.6/p
  gamma1[1:(p/2), 1:(p/2)] <- gamma_offdiag
  gamma1[(p/2+1):p, (p/2+1):p] <- gamma_offdiag
  # off-diagonal blocks: between valence
  gamma1[1:(p/2), (p/2+1):p] <- -gamma_offdiag
  gamma1[(p/2+1):p, 1:(p/2)] <- -gamma_offdiag
  diag(gamma1) <- 1

  # ----- Define Group differences -----

  # Run while, in case group differences stack up in a way that renders a VAR unstable / a Gaussian rank-deficient
  check_base <- FALSE
  counter <- 0
  while(!check_base) {

    ## Phi / Lagged effects
    # number of group differences
    n_GD_phi <- ceiling(p^2*0.15)
    # Fixed Effects
    v_GD_phi <- sample(c(rep(1,n_GD_phi), rep(0, p^2-n_GD_phi)), size=p^2, replace=FALSE)
    v_GD_phi_signs <- sample(c(-1, 1), size=p^2, replace=T)
    m_GD_phi <- matrix(v_GD_phi*v_GD_phi_signs*Delta, p, p)
    # RE Variances
    v_GD_phi_RE <- sample(c(rep(1,n_GD_phi), rep(0, p^2-n_GD_phi)), size=p^2, replace=FALSE)
    m_GD_phi_RE <- matrix(v_GD_phi_RE*Delta, p, p)

    ## Residual Partial correlations
    # Number of group differences
    n_GD_theta <- ceiling(0.15*p*(p-1)/2)
    # Fixed Effects
    v_GD_theta <- sample(c(rep(1,n_GD_theta), rep(0, p*(p-1)/2-n_GD_theta)), size=p*(p-1)/2, replace=FALSE)
    v_GD_theta_signs <- sample(c(-1, 1), size=p*(p-1)/2, replace=T)
    m_GD_theta <- matrix(0, p, p)
    m_GD_theta[upper.tri(m_GD_theta)] <- v_GD_theta*v_GD_theta_signs*Delta
    m_GD_theta <- m_GD_theta + t(m_GD_theta)
    # RE Variances
    v_GD_theta_RE <- sample(c(rep(1,n_GD_theta), rep(0, p*(p-1)/2-n_GD_theta)), size=p*(p-1)/2, replace=FALSE)
    m_GD_theta_RE <- matrix(0, p, p)
    m_GD_theta_RE[upper.tri(m_GD_theta_RE)] <- v_GD_theta_RE*Delta
    m_GD_theta_RE <- m_GD_theta_RE + t(m_GD_theta_RE)

    ## Partial correlations across random intercepts;
    #  (translates into "mean between network" (since all predictors are centered))
    v_GD_gam <- sample(c(rep(1,n_GD_theta), rep(0, p*(p-1)/2-n_GD_theta)), size=p*(p-1)/2, replace=FALSE)
    v_GD_gam_signs <- sample(c(-1, 1), size=p*(p-1)/2, replace=T)
    m_GD_gam <- matrix(0, p, p)
    m_GD_gam[upper.tri(m_GD_gam)] <- v_GD_gam_signs*v_GD_gam*Delta
    m_GD_gam <- m_GD_gam + t(m_GD_gam)


    # ----- Defining Group 2 -----
    ## phi
    phi2 <- phi1 + m_GD_phi # fixed effects
    phi2_RE_sd <- phi1_RE_sd + m_GD_phi_RE # random effects SD
    ## res pcor
    theta2 <- theta1 + m_GD_theta # fixed effects
    theta2_RE_sd <- theta1_RE_sd + m_GD_theta_RE # random effects SD
    ## between pcor (pcor of random intercepts)
    gamma2 <- gamma1 + m_GD_gam

    # ----- Check whether both groups have admissible distributions -----
    # Check Group 1
    G1_phi1_ch <- all(abs(eigen(phi1)$values) < 1)
    G1_theta1_ch <- all(eigen(theta1)$values > 0)
    G1_gam1_ch <- all(eigen(gamma1)$values > 0)
    # Check Group 2
    G2_phi1_ch <- all(abs(eigen(phi2)$values) < 1)
    G2_theta2_ch <- all(Re(eigen(theta2)$values) > 0)
    G2_gam1_ch <- all(Re(eigen(gamma2)$values) > 0)

    # Combine & evaluate checks
    v_checks <- c(G1_phi1_ch, G1_gam1_ch,
                  G1_theta1_ch, G2_theta2_ch,
                  G2_phi1_ch, G2_gam1_ch)
    if(all(v_checks))  {
      check_base <- TRUE
    } else {
      counter <- counter + 1
      if(counter > 1000) stop("Unable to generate Base model in data generation!")
    }

  } # end while: check base model


  # ----- Save True Models -----

  # Combine
  l_phi <- list(phi1, phi2)
  l_phi_RE_sd <- list(phi1_RE_sd, phi2_RE_sd)
  l_theta <- list(theta1, theta2)
  l_theta_RE_sd <- list(theta1_RE_sd, theta2_RE_sd)
  l_gamma <- list(gamma1, gamma2)

  # ----- Generate Data -----

  # Storage
  l_data <- list(list(),
                 list())

  # Get proportions for both groups
  v_props <- c(prop, 1-prop)

  # Progress bar
  if(pbar) pb <- txtProgressBar(0, Nsubj, style = 3)
  counter_dg <- 0

  for(j in 1:2) {  # groups

    N_j <- Nsubj*v_props[j] # number if subjects sampled in group j

    for(i in 1:N_j) { # individuals

      counter_ij <- 0
      check_ij <- FALSE
      while(!check_ij) {

        ## Draw VAR parameters
        RE_phi_i <- matrix(NA, p, p)
        for(s1 in 1:p) for(s2 in 1:p) RE_phi_i[s1, s2] <- rnorm(1, mean=0, sd=l_phi_RE_sd[[j]][s1, s2])
        phi_i <- l_phi[[j]] + RE_phi_i # add random effect
        check_ij_1 <- all(abs(eigen(phi_i)$value) < 1)

        ## Draw Residual pcors
        noise_theta <- matrix(0, p, p)
        for(s1 in 1:p) for(s2 in 1:p) noise_theta[s1, s2] <- rnorm(1, mean=0, sd=l_theta_RE_sd[[j]][s1, s2]) # this could be cleaning; I'm drawing >twice is much as I need
        noise_theta[lower.tri(noise_theta)] <- 0
        diag(noise_theta) <- 0
        noise_theta <- noise_theta + t(noise_theta)
        theta_i <- l_theta[[j]] + noise_theta
        check_ij_2 <- all(Re(eigen(theta_i)$value) > 0)

        # Convert partial cor matrix to cor matrix -> used below
        theta_i_cor <- pcor2cor(theta_i)
        check_ij_3 <- all(!is.na(theta_i_cor))

        if(all(c(check_ij_1, check_ij_2, check_ij_3))) {
          check_ij <- TRUE
        } else {
          counter_ij <- counter_ij + 1
          if(counter_ij > 1000) stop("Unable to generate REs that give stable/posdef matrices")
        }

      } # end while: check stable/pos-def

      ## Draw intercepts
      gamma_j_cors <- pcor2cor(l_gamma[[j]]) #* sigma^2 # convert to cor + scale to RE variance 0.05
      gamma_j_var <- gamma_j_cors*sigma # same random effects variance as for other parameters
      v_ints <- MASS::mvrnorm(n=1, mu=rep(0, p), Sigma = gamma_j_var)

      # Sample data
      l_data[[j]][[i]] <- simulateVAR(
        # pars = diag(p)*0.5,
        pars = phi_i,
        means = v_ints, # those are actually intercepts
        # means = rep(0, p),
        residuals = theta_i_cor,
        # residuals = diag(p),
        Nt = Nt,
        burnin = var_burnin)

      counter_dg <- counter_dg + 1
      if(pbar) setTxtProgressBar(pb, counter_dg)

    } # for: i

  } # for j

  data1 <- cbind(do.call(rbind, l_data[[1]]), rep(1:(Nsubj*v_props[1]), each=Nt))
  colnames(data1)[p+1] <- "id"
  data2 <- cbind(do.call(rbind, l_data[[2]]), rep((1:(Nsubj*v_props[2]))+1000, each=Nt))
  colnames(data2)[p+1] <- "id"

  # ----- Save -----

  outlist <- list("call" = list("p" = p,
                                "Nsubj" = Nsubj,
                                "prop" = prop,
                                "Delta" = Delta),
                  "true" = list("basemodel" = list("l_phi" = l_phi,
                                                   "phi2_RE_sd" = phi2_RE_sd,
                                                   "l_theta" = l_theta,
                                                   "l_gamma" = l_gamma),
                                "GroupDiffs" = list("m_GD_phi" = m_GD_phi,
                                                   "m_GD_phi_RE" = m_GD_phi_RE,
                                                   "m_GD_theta" = m_GD_theta,
                                                   "m_GD_theta_RE" = m_GD_theta_RE,
                                                   "m_GD_gam" = m_GD_gam)),
                  "data" = list(data1, data2),
                  "dataSamp"=list("counter_base"=as.numeric(counter),
                                  "counter_datasamp"=as.numeric(check_ij)))

  return(outlist)

} # eoF


