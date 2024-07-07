# jonashaslbeck@protonmail.com; April 8th, 2024

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
                    n_GD=NULL,
                    Nt = 100,
                    seed = 1,
                    pbar = FALSE,
                    maxresamp = 10000,
                    DeltaPars = c("phi", "phi_RE", "theta", "theta_RE", "gamma"),
                    Signs = c(-1, 1),
                    var_burnin = 500) {

  # ----- Fixed Parameters -----
  # Number of time points per subject
  # Nt <- 100
  # Random effects variances (for: phi, theta, random intercepts [gamma])
  sigma <- 0.05 # keep constant for now; this is a standard deviation
  sigma_int <- 0.5 # for intercepts: much bigger

  ## Phi eigen cutoff;
  # we require max abs eigenvalue to be this value;
  # we chose a value smaller than 1 to avoid near-instability
  phi_eig_c <- 0.95

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

  ## Partial correlations across intercepts (between person effects)
  gamma1 <- matrix(0, p, p) # all zero
  gamma_offdiag <- 0.6/p
  gamma1[1:(p/2), 1:(p/2)] <- gamma_offdiag
  gamma1[(p/2+1):p, (p/2+1):p] <- gamma_offdiag
  # off-diagonal blocks: between valence
  gamma1[1:(p/2), (p/2+1):p] <- - gamma_offdiag
  gamma1[(p/2+1):p, 1:(p/2)] <- - gamma_offdiag
  diag(gamma1) <- 1

  # ----- Define Group differences -----

  v_signs <- Signs

  # Run while, in case group differences stack up in a way that renders a VAR unstable / a Gaussian rank-deficient
  check_base <- FALSE
  counter <- 0
  while(!check_base) {

    ## Phi / Lagged effects
    # number of group differences (non-symmetric matrices)
    if(is.null(n_GD)) n_GD_phi <- ceiling(p^2*0.15) else n_GD_phi <- 1
    # Fixed Effects
    v_GD_phi <- sample(c(rep(1,n_GD_phi), rep(0, p^2-n_GD_phi)), size=p^2, replace=FALSE)
    v_GD_phi_signs <- sample(v_signs, size=p^2, replace=T)
    m_GD_phi <- matrix(v_GD_phi*v_GD_phi_signs*Delta, p, p)
    # RE Variances
    v_GD_phi_RE <- sample(c(rep(1,n_GD_phi), rep(0, p^2-n_GD_phi)), size=p^2, replace=FALSE)
    m_GD_phi_RE <- matrix(v_GD_phi_RE*Delta, p, p) # here we always add to avoid negative variances

    ## Residual Partial correlations
    # Number of group differences (symmetric matrices)
    if(is.null(n_GD)) n_GD_theta <- ceiling(0.15*p*(p-1)/2) else n_GD_theta <- 1

    # Fixed Effects
    v_GD_theta <- sample(c(rep(1,n_GD_theta), rep(0, p*(p-1)/2-n_GD_theta)), size=p*(p-1)/2, replace=FALSE)
    v_GD_theta_signs <- sample(v_signs, size=p*(p-1)/2, replace=T)
    m_GD_theta <- matrix(0, p, p)
    m_GD_theta[upper.tri(m_GD_theta)] <- v_GD_theta*v_GD_theta_signs*Delta
    m_GD_theta <- m_GD_theta + t(m_GD_theta)
    # RE Variances
    v_GD_theta_RE <- sample(c(rep(1,n_GD_theta), rep(0, p*(p-1)/2-n_GD_theta)), size=p*(p-1)/2, replace=FALSE)
    m_GD_theta_RE <- matrix(0, p, p)
    m_GD_theta_RE[upper.tri(m_GD_theta_RE)] <- v_GD_theta_RE*Delta # here we always add to avoid negative variances
    m_GD_theta_RE <- m_GD_theta_RE + t(m_GD_theta_RE)

    ## Partial correlations across random intercepts;
    #  (translates into "mean between network" (since all predictors are centered))
    v_GD_gam <- sample(c(rep(1,n_GD_theta), rep(0, p*(p-1)/2-n_GD_theta)), size=p*(p-1)/2, replace=FALSE)
    v_GD_gam_signs <- sample(v_signs, size=p*(p-1)/2, replace=T)
    m_GD_gam <- matrix(0, p, p)
    m_GD_gam[upper.tri(m_GD_gam)] <- v_GD_gam_signs*v_GD_gam*Delta
    m_GD_gam <- m_GD_gam + t(m_GD_gam)

    # ----- Defining Group 2 -----
    # If Deltas are not specified for a given parameter, we just set differences to zero before adding them

    ## phi, fixed effects
    if(!("phi" %in% DeltaPars)) m_GD_phi <- matrix(0, p, p)
    phi2 <- phi1 + m_GD_phi
    # phi, REs
    if(!("phi_RE" %in% DeltaPars)) m_GD_phi_RE <- matrix(0, p, p)
    phi2_RE_sd <- phi1_RE_sd + m_GD_phi_RE
    ## res pcor, fixed effects
    if(!("theta" %in% DeltaPars)) m_GD_theta <- matrix(0, p, p)
    theta2 <- theta1 + m_GD_theta
    # res pcor, REs
    if(!("theta_RE" %in% DeltaPars)) m_GD_theta_RE <- matrix(0, p, p)
    theta2_RE_sd <- theta1_RE_sd + m_GD_theta_RE
    ## between pcor (pcor of random intercepts)
    if(!("gamma" %in% DeltaPars)) m_GD_gam <- matrix(0, p, p)
    gamma2 <- gamma1 + m_GD_gam

    # ----- Check whether both groups have admissible distributions -----
    # Check Group 1
    phi_eigen1 <- eigen(phi1)$values
    G1_phi1_ch <- all(abs(phi_eigen1) < phi_eig_c)
    G1_phi1_ch2 <- all(Im(phi_eigen1)==0)
    G1_theta1_ch <- all(eigen(theta1)$values > 0)
    G1_gam_ch <- all(eigen(gamma1)$values > 0)
    # Check Group 2
    phi_eigen2 <- eigen(phi2)$values
    G2_phi1_ch <- all(abs(phi_eigen2) < phi_eig_c)
    G2_phi1_ch2 <- all(Im(phi_eigen2)==0)
    G2_theta2_ch <- all(Re(eigen(theta2)$values) > 0)
    G2_gam_ch <- all(Re(eigen(gamma2)$values) > 0)

    # Combine & evaluate checks
    v_checks <- c(G1_phi1_ch, G2_phi1_ch,
                  G1_phi1_ch2, G2_phi1_ch2,
                  G1_theta1_ch, G2_theta2_ch,
                  G1_gam_ch, G2_gam_ch)

    if(all(v_checks))  {
      check_base <- TRUE
    } else {
      counter <- counter + 1
      if(counter > maxresamp) stop("Unable to generate Base model in data generation!")
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
  l_means <- list()

  # Get proportions for both groups
  v_props <- c(prop, 1-prop)

  # Progress bar
  if(pbar) pb <- txtProgressBar(0, Nsubj, style = 3)
  counter_dg <- 0

  for(j in 1:2) {  # groups

    N_j <- Nsubj*v_props[j] # number if subjects sampled in group j
    # This also means that G2 always has fewer people in the prop=2 condition

    # Storage for person-means
    m_means_j <- matrix(NA, N_j, p)

    for(i in 1:N_j) { # individuals

      counter_ij <- 0
      check_ij <- FALSE
      while(!check_ij) {

        ## Draw VAR parameters
        RE_phi_i <- matrix(NA, p, p)
        for(s1 in 1:p) for(s2 in 1:p) RE_phi_i[s1, s2] <- rnorm(1, mean=0, sd=l_phi_RE_sd[[j]][s1, s2])
        phi_i <- l_phi[[j]] + RE_phi_i # add random effect
        check_ij_1 <- all(abs(eigen(phi_i)$value) < phi_eig_c)

        ## Draw Residual pcors
        noise_theta <- matrix(0, p, p)
        for(s1 in 1:p) for(s2 in 1:p) noise_theta[s1, s2] <- rnorm(1, mean=0, sd=l_theta_RE_sd[[j]][s1, s2]) # this could be cleaner; I'm drawing >twice is much as I need
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
          if(counter_ij > maxresamp) stop("Unable to generate REs that give stable/posdef matrices")
        }

      } # end while: check stable/pos-def

      ## Draw intercepts
      gamma_j_cors <- pcor2cor(l_gamma[[j]]) #* sigma^2 # convert to cor + scale to RE variance 0.5 (new!)
      gamma_j_var <- gamma_j_cors*sigma_int # same random effects variance as for other parameters
      v_ints <- MASS::mvrnorm(n=1, mu=rep(0, p), Sigma = gamma_j_var)

      # Compute model-implied means
      means_i <- solve(diag(p) - phi_i) %*% v_ints
      means_i <- as.numeric(means_i)
      # Save means
      m_means_j[i, ] <- means_i

      # browser()

      # Sample data
      l_data[[j]][[i]] <- simulateVAR(pars = phi_i,
                                      # means = means_i, # THOSE ARE MEANS !!!!
                                      means = v_ints,
                                      residuals = theta_i_cor,
                                      Nt = Nt,
                                      burnin = var_burnin)

      # cbind(apply(l_data[[j]][[i]], 2, mean), means_i) # check!

      counter_dg <- counter_dg + 1
      if(pbar) setTxtProgressBar(pb, counter_dg)

    } # for: i

    l_means[[j]] <- m_means_j

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
                  "true" = list("basemodel" = list("l_means" = l_means,
                                                   "l_phi" = l_phi,
                                                   "l_phi_RE_sd" = l_phi_RE_sd,
                                                   "l_theta" = l_theta,
                                                   "l_theta_RE_sd" = l_theta_RE_sd,
                                                   "l_gamma" = l_gamma),
                                "GroupDiffs" = list("m_GD_phi" = m_GD_phi,
                                                    "m_GD_phi_RE" = m_GD_phi_RE,
                                                    "m_GD_theta" = m_GD_theta,
                                                    "m_GD_theta_RE" = m_GD_theta_RE,
                                                    "m_GD_gam" = m_GD_gam)),
                  "data" = list(data1, data2),
                  "dataSamp" = list("counter_base" = as.numeric(counter),
                                    "counter_datasamp" = as.numeric(counter_ij)))

  return(outlist)

} # eoF

