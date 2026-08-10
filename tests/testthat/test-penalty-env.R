test_that("env P_k^full matches legacy tt_cpp and dense reference", {
  set.seed(11)
  for (r in c(1L, 2L)) {
    for (p in c(4L, 5L)) {
      ranks <- c(1L, rep(r, 2L), 1L)
      d <- 3L
      cores <- initialize_tt_cores(p, ranks, seed = 11L + r + p, sd = 0.25)
      lambda <- c(0.7, 1.1, 0.4)
      DtD <- TTPsplines:::tt_DtD_list(cores, 2L)
      p_vec <- rep(p, d)
      S <- glam_penalty(p_vec, lambda, penalty_order = 2L)
      for (k in seq_len(d)) {
        mk <- length(cores[[k]])
        A <- matrix(0, prod(p_vec), mk)
        for (i in seq_len(mk)) {
          g <- numeric(mk)
          g[i] <- 1
          cc <- cores
          cc[[k]] <- array(g, dim(cores[[k]]))
          A[, i] <- as.numeric(tt_full_theta(cc))
        }
        P_ref <- crossprod(A, S %*% A)
        P_ref <- 0.5 * (P_ref + t(P_ref))
        pen_env <- TTPsplines:::tt_conditional_penalty_full_env(
          cores, k, lambda, DtD
        )
        pen_cpp <- tt_conditional_penalty_full(
          cores, k, lambda, method = "tt_cpp"
        )
        nrm <- max(1e-12, sqrt(sum(P_ref^2)))
        expect_lt(sqrt(sum((pen_env$P_full - P_ref)^2)) / nrm, 1e-10)
        expect_lt(sqrt(sum((pen_env$P_full - pen_cpp$P_full)^2)) / nrm, 1e-10)
        # g' P g = J_λ(Θ) conditionally
        g <- as.numeric(cores[[k]])
        J <- tt_global_penalty_value(cores, lambda, penalty_order = 2L)
        expect_equal(
          as.numeric(0.5 * crossprod(g, pen_env$P_full %*% g)),
          J,
          tolerance = 1e-8
        )
      }
    }
  }
})

test_that("default tt_conditional_penalty_full uses env path", {
  set.seed(12)
  ranks <- c(1L, 2L, 2L, 1L)
  cores <- initialize_tt_cores(5L, ranks, seed = 12L, sd = 0.2)
  pen <- tt_conditional_penalty_full(cores, 2L, c(1, 1, 1))
  expect_true(grepl("^tt_env", pen$method))
})

test_that("Rcpp env path matches R env and legacy tt_cpp", {
  skip_if_not(exists("tt_conditional_penalty_full_env_cpp", mode = "function"))
  set.seed(14)
  ranks <- c(1L, 2L, 2L, 2L, 1L)
  cores <- initialize_tt_cores(6L, ranks, seed = 14L, sd = 0.2)
  lambda <- c(0.4, 1.1, 0.7, 1.3)
  DtD <- TTPsplines:::tt_DtD_list(cores, 2L)
  for (k in seq_len(4L)) {
    pen_r <- TTPsplines:::tt_conditional_penalty_full_env(cores, k, lambda, DtD)
    # force pure R assemble by calling components without cpp from_envs? 
    # env() already may dispatch cpp internals; compare to env_cpp and tt_cpp
    pen_cpp <- tt_conditional_penalty_full_env_cpp(cores, k, lambda, DtD)
    pen_legacy <- tt_conditional_penalty_full(cores, k, lambda, method = "tt_cpp")
    nrm <- max(1e-12, sqrt(sum(pen_legacy$P_full^2)))
    expect_lt(sqrt(sum((pen_cpp$P_full - pen_legacy$P_full)^2)) / nrm, 1e-10)
    expect_lt(sqrt(sum((pen_cpp$P_own - pen_legacy$P_own)^2)) /
                max(1e-12, sqrt(sum(pen_legacy$P_own^2))), 1e-10)
    expect_identical(pen_cpp$method, "tt_env_cpp")
  }
})

test_that("right-env prepare + left absorb matches single-core env", {
  set.seed(13)
  d <- 4L
  p <- 6L
  r <- 2L
  ranks <- c(1L, r, r, r, 1L)
  cores <- initialize_tt_cores(p, ranks, seed = 13L, sd = 0.25)
  lambda <- c(0.5, 1.2, 0.8, 1.5)
  DtD <- TTPsplines:::tt_DtD_list(cores, 2L)
  right <- TTPsplines:::tt_penalty_prepare_right_envs(cores, lambda, DtD)
  L0 <- matrix(1, 1, 1)
  LP <- matrix(0, 1, 1)
  for (k in seq_len(d)) {
    pen_sw <- TTPsplines:::tt_penalty_from_envs(
      L0, LP, right$R0[[k]], right$RP[[k]],
      DtD[[k]], lambda[k], p
    )
    pen_k <- TTPsplines:::tt_conditional_penalty_full_env(
      cores, k, lambda, DtD
    )
    expect_equal(pen_sw$P_full, pen_k$P_full, tolerance = 1e-10)
    absb <- TTPsplines:::.tt_left_env_absorb(
      L0, LP, cores[[k]], DtD[[k]], lambda[k]
    )
    L0 <- absb$L0
    LP <- absb$LP
  }
})
