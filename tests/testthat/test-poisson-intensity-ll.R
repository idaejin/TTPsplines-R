# Domain integral of a binned Poisson intensity. A Berman--Turner empty
# sample is a different sum: the score that chooses lambda uses the domain.

test_that("intensity score is the domain integral, not a Berman-Turner sample", {
  # 10 x 10 lattice, Delta = 0.01. Four occupied cells and 48 empty cells at
  # eta = 0; the other 48 empty cells have exp(eta) = 4.
  eta <- c(rep(0, 4), rep(0, 48), rep(log(4), 48))
  ll <- poisson_intensity_ll(rep(1, 4), rep(0, 4), eta,
                             n_cells = 100, delta = 0.01, p_te = 1)
  expect_equal(ll, -2.44)

  # Eight empty cells drawn entirely in the high region. As a domain sample
  # their mean estimates integral 4, not the array sum 2.44, and not the
  # Berman-Turner compensator 0.04 + 8 * 0.12 * 4 = 3.88.
  bad <- poisson_intensity_ll(rep(1, 4), rep(0, 4), rep(log(4), 8),
                              n_cells = 100, delta = 0.01, p_te = 1)
  expect_equal(bad, -4)
  expect_equal(0.04 + 0.12 * 8 * 4, 3.88)

  # Unit cube: Delta = 1/M, so M * Delta * mean(exp(eta)) is the mean.
  ll1 <- poisson_intensity_ll(1, 0, c(0, log(4)), n_cells = 2, delta = 0.5,
                              p_te = 0.2)
  expect_equal(ll1, -0.2 * mean(c(1, 4)))
})
