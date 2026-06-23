test_that("union inference keeps the observed path result", {
  set.seed(2022)
  n <- 30
  true_clusters <- rep(1:3, each=10)
  delta <- 4
  q <- 2
  mu <- rbind(c(delta/2,0), c(0,sqrt(3)*delta/2), c(-delta/2,0))
  sig <- 1
  X <- matrix(rnorm(n*q, sd=sig), n, q) + mu[true_clusters, ]

  old <- kmeans_inference(X, k=3, cluster_1=1, cluster_2=3,
                          sig=sig, iter.max=8, seed=2021)
  new <- kmeans_inference_union(X, k=3, cluster_1=1, cluster_2=3,
                                sig=sig, iter.max=8, seed=2021,
                                max_paths=75)

  expect_equal(new$pval_path, old$pval, tolerance=1e-10)
  expect_true(KmeansInference:::isSameIntervals(
    intervals::reduce(new$interval_path, check_valid=FALSE),
    intervals::reduce(old$final_interval, check_valid=FALSE)
  ))
  expect_true(all(new$regions$lower <= new$regions$upper))
  expect_true(any(new$regions$preserves_observed_path))
  expect_true(any(new$regions$preserves_tested_clusters))
})

test_that("union interval contains the path interval", {
  set.seed(2022)
  n <- 30
  true_clusters <- rep(1:3, each=10)
  delta <- 4
  q <- 2
  mu <- rbind(c(delta/2,0), c(0,sqrt(3)*delta/2), c(-delta/2,0))
  X <- matrix(rnorm(n*q), n, q) + mu[true_clusters, ]

  fit <- kmeans_inference_union(X, k=3, cluster_1=1, cluster_2=3,
                                sig=1, iter.max=8, seed=2021,
                                max_paths=75)
  missing <- suppressWarnings(intervals::interval_intersection(
    fit$interval_path,
    intervals::interval_complement(fit$interval_union)
  ))

  expect_equal(nrow(missing), 0)
  expect_true(fit$pval_union >= 0 && fit$pval_union <= 1)
  expect_true(fit$pval_path >= 0 && fit$pval_path <= 1)
})

test_that("region rows represent multiple components explicitly", {
  interval <- intervals::Intervals_full(matrix(c(0, 1, 2, 3),
                                               ncol=2, byrow=TRUE))
  rows <- KmeansInference:::.kmeans_interval_to_region_rows(
    interval=interval, path_id="path_test", anchor_phi=0.5,
    preserves_observed_path=TRUE, preserves_tested_clusters=FALSE
  )

  expect_equal(nrow(rows), 2)
  expect_equal(rows$set_id, c("path_test", "path_test"))
  expect_equal(rows$path_id, c("path_test", "path_test"))
  expect_equal(rows$lower, c(0, 2))
  expect_equal(rows$upper, c(1, 3))
})

test_that("general-covariance union matches the original path p-value", {
  set.seed(2022)
  n <- 30
  true_clusters <- rep(1:3, each=10)
  delta <- 4
  q <- 2
  mu <- rbind(c(delta/2,0), c(0,sqrt(3)*delta/2), c(-delta/2,0))
  X <- matrix(rnorm(n*q), n, q) + mu[true_clusters, ]
  SigInv <- solve(matrix(c(1.5, 0.4, 0.4, 0.8), 2, 2))

  # the original general-cov path emits a benign Intervals coercion warning
  old <- suppressWarnings(kmeans_inference(X, k=3, cluster_1=1, cluster_2=3,
                          iso=FALSE, SigInv=SigInv, iter.max=8, seed=2021))
  new <- kmeans_inference_union(X, k=3, cluster_1=1, cluster_2=3,
                                iso=FALSE, SigInv=SigInv, iter.max=8, seed=2021,
                                max_paths=120)

  # path p-value and statistic must reproduce the original general-cov test
  expect_equal(new$pval_path, old$pval, tolerance=1e-8)
  expect_equal(new$test_stat, old$test_stat, tolerance=1e-10)
  # union conditions on a weaker event -> never less powerful
  expect_true(new$pval_union <= new$pval_path + 1e-10)
  expect_true(new$pval_union >= 0 && new$pval_union <= 1)
  expect_true(is.null(new$sig) && !is.null(new$SigInv))
})

test_that("scan cutoff leaves the p-value unchanged but explores no more paths", {
  set.seed(2022)
  n <- 30
  true_clusters <- rep(1:3, each=10)
  delta <- 4
  q <- 2
  mu <- rbind(c(delta/2,0), c(0,sqrt(3)*delta/2), c(-delta/2,0))
  X <- matrix(rnorm(n*q), n, q) + mu[true_clusters, ]

  capped <- kmeans_inference_union(X, k=3, cluster_1=1, cluster_2=3,
                                   sig=1, iter.max=8, seed=2021, max_paths=200)
  full   <- kmeans_inference_union(X, k=3, cluster_1=1, cluster_2=3,
                                   sig=1, iter.max=8, seed=2021, max_paths=200,
                                   phi_max=Inf)

  # default cutoff must not change the result relative to scanning the whole line
  expect_equal(capped$pval_union, full$pval_union, tolerance=1e-6)
  expect_true(length(capped$paths) <= length(full$paths))
  expect_true(capped$exhaustive)
  # cutoff sits at or above the observed statistic, never below it
  expect_true(all(capped$regions$lower >= -1e-8))
})

test_that("anchor path intervals contain their anchors", {
  set.seed(2022)
  n <- 24
  true_clusters <- rep(1:3, each=8)
  delta <- 4
  q <- 2
  mu <- rbind(c(delta/2,0), c(0,sqrt(3)*delta/2), c(-delta/2,0))
  X <- matrix(rnorm(n*q), n, q) + mu[true_clusters, ]

  fit <- kmeans_inference_union(X, k=3, cluster_1=1, cluster_2=3,
                                sig=1, iter.max=8, seed=2021,
                                max_paths=60)

  for(path in fit$paths){
    expect_true(KmeansInference:::.kmeans_interval_contains(
      path$interval, path$anchor_phi, tol=1e-6
    ))
  }
})

test_that("tested cluster labels must be distinct", {
  set.seed(2022)
  true_clusters <- rep(1:3, each=8)
  mu <- rbind(c(2,0), c(0,2), c(-2,0))
  X <- matrix(rnorm(24*2), 24, 2) + mu[true_clusters, ]

  expect_error(
    kmeans_inference(X, k=3, cluster_1=1, cluster_2=1,
                     sig=1, iter.max=8, seed=2021),
    "must be different"
  )
  expect_error(
    kmeans_inference_union(X, k=3, cluster_1=1, cluster_2=1,
                           sig=1, iter.max=8, seed=2021),
    "must be different"
  )
})

test_that("general-covariance path summary reports a scalar naive p-value", {
  set.seed(2022)
  n <- 30
  true_clusters <- rep(1:3, each=10)
  delta <- 4
  q <- 2
  mu <- rbind(c(delta/2,0), c(0,sqrt(3)*delta/2), c(-delta/2,0))
  X <- matrix(rnorm(n*q), n, q) + mu[true_clusters, ]
  SigInv <- solve(matrix(c(1.5, 0.4, 0.4, 0.8), 2, 2))

  fit <- suppressWarnings(kmeans_inference(
    X, k=3, cluster_1=1, cluster_2=3,
    iso=FALSE, SigInv=SigInv, iter.max=8, seed=2021
  ))
  sm <- summary(fit)

  expect_true("p_naive" %in% names(sm))
  expect_equal(length(sm$p_naive), 1)
  expect_true(is.finite(sm$p_naive))
  expect_true(sm$p_naive >= 0 && sm$p_naive <= 1)
})
