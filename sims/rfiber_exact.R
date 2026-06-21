# ===========================================================================
# Exact R-fiber solver, MODULAR. Each function is pure and unit-tested in
# sims/test_rfiber_exact.R. Used to compute, EXACTLY, the theta-arcs of the
# R-fiber x(theta)=P2X+sqrt(T)(sin(theta) U0 + cos(theta) U1) on which a given
# k-means assignment structure holds, then map to r and apply a truncated-F.
#
# Key identity (verified empirically, residual ~1e-16): along the fiber, each
# assignment margin ||x_i - m_j||^2 - ||x_i - m_k||^2 is an exact QUADRATIC FORM
#   Q(a,b) = A a^2 + B b^2 + C a b + D a + E b + Fc,   (a,b)=(sin th, cos th),
# and with t=tan(theta/2) (theta in [0,pi/2] <=> t in [0,1]),
#   Q(sin th, cos th) * (1+t^2)^2 = N(t),  a QUARTIC in t (Module 5).
# ===========================================================================

## ---- Module 1: r <-> theta map (d = m-2; R = d tan^2 theta) ----------------
r_of_theta <- function(theta, d) d * tan(theta)^2
theta_of_r <- function(r, d) atan(sqrt(r / d))

## ---- Module 2: polynomial helpers -----------------------------------------
# Horner eval; coef in INCREASING degree order: coef[1] + coef[2] t + ...
poly_eval <- function(coef, t) {
  v <- numeric(length(t))
  for (k in rev(seq_along(coef))) v <- v * t + coef[k]
  v
}
# Real roots of a polynomial (increasing-degree coef) lying in [lo, hi].
poly_real_roots_in <- function(coef, lo, hi, imtol = 1e-8) {
  while (length(coef) > 1 && abs(coef[length(coef)]) < 1e-13) coef <- coef[-length(coef)]
  if (length(coef) <= 1) return(numeric(0))
  rt <- polyroot(coef)
  re <- Re(rt)[abs(Im(rt)) <= imtol * (1 + abs(Re(rt)))]
  sort(re[re >= lo - 1e-12 & re <= hi + 1e-12])
}

## ---- Module 3: assignment-margin as a quadratic form in (alpha,beta) -------
# x_i(theta)-m_j(theta) = g_j + alpha p_j + beta s_j (affine in (alpha,beta)).
# Returns the coeffs of Q = ||.||^2(j) - ||.||^2(k):
#   A a^2 + B b^2 + C a b + D a + E b + Fc.
margin_form <- function(g_j, p_j, s_j, g_k, p_k, s_k) {
  c(A  = sum(p_j^2)   - sum(p_k^2),
    B  = sum(s_j^2)   - sum(s_k^2),
    C  = 2 * (sum(p_j * s_j) - sum(p_k * s_k)),
    D  = 2 * (sum(g_j * p_j) - sum(g_k * p_k)),
    E  = 2 * (sum(g_j * s_j) - sum(g_k * s_k)),
    Fc = sum(g_j^2)   - sum(g_k^2))
}

## ---- Module 4: evaluate the form at (alpha,beta) ---------------------------
eval_form <- function(co, alpha, beta)
  co["A"] * alpha^2 + co["B"] * beta^2 + co["C"] * alpha * beta +
  co["D"] * alpha + co["E"] * beta + co["Fc"]

## ---- Module 5: form -> quartic N(t), and arcs where Q<=0 on [0,pi/2] -------
# N(t) = Q(sin th, cos th) * (1+t^2)^2, t = tan(th/2). Increasing-degree coeffs.
form_to_Nt <- function(co) {
  A <- co[["A"]]; B <- co[["B"]]; C <- co[["C"]]; D <- co[["D"]]; E <- co[["E"]]; Fc <- co[["Fc"]]
  c(B + E + Fc, 2 * C + 2 * D, 4 * A - 2 * B + 2 * Fc, -2 * C + 2 * D, B - E + Fc)
}
# theta-intervals (rows [lo,hi]) in [0,pi/2] on which Q(sin th, cos th) <= 0.
solve_form_arcs <- function(co, tol = 1e-10) {
  Nt  <- form_to_Nt(co)
  bps <- sort(unique(c(0, poly_real_roots_in(Nt, 0, 1), 1)))   # t-breakpoints
  out <- list()
  for (i in seq_len(length(bps) - 1)) {
    tmid <- (bps[i] + bps[i + 1]) / 2; thmid <- 2 * atan(tmid)
    if (eval_form(co, sin(thmid), cos(thmid)) <= tol)
      out[[length(out) + 1]] <- c(2 * atan(bps[i]), 2 * atan(bps[i + 1]))
  }
  if (length(out) == 0) return(matrix(numeric(0), ncol = 2))
  .merge_intervals(do.call(rbind, out))
}

## ---- Module 6: interval-set algebra (intersection / union / merge) ---------
.merge_intervals <- function(m, gap = 1e-12) {
  if (is.null(m) || nrow(m) == 0) return(matrix(numeric(0), ncol = 2))
  m <- m[order(m[, 1]), , drop = FALSE]; out <- m[1, , drop = FALSE]
  for (i in seq_len(nrow(m))[-1]) {
    if (m[i, 1] <= out[nrow(out), 2] + gap) out[nrow(out), 2] <- max(out[nrow(out), 2], m[i, 2])
    else out <- rbind(out, m[i, ])
  }
  out
}
intersect_intervals <- function(a, b) {
  if (nrow(a) == 0 || nrow(b) == 0) return(matrix(numeric(0), ncol = 2))
  out <- list()
  for (i in seq_len(nrow(a))) for (j in seq_len(nrow(b))) {
    lo <- max(a[i, 1], b[j, 1]); hi <- min(a[i, 2], b[j, 2])
    if (hi > lo) out[[length(out) + 1]] <- c(lo, hi)
  }
  if (length(out) == 0) return(matrix(numeric(0), ncol = 2))
  .merge_intervals(do.call(rbind, out))
}

## ---- Module 7: truncated-F upper p-value over r-intervals ------------------
trunc_F_p <- function(R_obs, ivl_r, df1, df2) {
  if (is.null(ivl_r) || nrow(ivl_r) == 0) return(NA_real_)
  pint <- function(a, b) pf(b, df1, df2) - pf(a, df1, df2)
  den <- sum(apply(ivl_r, 1, function(z) pint(z[1], z[2]))); if (den <= 0) return(NA_real_)
  num <- sum(apply(ivl_r, 1, function(z) { a <- max(z[1], R_obs); if (z[2] <= a) 0 else pint(a, z[2]) }))
  max(0, min(1, num / den))
}

## ---- Module 8: ASSEMBLER -- exact theta-arc of one Lloyd trajectory --------
# fiber: list(cmat=P2X, amat=sqrt(T)*U0, bmat=sqrt(T)*U1)  (each n x q), so
#   x_i(theta) = cmat[i,] + sin(theta)*amat[i,] + cos(theta)*bmat[i,].
# traj: list of assignment vectors (length n): traj[[1]] = initial assignment,
#   ..., traj[[length]] = final (e.g. fit$cluster[seq_len(fit$iter)]).
# init_idx: the k initial-centroid row indices. Returns theta-interval(s) on
# which this EXACT trajectory is realized (intersection of all assignment margins).
path_arc <- function(fiber, traj, init_idx, k) {
  cmat <- fiber$cmat; amat <- fiber$amat; bmat <- fiber$bmat; n <- nrow(cmat)
  arc <- matrix(c(0, pi / 2), 1)
  restrict <- function(arc, js, cbar, abar, bbar, assign) {
    for (i in seq_len(n)) {
      a <- assign[i]
      g_j <- cmat[i, ] - cbar[a, ]; p_j <- amat[i, ] - abar[a, ]; s_j <- bmat[i, ] - bbar[a, ]
      for (jp in seq_len(k)) if (jp != a) {
        g_k <- cmat[i, ] - cbar[jp, ]; p_k <- amat[i, ] - abar[jp, ]; s_k <- bmat[i, ] - bbar[jp, ]
        co  <- margin_form(g_j, p_j, s_j, g_k, p_k, s_k)
        arc <- intersect_intervals(arc, solve_form_arcs(co))
        if (nrow(arc) == 0) return(arc)
      }
    }
    arc
  }
  # initial assignment: centroids are the rows at init_idx
  arc <- restrict(arc, NULL, cmat[init_idx, , drop = FALSE], amat[init_idx, , drop = FALSE],
                  bmat[init_idx, , drop = FALSE], traj[[1]])
  if (nrow(arc) == 0) return(arc)
  # Lloyd iterations: centroids at step l+1 are means over the step-l assignment
  if (length(traj) > 1) for (l in seq_len(length(traj) - 1)) {
    last_cl <- traj[[l]]; cur_cl <- traj[[l + 1]]
    cbar <- matrix(0, k, ncol(cmat)); abar <- cbar; bbar <- cbar
    for (j in seq_len(k)) { mem <- which(last_cl == j); if (!length(mem)) next
      cbar[j, ] <- colMeans(cmat[mem, , drop = FALSE]); abar[j, ] <- colMeans(amat[mem, , drop = FALSE])
      bbar[j, ] <- colMeans(bmat[mem, , drop = FALSE]) }
    arc <- restrict(arc, NULL, cbar, abar, bbar, cur_cl)
    if (nrow(arc) == 0) return(arc)
  }
  arc
}
