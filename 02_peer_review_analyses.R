# ==============================================================================
# CELL 2 | ADDITIONAL ANALYSES FOR PEER REVIEW (S1–S11)
# Run after the completed Cell 1 in the same R session.
# All results are printed in the notebook.
# ==============================================================================

required_objects <- c(
  "CELL_1_COMPLETED", "show", "df", "X", "N", "SEED", "N_BOOT", "N_PERM",
  "items", "blocks", "constructs", "inhibitors", "creators", "paths", "inner",
  "table1", "table2_items", "table5", "fit_full", "scores_full", "pls_cross",
  "fit_pls", "path_values", "score_signs", "orient", "quantiles",
  "bootstrap_summary", "boot_coefs", "boot_paths", "boot_indices",
  "pairs", "htmt_values", "htmt_point", "missing_audit",
  "compare_groups", "university_mga"
)
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]
if (length(missing_objects)) {
  stop("Complete Cell 1 first. Missing objects: ",
       paste(missing_objects, collapse = ", "))
}
if (!isTRUE(CELL_1_COMPLETED)) stop("Cell 1 did not finish successfully.")

# S1. Missing-data audit.
show("S1. MISSING-DATA AUDIT", missing_audit, digits = 0)
cat("Imputed item responses: 0\n")

# S2. HTMT thresholds and percentile bootstrap intervals.
htmt_draws <- t(vapply(seq_len(N_BOOT), function(i) {
  htmt_values(X[boot_indices[, i], , drop = FALSE])
}, numeric(nrow(pairs))))
htmt_ci <- quantiles(htmt_draws, c(.025, .05, .95, .975))

htmt_inference <- data.frame(
  Construct_1 = pairs[, 1], Construct_2 = pairs[, 2], HTMT = htmt_point,
  Above_085 = htmt_point > .85, Above_090 = htmt_point > .90,
  CI90_low = htmt_ci[, 2], CI90_high = htmt_ci[, 3],
  CI95_low = htmt_ci[, 1], CI95_high = htmt_ci[, 4],
  CI90_upper_below_1 = htmt_ci[, 3] < 1
)
show("S2. HTMT THRESHOLDS AND PERCENTILE BOOTSTRAP INTERVALS", htmt_inference)

# S3. Item-level diagnostics for the two target pairs.
item_cor <- cor(X)
focus_pairs <- list(c("TO", "TI"), c("TC", "TIS"))

item_diagnostics <- do.call(rbind, lapply(focus_pairs, function(pair) {
  do.call(rbind, lapply(pair, function(own) {
    other <- setdiff(pair, own)
    own_items <- blocks[[own]]
    other_items <- blocks[[other]]
    between <- item_cor[own_items, other_items, drop = FALSE]
    strongest <- max.col(abs(between), ties.method = "first")

    data.frame(
      Pair = paste(pair, collapse = " - "),
      Item = own_items,
      Own_PLS_loading = pls_cross[own_items, own],
      Other_PLS_loading = pls_cross[own_items, other],
      Mean_abs_r_within = vapply(own_items, function(item) {
        mean(abs(item_cor[item, setdiff(own_items, item)]))
      }, numeric(1)),
      Mean_abs_r_between = rowMeans(abs(between)),
      Strongest_other_item = other_items[strongest],
      Strongest_other_item_r = between[cbind(seq_along(own_items), strongest)]
    )
  }))
}))
show("S3. ITEM DIAGNOSTICS FOR TARGET HTMT PAIRS", item_diagnostics)

# S4. Direct, specific indirect, total indirect, and total associations.
# Compute products and sums within each common bootstrap replicate.
indirect_map <- expand.grid(
  Mediator = creators, Predictor = inhibitors, stringsAsFactors = FALSE
)
association_labels <- rbind(
  data.frame(Type = "Direct", Predictor = inhibitors, Mediator = ""),
  data.frame(Type = "Specific indirect", Predictor = indirect_map$Predictor,
             Mediator = indirect_map$Mediator),
  data.frame(Type = "Total indirect", Predictor = inhibitors, Mediator = ""),
  data.frame(Type = "Total association", Predictor = inhibitors, Mediator = "")
)
association_values <- function(B) {
  a <- B[cbind(match(indirect_map$Mediator, rownames(B)),
               match(indirect_map$Predictor, colnames(B)))]
  specific <- unname(a * B["WP", indirect_map$Mediator])
  total_indirect <- vapply(inhibitors, function(k) {
    sum(specific[indirect_map$Predictor == k])
  }, numeric(1))
  direct <- unname(B["WP", inhibitors])
  c(direct, specific, total_indirect, direct + total_indirect)
}

association_draws <- t(vapply(seq_len(N_BOOT), function(i) {
  association_values(boot_coefs[, , i])
}, numeric(nrow(association_labels))))

association_results <- cbind(
  association_labels, Outcome = "WP",
  bootstrap_summary(association_values(fit_full$path_coefs), association_draws)
)
show("S4. DIRECT, INDIRECT, AND TOTAL ASSOCIATIONS WITH WORK PERFORMANCE",
     association_results)
cat("Inference uses the 95% percentile bootstrap interval.\n")
cat("The t-based columns are approximate descriptive statistics.\n")

# S5. Full-model R squared, inner VIF, and f squared.
equation_r2 <- function(scores, model_matrix) {
  outcomes <- rownames(model_matrix)[rowSums(model_matrix) > 0]
  data <- as.data.frame(scores)
  setNames(vapply(outcomes, function(outcome) {
    predictors <- colnames(model_matrix)[model_matrix[outcome, ] == 1]
    summary(lm(reformulate(predictors, outcome), data = data))$r.squared
  }, numeric(1)), outcomes)
}

r2_full <- equation_r2(scores_full, inner)
score_data <- as.data.frame(scores_full)
vif_f2 <- t(vapply(seq_len(nrow(paths)), function(i) {
  outcome <- paths$To[i]
  predictor <- paths$From[i]
  others <- setdiff(colnames(inner)[inner[outcome, ] == 1], predictor)
  auxiliary_r2 <- summary(lm(reformulate(others, predictor), data = score_data))$r.squared
  omitted_r2 <- summary(lm(reformulate(others, outcome), data = score_data))$r.squared
  full_r2 <- unname(r2_full[outcome])
  if (auxiliary_r2 >= 1 || full_r2 >= 1) stop("Perfect fit in equation for ", outcome)
  c(Inner_VIF = 1 / (1 - auxiliary_r2),
    f2 = (full_r2 - omitted_r2) / (1 - full_r2))
}, numeric(2)))

show("S5A. STRUCTURAL PATHS WITH INNER VIF AND F SQUARED", cbind(table5, vif_f2))
show("S5B. R SQUARED IN THE FULL MODEL", data.frame(
  Construct = names(r2_full),
  Predictors = unname(rowSums(inner)[names(r2_full)]),
  R2 = unname(r2_full)
))

# S6. Paired sensitivity analysis excluding techno-invasion.
keep <- setdiff(constructs, "TI")
inner_reduced <- inner[keep, keep, drop = FALSE]
blocks_reduced <- blocks[keep]
X_reduced <- X[, unlist(blocks_reduced, use.names = FALSE), drop = FALSE]

fit_reduced <- fit_pls(X_reduced, inner_reduced, blocks_reduced)
scores_reduced <- as.matrix(fit_reduced$scores[, keep, drop = FALSE])
reduced_signs <- score_signs(scores_full, scores_reduced)
aligned_scores <- sweep(scores_reduced, 2, reduced_signs, "*")
common <- paths$From != "TI" & paths$To != "TI"
common_paths <- paths[common, , drop = FALSE]

r2_reduced <- equation_r2(scores_reduced, inner_reduced)
r2_comparison <- data.frame(
  Construct = names(r2_full), R2_full = unname(r2_full),
  R2_without_TI = unname(r2_reduced[names(r2_full)])
)
r2_comparison$Delta_percentage_points <-
  100 * (r2_comparison$R2_without_TI - r2_comparison$R2_full)
show("S6A. R SQUARED IN THE FULL AND REDUCED MODELS", r2_comparison)
cat("The reduced model excludes the TI equation.\n")

reduced_point <- path_values(orient(fit_reduced$path_coefs, reduced_signs), common_paths)
reduced_draws <- matrix(NA_real_, N_BOOT, nrow(common_paths))

for (i in seq_len(N_BOOT)) {
  rows <- boot_indices[, i]
  model <- tryCatch(
    fit_pls(X_reduced[rows, , drop = FALSE], inner_reduced, blocks_reduced),
    error = function(e) {
      stop("Reduced-model bootstrap draw ", i, ": ", conditionMessage(e), call. = FALSE)
    }
  )
  signs <- score_signs(aligned_scores[rows, , drop = FALSE], as.matrix(model$scores))
  reduced_draws[i, ] <- path_values(orient(model$path_coefs, signs), common_paths)

  if (i %% 1000L == 0L) cat("Reduced-model bootstrap:", i, "/", N_BOOT, "\n")
}

reduced_ci <- quantiles(reduced_draws, c(.025, .975))
delta_ci <- quantiles(reduced_draws - boot_paths[, common, drop = FALSE], c(.025, .975))
sensitivity_paths <- data.frame(
  common_paths[, c("Hypothesis", "From", "To")],
  Beta_full = table5$Beta[common],
  Full_CI95_low = table5$CI95_low[common],
  Full_CI95_high = table5$CI95_high[common],
  Beta_without_TI = reduced_point,
  Reduced_CI95_low = reduced_ci[, 1],
  Reduced_CI95_high = reduced_ci[, 2],
  Delta_beta = reduced_point - table5$Beta[common],
  Delta_CI95_low = delta_ci[, 1],
  Delta_CI95_high = delta_ci[, 2]
)
show("S6B. PAIRED SENSITIVITY ANALYSIS FOR COMMON PATHS", sensitivity_paths)

# S7. University-type MICOM, using the same permutations as Table 6.
configural <- data.frame(
  Construct = constructs,
  Indicators_in_both_groups = vapply(blocks, paste, character(1), collapse = ", "),
  Mode = "A", Scheme = "centroid", Scaling = "Within group"
)
show("S7A. MICOM STEP 1 MODEL SPECIFICATION", configural)
cat("Verify equivalent survey administration using the study records.\n")
show("S7B. MICOM STEP 2 FOR UNIVERSITY TYPE", university_mga$step2, digits = 6)

mean_limits <- quantiles(university_mga$perm_mean, c(.025, .975))
var_limits <- quantiles(university_mga$perm_var, c(.025, .975))
observed <- university_mga$observed
pass2 <- university_mga$step2$Compositional_invariance_supported
mean_inside <- observed$Mean >= mean_limits[, 1] & observed$Mean <= mean_limits[, 2]
var_inside <- observed$LogVar >= var_limits[, 1] & observed$LogVar <= var_limits[, 2]

micom_step3 <- data.frame(
  Construct = constructs,
  Observed_mean_A_minus_B = observed$Mean,
  Mean_null_025 = mean_limits[, 1], Mean_null_975 = mean_limits[, 2],
  Observed_log_variance_ratio = observed$LogVar,
  LogVar_null_025 = var_limits[, 1], LogVar_null_975 = var_limits[, 2],
  Mean_within_limits = ifelse(pass2, mean_inside, NA),
  LogVar_within_limits = ifelse(pass2, var_inside, NA),
  Classification_given_step1 = ifelse(
    !pass2, "Step 2 not supported",
    ifelse(mean_inside & var_inside, "Full", "Partial")
  )
)
show("S7C. MICOM STEP 3 FOR UNIVERSITY TYPE", micom_step3)

# S8. First unrotated principal component diagnostic.
ev <- eigen(cor(X), symmetric = TRUE, only.values = TRUE)$values
show("S8. VARIANCE EXPLAINED BY THE FIRST UNROTATED PRINCIPAL COMPONENT",
     data.frame(First_component_variance_percent = 100 * ev[1] / sum(ev)), digits = 2)

# S9. Ordinal EFA and CFA on the same analytical sample.
extra_packages <- c("psych", "GPArotation", "lavaan")
to_install <- extra_packages[
  !vapply(extra_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(to_install)) {
  tryCatch(
    install.packages(to_install, repos = "https://cloud.r-project.org"),
    error = function(e) message("Package installation: ", conditionMessage(e))
  )
}
N_PARALLEL <- 50L

efa_status <- tryCatch({
  cat("\nRunning ordinal EFA diagnostics and parallel analysis...\n")
  r_raw <- psych::polychoric(X, smooth = FALSE, progress = FALSE)$rho
  dimnames(r_raw) <- list(items, items)
  min_eigen <- min(eigen(r_raw, symmetric = TRUE, only.values = TRUE)$values)
  smoothed <- min_eigen <= 1e-8
  r_used <- if (smoothed) psych::cor.smooth(r_raw) else r_raw
  kmo <- psych::KMO(r_used)
  bartlett <- psych::cortest.bartlett(r_used, n = N)

  set.seed(SEED)
  parallel <- psych::fa.parallel(
    X, fa = "fa", fm = "minres", cor = "poly",
    n.iter = N_PARALLEL, quant = .95, plot = FALSE
  )

  show("S9A. ORDINAL EFA FACTORABILITY AND RETENTION", data.frame(
    N = N, Overall_KMO = kmo$MSA, Bartlett_chisq = bartlett$chisq,
    Bartlett_df = bartlett$df, Bartlett_p_approx = bartlett$p.value,
    Raw_polychoric_min_eigen = min_eigen, Matrix_smoothed = smoothed,
    Parallel_factors = parallel$nfact, Parallel_repetitions = N_PARALLEL
  ))
  show("S9B. ITEM-LEVEL KMO",
       data.frame(Item = names(kmo$MSAi), MSA = unname(kmo$MSAi)))

  factor_counts <- unique(c(parallel$nfact, length(constructs)))
  factor_counts <- factor_counts[
    is.finite(factor_counts) & factor_counts >= 1 & factor_counts < ncol(X)
  ]

  for (k in factor_counts) {
    efa <- psych::fa(
      r_used, nfactors = k, n.obs = N, fm = "minres",
      rotate = if (k > 1L) "oblimin" else "none"
    )
    pattern <- as.matrix(unclass(efa$loadings))
    absolute <- abs(pattern)
    primary_index <- max.col(absolute, ties.method = "first")
    primary <- absolute[cbind(seq_len(nrow(pattern)), primary_index)]
    secondary <- if (k > 1L) {
      apply(absolute, 1, function(z) sort(z, decreasing = TRUE)[2])
    } else {
      rep(NA_real_, nrow(pattern))
    }

    diagnostics <- data.frame(
      Item = rownames(pattern),
      Intended_construct = table2_items$Construct[
        match(rownames(pattern), table2_items$Item)
      ],
      Primary_factor = colnames(pattern)[primary_index],
      Primary_abs_loading = primary, Secondary_abs_loading = secondary,
      Primary_below_040 = primary < .40,
      Cross_loading_at_least_030 = secondary >= .30,
      Loading_gap_below_020 = primary - secondary < .20
    )
    show(paste("S9C. EFA PATTERN LOADINGS | FACTORS =", k),
         cbind(diagnostics, as.data.frame(pattern)))

    if (!is.null(efa$Phi)) {
      show(paste("S9D. EFA FACTOR CORRELATIONS | FACTORS =", k),
           efa$Phi, row_names = TRUE)
    }
  }
  "Completed"
}, error = function(e) paste("Review required:", conditionMessage(e)))

cfa_status <- tryCatch({
  syntax <- paste(vapply(constructs, function(k) {
    paste(k, "=~", paste(blocks[[k]], collapse = " + "))
  }, character(1)), collapse = "\n")

  model <- lavaan::cfa(
    syntax, data = X, ordered = items,
    estimator = "WLSMV", std.lv = TRUE
  )
  if (!lavaan::lavInspect(model, "converged")) stop("Ordinal CFA did not converge.")
  admissible <- lavaan::lavInspect(model, "post.check")
  measures <- lavaan::fitMeasures(model, c(
    "chisq.scaled", "df.scaled", "pvalue.scaled", "cfi.scaled",
    "tli.scaled", "rmsea.scaled", "srmr", "cfi.robust", "tli.robust", "rmsea.robust"
  ))

  show("S9E. ORDINAL CFA FIT", data.frame(Index = names(measures), Value = unname(measures)))
  estimates <- lavaan::standardizedSolution(model)
  show("S9F. ORDINAL CFA STANDARDIZED INDICATOR LOADINGS",
       estimates[estimates$op == "=~", , drop = FALSE])

  factor_rows <- estimates$op == "~~" &
    estimates$lhs %in% constructs & estimates$rhs %in% constructs &
    estimates$lhs != estimates$rhs
  show("S9G. ORDINAL CFA STANDARDIZED FACTOR CORRELATIONS",
       estimates[factor_rows, , drop = FALSE])

  if (isTRUE(admissible)) "Completed; post-check passed" else "Inspect inadmissible solution"
}, error = function(e) paste("Review required:", conditionMessage(e)))

cat("\nEFA status:", efa_status, "\nCFA status:", cfa_status, "\n")
cat("EFA and CFA use the same analytical sample.\n")

# S10. Exploratory age-band and gender comparisons.
MIN_GROUP_N <- 20L
gender_levels <- c("Kad\u0131n", "Erkek")
exploratory_results <- list()

for (variable in c("Age", "Gender")) {
  counts <- table1[
    table1$Variable == variable, c("Category", "N", "Percent"), drop = FALSE
  ]
  counts$Eligible <- counts$N >= MIN_GROUP_N &
    !counts$Category %in% c("(Missing)", "")
  if (variable == "Gender") {
    counts$Eligible <- counts$Eligible & counts$Category %in% gender_levels
  }

  show(paste("S10A. EXPLORATORY GROUP COUNTS |", variable), counts, digits = 1)
  eligible <- counts$Category[counts$Eligible]
  excluded <- counts$Category[!counts$Eligible]

  if (length(excluded)) {
    cat("Excluded categories:", paste(excluded, collapse = ", "),
        "| Minimum group size:", MIN_GROUP_N, "\n")
  }
  if (length(eligible) < 2L) {
    cat("Skipped", variable, "because fewer than two groups qualified.\n")
    next
  }

  for (pair in combn(eligible, 2, simplify = FALSE)) {
    rows_a <- which(df[[variable]] == pair[1])
    rows_b <- which(df[[variable]] == pair[2])
    label <- paste(variable, pair[1], "vs", pair[2])
    cat("\nStarting:", label, "\n")

    result <- tryCatch(
      compare_groups(rows_a, rows_b, label),
      error = function(e) {
        cat("Comparison unavailable:", label, "|", conditionMessage(e), "\n")
        NULL
      }
    )
    if (is.null(result)) next

    show(paste("S10B. EXPLORATORY MICOM STEP 2 |", label), result$step2, digits = 6)
    exploratory_results[[label]] <- cbind(
      Variable = variable, Group_A = pair[1], Group_B = pair[2], result$paths
    )
  }
}

if (length(exploratory_results)) {
  exploratory_mga <- do.call(rbind, exploratory_results)
  exploratory_mga$Holm_p_all_age_gender_paths <-
    p.adjust(exploratory_mga$Permutation_p, method = "holm")

  show("S10C. EXPLORATORY AGE AND GENDER PATH COMPARISONS", exploratory_mga)
  cat("Holm correction covers all completed age and gender path comparisons.\n")
  cat("Interpretation is conditional on comparable administration and MICOM Step 2.\n")
} else {
  cat("No age or gender comparison completed successfully.\n")
}


# S11. Item-level sensitivity for the two high-HTMT pairs (Reviewer 3, comment 4a).
htmt_pair <- function(data, blk, a, b) {
  r <- abs(cor(data))
  aa <- r[blk[[a]], blk[[a]], drop = FALSE]
  bb <- r[blk[[b]], blk[[b]], drop = FALSE]
  mean(r[blk[[a]], blk[[b]], drop = FALSE]) /
    sqrt(mean(aa[upper.tri(aa)]) * mean(bb[upper.tri(bb)]))
}

drop_sets <- list("TO6", "TIS1", c("TO6", "TIS1"))
s11a <- do.call(rbind, lapply(drop_sets, function(drop) {
  blk <- lapply(blocks, setdiff, drop)
  Xd <- X[, unlist(blk, use.names = FALSE), drop = FALSE]
  fit_d <- fit_pls(Xd, inner, blk)
  sc <- as.matrix(fit_d$scores[, constructs, drop = FALSE])
  beta_d <- path_values(orient(fit_d$path_coefs, score_signs(scores_full, sc)))
  delta <- abs(beta_d - table5$Beta)
  data.frame(
    Dropped_items = paste(drop, collapse = " + "),
    HTMT_TO_TI = htmt_pair(Xd, blk, "TO", "TI"),
    HTMT_TC_TIS = htmt_pair(Xd, blk, "TC", "TIS"),
    Max_abs_delta_beta = max(delta),
    Path_with_max_delta = paste(paths$From, "->", paths$To)[which.max(delta)]
  )
}))
show("S11A. HTMT AND PATH CHANGES AFTER DROPPING CROSS-LOADING ITEMS", s11a)

# Full bootstrap for the model without both items, using the same resamples.
blk_both <- lapply(blocks, setdiff, c("TO6", "TIS1"))
X_both <- X[, unlist(blk_both, use.names = FALSE), drop = FALSE]
fit_both <- fit_pls(X_both, inner, blk_both)
scores_both <- as.matrix(fit_both$scores[, constructs, drop = FALSE])
signs_both <- score_signs(scores_full, scores_both)
aligned_both <- sweep(scores_both, 2, signs_both, "*")
beta_both <- path_values(orient(fit_both$path_coefs, signs_both))

draws_both <- matrix(NA_real_, N_BOOT, nrow(paths))
for (i in seq_len(N_BOOT)) {
  rows <- boot_indices[, i]
  model <- tryCatch(
    fit_pls(X_both[rows, , drop = FALSE], inner, blk_both),
    error = function(e) {
      stop("Item-sensitivity bootstrap draw ", i, ": ", conditionMessage(e), call. = FALSE)
    }
  )
  s <- score_signs(aligned_both[rows, , drop = FALSE],
                   as.matrix(model$scores[, constructs, drop = FALSE]))
  draws_both[i, ] <- path_values(orient(model$path_coefs, s))
  if (i %% 1000L == 0L) cat("Item-sensitivity bootstrap:", i, "/", N_BOOT, "\n")
}

ci_both <- quantiles(draws_both, c(.025, .975))
decision_both <- ifelse(
  !(ci_both[, 1] > 0 | ci_both[, 2] < 0), "Not supported: CI includes zero",
  ifelse(sign(beta_both) == paths$Expected_sign, "Supported",
         "Not supported: opposite direction")
)
s11b <- data.frame(
  paths[, c("Hypothesis", "From", "To")],
  Beta_full = table5$Beta, Beta_without_TO6_TIS1 = beta_both,
  CI95_low = ci_both[, 1], CI95_high = ci_both[, 2],
  Decision_full = table5$Decision, Decision_without_TO6_TIS1 = decision_both,
  Decision_changed = table5$Decision != decision_both
)
show("S11B. STRUCTURAL PATHS WITHOUT TO6 AND TIS1", s11b)
cat("Decisions changed:", sum(s11b$Decision_changed), "of", nrow(s11b), "\n")

cat("\nCell 2 completed.\n")
print(sessionInfo())