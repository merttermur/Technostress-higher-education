# ==============================================================================
# CELL 1 | MANUSCRIPT TABLES
# Input: original Anket.csv survey export
# Output: manuscript tables printed in the notebook
# ==============================================================================

CELL_1_COMPLETED <- FALSE
DATA_FILE <- if (file.exists("/content/Anket.csv")) "/content/Anket.csv" else "Anket.csv"
SEED <- 42L
N_BOOT <- 5000L
N_PERM <- 5000L
ALPHA <- 0.05

if (!file.exists(DATA_FILE)) stop("Upload Anket.csv before running Cell 1.")

if (!requireNamespace("plspm", quietly = TRUE)) {
  tryCatch(
    install.packages("plspm", repos = "https://cloud.r-project.org"),
    error = function(e) message("CRAN installation: ", conditionMessage(e))
  )
}
if (!requireNamespace("plspm", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("gastonstat/plspm", upgrade = "never")
}
if (!requireNamespace("plspm", quietly = TRUE)) stop("The plspm package is unavailable.")

show <- function(title, x, digits = 4L, row_names = FALSE) {
  cat("\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
  out <- as.data.frame(x, check.names = FALSE)
  numeric_columns <- vapply(out, is.numeric, logical(1))
  out[numeric_columns] <- lapply(out[numeric_columns], function(z) {
    formatted <- formatC(z, format = "f", digits = digits)
    tiny <- is.finite(z) & abs(z) > 0 & abs(z) < 10^(-digits)
    formatted[tiny] <- formatC(z[tiny], format = "e", digits = 3)
    formatted[is.na(z)] <- ""
    formatted
  })
  print(
    out, row.names = row_names, right = TRUE, quote = FALSE,
    na.print = "", max = max(1000L, nrow(out) * ncol(out) + 1L)
  )
  invisible(x)
}

# Read participant records from the original export.
raw_lines <- readLines(DATA_FILE, encoding = "UTF-8", warn = FALSE)
clean_lines <- trimws(sub("^\ufeff", "", raw_lines))
records <- clean_lines[grepl('^"?2026/', clean_lines)]
if (!length(records)) stop("No participant rows matched the original export layout.")

records <- sub(";+$", "", records)
records <- sub('^"', "", records)
records <- sub('"$', "", records)
records <- gsub('""', '"', records, fixed = TRUE)

df <- read.csv(
  text = paste(records, collapse = "\n"),
  header = FALSE, colClasses = "character", check.names = FALSE,
  quote = '"', comment.char = "", fill = FALSE, encoding = "UTF-8"
)

if (!identical(dim(df), c(408L, 47L))) {
  stop("Expected 408 rows and 47 columns. Parsed: ", paste(dim(df), collapse = " x "))
}

demo_names <- c(
  "Timestamp", "Gender", "Age", "Education", "Experience",
  "UniversityType", "Faculty", "AcademicTitle"
)
item_counts <- c(
  LF = 4L, TSP = 4L, IF = 5L, TO = 7L, TI = 2L,
  TC = 5L, TIS = 5L, TU = 3L, WP = 4L
)
constructs <- names(item_counts)
blocks <- setNames(lapply(constructs, function(k) {
  paste0(k, seq_len(item_counts[[k]]))
}), constructs)
items <- unlist(blocks, use.names = FALSE)

names(df) <- c(demo_names, items)
df[] <- lapply(df, trimws)

# Normalize university labels for matching and retain the original labels for display.
university_key <- enc2utf8(df$UniversityType)
university_key <- gsub("[\u200b\u200c\u200d\ufeff]", "", university_key, perl = TRUE)
university_key <- gsub("[\u00a0\u202f]", " ", university_key, perl = TRUE)
university_key <- gsub("[[:space:]]+", " ", trimws(university_key))
university_key <- chartr("\u00dc\u00fc\u0130\u0131", "UuIi", university_key)
university_key <- tolower(university_key)
university_key <- gsub("[\u0307\u0308]", "", university_key, perl = TRUE)

university_levels <- c("devlet universitesi", "vakif universitesi")
if (anyNA(university_key) || !all(university_key %in% university_levels)) {
  print(table(university_key, useNA = "ifany"))
  stop("Unexpected normalized university labels. Inspect the categories above.")
}

state <- which(university_key == university_levels[1])
foundation <- which(university_key == university_levels[2])
cat("University groups verified: State =", length(state),
    "| Foundation =", length(foundation), "\n")

if (length(state) != 178L || length(foundation) != 230L) {
  stop("University group counts differ from the expected 178 state and 230 foundation participants.")
}

parse_likert <- function(z) {
  valid <- !is.na(z) & grepl("^[1-5]($|[[:space:][:punct:]])", z)
  out <- rep(NA_integer_, length(z))
  out[valid] <- as.integer(substr(z[valid], 1, 1))
  out
}

X <- as.data.frame(lapply(df[items], parse_likert))
N <- nrow(X)

missing_audit <- data.frame(
  Item = items,
  Missing_or_blank = vapply(items, function(k) {
    sum(is.na(df[[k]]) | df[[k]] == "")
  }, integer(1)),
  Unrecognized_nonblank = vapply(items, function(k) {
    sum(!is.na(df[[k]]) & df[[k]] != "" & is.na(X[[k]]))
  }, integer(1))
)

if (anyNA(X)) {
  show("DATA IMPORT AUDIT", missing_audit, digits = 0)
  stop("Review missing or unrecognized item responses before fitting the model.")
}
if (any(vapply(X, sd, numeric(1)) <= 0)) stop("At least one indicator is constant.")

date_text <- sub("^(\\d{4}/\\d{1,2}/\\d{1,2}).*$", "\\1", df$Timestamp)
response_dates <- as.Date(date_text, format = "%Y/%m/%d")
cat("Participants:", N, "| Indicators:", ncol(X), "| Imputed item responses: 0\n")
cat("Dates parsed:", sum(!is.na(response_dates)), "/", N, "\n")
if (any(!is.na(response_dates))) {
  cat("Observed response dates:",
      paste(format(range(response_dates, na.rm = TRUE), "%Y-%m-%d"), collapse = " to "),
      "\n")
}

# Table 1. Participant characteristics.
demo_levels <- list(
  Age = c("27 ve alt\u0131", "28-35", "36-45", "46 ve \u00fcst\u00fc"),
  Education = c("Lisans", "Y\u00fcksek Lisans", "Doktora"),
  Experience = c(
    "1-5 Y\u0131l", "6-10 Y\u0131l", "11-15 Y\u0131l",
    "16-20 Y\u0131l", "21 ve \u00dczeri"
  ),
  UniversityType = c("Devlet \u00dcniversitesi", "Vak\u0131f \u00dcniversitesi")
)
table1 <- do.call(rbind, lapply(setdiff(demo_names, "Timestamp"), function(k) {
  z <- if (k == "UniversityType") {
    demo_levels$UniversityType[match(university_key, university_levels)]
  } else {
    df[[k]]
  }
  z[is.na(z) | z == ""] <- "(Missing)"
  levels <- unique(c(demo_levels[[k]], sort(unique(z))))
  counts <- table(factor(z, levels = levels))
  counts <- counts[counts > 0]
  data.frame(
    Variable = k, Category = names(counts), N = as.integer(counts),
    Percent = 100 * as.integer(counts) / N
  )
}))
show("TABLE 1. PARTICIPANT CHARACTERISTICS", table1, digits = 1)

# Define the full PLS path model.
inhibitors <- c("LF", "TSP", "IF")
creators <- c("TO", "TI", "TC", "TIS", "TU")

inner <- matrix(0, length(constructs), length(constructs),
                dimnames = list(constructs, constructs))
inner[creators, inhibitors] <- 1
inner["WP", c(inhibitors, creators)] <- 1

paths <- data.frame(
  Hypothesis = c(
    paste0("H1", letters[1:5]), paste0("H2", letters[1:5]),
    paste0("H3", letters[1:5]), paste0("H4", letters[1:5]),
    paste0("H5", letters[1:3])
  ),
  From = c(rep(inhibitors, each = 5), creators, inhibitors),
  To = c(rep(creators, 3), rep("WP", 8)),
  Expected_sign = c(rep(-1L, 20), rep(1L, 3))
)
stopifnot(nrow(paths) == sum(inner))

fit_pls <- function(data, model_matrix = inner, model_blocks = blocks) {
  fit <- plspm::plspm(
    Data = data, path_matrix = model_matrix,
    blocks = lapply(model_blocks, match, table = names(data)),
    modes = rep("A", length(model_blocks)),
    scheme = "centroid", scaled = TRUE, boot.val = FALSE, dataset = FALSE
  )
  if (!all(is.finite(fit$path_coefs)) || !all(is.finite(as.matrix(fit$scores)))) {
    stop("The fitted PLS model returned non-finite estimates.")
  }
  fit
}

path_values <- function(coefs, mapping = paths) {
  unname(coefs[cbind(
    match(mapping$To, rownames(coefs)), match(mapping$From, colnames(coefs))
  )])
}

score_signs <- function(reference, candidate) {
  keys <- colnames(candidate)
  r <- diag(cor(reference[, keys, drop = FALSE], candidate))
  if (any(!is.finite(r) | abs(r) < 1e-10)) stop("Could not orient construct scores.")
  setNames(sign(r), keys)
}

orient <- function(B, signs) {
  B * outer(signs[rownames(B)], signs[colnames(B)])
}

quantiles <- function(draws, probs) {
  t(apply(draws, 2, quantile, probs = probs, names = FALSE))
}

bootstrap_summary <- function(point, draws) {
  interval <- quantiles(draws, c(.025, .975))
  se <- apply(draws, 2, sd)
  t_approx <- rep(NA_real_, length(point))
  valid <- is.finite(se) & se > 0
  t_approx[valid] <- unname(point[valid]) / se[valid]
  data.frame(
    Beta = unname(point), Bootstrap_SE = se, t_approx = t_approx,
    p_t_approx = 2 * pt(-abs(t_approx), df = N - 1L),
    CI95_low = interval[, 1], CI95_high = interval[, 2],
    CI95_excludes_zero = interval[, 1] > 0 | interval[, 2] < 0
  )
}

fit_full <- fit_pls(X)
scores_full <- as.matrix(fit_full$scores[, constructs, drop = FALSE])
loadings_full <- setNames(fit_full$outer_model$loading,
                          as.character(fit_full$outer_model$name))
stopifnot(all(items %in% names(loadings_full)))

# Table 2. Construct and indicator statistics.
table2 <- do.call(rbind, lapply(constructs, function(k) {
  l <- loadings_full[blocks[[k]]]
  composite <- rowMeans(X[, blocks[[k]], drop = FALSE])
  data.frame(
    Construct = k, Items = length(l), Alpha = fit_full$unidim[k, "C.alpha"],
    CR = sum(l)^2 / (sum(l)^2 + sum(1 - l^2)), AVE = mean(l^2),
    Mean = mean(composite), SD = sd(composite)
  )
}))
table2_items <- data.frame(
  Construct = rep(constructs, item_counts), Item = items,
  Loading = unname(loadings_full[items]), Mean = colMeans(X),
  SD = vapply(X, sd, numeric(1))
)

show("TABLE 2A. RELIABILITY, CONVERGENT VALIDITY, AND CONSTRUCT DESCRIPTIVES", table2)
show("TABLE 2B. INDICATOR LOADINGS AND DESCRIPTIVES", table2_items)

# Table 3. Fornell-Larcker criterion.
fl <- cor(scores_full)
diag(fl) <- sqrt(table2$AVE[match(rownames(fl), table2$Construct)])
fl_lower <- fl
fl_lower[upper.tri(fl_lower)] <- NA_real_
show("TABLE 3. FORNELL-LARCKER MATRIX", fl_lower, row_names = TRUE)

max_absolute_r <- vapply(constructs, function(k) {
  max(abs(fl[k, setdiff(constructs, k)]))
}, numeric(1))
show("TABLE 3 CHECK. SQUARE ROOT OF AVE VS ABSOLUTE CORRELATIONS", data.frame(
  Construct = constructs, Sqrt_AVE = diag(fl), Max_absolute_r = max_absolute_r,
  Criterion_met = diag(fl) > max_absolute_r
))

# Table 4. HTMT based on absolute item correlations.
pairs <- t(combn(constructs, 2))
htmt_values <- function(data) {
  r <- abs(cor(data))
  vapply(seq_len(nrow(pairs)), function(i) {
    a <- blocks[[pairs[i, 1]]]
    b <- blocks[[pairs[i, 2]]]
    aa <- r[a, a, drop = FALSE]
    bb <- r[b, b, drop = FALSE]
    denom <- sqrt(mean(aa[upper.tri(aa)]) * mean(bb[upper.tri(bb)]))
    if (!is.finite(denom) || denom <= 0) stop("Invalid HTMT denominator.")
    mean(r[a, b, drop = FALSE]) / denom
  }, numeric(1))
}

htmt_point <- htmt_values(X)
htmt_matrix <- matrix(NA_real_, length(constructs), length(constructs),
                      dimnames = list(constructs, constructs))
htmt_matrix[pairs[, 2:1, drop = FALSE]] <- htmt_point
show("TABLE 4. HTMT USING ABSOLUTE ITEM CORRELATIONS", htmt_matrix, row_names = TRUE)

# Table 5. Bootstrap path estimates.
boot_indices <- matrix(NA_integer_, N, N_BOOT)
boot_coefs <- array(
  NA_real_, c(length(constructs), length(constructs), N_BOOT),
  dimnames = list(constructs, constructs, NULL)
)

set.seed(SEED)
for (b in seq_len(N_BOOT)) {
  rows <- sample.int(N, N, replace = TRUE)
  model <- tryCatch(
    fit_pls(X[rows, , drop = FALSE]),
    error = function(e) stop("Full-model bootstrap draw ", b, ": ", conditionMessage(e))
  )
  signs <- score_signs(scores_full[rows, , drop = FALSE], as.matrix(model$scores))
  boot_coefs[, , b] <- orient(model$path_coefs, signs)
  boot_indices[, b] <- rows
  if (b %% 1000L == 0L) cat("Full-model bootstrap:", b, "/", N_BOOT, "\n")
}

boot_paths <- t(vapply(seq_len(N_BOOT), function(b) {
  path_values(boot_coefs[, , b])
}, numeric(nrow(paths))))

table5 <- cbind(paths, bootstrap_summary(path_values(fit_full$path_coefs), boot_paths))
table5$Decision <- ifelse(
  !table5$CI95_excludes_zero, "Not supported: CI includes zero",
  ifelse(sign(table5$Beta) == table5$Expected_sign, "Supported",
         "Not supported: opposite direction")
)
show("TABLE 5. STRUCTURAL PATHS AND HYPOTHESIS DECISIONS", table5)
cat("Approximate p values use a two-sided t reference with N - 1 df.\n")
cat("Decisions use the 95% percentile bootstrap interval and expected direction.\n")
print(table(table5$Decision))

# Table 6. Two-group permutation function reused in Cell 2.
compare_groups <- function(rows_a, rows_b, label, progress_every = 1000L) {
  pool <- c(rows_a, rows_b)
  if (anyDuplicated(pool)) stop("Groups must contain disjoint observations.")

  pool_data <- X[pool, , drop = FALSE]
  pooled_x <- as.matrix(pool_data)
  anchors <- vapply(blocks, function(z) {
    rowMeans(scale(pooled_x[, z, drop = FALSE]))
  }, numeric(length(pool)))

  pooled_fit <- fit_pls(pool_data)
  pooled_scores <- as.matrix(pooled_fit$scores[, constructs, drop = FALSE])
  pooled_scores <- sweep(pooled_scores, 2, score_signs(anchors, pooled_scores), "*")

  group_fit <- function(group_rows) {
    data <- pool_data[group_rows, , drop = FALSE]
    item_sd <- vapply(data, sd, numeric(1))
    if (any(!is.finite(item_sd) | item_sd <= 0)) stop("Constant indicator in a group.")

    model <- fit_pls(data)
    weights <- setNames(model$outer_model$weight, as.character(model$outer_model$name))
    projected <- vapply(blocks, function(z) {
      drop(pooled_x[, z, drop = FALSE] %*% (weights[z] / item_sd[z]))
    }, numeric(length(pool)))
    signs <- score_signs(anchors, projected)

    list(
      Scores = sweep(projected, 2, signs, "*"),
      Beta = path_values(orient(model$path_coefs, signs))
    )
  }

  group_stats <- function(a, b) {
    fit_a <- group_fit(a)
    fit_b <- group_fit(b)
    s_a <- pooled_scores[a, , drop = FALSE]
    s_b <- pooled_scores[b, , drop = FALSE]
    list(
      Beta_A = fit_a$Beta, Beta_B = fit_b$Beta, Delta = fit_a$Beta - fit_b$Beta,
      Correlation = diag(cor(fit_a$Scores, fit_b$Scores)),
      Mean = colMeans(s_a) - colMeans(s_b),
      LogVar = log(apply(s_a, 2, var)) - log(apply(s_b, 2, var))
    )
  }

  ia <- seq_along(rows_a)
  ib <- seq.int(length(rows_a) + 1L, length(pool))
  observed <- group_stats(ia, ib)

  perm_c <- perm_mean <- perm_var <- matrix(NA_real_, N_PERM, length(constructs))
  perm_delta <- matrix(NA_real_, N_PERM, nrow(paths))

  set.seed(SEED)
  for (i in seq_len(N_PERM)) {
    shuffled <- sample.int(length(pool))
    draw <- tryCatch(
      group_stats(shuffled[ia], shuffled[ib]),
      error = function(e) stop(label, ", permutation ", i, ": ", conditionMessage(e))
    )
    perm_c[i, ] <- draw$Correlation
    perm_mean[i, ] <- draw$Mean
    perm_var[i, ] <- draw$LogVar
    perm_delta[i, ] <- draw$Delta
    if (i %% progress_every == 0L) cat(label, "permutations:", i, "/", N_PERM, "\n")
  }

  stopifnot(all(is.finite(c(perm_c, perm_mean, perm_var, perm_delta))))

  c05 <- apply(perm_c, 2, quantile, probs = ALPHA, names = FALSE)
  pass <- setNames(observed$Correlation >= c05, constructs)
  equation_pass <- vapply(paths$To, function(outcome) {
    predictors <- colnames(inner)[inner[outcome, ] == 1]
    all(pass[c(outcome, predictors)])
  }, logical(1))

  p_comp <- (1 + colSums(sweep(perm_c, 2, observed$Correlation, "<="))) / (N_PERM + 1)
  p_mga <- (1 + colSums(sweep(abs(perm_delta), 2, abs(observed$Delta), ">="))) /
    (N_PERM + 1)

  list(
    paths = data.frame(
      paths[, c("Hypothesis", "From", "To")],
      N_A = length(rows_a), N_B = length(rows_b),
      Beta_A = observed$Beta_A, Beta_B = observed$Beta_B,
      Delta_A_minus_B = observed$Delta,
      Permutation_p = p_mga, Equation_MICOM_step2 = equation_pass
    ),
    step2 = data.frame(
      Construct = constructs, Observed_c = observed$Correlation,
      Permutation_c_05 = c05, Lower_tail_p = p_comp,
      Compositional_invariance_supported = unname(pass)
    ),
    observed = observed, perm_mean = perm_mean, perm_var = perm_var
  )
}

university_mga <- compare_groups(state, foundation, "University type")
table6 <- university_mga$paths
names(table6)[match(
  c("N_A", "N_B", "Beta_A", "Beta_B", "Delta_A_minus_B"), names(table6)
)] <- c("N_State", "N_Foundation", "Beta_State", "Beta_Foundation",
        "Delta_State_minus_Foundation")
table6$Holm_p <- p.adjust(table6$Permutation_p, method = "holm")
show("TABLE 6. STATE VS FOUNDATION UNIVERSITY PATH COMPARISONS", table6)
cat("Path comparisons are conditional on configural and compositional invariance.\n")

# Table 7. PLS indicator cross-loadings.
pls_cross <- cor(X, scores_full)
table7 <- cbind(table2_items[, c("Construct", "Item")], as.data.frame(pls_cross))
show("TABLE 7. PLS INDICATOR CROSS-LOADINGS", table7)

CELL_1_COMPLETED <- TRUE
cat("\nCell 1 completed. Run Cell 2 in the same R session.\n")