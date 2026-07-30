# Statistical helpers shared by descriptive and prediction analyses.

cohen_kappa <- function(observed, expected, levels) {
  confusion <- table(
    factor(observed, levels = levels),
    factor(expected, levels = levels)
  )
  n <- sum(confusion)
  agreement <- sum(diag(confusion)) / n
  chance_agreement <- sum(rowSums(confusion) * colSums(confusion)) / n^2
  kappa <- if (chance_agreement < 1) {
    (agreement - chance_agreement) / (1 - chance_agreement)
  } else {
    NA_real_
  }

  list(
    kappa = kappa,
    observed = agreement,
    expected = chance_agreement,
    table = confusion
  )
}

cliffs_delta <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  ranks <- rank(c(x, y), ties.method = "average")
  u <- sum(ranks[seq_len(nx)]) - nx * (nx + 1) / 2
  2 * u / (nx * ny) - 1
}

format_p_value <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) format(p, scientific = TRUE, digits = 2) else sprintf("%.3f", p)
}

make_stratified_folds <- function(outcome, k, seed) {
  set.seed(seed)
  fold <- integer(length(outcome))
  for (level in levels(outcome)) {
    index <- sample(which(outcome == level))
    fold[index] <- rep(seq_len(k), length.out = length(index))
  }
  fold
}

operating_point <- function(outcome, probability, target_sensitivity = 0.95) {
  outcome <- as.character(outcome)
  thresholds <- sort(unique(probability), decreasing = TRUE)
  candidates <- data.table::rbindlist(lapply(thresholds, function(threshold) {
    positive <- probability >= threshold
    true_positive <- sum(positive & outcome == "TP")
    false_negative <- sum(!positive & outcome == "TP")
    true_negative <- sum(!positive & outcome == "FP")
    false_positive <- sum(positive & outcome == "FP")

    data.table::data.table(
      threshold = threshold,
      tp = true_positive,
      fn = false_negative,
      tn = true_negative,
      fp = false_positive,
      sensitivity = true_positive / (true_positive + false_negative),
      specificity = true_negative / (true_negative + false_positive)
    )
  }))

  candidates <- candidates[sensitivity >= target_sensitivity]
  if (!nrow(candidates)) stop("No empirical threshold reaches target sensitivity.")
  data.table::setorder(candidates, -specificity, -threshold)
  candidates[1L]
}

classification_metrics <- function(outcome, probability, target_sensitivity = 0.95) {
  outcome <- factor(outcome, levels = c("FP", "TP"))
  roc <- pROC::roc(
    response = outcome,
    predictor = probability,
    levels = c("FP", "TP"),
    direction = "<",
    quiet = TRUE
  )
  point <- operating_point(outcome, probability, target_sensitivity)

  list(
    auc = as.numeric(pROC::auc(roc)),
    specificity95 = point$specificity,
    brier = mean((probability - as.numeric(outcome == "TP"))^2),
    operating_point = point
  )
}

quantile_ci <- function(x) {
  unname(quantile(x, c(0.025, 0.975), na.rm = TRUE, type = 6))
}
