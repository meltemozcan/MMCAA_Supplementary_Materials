
CAI <- c("True Positive (TP)","False Positive (FP)","True Negative (TN)",
         "False Negative (FN)", "Proportion Selected (PS)", 
         "Success Ratio (SR)", "Sensitivity (SE)", "Specificity (SP)")

# Rounding that preserves trailing 0s after rounding to d decimal places
rnd <- function(vec, d) format(round(vec, d), nsmall = d)

### Model fit table functions ####
build_model_fit_table <- function(dig, path_config_fit, path_pinsearch_fit, suffix = "") {
  config_fit_stat <- readRDS(path_config_fit)
  partial_fit_stat <- readRDS(path_pinsearch_fit)
  tab_fit <- rbind("Configural invariance" = config_fit_stat, 
                   "Partial invariance" = partial_fit_stat)
  tab_fit_r <- cbind(rnd(tab_fit[ ,1], dig), tab_fit[, 2:3], rnd(tab_fit[,4:8], dig))
  tab <- cbind(c("Configural", "Partial"), tab_fit_r)
  saveRDS(tab, paste0("rds/fit_tab", suffix, ".rds"))
}
model_fit_kable <- function(
    pth, cap = "", labl = "", addft = "", cases = NULL) {
  if (length(pth) == 1) {
    tab <- readRDS(pth)
  } else {
    tabs <- lapply(pth, readRDS)
    tab <- do.call(rbind, tabs) # handle multiple paths
  }
  capt <- paste0("Model fit indices for the configural and partial invariance models", cap) 
  ft <- "CFI = Robust Comparative Fit Index. TLI = Robust Tucker-Lewis Index. RMSEA = Robust Root Mean Square Error of Approximation. SRMR = Square Root Mean Residual."
  clnames <- c("$\\chi^2$", "npar", "df", "p", "RMSEA", "CFI", "TLI", "SRMR")
  
  kbl <- tab %>% 
    kable(booktabs = TRUE, align = "c", col.names = c("Model", clnames), 
          linesep = "", escape = FALSE, row.names = FALSE, format = "latex", 
          label = labl, caption = capt) %>%  
    column_spec(2:9, "1.5cm") %>%
    kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    footnote(general = paste(ft, addft), threeparttable = TRUE, 
             footnote_as_chunk = TRUE) %>% landscape()
  # add pack_rows if multiple tables
  if (length(pth) > 1 && !is.null(cases)) {
    starts <- seq(1, by = 2, length.out = length(cases))
    ends <- starts + 1
    for (i in seq_along(cases)) {
      kbl <- pack_rows(kbl, cases[i], starts[i], ends[i])
    }
  }
  return(kbl)
}

### Model parameter estimate table functions ####
# pth is the pinsearch rds
build_parameter_table <- function(dig, pth, group_labs, suffix = "") {
  pinsearch_fit <- readRDS(pth)
  partial_fit <- pinsearch_fit[[1]]
  params <- unnest_list(lavInspect(partial_fit, "est"))
  loadings <- data.frame(lapply(params$lambda, FUN = function(x) 
    as.numeric(c(x[1:6, 1], x[7:13, 2], x[14:17, 3], x[18:19, 4]))))
  intercepts <- data.frame(params$nu)
  colnames(loadings) <- colnames(intercepts) <- group_labs
  tab <- cbind(loadings, intercepts)
  rn <- rownames(tab)
  rn[4] <- "i7r"
  rn[6] <- "i20r"
  tab <- cbind("F" = c(rep("F1", 6), rep("F2", 7), rep("F3", 4), rep("F4", 2)), 
                rn, rnd(tab, dig))
  saveRDS(tab, paste0("rds/parameter_tab", suffix, ".rds"))
}

parameter_kable <- function(dig, pth, cap = "", labl = "", addft = "",
                            groups = NULL, ng = NULL, cases = "") {
  if (length(pth) == 1) {
    tab <- readRDS(pth)
  } else {
    tabs <- lapply(pth, readRDS)
    # take the entire first df, and drop the first 2 columns of each following
    tab <- Reduce(function(x, y) cbind(x, y[, -(1:2), drop = FALSE]), tabs)
  }
  coln <- unlist(lapply(groups, function(x) rep(x, times = 2)))
  # make the second-level header (Loadings / Intercepts repeated for each case)
  header2 <- c(" " = 1,
               rep(setNames(rep(ng, 2), c("Loadings", "Intercepts")), length(cases)))
  header1 <- c(" " = 1, setNames(rep(2 * ng, length(cases)), cases))
  ft <- "" # placeholder
  capt <- paste0("Factor loading and measurement intercept estimates under partial invariance (", ng, " groups).", cap)
  kbl <- tab[-1] %>% 
    kable(booktabs = TRUE, align = "c", col.names = c("Item", coln),
          linesep = "", escape = FALSE, row.names = FALSE, label = labl,
          caption = capt) 
  
  # only add header1 if there is more than one case
  if (length(cases) == 1) {
    kbl <- kbl %>% add_header_above(header2)
  } else {
    kbl <- kbl %>% add_header_above(header2) %>% add_header_above(header1)
  }
  kbl <- kbl %>% 
    kableExtra::group_rows("1. Somatic Complaints", 1, 6, italic = T, bold = F, hline_after = T) %>% 
    kableExtra::group_rows("2. Negative Affect", 7, 13, italic = T, bold = F, hline_before = T, hline_after = T) %>%
    kableExtra::group_rows("3. (Lack of) Positive Affect", 14, 17, italic = T, bold = F, hline_before = T, hline_after = T) %>%
    kableExtra::group_rows("4. Interpersonal Relations", 18, 19, italic = T, bold = F, hline_before = T, hline_after = T) %>%
    kable_styling(latex_options = c("hold_position", "scale_down")) 
  
  if (addft != "") {
    ft <- "" # placeholder
   kbl <- kbl %>%
    footnote(general = paste(ft, addft), threeparttable = TRUE, 
             footnote_as_chunk = TRUE)
  }
  
  kbl
}
### CAI table functions ####
build_CAI_table <- function(dig, pth, ng, suffix = "") {
  cai <- readRDS(pth)
  cai <- rnd(cbind(cai$summary_mi, cai$summary), dig)
  cai <- cbind(cai[, 1:(2 * ng)],
               "ERW" = rep("-",8), cai[,(2 * ng + 1):(2 * ng + 1 + ng -2)])
  colnames(cai) <- NULL; rownames(cai) <- NULL
  reshaped_rows <- list()
  # loop through each row and interweave the columns
  for (i in 1:nrow(cai)) {
    row_combined <- rbind(cai[i, 1:ng, drop = TRUE], 
                          cai[i, (ng + 1):(2 * ng), drop = TRUE], 
                          cai[i, (2 * ng + 1):(3 * ng), drop = TRUE])
    reshaped_rows[[i]] <- row_combined
  }
  cai_tab <- do.call(rbind, reshaped_rows)
  rownames(cai_tab) <- rep(c("Strict", "Partial", "Partial (Ef)"), 8)  
  saveRDS(cai_tab, paste0("rds/CAI_tab", suffix, ".rds"))
} 

CAI_kable <- function(pth, groups, cap = "", labl = "", addft = "", cases = "", ng = NULL) {
   if (length(pth) == 1) {
    tab <- readRDS(pth)
  } else {
    tabs <- lapply(pth, readRDS)
    tab <- Reduce(function(x, y) cbind(x, y), tabs)
  }
  kbl <- tab %>% kable(
    booktabs = TRUE, align = "c", row.names = TRUE, col.names = c(" ", groups), 
    linesep = "", escape = FALSE, label = labl, # format = "latex",
    caption = paste0("Observed and expected CAI", cap, ".")) %>%
    kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    column_spec(1, "3cm") %>% column_spec(2:(ncol(tab)), "2.5cm") %>%
    footnote(
      general = paste0("Rows labeled \"Strict\" and \"Partial\" contain observed classification accuracy indices (CAI), and rows labeled \"Partial (Ef)\" contain the expected CAI for the focal group(s) if the underlying distributions matched the latent distribution of the reference group. ", addft),
      threeparttable = TRUE, footnote_as_chunk = TRUE)
  # 3 rows per CAI index
  row_ranges <- split(seq_len(nrow(tab)), rep(1:length(CAI), each = 3))
  for (i in seq_along(CAI)) {
    kbl <- kbl %>%
      pack_rows(CAI[i], row_ranges[[i]][1], row_ranges[[i]][3],
                hline_before = (i > 1), hline_after = TRUE, bold = FALSE, italic = TRUE)
  }
  # only add header1 if there is more than one case
  if (length(cases) == 1) {
    kbl <- kbl 
  } else {
    kbl <- kbl %>% 
      add_header_above(c(" " = 1, setNames(rep(ng, length(cases)), cases)))
  }
  kbl
}

###  rEf_h table functions ####
build_rEf_h_table <- function(pth, dig, ng, suffix ="") {
  cai <- readRDS(pth)
  cai <- cbind(cai$summary_mi, cai$summary)
  E_R <- cai[,(2 * ng + 1):(2 * ng + 1 + (ng - 2)), drop = FALSE]
  hs <- apply(E_R, MARGIN = 2,
              FUN = function(x) unbiasr::cohens_h(cai[, ng + 1], x))  
  rownames(hs) <- CAI
  colnames(hs) <-  paste0("$h_{\\text{", colnames(cai)[2:ng], "}}$")
  hs  
  saveRDS(round(hs, dig), paste0("rds/href_tab", suffix, ".rds"))
}


rEf_h_kable <- function(pth, groups, cap = "", labl = "", addft = "", cases = "", 
                        ng = NULL, landscape = FALSE) {
  if (length(pth) == 1) {
    tab <- readRDS(pth)
  } else {
    tabs <- lapply(pth, readRDS)
    tab <- Reduce(function(x, y) cbind(x, y), tabs)
  }
  kbl <- tab %>%
    kable(booktabs = TRUE, align = "c", row.names = TRUE, linesep = "", format = "latex",
          escape = FALSE, label =  labl,#"tbl-href", 
          caption = paste0("Cohen's $h$ effect sizes for discrepancies between observed and expected CAI", cap, ".")) %>%
     column_spec(1, "5cm") %>% column_spec(2:ng, "2cm") %>%
    kable_styling(latex_options = c("hold_position", "scale_down")) %>%
    kableExtra::footnote(
      general = paste0("Expected classification accuracy indices (CAI) are computed under partial invariance assuming that the underlying distribution for the focal group(s) match that of the reference group.", addft),
      threeparttable = TRUE, footnote_as_chunk = TRUE)
  # only add header1 if there is more than one case
  if (length(cases) > 1) {
    kbl <- kbl %>%
      add_header_above(c(" " = 1, setNames(rep((ng - 1), length(cases)), cases)))
  }
  if (landscape == TRUE) {kbl <- kbl %>% landscape() }
 
   kbl
}

### Latent variance-covariance and mean table functions ####
build_psi_alpha_table <- function(dig, pth_psi = NULL, pth_alpha = NULL,
                            est_psi = NULL, est_alpha = NULL,
                            long_rownames = TRUE, suffix = "") {
  if (!is.null(pth_psi)) psi <- readRDS(pth_psi)
  if (!is.null(pth_alpha)) alpha <- readRDS(pth_alpha)
  if (!is.null(est_psi)) psi <- est_psi
  if (!is.null(est_alpha)) alpha <- est_alpha
  
  psi_list <- lapply(psi, FUN = function(x) rnd(x, dig))
  
  
  rwnames <- c("1. Somatic Complaints", "2. Negative Affect",
               "3. (Lack of) Positive Affect", "4. Interpersonal Relations")
  
  group_tables <- lapply(seq_along(psi_list), function(i) {
    mat <- psi_list[[i]]
    mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
    mat[upper.tri(mat)] <- " "
    if (long_rownames) rownames(mat) <- rwnames
    mat
    })
  
  tab <- do.call(rbind, group_tables)
  tab <- cbind(tab, rnd(unlist(alpha), dig))
  colnames(tab) <- c("1", "2", "3", "4", "")
  tab
  saveRDS(tab, paste0("rds/psi_alpha_tab", suffix, ".rds"))
} 

psi_alpha_kable <- function(pth = NULL, groups, labl = "") {
  tab <- readRDS(pth)
  ng <- nrow(tab) / 4
  
  kbl <- tab %>% 
    kable(booktabs = TRUE, align = "c", label = labl,
          linesep = "", escape = FALSE, 
          caption = "Latent variance-covariance and mean estimates across groups.") %>%
    kable_styling(latex_options = c("hold_position", "scale_down"))
  
  for (i in 1:ng) {
    start_row <- (i - 1) * 4 + 1; end_row <- i * 4
    kbl <- kbl %>%
      pack_rows(groups[i], start_row, end_row, bold = FALSE, italic = TRUE,
                hline_before = (i > 1), hline_after  = TRUE)
  }
  kbl <- kbl %>% 
    add_header_above(c(" " = 1, "Latent variance-covariances" = 4, 
                       "Latent means" = 1))
  kbl
  }