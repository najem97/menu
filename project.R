# Sex‑Specific Differences in Cardiovascular Adverse Events for Common Anti‑arrhythmics (FAERS)
# -------------------------------------------------------------------------
# Re‑ordered & streamlined script based on recommended analytical flow
# Author: Mohammad Najm Dadam · Last update: 21.08.2025
# -------------------------------------------------------------------------

# 0.  House‑keeping ---------------------------------------------------------

required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "stringr", "janitor", "purrr", "ggplot2",
  "scales", "forcats", "ggrepel", "viridis", "epiR", "broom", "writexl", "exact2x2")

if (interactive()) {
  new_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  if (length(new_pkgs)) install.packages(new_pkgs)
}

purrr::walk(required_pkgs, library, character.only = TRUE)

# ---- I/O ----
file_path <- "D:/ORC/FDA - Events/Common drugs - AEs/8 Drugs.xlsx"   # <- edit

# 1 ── Helper: canonical drug map -------------------------------------------

canonical_drugs <- tibble::tribble(
  ~drug_simple,   ~pattern,
  "Amiodarone",  "amiodarone",
  "Sotalol",     "sotalol",
  "Dofetilide",  "dofetilide",
  "Flecainide",  "flecainide",
  "Adenosine",   "adenosine",
  "Verapamil",   "verapamil",
  "Propafenone", "propafenone",
  "Dronedarone", "dronedarone"
)

# 2 ── Helper: canonical CV term map ----------------------------------------
# "Atrial arrhythmia",             "atrial (fibrillation|flutter|tachycardia)",

cv_map <- tibble::tribble(
  ~canon,                           ~pattern,
  "Torsade de pointes",            "torsad",
  "Supraventricular arrhythmia",   "supraventricular( extrasystoles)?",
  "Ventricular arrhythmia",        "ventricular (tachyarrhythmia|arrhythmia|extrasystoles|tachycardia|fibrillation|flutter)|idioventricular",
  "Tachycardia",                   "(?<!induced\\s|foetal\\s)\\btachycardia\\b(?!\\s(induced|foetal|fetal))|tachyarrhythmia",
  "Bradycardia",                   "bradycardia|bradyarrhythmia",
  "Conduction block",              "block|atrioventricular dissociation",
  "Palpitations/Syncope",          "palpitations|presyncope|syncope",
  "Cardiac arrest/Sudden death",   "asystole|arrest|resuscitation|sudden cardiac death",
  "Cardiogenic shock",             "cardiogenic shock",
  "Myocardial infarction/ACS",     "myocardial infarction|acute coronary|myocardial ischaemia",
  "Angina",                        "angina",
  "Hypertension",                  "(?<!pulmonary\\s|pulmonary arterial\\s)(?<!essential\\s|portal\\s|labile\\s|white coat\\s|intracranial\\s)\\bhypertension\\b|blood pressure .*increased|hypertensive crisis",
  "Hypotension",                   "hypotension|blood pressure .*decreased|orthostatic",
  "QT prolongation",               "qt .*prolong",
  "QRS complex prolongation",      "qrs complex prolonged",
  "Oedema",                        "oedema peripheral|pulmonary oedema",
  "Pericarditis/Pericardial effusion", "pericarditis|pericardial effusion",
  "Cardiotoxicity",                "cardiotoxicity",
  "Potassium disturbance",         "hyperkalaemia|hypokalaemia"
)

clean_faers_with_audit <- function(path) {
  library(readxl)
  library(janitor)
  library(dplyr)
  library(stringr)
  library(tidyr)
  
  # (0) Read and clean raw data
  raw0 <- read_excel(path) |>
    clean_names() |>
    mutate(.orig_row = row_number(), drop_reason = "kept")
  message("Rows READ: ", nrow(raw0))
  
  # (a) Remove duplicates (by case_id, drug_name, event trio as in your code)
  raw0 <- raw0 |>
    mutate(.dup_flag = duplicated(paste(case_id, drug_name, event)))
  duplicated_reports <- raw0 |>
    filter(.dup_flag) |>
    mutate(drop_reason = "duplicate")
  raw1 <- raw0 |> filter(!.dup_flag)
  message("Duplicates removed: ", nrow(duplicated_reports))
  
  # Save raw after deduplication
  raw_dedup <- raw1
  
  # (b) Remove unknown/missing sex reports (but save them)
  raw1 <- raw1 |>
    mutate(gender_cln = toupper(str_trim(gender)))
  unknown_sex <- raw1 |>
    filter(!gender_cln %in% c("F", "M")) |>
    mutate(drop_reason = "unknown_sex")
  raw1 <- raw1 |>
    filter(gender_cln %in% c("F", "M")) |>
    mutate(sex = factor(gender_cln, levels = c("F", "M")))
  message("Unknown-sex rows DROPPED: ", nrow(unknown_sex))
  
  # (b2) Remove entire cases that include "drug ineffective" or "no adverse event"
  ineff_rows <- raw1 |>
    filter(str_detect(str_to_lower(event), "drug ineffective|no adverse event"))
  ineff_case_ids <- unique(ineff_rows$case_id)
  raw1 <- raw1 |>
    filter(!case_id %in% ineff_case_ids)
  message("Rows removed due to 'drug ineffective'/'no adverse event' (by case): ", nrow(ineff_rows))
  
  raw3 <- raw1 %>%
    mutate(drug_line = str_to_lower(coalesce(drug_name, ""))) %>%
    rowwise() %>%
    mutate(
      hit = { w <- which(str_detect(drug_line, canonical_drugs$pattern)); if (length(w)) w[1] else NA_integer_ },
      drug_simple = if (!is.na(hit)) canonical_drugs$drug_simple[hit] else NA_character_
    ) %>%
    ungroup() %>%
    mutate(drop_reason = if_else(is.na(drug_simple) & drop_reason == "kept",
                                 "non_target_drug", drop_reason)) %>%
    select(-drug_line, -hit)
  
  # ── NEW: Build df_desc as UNIQUE case×drug, merging all AE text per pair
  # We split any multi-AE cells, clean, and aggregate unique events back.
  df_desc <- raw3 %>%
    filter(drop_reason == "kept") %>%
    mutate(event = coalesce(event, "")) %>%
    separate_rows(event, sep = "[,;/]") %>%
    mutate(event = str_squish(event)) %>%
    group_by(case_id, drug_simple) %>%
    summarise(
      sex            = first(sex),
      age            = first(age_in_report),            # or dplyr::first(age[!is.na(age)], default = NA_real_)
      isr            = first(isr),            # adjust default type if needed
      date_received  = first(date_received),  # or: min(date_received, na.rm = TRUE) for earliest date
      events_all     = paste(sort(unique(event[event != ""])), collapse = "; "),
      n_unique_AEs   = n_distinct(event[event != ""]),
      n_rows_merged  = n(),
      .groups        = "drop"
    )
  message("df_desc rows (UNIQUE case×drug, merged AEs): ", nrow(df_desc))
  
  # (e) Explode multi-event strings for event-level analytics
  raw4 <- raw3 |>
    separate_rows(event, sep = "[,;/]") |>
    mutate(event = str_squish(event))
  message("Rows AFTER event explode: ", nrow(raw4))
  
  # (f) CV event labeling (expects global 'cv_map' with columns: pattern, canon)
  # IMPORTANT: Do NOT drop non-CV events; label them for denominator use in ROR/PRR.
  raw5 <- raw4 |>
    mutate(event_clean = str_to_lower(event) |>
             str_replace_all("[^[:alnum:] ]", " ") |>
             str_squish()) |>
    rowwise() |>
    mutate(canon_term = {
      idx <- which(str_detect(event_clean, cv_map$pattern))
      if (length(idx)) cv_map$canon[idx[1]] else NA_character_
    }) |>
    ungroup() |>
    mutate(
      event_group = if_else(!is.na(canon_term), "cv", "non_cv"),
      # preserve only target-drug rows for analytics; others go to 'dropped'
      drop_reason = if_else(is.na(drug_simple) & drop_reason == "kept", "non_target_drug", drop_reason)
    )
  
  # Final analytic tables
  df_events_all <- raw5 |> filter(drop_reason == "kept")  # case × drug × event, with event_group ∈ {cv, non_cv}
  df_cv_only    <- df_events_all |> filter(event_group == "cv")  # convenience view; do NOT use alone for ROR/PRR
  dropped       <- raw5 |> filter(drop_reason != "kept")
  
  message("FINAL df_events_all rows (case × drug × event): ", nrow(df_events_all))
  message("...of which CV events: ", nrow(df_cv_only))
  
  return(list(
    raw_data                 = raw0,
    raw_dedup                = raw_dedup,
    df_desc                  = df_desc,
    df_events_all            = df_events_all,
    df_cv_only               = df_cv_only,
    dropped                  = dropped,
    duplicated_reports       = duplicated_reports,
    unknown_sex              = unknown_sex,
    drug_ineffective_reports = ineff_rows
  ))
}

# ── Run cleaning -----------------------------------------------------------
cleaned      <- clean_faers_with_audit(file_path)
raw_data_tbl <- cleaned$raw_data
duplicated_tbl <- cleaned$duplicated_reports
df_desc      <- cleaned$df_desc   # descriptive level
df_cv        <- cleaned$df_cv_only     # event level
df_events_all = cleaned$df_events_all
dropped_tbl  <- cleaned$dropped
unknown_sex <- cleaned$unknown_sex
drug_ineffective_tbl <- cleaned$drug_ineffective_reports

# install.packages("writexl")  # run once
library(writexl)
dir.create("exports", showWarnings = FALSE)

# (optional) ensure proper date class for nice Excel dates
df_desc_export <- df_desc
if ("date_received" %in% names(df_desc_export)) {
  df_desc_export$date_received <- suppressWarnings(as.Date(df_desc_export$date_received))
}

write_xlsx(df_desc_export, "exports/df_desc.xlsx")

# 3 ── Descriptive epidemiology (df_desc) -----------------------------------
df_plot <- bind_rows(
  df_desc |> mutate(sex = as.character(sex)),  # known sex
  unknown_sex |> 
    mutate(sex = "Unknown", drug_simple = NA_character_) |> 
    separate_rows(drug_name, sep = "[,;/]") |> 
    mutate(drug_name = str_squish(drug_name)) |> 
    rowwise() |> 
    mutate(drug_simple = {
      idx <- which(str_detect(tolower(drug_name), canonical_drugs$pattern))
      if (length(idx)) canonical_drugs$drug_simple[idx[1]] else NA_character_
    }) |> 
    ungroup()
) |> 
  filter(!is.na(drug_simple))  # remove non-mappable drugs

# Count number of reports by drug and sex
plot_data <- df_plot |> 
  count(drug_simple, sex)

ggplot(plot_data, aes(x = drug_simple, y = n, fill = sex)) +
  geom_col(position = position_dodge(width = 0.7)) +
  scale_fill_manual(
    values = c(F = "#C51B7D", M = "#0087BD", Unknown = "gray70"),
    labels = c(F = "Female", M = "Male", Unknown = "Unknown")
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Distribution of FAERS Reports for Antiarrhythmic Drugs by Drug Type and Patient Sex",
    subtitle = "Includes all mapped drug reports with known sex, post-deduplication",
    x = NULL,
    y = "Number of Reports",
    fill = "Sex"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

age_tbl <- df_desc |> summarise(
  n_total = n(),
  n_age   = sum(!is.na(age)),
  pct_age = round(n_age / n_total * 100, 1),
  median  = median(age, na.rm = TRUE),
  q1      = quantile(age, .25, na.rm = TRUE),
  q3      = quantile(age, .75, na.rm = TRUE)
)

# 4 ── Sex‑specific AE profiles (df_cv) -------------------------------------
profile_df <- df_cv |> 
  count(drug_simple, sex, canon_term) |> 
  group_by(drug_simple, sex) |> 
  mutate(pct = n / sum(n) * 100) |> 
  ungroup()

# Fisher exact test Female vs Male within each drug for every AE ------------
fisher_wide <- df_cv |> 
  count(drug_simple, sex, canon_term) |> 
  pivot_wider(names_from = sex, values_from = n, values_fill = 0) |> 
  group_by(drug_simple) |> 
  mutate(tot_F = sum(F),
         tot_M = sum(M)) |> 
  ungroup() |> 
  rowwise() |> 
  mutate(
    ft = list(fisher.test(matrix(c(F, tot_F - F, M, tot_M - M), nrow = 2))),
    p = ft$p.value,
    ci = list(ft$conf.int)
  ) |> 
  unnest_wider(ci, names_sep = "_") |> 
  rename(l95 = ci_1, u95 = ci_2) |> 
  ungroup() |> 
  group_by(drug_simple) |> 
  mutate(p_adj = p.adjust(p, method = "BH")) |> 
  ungroup() |> 
  mutate(sig = case_when(
    p_adj < 0.001 ~ "***",
    p_adj < 0.01  ~ "**",
    p_adj < 0.05  ~ "*",
    TRUE          ~ ""
  )) |> 
  select(drug_simple, canon_term, p, p_adj, l95, u95, sig)

# replicate significance flag for both sexes so join works ------------------
fisher_long <- tidyr::crossing(fisher_wide, sex = c("F", "M"))

heat_df <- profile_df |> 
  left_join(
    fisher_long |> 
      mutate(
        `Fisher p` = p,
        `BH‑adjusted p` = p_adj,
        `Fisher CI (95%)` = sprintf("%.2f–%.2f", l95, u95),
        `Significance` = sig
      ) |> 
      select(drug_simple, canon_term, sex, `Fisher p`, `BH‑adjusted p`, `Fisher CI (95%)`, `Significance`),
    by = c("drug_simple", "canon_term", "sex")
  )

# 1) Build the y-axis order: F then M within each drug
y_order <- heat_df %>%
  distinct(drug_simple, sex) %>%
  mutate(sex = factor(sex, levels = c("F","M"))) %>%
  arrange(drug_simple, sex) %>%
  mutate(y_lab = paste(drug_simple, sex, sep = " – "),
         idx   = row_number())

y_levels <- y_order$y_lab

# 2) Positions for separators: after the LAST row of each drug (except the last drug)
boundary_df <- y_order %>%
  group_by(drug_simple) %>%
  summarise(last_idx = max(idx), .groups = "drop") %>%
  arrange(last_idx) %>%
  mutate(yintercept = last_idx + 0.5) %>%
  slice(1:(n()-1))     # drop the final boundary

# 3) Plot with separators between drugs
plot_df <- heat_df %>%
  mutate(y = factor(paste(drug_simple, sex, sep = " – "), levels = y_levels))

ggplot(plot_df, aes(canon_term, y, fill = pct)) +
  geom_tile(color = "grey92") +
  # draw lines BETWEEN drugs
  geom_hline(data = boundary_df, aes(yintercept = yintercept),
             inherit.aes = FALSE, colour = "black", size = 0.4) +
  geom_text(aes(label = Significance), colour = "white", size = 3) +
  scale_fill_viridis(name = "% within drug & sex") +
  labs(x = "Cardiovascular adverse event category", y = NULL,
       title = "Sex-specific distribution of cardiovascular AEs by drug and event type",
       subtitle = "Asterisks mark Female vs Male differences (BH–adjusted p < 0.05)") +
  theme_minimal(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(heat_df, aes(canon_term, paste(drug_simple, sex, sep = " – "), fill = pct)) +
  geom_tile(color = "grey92") +
  geom_text(aes(label = Significance), colour = "white", size = 3) +
  scale_fill_viridis(name = "% within drug & sex") +
  labs(x = "Cardiovascular adverse event category", y = NULL,
       title = "Sex-specific distribution of cardiovascular AEs by drug and event type",
       subtitle = "Asterisks mark Female vs Male differences (BH–adjusted p < 0.05)",
  ) +
  theme_minimal(12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5 ── Global Female vs Male OR forest plot ----------------------------------
# Build 2×2 tables: AE vs rest × sex (all drugs pooled)
orc_totals <- df_cv |> count(sex, name = "total")

ae_or_df <- df_cv |> 
  count(canon_term, sex, name = "a") |> 
  left_join(orc_totals, by = "sex") |> 
  mutate(b = total - a) |> 
  select(-total) |> 
  tidyr::pivot_wider(names_from = sex, values_from = c(a, b), values_fill = 0) |> 
  mutate(
    or   = (a_F * b_M) / (b_F * a_M),
    se   = sqrt(1 / a_F + 1 / b_F + 1 / a_M + 1 / b_M),
    l95  = exp(log(or) - 1.96 * se),
    u95  = exp(log(or) + 1.96 * se)
  )

sig_forest_df <- ae_or_df |> 
  filter((l95 > 1 & or > 1) | (u95 < 1 & or < 1)) |> 
  mutate(sex = factor(ifelse(or > 1, "Female", "Male"), levels = c("Female", "Male")))

forest_plot <- ggplot(sig_forest_df,
                      aes(x = or,
                          y = forcats::fct_reorder(canon_term, or),
                          colour = sex, shape = sex)) +
  geom_pointrange(aes(xmin = l95, xmax = u95),
                  position = position_dodge(width = 0.6),
                  size = 0.9) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 3)) +
  scale_colour_manual(values = c("Female" = "#C51B7D", "Male" = "#0087BD"),
                      labels = c("Female" = "Higher in females", "Male" = "Higher in males"),
                      name = NULL) +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17),
                     labels = c("Female" = "Higher in females", "Male" = "Higher in males"),
                     name = NULL) +
  labs(title = "Sex-specific Odds Ratios for Cardiovascular Adverse Events",
       x = "Odds ratio (Female vs Male, log scale)",
       y = NULL,
       caption = "Only AE categories whose 95 % CI excludes 1 are shown.  Points offset vertically for clarity.") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")
print(forest_plot)
# ──------------------------------------------------------------------------------------------------------

ae_or_df <- df_cv %>%
  filter(sex %in% c("F","M")) %>%
  distinct(case_id, sex, canon_term) %>%                        # report-level
  { base <- .
  totals <- base %>% distinct(case_id, sex) %>% count(sex, name = "total")
  base %>%
    count(canon_term, sex, name = "a") %>%
    left_join(totals, by = "sex") %>%
    mutate(b = pmax(total - a, 0L)) %>%
    select(-total) %>%
    pivot_wider(names_from = sex, values_from = c(a, b), values_fill = 0) %>%
    mutate(
      # Wald OR & CI with 0.5 CC if any zero cell
      cc = if_else(a_F==0 | b_F==0 | a_M==0 | b_M==0, 0.5, 0),
      aF = a_F + cc, bF = b_F + cc, aM = a_M + cc, bM = b_M + cc,
      OR  = (aF * bM) / (bF * aM),
      SE  = sqrt(1/aF + 1/bF + 1/aM + 1/bM),
      L95 = exp(log(OR) - 1.96*SE),
      U95 = exp(log(OR) + 1.96*SE),
      z   = log(OR) / SE,
      p_wald     = 2*pnorm(abs(z), lower.tail = FALSE),
      p_wald_adj = p.adjust(p_wald, method = "BH")
    ) %>%
    # mid-p Fisher p-value on raw integer cells (no CC)
    rowwise() %>%
    mutate(p_fisher_midp = fisher.exact(matrix(c(a_F, b_F, a_M, b_M), 2, byrow = TRUE),
                                        alternative = "two.sided", midp = TRUE)$p.value) %>%
    ungroup() %>%
    mutate(p_fisher_adj = p.adjust(p_fisher_midp, method = "BH"))
  }

# Optional quick comparison of calls between tests
method_summary <- ae_or_df %>%
  transmute(canon_term,
            sig_wald = p_wald_adj < 0.05,
            sig_fisher = p_fisher_adj < 0.05) %>%
  mutate(category = case_when(
    sig_wald & sig_fisher ~ "Both",
    sig_wald & !sig_fisher ~ "Wald only",
    !sig_wald & sig_fisher ~ "Fisher only",
    TRUE ~ "Neither"
  )) %>%
  count(category, name = "N")

print(method_summary)
  
# ──-----------------------------------------------------

# Per-drug, per-sex counts at report level
counts <- df_cv %>%
  filter(sex %in% c("F","M")) %>%
  distinct(case_id, sex, drug_simple, canon_term) %>%
  count(canon_term, drug_simple, sex, name = "a") %>%
  left_join(
    df_cv %>% filter(sex %in% c("F","M")) %>%
      distinct(case_id, sex, drug_simple) %>%
      count(drug_simple, sex, name = "tot"),
    by = c("drug_simple","sex")
  ) %>%
  mutate(a = coalesce(a, 0L), tot = coalesce(tot, 0L))

cmh_results <- counts %>%
  group_by(canon_term) %>% group_split() %>%
  map_df(function(df) {
    # summarise to one row per drug stratum with complete cells
    strata <- df %>%
      group_by(drug_simple) %>%
      summarise(
        a_F  = first(a[sex=="F"],  default = 0L),
        a_M  = first(a[sex=="M"],  default = 0L),
        tot_F= first(tot[sex=="F"],default = 0L),
        tot_M= first(tot[sex=="M"],default = 0L),
        .groups = "drop"
      ) %>%
      mutate(
        b_F = pmax(tot_F - a_F, 0L),
        b_M = pmax(tot_M - a_M, 0L)
      ) %>%
      # drop strata with zero total in both sexes (no info)
      filter(tot_F + tot_M > 0)
    
    if (nrow(strata) == 0) {
      return(tibble(canon_term = df$canon_term[1],
                    or_mh = NA_real_, l95_mh = NA_real_, u95_mh = NA_real_, p_cmh = NA_real_))
    }
    
    # build 2x2xK array safely
    K <- nrow(strata)
    arr <- array(0, dim = c(2, 2, K),
                 dimnames = list(c("F","M"), c("AE","Rest"), NULL))
    for (i in seq_len(K)) {
      arr[,,i] <- matrix(c(strata$a_F[i], strata$b_F[i],
                           strata$a_M[i], strata$b_M[i]),
                         nrow = 2, byrow = TRUE)
    }
    
    mt <- mantelhaen.test(arr, correct = FALSE)
    tibble(
      canon_term = df$canon_term[1],
      or_mh  = unname(mt$estimate),
      l95_mh = unname(mt$conf.int[1]),
      u95_mh = unname(mt$conf.int[2]),
      p_cmh  = mt$p.value
    )
  }) %>%
  mutate(
    p_cmh_adj = p.adjust(p_cmh, "BH"),
    direction = case_when(
      or_mh > 1 ~ "Higher in females",
      or_mh < 1 ~ "Higher in males",
      TRUE ~ "No difference"
    )
  )

cmh_sig <- cmh_results %>%
  filter(!is.na(or_mh), is.finite(or_mh), p_cmh_adj < 0.05) %>%
  arrange(desc(abs(log(or_mh))))

# quick tabulation
cmh_sig %>% count(direction)

# forest (top 25 by |log OR_MH|); OR_MH > 1 = female; < 1 = male
ggplot(cmh_sig %>% slice_head(n = 25),
       aes(or_mh, fct_reorder(canon_term, abs(log(or_mh))),
           color = direction, shape = direction)) +
  geom_pointrange(aes(xmin = l95_mh, xmax = u95_mh), size = 0.9) +
  geom_point(size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 3, 5)) +
  scale_color_manual(values = c("Higher in females"="#C51B7D","Higher in males"="#0087BD")) +
  scale_shape_manual(values = c("Higher in females"=16,"Higher in males"=17)) +
  labs(title = "Adjusted Sex Differences in Cardiovascular Adverse Events",
       subtitle = "Female–male odds ratios adjusted for drug; confidence intervals at 95%. BH-adjusted p < 0.05",
       x = "Common odds ratio (log scale)", y = NULL, color = NULL, shape = NULL) +
  theme_minimal(11) + theme(legend.position = "bottom")

pooled_wald <- df_cv %>%
  dplyr::filter(sex %in% c("F","M")) %>%
  dplyr::distinct(case_id, sex, canon_term) %>%
  dplyr::count(canon_term, sex, name = "a") %>%
  dplyr::left_join(dplyr::count(dplyr::distinct(df_cv, case_id, sex), sex, name="total"), by="sex") %>%
  dplyr::mutate(b = pmax(total - a, 0L)) %>%
  tidyr::pivot_wider(names_from = sex, values_from = c(a,b), values_fill = 0) %>%
  dplyr::mutate(cc = ifelse(a_F==0|b_F==0|a_M==0|b_M==0, 0.5, 0),
                aF=a_F+cc, bF=b_F+cc, aM=a_M+cc, bM=b_M+cc,
                OR=(aF*bM)/(bF*aM),
                SE=sqrt(1/aF+1/bF+1/aM+1/bM),
                p_wald=2*pnorm(abs(log(OR)/SE), lower.tail=FALSE),
                p_wald_adj=p.adjust(p_wald,"BH")) %>%
  dplyr::transmute(canon_term, OR_pooled=OR, p_wald_adj)

# Compare to CMH (drug-adjusted) results you computed
compare <- cmh_results %>%
  dplyr::select(canon_term, or_mh, p_cmh_adj) %>%
  dplyr::left_join(pooled_wald, by="canon_term") %>%
  dplyr::mutate(
    sig_pooled = p_wald_adj < 0.05,
    sig_cmh    = p_cmh_adj  < 0.05,
    category = dplyr::case_when(
      sig_pooled &  sig_cmh ~ "Both",
      sig_pooled & !sig_cmh ~ "Pooled only (drug confounding)",
      !sig_pooled &  sig_cmh ~ "CMH only",
      TRUE                  ~ "Neither"
    )
  )

table(compare$category)
  
# ── Build 2×2 per AE: (Female vs Male) × (AE vs rest), case-based, all drugs pooled
ae_or_df_new <- df_cv %>%
  filter(sex %in% c("F","M")) %>%
  distinct(case_id, sex, canon_term) %>%
  count(canon_term, sex, name = "a") %>%
  left_join(
    df_events_all %>% filter(sex %in% c("F","M")) %>%
      distinct(case_id, sex) %>% count(sex, name = "total"),
    by = "sex"
  ) %>%
  mutate(b = pmax(total - a, 0L)) %>%
  select(-total) %>%
  pivot_wider(names_from = sex, values_from = c(a, b), values_fill = 0) %>%
  # --- WALD log-OR, with 0.5 CC when needed (for OR/SE/CI only)
  mutate(
    cc   = if_else(a_F==0 | b_F==0 | a_M==0 | b_M==0, 0.5, 0),
    aF = a_F + cc, bF = b_F + cc, aM = a_M + cc, bM = b_M + cc,
    or   = (aF * bM) / (bF * aM),
    se   = sqrt(1/aF + 1/bF + 1/aM + 1/bM),
    l95  = exp(log(or) - 1.96 * se),
    u95  = exp(log(or) + 1.96 * se),
    z    = log(or) / se,
    p_wald      = 2 * pnorm(abs(z), lower.tail = FALSE),
    p_wald_adj  = p.adjust(p_wald, method = "BH")
  ) %>%
  # --- mid-p Fisher on the raw integer cells (no continuity correction)
  rowwise() %>%
  mutate(
    p_fisher_midp = fisher.exact(matrix(c(a_F, b_F, a_M, b_M), nrow = 2, byrow = TRUE),
                                 alternative = "two.sided", midp = TRUE)$p.value
  ) %>%
  ungroup() %>%
  mutate(p_fisher_adj = p.adjust(p_fisher_midp, method = "BH"))

# ── Show differences between methods
method_compare <- ae_or_df %>%
  transmute(
    canon_term, or, l95, u95,
    sig_wald   = p_wald_adj   < 0.05,
    sig_fisher = p_fisher_adj < 0.05
  ) %>%
  mutate(category = case_when(
    sig_wald & sig_fisher ~ "Both",
    sig_wald & !sig_fisher ~ "Wald only",
    !sig_wald & sig_fisher ~ "Fisher only",
    TRUE ~ "Neither"
  ))

summary_methods <- method_compare %>% count(category, name = "N")
print(summary_methods)

# ── Prepare plotting data: show significant AEs for each method (top 15 by |log(OR)|)
sig_long <- bind_rows(
  ae_or_df %>% filter(p_wald_adj < 0.05) %>%
    transmute(canon_term, OR, L95, U95, method = "Wald (BH)"),
  ae_or_df %>% filter(p_fisher_adj < 0.05) %>%
    transmute(canon_term, OR, L95, U95, method = "Mid-p Fisher (BH)")
) %>%
  group_by(method) %>%
  slice_max(order_by = abs(log(OR)), n = 15, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(dir = if_else(OR > 1, "Female", "Male"))

# ── Faceted forest: left = Wald, right = mid-p Fisher (same OR/CI; methods differ in p-values)
ggplot(sig_long, aes(x = OR, y = fct_reorder(canon_term, OR), colour = dir, shape = dir)) +
  geom_pointrange(aes(xmin = L95, xmax = U95), position = position_dodge(width = 0.6), size = 0.9) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(breaks = c(0.5, 1, 1.5, 3)) +
  scale_colour_manual(values = c(Female = "#C51B7D", Male = "#0087BD"),
                      labels = c(Female = "Higher in females", Male = "Higher in males"),
                      name = NULL) +
  scale_shape_manual(values = c(Female = 16, Male = 17),
                     labels = c(Female = "Higher in females", Male = "Higher in males"),
                     name = NULL) +
  facet_wrap(~ method, nrow = 1, scales = "free_y") +
  labs(title = "Female vs Male odds ratios for cardiovascular AEs",
       subtitle = "Left: Wald (BH). Right: mid-p Fisher (BH). Forests show OR & Wald 95% CI; methods differ in p-values.",
       x = "Odds ratio (Female / Male, log scale)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")


# ─────────────────────────────────────────────────────────────
# 1) POOLED DISPROPORTIONALITY — case×drug level with all events
#    Requires: df_events_all with columns case_id, drug_simple, sex, canon_term
# ─────────────────────────────────────────────────────────────
# ── INPUT: df_events_all with columns: case_id, drug_simple, canon_term
# One row = one AE occurrence (CV-only already or we filter for canon_term != NA below)

overall_disp <- df_events_all %>%
  filter(!is.na(canon_term)) %>%
  distinct(case_id, drug_simple, canon_term) %>%
  count(drug_simple, canon_term, name = "a") %>%
  left_join(
    df_events_all %>% distinct(case_id, drug_simple) %>%
      count(drug_simple, name = "drug_pairs"),
    by = "drug_simple"
  ) %>%
  left_join(
    df_events_all %>% filter(!is.na(canon_term)) %>%
      distinct(case_id, canon_term) %>%
      count(canon_term, name = "tot_a_canon"),
    by = "canon_term"
  ) %>%
  mutate(
    total_pairs      = nrow(df_events_all %>% distinct(case_id, drug_simple)),
    b                = drug_pairs - a,
    other_drug_pairs = total_pairs - drug_pairs,
    c                = tot_a_canon - a,
    d                = other_drug_pairs - c,
    # <-- FIXED: avoid integer overflow
    cc               = if_else((a == 0L) | (b == 0L) | (c == 0L) | (d == 0L), 0.5, 0.0),
    a2 = a + cc, b2 = b + cc, c2 = c + cc, d2 = d + cc,
    ror     = (a2 * d2) / (b2 * c2),
    se_ror  = sqrt(1/a2 + 1/b2 + 1/c2 + 1/d2),
    ror_l95 = exp(log(ror) - 1.96 * se_ror),
    ror_u95 = exp(log(ror) + 1.96 * se_ror),
    prr     = (a2/(a2 + b2)) / (c2/(c2 + d2)),
    se_prr  = sqrt(1/a2 - 1/(a2 + b2) + 1/c2 - 1/(c2 + d2)),
    prr_l95 = exp(log(prr) - 1.96 * se_prr),
    prr_u95 = exp(log(prr) + 1.96 * se_prr)
  ) %>%
  filter(a >= 5, ror_l95 > 1)

# ── Top-10 forest plot (by ROR)
pooled_top <- overall_disp %>%
  arrange(desc(ror)) %>%
  slice_head(n = 10) %>%
  mutate(
    Label          = paste(drug_simple, canon_term, sep = " – "),
    `ROR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", ror, ror_l95, ror_u95),
    `PRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", prr, prr_l95, prr_u95)
  )

ggplot(pooled_top, aes(ror, fct_reorder(Label, ror))) +
  geom_pointrange(aes(xmin = ror_l95, xmax = ror_u95), linewidth = 0.4) +
  geom_pointrange(aes(xmin = ror_l95, xmax = ror_u95), color = "#0F4D92") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#0F4D92") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10)) +
  labs(title = "Top-10 pooled disproportionality signals",
       x = "Reporting odds ratio (log scale)", y = NULL) +
  theme_minimal(base_size = 11)

# ── Exportable table
pooled_table <- overall_disp %>%
  mutate(
    `N reports (case×drug)` = a,
    `ROR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", ror, ror_l95, ror_u95),
    `PRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", prr, prr_l95, prr_u95)
  ) %>%
  arrange(desc(ror)) %>%
  select(Drug = drug_simple, AE_category = canon_term,
         `N reports (case×drug)`, `ROR (95% CI)`, `PRR (95% CI)`)

write_xlsx(list("Pooled signals" = pooled_table),
           path = "pooled_disproportionality_signals.xlsx")

# ─────────────────────────────────────────────────────────────
# 2) SEX-STRATIFIED DISPROPORTIONALITY — same logic within sex
# ─────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)
library(forcats)

sex_disp_strong <- df_events_all %>%
  filter(sex %in% c("F","M"), !is.na(canon_term)) %>%
  distinct(case_id, sex, drug_simple, canon_term) %>%
  count(sex, drug_simple, canon_term, name = "a") %>%
  left_join(
    df_events_all %>% filter(sex %in% c("F","M")) %>%
      distinct(case_id, sex, drug_simple) %>%
      count(sex, drug_simple, name = "drug_tot"),
    by = c("sex","drug_simple")
  ) %>%
  left_join(
    df_events_all %>% filter(sex %in% c("F","M"), !is.na(canon_term)) %>%
      distinct(case_id, sex, canon_term) %>%
      count(sex, canon_term, name = "ae_tot"),
    by = c("sex","canon_term")
  ) %>%
  left_join(
    df_events_all %>% filter(sex %in% c("F","M")) %>%
      distinct(case_id, sex) %>%
      count(sex, name = "sex_tot"),
    by = "sex"
  ) %>%
  mutate(
    b = drug_tot - a,
    c = ae_tot   - a,
    d = sex_tot  - (a + b + c),
    ror     = (a * d) / (b * c),
    se_ror  = sqrt(1/a + 1/b + 1/c + 1/d),
    ror_l95 = exp(log(ror) - 1.96 * se_ror),
    ror_u95 = exp(log(ror) + 1.96 * se_ror),
    prr     = (a/(a+b)) / (c/(c+d)),
    se_prr  = sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)),
    prr_l95 = exp(log(prr) - 1.96 * se_prr),
    prr_u95 = exp(log(prr) + 1.96 * se_prr)
  ) %>%
  filter(a >= 5, ror_l95 > 1)

strong_sex_signals_tbl <- sex_disp_strong %>%
  arrange(desc(ror)) %>%
  mutate(
    Sex = ifelse(sex == "F", "Female", "Male"),
    Drug = drug_simple, AE_category = canon_term,
    `N reports` = a,
    `ROR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", ror, ror_l95, ror_u95),
    `PRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", prr, prr_l95, prr_u95)
  ) %>%
  select(Sex, Drug, AE_category, `N reports`, `ROR (95% CI)`, `PRR (95% CI)`)

# Top-10 plot
sex_disp_strong %>%
  arrange(desc(ror)) %>% slice_head(n = 10) %>%
  mutate(Label = paste(ifelse(sex=="F","Female","Male"), drug_simple, canon_term, sep = " – ")) %>%
  ggplot(aes(ror, fct_reorder(Label, ror), colour = sex)) +
  geom_pointrange(aes(xmin = ror_l95, xmax = ror_u95)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(breaks = c(0.5,1,2,5,10,20)) +
  scale_colour_manual(values = c(F="#C51B7D", M="#0087BD"),
                      labels = c(F="Female", M="Male")) +
  labs(title = "Top-10 sex-specific disproportionality signals",
       x = "Reporting odds ratio (log scale)", y = NULL, colour = "Sex") +
  theme_minimal(11)


# Export
write_xlsx(list("Sex signals (F vs M)" = sex_table),
           "sex_disproportionality_signals_F_over_M.xlsx")

# ─────────────────────────────────────────────────────────────
# 3. RATIO OF RORs  (female vs male)  ── sig_ratio
# ─────────────────────────────────────────────────────────────
# Female vs Male ratio of ROR with continuity correction
ror_ratio_df <- sex_disp_full %>%
  select(sex, drug_simple, canon_term, a, b, c, d) %>%
  pivot_wider(names_from = sex, values_from = c(a, b, c, d), values_fill = 0) %>%
  mutate(
    # Haldane–Anscombe 0.5 if any zero in the 2×2 per sex
    ccF = if_else(a_F == 0 | b_F == 0 | c_F == 0 | d_F == 0, 0.5, 0),
    ccM = if_else(a_M == 0 | b_M == 0 | c_M == 0 | d_M == 0, 0.5, 0),
    aF = a_F + ccF, bF = b_F + ccF, cF = c_F + ccF, dF = d_F + ccF,
    aM = a_M + ccM, bM = b_M + ccM, cM = c_M + ccM, dM = d_M + ccM,
    ror_F    = (aF * dF) / (bF * cF),
    ror_M    = (aM * dM) / (bM * cM),
    log_diff = log(ror_F) - log(ror_M),
    se_diff  = sqrt(1/aF + 1/bF + 1/cF + 1/dF + 1/aM + 1/bM + 1/cM + 1/dM),
    ratio_ror = exp(log_diff),
    l95 = exp(log_diff - 1.96 * se_diff),
    u95 = exp(log_diff + 1.96 * se_diff),
    p   = 2 * pnorm(abs(log_diff / se_diff), lower.tail = FALSE),
    p_adj = p.adjust(p, method = "BH")
  )

# Significant sex differences (BH-adjusted), require some signal evidence
sig_ratio <- ror_ratio_df %>%
  filter(p_adj < 0.05, (a_F + a_M) >= 5) %>%
  arrange(desc(abs(log(ratio_ror)))) %>%
  mutate(label = paste(drug_simple, canon_term, sep = " – "))

# Plot
ggplot(sig_ratio, aes(ratio_ror, fct_reorder(label, ratio_ror))) +
  geom_pointrange(aes(xmin = l95, xmax = u95), colour = "darkorange") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(breaks = c(0.25, 0.5, 1, 2, 4)) +
  labs(title = "Sex differential (female / male) in disproportionality",
       x = "Female / Male ROR", y = NULL) +
  theme_minimal(11)



# Apply your “strong” filter (mirror pooled criteria; tweak if desired)
sex_disp_strong <- sex_disp %>%
  filter(a >= 5, ror_l95 > 1)

# ── 1. Strong pooled signals ----------------------------------------------
pooled_signals <- overall_disp %>%                 # already strong-filtered earlier
  mutate(signal_status = "Detected in pooled")     # keep all columns

# ── 2. Strong stratified signals (female / male) ---------------------------
female_signals <- sex_disp_strong %>% 
  filter(sex == "F") %>% 
  mutate(sex_specific = "Female")

male_signals   <- sex_disp_strong %>% 
  filter(sex == "M") %>% 
  mutate(sex_specific = "Male")

sex_signals <- bind_rows(female_signals, male_signals)

# ── 3. Hidden signals = sex-specific but NOT in pooled ---------------------
hidden_signals <- anti_join(
  sex_signals,
  pooled_signals %>% select(drug_simple, canon_term),   # just keys
  by = c("drug_simple", "canon_term")
) %>% 
  mutate(signal_status = "Hidden sex-specific signal")

# ── 4. Quick summary table (optional) --------------------------------------
summary_concord <- tibble::tibble(
  Category = c("Pooled only", "Sex-specific only", "Both"),
  N = c(
    nrow(anti_join(pooled_signals, sex_signals,
                   by = c("drug_simple","canon_term"))),
    nrow(hidden_signals),
    nrow(inner_join(pooled_signals, sex_signals,
                    by = c("drug_simple","canon_term")))
  )
)

print(summary_concord)
print(hidden_signals)

# ── Export pretty table of hidden signals ----------------------------------
tidy_hidden <- hidden_signals %>% 
  mutate(
    Sex         = ifelse(sex == "F", "Female", "Male"),
    Drug        = drug_simple,
    AE_category = canon_term,
    `N reports` = a,
    `ROR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", ror,  ror_l95,  ror_u95),
    `PRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", prr,  prr_l95,  prr_u95)
  ) %>% 
  select(Sex, Drug, AE_category, `N reports`, `ROR (95% CI)`, `PRR (95% CI)`)

write_xlsx(
  list("Hidden sex-specific signals" = tidy_hidden),
  path = "hidden_sex_specific_signals.xlsx"
)

# ── 1. Strong pooled signals ----------------------------------------------
pooled_signals <- overall_disp %>%                       # already strong-filtered
  mutate(signal_status = "Detected in pooled")           # keep all columns

# ── 2. Strong stratified signals (female / male) ---------------------------
female_signals <- sex_disp_strong %>% 
  filter(sex == "F") %>% 
  mutate(sex_specific = "Female")

male_signals   <- sex_disp_strong %>% 
  filter(sex == "M") %>% 
  mutate(sex_specific = "Male")

sex_signals <- bind_rows(female_signals, male_signals)   # combine

# ── 3. Hidden signals = sex-specific but NOT in pooled ---------------------
hidden_signals <- anti_join(
  sex_signals,
  pooled_signals %>% select(drug_simple, canon_term),     # just keys
  by = c("drug_simple", "canon_term")
) %>% 
  mutate(signal_status = "Hidden sex-specific signal")

# ── 4. Quick summary table (optional) --------------------------------------
summary_concord <- tibble::tibble(
  Category = c("Pooled only", "Sex-specific only", "Both"),
  N = c(
    nrow(anti_join(pooled_signals, sex_signals,
                   by = c("drug_simple","canon_term"))),
    nrow(hidden_signals),
    nrow(inner_join(pooled_signals, sex_signals,
                    by = c("drug_simple","canon_term")))
  )
)

print(summary_concord)
print(hidden_signals)

tidy_hidden <- hidden_signals %>% 
  mutate(
    Sex         = ifelse(sex == "F", "Female", "Male"),
    Drug        = drug_simple,
    AE_category = canon_term,
    `N reports` = a,
    `ROR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", ror,  ror_l95,  ror_u95),
    `PRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", prr,  prr_l95,  prr_u95)
  ) %>% 
  select(Sex, Drug, AE_category, `N reports`, `ROR (95% CI)`, `PRR (95% CI)`)
library(writexl)
write_xlsx(
  list("Hidden sex-specific signals" = tidy_hidden),
  path = "hidden_sex_specific_signals.xlsx"
)


hidden_signals <- sex_disp_strong %>%
  anti_join(
    overall_disp %>% distinct(drug_simple, canon_term),
    by = c("drug_simple", "canon_term")
  )

# Tidy view for export
tidy_hidden <- hidden_signals %>%
  mutate(
    Sex = case_when(
      sex %in% c("F","Female","FEMALE") ~ "Female",
      sex %in% c("M","Male","MALE")     ~ "Male",
      TRUE ~ as.character(sex)
    ),
    Drug        = drug_simple,
    AE_category = canon_term,
    `N reports` = a,
    `ROR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", ror, ror_l95, ror_u95),
    `PRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", prr, prr_l95, prr_u95)
  ) %>%
  select(Sex, Drug, AE_category, `N reports`, `ROR (95% CI)`, `PRR (95% CI)`) %>%
  arrange(desc(`N reports`), Drug, AE_category)

# Inspect in console
print(tidy_hidden, n = 30)

# Export
write_xlsx(
  list("Hidden sex-specific signals" = tidy_hidden),
  path = "hidden_sex_specific_signals.xlsx"
)


# ─────────────────────────────────────────────────────────────
# 9. Word cloud of canonical terms (top 40)
# ─────────────────────────────────────────────────────────────

top_n_cloud <- 40
wc_data <- ae_mapped |>
  count(canon_term, sort = TRUE) |>
  slice_head(n = top_n_cloud)

wc_plot <- ggplot(wc_data, aes(label = canon_term, size = n)) +
  geom_text_wordcloud(eccentricity = 0.6, rm_outside = TRUE) +
  scale_size_area(max_size = 18) +
  theme_minimal(base_size = 14) +
  labs(title = paste0("Top Cardiovascular Advers Events Categories"),
       size  = "Reports")
print(wc_plot)



# ─────────────────────────────────────────────────────────────
# 1. POOLED DISPROPORTIONALITY  ── overall_disp
# ─────────────────────────────────────────────────────────────
overall_disp <- df_cv %>%                                   # df_cv = case×drug×event
  count(drug_simple, canon_term, name = "a") %>%            # cell a
  left_join(count(df_cv, drug_simple,  name = "drug_tot"), by = "drug_simple") %>%
  left_join(count(df_cv, canon_term,  name = "ae_tot"),   by = "canon_term") %>%
  mutate(
    b = drug_tot - a,
    c = ae_tot   - a,
    d = nrow(df_cv) - (a + b + c),
    ror     = (a * d) / (b * c),
    se_ror  = sqrt(1/a + 1/b + 1/c + 1/d),
    ror_l95 = exp(log(ror) - 1.96 * se_ror),
    ror_u95 = exp(log(ror) + 1.96 * se_ror),
    prr     = (a/(a+b)) / (c/(c+d)),
    se_prr  = sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)),
    prr_l95 = exp(log(prr) - 1.96 * se_prr),
    prr_u95 = exp(log(prr) + 1.96 * se_prr)
  ) %>%
  filter(a >= 5, ror_l95 > 1)

# ── Top-10 forest plot & export table
pooled_top <- overall_disp %>%
  arrange(desc(ror)) %>% slice_head(n = 10) %>%
  mutate(Label          = paste(drug_simple, canon_term, sep = " – "),
         `ROR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", ror, ror_l95, ror_u95),
         `PRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", prr, prr_l95, prr_u95))

ggplot(pooled_top,
       aes(ror, fct_reorder(Label, ror))) +
  geom_pointrange(aes(xmin = ror_l95, xmax = ror_u95),
                  colour = "#0F4D92") +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "#0F4D92") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10)) +
  labs(title = "Top-10 pooled disproportionality signals",
       x = "Reporting odds ratio", y = NULL) +
  theme_minimal(11)

pooled_table <- overall_disp %>%
  mutate(`N reports`   = a,
         `ROR (95% CI)`= sprintf("%.2f (%.2f–%.2f)", ror, ror_l95, ror_u95),
         `PRR (95% CI)`= sprintf("%.2f (%.2f–%.2f)", prr, prr_l95, prr_u95)) %>%
  arrange(desc(ror)) %>%
  select(Drug = drug_simple, AE_category = canon_term,
         `N reports`, `ROR (95% CI)`, `PRR (95% CI)`)
library(writexl)
write_xlsx(list("Pooled signals" = pooled_table),
           "pooled_disproportionality_signals.xlsx")
getwd()
# ─────────────────────────────────────────────────────────────
# 2. SEX-STRATIFIED DISPROPORTIONALITY  ── within_sex_disp
# ─────────────────────────────────────────────────────────────
# --- stratified build mimicking overall_disp ------------------------------
library(dplyr)
library(ggplot2)
library(forcats)
library(writexl)

# case × drug × event (CV-only)
df_cv <- df_events_all %>%
  filter(sex %in% c("F","M"), !is.na(canon_term)) %>%
  distinct(case_id, sex, drug_simple, canon_term)

# Sex-stratified disproportionality within sex (mirrors your formulae)
sex_disp_strong <- df_cv %>%
  count(sex, drug_simple, canon_term, name = "a") %>%                                   # a
  left_join(count(df_cv, sex, drug_simple, name = "drug_tot"), by = c("sex","drug_simple")) %>%
  left_join(count(df_cv, sex, canon_term,  name = "ae_tot"),   by = c("sex","canon_term")) %>%
  left_join(count(df_cv, sex,              name = "sex_tot"),  by = "sex") %>%
  mutate(
    b = drug_tot - a,
    c = ae_tot   - a,
    d = sex_tot  - (a + b + c),
    ror     = (a * d) / (b * c),
    se_ror  = sqrt(1/a + 1/b + 1/c + 1/d),
    ror_l95 = exp(log(ror) - 1.96 * se_ror),
    ror_u95 = exp(log(ror) + 1.96 * se_ror),
    prr     = (a/(a+b)) / (c/(c+d)),
    se_prr  = sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)),
    prr_l95 = exp(log(prr) - 1.96 * se_prr),
    prr_u95 = exp(log(prr) + 1.96 * se_prr)
  ) %>%
  filter(a >= 5, ror_l95 > 1)

# Tidy table (ready to export)
strong_sex_signals_tbl <- sex_disp_strong %>%
  arrange(desc(ror)) %>%
  mutate(
    Sex = ifelse(sex == "F", "Female", "Male"),
    Drug = drug_simple, AE_category = canon_term,
    `N reports` = a,
    `ROR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", ror, ror_l95, ror_u95),
    `PRR (95% CI)` = sprintf("%.2f (%.2f–%.2f)", prr, prr_l95, prr_u95)
  ) %>%
  select(Sex, Drug, AE_category, `N reports`, `ROR (95% CI)`, `PRR (95% CI)`)

# writexl::write_xlsx(strong_sex_signals_tbl, "strong_sex_disproportionality_signals.xlsx")

# Top-10 plot
sig_sex_top <- sex_disp_strong %>% arrange(desc(ror)) %>% slice_head(n = 10) %>%
  mutate(Label = paste(ifelse(sex=="F","Female","Male"), drug_simple, canon_term, sep = " – "))

ggplot(sig_sex_top, aes(ror, fct_reorder(Label, ror), colour = sex)) +
  geom_pointrange(aes(xmin = ror_l95, xmax = ror_u95)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10(breaks = c(0.5,1,2,5,10,20)) +
  scale_colour_manual(values = c(F="#C51B7D", M="#0087BD"),
                      labels = c(F="Female", M="Male")) +
  labs(title = "Top-10 sex-specific disproportionality signals",
       x = "Reporting odds ratio (log scale)", y = NULL, colour = "Sex") +
  theme_minimal(11)

 
