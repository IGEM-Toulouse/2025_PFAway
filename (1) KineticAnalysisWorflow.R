# ================================================================
# Growth-OD pipeline — preprocessing, kinetics, time-normalization
# Dependencies: tidyverse, readxl, writexl, gcplyr
# ================================================================

# ---------------------- USER PARAMETERS -------------------------
# (1) Input file (CSV or XLSX) and, if XLSX, the sheet name
input_path   <- "path/to/data.xlsx"   # e.g. "data/growth_curves.csv"
input_sheet  <- "Sheet"                        # ignored for CSV

# (2) Column names in your input table
# Required columns (exact names or adapt here):
#   Time (in hours), OD, Date, Media, Strain, Molecule, Concentration, RepBio, RepTech
time_col <- "Time"
od_col   <- "OD"
group_cols <- c("Date","Media","Strain","Molecule","Concentration","RepBio","RepTech")

# (3) Time normalization step (in hours) 0.279722222 h ≈ 16.7833 min
step_h <- 0.279722222

# (4) OD blank (offset for gcplyr, usually 0 if already blank-corrected)
blank_od <- 0

# (5) Output directory
output_dir <- "path/to/output/file/"

# ------------------------- LIBRARIES ----------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(writexl)
  library(gcplyr)
})

# ---------------------- UTILITAIRES -----------------------------
stop_if_missing <- function(df, cols) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
}

read_any <- function(path, sheet = NULL) {
  if (grepl("\\.xlsx?$", path, ignore.case = TRUE)) {
    readxl::read_excel(path, sheet = sheet)
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    readr::read_csv(path, show_col_types = FALSE)
  } else {
    stop("Unsupported file type. Use .csv or .xlsx.")
  }
}


soft_clean_names <- function(x) {
  x <- gsub("\\s+", "_", x)          
  x <- gsub("[^A-Za-z0-9_]", "", x)   
  x
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------- LOAD --------------------------------
Kinetics <- read_any(input_path, sheet = input_sheet)
names(Kinetics) <- soft_clean_names(names(Kinetics)) 

name_map <- function(df, want) {
  current <- names(df)
  out <- want
  for (i in seq_along(want)) {
    hit <- which(tolower(current) == tolower(want[i]))
    if (length(hit) == 1) out[i] <- current[hit]
  }
  setNames(out, want)
}

nm <- name_map(Kinetics, c(time_col, od_col, group_cols))

stop_if_missing(Kinetics, unname(nm))

Kinetics <- Kinetics %>%
  rename(
    Time = all_of(nm[[time_col]]),
    OD   = all_of(nm[[od_col]])
  )

for (g in group_cols) {
  if (g != nm[[g]]) {
    Kinetics <- rename(Kinetics, !!g := all_of(nm[[g]]))
  }
}

Kinetics <- Kinetics %>%
  mutate(
    Time = as.numeric(Time),
    OD   = as.numeric(OD)
  )

Kinetics <- Kinetics %>%
  mutate(
    Condition = paste(Media, Molecule),
    ref = paste(Date, Media, Strain, Molecule, Concentration, RepBio, RepTech, sep = "_")
  )

# ----------------- GROWTH RATE & LAG CALCULATION -----------------
Kinetics <- Kinetics %>%
  arrange(ref, Time) %>%
  group_by(across(all_of(c(group_cols, "Condition", "ref")))) %>%
  mutate(
    mu7 = calc_deriv(
      y = OD, x = Time,
      percapita = TRUE, blank = blank_od,
      window_width_n = 7, trans_y = "log"
    )
  ) %>%
  ungroup()

get_y0 <- function(t, y) {
  if (any(t == 0, na.rm = TRUE)) y[which(t == 0)[1]] else y[which.min(abs(t - min(t, na.rm = TRUE)))]
}

Synthesis <- Kinetics %>%
  arrange(ref, Time) %>%
  group_by(across(all_of(c(group_cols, "Condition", "ref")))) %>%
  summarise(
    lag7   = { t <- Time; y <- OD; mu <- mu7; lag_time(y = y, x = t, deriv = mu, y0 = get_y0(t, y), blank = blank_od) },
    muMAX7 = suppressWarnings(max_gc(mu7)),
    .groups = "drop_last"
  ) %>%
  ungroup()

SynthesisTech <- Synthesis %>%
  group_by(across(all_of(c("Date","Media","Strain","Molecule","Concentration","RepBio","Condition")))) %>%
  summarise(
    lag7   = mean(lag7,   na.rm = TRUE),
    muMAX7 = mean(muMAX7, na.rm = TRUE),
    .groups = "drop"
  )

SynthesisBio <- Synthesis %>%
  group_by(across(all_of(c("Media","Strain","Molecule","Concentration","Condition")))) %>%
  summarise(
    lag7   = mean(lag7,   na.rm = TRUE),
    muMAX7 = mean(muMAX7, na.rm = TRUE),
    .groups = "drop"
  )

# ---------------------- TIME NORMALIZATION -----------------------
per_rep <- Kinetics %>%
  group_by(across(all_of(c(group_cols, "Condition", "ref")))) %>%
  arrange(Time, .by_group = TRUE) %>%
  mutate(Time0 = Time - first(Time)) %>%
  summarise(t_max = max(Time0, na.rm = TRUE), .groups = "drop")

overlap <- per_rep %>%
  group_by(across(all_of(c("Media","Strain","Molecule","Concentration","Condition")))) %>%
  summarise(t_end = min(t_max, na.rm = TRUE), .groups = "drop")

Kinetics_interp <- Kinetics %>%
  group_by(across(all_of(c(group_cols, "Condition", "ref")))) %>%
  arrange(Time, .by_group = TRUE) %>%
  mutate(Time0 = Time - first(Time)) %>%
  ungroup() %>%
  inner_join(overlap, by = c("Media","Strain","Molecule","Concentration","Condition")) %>%
  group_by(across(all_of(c(group_cols, "Condition", "ref")))) %>%
  group_modify(~{
    t_end <- unique(.x$t_end)
    grid  <- seq(0, t_end, by = step_h)
    y     <- approx(x = .x$Time0, y = .x$OD, xout = grid, method = "linear", rule = 1)$y
    tibble(Time = grid, OD = y)
  }) %>%
  ungroup()

Kinetics_mean <- Kinetics_interp %>%
  group_by(Media, Strain, Molecule, Concentration, Condition, Time) %>%
  summarise(
    OD_m  = mean(OD, na.rm = TRUE),
    OD_sd = sd(OD,   na.rm = TRUE),
    n     = dplyr::n(),
    .groups = "drop"
  )

# --------------------------- SAVE --------------------------------
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write_xlsx(Kinetics,        file.path(output_dir, "Kinetics_raw.xlsx"))
write_xlsx(Kinetics_interp, file.path(output_dir, "Kinetics_norm.xlsx"))
write_xlsx(Kinetics_mean,   file.path(output_dir, "Kinetics_norm_mean.xlsx"))
write_xlsx(Synthesis,     file.path(output_dir, "gr_lag_by_replicate.xlsx"))
write_xlsx(SynthesisBio,  file.path(output_dir, "gr_lag_bio_mean.xlsx"))

message("Done. Files written to: ", normalizePath(output_dir))

