# app.R — Kinetic Viewer (generic / “nude”)
suppressPackageStartupMessages({
    library(shiny); library(shinydashboard); library(ggplot2)
    library(dplyr); library(tidyr); library(readxl); library(readr); library(stringr)
})

# ========================== USER SETTINGS ==========================
data_dir <- "data"
files <- list(
    interp     = "Kinetics_norm.xlsx",        # Time, OD, Condition, Concentration, (Strain), ref, RepBio, RepTech
    interp_sum = "gr_lag_by_replicate.xlsx",# lag7, muMAX7 per replicate (+ Condition, Concentration, Strain)
    mean       = "Kinetics_norm_mean.xlsx",    # Time, OD_m, Condition, Concentration, (Strain)
    mean_sum   = "gr_lag_bio_mean.xlsx"     # lag7, muMAX7 (bio mean) (+ Condition, Concentration, Strain)
)

# ============================ HELPERS ==============================
`%||%` <- function(a, b) if (!is.null(a)) a else b

safe_read_xlsx <- function(path) if (file.exists(path)) readxl::read_xlsx(path) else tibble::tibble()

ensure_cols <- function(df, cols) { miss <- setdiff(cols, names(df)); if (length(miss)) df[miss] <- NA; df }
filter_by_condition <- function(df, cond) {
    if (is.null(cond) || !length(cond) || !("Condition" %in% names(df))) return(df)
    df %>% mutate(Condition = trimws(as.character(Condition))) %>% filter(Condition == trimws(as.character(cond)))
}
ordered_conc_levels <- function(...) {
    x <- c(...); x <- as.character(x); if (!length(x)) return(character())
    nums <- suppressWarnings(readr::parse_number(x)); ord <- order(nums, x, na.last = TRUE); unique(x[ord])
}
apply_conc_order <- function(df, levels_vec) {
    if (!("Concentration" %in% names(df))) return(df)
    df %>% mutate(Concentration = factor(as.character(Concentration), levels = levels_vec))
}
build_conc_palette <- function(levels, from = "#BBDFF6", to = "#1F78B4", mid = NULL, reverse = FALSE) {
    lev <- as.character(levels); if (!length(lev)) return(NULL)
    ramp <- if (is.null(mid)) grDevices::colorRampPalette(c(from, to)) else grDevices::colorRampPalette(c(from, mid, to))
    cols <- ramp(length(lev)); if (reverse) cols <- rev(cols); stats::setNames(cols, lev)
}
pick_metric_col <- function(d1, d2, candidates) {
    common <- intersect(names(d1), names(d2)); for (c in candidates) if (c %in% common) return(c)
    alln <- union(names(d1), names(d2));       for (c in candidates) if (c %in% alln)  return(c)
    candidates[[1]]
}

# ============================== DATA ===============================
ds <- list(
    interp     = safe_read_xlsx(file.path(data_dir, files$interp)),
    interp_sum = safe_read_xlsx(file.path(data_dir, files$interp_sum)),
    mean       = safe_read_xlsx(file.path(data_dir, files$mean)),
    mean_sum   = safe_read_xlsx(file.path(data_dir, files$mean_sum))
)

get_conditions <- function() {
    bind_rows(
        ds$interp %>% dplyr::select(any_of("Condition")),
        ds$mean   %>% dplyr::select(any_of("Condition"))
    ) %>% dplyr::distinct() %>% dplyr::pull(Condition) %>% as.character() %>% sort()
}

# ============================= THEME ===============================
plot_theme_base <- theme_classic(base_size = 14) +
    theme(strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          legend.position = "bottom")

LABEL_KIN_X  <- "Time (h)"
LABEL_KIN_Y  <- expression(OD[600*nm])
LABEL_MU_Y   <- expression(Growth~Rate~(h^{-1}))
LABEL_LAG_Y  <- "Lag Time (h)"

KIN_LINE_REP_ALPHA <- 0.30
KIN_LINE_REP_SIZE  <- 0.5
KIN_LINE_MEAN_SIZE <- 1.0
PT_REP_ALPHA <- 0.35
PT_REP_SIZE  <- 1.8
PT_MEAN_SIZE <- 3.0

# ============================== PLOTS ==============================
# Kinetic (generic) — two display modes
kinetic_plot_generic <- function(d_interp, d_mean, mode = c("facetConc_colorStrain","facetStrain_colorConc")) {
    mode <- match.arg(mode)
    d_interp <- ensure_cols(d_interp, c("Time","OD","Condition","Concentration","ref","RepBio","RepTech","Strain"))
    d_mean   <- ensure_cols(d_mean,   c("Time","OD_m","OD","Condition","Concentration","Strain"))
    
    validate(need(nrow(d_interp) + nrow(d_mean) > 0, "No data to plot (Suivi_norm*)."))
    y_mean <- if ("OD_m" %in% names(d_mean)) "OD_m" else if ("OD" %in% names(d_mean)) "OD" else NULL
    validate(need(!is.null(y_mean), "Mean table needs 'OD_m' or 'OD'."))
    
    conc_levels <- ordered_conc_levels(d_interp$Concentration, d_mean$Concentration)
    d_interp <- apply_conc_order(d_interp, conc_levels)
    d_mean   <- apply_conc_order(d_mean,   conc_levels)
    
    # choose facet/color variables
    if (mode == "facetConc_colorStrain") { facet_var <- "Concentration"; color_var <- "Strain"; show_reps <- TRUE
    } else                                { facet_var <- "Strain";        color_var <- "Concentration"; show_reps <- FALSE }
    
    # palettes
    use_manual_palette <- FALSE; pal_vals <- NULL
    if (color_var == "Concentration") {
        pal_vals <- build_conc_palette(conc_levels); use_manual_palette <- TRUE
    } else if (color_var == "Strain") {
        lev <- sort(unique(as.character(c(d_interp$Strain, d_mean$Strain))))
        lev <- lev[nzchar(lev)]
        if (!length(lev)) { color_var <- "Concentration"; pal_vals <- build_conc_palette(conc_levels); use_manual_palette <- TRUE }
    }
    
    p <- ggplot()
    if (show_reps && nrow(d_interp)) {
        p <- p + geom_line(
            data = d_interp,
            aes(x = Time, y = OD, color = .data[[color_var]],
                group = interaction(Condition, Strain, Concentration, ref, RepBio, RepTech)),
            alpha = KIN_LINE_REP_ALPHA, linewidth = KIN_LINE_REP_SIZE
        )
    }
    if (nrow(d_mean)) {
        p <- p + geom_line(
            data = d_mean,
            aes(x = Time, y = .data[[y_mean]], color = .data[[color_var]],
                group = interaction(Condition, Strain, Concentration)),
            linewidth = KIN_LINE_MEAN_SIZE
        )
    }
    p <- p + labs(x = LABEL_KIN_X, y = LABEL_KIN_Y, color = color_var) + plot_theme_base
    if (use_manual_palette && !is.null(pal_vals)) p <- p + scale_color_manual(values = pal_vals, drop = FALSE, name = color_var)
    
    if (facet_var %in% names(d_mean) || facet_var %in% names(d_interp)) {
        p <- p + facet_wrap(reformulate(facet_var), ncol = 3, drop = FALSE)
    }
    p
}

# Points (generic) — Growth / Lag
points_plot_generic <- function(d_rep, d_mean, y_candidates, y_lab) {
    y_col <- pick_metric_col(d_rep, d_mean, y_candidates)
    
    d_rep  <- ensure_cols(d_rep,  c("Condition","Concentration","Strain", y_col))
    d_mean <- ensure_cols(d_mean, c("Condition","Concentration","Strain", y_col))
    # numeric
    if (!is.numeric(d_rep[[y_col]]))  d_rep[[y_col]]  <- suppressWarnings(readr::parse_number(as.character(d_rep[[y_col]])))
    if (!is.numeric(d_mean[[y_col]])) d_mean[[y_col]] <- suppressWarnings(readr::parse_number(as.character(d_mean[[y_col]])))
    
    validate(need(sum(is.finite(d_rep[[y_col]]), na.rm = TRUE) + sum(is.finite(d_mean[[y_col]]), na.rm = TRUE) > 0,
                  "No numeric values to plot."))
    
    conc_levels <- ordered_conc_levels(d_rep$Concentration, d_mean$Concentration)
    d_rep  <- apply_conc_order(d_rep, conc_levels)
    d_mean <- apply_conc_order(d_mean, conc_levels)
    
    # color by Strain if present, else by Concentration
    color_var <- if (any(!is.na(c(d_rep$Strain, d_mean$Strain)))) "Strain" else "Concentration"
    pal_vals <- if (color_var == "Concentration") build_conc_palette(conc_levels) else NULL
    
    pd_rep  <- position_jitterdodge(jitter.width = 0.08, jitter.height = 0, dodge.width = 0.5, seed = 1)
    pd_mean <- position_dodge(width = 0.5)
    
    p <- ggplot()
    if (nrow(d_rep))
        p <- p + geom_point(
            data = d_rep %>% filter(is.finite(.data[[y_col]])),
            aes(x = Concentration, y = .data[[y_col]], color = .data[[color_var]], group = .data[[color_var]]),
            alpha = PT_REP_ALPHA, size = PT_REP_SIZE, position = pd_rep
        )
    if (nrow(d_mean))
        p <- p + geom_point(
            data = d_mean %>% filter(is.finite(.data[[y_col]])),
            aes(x = Concentration, y = .data[[y_col]], color = .data[[color_var]], group = .data[[color_var]]),
            size = PT_MEAN_SIZE, position = pd_mean
        )
    
    p <- p + scale_x_discrete(drop = FALSE) + labs(x = "Concentration", y = y_lab, color = color_var) + plot_theme_base
    if (!is.null(pal_vals)) p <- p + scale_color_manual(values = pal_vals, drop = FALSE, name = color_var)
    p
}

# ================================ UI =================================
header <- dashboardHeader(title = "Kinetic Viewer")

sidebar <- dashboardSidebar(
    width = 260,
    sidebarMenu(
        id = "sidebar",
        menuItem("Kinetic",     tabName = "kinetic", icon = icon("wave-square")),
        menuItem("Growth Rate", tabName = "mu",      icon = icon("chart-line")),
        menuItem("Lag Time",    tabName = "lag",     icon = icon("clock"))
    )
)

body <- dashboardBody(
    tabItems(
        tabItem(tabName = "kinetic",
                fluidRow(
                    box(title = "Options", width = 3, solidHeader = TRUE,
                        selectizeInput("cond", "Condition", choices = NULL, multiple = FALSE,
                                       options = list(placeholder = "Choose a condition…")),
                        radioButtons("view_mode", "Representation",
                                     choices = c("(A) Facet by concentration, color by strain" = "facetConc_colorStrain",
                                                 "(B) Facet by strain (means only), color by concentration" = "facetStrain_colorConc"),
                                     selected = "facetConc_colorStrain"),
                        helpText(HTML("Thick line = biological mean<br>Thin line = replicate (only in mode A)"))
                    ),
                    box(title = "Kinetic", width = 9, solidHeader = TRUE,
                        plotOutput("plot_kinetic", height = 620))
                )
        ),
        tabItem(tabName = "mu",
                fluidRow(
                    box(title = "Options", width = 3, solidHeader = TRUE,
                        selectizeInput("cond_gr", "Condition", choices = NULL, multiple = FALSE,
                                       options = list(placeholder = "Choose a condition…")),
                        helpText(HTML("Opaque point = biological mean<br>Transparent point = replicate"))
                    ),
                    box(title = "Growth Rate", width = 9, solidHeader = TRUE,
                        plotOutput("plot_mu", height = 560))
                )
        ),
        tabItem(tabName = "lag",
                fluidRow(
                    box(title = "Options", width = 3, solidHeader = TRUE,
                        selectizeInput("cond_lag", "Condition", choices = NULL, multiple = FALSE,
                                       options = list(placeholder = "Choose a condition…")),
                        helpText(HTML("Opaque point = biological mean<br>Transparent point = replicate"))
                    ),
                    box(title = "Lag Time", width = 9, solidHeader = TRUE,
                        plotOutput("plot_lag", height = 560))
                )
        )
    )
)

ui <- dashboardPage(skin = "black", header, sidebar, body)

# ================================= SERVER ===============================
server <- function(input, output, session) {
    observe({
        conds <- get_conditions()
        updateSelectizeInput(session, "cond",     choices = conds, selected = conds[1] %||% NULL, server = TRUE)
        updateSelectizeInput(session, "cond_gr",  choices = conds, selected = conds[1] %||% NULL, server = TRUE)
        updateSelectizeInput(session, "cond_lag", choices = conds, selected = conds[1] %||% NULL, server = TRUE)
    })
    
    # Kinetic
    output$plot_kinetic <- renderPlot({
        d_interp <- ds$interp; d_mean <- ds$mean
        if (isTruthy(input$cond)) {
            d_interp <- filter_by_condition(d_interp, input$cond)
            d_mean   <- filter_by_condition(d_mean,   input$cond)
        }
        kinetic_plot_generic(d_interp, d_mean, mode = input$view_mode %||% "facetConc_colorStrain")
    }, res = 96) %>% bindCache(input$cond, input$view_mode)
    
    # Growth rate
    output$plot_mu <- renderPlot({
        d_rep  <- ds$interp_sum; d_mean <- ds$mean_sum
        if (isTruthy(input$cond_gr)) {
            d_rep  <- filter_by_condition(d_rep,  input$cond_gr)
            d_mean <- filter_by_condition(d_mean, input$cond_gr)
        }
        points_plot_generic(
            d_rep, d_mean,
            y_candidates = c("muMAX7","muMAX_m3","GrowthRate","mu","mu_max","muMAX"),
            y_lab = LABEL_MU_Y
        )
    }, res = 96) %>% bindCache(input$cond_gr)
    
    # Lag time
    output$plot_lag <- renderPlot({
        d_rep  <- ds$interp_sum; d_mean <- ds$mean_sum
        if (isTruthy(input$cond_lag)) {
            d_rep  <- filter_by_condition(d_rep,  input$cond_lag)
            d_mean <- filter_by_condition(d_mean, input$cond_lag)
        }
        points_plot_generic(
            d_rep, d_mean,
            y_candidates = c("lag7","Lag","LagTime","lag_time","LT"),
            y_lab = LABEL_LAG_Y
        )
    }, res = 96) %>% bindCache(input$cond_lag)
}

shinyApp(ui, server)
