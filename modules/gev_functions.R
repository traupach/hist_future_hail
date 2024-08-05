###
# Helper functions for GEV analysis of hail sizes and wind speed in R.
#
# Tim Raupach <t.raupach@unsw.edu.au>
###

# Quietly load libraries.
suppressMessages(library(arrow))
suppressMessages(library(extRemes))
suppressMessages(library(lubridate))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(stringr))
suppressMessages(library(tables))

# Labellers for plot elements.
default_labels <- labeller(
    epoch = c(historical = "Historical", ssp245 = "Future"),
    variable = c(
        hailcast_diam_max = "Maximum hail size",
        wind_10m = "Max. 10 m wind with hail"
    ),
    parameter = c(
        shape = "Shape",
        location = "Location",
        scale = "Scale"
    ),
    domain = c(
        "Sydney + Canberra" = "Sydney/Canberra"
    ),
    .multi_line = FALSE
)

default_labels_ml <- labeller(
    epoch = c(historical = "Historical", ssp245 = "Future"),
    variable = c(
        hailcast_diam_max = "Maximum hail size",
        wind_10m = "Max. 10 m wind with hail"
    ),
    parameter = c(
        shape = "Shape",
        location = "Location",
        scale = "Scale"
    ),
    domain = c(
        "Sydney + Canberra" = "Sydney/\nCanberra"
    ),
    .multi_line = TRUE
)

default_fontsize <- 18 # Font size for plots.

# Read and concatenate all feather files in `results-dir`.
read_feathers <- function(results_dir) {
    all_dat <- list()

    files <- list.files(results_dir, pattern = "*.feather", full.names = TRUE)
    for (file in files) {
        print(file)
        dat <- read_feather(file)
        all_dat <- append(all_dat, list(dat))
    }

    all_dat <- bind_rows(all_dat)
    return(all_dat)
}

# Plot a timeseries of data from `dat`.
plot_ts <- function(dat, var, ylabel, xlabel = "Year", file = NA,
                    width = 12, height = 10, fontsize = default_fontsize, labels = default_labels_ml) {
    g <- ggplot(dat, aes(x = time, y = .data[[var]])) +
        geom_point(shape = 1) +
        facet_grid(domain ~ epoch, scales = "free_x", labeller = labels) +
        theme_bw(fontsize) +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        labs(x = xlabel, y = ylabel) +
        theme(
            panel.spacing.x = unit(2, "lines"),
            panel.spacing.y = unit(1, "lines")
        )
    print(g)

    if (!is.na(file)) {
        ggsave(file, dpi = 300, width = width, height = height)
    }
}

# Plot parameters of GEV fits returned by fit_gevs.
plot_params <- function(gev_fits, fontsize = default_fontsize, dodge = 0.3, labels = default_labels_ml, file = NULL,
                        width = 12, height = 8) {
    domain <- low <- high <- epoch <- est <- NULL

    g <- ggplot(gev_fits$params) +
        ggh4x::facet_grid2(variable ~ parameter, scale = "free", labeller = labels, independent = "y") +
        theme_bw(fontsize) +
        geom_errorbar(aes(x = domain, ymin = low, ymax = high, colour = epoch),
            stat = "identity", width = 0.25, linewidth = 1, position = position_dodge(dodge)
        ) +
        geom_point(aes(x = domain, y = est, colour = epoch), position = position_dodge(dodge), size = 3) +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        labs(x = "Domain", y = "Parameter value") +
        scale_colour_discrete(name = "Epoch", breaks = c("historical", "ssp245"), labels = c("Historical", "Future")) +
        theme(axis.text.x = element_text(angle = 25, vjust = 1, hjust = 1))
    print(g)

    if (!is.na(file)) {
        ggsave(file, dpi = 300, width = width, height = height)
    }
}

# Plot KS test p value distributions for the various fits in a fitted set of GEVs.
plot_ks_fits <- function(gev_fits, file = NA, fontsize = default_fontsize, labels = default_labels) {
    domain <- scenario <- pval <- NULL

    g <- ggplot(gev_fits$ks_fits) +
        geom_boxplot(aes(x = domain, fill = scenario, y = pval), colour = "black", alpha = 0.75, width = 0.5) +
        theme_bw(fontsize) +
        scale_fill_discrete(name = "Scenario") +
        geom_hline(yintercept = 0.05, colour = "red") +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        labs(y = "p value", x = "Domain") +
        facet_wrap(~variable, nrow = 2, labeller = labels)
    print(g)

    if (!is.na(file)) {
        ggsave(file, dpi = 300, width = 12, height = 6)
    }
}

# Plot qq plots of GEV fitted functions.s
plot_quantiles <- function(gev_fits, var, unit, labels = default_labels_ml, fontsize = default_fontsize,
                           width = 12, height = 5, file = NA) {
    variable <- model <- empirical <- NULL

    vals <- gev_fits$quantiles %>% filter(variable == var)
    mins <- min(min(vals$model), min(vals$empirical))
    maxs <- max(max(vals$model), max(vals$empirical))

    g <- ggplot(gev_fits$quantiles %>% filter(variable == var), aes(x = model, y = empirical)) +
        # facet_wrap(epoch ~ domain, nrow = 2, labeller = labels) +
        facet_grid(epoch ~ domain, labeller = labels) +
        geom_point(shape = 1) +
        theme_bw(fontsize) +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        labs(
            x = parse(text = paste("Model~quantile~group('[',", unit, ",']')", sep = "")),
            y = parse(text = paste("Empirical~quantile~group('[',", unit, ",']')", sep = ""))
        ) +
        geom_abline(slope = 1, intercept = 0) +
        coord_fixed(xlim = c(mins, maxs), ylim = c(mins, maxs))
    print(g)

    if (!is.na(file)) {
        ggsave(file, dpi = 300, width = width, height = height)
    }
}

# Plot return level curves for fitted GEVs.
plot_return_levels <- function(gev_fits, var, varname, file = NA, width = 12, height = 3,
                               fontsize = default_fontsize, labels = default_labels_ml) {
    variable <- epoch <- low <- high <- est <- NULL

    g <- gev_fits$return_levels %>%
        filter(variable == var) %>%
        ggplot(aes(x = period, y = est)) +
        geom_ribbon(aes(fill = epoch, ymin = low, ymax = high), linewidth = 0.5, alpha = 0.2) +
        geom_line(aes(colour = epoch), linewidth = 1) +
        facet_wrap(~domain, nrow = 1, labeller = labels) +
        theme_bw(fontsize) +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        labs(y = parse(text = varname), x = "Return period [hail days]") +
        scale_fill_discrete(name = "Epoch", breaks = c("historical", "ssp245"), labels = c("Historical", "Future")) +
        scale_colour_discrete(name = "Epoch", breaks = c("historical", "ssp245"), labels = c("Historical", "Future")) +
        scale_x_log10()
    print(g)

    if (!is.na(file)) {
        ggsave(file, dpi = 300, width = width, height = height)
    }
}

# Plot hail probabilities for comparison of GEV fits.
plot_hail_probs <- function(gev_fits, file = NA, width = 12, height = 3,
                            fontsize = default_fontsize, labels = default_labels_ml) {
    diam <- epoch <- p <- NULL

    g <- ggplot(gev_fits$hail_probs) +
        geom_point(aes(x = diam, y = p, colour = epoch), shape = 1, size = 3, stroke = 2) +
        facet_wrap(~domain, nrow = 1, labeller = labels) +
        theme_bw(fontsize) +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        scale_colour_discrete(name = "Epoch", breaks = c("historical", "ssp245"), labels = c("Historical", "Future")) +
        labs(x = "Hail diameter [mm]", y = "Probability [%]")
    print(g)

    if (!is.na(file)) {
        ggsave(file, dpi = 300, width = width, height = height)
    }
}

# Calculate ks tests to show differences between historical and ssp245 model for each domain.
ks_tests <- function(gevs, domains, variables, ks_iterations) {
    ks_change <- list()
    for (d in domains) {
        for (v in variables) {
            ks <- vector()
            for (i in seq(1, ks_iterations)) {
                ks <- append(ks, ks.test(
                    rextRemes(gevs[[d]][[v]][["historical"]], 1000),
                    rextRemes(gevs[[d]][[v]][["ssp245"]], 1000)
                )$p.value)
            }
            ks <- tibble(pval = ks, domain = d, scenario = "historical model vs ssp245 model", variable = v)
            ks_change <- append(ks_change, list(ks))
        }
    }

    return(ks_change)
}

# Collect together model parameters and confidence intervals.
collect_params <- function(gevs, domains, variables, epochs) {
    params <- list()
    for (d in domains) {
        for (v in variables) {
            for (e in epochs) {
                cis <- ci.fevd(gevs[[d]][[v]][[e]], type = "parameter")
                vars <- rownames(cis)
                for (i in seq(1, length(vars))) {
                    params <- append(params, list(tibble(
                        parameter = vars[i], low = cis[i, 1],
                        est = cis[i, 2], high = cis[i, 3],
                        domain = d, epoch = e, variable = v
                    )))
                }
            }
        }
    }
    params <- bind_rows(params)
    return(params)
}

probabilities_table <- function(gev_fits, var, units) {
    diam <- windspeed <- NULL
    probs <- gev_fits[[var]]
    probs$domain <- factor(probs$domain)
    probs$epoch <- factor(probs$epoch, levels = c("historical", "ssp245"), labels = c("Historical", "Future"))

    if (var == "hail_probs") {
        probs <- rename(probs, by_var = diam)
    } else if (var == "wind_probs") {
        probs <- rename(probs, by_var = windspeed)
    } else {
        stop("Unknown value for var.")
    }

    probs$by_var <- factor(probs$by_var,
        levels = unique(probs$by_var),
        labels = paste(unique(probs$by_var), units)
    )

    tab <- tabular(
        Heading("Domain") * domain *
            Heading("Epoch") * epoch ~ Heading("Probability [\\%]") * by_var *
            Heading() * p * Heading() * identity * Format(digits = 1),
        data = probs,
    )
    print(tab)
}

# Do GEV fits per domain, variable and epoch. Collect quantiles for qq plots, do
# ks tests to compare distributions.
#
# Arguments:
#   all_dat: The data to fit to.
#   epochs: Epochs to calculate for.
#   prob_diams: Hail sizes to find probabilities for.
#   p: Quantiles to use in qqplots.
#   return_periods: Return periods to use (years of hail days).
#   ks.iterations: Iterations to use for ks tests.
#
# Returns: gevs: GEV objects,
#   return_levels: return levels summary,
#   params: GEV parameters,
#   hail_probs: hail probabilies for given sizes,
#   wind_probs: wind probabilies for given speeds,
#   ks_fits: KS test results,
#   quantiles: Quantiles of empirical vs modeled amounts.
fit_gevs <- function(all_dat,
                     epochs = c("historical", "ssp245"),
                     prob_diams = c(20, 50, 100),
                     prob_windspeeds = c(27.78), # 100 km/h
                     p = seq(1, 99) / 100,
                     return_periods = seq(2, 100),
                     ks_iterations = 100,
                     variables = c("hailcast_diam_max", "wind_10m")) {
    # Set tidyverse components to NULL to avoid lintr complaints.
    epoch <- domain <- NULL

    # List of domains from data.
    domains <- distinct(all_dat, domain)$domain

    # gevs contains [domain][variable][epoch]
    gevs <- list(list(list()))
    quantiles <- list()
    return_levels <- list()
    hail_probs <- list()
    wind_probs <- list()
    ks_fits <- list()

    for (d in domains) {
        for (v in variables) {
            for (e in epochs) {
                print(paste("Fitting for", v, "in", d, "for", e))

                dat <- filter(all_dat, domain == d, epoch == e)[[v]]
                gev <- fevd(dat, type = "GEV", time.units = "days")
                gevs[[d]][[v]][[e]] <- gev

                # Calculate model and empirical quantiles for a qqplot.
                qs <- tibble(p = p)
                qs["model"] <- qevd(
                    p = p,
                    loc = gev$results$par[["location"]],
                    scale = gev$results$par[["scale"]],
                    shape = gev$results$par[["shape"]]
                )
                qs["empirical"] <- quantile(p = p, dat)
                qs["domain"] <- d
                qs["epoch"] <- e
                qs["variable"] <- v
                quantiles <- append(quantiles, list(qs))

                # Calculate return levels for given periods.
                ret_level <- return.level(gev, return.period = return_periods, do.ci = TRUE)
                return_levels <- append(return_levels, list(tibble(
                    low = ret_level[, 1], est = ret_level[, 2],
                    high = ret_level[, 3], domain = d, epoch = e,
                    variable = v,
                    period = return_periods
                )))

                # Calculate return periods for given levels.
                if (v == "hailcast_diam_max") {
                    hail_probs <- append(hail_probs, list(tibble(
                        domain = d, epoch = e, diam = prob_diams,
                        p = pextRemes(gev, q = prob_diams, lower.tail = FALSE) * 100
                    )))
                }
                if (v == "wind_10m") {
                    wind_probs <- append(wind_probs, list(tibble(
                        domain = d, epoch = e, windspeed = prob_windspeeds,
                        p = pextRemes(gev, q = prob_windspeeds, lower.tail = FALSE) * 100
                    )))
                }

                # Calculate KS tests for model fits.
                ks <- vector()
                for (i in seq(1, ks_iterations)) {
                    ks <- append(ks, ks.test(dat, rextRemes(gev, 1000))$p.value)
                }
                ks <- tibble(pval = ks, domain = d, scenario = paste(e, "model vs empirical"), variable = v)
                ks_fits <- append(ks_fits, list(ks))
            }
        }
    }

    quantiles <- bind_rows(quantiles)
    return_levels <- bind_rows(return_levels)
    ks_fits <- bind_rows(ks_fits)
    hail_probs <- bind_rows(hail_probs)
    wind_probs <- bind_rows(wind_probs)

    ks_fits <- bind_rows(ks_fits, ks_tests( # nolint
        gevs = gevs, domains = domains,
        variables = variables, ks_iterations = ks_iterations
    ))
    params <- collect_params( # nolint
        gevs = gevs, domains = domains,
        variables = variables, epochs = epochs
    )

    # Set plot orders.
    hail_probs$diam <- factor(hail_probs$diam, levels = prob_diams)
    ks_fits$scenario <- factor(ks_fits$scenario,
        levels = c(
            "historical model vs empirical",
            "ssp245 model vs empirical",
            "historical model vs ssp245 model"
        ),
        labels = c("Model vs empirical: historical",
            "ssp245 model vs empirical" = "Model vs empirical: future",
            "historical model vs ssp245 model" = "Historical model vs future model"
        )
    )

    return(list(
        gevs = gevs, return_levels = return_levels, params = params,
        hail_probs = hail_probs, wind_probs = wind_probs,
        ks_fits = ks_fits, quantiles = quantiles
    ))
}
