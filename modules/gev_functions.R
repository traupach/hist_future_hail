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
suppressMessages(library(purrr))
suppressMessages(library(broom))
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(psych))
suppressMessages(library(tibble))
suppressMessages(library(gridExtra))

letters <- c(
    "a", "b", "c", "d", "e", " f", "g", "h", "i", "j", "k", "l", "m",
    "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"
)

# Labellers for plot elements.
default_labels <- labeller(
    epoch = c(historical = "Historical", ssp245 = "Future"),
    variable = c(
        hailcast_diam_max = "Max. hail size",
        wind_10m = "Max. 10 m wind"
    ),
    parameter = c(
        shape = "Shape",
        location = "Location",
        scale = "Scale"
    ),
    .multi_line = FALSE
)

default_labels_units <- labeller(
    epoch = c(historical = "Historical", ssp245 = "Future"),
    variable = c(
        hailcast_diam_max = "Max.~hail~size~group('[',mm,']')",
        wind_10m = "Max.~wind~group('[',km~h^{-1},']')"
    ),
    domain = c(
        "Adelaide" = "Adelaide",
        "Melbourne" = "Melbourne",
        "Brisbane" = "Brisbane",
        "Kalgoorlie" = "Kalgoorlie",
        "Perth" = "Perth",
        "Sydney/Canberra" = "'Sydney/\nCanberra'"
    ),
    .default = label_parsed
)

default_labels_ml <- labeller(
    epoch = c(historical = "Historical", ssp245 = "Future"),
    variable = c(
        hailcast_diam_max = "Max. hail size",
        wind_10m = "Max. 10 m wind"
    ),
    parameter = c(
        shape = "Shape",
        location = "Location",
        scale = "Scale"
    ),
    domain = c(
        "Sydney/Canberra" = "Sydney/\nCanberra"
    ),
    .multi_line = TRUE
)

default_fontsize <- 18 # Font size for plots.

# Read and concatenate all feather files in `results-dir`.
read_feathers <- function(results_dir, pattern = "*.feather", remove_leaps = TRUE) {
    all_dat <- list()

    files <- list.files(results_dir, pattern = pattern, full.names = TRUE)
    for (file in files) {
        print(file)
        dat <- read_feather(file)
        all_dat <- append(all_dat, list(dat))
    }

    all_dat <- bind_rows(all_dat)

    if (remove_leaps == TRUE) {
        all_dat <- all_dat %>% filter(!(month(time) == 2 & day(time) == 29)) # Remove 29th of February.
    }

    # Replace Sydney + Canberra with Sydney/Canberra.
    all_dat <- all_dat %>% mutate(domain = case_when(domain == "Sydney + Canberra" ~ "Sydney/Canberra", TRUE ~ domain))

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
                        width = 12, height = 6) {
    domain <- low <- high <- epoch <- est <- parameter <- variable <- label <- NULL

    letter_labels <- gev_fits$params %>%
        select(parameter, variable) %>%
        unique() %>%
        group_by(variable, parameter) %>%
        mutate(label = paste("bold(", letters[cur_group_id()], ")", sep = "")) %>%
        ungroup()

    g <- ggplot(gev_fits$params) +
        ggh4x::facet_grid2(variable ~ parameter, scale = "free", labeller = labels, independent = "y") +
        theme_bw(fontsize) +
        geom_errorbar(aes(x = domain, ymin = low, ymax = high, colour = epoch),
            stat = "identity", width = 0.25, linewidth = 1, position = position_dodge(dodge)
        ) +
        geom_label(aes(x = Inf, y = -Inf, label = label),
            data = letter_labels, hjust = "right", vjust = "bottom", size = 5.5,
            parse = TRUE, label.size = 0
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

plot_densities <- function(gev_fits, variable, label, ev_types, gp_thresholds,
                           epochs = c("historical", "ssp245"),
                           x = seq(0, 150), fontsize = default_fontsize, 
                           labels = default_labels_ml,
                           file = NA, width = 12, height = 5) {
    datagrabber <- Data <- NULL # nolint

    domains <- unique(gev_fits$params$domain)
    densities <- tibble()

    for (domain in domains) {
        for (epoch in epochs) {
            d <- gev_fits$gev[[domain]][[variable]][[epoch]]


            exceedances = datagrabber(d)
            if (ev_types[[variable]] == "GEV") {
                threshold = 0
                mod <- devd(
                    x = x,
                    scale = d$results$par[["scale"]],
                    shape = d$results$par[["shape"]],
                    loc = d$results$par[["location"]],
                    type = "GEV"
                )
            } else if (ev_types[[variable]] == "GP") {
                threshold = gp_thresholds[[variable]]
                exceedances = exceedances[exceedances > threshold]
                mod <- devd(
                    x = x,
                    scale = d$results$par[["scale"]],
                    shape = d$results$par[["shape"]],
                    type = "GP"
                )
            }

            res_mod <- tibble(
                epoch = epoch, variable = variable, domain = domain,
                x = x + threshold, density = mod, Data = "EVD model"
            )

            emp <- density(exceedances)
            res_emp <- tibble(
                epoch = epoch, variable = variable, domain = domain,
                x = emp$x, density = emp$y, Data = "WRF simulations"
            )

            densities <- rbind(densities, res_mod, res_emp)
        }
    }

    g <- ggplot(densities, aes(x = x, y = density)) +
        facet_wrap(~domain, ncol = 3) +
        geom_line(aes(colour = epoch, linetype = Data), linewidth = 1) +
        theme_bw(fontsize) +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        labs(x = parse(text = label), y = "Density") +
        scale_colour_discrete(name = "Epoch", breaks = c("historical", "ssp245"), labels = c("Historical", "Future")) +
        scale_x_continuous()
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

    g <- ggplot(gev_fits$quantiles %>% filter(variable == var), aes(x = empirical, y = model)) +
        facet_grid(epoch ~ domain, labeller = labels) +
        geom_point(shape = 1) +
        theme_bw(fontsize) +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        labs(
            x = parse(text = paste("WRF~simulations~quantile~group('[',", unit, ",']')", sep = "")),
            y = parse(text = paste("EVD~model~quantile~group('[',", unit, ",']')", sep = ""))
        ) +
        geom_abline(slope = 1, intercept = 0) +
        coord_fixed(xlim = c(mins, maxs), ylim = c(mins, maxs))
    print(g)

    if (!is.na(file)) {
        ggsave(file, dpi = 300, width = width, height = height)
    }
}

# Plot return level curves for fitted GEVs.
plot_return_levels <- function(gev_fits, var, file = NA, width = 12, height = 6.5,
                               fontsize = default_fontsize, labels = default_labels_units) {
    variable <- epoch <- low <- high <- est <- domain <- label <- NULL

    letter_labels <- gev_fits$return_levels %>%
        select(variable, domain) %>%
        unique() %>%
        group_by(variable, domain) %>%
        mutate(label = paste("bold(", letters[cur_group_id()], ")", sep = "")) %>%
        ungroup()

    g <- gev_fits$return_levels %>%
        ggplot(aes(x = period, y = est)) +
        geom_ribbon(aes(fill = epoch, ymin = low, ymax = high), linewidth = 0.5, alpha = 0.2) +
        geom_line(aes(colour = epoch), linewidth = 1) +
        facet_grid(variable ~ domain, scale = "free_y", labeller = labels) +
        theme_bw(fontsize) +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        labs(y = "Extreme value", x = "Return period [years]") +
        scale_fill_discrete(name = "Epoch", breaks = c("historical", "ssp245"), labels = c("Historical", "Future")) +
        scale_colour_discrete(name = "Epoch", breaks = c("historical", "ssp245"), labels = c("Historical", "Future")) +
        geom_label(aes(x = Inf, y = -Inf, label = label),
            data = letter_labels, hjust = "right", vjust = "bottom",
            label.size = 0, size = 5.5, parse = TRUE
        )
    print(g)

    if (!is.na(file)) {
        ggsave(file, dpi = 300, width = width, height = height)
    }
}

# Plot hail probabilities for comparison of GEV fits.
plot_probs <- function(gev_fits, file = NA, width = 12, height = 6,
                       fontsize = default_fontsize, labels = default_labels_ml) {
    diam <- epoch <- p <- thresh <- windspeed <- NULL
    variable <- domain <- label <- NULL

    probs <- rbind(
        gev_fits$hail_probs %>%
            rename(thresh = diam) %>%
            mutate(variable = "hailcast_diam_max", thresh = factor(paste(thresh, "mm hail"),
                levels = c(
                    "50 mm hail",
                    "100 mm hail"
                )
            )),
        gev_fits$wind_probs %>%
            rename(thresh = windspeed) %>%
            mutate(variable = "wind_10m", thresh = factor(paste(thresh, "km/h wind"),
                levels = c(
                    "70 km/h wind",
                    "90 km/h wind"
                ),
                labels = c(
                    "70 km/h wind",
                    "90 km/h wind"
                )
            ))
    )

    letter_labels <- probs %>%
        select(variable, domain) %>%
        unique() %>%
        group_by(variable, domain) %>%
        mutate(label = paste("bold(", letters[cur_group_id()], ")", sep = "")) %>%
        ungroup()

    g <- ggplot(probs) +
        geom_point(aes(x = factor(thresh), y = p, colour = epoch), shape = 1, size = 3, stroke = 2) +
        facet_grid(variable ~ domain, labeller = labels, scales = "free") +
        theme_bw(fontsize) +
        theme(strip.background = element_blank(), strip.text = element_text(size = fontsize)) +
        scale_colour_discrete(name = "Epoch", breaks = c("historical", "ssp245"), labels = c("Historical", "Future")) +
        labs(x = "Damage threshold", y = "Probability [%]") +
        theme(axis.text.x = element_text(angle = 38, vjust = 1, hjust = 1)) +
        geom_label(aes(x = -Inf, y = -Inf, label = label),
            data = letter_labels, hjust = "left", vjust = "bottom",
            label.size = 0, size = 5.5, parse = TRUE
        )
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
            ks <- tibble(pval = ks, domain = d, scenario = "historical EVD vs ssp245 EVD", variable = v)
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

probabilities_table <- function(gev_fits, out_file,
                                vars = list(hail_probs = "diam", wind_probs = "windspeed"),
                                labels = list(hail_probs = "hail", wind_probs = "wind"),
                                units = list(hail_probs = "mm", wind_probs = "m s$^{-1}$")) {
    thresh <- NULL

    probs <- list()
    for (var in names(vars)) {
        p <- rename(gev_fits[[var]], "thresh" = vars[[var]])
        p <- p %>% mutate(thresh = paste(as.character(thresh), units[[var]], labels[[var]]))
        probs <- c(probs, list(p))
    }

    probs <- reduce(probs, full_join, by = c("domain", "epoch", "thresh", "p"))
    probs$domain <- factor(probs$domain)
    probs$epoch <- factor(probs$epoch, levels = c("historical", "ssp245"), labels = c("Historical", "Future"))
    probs$thresh <- factor(probs$thresh,
        levels = unique(probs$thresh),
    )

    probs <- mutate(probs, p = round(p, 2))

    tab <- tabular(
        Heading("Domain") * domain *
            Heading("Epoch") * epoch ~ Heading("Probability [\\%]") * thresh *
            Heading() * Format(digits = 3) * p * Heading() * identity * Format(digits = 1),
        data = probs,
    )
    print(toLatex(tab, file = out_file))
}

# Do GEV fits per domain, variable and epoch. Collect quantiles for qq plots, do
# ks tests to compare distributions.
#
# Arguments:
#   all_dat: The data to fit to.
#   epochs: Epochs to calculate for.
#   prob_diams: Hail sizes to find probabilities for.
#   p: Quantiles to use in qqplots.
#   return_periods: Return periods to use (years/seasons).
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
                     gp_thresholds,
                     ev_types,
                     span = 20,
                     expected_per_year = 151,
                     epochs = c("historical", "ssp245"),
                     prob_diams = c(50, 100), # mm
                     prob_windspeeds = c(70, 90), # km/h
                     p = seq(1, 99) / 100,
                     return_periods = seq(2, 20, by = 1),
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
                stopifnot(length(dat) / span == expected_per_year)

                type = ev_types[[v]]
                thresh = NA

                if (type == "GP") {
                    # If using GP distribution, subset data to those above threshold; no location parameter.
                    thresh = gp_thresholds[[v]]
                    loc = NA
                    gev <- fevd(dat, type = "GP", threshold = thresh, span = span)
                    dat = dat[dat > thresh]
                    count_per_year = 1
                } else if (type == "GEV") {
                    # If using GEV distribution, no threshold, subset to non-zero values; no span.
                    dat = dat[dat != 0]
                    gev <- fevd(dat, type = "GEV")
                    thresh = NA
                    count_per_year = tibble(v = dat > 0) %>%
                        mutate(year = (row_number() - 1) %/% 151) %>%
                        group_by(year) %>%
                        summarise(n = sum(v)) %>%
                        ungroup() %>%
                        summarize(n = mean(n))
                    count_per_year = count_per_year$n
                }
                gevs[[d]][[v]][[e]] <- gev

                if (type == "GEV") {
                    loc = gev$results$par[["location"]]
                }

                # Calculate model and empirical quantiles for a qqplot.
                qs <- tibble(p = p)
                qs["model"] <- qevd(
                    p = p,
                    threshold = thresh,
                    type = type,
                    scale = gev$results$par[["scale"]],
                    shape = gev$results$par[["shape"]],
                    loc = loc
                )
                qs["empirical"] <- quantile(p = p, dat)
                qs["domain"] <- d
                qs["epoch"] <- e
                qs["variable"] <- v
                quantiles <- append(quantiles, list(qs))

                # Calculate KS tests for model fits.
                ks <- vector()
                for (i in seq(1, ks_iterations)) {
                    ks <- append(ks, ks.test(dat, rextRemes(gev, 1000))$p.value)
                }
                ks <- tibble(pval = ks, domain = d, scenario = paste(e, "EVD vs WRF"), variable = v)
                ks_fits <- append(ks_fits, list(ks))

                # Calculate return levels for given periods.
                ret_level <- return.level(gev,
                    return.period = return_periods * count_per_year,
                    do.ci = TRUE
                )
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
            "historical EVD vs WRF",
            "ssp245 EVD vs WRF",
            "historical EVD vs ssp245 EVD"
        ),
        labels = c(
            "historical EVD vs WRF" = "EVD vs WRF: historical",
            "ssp245 EVD vs WRF" = "EVD vs WRF: future",
            "historical EVD vs ssp245 EVD" = "Historical EVD vs future EVD"
        )
    )

    return(list(
        gevs = gevs, return_levels = return_levels, params = params,
        hail_probs = hail_probs, wind_probs = wind_probs,
        ks_fits = ks_fits, quantiles = quantiles
    ))
}

hail_day_changes <- function(dat, out_file, fontsize = default_fontsize, plot_file = NA, width = 12, height = 3) {
    domain <- epoch <- season <- historical <- ssp245 <- change_from <- change_to <- NULL
    estimate <- estimate1 <- estimate2 <- p.value <- historic <- rel_change <- sig <- NULL # nolint
    conf.low <- conf.high <- mean_ssp245 <- sd_historical <- sd_ssp245 <- NULL # nolint

    hail_days <- dat %>%
        mutate(season = year(time - ddays(60))) %>%
        group_by(domain, epoch, season) %>%
        count() %>%
        group_by(epoch) %>%
        mutate(year = season - min(season) + 1) %>%
        ungroup()

    stats <- hail_days %>%
        group_by(domain, epoch) %>%
        reframe(mean = mean(n), sd = round(sd(n), 1)) %>%
        pivot_wider(names_from = "epoch", values_from = c("mean", "sd"))

    t_test <- hail_days %>%
        select(!season) %>%
        pivot_wider(names_from = "epoch", values_from = "n") %>%
        group_by(domain) %>%
        summarise(tidy(t.test(x = ssp245, y = historical)))

    t_test_res <- t_test %>%
        reframe(domain,
            historic = estimate2,
            rel_change = estimate / abs(estimate2) * 100,
            abs_change = estimate,
            change_from = conf.low / abs(estimate2) * 100, # nolint
            change_to = conf.high / abs(estimate2) * 100, # nolint
            sig_010 = p.value < 0.1, # nolint
            sig_005 = p.value < 0.05, # nolint
            sig_001 = p.value < 0.01 # nolint
        ) %>%
        arrange(desc(historic))

    t_test_res <- full_join(t_test_res, stats, by = "domain")
    stopifnot(all(t_test_res$historic == t_test_res$mean_historical))

    t_test_disp <- t_test_res %>%
        mutate(rel_change = paste(as.character(round(rel_change, 0)), "\\%", sep = "")) %>%
        mutate(change_range = paste("(", as.character(round(change_from, 0)),
            " to ", as.character(round(change_to, 0)), "\\%)",
            sep = ""
        )) %>%
        mutate(sig = case_when(sig_010 == TRUE ~ "\\ast{}", TRUE ~ "")) %>%
        mutate(sig = case_when(sig_005 == TRUE ~ paste(sig, "\\!\\ast{}", sep = ""), TRUE ~ sig)) %>%
        mutate(sig = case_when(sig_001 == TRUE ~ paste(sig, "\\!\\!\\ast{}", sep = ""), TRUE ~ sig)) %>%
        mutate(sig = paste("$", sig, "$", sep = "")) %>%
        select(!starts_with("sig_"))

    # Plot boxplot of changes by domain.
    p <- ggplot(hail_days, aes(x = domain, y = n)) +
        geom_boxplot(aes(fill = epoch), width = 0.75, alpha = 0.75) +
        theme_bw(fontsize) +
        scale_fill_discrete(name = "Epoch", breaks = c("historical", "ssp245"), labels = c("Historical", "Future")) +
        labs(x = "Domain", y = "Seasonal hail days")
    print(p)

    if (!is.na(plot_file)) {
        ggsave(plot_file, dpi = 300, width = width, height = height)
    }

    print(t_test)

    seasonal_hail_days <- t_test %>%
        select(domain, estimate1, estimate2) %>%
        rename(historical = estimate2, ssp245 = estimate1) %>%
        pivot_longer(cols = c("ssp245", "historical"), names_to = "epoch", values_to = "frequency")

    return(list(seasonal_hail_days = seasonal_hail_days, t_test_disp = t_test_disp))
}

ingredients_changes <- function(ings, plot = FALSE) {
    domain <- historical <- ssp245 <- estimate2 <- estimate <- v <- variable <- NULL
    conf.low <- conf.high <- p.value <- NULL # nolint
    epoch <- rel_change <- change_from <- change_to <- sig <- NULL

    vars <- c(
        "hailcast_diam_max", "wind_10m", "mixed_100_cape", "mixed_100_cin", "mixed_100_lifted_index",
        "lapse_rate_700_500", "temp_500", "freezing_level", "melting_level", "shear_magnitude"
    )

    t_results <- tibble()
    for (var in vars) {
        d <- ings %>% select("domain", "time", "epoch", v = all_of(var))

        t_res <- d %>%
            pivot_wider(names_from = epoch, values_from = v) %>%
            select(domain, historical, ssp245) %>%
            group_by(domain) %>%
            summarise(tidy(t.test(x = ssp245, y = historical))) %>%
            mutate(variable = var)

        t_results <- rbind(t_results, t_res)
    }

    t_test_res <- t_results %>%
        reframe(domain, variable,
            historic = estimate2,
            rel_change = estimate / abs(estimate2) * 100,
            abs_change = estimate,
            change_from = conf.low / abs(estimate2) * 100,
            change_to = conf.high / abs(estimate2) * 100,
            sig_010 = p.value < 0.1,
            sig_005 = p.value < 0.05,
            sig_001 = p.value < 0.01
        )

    t_test_disp <- t_test_res %>%
        mutate(rel_change = paste(as.character(round(rel_change, 0)), "\\%", sep = "")) %>%
        mutate(change_range = paste("(", as.character(round(change_from, 0)),
            " to ", as.character(round(change_to, 0)), "\\%)",
            sep = ""
        )) %>%
        mutate(sig = case_when(sig_010 == TRUE ~ "\\ast{}", TRUE ~ "")) %>%
        mutate(sig = case_when(sig_005 == TRUE ~ paste(sig, "\\!\\ast{}", sep = ""), TRUE ~ sig)) %>%
        mutate(sig = case_when(sig_001 == TRUE ~ paste(sig, "\\!\\!\\ast{}", sep = ""), TRUE ~ sig)) %>%
        mutate(sig = paste("$", sig, "$", sep = "")) %>%
        select(!starts_with("sig_"))

    t_test_disp <- t_test_disp %>%
        mutate(variable = replace(variable, variable == "hailcast_diam_max", "Hail size")) %>%
        mutate(variable = replace(variable, variable == "mixed_100_cape", "CAPE")) %>%
        mutate(variable = replace(variable, variable == "mixed_100_cin", "CIN")) %>%
        mutate(variable = replace(variable, variable == "mixed_100_lifted_index", "LI")) %>%
        mutate(variable = replace(variable, variable == "wind_10m", "Wind")) %>%
        mutate(variable = replace(variable, variable == "lapse_rate_700_500", "LR")) %>%
        mutate(variable = replace(variable, variable == "temp_500", "T500")) %>%
        mutate(variable = replace(variable, variable == "freezing_level", "FLH")) %>%
        mutate(variable = replace(variable, variable == "melting_level", "MLH")) %>%
        mutate(variable = replace(variable, variable == "shear_magnitude", "S06"))

    # Error-bar plot.
    if (plot) {
        g <- t_test_res %>% ggplot(aes(x = variable, y = rel_change)) +
            geom_point(aes(color = domain), position = position_dodge(0.5)) +
            geom_errorbar(aes(ymax = change_to, ymin = change_from, color = domain),
                width = 0.5,
                linewidth = 1, position = position_dodge(0.5)
            ) +
            theme_bw()
        print(g)
    }

    return(t_test_disp)
}

domain_correlation_plot <- function(means, fontsize = default_fontsize, plot_file = NA, width = 12, height = 10) {
    v <- domain <- epoch <- season <- rowname <- name <- value <- from <- to <- sig <- label <- NULL

    vars <- c(
        "hailcast_diam_max", "wind_10m", "shear_magnitude", "mixed_100_cape", "mixed_100_lifted_index",
        "mixed_100_cin", "lapse_rate_700_500", "temp_500", "freezing_level", "melting_level"
    )

    domains <- (means %>% select(domain) %>% unique())[["domain"]]
    mask <- lower.tri(matrix(nrow = length(domains), ncol = length(domains))) # nolint

    seasonal <- means %>%
        mutate(season = year(time - ddays(60))) %>%
        group_by(domain, epoch, season) %>%
        summarise_at(vars, mean) %>%
        arrange(domain, epoch, season) %>%
        ungroup()

    counts <- means %>%
        mutate(season = year(time - ddays(60))) %>%
        group_by(domain, epoch, season) %>%
        count() %>%
        rename(hail_days = n) %>%
        ungroup()

    seasonal <- full_join(seasonal, counts, by = c("domain", "epoch", "season"))
    corrs <- tibble()
    ps <- tibble()

    for (var in c(vars, "hail_days")) {
        d <- seasonal %>% select(domain, epoch, season, v = all_of(var))

        r <- d %>%
            pivot_wider(names_from = domain, values_from = v) %>%
            arrange(epoch, season) %>%
            select(!season) %>%
            group_by(epoch) %>%
            group_map(~ rownames_to_column(cbind(.y, var, replace(corr.test(.x, ci = FALSE)$r, !mask, NA))))

        p <- d %>%
            pivot_wider(names_from = domain, values_from = v) %>%
            arrange(epoch, season) %>%
            select(!season) %>%
            group_by(epoch) %>%
            group_map(~ rownames_to_column(cbind(.y, var, replace(corr.test(.x, ci = FALSE)$p, !mask, NA))))

        corrs <- rbind(corrs, bind_rows(r))
        ps <- rbind(ps, bind_rows(p))
    }

    corrs <- corrs %>%
        pivot_longer(cols = all_of(domains)) %>%
        rename(from = rowname, to = name, r = value) %>%
        drop_na()
    ps <- ps %>%
        pivot_longer(cols = all_of(domains)) %>%
        rename(from = rowname, to = name, p = value) %>%
        drop_na()

    corrs <- full_join(corrs, ps, by = c("from", "epoch", "var", "to"))

    corrs <- corrs %>% mutate(sig = case_when(p < 0.1 ~ "*", TRUE ~ ""))
    corrs <- corrs %>% mutate(sig = case_when(p < 0.05 ~ "**", TRUE ~ sig))
    corrs <- corrs %>% mutate(sig = case_when(p < 0.01 ~ "***", TRUE ~ sig))

    ing_labels <- labeller(
        epoch = c(historical = "Hist.", ssp245 = "Fut."),
        var = c(
            freezing_level = "FLH",
            hailcast_diam_max = "Hail size",
            lapse_rate_700_500 = "LR",
            melting_level = "MLH",
            wind_10m = "Wind",
            mixed_100_cape = "CAPE",
            mixed_100_cin = "CIN",
            mixed_100_lifted_index = "LI",
            shear_magnitude = "S06",
            temp_500 = "T500",
            hail_days = "Days"
        ),
        .multi_line = FALSE
    )

    grobs <- list()
    i <- 1
    for (p in list(c(1:6), c(7:11))) {
        vs <- c("hail_days", vars)[p]
        dat <- corrs %>%
            filter(var %in% vs) %>%
            mutate(var = factor(var, levels = vs))

        letter_labels <- dat %>%
            select(epoch, var) %>%
            unique() %>%
            group_by(epoch, var) %>%
            mutate(label = paste("bold(", letters[cur_group_id() + (min(p) * 2) - 2], ")", sep = "")) %>%
            ungroup()

        grobs[[i]] <- ggplot(dat, aes(x = from, y = to)) +
            geom_tile(aes(fill = r)) +
            facet_grid(epoch ~ var, labeller = ing_labels) +
            theme_bw(fontsize) +
            theme(
                strip.background = element_blank(),
                strip.text = element_text(size = fontsize)
            ) +
            scale_fill_gradientn(
                colours = c("darkblue", "white", "darkred"),
                na.value = "white", limits = c(-1, 1),
                name = "Pearson's r"
            ) +
            geom_text(aes(label = sig), vjust = 1, hjust = 0.5) +
            coord_equal() +
            labs(x = "", y = "") +
            theme(panel.spacing = unit(0.15, "lines")) +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
            geom_label(aes(x = -Inf, y = Inf, label = label),
                data = letter_labels,
                hjust = "left", vjust = "top",
                label.size = 0, size = 5.5, parse = TRUE
            )

        if (i == 1) {
            grobs[[i]] <- grobs[[i]] + guides(fill = "none")
            grobs[[i]] <- grobs[[i]] + theme(axis.text.x = element_blank())
        }

        i <- i + 1
    }

    if (!is.na(plot_file)) {
        ggsave(plot_file, arrangeGrob(grobs = grobs, nrow = 2, heights = c(1, 1.3)),
            dpi = 300, width = width, height = height
        )
    }
    grid.arrange(grobs = grobs, nrow = 2, heights = c(1, 1.285))
}

periods_for_thresholds <- function(var, thresh) {
    res <- tibble()
    for (t in thresh) {
        res <- rbind(res, gev_fits$return_levels %>%
            filter(est > t, variable == var) %>%
            group_by(domain, epoch) %>%
            filter(est == min(est)) %>%
            mutate(threshold = t, period = round(period, 0)) %>%
            select(threshold, variable, domain, epoch, period))
    }

    return(res)
}
