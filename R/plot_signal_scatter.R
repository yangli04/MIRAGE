#' Scatter plot of treatment vs control mutation rate from MIRAGE results
#'
#' \code{plot_signal_scatter} plots per-site control vs treatment mutation
#' (or conversion) rates from a MIRAGE result, optionally colored by
#' significance call or by inferred site-level signal \code{beta_est}. It is a
#' one-line diagnostic for the output of
#' \code{\link{estimate_inference_with_empirical}}.
#'
#' @param res The list returned by \code{\link{estimate_inference_with_empirical}}.
#' @param color_by How to color points. \code{"significance"} marks sites with
#'   \code{res[[fdr_col]] < fdr_cutoff} as "Identified". \code{"beta_est"}
#'   shades by the estimated site-level signal. \code{"none"} uses one color.
#'   Default \code{"significance"}.
#' @param fdr_col,fdr_cutoff Column and cutoff used when
#'   \code{color_by = "significance"}. Default \code{"lrt_fdr"} and 0.05.
#' @param sites Which sites to plot: \code{"both"} (default), \code{"homo"},
#'   or \code{"heter"}.
#' @param xlab,ylab Axis labels.
#' @param point_size,point_alpha Aesthetic controls for \code{geom_point}.
#' @param rasterize If \code{TRUE} and \pkg{ggrastr} is installed, rasterize
#'   the point layer at \code{rasterize_dpi}. Useful for millions of points.
#' @param rasterize_dpi DPI for the rasterized point layer.
#'
#' @return A \code{ggplot} with x and y on a shared scale and aspect ratio 1.
#'
#' @examples
#' count_table <- read.delim(system.file("extdata",
#'   "pos.read.count.example.parclip.txt", package = "MIRAGE"))
#' res <- with(count_table, estimate_inference_with_empirical(
#'   cbind(pos, motif, type), treated_fixed_count, control_fixed_count,
#'   treated_depth, control_depth, bg.method = "lrt", bg.target = "both",
#'   seed = 123))
#' plot_signal_scatter(res)                               # by significance
#' plot_signal_scatter(res, color_by = "beta_est")        # by inferred signal
#'
#' @export

plot_signal_scatter <- function(res,
                                color_by = c("significance", "beta_est", "none"),
                                fdr_col = "lrt_fdr",
                                fdr_cutoff = 0.05,
                                sites = c("both", "homo", "heter"),
                                xlab = "Control mutation rate (%)",
                                ylab = "Treatment mutation rate (%)",
                                point_size = 0.6,
                                point_alpha = 0.4,
                                rasterize = FALSE,
                                rasterize_dpi = 300) {

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_signal_scatter() requires the ggplot2 package.", call. = FALSE)

  color_by <- match.arg(color_by)
  sites <- match.arg(sites)

  needed <- c("treatment_X_rate", "control_X_rate",
              if (color_by == "significance") fdr_col,
              if (color_by == "beta_est") "beta_est")

  pick <- function(df) df[, needed, drop = FALSE]
  parts <- list()
  if (sites != "heter") parts$homo  <- pick(res$homosites)
  if (sites != "homo")  parts$heter <- pick(res$hetersites)
  df <- do.call(rbind, parts)
  row.names(df) <- NULL
  df$x <- df$control_X_rate * 100
  df$y <- df$treatment_X_rate * 100
  lim <- c(0, max(c(df$x, df$y, 1), na.rm = TRUE))

  point_layer <- ggplot2::geom_point(alpha = point_alpha, size = point_size,
                                     stroke = 0)
  if (isTRUE(rasterize)) {
    if (requireNamespace("ggrastr", quietly = TRUE))
      point_layer <- ggrastr::rasterise(point_layer, dpi = rasterize_dpi)
    else
      warning("ggrastr is not installed; falling back to vector points.",
              call. = FALSE)
  }

  if (color_by == "significance") {
    df$status <- factor(
      ifelse(!is.na(df[[fdr_col]]) & df[[fdr_col]] < fdr_cutoff,
             "Identified", "Background"),
      levels = c("Background", "Identified"))
    df <- df[order(df$status), , drop = FALSE]
    sub <- sprintf("Identified: %s    Background: %s",
                   format(sum(df$status == "Identified"), big.mark = ","),
                   format(sum(df$status == "Background"), big.mark = ","))
    mapping <- ggplot2::aes(x = .data$x, y = .data$y, color = .data$status)
    color_scale <- ggplot2::scale_color_manual(
      name = NULL,
      values = c(Background = "#CCCCCC", Identified = "#D55E00"))
    color_guide <- ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 3, alpha = 0.9)))
  } else if (color_by == "beta_est") {
    df <- df[order(df$beta_est, na.last = FALSE), , drop = FALSE]
    sub <- sprintf("Sites: %s", format(nrow(df), big.mark = ","))
    mapping <- ggplot2::aes(x = .data$x, y = .data$y, color = .data$beta_est)
    color_scale <- ggplot2::scale_color_viridis_c(
      name = expression(hat(beta)), option = "magma", direction = -1)
    color_guide <- NULL
  } else {
    sub <- sprintf("Sites: %s", format(nrow(df), big.mark = ","))
    mapping <- ggplot2::aes(x = .data$x, y = .data$y)
    df$.col <- "Sites"
    mapping <- ggplot2::aes(x = .data$x, y = .data$y, color = .data$.col)
    color_scale <- ggplot2::scale_color_manual(values = c(Sites = "#666666"),
                                               guide = "none")
    color_guide <- NULL
  }

  ggplot2::ggplot(df, mapping) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                         color = "gray50", linewidth = 0.4) +
    point_layer +
    color_scale +
    color_guide +
    ggplot2::scale_x_continuous(limits = lim) +
    ggplot2::scale_y_continuous(limits = lim) +
    ggplot2::labs(x = xlab, y = ylab, subtitle = sub) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(color = "gray90",
                                               linetype = "dashed",
                                               linewidth = 0.2),
      text = ggplot2::element_text(family = "Helvetica"),
      axis.title = ggplot2::element_text(size = 14),
      plot.subtitle = ggplot2::element_text(size = 10),
      legend.position = "bottom",
      aspect.ratio = 1)
}
