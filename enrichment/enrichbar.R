enrichbar <- function(res, x, y, filter_col = "NES", colors = NULL, 
                      top = 15, text_size = 3,
                      left = NULL, right = NULL, group_label = c("up", "down"),
                      fill_title = "group", step_len = 2, file = NULL, width = 7,
                      height = 7, scale = 1, offset = 0.01, 
                      legend.position = "top", ...) {

  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(dplyr))

  # add direction info
  up <- res %>%
    as.data.frame() %>%
    dplyr::filter(!!rlang::sym(filter_col) > 0) %>%
    dplyr::arrange(!!rlang::sym(x)) %>%
    dplyr::slice_head(n = top)
  up$group <- group_label[[1]]

  down <- res %>%
    as.data.frame() %>%
    dplyr::filter(!!rlang::sym(filter_col) < 0) %>%
    dplyr::arrange(!!rlang::sym(x)) %>%
    dplyr::slice_head(n = top)
  down$group <- group_label[[2]]

  if (is.null(colors)) {
    colors <- c("tomato", "steelblue")
  }

  # prepare plot data
  df_to_plot <- rbind(up, down) %>%
    dplyr::mutate(
      logp = ifelse(
        group == group_label[[1]],
        -log10(!!rlang::sym(x)),
        log10(!!rlang::sym(x))
      )
    ) %>%
    dplyr::arrange(desc(logp)) %>%
    dplyr::mutate(y = rev(seq_len(nrow(.))))

  # custom x axis
  x_max <- ceiling(max(abs(df_to_plot$logp))/step_len) * step_len
  x_min <- -x_max
  breaks <- c(
    -seq(x_max, 0, -step_len), 
    seq(step_len, x_max, step_len)
  )

  # if (!is.null(left)) {
  #   x_min <- min(x_min, left)
  # }

  # if (!is.null(right)) {
  #   x_max <- max(x_max, right)
  # }

  # determine step
  step <- step_len

  xlab <- switch(as.character(x),
    "pvalue" = "-log10(p value)",
    "p.adjust" = "-log10(p adjust)",
    "qvalue" = "-log10(q value)",
    x
  )

  # plot
  p <- df_to_plot %>%
    ggplot(aes(logp, reorder(!!rlang::sym(y), logp))) +
    geom_bar(aes(fill = group), stat = "identity") +
    geom_text(
      data = df_to_plot[df_to_plot$group == group_label[[1]], ],
      aes(-offset, !!rlang::sym(y), label = !!rlang::sym(y), hjust = 1),
      size = text_size
    ) +
    geom_text(
      data = df_to_plot[df_to_plot$group == group_label[[2]], ],
      aes(offset, !!rlang::sym(y), label = !!rlang::sym(y), hjust = 0),
      size = text_size
    ) +
    scale_fill_manual(
      name = fill_title,
      values = setNames(colors, group_label)
    ) +
    scale_y_discrete(name = NULL, labels = NULL) +
    scale_x_continuous(
      name = xlab,
      limits = c(x_min, x_max),
      breaks = breaks, 
      labels = breaks
    ) +
    theme_bw() +
    theme(
      axis.ticks.y = element_blank(),
      panel.background = element_blank(),
      axis.line.x = element_line(),
      panel.grid = element_blank(), 
      legend.position = legend.position,
      ...
    )

  # save or print
  if (!is.null(file)) {
    ggsave(file = file, width = width, height = height, scale = scale, plot = p)
  } else {
    print(p)
  }
}