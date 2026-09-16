###############################################################################
# 11 - PRISMA flow diagram (Fig. S1) from data/prisma_counts.csv
# Run: Rscript 00_run_all.R only11
###############################################################################
suppressPackageStartupMessages(library(ggplot2))
pc <- read.csv(file.path(PROJECT_DIR, "data", "prisma_counts.csv"), stringsAsFactors = FALSE)
N <- function(k) { v <- pc$n[pc$key == k]; if (!length(v) || is.na(v) || v == "") "n = ..." else paste0("n = ", v) }
L <- function(k) pc$label[pc$key == k]
wrap <- function(s, n) paste(strwrap(s, n), collapse = "\n")
main <- function(y, key, h = 0.75, extra = NULL) list(
  annotate("rect", xmin = 0.4, xmax = 5.4, ymin = y - h / 2, ymax = y + h / 2, fill = "#EAF2F8", colour = "#2F5F8A", linewidth = 0.4),
  annotate("text", x = 2.9, y = y, label = paste0(wrap(L(key), 52), " (", N(key), ")", if (!is.null(extra)) paste0("\n", extra) else ""),
           size = 2.7, lineheight = 0.95))
side <- function(y, key) list(
  annotate("rect", xmin = 6.2, xmax = 9.8, ymin = y - 0.3, ymax = y + 0.3, fill = "white", colour = "#2F5F8A", linewidth = 0.4),
  annotate("text", x = 8.0, y = y, label = paste0(wrap(L(key), 38), " (", N(key), ")"), size = 2.6, lineheight = 0.95),
  annotate("segment", x = 2.9, xend = 6.2, y = y, yend = y, linewidth = 0.35, arrow = arrow(length = unit(1.5, "mm"), type = "closed")))
down <- function(y1, y2) annotate("segment", x = 2.9, xend = 2.9, y = y1, yend = y2, linewidth = 0.4,
                                  arrow = arrow(length = unit(1.6, "mm"), type = "closed"))
ys <- c(identified = 10, non_duplicated = 8.4, available = 6.8, in_scope = 5.2, kept = 3.4, included = 1.4)
g <- ggplot() + xlim(0, 10) + ylim(0.6, 10.6) + theme_void() +
  down(ys[1] - 0.45, ys[2] + 0.38) + down(ys[2] - 0.38, ys[3] + 0.38) + down(ys[3] - 0.38, ys[4] + 0.38) +
  down(ys[4] - 0.38, ys[5] + 0.38) + down(ys[5] - 0.38, ys[6] + 0.55) +
  side(9.2, "duplicates") + side(7.6, "not_accessible") + side(6.0, "out_of_scope") +
  side(4.55, "no_control") + side(3.95, "other_exclusions") + side(2.65, "excluded_analysis") + side(2.05, "excluded_revision") +
  main(ys[["identified"]], "identified", h = 0.9) + main(ys[["non_duplicated"]], "non_duplicated") +
  main(ys[["available"]], "available") + main(ys[["in_scope"]], "in_scope") + main(ys[["kept"]], "kept") +
  main(ys[["included"]], "included", h = 1.0, extra = paste0(L("included_rr"), ": ", N("included_rr"), "\n", L("included_sr"), ": ", N("included_sr")))
dir.create(file.path(OUT_DIR, "figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(OUT_DIR, "figures", "FigS1_PRISMA.png"), g, width = 174, height = 150, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(OUT_DIR, "figures", "FigS1_PRISMA.pdf"), g, width = 174, height = 150, units = "mm")
chk <- function(a, b, c) as.numeric(pc$n[pc$key == a]) - as.numeric(pc$n[pc$key == b]) == as.numeric(pc$n[pc$key == c])
message("PRISMA arithmetic: identified-duplicates=non_duplicated: ", chk("identified", "duplicates", "non_duplicated"))
