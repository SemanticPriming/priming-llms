# Run from the repository root: Rscript presentation/plot-benchmark.R
winner <- read.csv("outputs/discovery_grid.csv")
winner <- winner[which.max(winner$r_dS), ]
stopifnot(winner$model == "allenai/OLMo-1B-hf", winner$format == "after_word")

human <- read.csv("data/priming_data/en_priming.csv")
model <- read.csv("outputs/surprisal/allenai_OLMo-1B-hf__after_word.csv")
split <- read.csv("outputs/split.csv")
items <- merge(model, human, by.x = "target", by.y = "target_word_unique")
test <- subset(merge(items, split, by = "target"), half == "confirmation")
reported <- subset(read.csv("outputs/winner_summary.csv"), half == "confirmation")
stopifnot(nrow(test) == reported$n,
          isTRUE(all.equal(cor(test$dS, test$avgZ_prime, method = "spearman"),
                           reported$spearman)))

ink <- "#17344A"
teal <- "#087F78"
svg("presentation/held-out-items.svg", width = 8.2, height = 5.4, pointsize = 15)
par(mar = c(4.4, 5.6, 0.6, 0.9), fg = ink, col.axis = ink, col.lab = ink,
    las = 1, bty = "l", mgp = c(2.65, 0.6, 0), family = "sans")
plot(test$dS, test$avgZ_prime, pch = 19, cex = 0.8, ylim = c(min(test$avgZ_prime), .7),
     col = adjustcolor(teal, alpha.f = 0.42),
     xlab = "Model priming (surprisal difference, nats)",
     ylab = "Human priming (standardized RT difference)")
labelled <- list(stupid = 4, fingerprint = 2, sink = 3, that = 3, scandal = 4,
                 july = 2, spaghetti = 1, fly = 4)
stopifnot(all(names(labelled) %in% test$target))
for (w in names(labelled)) {
  r <- test[test$target == w, ]
  points(r$dS, r$avgZ_prime, pch = 21, cex = 1.1, col = ink, lwd = 1.4)
  text(r$dS, r$avgZ_prime, paste0(r$cue_related, " \u2013 ", w), pos = labelled[[w]],
       cex = .72, col = ink, offset = .5)
}
dev.off()

formats <- read.csv("outputs/confirmation_by_format.csv", check.names = FALSE)
selected <- match(c("bare", "bare+k2", "bare+k5", "bare+k10", "after_word"), formats$format)
stopifnot(!anyNA(selected))
averages <- vapply(list(selected[1], selected[2:4], selected[5]),
                   function(rows) mean(as.matrix(formats[rows, -1])), numeric(1))
print(averages)
svg("presentation/prompt-comparison.svg", width = 9.6, height = 4.2, pointsize = 16)
par(mar = c(3.8, 9.8, 0.6, 1.3), fg = ink, col.axis = ink, col.lab = ink,
    las = 1, mgp = c(2.45, 0.55, 0), family = "sans")
positions <- barplot(rev(averages), horiz = TRUE, xlim = c(0, .265),
                     names.arg = rev(c("Bare pair", "Bare + examples", "After the word…")),
                     col = c(teal, "#7695A4", "#7695A4"), border = NA, axes = FALSE)
axis(1, at = seq(0, .25, .05))
mtext("Mean held-out Spearman correlation across 13 models", side = 1, line = 2.45)
text(rev(averages) + .005, positions, sprintf("%.2f", rev(averages)),
     adj = c(0, .5), col = ink)
dev.off()

# Heat map: discovery-half Spearman for every model x prompt setup (13 x 14 = 182 cells)
grid <- read.csv("outputs/discovery_grid.csv")
m <- tapply(grid$r_dS, list(grid$format, grid$model), identity)
stopifnot(!anyNA(m), all(dim(m) == c(14, 13)))
colnames(m) <- sub("-hf$", "", sub("^.*/", "", colnames(m)))
colnames(m) <- c("bert-base-uncased" = "BERT-base", "roberta-base" = "RoBERTa-base",
                 "roberta-large" = "RoBERTa-large", gpt2 = "GPT-2",
                 "Qwen2.5-0.5B" = "Qwen2.5-0.5B", "Qwen2.5-0.5B-Instruct" = "Qwen2.5-0.5B-Inst.",
                 "OLMo-1B" = "OLMo-1B", "SmolLM2-360M" = "SmolLM2-360M",
                 "pythia-70m" = "Pythia-70M", "pythia-160m" = "Pythia-160M",
                 "pythia-410m" = "Pythia-410M", "pythia-1b" = "Pythia-1B",
                 "pythia-1.4b" = "Pythia-1.4B")[colnames(m)]
m <- m[order(rowMeans(m)), order(colMeans(m))]
row_labels <- c(after_word = "After the word…", given_word = "Given the word…",
                first_word = "If the first word is…", pair_begins = "The pair begins with…",
                and = "X and Y", ldt = "Word 1 / Word 2", list = "Words: X, Y",
                "bare+k10" = "Bare + 10 examples", "bare+k5" = "Bare + 5 examples",
                "bare+k2" = "Bare + 2 examples", comma = "X, Y", newline = "Line break",
                assoc = "X makes me think of Y", bare = "Bare pair")
pal <- colorRampPalette(c("#FFFFFF", "#BFE6E1", teal, "#0B3C49"))(100)
lim <- .36
svg("presentation/prompt-heatmap.svg", width = 10.4, height = 6.6, pointsize = 17)
par(mar = c(8.4, 12.2, 0.3, 0.4), fg = ink, col.axis = ink, family = "sans")
image(seq_len(ncol(m)), seq_len(nrow(m)), t(pmin(pmax(m, 0), lim)), col = pal,
      zlim = c(0, lim), axes = FALSE, xlab = "", ylab = "")
abline(h = seq(.5, nrow(m) + .5), v = seq(.5, ncol(m) + .5), col = "white", lwd = 1.5)
axis(2, at = seq_len(nrow(m)), labels = row_labels[rownames(m)], las = 1, tick = FALSE, line = -.6)
axis(1, at = seq_len(ncol(m)), labels = colnames(m), las = 2, tick = FALSE, line = -.6)
cells <- expand.grid(col = seq_len(ncol(m)), row = seq_len(nrow(m)))
vals <- m[cbind(cells$row, cells$col)]
text(cells$col, cells$row, sub("^-?0\\.", ".", sub("^-0.00$", "0.00", sprintf("%.2f", vals))), cex = .95,
     col = ifelse(vals > .20, "white", ink))
# Box the top setup on the discovery half
best <- which(m == max(m), arr.ind = TRUE)
wr <- best[1, "row"]; wc <- best[1, "col"]
stopifnot(rownames(m)[wr] == "after_word", colnames(m)[wc] == "OLMo-1B")
rect(wc - .5, wr - .5, wc + .5, wr + .5, border = "#F2A900", lwd = 5)
dev.off()

# Side-by-side distributions in raw units, selected setup, all 1,000 items
both <- items[is.finite(items$dS) & is.finite(items$avgZ_prime), ]
print(c(n = nrow(both), pos_human = mean(both$avgZ_prime > 0), pos_model = mean(both$dS > 0)))
amber <- "#F2A900"
svg("presentation/prompt-distributions.svg", width = 10.4, height = 4.4, pointsize = 17)
par(mfrow = c(1, 2), mar = c(4.4, 4.6, 2.2, 0.8), fg = ink, col.axis = ink, col.lab = ink,
    col.main = ink, las = 1, bty = "l", mgp = c(2.7, .6, 0), family = "sans")
panel <- function(x, col, main, xlab) {
  hist(x, breaks = 30, col = adjustcolor(col, .75), border = "white", main = main,
       xlab = xlab, ylab = "Items", font.main = 1)
  abline(v = 0, lty = 2, col = ink)
}
panel(both$avgZ_prime, teal, "Human", "Priming (standardized RT difference)")
panel(both$dS, amber, "Model", "Priming (surprisal difference, nats)")
dev.off()

# Semantic Priming Project: the eight discovery-selected setups x cue type x delay (all targets)
spp <- subset(read.csv("outputs/spp_validation.csv"), subset == "all targets")
spp$col <- paste(spp$prime, spp$soa)
sm <- tapply(spp$r, list(spp$config, spp$col), identity)
sm <- sm[, c("first short", "first long", "other short", "other long")]
stopifnot(!anyNA(sm), nrow(sm) == 8)
print(round(colMeans(sm), 3))
model_names <- c("EleutherAI/pythia-1b" = "Pythia-1B", "HuggingFaceTB/SmolLM2-360M" = "SmolLM2-360M",
                 "Qwen/Qwen2.5-0.5B-Instruct" = "Qwen2.5-0.5B-Inst.", "Qwen/Qwen2.5-0.5B" = "Qwen2.5-0.5B",
                 "allenai/OLMo-1B-hf" = "OLMo-1B")
frame_names <- c(ldt = "Word 1 / Word 2", newline = "Line break", after_word = "After the word…",
                 pair_begins = "The pair begins with…", and = "X and Y")
parts <- strsplit(rownames(sm), " \\| ")
rownames(sm) <- vapply(parts, function(p) paste0(model_names[p[1]], " · ", frame_names[p[2]]), "")
sm <- sm[order(rowMeans(sm), decreasing = TRUE), ]
sm <- rbind(sm, Mean = colMeans(sm))
sm <- sm[nrow(sm):1, ]   # image() draws row 1 at the bottom, so Mean ends up last (bottom)
nr <- nrow(sm); nc <- ncol(sm)
svg("presentation/spp-heatmap.svg", width = 9.4, height = 5.6, pointsize = 17)
par(mar = c(5.6, 15.2, 0.4, 0.4), fg = ink, col.axis = ink, family = "sans")
image(seq_len(nc), seq_len(nr), t(pmin(pmax(sm, 0), lim)), col = pal, zlim = c(0, lim),
      axes = FALSE, xlab = "", ylab = "")
abline(h = seq(.5, nr + .5), v = seq(.5, nc + .5), col = "white", lwd = 1.5)
abline(h = 1.5, col = "white", lwd = 6)
axis(2, at = seq_len(nr), labels = rownames(sm), las = 1, tick = FALSE, line = -.6)
mtext(rep(c("Short", "Long"), 2), side = 1, at = 1:4, line = 0.5)
mtext(c("First associate", "Other associate"), side = 1, at = c(1.5, 3.5), line = 2.4, font = 2)
cells <- expand.grid(col = seq_len(nc), row = seq_len(nr))
vals <- sm[cbind(cells$row, cells$col)]
text(cells$col, cells$row, sub("^-?0\\.", ".", sprintf("%.2f", vals)), cex = .95,
     col = ifelse(vals > .20, "white", ink), font = ifelse(cells$row == 1, 2, 1))
dev.off()

# English vs German: best English model, best German model and mGPT, each on English and German words
# (confirmation half). Real values come from outputs/crosslang_final.csv (rows: model, columns: EN, DE).
# Until that file exists, the figure is drawn with PLACEHOLDER numbers and stamped as such.
final <- "outputs/crosslang_final.csv"
placeholder <- !file.exists(final)
if (placeholder) {
  cm <- rbind("Best English model · OLMo-1B" = c(.24, .06),
              "Best German model · LLäMmlein 1B" = c(.09, .12),
              "Multilingual · mGPT" = c(.17, .09))
} else {
  cm <- as.matrix(read.csv(final, row.names = 1, check.names = FALSE))
}
colnames(cm) <- c("EN", "DE")
cm <- cm[nrow(cm):1, , drop = FALSE]   # first row at the top
svg("presentation/language-heatmap.svg", width = 9.6, height = 3.9, pointsize = 17)
par(mar = c(3.6, 16.4, 0.4, 0.4), fg = ink, col.axis = ink, family = "sans")
image(1:2, seq_len(nrow(cm)), t(pmin(pmax(cm, 0), lim)), col = pal, zlim = c(0, lim),
      axes = FALSE, xlab = "", ylab = "")
abline(h = seq(.5, nrow(cm) + .5), v = seq(.5, 2.5), col = "white", lwd = 1.5)
axis(2, at = seq_len(nrow(cm)), labels = rownames(cm), las = 1, tick = FALSE, line = -.6)
mtext(c("English words", "German words"), side = 1, at = 1:2, line = 0.4, font = 2)
cells <- expand.grid(col = 1:2, row = seq_len(nrow(cm)))
vals <- cm[cbind(cells$row, cells$col)]
text(cells$col, cells$row, sub("^-?0\\.", ".", sprintf("%.2f", vals)), cex = 1.1,
     col = ifelse(vals > .20, "white", ink))
# box the language-matched cells (English model on English words, German model on German words)
rect(.5, nrow(cm) - .5, 1.5, nrow(cm) + .5, border = "#F2A900", lwd = 4)
rect(1.5, nrow(cm) - 1.5, 2.5, nrow(cm) - .5, border = "#F2A900", lwd = 4)
if (placeholder) text(1.5, (nrow(cm) + 1) / 2, "FAKE NUMBERS", cex = 1.9,
                      col = adjustcolor("#C0392B", .75), srt = 10, font = 2)
dev.off()

# Ukrainian pilot: item priming scores and example norming ratings in one figure (aggregate values
# only; needs the private ../spaml2-private repository next to this one). Priming = unrelated-cue
# z RT minus related-cue z RT, from the current pipeline output as it is (17 targets currently appear
# twice; to be fixed in the processing).
uk_dir <- "../spaml2-private/formatted"
if (file.exists(file.path(uk_dir, "priming", "uk_summary.csv"))) {
  usum <- read.csv(file.path(uk_dir, "priming", "uk_summary.csv"))
  rel <- usum[usum$type == "related", c("target", "M_Z")]; names(rel)[2] <- "rel_z"
  unr <- usum[usum$type == "unrelated", c("target", "M_Z")]; names(unr)[2] <- "unrel_z"
  ukp <- merge(rel, unr, by = "target")
  ukp$priming <- ukp$unrel_z - ukp$rel_z
  print(c(n = nrow(ukp), mean = mean(ukp$priming), sd = sd(ukp$priming), pos = mean(ukp$priming > 0)))

  # example ratings: eight words, six norms, colored within each norm across all rated words
  rt <- read.csv(file.path(uk_dir, "clean_ratings.csv"))
  dm <- read.csv(file.path(uk_dir, "clean_demos.csv"))[, c("file", "participant_uuid", "task")]
  rt <- merge(rt[!(rt$unknown %in% c(TRUE, "TRUE", 1)), c("file", "participant_uuid", "word", "score")], dm,
              by = c("file", "participant_uuid"))
  wm <- tapply(rt$score, list(rt$word, rt$task), mean)
  wm <- wm[complete.cases(wm), ]
  pick <- c("акула" = "shark", "бекон" = "bacon", "бетмен" = "Batman", "вампір" = "vampire",
            "дракон" = "dragon", "шоколад" = "chocolate", "смерть" = "death", "любов" = "love")
  stopifnot(all(names(pick) %in% rownames(wm)))
  cols <- c(aoa = "AoA", familiar = "Fam", image = "Img", concrete = "Conc", valence = "Val", arousal = "Aro")
  ex <- wm[names(pick), names(cols)]
  scaled <- sapply(names(cols), function(k) {
    lo <- quantile(wm[, k], .02); hi <- quantile(wm[, k], .98)
    pmin(pmax((ex[, k] - lo) / (hi - lo), 0), 1)
  })
  ex <- ex[nrow(ex):1, ]; scaled <- scaled[nrow(scaled):1, ]

  svg("presentation/ukrainian-pilot.svg", width = 12.4, height = 4.6, pointsize = 16)
  layout(matrix(1:2, 1), widths = c(1, 1.5))
  par(mar = c(4.6, 4.8, 2.4, 0.8), fg = ink, col.axis = ink, col.lab = ink, col.main = ink,
      las = 1, bty = "l", mgp = c(3, .6, 0), family = "sans")
  hist(ukp$priming, breaks = 30, col = adjustcolor(teal, .75), border = "white",
       main = "Priming scores", font.main = 1,
       xlab = "Priming (standardized RT difference)", ylab = "Pairs")
  abline(v = 0, lty = 2, col = ink)
  par(mar = c(4.6, 9.8, 2.4, 0.4), bty = "n")
  image(seq_len(ncol(ex)), seq_len(nrow(ex)), t(scaled), col = pal, zlim = c(0, 1),
        axes = FALSE, xlab = "", ylab = "", main = "Example rating norms", font.main = 1)
  abline(h = seq(.5, nrow(ex) + .5), v = seq(.5, ncol(ex) + .5), col = "white", lwd = 1.5)
  axis(2, at = seq_len(nrow(ex)), labels = paste0(rownames(ex), " · ", pick[rownames(ex)]),
       las = 1, tick = FALSE, line = -.6)
  axis(1, at = seq_len(ncol(ex)), labels = cols, las = 1, tick = FALSE, line = -.4)
  cells <- expand.grid(col = seq_len(ncol(ex)), row = seq_len(nrow(ex)))
  text(cells$col, cells$row, sprintf("%.1f", ex[cbind(cells$row, cells$col)]), cex = .9,
       col = ifelse(scaled[cbind(cells$row, cells$col)] > .55, "white", ink))
  dev.off()
}
