library(dplyr)
library(purrr)
library(stringr)
library(tidyr)
library(quanteda)
library(tibble)
library(pangoling)

# same length requirement colthearts N
calc_neighbors <- function(DF, subs_DF) {
  # --- Prep lexicon ---
  lex <- subs_DF %>%
    mutate(
      unigram = str_to_lower(str_trim(unigram)),
      len = nchar(unigram)
    ) %>%
    group_by(unigram) %>%
    summarise(freq = sum(unigram_freq), .groups = "drop")
  
  lex_words <- lex$unigram
  lex_len   <- setNames(nchar(lex_words), lex_words)
  freq_lookup <- setNames(lex$freq, lex$unigram)
  lex_words <- lex_words[nchar(lex_words) > 0]
  
  # hash-like env for fast lookup
  lex_env <- new.env(parent = emptyenv())
  invisible(lapply(lex_words, function(w) assign(w, TRUE, envir = lex_env)))
  in_lex <- function(x) !vapply(mget(x, ifnotfound = NA, envir = lex_env), is.na, logical(1))
  
  # neighbor generator
  gen_sub_neighbors <- function(w) {
    ch <- str_split(w, "", simplify = TRUE)
    L  <- ncol(ch)
    out <- character(L * (26 - 1))
    k <- 1L
    for (i in seq_len(L)) {
      orig <- ch[1, i]
      for (a in letters) {
        if (a != orig) {
          ch[1, i] <- a
          out[k] <- str_c(ch, collapse = "")
          k <- k + 1L
        }
      }
      ch[1, i] <- orig
    }
    unique(out)
  }
  
  # wrapper for a single word
  calc_for_word <- function(w, L) {
    if (is.na(w) || w == "" || L <= 0) return(c(N = 0L, N_freq = 0))
    if (!any(lex_len == L)) return(c(N = 0L, N_freq = 0))
    cand <- gen_sub_neighbors(w)
    keep <- in_lex(cand)
    nbrs <- cand[keep & cand != w]
    N <- length(nbrs)
    N_freq <- if (N == 0) 0 else sum(freq_lookup[nbrs], na.rm = TRUE)
    c(N = N, N_freq = N_freq)
  }
  
  DF %>%
    mutate(
      word = str_to_lower(str_trim(word)),
      length = if ("length" %in% names(DF)) length else nchar(word)
    ) %>%
    mutate(
      res = map2(word, length, calc_for_word)
    ) %>%
    mutate(
      orthographic_neighbors = map_int(res, 1),
      orthographic_neighbors_freq = map_dbl(res, 2)
    ) %>%
    select(-res)
}

# Helper: extract adjacent letter bigrams from a single word
.bigrams_of <- function(w) {
  w <- str_to_lower(str_trim(w))
  if (is.na(w) || nchar(w) < 2) return(character(0))
  starts <- 1:(nchar(w) - 1)
  str_sub(w, starts, starts + 1)
}

# Main: add bigram frequency features to DF using subs_DF as the corpus
add_bigram_frequency <- function(DF, subs_DF, weight = c("unigram_freq", "count")) {
  weight <- match.arg(weight)
  
  # --- Prep lexicon (lowercase, trim, drop blanks) ---
  lex <- subs_DF %>%
    transmute(
      unigram = str_to_lower(str_trim(unigram)),
      freq = as.numeric(unigram_freq)
    ) %>%
    filter(!is.na(unigram), unigram != "", !is.na(freq))
  
  # --- Build bigram table from the lexicon ---
  # Each word contributes its internal bigrams; we weight by its frequency
  bigram_tbl <- lex %>%
    mutate(bigrams = map(unigram, .bigrams_of)) %>%
    filter(map_int(bigrams, length) > 0) %>%
    unnest_longer(bigrams, values_to = "bigram") %>%
    {
      if (weight == "unigram_freq") {
        group_by(., bigram) %>% summarise(bigram_freq = sum(freq), .groups = "drop")
      } else {
        # Count how many *words* contain the bigram (binary per word)
        distinct(., unigram, bigram) %>%
          count(bigram, name = "bigram_freq")
      }
    }
  
  # --- Score each target word in DF ---
  DF %>%
    mutate(
      word = str_to_lower(str_trim(word)),
      length = if ("length" %in% names(.)) length else nchar(word),
      bigrams = map(word, .bigrams_of)
    ) %>%
    left_join(
      # sum/mean across that word's bigrams
      tibble(word = .$word, bigrams = .$bigrams) %>%
        unnest_longer(bigrams, values_to = "bigram", keep_empty = TRUE) %>%
        left_join(bigram_tbl, by = "bigram") %>%
        group_by(word) %>%
        summarise(
          bigram_freq_sum  = sum(replace_na(bigram_freq, 0)),
          bigram_freq_mean = if (n() > 0) mean(replace_na(bigram_freq, 0)) else 0,
          .groups = "drop"
        ),
      by = "word"
    ) %>%
    mutate(
      bigram_freq_sum  = replace_na(bigram_freq_sum, 0),
      bigram_freq_mean = replace_na(bigram_freq_mean, 0)
    ) %>%
    select(-bigrams)
}

add_orthographic_levenshtein <- function(DF, subs_DF, k = 20, max_len_diff = 3) {
  # --- Prep lexicon ---
  lex <- subs_DF %>%
    transmute(
      unigram = str_to_lower(str_trim(unigram))
    ) %>%
    filter(!is.na(unigram), unigram != "") %>%
    distinct() %>%
    mutate(len = nchar(unigram))
  
  # quick split by length for fast lookup
  lex_by_len <- split(lex$unigram, lex$len)
  all_lens <- as.integer(names(lex_by_len))
  
  # helper to fetch candidate set within ± max_len_diff
  get_candidates <- function(L) {
    lens <- all_lens[abs(all_lens - L) <= max_len_diff]
    if (length(lens) == 0) character(0) else unlist(lex_by_len[as.character(lens)], use.names = FALSE)
  }
  
  # per-word computation
  compute_OLDk <- function(w) {
    w <- str_to_lower(str_trim(w))
    if (is.na(w) || w == "") return(c(OLDk = NA_real_, lev_min = NA_real_, lev1_count = 0))
    L <- nchar(w)
    cands <- get_candidates(L)
    # exclude exact self if present
    cands <- cands[cands != w]
    if (length(cands) == 0) return(c(OLDk = NA_real_, lev_min = NA_real_, lev1_count = 0))
    
    # distances via adist (Levenshtein; cost 1 for ins/del/sub)
    d <- as.vector(adist(w, cands))
    # sort once
    ord <- order(d)
    d_sorted <- d[ord]
    # OLDk: mean of the k smallest (or all if fewer than k)
    kk <- min(k, length(d_sorted))
    OLDk <- mean(d_sorted[seq_len(kk)])
    lev_min <- d_sorted[1]
    lev1_count <- sum(d == 1L)
    c(OLDk = OLDk, lev_min = lev_min, lev1_count = lev1_count)
  }
  
  DF %>%
    mutate(word = str_to_lower(str_trim(word))) %>%
    mutate(tmp = map(word, compute_OLDk)) %>%
    mutate(
      OLDk = map_dbl(tmp, 1),
      lev_min = map_dbl(tmp, 2),
      lev1_count = map_int(tmp, 3)
    ) %>%
    select(-tmp)
}

# deps: dplyr, purrr, stringr
add_phonographic_levenshtein <- function(DF, k = 20, max_len_diff = 3,
                                         word_col = "word", phon_col = "phon") {
  # --- Prep phoneme lexicon from DF itself ---
  lex <- DF %>%
    transmute(
      word = stringr::str_to_lower(stringr::str_trim(.data[[word_col]])),
      phon = stringr::str_trim(.data[[phon_col]])
    ) %>%
    dplyr::filter(!is.na(phon), phon != "") %>%
    dplyr::distinct() %>%
    dplyr::mutate(len = lengths(stringr::str_split(phon, "\\s+")))
  
  # split by phoneme length for quick candidate lookup
  lex_by_len <- split(lex$phon, lex$len)
  all_lens <- as.integer(names(lex_by_len))
  get_candidates <- function(L) {
    lens <- all_lens[abs(all_lens - L) <= max_len_diff]
    if (length(lens) == 0) character(0) else unlist(lex_by_len[as.character(lens)], use.names = FALSE)
  }
  
  # distance calc on space-separated phoneme strings
  compute_PLDk <- function(w) {
    w <- stringr::str_to_lower(stringr::str_trim(w))
    phon_w <- lex$phon[match(w, lex$word)]
    if (is.na(phon_w)) return(c(PLDk = NA_real_, lev_min = NA_real_, lev1_count = NA_integer_))
    
    L <- lengths(stringr::str_split(phon_w, "\\s+"))
    cands <- get_candidates(L)
    cands <- cands[cands != phon_w]
    if (length(cands) == 0) return(c(PLDk = NA_real_, lev_min = NA_real_, lev1_count = 0))
    
    # Levenshtein on the phoneme strings (simple & fast)
    d <- as.vector(adist(phon_w, cands))
    ord <- order(d)
    d_sorted <- d[ord]
    kk <- min(k, length(d_sorted))
    
    PLDk <- mean(d_sorted[seq_len(kk)])
    lev_min <- d_sorted[1]
    lev1_count <- sum(d == 1L)
    c(PLDk = PLDk, lev_min = lev_min, lev1_count = lev1_count)
  }
  
  DF %>%
    dplyr::mutate(!!word_col := stringr::str_to_lower(stringr::str_trim(.data[[word_col]]))) %>%
    dplyr::mutate(`.tmp_pld` = purrr::map(.data[[word_col]], compute_PLDk)) %>%
    dplyr::mutate(
      PLDk = purrr::map_dbl(`.tmp_pld`, 1),
      phon_lev_min = purrr::map_dbl(`.tmp_pld`, 2),
      phon_lev1_count = purrr::map_int(`.tmp_pld`, 3)
    ) %>%
    dplyr::select(-`.tmp_pld`)
}

# ---------- 2) Hoffman SemD per word + prop across contexts ----------
# SemD(word) = -log( mean_{(i<j) in contexts(word)} cosine(ctx_i, ctx_j) ), cap contexts at 2000
# A: m x k context embeddings for a single word's contexts
# returns mean of off-diagonal cosines without storing the full m x m Gram
.fast_mean_pairwise_cosine <- function(A, block_rows = 1500) {
  m <- nrow(A)
  if (m <= 1) return(NA_real_)
  # row-normalize
  rn <- sqrt(rowSums(A * A))
  A  <- A / rn
  
  # small m → just use tcrossprod once
  if (m <= block_rows) {
    G <- tcrossprod(A)                # m x m
    return(mean(G[upper.tri(G)], na.rm = TRUE))
  }
  
  # large m → accumulate upper-tri in blocks
  # partition rows into blocks; sum upper-tri entries across blocks
  idx <- split(seq_len(m), ceiling(seq_len(m) / block_rows))
  total_sum <- 0
  total_pairs <- 0
  
  for (bi in seq_along(idx)) {
    Ai <- A[idx[[bi]], , drop = FALSE]
    # diagonal block: take upper triangle within this block
    Gi <- tcrossprod(Ai)
    if (nrow(Gi) > 1) {
      ut <- Gi[upper.tri(Gi)]
      total_sum   <- total_sum   + sum(ut, na.rm = TRUE)
      total_pairs <- total_pairs + length(ut)
    }
    # off-diagonal blocks: only one triangle (i<j)
    if (bi < length(idx)) {
      for (bj in (bi + 1L):length(idx)) {
        Aj <- A[idx[[bj]], , drop = FALSE]
        Gij <- Ai %*% t(Aj)           # |Ai| x |Aj|
        total_sum   <- total_sum   + sum(Gij, na.rm = TRUE)
        total_pairs <- total_pairs + length(Gij)
      }
    }
  }
  if (total_pairs == 0) return(NA_real_)
  total_sum / total_pairs
}

# space: list returned by build_context_space(...): $X (contexts x terms, sparse), $context_vectors (contexts x k), $vocab
# terms: optional vector to evaluate (else full vocab)
# max_contexts: cap per word (Hoffman capped at 2,000)
# n_workers: parallel workers (set >1 on Linux/Mac; on Windows this code falls back to PSOCK)
compute_semd_hoffman_parallel <- function(
    space, terms = NULL, max_contexts = 2000, n_workers = parallel::detectCores(logical = TRUE) - 1
) {
  X     <- space$X                 # dgCMatrix, contexts x terms (weighted, nonzero pattern OK)
  Cvec  <- space$context_vectors   # contexts x k (dense)
  vocab <- space$vocab
  ndocs <- nrow(X)
  
  # which terms to process
  if (is.null(terms)) {
    terms <- vocab
  }
  map <- match(tolower(terms), tolower(vocab), nomatch = 0L)
  
  # work items (skip 0 matches up front but we’ll keep NA rows to preserve order)
  jobs <- which(map > 0L)
  
  # progress bar
  pb <- utils::txtProgressBar(min = 0, max = length(jobs), style = 3)
  pb_env <- environment()
  .tick <- function() utils::setTxtProgressBar(pb, get("i", envir = pb_env))
  
  # worker function (pure)
  worker <- function(j) {
    col_j <- map[j]
    # context indices where term occurs (nonzero in column col_j)
    idx <- which(X[, col_j] != 0)
    prop <- length(idx) / ndocs
    if (length(idx) == 0) return(list(term = terms[j], semd = NA_real_, prop = prop))
    
    if (length(idx) > max_contexts) {
      # deterministic sample for reproducibility across workers
      set.seed(1009 + j)
      idx <- sample(idx, max_contexts)
    }
    A <- Cvec[idx, , drop = FALSE]
    mean_cos <- .fast_mean_pairwise_cosine(A, block_rows = 1500)
    semd <- if (is.na(mean_cos) || mean_cos <= 0) NA_real_ else -log(mean_cos)
    list(term = terms[j], semd = semd, prop = prop)
  }
  
  results <- vector("list", length(terms))
  # serial fill for missing terms (kept NA)
  for (j in setdiff(seq_along(terms), jobs)) {
    results[[j]] <- list(term = terms[j], semd = NA_real_, prop = NA_real_)
  }
  
  # parallel over jobs
  os_is_windows <- tolower(Sys.info()[["sysname"]]) == "windows"
  i <- 0L
  if (!os_is_windows && n_workers > 1) {
    out <- parallel::mclapply(jobs, function(j) {
      res <- worker(j); res
    }, mc.cores = max(1L, n_workers), mc.preschedule = TRUE)
    # collect with progress (coarse)
    for (k in seq_along(jobs)) {
      results[[jobs[k]]] <- out[[k]]
      i <- i + 1L; .tick()
    }
  } else {
    # Windows or single-core: PSOCK with parLapply or plain lapply
    if (n_workers > 1 && os_is_windows) {
      cl <- parallel::makeCluster(n_workers)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      # export needed objects to workers
      parallel::clusterExport(cl, varlist = c("worker", "map", "X", "Cvec", "ndocs", "terms", ".fast_mean_pairwise_cosine"), envir = environment())
      out <- parallel::parLapply(cl, jobs, function(j) worker(j))
      for (k in seq_along(jobs)) {
        results[[jobs[k]]] <- out[[k]]
        i <- i + 1L; .tick()
      }
    } else {
      for (j in jobs) {
        results[[j]] <- worker(j)
        i <- i + 1L; .tick()
      }
    }
  }
  close(pb)
  
  # bind
  DT <- data.table::rbindlist(lapply(results, as.data.frame))
  data.table::setDT(DT)[, .(term, semd, prop)]
}

# Compute surprisal of the *last word* in each sentence, conditioned on its left context.
# log.p = 0.5 -> bits. Use 1 for nats.
last_word_surprisal <- function(sentences,
                                log.p = 0.5,
                                drop_trailing_punct = TRUE,
                                model = "gpt2") {
  word_lists <- str_split(sentences, "\\s+")
  df <- map2_dfr(word_lists, seq_along(word_lists),
                 ~ tibble(sent_n = paste0("s", .y), word = .x))
  
  if (drop_trailing_punct) {
    df <- df %>%
      mutate(word = str_replace(word, "[[:punct:]]+$", "")) %>%
      filter(word != "")
  }
  
  scored <- df %>%
    mutate(surprisal = causal_words_pred(
      word, by = sent_n, log.p = log.p, model = model
    ))
  
  scored %>%
    group_by(sent_n) %>%
    slice_tail(n = 1) %>%
    ungroup() %>%
    transmute(sent_n, last_word = word, surprisal)
}

# Priming effect for pairs: delta = surprisal_B - surprisal_A
# Positive delta => B is *more* surprising on its final word (A "primes" relative to B).
priming_delta <- function(sentence_A, sentence_B, log.p = 0.5, model = "gpt2") {
  stopifnot(length(sentence_A) == length(sentence_B))
  all_sent <- c(sentence_A, sentence_B)
  
  # give deterministic names to last_word_surprisal output
  names(all_sent) <- paste0("s", seq_along(all_sent))
  
  lw <- last_word_surprisal(all_sent, log.p = log.p, model = model)
  
  n <- length(sentence_A)
  
  tibble(
    pair_id = seq_len(n),
    sentence_A, sentence_B
  ) |>
    mutate(
      sent_n_A = paste0("s", seq_len(n)),          # s1..sn  -> A
      sent_n_B = paste0("s", n + seq_len(n))       # s(n+1)..s(2n) -> B
    ) |>
    left_join(lw, by = c("sent_n_A" = "sent_n")) |>
    rename(last_word_A = last_word, surprisal_A = surprisal) |>
    left_join(lw, by = c("sent_n_B" = "sent_n")) |>
    rename(last_word_B = last_word, surprisal_B = surprisal) |>
    mutate(priming_delta = surprisal_B - surprisal_A) |>
    select(-sent_n_A, -sent_n_B)
}

# find exact sub-sequence indices
find_subseq <- function(haystack_ids, needle_ids) {
  H <- as.integer(haystack_ids); N <- as.integer(needle_ids)
  hits <- integer(0)
  if (length(N) == 0 || length(H) < length(N)) return(hits)
  for (i in seq_len(length(H) - length(N) + 1L)) {
    if (all(H[i:(i + length(N) - 1L)] == N)) hits <- c(hits, i)
  }
  hits
}

get_embeddings <- function(word, sentence = NULL, tokzr, model,
                           occurrence = c("last", "first", "all"),
                           layer = NULL) {
  occurrence <- match.arg(occurrence)
  
  # --- Static (average over subtokens) ---
  word_ids   <- tokzr$encode(word, add_special_tokens = FALSE)
  emb_mat    <- model$transformer$wte$weight$detach()$cpu()$numpy()
  static_avg <- colMeans(emb_mat[as.integer(word_ids) + 1L, , drop = FALSE])
  
  # --- Contextual (optional) ---
  contextual_avg <- NULL
  if (!is.null(sentence)) {
    input_ids <- tokzr$encode(sentence, return_tensors = "pt")
    out <- model$transformer(input_ids, output_hidden_states = TRUE)
    hstates <- out$hidden_states
    target_layer <- if (is.null(layer)) length(hstates) else layer
    chosen <- hstates[[target_layer]]
    chosen_np <- chosen$detach()$cpu()$numpy()[1,,]  # [seq_len, hidden_dim]
    
    # Convert tensor -> numpy BEFORE indexing
    ids_np   <- input_ids$cpu()$numpy()
    sent_ids <- as.integer(reticulate::py_to_r(ids_np)[1, ])
    
    # --- NEW: try both encodings of the target word ---
    ids_space   <- tokzr$encode(paste0(" ", word), add_special_tokens = FALSE)
    ids_nospace <- tokzr$encode(word,              add_special_tokens = FALSE)
    
    starts_space   <- find_subseq(sent_ids, ids_space)
    starts_nospace <- find_subseq(sent_ids, ids_nospace)
    
    if (length(starts_space) > 0) {
      word_ids <- ids_space
      starts   <- starts_space
    } else if (length(starts_nospace) > 0) {
      word_ids <- ids_nospace
      starts   <- starts_nospace
    } else {
      starts <- integer(0)  # no match found
    }
    
    if (length(starts) > 0) {
      spans <- lapply(starts, function(s) s:(s + length(word_ids) - 1L))
      if (occurrence == "first") spans <- spans[1]
      if (occurrence == "last")  spans <- spans[length(spans)]
      mats <- lapply(spans, function(idx) chosen_np[idx, , drop = FALSE])
      contextual_avg <- Reduce(`+`, lapply(mats, colMeans)) / length(mats)
    }
  }
  
  list(
    static = static_avg,
    contextual = contextual_avg,
    layer_used = if (is.null(layer)) "last" else layer
  )
}

# --- assumes you already have: hf$tokzr, hf$model, find_subseq(), get_embeddings() ---

embed_dataframe <- function(df,
                            tokzr, model,
                            layer = NULL,
                            occurrence = "last",
                            word_col = "target_word_unique",
                            sent_rel_col = "sentence_related",
                            sent_unrel_col = "sentence_unrelated",
                            quiet = TRUE) {
  stopifnot(all(c(word_col, sent_rel_col, sent_unrel_col) %in% names(df)))
  
  # Determine hidden size from the model's embedding matrix
  emb_mat <- model$transformer$wte$weight$detach()$cpu()$numpy()
  hidden_size <- ncol(emb_mat)
  
  n <- nrow(df)
  # Preallocate (faster & predictable)
  static_mat   <- matrix(NA_real_, nrow = n, ncol = hidden_size)
  ctx_rel_mat  <- matrix(NA_real_, nrow = n, ncol = hidden_size)
  ctx_unrel_mat<- matrix(NA_real_, nrow = n, ncol = hidden_size)
  
  colnames(static_mat)    <- paste0("d", seq_len(hidden_size))
  colnames(ctx_rel_mat)   <- paste0("d", seq_len(hidden_size))
  colnames(ctx_unrel_mat) <- paste0("d", seq_len(hidden_size))
  
  for (i in seq_len(n)) {
    w   <- as.character(df[[word_col]][i])
    s_r <- as.character(df[[sent_rel_col]][i])
    s_u <- as.character(df[[sent_unrel_col]][i])
    
    # Static once per row
    emb_static <- try(
      get_embeddings(w, sentence = NULL,
                     tokzr = tokzr, model = model),
      silent = TRUE
    )
    if (!inherits(emb_static, "try-error") && !is.null(emb_static$static)) {
      static_mat[i, ] <- as.numeric(emb_static$static)
    }
    
    # Contextual in related
    emb_rel <- try(
      get_embeddings(w, sentence = s_r,
                     tokzr = tokzr, model = model,
                     layer = layer, occurrence = occurrence),
      silent = TRUE
    )
    if (!inherits(emb_rel, "try-error") && !is.null(emb_rel$contextual)) {
      ctx_rel_mat[i, ] <- as.numeric(emb_rel$contextual)
    }
    
    # Contextual in unrelated
    emb_unrel <- try(
      get_embeddings(w, sentence = s_u,
                     tokzr = tokzr, model = model,
                     layer = layer, occurrence = occurrence),
      silent = TRUE
    )
    if (!inherits(emb_unrel, "try-error") && !is.null(emb_unrel$contextual)) {
      ctx_unrel_mat[i, ] <- as.numeric(emb_unrel$contextual)
    }
    
    if (!quiet && i %% 50 == 0) message("Processed ", i, " / ", n)
  }
  
  # Return a tidy list; keep original ids for alignment
  rownames(static_mat)    <- rownames(df)
  rownames(ctx_rel_mat)   <- rownames(df)
  rownames(ctx_unrel_mat) <- rownames(df)
  
  list(
    static            = static_mat,
    context_related   = ctx_rel_mat,
    context_unrelated = ctx_unrel_mat,
    layer_used        = if (is.null(layer)) "last" else layer,
    occurrence        = occurrence
  )
}