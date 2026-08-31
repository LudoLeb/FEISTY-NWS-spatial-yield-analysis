# =============================================================================
# MONITOR THE COMPLETE ONE-DEGREE PRODUCTION RUN
# =============================================================================
# Run this script from the Global project root while script 06 is running.
# It reads checkpoints and logs only; it does not modify FEISTY results.
# The generated HTML page refreshes automatically to show progress.
# =============================================================================

# ---- 1. Paths and expected number of chunks ---------------------------------

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
checkpoint_dir <- file.path(project_dir, "outputs", "checkpoints", "full_grid")
log_file <- file.path(project_dir, "outputs", "logs", "global_full_grid.log")
grid_file <- file.path(project_dir, "outputs", "tables", "global_full_grid_cells.csv")
final_file <- file.path(project_dir, "outputs", "tables", "global_full_grid_F040_SSB_constraint.csv")
complete_marker <- file.path(project_dir, "outputs", "checkpoints", "full_grid", "COMPLETE.txt")
progress_file <- file.path(project_dir, "outputs", "global_full_grid_progress.html")
visible_progress_file <- Sys.getenv(
  "FEISTY_GLOBAL_VISIBLE_PROGRESS",
  unset = progress_file
)

chunk_size <- 70L
total_chunks <- if (file.exists(grid_file)) {
  ceiling(nrow(read.csv(grid_file, check.names = FALSE)) / chunk_size)
} else {
  586L
}

# ---- 2. Helper functions -----------------------------------------------------

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

format_remaining <- function(seconds) {
  if (!is.finite(seconds) || seconds <= 0) return("estimation en cours")
  hours <- seconds / 3600
  if (hours >= 24) {
    sprintf("environ %.1f jours restants", hours / 24)
  } else if (hours >= 1) {
    sprintf("environ %.1f heures restantes", hours)
  } else {
    sprintf("environ %d minutes restantes", round(seconds / 60))
  }
}

# ---- 3. Build the progress page ---------------------------------------------

write_progress <- function() {
  checkpoints <- list.files(checkpoint_dir, pattern = "^chunk_[0-9]+\\.rds$", full.names = TRUE)
  completed <- length(checkpoints)
  percent <- min(100, 100 * completed / total_chunks)

  log_lines <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character()
  useful_lines <- grep(
    "^(Chunk|Running chunk|Full-grid run completed)",
    log_lines,
    value = TRUE
  )
  latest <- if (length(useful_lines)) tail(useful_lines, 1) else "Initialisation du calcul"

  remaining <- "estimation en cours"
  if (length(checkpoints) >= 3L && completed < total_chunks) {
    times <- sort(as.numeric(file.info(checkpoints)$mtime))
    intervals <- diff(times)
    typical <- median(intervals[intervals > 0], na.rm = TRUE)
    remaining <- format_remaining((total_chunks - completed) * typical)
  }

  complete <- completed >= total_chunks && file.exists(final_file) && file.exists(complete_marker)
  status <- if (complete) {
    "Calcul terminé et validé"
  } else if (completed >= total_chunks) {
    "Consolidation et validation finales"
  } else if (completed == total_chunks - 1L) {
    "Réparation du dernier lot incomplet"
  } else {
    "Calcul en cours"
  }
  if (!complete && completed >= total_chunks - 1L) {
    latest <- "Recalcul du lot incomplet puis reconstruction des résultats finaux"
  }
  status_colour <- if (complete) "#18864b" else "#1769aa"
  refresh <- if (complete) "" else '<meta http-equiv="refresh" content="30">'
  updated <- format(Sys.time(), "%d/%m/%Y à %H:%M:%S")

  html <- sprintf(
'<!doctype html>
<html lang="fr">
<head>
  <meta charset="utf-8">
  %s
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Progression FEISTY — grille mondiale complète</title>
  <style>
    body { margin: 0; background: #f3f6f8; color: #24313a; font-family: Arial, Helvetica, sans-serif; }
    main { max-width: 760px; margin: 8vh auto; padding: 0 24px; }
    .card { background: white; border-radius: 14px; padding: 34px; box-shadow: 0 8px 30px rgba(25,45,60,.10); }
    h1 { margin: 0 0 6px; font-size: 28px; }
    .subtitle { color: #667680; margin-bottom: 10px; }
    .status { color: %s; font-weight: 700; margin-bottom: 28px; }
    .numbers { display: flex; justify-content: space-between; align-items: end; margin-bottom: 10px; }
    .percent { font-size: 38px; font-weight: 700; }
    .chunks { color: #5d6b75; text-align: right; }
    .track { height: 22px; background: #e2e9ed; border-radius: 999px; overflow: hidden; }
    .bar { width: %.4f%%; height: 100%%; background: %s; border-radius: 999px; transition: width .4s ease; }
    .remaining { margin-top: 15px; font-weight: 700; color: #44545f; }
    .latest { margin-top: 20px; padding: 15px 17px; background: #f4f7f9; border-radius: 9px; }
    .note { margin-top: 20px; color: #667680; font-size: 14px; line-height: 1.5; }
  </style>
</head>
<body>
  <main>
    <div class="card">
      <h1>Simulation mondiale FEISTY</h1>
      <div class="subtitle">Grille complète de 41 008 cellules — 2015 à 2087</div>
      <div class="status">%s</div>
      <div class="numbers">
        <div class="percent">%.1f %%</div>
        <div class="chunks">%d lots terminés sur %d</div>
      </div>
      <div class="track"><div class="bar"></div></div>
      <div class="remaining">%s</div>
      <div class="latest">%s</div>
      <div class="note">Dernière mise à jour : %s<br>Cette page se rafraîchit automatiquement toutes les 30 secondes.</div>
    </div>
  </main>
</body>
</html>',
    refresh,
    status_colour,
    percent,
    status_colour,
    status,
    percent,
    completed,
    total_chunks,
    remaining,
    html_escape(latest),
    updated
  )

  writeLines(html, progress_file, useBytes = TRUE)
  if (dir.exists(dirname(visible_progress_file))) {
    writeLines(html, visible_progress_file, useBytes = TRUE)
  }
  complete
}

repeat {
  if (write_progress()) break
  Sys.sleep(30)
}
