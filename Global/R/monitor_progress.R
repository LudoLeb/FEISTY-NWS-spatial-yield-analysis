# =============================================================================
# MONITOR THE INTERMEDIATE FIVE-DEGREE PRODUCTION RUN
# =============================================================================
# Run this script from the Global project root while script 05 is running.
# It reads checkpoints and the log, then writes an auto-refreshing HTML page.
# =============================================================================

# ---- 1. Paths and expected number of chunks ---------------------------------

project_dir <- normalizePath(".", winslash = "/", mustWork = TRUE)
checkpoint_dir <- file.path(project_dir, "outputs", "checkpoints", "production_5deg")
log_file <- file.path(project_dir, "outputs", "logs", "global_production_5deg.log")
selected_file <- file.path(project_dir, "outputs", "tables", "global_production_5deg_selected_cells.csv")
progress_file <- file.path(project_dir, "outputs", "global_simulation_progress.html")
final_file <- file.path(project_dir, "outputs", "tables", "global_production_5deg_F040_SSB_constraint.csv")

chunk_size <- 40L
total_chunks <- if (file.exists(selected_file)) {
  ceiling(nrow(read.csv(selected_file, check.names = FALSE)) / chunk_size)
} else {
  55L
}

# ---- 2. Helper function ------------------------------------------------------

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

# ---- 3. Build the progress page ---------------------------------------------

write_progress <- function() {
  checkpoints <- list.files(checkpoint_dir, pattern = "^chunk_[0-9]+\\.rds$", full.names = TRUE)
  completed <- length(checkpoints)
  percent <- min(100, 100 * completed / total_chunks)

  log_lines <- if (file.exists(log_file)) readLines(log_file, warn = FALSE) else character()
  useful_lines <- grep("^(Chunk|Running chunk|Production run completed)", log_lines, value = TRUE)
  latest <- if (length(useful_lines)) tail(useful_lines, 1) else "Initialisation du calcul"

  complete <- completed >= total_chunks && file.exists(final_file)
  status <- if (complete) "Calcul terminé" else "Calcul en cours"
  status_colour <- if (complete) "#18864b" else "#1769aa"
  refresh <- if (complete) "" else '<meta http-equiv="refresh" content="20">'
  updated <- format(Sys.time(), "%d/%m/%Y à %H:%M:%S")

  html <- sprintf(
'<!doctype html>
<html lang="fr">
<head>
  <meta charset="utf-8">
  %s
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Progression FEISTY global</title>
  <style>
    body { margin: 0; background: #f3f6f8; color: #24313a; font-family: Arial, Helvetica, sans-serif; }
    main { max-width: 720px; margin: 8vh auto; padding: 0 24px; }
    .card { background: white; border-radius: 14px; padding: 34px; box-shadow: 0 8px 30px rgba(25,45,60,.10); }
    h1 { margin: 0 0 8px; font-size: 28px; }
    .status { color: %s; font-weight: 700; margin-bottom: 28px; }
    .numbers { display: flex; justify-content: space-between; align-items: end; margin-bottom: 10px; }
    .percent { font-size: 38px; font-weight: 700; }
    .chunks { color: #5d6b75; }
    .track { height: 22px; background: #e2e9ed; border-radius: 999px; overflow: hidden; }
    .bar { width: %.3f%%; height: 100%%; background: %s; border-radius: 999px; transition: width .4s ease; }
    .latest { margin-top: 25px; padding: 15px 17px; background: #f4f7f9; border-radius: 9px; }
    .note { margin-top: 20px; color: #667680; font-size: 14px; line-height: 1.5; }
  </style>
</head>
<body>
  <main>
    <div class="card">
      <h1>Simulation globale FEISTY</h1>
      <div class="status">%s</div>
      <div class="numbers">
        <div class="percent">%.1f %%</div>
        <div class="chunks">%d lots terminés sur %d</div>
      </div>
      <div class="track"><div class="bar"></div></div>
      <div class="latest">%s</div>
      <div class="note">Dernière mise à jour : %s<br>Cette page se rafraîchit automatiquement toutes les 20 secondes.</div>
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
    html_escape(latest),
    updated
  )

  writeLines(html, progress_file, useBytes = TRUE)
  complete
}

repeat {
  if (write_progress()) break
  Sys.sleep(20)
}
