# Generated from function body. Editing this file has no effect.
run_vina <- function (receptor = NULL, ligand = NULL, center = NULL, size = NULL, 
    config_paths = NULL, exhaustiveness = 8, output_dir = "docked", 
    vina_path = NULL, seed = 12345, cpu = NULL) 
{
    if (!is.null(cpu)) {
        if (!is.numeric(cpu) || cpu != as.integer(cpu) || cpu <= 
            0) {
            warning(sprintf("Provided CPU value '%s' is not a valid positive integer, using auto-detection", 
                cpu))
            cpu <- NULL
        }
    }
    if (is.null(vina_path)) {
        vina_path <- Sys.getenv("GLUEDOCK_VINA_PATH", unset = NA)
        if (is.na(vina_path)) {
            stop("vina_path parameter is required and GLUEDOCK_VINA_PATH environment variable not set. Please provide vina_path parameter or set GLUEDOCK_VINA_PATH environment variable.")
        }
    }
    if (!file.exists(vina_path)) {
        stop("AutoDock Vina executable does not exist: ", vina_path)
    }
    if (is.null(config_paths)) {
        if (is.null(receptor) || is.null(ligand) || is.null(center) || 
            is.null(size)) {
            stop("In direct mode, receptor, ligand, center, and size must be provided")
        }
        if (!file.exists(receptor)) {
            stop("Receptor file does not exist: ", receptor)
        }
        if (!file.exists(ligand)) {
            stop("Ligand file does not exist: ", ligand)
        }
        if (length(center) != 3) {
            stop("Center must be a numeric vector of length 3")
        }
        if (length(size) != 3) {
            stop("Size must be a numeric vector of length 3")
        }
        receptor <- normalizePath(receptor, winslash = "/", mustWork = FALSE)
        ligand <- normalizePath(ligand, winslash = "/", mustWork = FALSE)
        if (!is.null(seed)) {
            if (!is.numeric(seed) || seed != as.integer(seed) || 
                seed <= 0) {
                warning(sprintf("Provided seed value '%s' is not a valid positive integer, using default value %d", 
                  seed, 12345))
                seed <- 12345
            }
        }
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        rec_id <- tolower(tools::file_path_sans_ext(basename(receptor)))
        lig_id <- tolower(tools::file_path_sans_ext(basename(ligand)))
        out_file <- file.path(output_dir, paste0(rec_id, "_", 
            lig_id, ".pdbqt"))
        cmd <- sprintf("%s --receptor %s --ligand %s --center_x %f --center_y %f --center_z %f --size_x %f --size_y %f --size_z %f --exhaustiveness %d --out %s --seed %d", 
            shQuote(vina_path), shQuote(receptor), shQuote(ligand), 
            center[1], center[2], center[3], size[1], size[2], 
            size[3], exhaustiveness, shQuote(out_file), seed)
        if (!is.null(cpu)) {
            cmd <- paste0(cmd, sprintf(" --cpu %d", as.integer(cpu)))
        }
        message("Running direct docking command:\n  ", cmd)
        tryCatch(system(cmd, intern = TRUE), error = function(e) {
            warning("Error executing vina: ", e$message)
        })
    }
    else {
        config_files <- character(0)
        for (path in config_paths) {
            if (dir.exists(path)) {
                files <- list.files(path = path, pattern = ".*_config\\.txt$", 
                  full.names = TRUE)
                config_files <- c(config_files, files)
            }
            else if (file.exists(path)) {
                config_files <- c(config_files, path)
            }
            else {
                warning("Config path does not exist: ", path)
            }
        }
        if (length(config_files) == 0) {
            stop("No valid config files found")
        }
        if (!is.null(seed)) {
            if (!is.numeric(seed) || seed != as.integer(seed) || 
                seed <= 0) {
                warning(sprintf("Provided seed value '%s' is not a valid positive integer, using default value %d", 
                  seed, 12345))
                seed <- 12345
            }
        }
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }
        results <- list()
        for (i in seq_along(config_files)) {
            config_file <- config_files[i]
            config_base <- tools::file_path_sans_ext(basename(config_file))
            config_base <- sub("_config$", "", config_base)
            out_file <- file.path(output_dir, paste0(config_base, 
                ".pdbqt"))
            cmd <- sprintf("%s --config %s --out %s --seed %d", 
                shQuote(vina_path), shQuote(config_file), shQuote(out_file), 
                seed)
            if (!is.null(cpu)) {
                cmd <- paste0(cmd, sprintf(" --cpu %d", as.integer(cpu)))
            }
            message(sprintf("[#%d/%d] Running config docking command:\n  %s", 
                i, length(config_files), cmd))
            tryCatch(system(cmd, intern = TRUE), error = function(e) {
                warning("Error executing vina with config ", 
                  config_file, ": ", e$message)
            })
        }
        message(sprintf("Processed %d config files", length(config_files)))
    }
}
