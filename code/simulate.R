## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Simulate ODE system with posteriors from STAN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

## Defines simulate_model(model_name, n_draws, internal_cores) which:
##   * loads posterior_draws_<model_name>.RDA,
##   * (optionally) subsamples to `n_draws` posterior draws,
##   * extracts the ODE body from stan_model_<model_name>.stan,
##   * integrates the ODE over a grid of initial conditions and posterior draws,
##   * writes tmp/ODE_kelp_<level>_<model_name>.RDA for each entry of A.level.
##
## When invoked from inside parallel::parLapply (parallel_models = TRUE in
## RunMe.R), pass `internal_cores = 1L` to avoid nested cluster oversubscription.
##
## Skip-if-unchanged: when `reuse_existing` is TRUE (default driven by the
## `reuse_existing_sims` global), the function builds a signature from the .stan
## source, posterior_draws_<name>.RDA, simulate.R, A.level, and n_draws, and
## compares against tmp/ODE_kelp_<name>.hash. If the signature matches and all
## per-level ODE_kelp_<level>_<name>.RDA outputs exist, the simulation is skipped.

## Muffle the cosmetic "closing unused connection" warnings emitted by R's GC
## when PSOCK cluster sockets get finalized. Legitimate warnings still surface.
.muffle_connection_warnings <- function(expr) {
  withCallingHandlers(expr,
    warning = function(w) {
      if (grepl("closing unused connection", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    }
  )
}

## Build R-callable resourceLoss_<name>(S0, A0, F0, params) from the Stan source
## by extracting the body between the ## ODE_BODY_START ## / END markers.
## Reads the .stan source from `models_dir` and writes the generated .R helper
## to `out_dir`.
make_resourceLoss <- function(model_name, models_dir, out_dir) {
  stan_lines <- readLines(paste0(models_dir, "/stan_model_", model_name, ".stan"))

  start <- which(grepl("## ODE_BODY_START ##", stan_lines)) + 1
  end   <- which(grepl("## ODE_BODY_END ##",   stan_lines)) - 1

  body <- stan_lines[start:end]
  body <- gsub("//", "#",  body)                    # Stan comments -> R comments
  body <- gsub(";(\\s*#|\\s*$)", "\\1", body)       # remove semicolons before # or at end of line

  writeLines(c(
    paste0("resourceLoss_", model_name, " <- function(S0, A0, F0, params) {"),
    "  with(as.list(c(S0, A0, F0)), {",
    body,
    "    return(list(c(dS_dt, dA_dt, dF_dt)))",
    "  })",
    "}"
  ), paste0(out_dir, "/resourceLoss_", model_name, ".R"))
}

simulate_model <- function(model_name, n_draws = NULL, internal_cores = n_cores,
                           reuse_existing = if (exists("reuse_existing_sims", inherits = TRUE))
                                              reuse_existing_sims else FALSE) {

  ## Reuse cached simulation if inputs unchanged and all per-level outputs exist.
  stan_file  <- paste0(models, "/stan_model_",      model_name, ".stan")
  draws_file <- paste0(tmp,    "/posterior_draws_", model_name, ".RDA")
  sim_R_file <- paste0(code,   "/simulate.R")
  hash_file  <- paste0(tmp,    "/ODE_kelp_",        model_name, ".hash")
  out_files  <- paste0(tmp,    "/ODE_kelp_", names(A.level), "_", model_name, ".RDA")

  current_sig <- c(
    paste0("stan=",       unname(tools::md5sum(stan_file))),
    paste0("draws=",      unname(tools::md5sum(draws_file))),
    paste0("simulate_R=", unname(tools::md5sum(sim_R_file))),
    paste0("A.level=",    paste(names(A.level), A.level, sep = ":", collapse = ",")),
    paste0("n_draws=",    if (is.null(n_draws)) "NULL" else as.character(n_draws))
  )

  if (reuse_existing && all(file.exists(out_files)) && file.exists(hash_file)) {
    saved_sig <- readLines(hash_file, warn = FALSE)
    if (identical(saved_sig, current_sig)) {
      message("[", model_name, "] reusing cached simulation (inputs unchanged).")
      return(invisible(NULL))
    }
  }

  ## load posts_df_raw
  load(paste0(tmp, "/posterior_draws_", model_name, ".RDA"))

  ## optional subsampling of posterior draws
  if (!is.null(n_draws) && n_draws < nrow(posts_df_raw)) {
    keep         <- sample.int(nrow(posts_df_raw), size = n_draws)
    posts_df_raw <- posts_df_raw[keep, , drop = FALSE]
  }

  ## build ODE function from Stan source (generated helper lives in tmp/)
  make_resourceLoss(model_name, models, tmp)
  source(paste0(tmp, "/resourceLoss_", model_name, ".R"))
  resourceLoss <- get(paste0("resourceLoss_", model_name))

  ## cluster setup (skip when running serially, e.g. inside mclapply over models)
  if (internal_cores > 1L) {
    cat("[", model_name, "] initiating parallel simulation on ",
        internal_cores, " cores.\n", sep = "")
    cl <- makeCluster(internal_cores)
    on.exit(.muffle_connection_warnings(stopCluster(cl)), add = TRUE)
    invisible(clusterEvalQ(cl, library(deSolve)))
  } else {
    cat("[", model_name, "] running simulation serially.\n", sep = "")
    cl <- NULL
  }

  for (AL in 1:length(A.level)) { # initial kelp abundance (low and high)

    ## create sequences of initial conditions
    len_init <- 80 # 100
    A0 <- seq(A.level[AL], A.level[AL], length.out = len_init)    # Kelp
    S0 <- seq(  1, 300, length.out = len_init)                    # Drift
    F0 <- seq(  0,   0, length.out = len_init)                    # Stomach Fullness
    U  <- 20                                                      # Urchins

    ## set time points (in hrs) for the three Periods
    P1 <- 44 # 44 hrs
    P2 <- 45 # 89 hrs
    P3 <- 45 # 134 hrs

    ## time sequences to pass to ode()
    t.list_P1 <- seq(1, P1, by = 1)
    t.list_P2 <- seq(1, P2, by = 1)
    t.list_P3 <- seq(1, P3, by = 1)

    ## set initial conditions
    init_P1 = c(S = S0[1],
                A = A0[1],
                F = F0[1])

    ## concatenate into list of lists (mapply() doesn't like 'F', so we rename after)
    inits_P1 <- mapply(c,
                       S = S0, A = A0, ff = F0,
                       SIMPLIFY = FALSE)
    inits_P1 <- lapply(inits_P1, function(x){ names(x) <- sub("ff", "F", names(x)); x})

    ## derive ODE parameter columns automatically from posts_df_raw
    non_ode_parms <- c("sigma")   # Stan params not used in the ODE
    ode_parm_cols <- setdiff(
      names(posts_df_raw)[!startsWith(names(posts_df_raw), ".")],
      non_ode_parms
    )

    ## first set of params (single draw, for test ODE run).
    ## U is appended so the ODE body (which references U) can find it inside
    ## with(as.list(c(S0, A0, F0))) -- the old simulate.R relied on U being a
    ## global, which no longer holds now that the body lives inside a function.
    parm_list <- c(setNames(as.numeric(posts_df_raw[1, ode_parm_cols]), ode_parm_cols),
                   U = U)

    ## run single ODE
    out_P1 <- ode(init_P1,
                  times = t.list_P1,
                  func = resourceLoss,
                  parms = parm_list)

    ## full list of params (all draws, for parallel run)
    full_parm_list <- lapply(seq_len(nrow(posts_df_raw)), function(i)
      c(setNames(as.numeric(posts_df_raw[i, ode_parm_cols]), ode_parm_cols),
        U = U)
    )

    ## t.lists for restocking model
    t.list_P1_restock <- seq(1, P1, by = 1)
    t.list_P2_restock <- seq(P1 + 1, P1 + P2, by = 1)
    t.list_P3_restock <- seq(P1 + P2 + 1, P1 + P2 + P3, by = 1)

    ## flatten parameter x initial-condition grid for better parallel utilisation
    print(paste('Kelp level', AL, 'of', length(A.level)))
    n_parms <- length(full_parm_list)
    n_inits <- length(inits_P1)
    combos  <- expand.grid(init_idx = seq_len(n_inits), parm_idx = seq_len(n_parms))

    if (!is.null(cl)) {
      clusterExport(cl,
        c("full_parm_list", "inits_P1", "combos",
          "t.list_P1_restock", "t.list_P2_restock", "t.list_P3_restock",
          "resourceLoss", "P1", "P2", "U"),
        envir = environment())
    }

    flat_results <- .muffle_connection_warnings(
      pblapply(seq_len(nrow(combos)), function(i) {
        x <- full_parm_list[[combos$parm_idx[i]]]
        y <- inits_P1[[combos$init_idx[i]]]

        p1 <- ode(y,
                  times = t.list_P1_restock,
                  func  = resourceLoss,
                  parms = x)

        p2 <- ode(c(y[1], y[2], p1[P1, 4]),
                  times = t.list_P2_restock,
                  func  = resourceLoss,
                  parms = x)

        p3 <- ode(c(y[1], y[2], p2[P2, 4]),
                  times = t.list_P3_restock,
                  func  = resourceLoss,
                  parms = x)

        rbind(p1, p2, p3)
      }, cl = cl)
    )

    ## re-nest into original outs_parms[[parm_idx]][[init_idx]] structure
    outs_parms <- split(flat_results, combos$parm_idx)

    save(outs_parms, P1, P2, P3, S0, A0, F0,
         file = paste0(tmp, "/ODE_kelp_", names(A.level[AL]), '_', model_name, ".RDA"))

    ## drop heavy iteration-local objects so the pblapply closure on the next
    ## AL iteration does not pick them up via lexical scoping (each worker
    ## dispatch would otherwise serialize the entire previous-iteration result).
    rm(outs_parms, flat_results)

  } # end AL (kelp low or high) loop

  writeLines(current_sig, hash_file)

  invisible(NULL)
}

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
