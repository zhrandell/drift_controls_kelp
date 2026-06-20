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
## When invoked from inside an outer future_lapply (parallel_models = TRUE in
## RunMe.R), pass `internal_cores = 1L` to avoid nested cluster oversubscription.
##
## Skip-if-unchanged: when `reuse_existing` is TRUE (default driven by the
## `reuse_existing_sims` global), the function builds a signature from the .stan
## source, posterior_draws_<name>.RDA, simulate.R, A.level, and n_draws, and
## compares against tmp/ODE_kelp_<name>.hash. If the signature matches and all
## per-level ODE_kelp_<level>_<name>.RDA outputs exist, the simulation is skipped.

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

  ## Sever the closure's link to simulate_model's frame. Without this,
  ## future_lapply's auto-globals would drag every local in scope (including
  ## posts_df_raw) across to each worker via resourceLoss's enclosing env.
  environment(resourceLoss) <- globalenv()

  ## Switch the future plan only for the duration of this function. When
  ## simulate_model is called from an outer future_lapply (parallel_models =
  ## TRUE in RunMe.R), RunMe.R passes internal_cores = 1L so we don't swap
  ## the plan -- the inner future_lapply runs sequentially under future's
  ## nested-future-safety default, avoiding cluster-of-clusters explosions.
  ##
  ## Cluster setup can fail with "sh: fork: Resource temporarily unavailable"
  ## when internal_cores exceeds the per-user process limit (ulimit -u on
  ## macOS/Linux). When that happens, halve the worker count and retry, and
  ## fall back to sequential execution if even a single worker can't be
  ## spawned. Use a short connectTimeout so the retries don't each wait the
  ## parallelly default of 125 s.
  if (internal_cores > 1L) {
    orig_plan <- plan()
    on.exit(plan(orig_plan), add = TRUE)
    old_timeout <- getOption("parallelly.makeNodePSOCK.connectTimeout", 125)
    options(parallelly.makeNodePSOCK.connectTimeout = 20)
    on.exit(options(parallelly.makeNodePSOCK.connectTimeout = old_timeout),
            add = TRUE)

    workers_try <- internal_cores
    cluster_up  <- FALSE
    while (workers_try > 1L) {
      cat("[", model_name, "] initiating parallel simulation on ",
          workers_try, " cores.\n", sep = "")
      cluster_up <- tryCatch({
        plan(multisession, workers = workers_try)
        TRUE
      }, error = function(e) {
        message("[", model_name, "] cluster setup failed with ",
                workers_try, " workers: ", conditionMessage(e))
        plan(sequential)   # clean up any half-spawned workers
        FALSE
      })
      if (cluster_up) break
      workers_try <- workers_try %/% 2L
    }

    if (!cluster_up) {
      cat("[", model_name, "] running simulation serially ",
          "after cluster setup failures.\n", sep = "")
    }
  } else {
    cat("[", model_name, "] running simulation serially.\n", sep = "")
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
    ## S_present / A_present are the presence flags read by the ODE body's
    ## preference-branching block; simulations always have both resources
    ## present, so both flags are 1 (matching the Low/High treatment branch).
    parm_list <- c(setNames(as.numeric(posts_df_raw[1, ode_parm_cols]), ode_parm_cols),
                   U = U, S_present = 1, A_present = 1)

    ## run single ODE
    out_P1 <- ode(init_P1,
                  times = t.list_P1,
                  func = resourceLoss,
                  parms = parm_list)

    ## full list of params (all draws, for parallel run)
    full_parm_list <- lapply(seq_len(nrow(posts_df_raw)), function(i)
      c(setNames(as.numeric(posts_df_raw[i, ode_parm_cols]), ode_parm_cols),
        U = U, S_present = 1, A_present = 1)
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

    ## future_lapply auto-exports the closed-over globals (full_parm_list,
    ## inits_P1, combos, t.list_P1_restock/P2/P3, resourceLoss, P1, P2) via
    ## future.globals = TRUE. We tick the progressor at most n_ticks times
    ## across the whole combo loop instead of once per ode trio -- otherwise
    ## ~16k progression signals per kelp level saturate the IPC channel from
    ## the multisession workers back to the master and freeze the R event loop.
    ## The inner bar is visible when this stage runs under plan(sequential);
    ## under outer-parallel mode it executes inside a worker and is not
    ## relayed to the master console.
    n_ticks   <- 100L
    tick_step <- max(1L, floor(nrow(combos) / n_ticks))
    flat_results <- with_progress({
      p <- progressor(steps = n_ticks)
      future_lapply(seq_len(nrow(combos)), function(i) {
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

        if (i %% tick_step == 0L) p()
        rbind(p1, p2, p3)
      }, future.globals = TRUE, future.packages = "deSolve", future.seed = TRUE)
    })

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
