#M_var_dir <- "resampled"
#batch_cal <- "Candidate_models"
#out_dir <- "Candidate_Models"
#reg_mult <- c(seq(0.1, 1, 0.1), 2:5)
#f_clas <- "all"
#args <- NULL
#maxent_path <- "/home/jccooper/Dropbox/maxent"
#wait <- FALSE
#run <- TRUE
#occ.joint = "occ_joint.csv"
#occ.tra = "occ_tra.csv"
#M.var.dir = M_var_dir
#batch = batch_cal
#out.dir = out_dir
#reg.mult = reg_mult
#f.clas = f_clas
#args = args
#maxent.path = maxent_path
#wait = wait
#run = run
#max.memory = 1000

jcc_kuenm_cal <- function (occ.joint, occ.tra, M.var.dir, batch, out.dir, max.memory = 1000, 
          reg.mult, f.clas = "all", args = NULL, maxent.path, wait = TRUE, 
          run = TRUE) 
{
  if (!file.exists(occ.joint)) {
    stop(paste(occ.joint, "does not exist in the working directory, check file name", 
               "\nor extension, example: species_joint.csv"))
  }
  if (!file.exists(occ.tra)) {
    stop(paste(occ.tra, "does not exist in the working directory, check file name", 
               "\nor extension, example: species_train.csv"))
  }
  if (missing(M.var.dir)) {
    stop("Argument M.var.dir is not defined.")
  }
  if (!dir.exists(M.var.dir)) {
    stop(paste(M.var.dir, "does not exist in the working directory, check folder name", 
               "\nor its existence."))
  }
  if (length(list.dirs(M.var.dir, recursive = FALSE)) == 0) {
    stop(paste(M.var.dir, "does not contain any subdirectory with environmental variables,", 
               "\neach set of variables must be in a subdirectory inside", 
               paste(M.var.dir, ".", sep = "")))
  }
  if (missing(reg.mult)) {
    warning(paste("Argument reg.mult is not defined, a default set will be used:", 
                  "\n0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 2.0 3.0 4.0 5.0 6.0 8.0 10.0"))
    reg.mult <- c(seq(0.1, 1, 0.1), seq(2, 6, 1), 8, 10)
  }
  if (class(reg.mult) != "numeric") {
    stop("Argument reg.mult must be numeric.")
  }
  if (missing(batch)) {
    warning(paste("Argument batch is not defined, the default name candidate_models", 
                  "\nwill be used."))
    batch <- "candidate_models"
  }
  if (missing(out.dir)) {
    warning(paste("Argument out.dir is not defined, the default name Candidate_Models", 
                  "\nwill be used."))
    out.dir <- "Candidate_Models"
  }
  if (missing(maxent.path)) {
    stop(paste("Argument maxent.path is not defined, it is necessary for executing", 
               "\nthe Maxent software."))
  }
  if (.Platform$OS.type == "unix") {
    sl <- "/"
    dl <- "/"
  } else {
    sl <- "\\"
    dl <- "\\\\"
  }
  # m <- dir(M.var.dir)
  m <- list.dirs(M_var_dir)[-1]
  ms <- paste(getwd(), m, sep = "/")
  env <- paste("environmentallayers=", ms, sep="")
  oc <- occ.joint
  #samp <- paste("samplesfile=", gsub("/", dl, paste("\"", paste(getwd(), 
  #                                                              oc, sep = sl), "\"", sep = "")), sep = "")
  
  samp <- paste0("samplesfile=",getwd(),sl,oc)
                                                                
  occ <- occ.tra
  #samp1 <- paste("samplesfile=", gsub("/", dl, paste("\"", 
  #                                                   paste(getwd(), occ, sep = sl), "\"", sep = "")), sep = "")
  samp1 <- paste0("samplesfile=",getwd(),sl,occ)
  fea <- feature_classes(f.clas)
  dir.create(out.dir)
  out.dir <- paste0(getwd(),sl,out.dir)
  # fix from Can Elverici
  #ram <- paste("-mx", max.memory, "m", sep = "")
  ram <- paste("-Djava.awt.headless=false -mx", max.memory, "m", sep = "")
  #in.comm <- paste("java", ram, paste("-jar", gsub("/", dl, 
  #                                                 paste("\"", paste(maxent.path, "maxent.jar", sep = sl), 
  #                                                       "\"", sep = ""))), sep = " ")
  in.comm <- paste("java", ram, paste("-jar",
                                      paste(maxent.path,"maxent.jar",sep = sl)), sep = " ")
  a.fea <- "autofeature=false"
  fin.com <- "extrapolate=false doclamp=false replicates=1 replicatetype=Crossvalidate responsecurves=false jackknife=false plots=false pictures=false outputformat=raw warnings=false visible=false redoifexists autorun\n"
  fin.com1 <- "extrapolate=false doclamp=false replicates=1 replicatetype=Crossvalidate responsecurves=false jackknife=false plots=false pictures=false outputformat=logistic warnings=false visible=false redoifexists autorun\n"
  if (.Platform$OS.type == "unix") {
    cat("\nCreating directories and maxent batch file, please wait...\n")
    sink(paste(batch, ".sh", sep = ""))
    cat("#! /bin/csh\n")
  } else {
    pb <- winProgressBar(title = "Progress bar", min = 0, 
                         max = length(reg.mult), width = 300)
    sink(paste(batch, ".bat", sep = ""))
  }
  for (i in 1:length(reg.mult)) {
    Sys.sleep(0.1)
    if (.Platform$OS.type == "unix") {} else {
      setWinProgressBar(pb, i, title = paste(round(i/length(reg.mult) * 
                                                     100, 0), "% finished"))
    }
    for (j in 1:length(fea)) {
      for (k in 1:length(ms)) {
        subfol <- paste("outputdirectory=",out.dir,sl,
                        paste("M",reg.mult[i], "F", names(fea)[j],
                              m[k], "all", sep = "_"),
                        sep = "")
        dir.create(paste(out.dir, sl, paste("M", reg.mult[i], 
                                            "F", names(fea)[j], m[k], "all", sep = "_"), 
                         sep = ""),recursive = T,showWarnings = FALSE)
        reg.m <- paste("betamultiplier=", reg.mult[i], 
                       sep = "")
        cat(paste(in.comm, env[k], samp, subfol, reg.m, 
                  a.fea, fea[j], args, fin.com, sep = " "))
        subfol1 <- paste0("outputdirectory=",out.dir, sl, 
                               paste("M", reg.mult[i], "F", names(fea)[j],
                                     m[k], "cal", sep = "_"), sep = "")
        dir.create(paste(out.dir, sl, paste("M", reg.mult[i], 
                                            "F", names(fea)[j], m[k], "cal", sep = "_"), 
                         sep = ""),showWarnings = FALSE)
        cat(paste(in.comm, env[k], samp1, subfol1, reg.m, 
                  a.fea, fea[j], args, fin.com1, sep = " "))
      }
    }
  }
  sink()
  if (.Platform$OS.type != "unix") {
    suppressMessages(close(pb))
  }
  cat("\nIf asked and run = TRUE, allow running as administrator.")
  if (run == TRUE) {
    if (.Platform$OS.type == "unix") {
      batfile_path <- file.path(getwd(), paste(batch, ".sh", 
                                               sep = ""))
      r_wd <- getwd()
      setwd(maxent.path)
      system(paste("bash", batfile_path), wait = wait)
    } else {
      batfile_path <- file.path(getwd(), paste(batch, ".bat", 
                                               sep = ""))
      r_wd <- getwd()
      setwd(maxent.path)
      system2(batfile_path, wait = wait, invisible = FALSE)
    }
    setwd(r_wd)
  }
  cat("\nProcess finished\n")
  if (.Platform$OS.type == "unix") {
    cat(paste("A maxent shell script for creating", i * j * 
                k, "calibration models has been written", sep = " "))
  }
  else {
    cat(paste("A maxent batch file for creating", i * j * 
                k, "calibration models has been written", sep = " "))
  }
  cat(paste("\nCheck your working directory!!!", getwd(), sep = "    "))
}
