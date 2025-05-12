rifi_preprocess_no_TI <-
function(inp,
         cores,
         FUN_filter = function(x) {
           FALSE
         },
         bg = 0, # used for determination which model so STD or ABG should be used for fit
         rm_FLT = FALSE,
         thrsh_check = 0, # used for filtering
         dista = 300,
         run_PDD = FALSE,
         pen_PDD = 2,
         pen_out_PDD = 1,
         thrsh_PDD = 0.001,
         pen_TI = 10,
         thrsh_TI = 0.5,
         add = 1000) {
  
  message("running check_input...")
  inp <- check_input(inp = inp, thrsh = thrsh_check)
  inp_save <- inp
  inp <- tryCatch({
    for(i in seq_along(metadata(inp)$replicate)){
      tmp_inp <- inp[,colData(inp)$replicate == i]
      logi <- apply(assay(tmp_inp), 1, FUN = FUN_filter)
      rows <- which(logi)
      inp <- encode_FLT(obj = inp, rows = rows, rep = i)
    }
    inp
  },
  error = function(e) {
    writeLines(
      paste(
        "An unknown error has appeared!\n",
        e,
        "An emergency output was returned!\n The given filtration
            malfunctioned!"
      )
    )
    inp <- inp_save
    return(inp)
  }
  )
  
  message("running make_df...")
  inp <- make_df(
    inp = inp,
    cores = cores,
    bg = bg,
    rm_FLT = rm_FLT
  )
  message("running segment_pos...")
  inp <- segment_pos(inp = inp, dista = dista)
  
  if (run_PDD == TRUE) {
    message("running finding_PDD...")
    tryCatch({
      inp <-
        finding_PDD(
          inp = inp,
          pen = pen_PDD,
          pen_out = pen_out_PDD,
          thrsh = thrsh_PDD,
          cores = cores
        )
    },
    error = function(e) {
      writeLines(
        paste(
          "An unknown error has appeared!\n",
          e,
          "An emergency output was returned!\n Please rerun finding_PDD
              manually!"
        )
      )
    }
    )
  }
  metadata(inp)$sessioninfo<-sessionInfo()
  metadata(inp)$date_of_analysis<-date()
  res <- inp
  res
}
