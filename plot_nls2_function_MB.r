plot_nls2_function_MB <- function(inp, pdf_name) {
  # Preprocess the input and order it
  inp <- inp_order(inp)
  assay(inp)[decode_FLT(inp)] <- NA  # Remove flagged data points
  
  # Extract metadata
  time <- metadata(inp)$timepoints
  fit_STD <- metadata(inp)$fit_STD
  fit_TI <- metadata(inp)$fit_TI
  
  # Create a PDF to save the plots
  pdf(file = paste0("~/Desktop/MB_fit_data/", pdf_name, "_fit_nls2_MB.pdf"),
      width = 7, height = 7)
  par(mfrow = c(2, 2))  # Layout for multiple plots per page
  
  # Iterate through rows in the input
  for (i in seq_len(nrow(inp))) {
    tmp_inp <- inp[i, ]
    # Skip if all data is NA
    if (all(is.na(assay(tmp_inp)))) {
      next
    }
    row_max <- max(assay(tmp_inp), na.rm = TRUE)  # Maximum intensity
    ID <- rowRanges(tmp_inp)$ID
    seg_ID <- rowRanges(tmp_inp)$seg_ID
    
    # Initialize parameters for the legend
    legend_text <- c()
    legend_colors <- c()
    legend_shapes <- c()
    
    # Plot replicates explicitly
    if (!all(is.na(assay(tmp_inp)))) {
      plot(time, assay(tmp_inp), pch = 16, xlab = "Time [min]",
           ylab = "Intensity [A.U.]",
           main = paste0("ID: ", ID,
                         " position: ", rowRanges(tmp_inp)$position,
                         " ", decode(strand(tmp_inp)),
                         "\n", seg_ID
           ),
           ylim = c(0, row_max * 1.2),
           col = colData(inp)$replicate + 1  # Color by replicate
      )
      
      # Add replicates to the legend
      replicates <- unique(colData(inp)$replicate)
      legend_text <- c(legend_text, paste0("Replicate ", replicates))
      legend_colors <- c(legend_colors, replicates + 1)
      legend_shapes <- c(legend_shapes, rep(16, length(replicates)))
    }
    
    # Add model fits (if available) to the legend and plot
    if (ID %in% fit_STD$ID) {
      f <- fit_STD[fit_STD$ID == ID, ]
      if (!any(is.na(f))) {
        curve((f$k / f$decay * x / x + f$bg) * row_max,
              from = 0,
              to = f$delay,
              type = "l",
              add = TRUE,
              pch = 7
        )
        curve(((f$k / f$decay) * exp(-f$decay * (x - f$delay)) + f$bg) * row_max,
              from = f$delay,
              to = max(time),
              type = "l",
              add = TRUE,
              pch = 7
        )
        legend_text <- c(legend_text,
                         paste0("delay = ", round(f$delay, 1)),
                         paste0("halflife = ", round(log(2) / f$decay, 1)),
                         paste0("background = ", round(f$bg * row_max, 1)),
                         paste0("RSS = ", round(f$RSS, 2)),
                         paste0("MAE = ", round(f$MAE, 2)),
                         paste0("R^2 = ", round(f$R_2, 2))
        )
        legend_colors <- c(legend_colors, NA)
        legend_shapes <- c(legend_shapes, NA)
      }
    }
    
    if (ID %in% fit_TI$ID) {
      f <- fit_TI[fit_TI$ID == ID, ]
      if (!any(is.na(f))) {
        curve(((f$k / f$decay - f$ti / f$decay) * x / x + f$bg) * row_max,
              from = 0,
              to = f$ti_delay,
              type = "l",
              add = TRUE,
              pch = 7
        )
        curve(((f$k / f$decay - f$ti / f$decay *
                  exp(-f$decay * (x - f$ti_delay))) + f$bg) * row_max,
              from = f$ti_delay,
              to = f$ti_delay + f$rest_delay,
              type = "l",
              add = TRUE,
              pch = 7
        )
        curve(((f$k / f$decay - f$ti / f$decay *
                  exp(-f$decay * f$rest_delay)) * exp(-f$decay * (
                    x - (f$ti_delay + f$rest_delay))) + f$bg) * row_max,
              from = f$ti_delay + f$rest_delay,
              to = max(time),
              type = "l",
              add = TRUE,
              pch = 7
        )
        legend_text <- c(legend_text,
                         paste0("ti_delay = ", round(f$ti_delay, 1)),
                         paste0("rest_delay = ", round(f$rest_delay, 1)),
                         paste0("halflife = ", round(log(2) / f$decay, 1)),
                         paste0("ti = ", round(f$ti, 1)),
                         paste0("term_prob = ", round(f$ti / f$k, 1)),
                         paste0("background = ", round(f$bg * row_max, 1))
        )
        legend_colors <- c(legend_colors, NA)
        legend_shapes <- c(legend_shapes, NA)
      }
    }
    
    # Add a legend summarizing the data
    legend(
      "topright",
      legend = legend_text,
      col = legend_colors,
      pch = legend_shapes,
      bty = "n",
      cex = 0.8
    )
  }
  
  # Close the PDF device
  dev.off()
}
