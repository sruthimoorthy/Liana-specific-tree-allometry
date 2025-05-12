library(ITSMe)


# 1) Define the folder containing .txt point-cloud files
pc_folder <- ""

# 2) Get a list of all .txt files in that folder
pc_files <- list.files(pc_folder, pattern = "\\.txt$", full.names = TRUE)

# 3) Create a data frame to store results
results <- data.frame(
  filename       = character(),
  dbh            = numeric(),
  residual_dbh   = numeric(),
  fdbh           = numeric(),
  dab            = numeric(),
  residual_dab   = numeric(),
  fdab           = numeric(),
  buttress_h    = numeric(),
  height = numeric(),
  crown_area = numeric(),
  crown_vol = numeric(),
  crown_depth = numeric(),
  stringsAsFactors = FALSE
)

# 4) Loop over each file, compute metrics
for (f in pc_files) {
  print(f)
  # Read the point cloud
  pc_tree <- read_tree_pc(path = f)
  
  # Compute DBH metrics
  out_dbh <- dbh_pc(pc = pc_tree, plot = FALSE)
  # out_dbh has components dbh, R2 (residual), and fdbh
  
  # Compute DAB metrics
  out_dab <- dab_pc(pc = pc_tree, plot = FALSE)
  buttress_h <- out_dab$h
  # out_dab has components dab, R2 (residual), and fdab
  H_out <- tree_height_pc(pc = pc_tree)
  if(buttress_h > 1.5) {
    C_out <- classify_crown_pc(pc = pc_tree, minheight = buttress_h, buttress = TRUE,plot = FALSE)
  }
  else {
    C_out <- classify_crown_pc(pc = pc_tree,plot = FALSE)
    
  }

  crown_pc <- C_out$crownpoints
  crownDepth <- tree_height_pc(pc = crown_pc)
  
  pca_out <- projected_area_pc(pc = crown_pc, plot = FALSE)
  cvol_out <- alpha_volume_pc(pc = crown_pc, plot = FALSE)
  
  # 5) Append one row of results
  row_to_add <- data.frame(
    filename      = basename(f),       # just the file name, not full path
    dbh           = out_dbh$dbh,
    residual_dbh  = out_dbh$R2,
    fdbh          = out_dbh$fdbh,
    dab           = out_dab$dab,
    residual_dab  = out_dab$R2,
    fdab          = out_dab$fdab,
    buttress_h    = out_dab$h,
    height        = H_out$h,
    crown_area    = pca_out,
    crown_vol     = cvol_out,
    crown_depth   = crownDepth$h,
    stringsAsFactors = FALSE
  )
  
  results <- rbind(results, row_to_add)
}

# 6) Write all results to a CSV
write.csv(results, file = "ITSMe_tree_metrics_summary.csv", row.names = FALSE)

cat("Done! Results saved to 'tree_metrics_summary.csv'.\n")
