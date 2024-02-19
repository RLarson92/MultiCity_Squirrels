wide_to_stacked <- function(input_df, surveys_per_bout, species){
  nbouts <- ncol(input_df) / surveys_per_bout
  inds <- split(1:(nbouts*surveys_per_bout), rep(1:nbouts, each=surveys_per_bout))
  split_df <- lapply(1:nbouts, function(i){
    out <- input_df[,inds[[i]]]
    out$Site <- FS_occu$Site
    out$City <- FS_occu$City
    out$Season <- i
    out$Species <- species
    names(out)[1:4] <- paste0("week",1:4)
    out
  })
  stack_df <- do.call("rbind", split_df)
  stack_df
}