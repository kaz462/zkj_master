plot_map <- function(counts){
  # define color buckets
  colors = c("#F0F8FF", "#87CEFA", "#00BFFF", "#1E90FF", "#0000FF", "#00008B")
  colorBuckets <- as.numeric(cut(counts, c(0, 1, 5, 10, 15, 20, 1000)))
  leg.txt <- c("0", "1-5", "6-10", "11-15", "16-20", ">20")
  # draw map
  map("county", "michigan", col = colors[colorBuckets], fill = T, resolution = 0,
      lty = 1)
  map.axes()	
  title("Car crash counts")
  legend(-89.8,44.5, leg.txt, horiz = F, fill = colors)
}


## plot real data in the map
plot_map(counts = yy)

## read from simulated data
test_data <- readRDS(sprintf("simulated_a5_alpha1/sim_%d.rds", i))
ysim <- test_data$y
y1996 <- matrix(ysim, 83, 14)[,14]
## plot simulated data in the map
plot_map(counts = y1996)
