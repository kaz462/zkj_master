## leaflet colorNumeric
map()	# low resolution map of the world
map(wrap = c(200,360), fill = TRUE, col = 2) # pacific-centered map of the world
map(wrap = c(0, 360, NA), fill = TRUE, col = 2) # idem, without Antarctica
map('usa')	# national boundaries
map('county', 'new jersey')	# county map of New Jersey
map('state', region = c('new york', 'new jersey', 'penn'))	# map of three states
map("state", ".*dakota", myborder = 0)	# map of the dakotas
map.axes()				# show the effect of myborder = 0
if(require(mapproj))
  map('state', proj = 'bonne', param = 45)	# Bonne equal-area projection of states

# names of the San Juan islands in Washington state
map('county', 'washington,san', names = TRUE, plot = FALSE)

# national boundaries in one linetype, states in another
# (figure 5 in the reference)
map("state", interior = FALSE)
map("state", boundary = FALSE, lty = 2, add = TRUE)

# plot the ozone data on a base map
# (figure 4 in the reference)
data(ozone)
map("state", xlim = range(ozone$x), ylim = range(ozone$y))
text(ozone$x, ozone$y, ozone$median)
box()
if(require(mapproj)) {	# mapproj is used for  projection="polyconic"
  # color US county map by 2009 unemployment rate
  # match counties to map using FIPS county codes
  # Based on J's solution to the "Choropleth Challenge"
  # http://blog.revolutionanalytics.com/2009/11/choropleth-challenge-result.html
  
  # load data
  # unemp includes data for some counties not on the "lower 48 states" county
  # map, such as those in Alaska, Hawaii, Puerto Rico, and some tiny Virginia
  #  cities
  data(unemp)
  data(county.fips)

  # define color buckets
  colors = c("#F0F8FF", "#87CEFA", "#00BFFF", "#1E90FF", "#0000FF", "#00008B")
  unemp$colorBuckets <- as.numeric(cut(unemp$unemp, c(0, 2, 4, 6, 8, 10, 100)))
  leg.txt <- c("0", "1-5", "6-10", "11-15", "16-20", ">20")
  
  # align data with map definitions by (partial) matching state,county
  # names, which include multiple polygons for some counties
  cnty.fips <- county.fips$fips[match(map("county", plot=FALSE)$names,
                                      county.fips$polyname)]
  colorsmatched <- unemp$colorBuckets [match(cnty.fips, unemp$fips)]
  
  # draw map
  map("county", col = colors[colorsmatched], fill = TRUE, resolution = 0,
      lty = 0, projection = "polyconic")
  map("state", col = "white", fill = FALSE, add = TRUE, lty = 1, lwd = 0.2,
      projection="polyconic")
  title("unemployment by county, 2009")

  
  
  
  colorsmatched <- unemp$colorBuckets [match(cnty.fips, unemp$fips)]
  # draw map
  map("county", "michigan", col = colors[colorsmatched], fill = T, resolution = 0,
          lty = 1)# 
  map.axes()	
  title("unemployment by county, 2009")
  legend(-89.8,44.5, leg.txt, horiz = F, fill = colors)
  

  
  
  
  
  
  
  
  
  
  
  
  
  
  test_data <- readRDS(sprintf("simulated_a5_alpha1/sim_%d.rds", i))
  ysim <- test_data$y
  y1996 <- matrix(ysim, 83, 14)[,14]

  # define color buckets
  colors = c("#F0F8FF", "#87CEFA", "#00BFFF", "#1E90FF", "#0000FF", "#00008B")
  colorBuckets <- as.numeric(cut(y1996, c(0, 1, 5, 10, 15, 20, 1000)))
  leg.txt <- c("0", "1-5", "6-10", "11-15", "16-20", ">20")
  # draw map
  map("county", "michigan", col = colors[colorBuckets], fill = T, resolution = 0,
      lty = 1)# 
  map.axes()	
  title("unemployment by county, 2009")
  legend(-89.8,44.5, leg.txt, horiz = F, fill = colors)

    
  # define color buckets
  colors = c("#F0F8FF", "#87CEFA", "#00BFFF", "#1E90FF", "#0000FF", "#00008B")
  colorBuckets <- as.numeric(cut(y1996, c(0, 1, 5, 10, 15, 20, 1000)))
  leg.txt <- c("0", "1-5", "6-10", "11-15", "16-20", ">20")
  # draw map
  map("county", "michigan", col = colors[colorBuckets], fill = T, resolution = 0,
      lty = 1)# 
  map.axes()	
  title("unemployment by county, 2009")
  legend(-89.8,44.5, leg.txt, horiz = F, fill = colors)
  
  
  
  
  
    