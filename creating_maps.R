# draw your own maps in R

packages <- c("sf", "ggplot2", "dplyr", "osmextract")
for (pkg in packages) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
invisible(lapply(packages, library, character.only = TRUE))

# --- 1. Parameters (unchanged) ------------------------------------------
city      <- "kiel"
center_ll <- c(10.170455259879034, 54.36673115437426)
radius    <- 10000
crs_m     <- 32600 + floor((center_ll[1] + 180) / 6) + 1

center   <- st_sfc(st_point(center_ll), crs = 4326)
center_m <- st_transform(center, crs_m)
box_m    <- st_buffer(center_m, radius, endCapStyle = "SQUARE")
bb_m     <- st_bbox(box_m)
bb_ll    <- st_bbox(st_transform(box_m, 4326))
box_ll   <- st_as_sfc(bb_ll)

# --- 2. Region file: downloaded once, then reused ------------------------
# --- 2. Storage: region file + processed layers in geodata/ ---------------
geo_dir <- "geodata"
pbf_dir <- file.path(geo_dir, "pbf")
dir.create(pbf_dir, recursive = TRUE, showWarnings = FALSE)

region   <- "Schleswig-Holstein"
area_tag <- sprintf("%s_%.4f_%.4f_r%d", city, center_ll[1], center_ll[2], radius)

# 1) Make the .pbf persistent: it is downloaded once into geodata/pbf, then reused
oe_download_directory_old <- oe_download_directory()
Sys.setenv(OSMEXT_DOWNLOAD_DIRECTORY = normalizePath(pbf_dir))

# 2) Cache each processed layer; the filename encodes area, layer and query
read_osm <- function(layer_name, layer, where, extra_tags = NULL, force = FALSE) {
  tag  <- substr(rlang::hash(list(region, layer, where, extra_tags, crs_m, as.numeric(bb_ll))), 1, 8)
  file <- file.path(geo_dir, sprintf("%s_%s_%s.rds", area_tag, layer_name, tag))
  
  if (!force && file.exists(file)) {
    message("Loading from cache: ", basename(file))
    return(readRDS(file))
  }
  
  message("Reading from OSM extract: ", layer_name)
  res <- oe_get(
    region,
    layer              = layer,
    download_directory = pbf_dir,
    query              = sprintf("SELECT osm_id, geometry FROM %s WHERE %s", layer, where),
    extra_tags         = extra_tags,
    boundary           = box_ll,
    boundary_type      = "clipsrc",
    quiet              = TRUE
  ) |> st_transform(crs_m)
  
  saveRDS(res, file)
  res
}

big <- read_osm("big", "lines",
                "highway IN ('motorway','trunk','primary','secondary','motorway_link','trunk_link')")

medium <- read_osm("medium", "lines",
                   "highway IN ('tertiary','tertiary_link','unclassified','residential')")

water_poly <- read_osm("water", "multipolygons", "\"natural\" IN ('water','bay')") |>
  st_make_valid()

coast     <- read_osm("coast",     "lines",  "\"natural\" = 'coastline'", extra_tags = "natural")
ferry     <- read_osm("ferry",     "lines",  "route = 'ferry'",           extra_tags = "route")
terminals <- read_osm("terminals", "points", "amenity = 'ferry_terminal'", extra_tags = "amenity")

# --- Views ---------------------------------------------------------------
# By hand: center (lon/lat) and half-width in metres
view_center <- function(center_ll, half) {
  c_m <- st_coordinates(st_transform(st_sfc(st_point(center_ll), crs = 4326), crs_m))
  c(xmin = c_m[1] - half, xmax = c_m[1] + half,
    ymin = c_m[2] - half, ymax = c_m[2] + half)
}

# Automatic: square view that keeps every stop inside the un-faded zone
view_fit <- function(pts, min_half = 1500) {
  xy  <- st_coordinates(pts)
  bb  <- st_bbox(pts)
  cx  <- mean(c(bb["xmin"], bb["xmax"]))
  cy  <- mean(c(bb["ymin"], bb["ymax"]))
  rmax <- max(sqrt((xy[, 1] - cx)^2 + (xy[, 2] - cy)^2))
  half <- max(rmax / fade_start, min_half)
  c(xmin = cx - half, xmax = cx + half, ymin = cy - half, ymax = cy + half)
}

# --- Base map for a given view -------------------------------------------
fade_start <- 0.7
fade_end   <- 0.95
n_steps    <- 50

# White overlay: stacked "square minus circle" layers whose alphas add up to a smooth fade
fade_start <- 0.7
fade_end   <- 0.95
n_steps    <- 50

view_box <- function(view) {
  st_bbox(c(xmin = unname(view["xmin"]), xmax = unname(view["xmax"]),
            ymin = unname(view["ymin"]), ymax = unname(view["ymax"])),
          crs = st_crs(crs_m))
}

# White overlay lying exactly inside the view box (square minus growing circle)
fade_mask <- function(view) {
  vb   <- view_box(view)
  sq   <- st_as_sfc(vb)
  cx   <- unname((view["xmin"] + view["xmax"]) / 2)
  cy   <- unname((view["ymin"] + view["ymax"]) / 2)
  half <- unname((view["xmax"] - view["xmin"]) / 2)
  ctr  <- st_sfc(st_point(c(cx, cy)), crs = crs_m)
  
  radii  <- seq(fade_start, fade_end, length.out = n_steps) * half
  target <- seq_len(n_steps) / n_steps
  prev   <- c(0, target[-n_steps])
  step_a <- (target - prev) / (1 - prev)
  
  layers <- lapply(radii, function(r) st_difference(sq, st_buffer(ctr, r, nQuadSegs = 90)))
  st_sf(alpha = step_a, geometry = do.call(c, layers))
}

base_map <- function(view) {
  vb   <- view_box(view)
  half <- unname((view["xmax"] - view["xmin"]) / 2)
  k    <- min(2, max(0.6, sqrt(4500 / half)))
  
  # Clip every layer to the view; skip empty layers safely
  crop <- function(x) {
    if (nrow(x) == 0) return(x)
    suppressWarnings(st_crop(x, vb))
  }
  
  ggplot() +
    geom_sf(data = crop(water_poly), fill = "grey80", color = NA) +
    geom_sf(data = crop(coast),      color = "grey80", linewidth = 0.4 * k) +
    geom_sf(data = crop(ferry),      color = "white",  linewidth = 0.4 * k, linetype = "dashed") +
    geom_sf(data = crop(terminals),  color = "white",  size = 1 * k) +
    geom_sf(data = crop(medium),     color = "grey30", linewidth = 0.25 * k) +
    geom_sf(data = crop(big),        color = "black",  linewidth = 0.55 * k) +
    geom_sf(data = fade_mask(view), aes(alpha = alpha), fill = "white", color = NA) +
    scale_alpha_identity() +
    coord_sf(xlim = unname(c(view["xmin"], view["xmax"])),
             ylim = unname(c(view["ymin"], view["ymax"])),
             expand = FALSE, datum = NA) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))
}

area_tag <- sprintf("%s_%.4f_%.4f_r%d", city, center_ll[1], center_ll[2], radius)

transform_coord_to_point <- function(coords, name) {
  st_sfc(st_point(coords), crs = 4326) |>
    st_transform(crs_m) |>
    st_sf(name = name, geometry = _)
}

transform_stops <- function(stops) {
  do.call(rbind, lapply(names(stops), function(n) transform_coord_to_point(stops[[n]], n)))
}

curve_between <- function(p1, p2, bend = 0.25, n = 60) {
  d      <- p2 - p1
  len    <- sqrt(sum(d^2))
  normal <- c(-d[2], d[1]) / len
  ctrl   <- (p1 + p2) / 2 + bend * len * normal
  t      <- seq(0, 1, length.out = n)
  m      <- outer((1 - t)^2, p1) + outer(2 * (1 - t) * t, ctrl) + outer(t^2, p2)
  st_linestring(m)
}


# --- One call per map: builds and saves every step -----------------------
make_route_map <- function(stops, map_name, view = NULL, bends = c(0.25, -0.25), color = 'red') {
  pts <- do.call(rbind, lapply(names(stops), function(n) transform_coord_to_point(stops[[n]], n)))
  xy  <- st_coordinates(pts)
  if (is.null(view)) view <- view_fit(pts)
  
  n_seg <- nrow(xy) - 1
  p <- base_map(view) + geom_sf(data = pts, color = color, size = 3)
  
  if (n_seg >= 1) {
    bends <- rep(bends, length.out = n_seg)
    route <- st_sfc(lapply(seq_len(n_seg), function(i)
      curve_between(xy[i, ], xy[i + 1, ], bend = bends[i])),
      crs = crs_m)
    p <- p + geom_sf(data = route, color = color, linewidth = 0.6, linetype = "dashed")
  }
  
  dir.create("maps", showWarnings = FALSE)
  ggsave(sprintf("maps/%s_%s.png", area_tag, map_name),
         p, width = 20, height = 20, units = "cm", dpi = 300)
  invisible(p)
}





stops_1 <- list(
  Bluecher = c(10.13754302323587, 54.3420781340728),
  Garden   = c(10.146973985434174, 54.330622524572235),
  Zoo      = c(10.144266894677251, 54.32821667644662),
  Momiji   = c(10.140664642108897, 54.32782974201529),
  Verace   = c(10.138220452762393, 54.324122914789704),
  Kino     = c(10.133571854610285, 54.32835394686505)
)
stops_2 <- list(
  Ferry    = c(10.133925969114692, 54.3152768222387),
  Mediadom = c(10.181361756045352, 54.332533985101136),
  Laboe    = c(10.217934389473765, 54.40438748591363),
  UBoot    = c(10.228734460059481, 54.41242049300569),
  Ocean    = c(10.216033569716602, 54.41196127688331),
  Pincho   = c(10.13711392575389, 54.32221485770658)
)
stops_3 <- list(
  Hammer   = c(10.056817705318851, 54.29652784400112),
  Hasseldi = c(10.082473267332333, 54.32167911059998),
  Suchsdor = c(10.066450104068341, 54.35974601904042),
  Tannenbe = c(10.115726647994343, 54.36466797873711),
  Uhlenkro = c(10.096297610674949, 54.30819410895936),
  Sri_Thai = c(10.13465818279696,  54.35867727745209)
)
stops_4 <- list(
  Trauma  = c(10.116871323058211, 54.334929027559966),
  Mum_Dad = c(10.13196606995149,  54.31915609906696)
)
stops_5 <- list(
  Doener = c(10.128668269951525, 54.319085341223065),
  Kamel  = c(10.13455925461032,  54.32840312887785),
  Bunker = c(10.13213295223185,  54.32191474691611),
  Tuch   = c(10.134149252762485, 54.32784552437226)
)

p1 <- make_route_map(stops_1, "route1", view = view_center(c(10.1394, 54.3233), half = 3500), color = "red")
p2 <- make_route_map(stops_2, "route2", color = 'blue')                       # fitted to the stops
p3 <- make_route_map(stops_3, "route3", color = 'seagreen')
p4 <- make_route_map(stops_4, "route4", view = view_fit(transform_stops(stops_4), min_half = 2500), color = 'sienna')
p5 <- make_route_map(stops_5, "route5", view = view_center(c(10.132, 54.3235), half = 1500), color = 'black', bends = c(0.25, 0.25, 0.25))

p1
p2
p3
p4
p5
