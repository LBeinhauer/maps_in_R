# draw your own maps in R

packages <- c("sf", "ggplot2", "rnaturalearth", "rnaturalearthdata", "osmdata", "dplyr")

# Check and install missing packages
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))




# --- 1. Parameters ------------------------------------------------------
city      <- "kiel"
center_ll <- c(10.1394, 54.3233)          # lon, lat
radius    <- 4500                         # metres

# Metric CRS chosen from the longitude (UTM zone), so other cities work too
crs_m <- 32600 + floor((center_ll[1] + 180) / 6) + 1

center   <- st_sfc(st_point(center_ll), crs = 4326)
center_m <- st_transform(center, crs_m)
box_m    <- st_buffer(center_m, radius, endCapStyle = "SQUARE")
bb_m     <- st_bbox(box_m)
bb_ll    <- st_bbox(st_transform(box_m, 4326))

# --- 2. Cache helper ----------------------------------------------------
cache_dir <- "geodata"
dir.create(cache_dir, showWarnings = FALSE)

# e.g. "kiel_10.1394_54.3233_r4500"
area_tag <- sprintf("%s_%.4f_%.4f_r%d", city, center_ll[1], center_ll[2], radius)

cached <- function(layer, query_spec, fn, force = FALSE) {
  tag  <- substr(rlang::hash(query_spec), 1, 8)        # changes if the query changes
  file <- file.path(cache_dir, sprintf("%s_%s_%s.rds", area_tag, layer, tag))
  if (!force && file.exists(file)) {
    message("Loading from cache: ", basename(file))
    return(readRDS(file))
  }
  message("Downloading: ", layer)
  res <- fn()                                          # an error here means nothing gets saved
  saveRDS(res, file)
  res
}

# --- 3. Download (or load) the OSM data ---------------------------------
get_lines <- function(layer, key, values, force = FALSE) {
  cached(layer, list(key, values, crs_m), force = force, fn = function() {
    opq(bbox = bb_ll, timeout = 120) |>
      add_osm_feature(key, values) |>
      osmdata_sf() |>
      (\(x) x$osm_lines)() |>
      st_transform(crs_m)
  })
}

big    <- get_lines("big",    "highway", c("motorway", "trunk", "primary", "secondary", "motorway_link", "trunk_link"))
medium <- get_lines("medium", "highway", c("tertiary", "tertiary_link", "unclassified", "residential"))
# small <- get_lines("small", "highway", c("living_street", "service", "pedestrian", "footway", "path"))

water_poly <- cached("water", "natural=water,bay;crs", function() {
  w <- opq(bbox = bb_ll, timeout = 120) |>
    add_osm_features(list("natural" = "water", "natural" = "bay")) |>
    osmdata_sf()
  bind_rows(
    w$osm_polygons      |> select(geometry),
    w$osm_multipolygons |> select(geometry)
  ) |> st_transform(crs_m) |> st_make_valid()
})

ferry <- cached("ferry", "route=ferry;crs", function() {
  f <- opq(bbox = bb_ll, timeout = 120) |>
    add_osm_feature("route", "ferry") |>
    osmdata_sf()
  bind_rows(
    f$osm_lines      |> select(geometry),
    f$osm_multilines |> select(geometry)
  ) |> st_transform(crs_m)
})

terminals <- cached("terminals", "amenity=ferry_terminal;crs", function() {
  t <- opq(bbox = bb_ll, timeout = 120) |>
    add_osm_feature("amenity", "ferry_terminal") |>
    osmdata_sf()
  t$osm_points |> select(geometry) |> st_transform(crs_m)
})
# --- 3. Radial fade overlay (white, transparent in the middle) ----------
n  <- 800
g  <- expand.grid(x = seq(-1, 1, length.out = n), y = seq(-1, 1, length.out = n))
d  <- sqrt(g$x^2 + g$y^2)
a <- pmin(1, pmax(0, (d - 0.7) / (0.95 - 0.5)))
fade <- matrix(rgb(1, 1, 1, a), nrow = n)

# 2. Let the overlay extend slightly past the map frame, so no edge pixels are left uncovered
pad <- 0.01 * radius * 2

# --- 4. Plot -------------------------------------------------------------
p <- ggplot() +
  geom_sf(data = water_poly, fill = "grey80", color = NA) +
  # geom_sf(data = small,  color = "grey35", linewidth = 0.10) +
  geom_sf(data = ferry, color = "white", linewidth = 0.4, linetype = "dashed") +
  geom_sf(data = terminals, color = "white", size = 1) +
  geom_sf(data = medium, color = "grey15", linewidth = 0.25) +
  geom_sf(data = big,    color = "black",  linewidth = 0.55) +
  annotation_raster(fade,
                    bb_m["xmin"] - pad, bb_m["xmax"] + pad,
                    bb_m["ymin"] - pad, bb_m["ymax"] + pad,
                    interpolate = TRUE) +
  # annotate("text", x = mean(c(bb_m["xmin"], bb_m["xmax"])), y = bb_m["ymin"] + 0.07 * radius * 2,
  #          label = "KIEL", size = 18, fontface = "bold", family = "serif") +
  # annotate("text", x = mean(c(bb_m["xmin"], bb_m["xmax"])), y = bb_m["ymin"] + 0.035 * radius * 2,
  #          label = "54.32°N  10.14°E", size = 5, color = "grey30") +
  coord_sf(xlim = c(bb_m["xmin"], bb_m["xmax"]),
           ylim = c(bb_m["ymin"], bb_m["ymax"]),
           expand = FALSE, datum = NA) +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = NA))

p

# additional functions

transform_coord_to_point <- function(coords, name){
  st_sfc(st_point(coords), crs = 4326) |>
    st_transform(crs_m) |>
    st_sf(name = name, geometry = _) 
}

curve_between <- function(p1, p2, bend = 0.25, n = 60) {
  d      <- p2 - p1
  len    <- sqrt(sum(d^2))
  normal <- c(-d[2], d[1]) / len                          # perpendicular to the segment
  ctrl   <- (p1 + p2) / 2 + bend * len * normal           # control point, offset sideways
  t      <- seq(0, 1, length.out = n)
  m      <- outer((1 - t)^2, p1) + outer(2 * (1 - t) * t, ctrl) + outer(t^2, p2)
  st_linestring(m)
}

bends <- c(0.25, -0.25, 0.25, -0.25, 0.25, -0.25)


stops <- list(
  Bluecher = c(10.13754302323587, 54.3420781340728),
  Garden   = c(10.146973985434174, 54.330622524572235),
  Zoo      = c(10.144266894677251, 54.32821667644662),
  Momiji   = c(10.140664642108897, 54.32782974201529),
  Verace   = c(10.138220452762393, 54.324122914789704),
  Kino     = c(10.133571854610285, 54.32835394686505)   # lon/lat order fixed
)

pts <- do.call(rbind, lapply(names(stops), function(n) {
  transform_coord_to_point(stops[[n]], n)
}))
xy <- st_coordinates(pts)   

route <- st_sfc(
  lapply(seq_len(nrow(xy) - 1), function(i) {
    curve_between(xy[i, ], xy[i + 1, ], bend = bends[i])
  }),
  crs = crs_m
)

# --- 5. plot stops ----------------------------------

p2 <- p +
  #geom_sf(data = route[1], color = "red", linewidth = 0.6, linetype = "dashed") +
  geom_sf(data = pts[1,],   color = "red", size = 3) +
  coord_sf(xlim = c(bb_m["xmin"], bb_m["xmax"]),
           ylim = c(bb_m["ymin"], bb_m["ymax"]),
           expand = FALSE, datum = NA)
p2

p3 <- p +
  geom_sf(data = route[1,], color = "red", linewidth = 0.6, linetype = "dashed") +
  geom_sf(data = pts[1:2,],   color = "red", size = 3) +
  coord_sf(xlim = c(bb_m["xmin"], bb_m["xmax"]),
           ylim = c(bb_m["ymin"], bb_m["ymax"]),
           expand = FALSE, datum = NA)
p3

p4 <- p +
  geom_sf(data = route[1:2,], color = "red", linewidth = 0.6, linetype = "dashed") +
  geom_sf(data = pts[1:3,],   color = "red", size = 3) +
  coord_sf(xlim = c(bb_m["xmin"], bb_m["xmax"]),
           ylim = c(bb_m["ymin"], bb_m["ymax"]),
           expand = FALSE, datum = NA)
p4

p5 <- p +
  geom_sf(data = route[1:3,], color = "red", linewidth = 0.6, linetype = "dashed") +
  geom_sf(data = pts[1:4,],   color = "red", size = 3) +
  coord_sf(xlim = c(bb_m["xmin"], bb_m["xmax"]),
           ylim = c(bb_m["ymin"], bb_m["ymax"]),
           expand = FALSE, datum = NA)
p5

p6 <- p +
  geom_sf(data = route[1:4,], color = "red", linewidth = 0.6, linetype = "dashed") +
  geom_sf(data = pts[1:5,],   color = "red", size = 3) +
  coord_sf(xlim = c(bb_m["xmin"], bb_m["xmax"]),
           ylim = c(bb_m["ymin"], bb_m["ymax"]),
           expand = FALSE, datum = NA)
p6

p7 <- p +
  geom_sf(data = route[1:5,], color = "red", linewidth = 0.6, linetype = "dashed") +
  geom_sf(data = pts[1:6,],   color = "red", size = 3) +
  coord_sf(xlim = c(bb_m["xmin"], bb_m["xmax"]),
           ylim = c(bb_m["ymin"], bb_m["ymax"]),
           expand = FALSE, datum = NA)
p7

ggsave(sprintf("maps/%s_route1_01_Blücher.png", area_tag), p2, width = 20, height = 20, units = "cm", dpi = 300)
ggsave(sprintf("maps/%s_route1_02_Garden.png", area_tag), p3, width = 20, height = 20, units = "cm", dpi = 300)
ggsave(sprintf("maps/%s_route1_03_Zoo.png", area_tag), p4, width = 20, height = 20, units = "cm", dpi = 300)
ggsave(sprintf("maps/%s_route1_04_Moniji.png", area_tag), p5, width = 20, height = 20, units = "cm", dpi = 300)
ggsave(sprintf("maps/%s_route1_05_Verace.png", area_tag), p6, width = 20, height = 20, units = "cm", dpi = 300)
ggsave(sprintf("maps/%s_route1_06_Kino.png", area_tag), p7, width = 20, height = 20, units = "cm", dpi = 300)

