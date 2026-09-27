# Land cells per zoogeographic realm.

suppressPackageStartupMessages({library(terra); library(sf)})
sf_use_s2(FALSE)

zip_url <- "https://macroecology.ku.dk/resources/wallace/cmec_regions___realms.zip"
work    <- file.path(tempdir(), "wallace")
dir.create(work, showWarnings = FALSE, recursive = TRUE)
zipf    <- file.path(work, "cmec.zip")

if (!file.exists(zipf)) download.file(zip_url, zipf, mode = "wb", quiet = TRUE)
unzip(zipf, exdir = work)

kmz <- list.files(work, pattern = "Wallace\\.kmz$", recursive = TRUE, full.names = TRUE)[1]
unzip(kmz, exdir = file.path(work, "kml"))
kml <- file.path(work, "kml", "doc.kml")

realms <- st_make_valid(st_zm(st_read(kml, "Realms", quiet = TRUE)))

v     <- vect(realms["Name"])
v$id  <- seq_len(nrow(v))
grid  <- rast(xmin = -180, xmax = 180, ymin = -90, ymax = 90,
              resolution = 1/24, crs = "EPSG:4326")
r     <- rasterize(v, grid, field = "id")

f   <- as.data.frame(freq(r))
out <- data.frame(Realm      = as.character(realms$Name)[f$value],
                  Land_cells = f$count)
out <- out[order(-out$Land_cells), ]

write.csv(out, "output/region_comparison/Table_S_Realm_Land_Cells.csv", row.names = FALSE)
print(out, row.names = FALSE)
cat("\nTotal land cells at 2.5 arcmin:", format(sum(out$Land_cells), big.mark = ","), "\n")
