#!/usr/bin/env Rscript

# make_seq_som_hex_logo.R
# Creates a hexagon-style package logo/sticker for SEQ SOM using an existing PNG image.
#
# Usage from terminal:
#   Rscript make_seq_som_hex_logo.R path/to/input_logo.png path/to/seq_som_hex_logo.png
#
# Example:
#   Rscript make_seq_som_hex_logo.R a_clean_vector_style_logo_graphic_on_a_white_back.png seq_som_hex_logo.png

# -----------------------------
# 1. Install/load dependencies
# -----------------------------
packages <- c("hexSticker", "magick", "ggplot2", "sysfonts", "showtext")

missing_packages <- packages[!packages %in% rownames(installed.packages())]
if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(hexSticker)
  library(magick)
  library(ggplot2)
  library(sysfonts)
  library(showtext)
})

showtext_auto()

# -----------------------------
# 2. Input/output paths
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

input_logo <- commandArgs(trailingOnly = TRUE)[1]
output_logo <- commandArgs(trailingOnly = TRUE)[2]

# -----------------------------
# 3. Prepare image
# -----------------------------
# This trims extra whitespace and makes the image easier to fit inside the hex.
prepared_png <- tempfile(fileext = ".png")

img <- image_read(input_logo) |>
  image_trim(fuzz = 8) |>
  image_resize("1200x1200")

image_write(img, prepared_png)

sticker(
  subplot = prepared_png,
  package = "",
  p_size = 10,
  p_color = "#0B2545",
  s_x = 1,
  s_y = 1.0,
  s_width = 0.58,
  s_height = 0.58,
  h_fill = "white",
  h_color = "#0B2545",
  filename = output_logo,
  dpi = 600
)

message("Saved hex logo to: ", normalizePath(output_logo))