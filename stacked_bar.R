
# ===== OrthoFinder shared vs lineage-specific barplot (fixed) =====
library(readr); library(dplyr); library(tidyr); library(stringr); library(ggplot2); library(forcats) ; library(scales)

# ----- Edit paths -----
og_count_file   <- "Orthogroups.GeneCount.tsv"
unassigned_file <- "Orthogroups_UnassignedGenes.tsv"
out_png <- "orthofinder_shared_vs_specific_barplot.png"
dpi_out <- 300


# If NULL, use all species columns in the file.
FOCAL_SPECIES <- NULL

# ---------- 1) Read GeneCount & clean ----------
og <- read_tsv(og_count_file, show_col_types = FALSE)

# Identify non-species columns commonly present
non_species_cols <- c("Orthogroup", "Total", "Average per species", "Average_per_species", "Average")
# Keep species columns only
species_cols_all <- setdiff(names(og), non_species_cols)

# Drop summary rows if present
og <- og %>% filter(!(Orthogroup %in% c("Total","Average per species","Average_per_species","Average")))

# Coerce counts to numeric
og <- og %>% mutate(across(all_of(species_cols_all), ~ suppressWarnings(as.numeric(.))))

# Choose focal species set
species_cols <- if (is.null(FOCAL_SPECIES)) species_cols_all else {
  missing <- setdiff(FOCAL_SPECIES, species_cols_all)
  if (length(missing) > 0) stop("FOCAL_SPECIES not found: ", paste(missing, collapse=", "))
  FOCAL_SPECIES
}

# Build presence matrix for focal species
mat <- as.matrix(og %>% select(all_of(species_cols)))
mat[is.na(mat)] <- 0

n_present <- rowSums(mat > 0)
any_dup   <- apply(mat, 1, function(x) any(x > 1))
all_single_copy <- apply(mat, 1, function(x) all(x == 1))  # strictly 1 in each focal species

# Categories (over focal species):
idx_111            <- which(all_single_copy & (length(species_cols) > 0))
idx_shared_multicp <- which((n_present >= 2) & any_dup)
idx_species_spec   <- which(n_present == 1)

# ---------- 2) Tally counts per species & category ----------
sum_counts <- function(idx, sp) if (length(idx)) sum(og[[sp]][idx], na.rm=TRUE) else 0

# 1:1:1 → one gene per OG per focal species
df_111 <- tibble(
  species  = species_cols,
  count    = length(idx_111),   # exactly one per species per OG
  category = "1:1:1"
)

# N:N:N (shared multicopy) → sum actual counts per species across those OGs
# mat = counts matrix for focal species; og = GeneCount table with those columns
mat <- as.matrix(og[species_cols]); mat[is.na(mat)] <- 0
n_present       <- rowSums(mat > 0)
all_single_copy <- apply(mat, 1, function(x) all(x == 1))

idx_111    <- which(all_single_copy)              # single-copy in ALL focal species
idx_spec   <- which(n_present == 1)               # species-specific
idx_shared <- setdiff(which(n_present >= 2), idx_111)  # shared but NOT 1:1:1  (this is your N:N:N)

sum_counts <- function(idx, sp) if (length(idx)) sum(og[[sp]][idx], na.rm=TRUE) else 0

df_111 <- tibble(species = species_cols, count = length(idx_111), category = "1:1:1")
df_NNN <- tibble(species = species_cols,
                 count   = vapply(species_cols, function(sp) sum_counts(idx_shared, sp), numeric(1)),
                 category= "N:N:N")
# species-specific
spec_counts <- sapply(species_cols, function(sp){
  idx_sp <- which((mat[, sp] > 0) & (n_present == 1))
  sum(og[[sp]][idx_sp], na.rm=TRUE)
})
df_spec <- tibble(species = species_cols, count = as.numeric(spec_counts), category = "species-specific")

# ---------- 3) Unassigned ----------
unassigned <- read_tsv(unassigned_file, show_col_types = FALSE)
# Keep overlapping species only
common_species <- intersect(species_cols, names(unassigned))
count_unassigned <- function(cell) {
  if (is.na(cell) || !nzchar(cell)) return(0L)
  sum(nzchar(str_split(cell, ",\\s*")[[1]]))
}
ua_counts <- sapply(species_cols, function(sp) {
  if (sp %in% common_species) sum(vapply(unassigned[[sp]], count_unassigned, integer(1))) else 0L
})
df_ua <- tibble(species = species_cols, count = as.numeric(ua_counts), category = "unassigned")

# ---------- 4) Combine & plot ----------
df_plot <- bind_rows(df_111, df_NNN, df_spec, df_ua) %>%
  group_by(species, category) %>%
  summarise(count = sum(count, na.rm=TRUE), .groups="drop")
unique(df_plot$species)

# ---------- Species mapping ----------
# Map long species names → cleaner labels
species_map <- c(
  "01_hilli_modified"          = "HIL",
  "02_quadrimaculata_modified" = "QUAD",
  "03_stygia_modified"         = "STY",
  "04_vicina_modified"         = "VIC",
  "05_dros_mel_modified"       = "DROS"
)

# Apply mapping & ordering
df_plot <- df_plot %>%
  mutate(
    species = recode(species, !!!species_map),
    species = factor(species, levels = c("HIL", "QUAD", "STY", "VIC", "DROS"))
  )

# ---------- Category order ----------
level_order <- c("1:1:1", "N:N:N", "species-specific", "unassigned")

# ---------- Colours ----------
cols <- c(
  "1:1:1"            = "#d73027",
  "N:N:N"            = "#4575b4",
  "species-specific" = "#4daf4a",
  "unassigned"       = "#7f7f7f"
)

# ---------- Normal barplot ----------
p <- ggplot(df_plot, aes(x = count, y = species, fill = category)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = cols, limits = level_order, drop = FALSE, name = NULL) +
  scale_x_continuous(
    limits = c(0, 16000),
    breaks = seq(0, 16000, 2000),
    labels = scales::comma,
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(x = "Number of Genes", y = NULL, title = " ") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save outputs
ggsave("orthofinder_shared_vs_specific_barplot.png", p,
       width = 11, height = 6.5, dpi = 600)

ggsave("orthofinder_shared_vs_specific_barplot.pdf", p,
       width = 11, height = 6.5, device = cairo_pdf)

ggsave("orthofinder_shared_vs_specific_barplot.tiff", p,
       width = 11, height = 6.5, dpi = 600, compression = "lzw")

p
names(og)



# ---------- Reversed-axis barplot (mirror version) ----------
p2 <- ggplot(df_plot, aes(x = count, y = species, fill = category)) +
  geom_col(width = 0.8, position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values = cols, limits = level_order, drop = FALSE, name = NULL) +
  scale_x_continuous(
    limits = c(16000, 0),              
    breaks = seq(0, 16000, 2000),
    labels = scales::comma,
    expand = expansion(mult = c(0.02, 0)),
    trans = "reverse"
  ) +
  scale_y_discrete(position = "right") +
  labs(x = "Number of Genes", y = NULL, title = " ") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(size = 13, face = "bold"),
    axis.title.x = element_text(size = 14, face = "bold"),
    legend.position = "right",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Save reversed plots
ggsave("orthofinder_shared_vs_specific_barplot_flipd.png", p2,
       width = 11, height = 6.5, dpi = 600)

ggsave("orthofinder_shared_vs_specific_barplot_flipd.pdf", p2,
       width = 11, height = 6.5, device = cairo_pdf)

p2
