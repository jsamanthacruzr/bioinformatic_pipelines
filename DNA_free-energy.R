# ==============================================================================
# Calculate DNA Duplex Free Energy and Generate Profiles
# Description: Calculates dinucleotide free energy profiles (kJ/mol) across 
#              promoter regions for Upregulated, Downregulated, and Unchanged groups.
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Load Required Libraries
# ------------------------------------------------------------------------------
library(seqinr)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# 2. Environment Setup & Relative Paths
# ------------------------------------------------------------------------------
# Set working directory relative to project root
# setwd("./")

output_dir <- "./results/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 3. Free Energy Tables
# ------------------------------------------------------------------------------

# Free energy values in kcal/mol
free_energy_cal <- data.frame(
  dinucleotide = c("aa", "tt", "at", "ag", "ac", "ca", "ct", "cg", "cc", "gg", "ta", "tg", "tc", "ga", "gt", "gc"),
  energy       = c(-1, -1, -0.88, -1.28, -1.44, -1.45, -1.28, -2.17, -1.42, -1.42, -0.58, -1.45, -1.3, -1.3, -1.44, -2.24)
)

# Free energy values in Joules (kJ/mol)
free_energy_jules <- data.frame(
  dinucleotide = c("aa", "tt", "at", "ag", "ac", "ca", "ct", "cg", "cc", "gg", "ta", "tg", "tc", "ga", "gt", "gc"),
  energy       = c(-4.184, -4.184, -3.68192, -5.35552, -6.02496, -6.0668, -5.35552, -9.07928, -5.94128, -5.94128, -2.42672, -6.0668, -5.4392, -5.4392, -6.02496, -9.37216)
)

# ------------------------------------------------------------------------------
# 4. Load Promoter FASTA Sequences
# ------------------------------------------------------------------------------
promoters_upregulated <- read.fasta("UPREGULATED_1uM_promoter_sequences-70bp-3.fasta", seqtype = "DNA")
promoters_silenced    <- read.fasta("DOWNREGULATED_1uM_promoter_sequences-70bp-3.fasta", seqtype = "DNA")
promoters_unchanged   <- read.fasta("NOCHANGE_1uM_promoter_sequences-70bp-3.fasta", seqtype = "DNA")

# ------------------------------------------------------------------------------
# 5. Functions for Free Energy Calculation
# ------------------------------------------------------------------------------

# Function to calculate free energy for a single sequence
calc_free_energy <- function(sequence, free_energy_table) {
  energies <- sapply(1:(nchar(sequence) - 1), function(i) {
    dinucleotide <- substr(sequence, i, i + 1)
    energy       <- free_energy_table$energy[free_energy_table$dinucleotide == dinucleotide]
    return(ifelse(length(energy) > 0, energy, NA)) 
  })
  return(energies)
}

# Function to process each experimental group to calculate free energy
calc_energy_group <- function(promoters, group_name, free_energy_table) {
  do.call(rbind, lapply(1:length(promoters), function(i) {
    sequence <- paste(promoters[[i]], collapse = "")
    energies <- calc_free_energy(sequence, free_energy_table)
    data.frame(
      Group    = group_name,
      Position = 1:length(energies),
      Energy   = energies
    )
  }))
}

# ------------------------------------------------------------------------------
# 6. Compute Free Energies & Bind Datasets
# ------------------------------------------------------------------------------
energy_upregulated <- calc_energy_group(promoters_upregulated, "Upregulated", free_energy_jules)
energy_silenced    <- calc_energy_group(promoters_silenced, "Downregulated", free_energy_jules)
energy_unchanged   <- calc_energy_group(promoters_unchanged, "Unchanged", free_energy_jules)

promoters_energies_df <- rbind(
  energy_upregulated,
  energy_silenced,
  energy_unchanged
)

# ------------------------------------------------------------------------------
# 7. Coordinate Adjustment & Confidence Intervals (95% CI)
# ------------------------------------------------------------------------------

# Adjust positions so that the range corresponds to position relative to TSS (-30 bp offset)
promoters_energies_df_2 <- promoters_energies_df %>%
  mutate(Position_Adjusted = Position - 30)

# Calculate summary statistics by group and position (Mean, SD, SE, 95% CI)
energy_stats_grouped_2 <- promoters_energies_df_2 %>%
  group_by(Group, Position_Adjusted) %>%
  summarise(
    Mean = mean(Energy, na.rm = TRUE),
    SD   = sd(Energy, na.rm = TRUE),
    n    = n(),
    .groups = "drop"
  ) %>%
  mutate(
    SE       = SD / sqrt(n),  
    CI_Lower = Mean - qt(0.975, df = n - 1) * SE,  
    CI_Upper = Mean + qt(0.975, df = n - 1) * SE  
  )

# ------------------------------------------------------------------------------
# 8. Visualization (ggplot2)
# ------------------------------------------------------------------------------
group_colors <- c(
  "Upregulated"   = "#C89092", 
  "Downregulated" = "#2A5783", 
  "Unchanged"     = "#4D4D4D"
)

plot <- ggplot(energy_stats_grouped_2, aes(x = Position_Adjusted, y = Mean, color = Group, group = Group)) +
  geom_line(linewidth = 1.2) +  
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper, fill = Group), alpha = 0.3, color = NA) +
  labs(
    title = "1h 1 µM TPL (± 95% CI)",
    x     = "Nucleotide position",
    y     = "DNA duplex free energy (kJ/mol)"
  ) +
  theme_minimal(base_size = 15) + 
  theme(
    plot.title       = element_text(hjust = 0.5, size = 18, face = "bold", family = "sans"),  
    axis.title       = element_text(size = 18, family = "sans", face = "bold"), 
    axis.text        = element_text(size = 14, family = "sans"),  
    legend.title     = element_text(size = 18, family = "sans", face = "bold"),  
    legend.text      = element_text(size = 14, family = "sans"),  
    legend.position  = "top", 
    panel.grid.major = element_line(color = "gray90", linewidth = 0.5),  
    panel.grid.minor = element_blank(),  
    plot.background  = element_rect(fill = "white", color = NA) 
  ) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  coord_cartesian(xlim = c(-30, 40))

# ------------------------------------------------------------------------------
# 9. Save Figures
# ------------------------------------------------------------------------------
ggsave(file.path(output_dir, "Promoters_1h1uM-v2.svg"), plot = plot, device = "svg")