# Load utility functions
source("utils.R")

#### ----------------------------
#### 1. Load and Clean BIN Data
#### ----------------------------

peces_bin <- read.csv("data/Actinopterygii_BINs_country.csv", header = TRUE)

# Replace empty strings with NA
peces_bin[] <- lapply(peces_bin, function(i) {
  i[grepl("^\\s*$", as.character(i))] <- NA
  return(i)
})

peces_bin$basin <- NA  # Initialize basin column

#### ----------------------------
#### 2. Fill Basin Information by Coordinates
#### ----------------------------

# Load coordinates and initialize basin
coord_location <- read.csv("data/amz_basin.txt", sep = "\t")
coord_location$basin <- NA
coord_location <- unique(coord_location)

# Define basin equivalences
equivalence <- data.frame(
  "gis_code" = c("6040280410", "6040280400", "6040285990", "6040262110",
                 "6040262220", "6040239720", "6040239820", "6040007950",
                 "6040344660", "6040344540", "6040345180", "6040345000",
                 "6040476790", "6040476800", "6040527670", "6040527680"),
  "basin" = c("Upper Amazon", "Negro", "Madeira", "Trombetas", "Tapajos",
              "Paru/Jari", "Xingu", rep("Tocantins", 9))
)

# Assign basins to coordinates
for (i in 1:nrow(equivalence)) {
  rows <- coord_location$name %in% equivalence$gis_code[i]
  coord_location$basin[rows] <- equivalence$basin[i]
}

# Remove rows without basin information
coord_location <- coord_location[!is.na(coord_location$basin), ]

# Match basin information to BIN data
for (i in 1:nrow(coord_location)) {
  pos <- which(coord_location$lat[i] == peces_bin$lat & coord_location$lon[i] == peces_bin$lon)
  peces_bin$basin[pos] <- coord_location$basin[i]
}

#### ----------------------------
#### 3. Filter for Manual Search from GenBank Data
#### ----------------------------

# Conditions: No coordinates and not marine species

# Load GenBank data
genbank_data <- read.csv("data/genbank_data.csv")
genbank_data$basin <- NA

# Filter GenBank data by removing entries with coordinates
lat_lookup <- setNames(peces_bin$lat, peces_bin$genbank_accession)
coord <- lat_lookup[genbank_data$gb_code]
keep <- which(is.na(coord))
genbank_data <- genbank_data[keep, ]

# Load habitats and filter out marine samples
habitats <- read.csv("data/habitats.csv", header = TRUE)
marin_sp <- habitats$species[habitats$habitat == "Marine"]

sp_lookup <- setNames(peces_bin$species_name, peces_bin$genbank_accession)
sps <- sp_lookup[genbank_data$gb_code]
sps <- gsub(" cf. ", " ", sps)
sps <- sub("^([A-Za-z]+ [A-Za-z]+).*", "\\1", sps)

genbank_data <- genbank_data[!sps %in% marin_sp, ]

# Display papers to manually review for basins
length(unique(genbank_data$paper))
head(unique(genbank_data$paper))

# From now on, these papers will be manually reviewed to determine the
# sample sub-basin.

# After filling in the file (which can be read with the following codes)
# file <- "data/gb_papers_filled.csv"
# genbank_data <- read.csv(file)
# the main database will be filled

# Add manually determined basins to BIN data
basins <- c("Upper Amazon", "Negro", "Madeira", "Trombetas", "Paru/Jari", "Tapajos", "Xingu", "Tocantins")
selected_codes <- genbank_data$gb_code[genbank_data$basin %in% basins]

for (code in selected_codes) {
  basin <- genbank_data$basin[genbank_data$gb_code == code]
  peces_bin$basin[peces_bin$genbank_accession %in% code] <- basin
}

#### ----------------------------
#### 4. Retain Only Amazon Basin Samples
#### ----------------------------

peces_bin <- peces_bin[peces_bin$basin %in% unique(basins), ]
write.csv(peces_bin, "amz_BOLD_bins.csv", row.names = FALSE)

#### ----------------------------
#### 5. Create Sequence File for Alignment
#### ----------------------------

sequences <- as.list(peces_bin$nucleotides)
names(sequences) <- peces_bin$sequenceID

write.fasta(sequences, "amz_BOLD_seqs.fasta")
