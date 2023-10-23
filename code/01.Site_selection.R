

### Set the working directory to the repository root

# Load required packages
require(terra)
require(tidyverse)

# Set my ggplot theme
MyTheme = function() {
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(colour = 'black', fill = NA, linewidth = 0.5),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size = rel(1)),
    axis.ticks = element_blank()
  )
}

### Step 1: Combine multiple metadata from the Australian Microbiome Initiative.

# I have all currently available ITS data from the Australian Microbiome
# Initiative dataset (23/09/2023), which contains two metadata files. One for
# 'old' data from the Biomes of Australian Soil Environments Initiative dataset,
# and another for 'new' data from the current, Australian Microbiome Initiative.
# Additionally, I have an ITS datasets that is not yet available the Australian
# Microbiome initiative data portal. The I will first combine the metadata into
# a single data frame.

data_aus_microbiome <-
  read.csv("data/AusMicrobiome/package_metadata_bpa_e1f6092f_20230923T1206_amdb-genomics-amplicon.csv") %>%
  # Remove data without longitude or latitude values.
  filter(longitude != 0) %>%
  filter(latitude != 0) %>%
  # Australian Microbiome collects soil from two depths: 0-10cm and 10-30cm.
  # Here I retain samples from the top soil layer, 0-10cm.
  filter(depth == 0) %>%
  # Select variables of interest.
  select(sample_id, reads, flow_id, date = utc_date_sampled, longitude,
         latitude, fire, ammonium = ammonium_nitrogen_wt, clay, conductivity,
         elev, exc_aluminium, exc_calcium, exc_magnesium, exc_potassium,
         exc_sodium, nitrate = nitrate_nitrogen, organic_carbon, ph = ph,
         phosphorus_colwell, potassium_colwell, sand, silt, sulphur, texture,
         vegetation_type)
glimpse(data_aus_microbiome)

data_base <-
  read.csv("data/AusMicrobiome/package_metadata_bpa_e1f6092f_20230923T1206_base-genomics-amplicon.csv") %>%
  filter(longitude != 0) %>%
  filter(latitude != 0) %>%
  filter(depth == 0) %>%
  select(sample_id, reads, flow_id, date = utc_date_sampled, longitude,
         latitude, fire, ammonium = ammonium_nitrogen_wt, clay, conductivity,
         elev, exc_aluminium, exc_calcium, exc_magnesium, exc_potassium,
         exc_sodium, nitrate = nitrate_nitrogen, organic_carbon, ph = ph,
         phosphorus_colwell, potassium_colwell, sand, silt, sulphur, texture,
         vegetation_type)
glimpse(data_base)

data_vic_microbiome <-
  read.csv('data/AusMicrobiome/VicMicrobiome.csv') %>%
  filter(Longitude_decimal_degrees != 0) %>%
  filter(Latitude_decimal_degrees != 0) %>%
  # Filter to retain samples from the 'Upper' soil layer, 0-10 cm.
  filter(grepl('Upper' , Depth)) %>%
  # Number of reads per sample are with these data so I will leave blank for
  # now.
  mutate(reads = '0') %>%
  # 'flow_id' refers to the sequencing run. 'flow_id' is not avaiabe for these
  # data. Assuming all samples were sequenced on the same run, I will assign
  # 'flow_id' to 'VicMicro' for now.
  mutate(flow_id = 'VicMicro') %>%
  # Combine fire data into one column
  unite(fire, c('Fire',  'Fire_intensity_if_known')) %>%
  # Select variables of intrest
  select(sample_id = Sample_ID, reads, flow_id, date =
         UTC_date_sampled_YYYY.MM.DD, longitude = Longitude_decimal_degrees,
         latitude = Latitude_decimal_degrees, fire, ammonium =
         Ammonium_nitrogen_mg_kg, clay = Clay_percent_2_microm, conductivity =
         Conductivity_ds_m, elev = Elevation_m, exc_aluminium =
         Exc_aluminium_meq_100g, exc_calcium = Exc_calcium_meq_100g,
         exc_magnesium = Exc_magnesium_meq_100g, exc_potassium =
         Exc_potassium_meq_100g, exc_sodium = Exc_sodium_meq_100g, nitrate =
         Nitrate_nitrogen_mg_kg, organic_carbon = Organic_carbon_percent, ph =
         pH_solid_CaCl2, phosphorus_colwell = Phosphorus_Colwell_mg_kg,
         potassium_colwell = Potassium_Colwell_mg_kg, sand = Sand_percent,
         silt = Silt_percent_2.20_microm, sulphur = Sulphur_mg_kg, texture =
         Texture, vegetation_type = Vegetation_type_controlled_vocab_4)
glimpse(data_vic_microbiome)

# Join data frames
# There are some samples that have been sequenced on multiple runs. I will
# combine these samples into one row, by multiplying the number of reads across
# the different runs, and then grabbing a single value for environemtal
# variable.
data <- rbind(data_aus_microbiome, data_base, data_vic_microbiome) %>%
  mutate(reads = as.numeric(reads)) %>%
  group_by(sample_id) %>%
  summarise(
    reads = sum(reads),
    flow_id = first(flow_id),
    date = first(date),
    longitude = first(longitude),
    latitude = first(latitude),
    fire = first(fire),
    ammonium = first(ammonium),
    clay = first(clay),
    conductivity = first(conductivity),
    elev = first(elev),
    exc_aluminium = first(exc_aluminium),
    exc_calcium = first(exc_calcium),
    exc_magnesium = first(exc_magnesium),
    exc_potassium = first(exc_potassium),
    exc_sodium = first(exc_sodium),
    nitrate = first(nitrate),
    organic_carbon = first(organic_carbon),
    ph = first(ph),
    phosphorus_colwell = first(phosphorus_colwell),
    potassium_colwell = first(potassium_colwell),
    sand = first(sand),
    silt = first(silt), 
    sulphur = first(sulphur),
    texture = first(texture),
    vegetation_type = first(vegetation_type)
  ) %>%
  ungroup() %>%
  mutate(sample_id = str_replace(sample_id, "102.100.100/", "s")) %>%
  mutate(sample_id = str_replace(sample_id, "102.100.100_", "s"))
# There are 2,253 plots with longitude and latitude values.
glimpse(data)

### Step 2: Organise map layers

# Global map
world_map <- rnaturalearth::ne_countries(scale = 'medium', type = 'map_units',
                                         returnclass = 'sf')
# Map of Australia
aus_map <- world_map[world_map$name == 'Australia',]
# Basic Australia layer
ausplot <-
  ggplot() +
  geom_sf(data = aus_map) +
  xlab('Longitude') +
  ylab('Latitude') +
  scale_x_continuous(limits = c(113, 154)) +
  scale_y_continuous(limits = c(-43.5, -10)) +
  MyTheme()

### Step 3: Extract woody vegetation cover data for each plot.

# To create a 'forest plots' dataset, I will use a woody vegetation cover
# threshold of 30%. Altough not required for this project, I will also create
# a 'woodland plots' dataset using a threshold of 10-30%.

# Download the woody vegetation cover data from https://www.tern.org.au/
URL <- "https://dap.tern.org.au/thredds/fileServer/landscapes/remote_sensing/landsat/persistent_green/v2/mosaics/lztmre_aus_y20002011_dm7a2_d20050630_r500m.tif"
destfile <- "data/woody_vegetation_cover_5s.tif"
download.file(URL, destfile)

# Import the woody vegetation cover data:
# Woody veg cover is standardised to between 100 and 200
# Here I will convert the standardised woody veg values to percent cover
woody_veg_cover_rast <- rast(destfile) %>%
  tidyterra::mutate(woody_vegetation_cover_5s =
                     (woody_vegetation_cover_5s * 0.01 - 1) * 100)
woody_veg_cover_rast

# Project plot coordinates:
# Longitude and latitude as points
longlat <- data %>%
  select(x = longitude, y = latitude) %>%
  as.matrix()
# Longitude and latitude as a spatial vector
xy <- vect(longlat, crs = '+proj=longlat +datum=WGS84')
xy
# Project coordinates according to woody_veg_cover_rast
coords <- project(xy, woody_veg_cover_rast)
coords

# Check the projection
plot(woody_veg_cover_rast)
points(coords)

# Extract woody veg cover values for each plot
woody_veg_values <-
  terra::extract(woody_veg_cover_rast, coords) %>%
  select(woody_veg_cover = woody_vegetation_cover_5s)

# Add woody veg cover values to data
# NA values arise from plots on islands and Antarctica that are not covered by 
# the woody veg cover map.
data <- data %>% cbind(woody_veg_values)
glimpse(data)

### Step 4: Subset forest and woodland plots
unique(data$vegetation_type)
# There are 546 forest plots
data_forest_plots <- data %>%
  filter(woody_veg_cover >= 30) %>%
  filter(vegetation_type != "Grassland" &
         vegetation_type != "Cropland" &
         vegetation_type != "Dune" &
         vegetation_type != "Other (Beach)")
glimpse(data_forest_plots)
write_csv(data_forest_plots, 'data/plots_forest.csv')

# There are 1,863 woodland plots
data_woodland_plots <- data %>%
  filter(woody_veg_cover >= 10 | woody_veg_cover <= 30) %>%
    filter(vegetation_type != "Cropland" &
           vegetation_type != "Dune" &
           vegetation_type != "Other (Beach)")
glimpse(data_woodland_plots)
write_csv(data_woodland_plots, 'data/plots_woodland.csv')

### Step 5: Create a map of forest plots

# Longitude and latitude for forest plots
longlat <- data_forest_plots %>%
  select(x = longitude, y = latitude) %>%
  as.matrix()
# Longitude and latitude to spatial vector
xy <- vect(longlat, crs = '+proj=longlat +datum=WGS84')
xy
# Project according to woody_veg_cover_rast
coords <- project(xy, woody_veg_cover_rast)

# Quick plot of forests sites
plot(woody_veg_cover_rast)
points(coords)
# Quick ggplot sites
ausplot +
  geom_point(data = data_forest_plots, 
             aes(x = longitude, y = latitude),
             shape = 21,
             fill = 'black')
# Snazz up your ggplot
forest_ggplots = ggplot() +
  tidyterra::geom_spatraster(data = woody_veg_cover_rast) +
  colorspace::scale_fill_continuous_sequential(
    na.value = NA,
    palette = 'Greens 2',
    name = 'Woody vegetation\ncover (%)') +
  geom_point(data = data_forest_plots, 
             aes(x = longitude, y = latitude),
             shape = 21,
             fill = 'black',
             alpha = 0.7) +
  xlab('Longitude') +
  ylab(NULL) +
  scale_x_continuous(limits = c(113, 154)) +
  scale_y_continuous(limits = c(-43.5, -10),
                     breaks = c(-40, -30, -20, -10)) +
  theme(axis.text.y=element_blank()) +
  MyTheme()

# Subset plots to Australian Temperate Broadleaf and Mixed Forests ecoregions
# Australian bioregions

# Download the data
URL <- "http://www.environment.gov.au/fed/catalog/search/resource/details.page?uuid=%7B4A2321F0-DD57-454E-BE34-6FD4BDE64703%7D/ibra7_regions.zip"
destfile <- "data/ibra7_regions.zip"
download.file(URL, destfile)
unzip(destfile, exdir = "data")
file.remove(destfile)

bioregions_vect = vect('data/IBRA7_regions/ibra7_regions.shp')
bioregions_vect

temp_broad_leaf_bioregions <- c('NSW South Western Slopes',
                                'South East Coastal Plain', 
                                'South East Corner',
                                'South Eastern Highlands',
                                'Southern Volcanic Plain',
                                'Victorian Midlands',
                                'Sydney Basin',
                                'NSW North Coast',
                                'Nandewar',
                                'New England Tablelands',
                                'South Eastern Queensland',
                                'Ben Lomond',
                                'Furneaux',
                                'King',
                                'Tasmanian Central Highlands',
                                'Tasmanian Northern Midlands',
                                'Tasmanian Northern Slopes',
                                'Tasmanian Southern Ranges',
                                'Tasmanian South East',
                                'Tasmanian West')

temp_broad_leaf_vect <- 
  terra::subset(bioregions_vect,
                bioregions_vect$REG_NAME_7 %in% temp_broad_leaf_bioregions)

# Project point data according to the temp_broad_leaf_vect
# Longitude and latitude for forest plots
longlat <- data_forest_plots %>%
  select(x = longitude, y = latitude) %>%
  as.matrix()
# Longitude and latitude to spatial vector
xy <- vect(longlat, crs = '+proj=longlat +datum=WGS84')
xy
# Project according to temp_broad_leaf_vect
coords <- project(xy, temp_broad_leaf_vect)

# Filter forest plots within the temp_broad_leaf_vect
temperete_forest_plots = terra::extract(temp_broad_leaf_vect, coords) %>%
  select(ecoregion = REG_NAME_7) %>%
  bind_cols(data_forest_plots,.) %>%
  drop_na(ecoregion)
glimpse(temperete_forest_plots)

# Create a temperate forest plot template
temperate_forest_ggplots <-
  ausplot +
  tidyterra::geom_spatvector(data = temp_broad_leaf_vect, 
                             aes(fill = REG_NAME_7),
                             fill = 'dark green', alpha = 0.3, colour = NA) +
  scale_x_continuous(limits = c(min(temperete_forest_plots$longitude) - 1,
                                max(temperete_forest_plots$longitude) + 0.5),
                     breaks = c(145, 150, 155)) +
  scale_y_continuous(limits = c(min(temperete_forest_plots$latitude) - 0.25,
                                max(temperete_forest_plots$latitude) + 3),
                     breaks = c(-25, -30, -35, -40))

### Step 6: Select plots to request from the Ausralian Microbiome Initiative

# Subset plots to those that have nitrogen data.
temperete_forest_plots_withN = temperete_forest_plots %>%
  mutate(mineral_nitrogen = nitrate+ammonium) %>%
  filter(!nitrate == 0) %>%
  filter(!ammonium == 0)
glimpse(temperete_forest_plots_withN)
# There are 296 samples with N data

# Graph the nitrogen distribution across samples
temperete_forest_plots_withN %>% 
  select(sample_id, nitrate, ammonium, mineral_nitrogen) %>%
  pivot_longer(-sample_id) %>%
  # Define plot order: arrange by ascending values and group by name
  arrange(value) %>%
  arrange(name) %>%
  ggplot(aes(x = 1:nrow(.), y = log1p(value))) +
  geom_point(shape = 21) +
  facet_wrap(~name, scales = 'free', ncol = 3) +
  labs(x = 'Sample ordered by incresing value',
       y = 'Logarithum of soil\nmineral nitrogen (µmol/L)') +
  theme(axis.text.x=element_blank()) +
  scale_y_continuous(labels = ~ format(.x, scientific = FALSE),
                     limits = c(0, 7),
                     breaks = c(0, 0.5, 1, 1.5, 2, 2.5, 3,
                                3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7)) +
  theme_bw()

# Remove nitrogen outliers
temperete_forest_plots_withN_noOutliers = temperete_forest_plots_withN %>%
  filter(log1p(ammonium) > 1) %>%
  filter(log1p(ammonium) < 4.5) %>%
  filter(log1p(nitrate) > 0.1) %>%
  filter(log1p(nitrate) < 4.5)
glimpse(temperete_forest_plots_withN_noOutliers)
# There are 214 samples after removing outliers

# Graph the nitrogen distribution across samples after removing outliers
temperete_forest_plots_withN_noOutliers %>%
  select(sample_id, nitrate, ammonium, mineral_nitrogen) %>%
  pivot_longer(-sample_id) %>%
  # Define plot order: arrange by ascending values and group by name
  arrange(value) %>%
  arrange(name) %>%
  ggplot(aes(x = 1:nrow(.), y = log1p(value))) +
  geom_point(shape = 21) +
  facet_wrap(~name, scales = 'free', ncol = 3) +
  labs(x = 'Sample ordered by incresing value',
       y = 'Logarithum of soil\nmineral nitrogen (µmol/L)') +
  theme(axis.text.x=element_blank()) +
  scale_y_continuous(labels = ~ format(.x, scientific = FALSE),
                     limits = c(0, 5.5),
                     breaks = c(0, 0.5, 1, 1.5, 2, 2.5,
                                3,3.5, 4, 4.5, 5, 5.5)) +
  theme_bw()

# On visual inspection in Google Earth, there are some obvious non-forest plots
# that I will remove
filter_not_forest <- c("s12444", "s12493", "s12574", "s12618", "s12620", 
                       "s12622", "s12624", "s14187", "s15955", "s15969",
                       "s19469", "s19471", "s19473", "s19513", "s401191",
                       "s401548", "s401549", "s401550", "s401551", "s401553",
                       "s62122", "s62124", "s62136", "s62138", "s62176",
                       "s7851", "s7859", "s7861", "s9573", "s9576")

# These are the plots that I will request for re-sequencing from the Australian
# Microbiome Initiative: There are 183 samples from forests that have a nice
# nitrogen distribution
my_plots <- temperete_forest_plots_withN_noOutliers %>%
  filter(!sample_id %in% filter_not_forest)
glimpse(my_plots)

# Plot the nitrogen distribuiton
my_plots %>%
  select(sample_id, nitrate, ammonium, mineral_nitrogen) %>%
  pivot_longer(-sample_id) %>%
  # Define plot order: arrange by ascending values and group by name
  arrange(value) %>%
  arrange(name) %>%
  ggplot(aes(x = 1:nrow(.), y = log1p(value))) +
  geom_point(shape = 21) +
  facet_wrap(~name, scales = 'free', ncol = 3) +
  labs(x = 'Sample ordered by incresing value',
       y = 'Logarithum of soil\nmineral nitrogen (µmol/L)') +
  theme(axis.text.x=element_blank()) +
  scale_y_continuous(labels = ~ format(.x, scientific = FALSE),
                     limits = c(0, 5.5),
                     breaks = c(0, 0.5, 1, 1.5, 2, 2.5,
                                3,3.5, 4, 4.5, 5, 5.5)) +
  theme_bw()

# Plot the geographic distribution
temperate_forest_ggplots +
  geom_point(data = my_plots, 
             aes(x = longitude, y = latitude),
             shape = 21,
             fill = 'black',
             alpha = 0.6)

# Save the metadata file of the plots that I have selected
write_csv(my_plots, "data/MyPlots.csv")

# Create a sample list to send to the Australian Microbiome Initiative to
# request samples.
my_plots %>%
  select(sample_id) %>%
  mutate(sample_id = str_replace(sample_id, "s", "102.100.100/")) %>%
  write_csv("data/sample_list.csv")
  
