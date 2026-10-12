# FireWater
Code and data to reproduce results in "Ecophysiological adaptations to drought and fire explain tree species occurrences across the western Unites States" by Daniel C. Laughlin, Robert K. Shriver, Karin Kralicek, Andrew Siefert, Sylvain Delzon, William M. Hammond

Abstract

Promoting adequate natural regeneration after fire and preventing the permanent conversion of forest into non-forest habitat is perhaps one of the most pressing challenges for forest managers in western North America. It is well known that climate drives variation in fire regimes, but we still do not understand the relationships between drought and fire adaptive traits in plants. We compiled published trait data for the most common 109 tree species in the western USA to quantify the relationships between fire and drought traits, and to quantify the effects of traits on occurrence probabilities along gradients of water limitation and fire return interval. Resprouters had significantly more vulnerable xylem than non-resprouters. Bark thickness was positively correlated with maximum rooting depth. Species with resistant xylem were more likely to occur in arid climates and in sites with deep water tables, and species with deep roots were more likely to occur in arid climates with a dry season. Species with thick bark were more likely to occur in climates with a dry season and in climates with a frequent fire regime, and resprouting species tended to occur in more frequent fire regimes. Increasing frequencies of drought and fire might select for species with more resistant xylem, deeper roots, thicker bark, and resprouting ability. If true, then the species pool of western USA trees could be missing phenotypes that can maintain forest and woodland structure in a drier, more fire-prone future.


---

## Overview

These data link tree species occurrences on forest inventory plots in the western United States to plot-level climate and environmental conditions, species functional traits, and a species-level phylogeny. The study area spans roughly 31.8°–49.0° N and 124.7°–104.3° W, covering 16 ecological provinces.

At a glance:

| Quantity | Value |
|---|---|
| Plots (presence/absence) | 28,224 |
| Species (presence/absence) | 106 |
| Species (traits, species list, phylogeny) | 109 |
| Ecological provinces | 16 |


## File inventory

| File | Format | Rows × columns | Description |
|---|---|---|---|
| `pres_abs_tree.rds` | R data frame | 963,489 × 4 | Species presence/absence by plot (long format) |
| `plot_climate_data.rds` | R data frame | 28,728 × 52 | Monthly climate normals for each plot |
| `plot_environment_data.rds` | R data frame | 28,138 × 7 | Non-climate environmental variables for each plot |
| `species_trait_data.csv` | CSV | 109 × 19 | Functional traits for each species |
| `phylo_species_list.csv` | CSV | 109 × 3 | Species taxonomy (genus, family) |
| `phylogeny.tre` | Newick | 109 tips | Time-calibrated phylogeny of the study species |
| `provinces.csv` | CSV | 16 × 3 | Lookup table of ecological province codes |

`.rds` files can be read in R with `readRDS("file.rds")`, or in Python with the `pyreadr` package.

## Linking the files

```
provinces.csv ──province──┐
                          ▼
plot_climate_data ──plt_cn──► pres_abs_tree ◄──species── species_trait_data
plot_environment_data ─plt_cn─┘       ▲                  phylo_species_list
                                      └──species──────── phylogeny.tre (tip labels)
```

Plot-level files join on `plt_cn`. Species-level files join on `species` (binomial, space-separated). Tip labels in `phylogeny.tre` use underscores instead of spaces (e.g., `Abies_amabilis`), so replace `_` with a space before joining.

---

## File descriptions

### `pres_abs_tree.rds`

Presence/absence of tree species on each plot, in long format (one row per plot × species combination). Each plot belongs to exactly one province.

The data are **not** a full plots × species matrix: each plot has rows only for the species considered in its province (between 8 and 60 species per province).

| Column | Type | Description |
|---|---|---|
| `province` | character | Ecological province code; see `provinces.csv` |
| `plt_cn` | character | Unique plot identifier |
| `species` | character | Species binomial |
| `pres` | numeric | 1 = present, 0 = absent |

### `plot_climate_data.rds`

Climate normals for each plot.

| Column | Units | Description |
|---|---|---|
| `plt_cn` | — | Unique plot identifier |
| `lon`, `lat` | decimal degrees | Plot coordinates [FIA fuzzed/swapped coordinates] |
| `PPT01`–`PPT12` | mm | Total precipitation, January–December |
| `Tave01`–`Tave12` | °C | Mean temperature, January–December |
| `Tmax01`–`Tmax12` | °C | Mean maximum temperature, January–December |
| `Tmin01`–`Tmin12` | °C | Mean minimum temperature, January–December |
| `pptcv` | unitless | Coefficient of variation of precipitation |


### `plot_environment_data.rds`

Non-climate environmental variables for each plot.

| Column | Units | Description |
|---|---|---|
| `plt_cn` | — | Unique plot identifier |
| `lon`, `lat` | decimal degrees | Plot coordinates |
| `arid` | [TODO] | Aridity index [P/PET, where lower values are drier; observed range 0.03–4.09] |
| `wtd` | [TODO: m?] | Water table depth [Values are ≤ 0, ranging from −704 to 0] |
| `fri` | [TODO: years?] | Fire return interval [observed range 1–213] |
| `extent` |


### `species_trait_data.csv`

Functional traits for each species. Missing values are coded as `NA`.

| Column | Units | Description | n (non-NA) |
|---|---|---|---|
| `species` | — | Species binomial | 109 |
| `p50_mpa` | MPa | Xylem water potential at 50% loss of hydraulic conductivity from Xylem Functional Traits database | 69 |
| `rd_m` | m | Maximum rooting depth from RSIP database | 44 |
| `resprout` | binary | Able to resprout (1) or not (0) from USDA Plants database | 109 |
| `bark_thickness_mm` | mm | Bark thickness at 20 cm Diameter at Breast Height from FVS Regional Variants | 82 |

### `phylo_species_list.csv`

Taxonomy for the species in the phylogeny: 109 species in 35 genera and 15 families. Pinaceae (36 species), Cupressaceae (18) and Fagaceae (18) are the most represented.

| Column | Description |
|---|---|
| `species` | Species binomial |
| `genus` | Genus |
| `family` | Family |

[TODO: taxonomic authority used.]

### `phylogeny.tre`

Rooted, ultrametric phylogeny of the 109 study species in Newick format, with branch lengths. Internal node labels (e.g., `mrcaott248ott650`) are Open Tree of Life identifiers, and the root is labeled `Spermatophyta`.

Read in R with `ape::read.tree("phylogeny.tre")`.

### `provinces.csv`

Lookup table for the ecological province codes used in `pres_abs_tree.rds`.

| Column | Description |
|---|---|
| `province` | Province code (codes starting with `M` are mountain provinces) |
| `province_name` | Full province name |
| `region` | Broader region: `IM` = Intermountain, `PNW` = Pacific Northwest, `PSW` = Pacific Southwest, `RM` = Rocky Mountains  |

---
