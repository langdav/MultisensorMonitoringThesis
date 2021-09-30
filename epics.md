## Aim: Budburst on Single Tree level
- can we detect the Budburst on single tree level (e.g. '70% of a tree sprouted (ausgetrieben))
- if so, what data can provide that information (SRS; Drone; Sensors)

## Epics

**Epic 1: Collecting data**

- [x] Story 1: Satellite Remote Sensing (SRS)
  - [x] <del>Landsat (30m)</del>, Sentinel-2 (10m)
    - [x] optional: get script from *Nico*
    - <del> [] familiarise with [Force](https://force-eo.readthedocs.io/en/latest/howto/index.html) to download the needed data </del>
      -  <del> [] Tutorial "The Datacube" </del>
      - <del> [] Tutorial "Coregistration" </del>
    - <del> [] optional: use gpkg to get the extent of the Uniwald </del>
  - [x] Sentinel-1 [e.g. here](https://scihub.copernicus.eu/)
    - downloaded via [Earthdata Search](https://search.earthdata.nasa.gov/search), using the 'mof.shp' files
    - [x] clip to extent of MOF
    - [x] preprocessing
    - <del> [] optional: get data from *Nico* </del>
  - [x] Planetscope (3.7m an nadir)
    - [x] optional: get data from *Nico*
    - [x] new account on [planet.com](https://www.planet.com/) (students account)
    - <del> [] load data with cloud cover of ~ 80% </del>
    - <del> [] get data for 2021 </del>

- [] Story 2: Drone :alarm_clock: *when processed*
  - [] get data from *Nico*, as soon as its processed
  - [] Dense Point Clouds
  - [] Orthophoto

- [] Story 3: TreeSensors
  - [x] TreeSense Magnitude :alarm_clock: 25.06.2021
    - [x] download data manually on [treesense.net](https://login.treesense.net/home/) (downloaded: 21.06.2012, 10:00)
    - [x] tidy data (remove duplicate entries; replace "," with "." and convert temperature and magnitude to numeric)
    - [] find out, what the impedance tells us
        - The electrical impedance in the xylem part represents the tree's heart rate. It indicates the well-being of the tree and provides a large variety of insights. (Source: treesense.net); *very 'informative', thanks*
  - [] Sapflow :alarm_clock: *when available*
  - [] Wood Water Content (VWWC) :alarm_clock: *when available*
  - [] Spectral Sensors :alarm_clock: *when available*
 
- [] Story 4: PhenoCams
  - [x] look for budbursts in the cams images
  - [x] ask *Martin*, which cam was placed in which tree
  - [x] search for sources of error; find potential improvement possibilities
  - [] add error/improvement notes to the PhenoCam repo
  
- [] Story 5: Budburst phases

**Epic 2: Paperwork**
 
- [x] Story 1: Comparison of own (planned) work with paper
  - [x] read 'Assessing Forest Phenology: A Multi-Scale Comparison of Near-Surface (UAV, Spectral Reﬂectance Sensor, PhenoCam) and Satellite (MODIS, Sentinel-2) Remote Sensing' (Thapa et al. 2021)
      - comparison
          - which methods did they use
              - Phenology of whole forest
              - UAV (multispectral camera; green, red, red edge, near-infrared; orthoimages)
              - Satellite Remote Sensing (MODIS, Sentinel-2)
              - Spectral Reflectance Sensors
              - PhenoCam for In Situ data; 1 Cam, facing south, mounted on a tower at height 7.25 m
              - Normalized Difference Vegetation Index (NDVI), Green Chromatic Coordinate (GCC), Normalized Difference of Green & Red (VIgreen)
          - which methods will I use
              - Budburst on single tree level
              - UAV (multispectral camera; ...; orthoimages; Denseclouds)
              - LiDAR (if applicable; I'd like to use it)
              - Satellite Remote Sensing (optical: Landsat (30m), Sentinel-2 (10m), Planetscope (3m); radar: Sentinel-1 (20m))
              - PhenoCams (more or less): Budburst times determined manually
              - TreeSense Sensors (Impedance, Sapflow, Wood Water Content)
              - Other sensors (Spectral information, temperature, ...)
              - Normalized Difference Vegetation Index (NDVI), Green Chromatic Coordinate (GCC), Enhanced Vegetation Index (EVI)

**Epic 3: Data Processing**

- [] Story 1: Satellite Remote Sensing (SRS) data
  - [] NDVI
  - [] EVI (Enhanced Vegetation Index)

- [] Story 2: Drone data
  - Dense Point Clouds
    - [] extract single trees
    - [] quantify green pixels on per tree basis
  
- [] Story 3: Tree Senors
  - [] higher sapflow after sprouting? (hypothesis: bigger leaf area after sprouting leads to higher sapflow as respiration area is bigger)
  
- <del> [] Story 4: Processing Sentinel-1 data like it was done in Frison et al., 2018 </del>

- [] Story 5: Processing Sentinel-1 data
  - [x] derive various Sentinel-1 indices
  - [x] create plots for visual comparison
      - comparison plots of all 5 trees per Sentinel-1 indice
      - inclusion of budburst phases
 
## (Potential) Sources of errors

- Satellite Remote Sensing
  - clouds (especially Planetscope)
  - data inavailability
  - resolution
 
- Drone
  - lighting (distorted colours)

## Backlog:

- [] familiarise with [Force](https://force-eo.readthedocs.io/en/latest/howto/index.html) to download the needed data </del>
    - [] Tutorial "The Datacube" </del>
    - [] Tutorial "Coregistration"

## Troubleshooting:
