## Aim: Budburst on Single Tree level
- can we detect the Budburst on single tree level (e.g. '70% of a tree sprouted (ausgetrieben))
- if so, what data can provide that information (SRS; Drone; Sensors)

## Epics

**Epic 1: Collect data**

- [] Story 1: Satellite Remote Sensing (SRS) :alarm_clock: 25.06.2021
  - [x] Landsat (30m), Sentinel-2 (10m)
    - [x] optional: get data from *Nico*
    - <del> [] familiarise with [Force](https://force-eo.readthedocs.io/en/latest/howto/index.html) to download the needed data </del>
      -  <del> [] Tutorial "The Datacube" </del>
      - <del> [] Tutorial "Coregistration" </del>
    - </del> [] optional: use gpkg to get the extent of the Uniwald </del>
  - [] Sentinel-1 [e.g. here](https://scihub.copernicus.eu/)
    - [] optional: get data from *Nico*
  - [x] Planetscope (3m)
    - [x] optional: get data from *Nico*
    - [x] new account on [planet.com](https://www.planet.com/) (students account)
    - <del> [] load data with cloud cover of ~ 80% </del>
    - <del> [] get data for 2021 </del>

- [] Story 2: Drone :alarm_clock: *when processed*
  - [] get data from *Nico*, as soon as its processed
  - [] Dense Point Clouds
  - [] Ortophoto

- [] Story 3: TreeSensors
  - [] TreeSense Magnitude :alarm_clock: 25.06.2021
    - [x] download data manually on [treesense.net](https://login.treesense.net/home/) (downloaded: 21.06.2012, 10:00)
    - [x] tidy data (remove duplicate entries; replace "," with "." and convert temperature and magnitude to numeric)
    - [] find out, what the impedance tells us
        - The electrical impedance in the xylem part represents the tree's heart rate. It indicates the well-being of the tree and provides a large variety of insights. (Source: treesense.net); *very 'informative', thanks*
 - [] Sapflow :alarm_clock: *when available*
 - [] Wood Water Content (VWWC) :alarm_clock: *when available*
 - [] Spectral Sensors :alarm_clock: *when available*
 
- [] Story 4: PhenoCams :alarm_clock: 25.06.2021
  - [x] look for budbursts in the cams images
  - [x] ask *Martin*, which cam was placed in which tree
  - [x] search for sources of error; find potential improvement possibilities
  - [] add error/improvement notes to the PhenoCam repo

**Epic 2: Paperwork**
 
- [] Story 1: Comparison of own (planned) work with paper
  - [] read 'Assessing Forest Phenology: A Multi-Scale Comparison of Near-Surface (UAV, Spectral Reﬂectance Sensor, PhenoCam) and Satellite (MODIS, Sentinel-2) Remote Sensing' (Thapa et al. 2021)
  - [] comparison
   - [] which methods did they use
      - forest
      - PhenoCams for In Situ data
   - [] which methods will I use
      - single tree
      - Sensors + manually taken Budburst times as In Situ data
      - validation via phenological observation
      - <del> PhenoCams </del>

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
