Predictive mapping of soil types using legacy soil observations
================
Tomislav Hengl (OpenGeoHub), Robert Minarik (OpenGeoHub)
2023-04-16

[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.7820796.svg)](https://doi.org/10.5281/zenodo.7820796)

## Soil types

[Soil type](https://en.wikipedia.org/wiki/Soil_type) is a result of
technical classification of a soil site. Soil type typically reflects
commons soil properties, vertical soil stratification / soil horizons,
diagnostic features, sometimes also soil-forming factors such as climate
and hydrological conditions. Unlike biology and species taxonomy, soil
types are abstract and often subjective, i.e. are not always easy to
validate using laboratory measurements. In that sense, soil
classification can be best compared to classification of climate
(e.g. [Köppen climate
classification](https://en.wikipedia.org/wiki/K%C3%B6ppen_climate_classification))
and [biomes](https://en.wikipedia.org/wiki/Biome), hence they can be
considered fuzzy geographical features (Lodwick, 2007).

Knowing the soil type of your land can be fairly useful as a soil type
conveys complex multivariate information and are ideal summaries of
soil. Soil types are used as input to crop-yield modeling, land use
management planning and landscape modeling of natural hazards and
similar. Soils can be classified using the classification manuals e.g. 
FAO’s IUSS’s [World Reference Base (WRB)
manual](https://www.fao.org/soils-portal/data-hub/soil-classification/world-reference-base/en/)
and/or [USDA Keys to Soil
Taxonomy](https://www.nrcs.usda.gov/resources/guides-and-instructions/keys-to-soil-taxonomy)
(KST). What somewhat adds to complexity of using soil classification is
that many countries have their own classification systems (in a local
language), but often national soil types can be **correlated to some
international system** e.g. WRB or KST (Krasilnikov, Arnold, Marti, &
Shoba, 2009; Verweij, 2017), so that even though different countries use
somewhat different concepts, many soil classification systems are really
comparable (Michéli, Láng, Owens, McBratney, & Hempel, 2016).
Correlation is, nevertheless, not always ideal and it is good to
consider that one soil type in a local system could be translated to at
least 2 or 3 soil types in the international (target) system.

In this computational notebook we show how to generate a global map of
soil types, how to access the [predictions we
produced](https://doi.org/10.5281/zenodo.7820796) (GB of data), and how
to further learn about your own soil. These predictions will be
continuously updated, so if you see an artifact or an issue, [please
report](https://github.com/OpenGeoHub/SoilTypeMapping/issues).

## Mapping soil types using legacy soil observations

Soil types can be mapped from (legacy) soil observations and
measurements (O&M), especially by using those O&M’s that include full or
partial soil profile descriptions (Hengl et al., 2017; Poggio et al.,
2021). Here we are only interested in the global compilations of soil
observations that come with harmonized *analysis-ready* soil type
designations and which are geolocated. For example the World Soil
Information Service (WoSIS) soil profile database (available via:
<https://www.isric.org/explore/wosis>) contains over 30,000 profiles
with soil classification (Batjes, Ribeiro, & Van Oostrum, 2020) and can
be considered one of the most consistent global soil profile
compilations.

To produce the soil type predictions at 1 km spatial resolution, we have
hence decided to use a combination of the legacy soil profiles,
surface-cover observations and simulated points from the HWSDv2. This
gave us about [70,000 training
points](https://opengeohub.github.io/SoilSamples/wrb-soil-types.html) in
total (see figure below). Although there are still some potential issues
with some countries being over-represented, the final training data set
can be considered spatially complete and representative.

<div class="figure">

<img src="img/preview_luvic.chernozems.png" alt="Training points used for global soil type mapping." width="100%" />
<p class="caption">
Training points used for global soil type mapping.
</p>

</div>

As covariate layers for soil type mapping, we use the common global,
(primarily remote sensing based) global layers including:

- Digital Terrain Model parameters based on the MERIT DTM and
  [MERITHydro](http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_Hydro/),
- [Global Lithological DB GLiM and landform
  units](https://doi.org/10.5281/zenodo.1447198),
- [CHELSA Climate bioclimatic
  layers](https://chelsa-climate.org/bioclim/),
- Long-term MODIS time-series day-time and night-time [Land Surface
  Temperature (LST)](https://doi.org/10.5281/zenodo.1420114),
- Long-term MODIS monthly EVI,
- Long-term [monthly snow probability
  images](https://doi.org/10.5281/zenodo.5774953),
- Long-term NASA’s [NEO Water vapor
  images](https://neo.gsfc.nasa.gov/view.php?datasetId=MODAL2_M_SKY_WV),
- Long-term [Cloud fraction images from
  EarthEnv](https://www.earthenv.org/cloud),

After we have imported and harmonized all training points, we run
spatial overlay using the terra package (Hijmans, 2019), which can also
be run in parallel by using e.g.:

``` r
ov.stat = parallel::mclapply(tif.lst, function(i){
        terra::extract(terra::rast(i), 
        terra::vect(as.matrix(tr.pntsF[,c("longitude", "latitude")]), crs="EPSG:4326"))}, 
        mc.cores=80)
ov.stat = dplyr::bind_cols(lapply(ov.stat, function(i){i[,2]}))
names(ov.stat) = make.names(tools::file_path_sans_ext(basename(tif.lst)))
ov.stat$row.id = tr.pntsF$row.id
wrb.rm = plyr::join_all(list(tr.pntsF, ov.stat))
```

A snapshot of the analysis-ready classification matrix can be loaded
from:

``` r
wrb.rm = readRDS("./data/WRB_global_cm_v20230412.rds")
head(wrb.rm[,1:10])
```

    ## # A tibble: 6 x 10
    ## # Groups:   h_wrb4 [4]
    ##   profile… latitu… longit… h_wrb4 source… row.id monthl… monthl… monthl… monthl…
    ##   <chr>      <dbl>   <dbl> <chr>  <chr>    <int>   <dbl>   <dbl>   <dbl>   <dbl>
    ## 1 47431       6.59    2.15 Eutri… {AF-Af…      1   4742.    269.   5476.    89.2
    ## 2 47478       6.59    2.19 Eutri… {AF-Af…      2   4892     254.   5132.   122. 
    ## 3 47503       6.59    2.23 Eutri… {AF-Af…      3   4477.    264.   4917.    16.9
    ## 4 47596       6.87    2.35 Hapli… {AF-Af…      4   5704.    444.   5332.    64.9
    ## 5 52678       1.07   34.9  Hapli… {AF-Af…      5   3474.    801.   4719.   377. 
    ## 6 52686       0.94   35.0  Geric… {AF-Af…      6   2952.    881.   4728.   475.

This contains the target variable `h_wrb4`, coordinates of the points +
values of some 130 covariate layers listed above. The most frequent WRB
soil types in the world seem to be:

``` r
summary(as.factor(wrb.rm$h_wrb4), maxsum=10)
```

    ##   Lithic Leptosols   Eutric Cambisols   Chromic Luvisols Calcaric Cambisols 
    ##               6128               5412               4126               3981 
    ##    Haplic Luvisols  Eutric Stagnosols    Haplic Acrisols   Haplic Calcisols 
    ##               3217               2879               2511               2487 
    ##     Shifting sands            (Other) 
    ##               2156              56927

We can fit a prediction model using this data by running e.g. the
randomForestSRC package (Ishwaran & Kogalur, 2022):

``` r
vs.wrb = readRDS("./data/topvars_wrb4.rds")
dfs = wrb.rm[,c("h_wrb4", make.names(vs.wrb$topvars))]
m.test = randomForestSRC::rfsrc(h_wrb4 ~ ., data=dfs, mtry=88, importance=TRUE, ntree=85)
```

To generate predictions, we can use the sample data for the 200 x 200 km
tile covering part of Hungary:

``` r
g1km = readRDS("./tiles/T9820/data_T9820_1km.rds")
pred = predict(m.test, newdata=g1km@data[,m.test$xvar.names], na.action="na.impute")
probs = as.data.frame(pred$predicted)
```

The variable importance shows that vegetation cover and climate are, in
general, the main explanatory factors for mapping soil types.
Surprisingly landform type and/or lithological unit do not come in the
top 20 most important covariates.

``` r
library(ggplot2)
feat_imp_df <- m.test$importance %>% 
   data.frame() %>% 
   mutate(feature = row.names(.)) 
# plot dataframe
feat_imp_df$relative_importance = 100*feat_imp_df$all/sum(feat_imp_df$all)
feat_imp_df = feat_imp_df[order(feat_imp_df$relative_importance, decreasing = TRUE),]
feat_imp_df$variable = paste0(1:nrow(feat_imp_df), ". ", sapply(row.names(feat_imp_df), function(i){paste(strsplit(i, "_")[[1]][1:4], collapse="_")}))
ggplot(data = feat_imp_df[1:20,], aes(x = reorder(variable, relative_importance), y = relative_importance)) +
   geom_bar(fill = "steelblue",
            stat = "identity") +
   coord_flip() +
   labs(title = "Variable importance",
        x = NULL,
        y = NULL) +
   theme_bw() + theme(text = element_text(size=15))
```

<div class="figure">

<img src="img/global_wrb_variable_importance.png" alt="Variable importance plot for mapping soil types." width="80%" />
<p class="caption">
Variable importance plot for mapping soil types.
</p>

</div>

For predicting the class probabilities the most important metric is most
likely the
[log-loss](https://stats.stackexchange.com/questions/276067/whats-considered-a-good-log-loss)
i.e. the average log loss for all classes.

## How to use these predictions?

You can access all predictions of soil types produced using the training
data and function listed above directly from
[Zenodo](https://doi.org/10.5281/zenodo.7820796). All TIF files are
provided as [COGs](https://www.cogeo.org/), which means that you can
open them directly in QGIS or similar.

You can also query soil types for any longitude latitude using the
function below (no large data download is needed!):

``` r
extract_xy = function(x, y, mc.cores=10){
  cogs = read.csv("./layers/WRB2020_maps_1km.csv")
  out = parallel::mclapply(paste0("/vsicurl/", cogs$filename.lst), function(i){terra::extract(terra::rast(i), 
        terra::vect(matrix(c(x, y), ncol = 2), crs="EPSG:4326"))}, mc.cores=mc.cores)
  out = dplyr::bind_cols(lapply(out, function(i){i[,2]}))
  names(out) = make.names(cogs$d.lst)
  return(out)
}
```

Test it on a location in Hungary:

``` r
out = extract_xy(x=19.2045, y=46.2251)
```

``` r
out[,which(rank(t(out), ties.method = "random") %in% c(ncol(out)-c(2,1,0)))]
```

    ## # A tibble: 1 x 3
    ##   Haplic.Chernozems_p_ Luvic.Chernozems_p_ Petric.Calcisols_p_
    ##                  <int>               <int>               <int>
    ## 1                   33                  37                   4

This shows that the two most probable soil types at this location are
`Haplic.Chernozems` or `Luvic.Chernozems`. You can read more about these
soils by using the [WRB
documents](https://www.fao.org/soils-portal/data-hub/soil-classification/world-reference-base/en/).

Load this function and test predicting soil types for an arbitrary
longitude and latitude. Find photographs of the soil type you get using
some [visual databases](http://www.photosoil.tsu.ru/en/soils?) and
quickly confirm if the predictions and your field observation match.

*Note*: we are not mapping ALL soil types that have ever been classified
on the field. Original list of soil types have been subset to classes
that appear at least 10 times and at least in 2 countries. If you notice
an error or artifact please report via the [Github
repository](https://github.com/OpenGeoHub/SoilTypeMapping/issues).

*Disclaimer*: Use at own risk. These are initial results with limited
accuracy and possible issues with quality of training points, location
errors and harmonization issues. Update of the predictions takes about
4–5 hrs and will be regularly run provided that new training points are
available.

To cite these maps please use:

    @dataset{hengl_t_2023_7820797,
      author       = {Hengl, T. and Minarik, R.},
      title        = {{Global distribution of predicted soil types at 1 
                       km resolution based on the WRB 2022 classification}},
      year         = 2023,
      publisher    = {OpenGeoHub foundation},
      address      = {Wageningen},
      version      = {v0.1},
      doi          = {10.5281/zenodo.7820797},
      url          = {https://doi.org/10.5281/zenodo.7820797}
    }

## How to contribute?

Help us make better soil maps of the world: contribute soil observations
& measurements. If you are a soil point data producer, [register your
data](https://opengeohub.github.io/SoilSamples/) and then send a note so
we can import and add it to the training points.

Test using these maps locally and send us screenshots of potential
issues you discover.

The following key improvements are planned in the next release:

- [Training
  points](https://opengeohub.github.io/SoilSamples/wrb-soil-types.html)
  will be quality controlled and potential artifacts removed,
- Prediction errors will be provided per pixel for each class
  probability,
- Finer-resolution covariate layers (up to 100 m spatial resolution) can
  be added to help increase prediction accuracy,
- [Ensemble
  ML](https://medium.com/nerd-for-tech/extrapolation-is-tough-for-trees-tree-based-learners-combining-learners-of-different-type-makes-659187a6f58d)
  i.e. using 3–4 learners instead of a single learner could further help
  increase accuracy and reduce extrapolation problems,

## Acknowledgments

**[EarthMonitor.org](https://EarthMonitor.org/)** project has received
funding from the European Union’s Horizon Europe research an innovation
programme under grant agreement
**[No. 101059548](https://cordis.europa.eu/project/id/101059548)**.

**[AI4SoilHealth.eu](https://AI4SoilHealth.eu/)** project has received
funding from the European Union’s Horizon Europe research an innovation
programme under grant agreement
**[No. 101086179](https://cordis.europa.eu/project/id/101086179)**.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
line-spacing="2">

<div id="ref-batjes2020standardised" class="csl-entry">

Batjes, N. H., Ribeiro, E., & Van Oostrum, A. (2020). Standardised soil
profile data to support global mapping and modelling (WoSIS snapshot
2019). *Earth System Science Data*, *12*(1), 299–320.
doi:[10.5194/essd-12-299-2020](https://doi.org/10.5194/essd-12-299-2020)

</div>

<div id="ref-hengl2017soilgrids250m" class="csl-entry">

Hengl, T., Mendes de Jesus, J., Heuvelink, G. B., Ruiperez Gonzalez, M.,
Kilibarda, M., Blagotić, A., et al.others. (2017). SoilGrids250m: Global
gridded soil information based on machine learning. *PLoS One*, *12*(2),
e0169748.
doi:[10.1371/journal.pone.0169748](https://doi.org/10.1371/journal.pone.0169748)

</div>

<div id="ref-hijmans2019spatial" class="csl-entry">

Hijmans, R. J. (2019). *<span class="nocase">Spatial data in R</span>*.
Davis, CA: United States: GFC for the Innovation Lab for Collaborative
Research on Sustainable Intensification. Retrieved from
<https://rspatial.org/>

</div>

<div id="ref-Ishwaran2021" class="csl-entry">

Ishwaran, H., & Kogalur, U. B. (2022). *<span class="nocase">Fast
Unified Random Forests for Survival, Regression, and Classification
(RF-SRC)</span>*. CRAN. Retrieved from
<https://cran.r-project.org/package=randomForestSRC>

</div>

<div id="ref-krasilnikov2009handbook" class="csl-entry">

Krasilnikov, P., Arnold, R., Marti, J. J. I., & Shoba, S. (2009). *A
handbook of soil terminology, correlation and classification*. Taylor &
Francis Group.

</div>

<div id="ref-lodwick2007fuzzy" class="csl-entry">

Lodwick, W. (2007). *Fuzzy surfaces in GIS and geographical analysis:
Theory, analytical methods, algorithms and applications*. CRC Press.

</div>

<div id="ref-micheli2016testing" class="csl-entry">

Michéli, E., Láng, V., Owens, P. R., McBratney, A., & Hempel, J. (2016).
Testing the pedometric evaluation of taxonomic units on soil taxonomy—a
step in advancing towards a universal soil classification system.
*Geoderma*, *264*, 340–349.
doi:[10.1016/j.geoderma.2015.09.008](https://doi.org/10.1016/j.geoderma.2015.09.008)

</div>

<div id="ref-poggio2021soilgrids" class="csl-entry">

Poggio, L., De Sousa, L. M., Batjes, N. H., Heuvelink, G., Kempen, B.,
Ribeiro, E., & Rossiter, D. (2021). <span class="nocase">SoilGrids 2.0:
producing soil information for the globe with quantified spatial
uncertainty</span>. *Soil*, *7*(1), 217–240.
doi:[10.5194/soil-7-217-2021](https://doi.org/10.5194/soil-7-217-2021)

</div>

<div id="ref-Verweij2017" class="csl-entry">

Verweij, S. (2017). *Exploring the use of multiple covariates and
machine learning in disaggregating complex soil maps*. Wageningen:
Wageningen University, Soil Geography; Landscape. Retrieved from
<https://edepot.wur.nl/425324>

</div>

</div>
