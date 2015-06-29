### Oligoprofile heatmap

Reproduce and modify the oligo heatmap of the [run.profiling.script](profiling). Plots heatmap of the oligo profiles together with phylotype groups and sample clusters. Reload and plot preprocessed data from the [run.profiling.script](profiling) by using the [read.profiling](reading) function, assuming you have stored the output in directory "datadirectory/". Start by reading the data:


```r
library(microbiome)
library(HITChipDB)

# Define data directory (here: simulated data directory)
data.directory <- system.file("extdata", package = "microbiome")

# Read Oligo level data in original domain
probedata <- read_hitchip(data.directory, "frpa")$probedata
```

```
## Loading pre-calculated RPA preprocessing parameters
```

```r
# Read Oligo-phylogeny mapping table (two methods):
taxonomy <- GetPhylogeny("HITChip", "filtered")
```

Save the image in a file (too large to be opened directly in
R). To prevent ordering of the rows, use hclust.method = NULL in the
function call:


```r
library(microbiome)
ppcm <- 150
png(filename = "oligoprofileClustering.png", 
     width = max(trunc(ppcm*21), 
                 trunc(ppcm*21*ncol(probedata)/70)), 
     height = trunc(ppcm*29.7))
tmp <- PlotPhylochipHeatmap(probedata, 
		taxonomy, level = "L1", metric = "pearson")
dev.off()
```

