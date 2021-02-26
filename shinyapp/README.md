Run the `shiny::runApp()` in this folder to launch the Shiny App locally.

This will require R with the following packages:
- shiny
- dplyr
- DT
- ggplot2
- shinydashboard

A docker with these dependencies is available at `quay.io/jmonlong/genevar` (more details below).

The [testdata](testdata) was used while preparing the real data.

## Live version on shinyapps.io

The demo app is set up on https://www.shinyapps.io/ at https://jmonlong.shinyapps.io/GeneVar/

To set it up, run the following in this folder:

```r
rsconnect::setAccountInfo(name='jmonlong',
			  token='<TOKEN>',
			  secret='<SECRET>')
## for the bioconductor repos
library(BiocManager)
options(repos = BiocManager::repositories())
# setReposirories() # for more recent R versions
## deploy on shinyapps
library(rsconnect)
deployApp(appName='GeneVar')
```

`<TOKEN>` and `<SECRET>` provided on the shinyapps.io "Tokens" account tab.

## Running it locally with docker

[`run-app.R`](run-app.R) launches the app on the 3457 port, so the docker command makes sure to link this port.
The app is then available by opening a web browser at [http://0.0.0.0:3457](http://0.0.0.0:3457)

```
docker run -it -v `pwd`:/app -w /app -p 3457:3457 jmonlong-r-genevar Rscript run-app.R
```

The docker container was created using the [Dockerfile](Dockerfile).
It was pushed to `quay.io/jmonlong/genevar`

This docker container can also be used to deploy the app on shinyapps.io
