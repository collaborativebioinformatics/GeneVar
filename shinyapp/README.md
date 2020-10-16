Run the `shiny::runApp()` in this folder to launch the Shiny App locally.

This will require R with the following packages:
- shiny
- dplyr
- DT
- ggplot2
- shinydashboard

A docker with these dependencies will be available soon.

The [testdata](testdata) was used while preparing the real data.

## Live version

The demo app is set up on https://www.shinyapps.io/ at https://jmonlong.shinyapps.io/GeneVar/

To set it up, run the following in this folder:

```r
library(rsconnect)
rsconnect::setAccountInfo(name='jmonlong',
			  token='<TOKEN>',
			  secret='<SECRET>')
deployApp(appName='GeneVar')
```

`<TOKEN>` and `<SECRET>` provided on the shinyapps.io "Tokens" account tab.
