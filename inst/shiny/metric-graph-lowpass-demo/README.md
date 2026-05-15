# Metric Graph Low-Pass Demo

Interactive Shiny test harness for `fit.metric.graph.lowpass()` on 1D
Gaussian-mixture synthetic data over `[0, 1]`.

Run from a package checkout after loading the development package:

```r
pkgload::load_all(".")
shiny::runApp("inst/shiny/metric-graph-lowpass-demo")
```

Or run from an installed package:

```r
shiny::runApp(system.file("shiny/metric-graph-lowpass-demo", package = "gflow"))
```

The app separates data generation from model fitting. Mixture, sampling, and
noise controls update the preview immediately; `fit.metric.graph.lowpass()` is
called only when the `Fit metric low-pass` button is clicked.
