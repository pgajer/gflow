# Collaborator installation

Use the [complete-help installation route](../INSTALL.md):

```sh
git clone https://github.com/pgajer/gflow.git
cd gflow
make install-user
```

This generates help and all four vignettes before installing the source archive.
It requires R, a C++17 compiler, GNU make, and Pandoc. The helper installs missing
core and documentation-building R dependencies in the active libraries.
Published dgraphs 0.2.0 is the tested graph dependency. No sibling checkout or
optional `grip`, `malo`, or viewer installation is needed for the introduction.

In R, open `help("gflow-package", package = "gflow")`, then
`vignette("function-guide", package = "gflow")`. The default build supports
serial operation and uses OpenMP when available. See [INSTALL.md](../INSTALL.md)
for library selection and optional compiler configuration.
