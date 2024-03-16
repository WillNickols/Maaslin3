

# MaAsLin 3 #

MaAsLin 3 is the next generation of MaAsLin (Microbiome Multivariable Association with Linear Models). This repository contains the MaAsLin 3 code an an early shell of the Bioconductor package.

The following packages are required dependencies:
```
optparse
logging
data.table
dplyr
pbapply
lmerTest
parallel
lme4
plyr
TcGSA
```

To load the Maaslin3 function, run the following atfter setting `Maaslin3_path` to be the path to `../Maaslin3/R/`:
```
for (R_file in dir(Maaslin3_path, pattern = "*.R$")) {
  source(file.path(Maaslin3_path, R_file))
}
```