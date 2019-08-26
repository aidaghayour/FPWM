
# Introduction
This model, helps user to have a better understanding of background scenario in sequence preference of TF, and being more precise, does not suffer from present models' drawbacks.
In detail, we aim to build an environment of multiple sequence logos forming a graph. Paths are labeled with weights implying the probability rate of following sub-sequence.


## Current functionalities

1)Building object of TFregulomeR() data, suitable for plotting using ObjectGenerator()

2)Storing a local .txt file of the form Forked-TRANSFAC using StoreFTRANSFACFile()

3)Constructing a class object for plotting with the help of ReadFTRANSFACFile()

4)Plotting a F-PWM object using FPWMPlotter()

5)Storing a bulk data file with concatenating  multiple sets of FPWMs with the help of StoreMultiTRANSFACFile()

6)Plotting local PDF of FPWMs by reading a bulk data file with the help of PlotMultiFTRANSFACFile()

7)Graphical studying the intersection of TF peaks in cell specific enviroment regarding the intersection score using Barchartandseqlogo()

## Release notes
#### This repository is FPWM Beta release! 


#### For development release, please visit [FPWM](https://github.com/aidaghayour/FPWM)


## Installation

#### Prerequisite packages (Will be installed automatically)

1) Required packages: the packages below are the basic prerequisite packages for *FPWM* functionalities

    - [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) (>= 3.0.0)
    - [ggseqlogo](https://cran.r-project.org/web/packages/ggseqlogo/index.html) (>= 0.1)
    - [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) (>= 2.3)
    - [ggplotify](https://cran.r-project.org/web/packages/ggplotify/index.html) (>= 3.4.0)
    - [grid](https://cran.r-project.org/web/packages/ggplotify/index.html) (>= 3.4.0)
    - [lattice](https://cran.r-project.org/web/packages/lattice/index.html) (>= 2)
    - [stringr](https://cran.r-project.org/web/packages/stringr/index.html) (>= 1.1.7)
    - [gridGraphics](https://cran.r-project.org/web/packages/gridGraphics/index.html) (>= 1.6)
    - [base2grob](https://cran.r-project.org/web/packages/base2grob/index.html) (>= 3.4)
    - [ggseqlogo](https://github.com/omarwagih/ggseqlogo)
    - [cowplot](https://cran.r-project.org/web/packages/cowplot/index.html) (≥ 3.5.0)
    - [readr](https://cran.r-project.org/web/packages/readr/index.html) (≥ 1.3.1)
    - [stringr](https://www.rdocumentation.org/packages/stringr/versions/1.4.0) (≥ 1.4.0)
    - [TFregulomeR](https://github.com/benoukraflab/TFregulomeR)

#### Install

In R console,

```r
# if you have not installed "devtools" package
install.packages("devtools")
devtools::install_github("aidaghayour/FPWM")
```
The step above will automatically install the required packages. Note that it needs R version of >= 3.5.2.



## Documentation
You can check detailed package instructions in [Vignettes](https://aidaghayour.github.io/)



