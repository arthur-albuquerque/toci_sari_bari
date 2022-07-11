This repository contains code and data to reproduce analysis in the research article (now accepted in the Clinical Microbiology and Infection journal):

## **Effect of Tocilizumab, Sarilumab, and Baricitinib on Mortality Among Patients Hospitalized for COVID-19 Treated with Corticosteroids: A Systematic Review and Meta-Analysis**

PROSPERO: CRD42022297413; [Preregistered Analysis Plan in OSF](https://osf.io/2kstc/)

Authors: Arthur M. Albuquerque, Igor Eckert BS, Lucas Tramujas MD, Guillaume Butler-Laporte MD, Emily G. McDonald MD MSc, James M. Brophy MD PhD, Todd C. Lee MD MPH

## Project organization

```
├── preregistration  <- files regarding our preregistered analysis plan in OSF
├── analysis         <- code for Bayesian model fitting
├── data             <- extracted data from included studies
├── functions        <- Custom function to automate model diagnostics
└── output           <- all files on model outputs, figures, tables, and Appendix Material
    ├── data            <- Modified data files
    ├── figures         <- Code and .pdf files for main figures
    ├── fits            <- Bayesian (brms) models (.rds files)
    └── tables          <- Code and .html files for tables
```

## Getting started

All analyses were conducted in R (R Environment version 4.1.2). 

1.  To download all files and reproduce our analyses, clone this repository using Git's integration with RStudio. Here is a tutorial article in case you are not familiar with this process:

    *Vuorre M, Curley JP. Curating Research Assets: A Tutorial on the Git Version Control System. Advances in Methods and Practices in Psychological Science 2018;1:219–36.* [https://doi.org/10.1177%2F2515245918754826](https://doi.org/10.1177%2F2515245918754826)

       After cloning this repository, open the `toci_sari_bari.Rproj` file and you will be able to run all files.

2. We used [Stan](http://mc-stan.org) through the R package [brms](https://paul-buerkner.github.io/brms/index.html) in this folder. Instructions to download RStan and brms can be found [here](https://mc-stan.org/users/interfaces/rstan.html) and [here](https://github.com/paul-buerkner/brms#how-do-i-install-brms), respectively.

Lastly, we use also used [CmdStanR](https://mc-stan.org/cmdstanr/), a lightweight interface to Stan for R users, to fit our models. Here is the code to install it:

```
> repos <- c(CRAN = "https://cloud.r-project.org", cmdstanr = "https://mc-stan.org/r-packages/")
> options(repos = repos)
> install.packages("cmdstanr")
```

3. We used the R package [{renv}](https://rstudio.github.io/renv/) to make this R project as reproducible as possible. In summary, {renv} guarantees that all required R packages for this project are downloaded to your computer in their necessary versions. Please check their ["Get Started" vignette](https://rstudio.github.io/renv/articles/renv.html) in case you would like to learn more about it.
