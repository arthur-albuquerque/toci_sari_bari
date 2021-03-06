# Load ----

# Load packages

# Ensures the package "pacman" is installed
if (!require("pacman")) install.packages("pacman")

pacman::p_load(here,
               rio,
               metafor,
               dplyr,
               brms,
               tidybayes)

# Ensure that packages are installed according to versions in renv's lockfile
renv::restore()

# Load data
d_logOR = readRDS(
  here::here("output", "data", "effect_sizes.Rds")
)

# Load meta-analyses
ma_toci = readRDS(here::here("output", "fits", "ma_toci.Rds"))
ma_sari = readRDS(here::here("output", "fits", "ma_sari.Rds"))
ma_bari = readRDS(here::here("output", "fits", "ma_bari.Rds"))
ma_bari_sens = readRDS(here::here("output", "fits", "ma_bari_sens.Rds"))

# Helper functions to add model information to forest plot

tau_text <- function(model) {
  
  posterior = 
    model |> 
    tidybayes::tidy_draws() |> 
    ggdist::median_hdi(sd_study__Intercept)
  
  bquote(paste(tau,
               .(" = "),
               .(formatC(posterior$sd_study__Intercept, digits=2, format="f")),
               .(" [95% CrI: "),
               .(formatC(posterior$.lower, digits=2, format="f")),
               ", ",
               .(formatC(posterior$.upper, digits=2, format="f")),
               "]"
  )
  )
}

PI_text <- function(pred) {
  
  bquote(paste("PI: [",
               .(formatC(exp(pred$ymin), digits=2, format="f")),
               .(", "),
               .(formatC(exp(pred$ymax), digits=2, format="f")),
               "]"
  )
  )
}

# Figure 1 ----

# Effect sizes

d_toci = d_logOR |> 
  dplyr::filter(treatment == "tocilizumab")

d_sari_bari = d_logOR |>
  dplyr::filter(treatment %in% c("sarilumab","baricitinib"),
                study != "RECOVERY Bari") |> 
  arrange(desc(study), treatment)
  

# Bayesian meta-analyses overall results (mu), median and 95% HDI

toci_results = ma_toci |> brms::fixef(summary = F) |> ggdist::median_hdi()
sari_results = ma_sari |> brms::fixef(summary = F) |> ggdist::median_hdi()
bari_results = ma_bari |> brms::fixef(summary = F) |> ggdist::median_hdi()

# Posterior predictive distributions (predictive intervals [PI])

nd = data.frame(study = "new", sei = 0)

predictions = function(model){
  
  set.seed(123)
  
  # Explanations: 
  # https://twitter.com/bmwiernik/status/1473306749906169858
  # https://twitter.com/IsabellaGhement/status/1458539033131302913/photo/1
  # https://www.andrewheiss.com/blog/2021/11/10/ame-bayes-re-guide/
  
  brms::posterior_predict(object = model,
                          newdata = nd,
                          re_formula = NULL,
                          allow_new_levels = TRUE,
                          sample_new_levels = "gaussian") |> 
    data.frame() |>
    ggdist::median_hdi() |> 
    dplyr::rename(y = 1,
                  ymin = .lower,
                  ymax = .upper)
}

toci_pred = predictions(ma_toci)
sari_pred = predictions(ma_sari)
bari_pred = predictions(ma_bari)

# Tocilizumab

forest_fun_toci = function(){
  
  metafor::forest(d_toci$yi, d_toci$vi,
                  cex=0.75, # text size
                  ylim=c(0, 22.5), # Y length
                  rows=c(3:18),
                  alim=log(c(0.25, 4)),
                  xlim=c(-8, 3.5), 
                  at=log(c(0.25, 0.5, 1, 2, 4)),
                  slab=d_toci$study,
                  ilab=cbind(paste(d_toci$trt_events, "/", d_toci$trt_total),
                             paste(d_toci$control_events, "/", d_toci$control_total)),
                  ilab.xpos=c(-4,-2.5),
                  header="Study",
                  atransf=exp,
                  order=d_toci$vi,
                  pch = 19) # Circles 
  

abline(h=2.2) # horizontal line

metafor::addpoly(x = toci_results$y,
                 ci.lb = toci_results$ymin,
                 ci.ub = toci_results$ymax,
                 
                 mlab = "Total [95% CrI]", # label
                 row=1, # location
                 atransf=exp,
                 efac = 1.5 # polygon size
                 )

text(c(-4,-2.5),
     1,
     c(paste(sum(d_toci$trt_events), "/", sum(d_toci$trt_total)),
       paste(sum(d_toci$control_events), "/", sum(d_toci$control_total))),
     font = 2, cex = 0.75)


# Prediction interval
colp = "#6b58a6"
coll = "#a7a9ac"

metafor::addpoly(x = toci_pred$y,
                 ci.lb = toci_pred$ymin,
                 ci.ub = toci_pred$ymax,
                 
                 pi.lb = toci_pred$ymin, # prediction interval
                 pi.ub = toci_pred$ymax, # prediction interval
                 annotate=FALSE,
                 
                 mlab = " ", # label
                 row=-0.25, # location
                 atransf=exp,
                 efac = 0.5, # polygon size
                 col=colp, border=colp
)

# Headers

text(mean(c(-4,-2.5)) - 0.1,
     22.5,
     "No of events / total",
     font=2,
     cex = 0.7)

# Horizontal line
segments(-4 - 0.6,
         22,
         -2.5 + 0.35,
         22)

text(c(-4,-2.5), 21.5, c("Experimental","Control"), font = 2, cex = 0.7)

text(c(-0.6,0.6), # X axis
     19.5, # Y axis
     cex = 0.7, # size
     c("Favors\nExperimental","Favors\nControl"),
     pos=c(2,4), # Right + Left aligned
     offset=-1)

text(-7.1,19.5,
     "Tocilizumab",
     font=4, # bold
     cex = 0.9)

text(-7.95, -0.25,
     tau_text(ma_toci),
     pos=4,
     cex=0.8)

text(1.55, -0.25,
     PI_text(toci_pred),
     pos=4,
     cex = 0.85)
}

# Baricitinib + Sarilumab

forest_fun_bari_sari = function(){

  metafor::forest(d_sari_bari$yi, d_sari_bari$vi,
                  cex=0.75, # text size
                  ylim=c(0, 23.5), # Y length
                  rows=c(3:10, 16:19),
                  alim=log(c(0.25, 4)),
                  xlim=c(-8, 3.5), 
                  at=log(c(0.25, 0.5, 1, 2, 4)),
                  slab=d_sari_bari$study,
                  ilab=cbind(paste(d_sari_bari$trt_events, "/", d_sari_bari$trt_total),
                             paste(d_sari_bari$control_events, "/", d_sari_bari$control_total)),
                  ilab.xpos=c(-4,-2.5),
                  header="Study",
                  atransf=exp,
                  order=d_sari_bari$treatment,
                  pch = 19) # Circles 
  
  abline(h=15) # horizontal line
  
  metafor::addpoly(x = bari_results$y,
                   ci.lb = bari_results$ymin,
                   ci.ub = bari_results$ymax,
                   
                   mlab = "Total [95% CrI]",
                   row=14,
                   efac = 1, # polygon size
                   atransf=exp)
  
  bari = dplyr::filter(d_sari_bari,
                       treatment == "baricitinib")
  
  text(c(-4,-2.5),
       14,
       c(paste(sum(bari$trt_events), "/", sum(bari$trt_total)),
         paste(sum(bari$control_events), "/", sum(bari$control_total))),
       font = 2, cex = 0.75)
  
  # Prediction interval
  colp = "#6b58a6"
  coll = "#a7a9ac"
  
  metafor::addpoly(x = bari_pred$y,
                   ci.lb = bari_pred$ymin,
                   ci.ub = bari_pred$ymax,
                   
                   pi.lb = bari_pred$ymin, # prediction interval
                   pi.ub = bari_pred$ymax, # prediction interval
                   annotate=FALSE,
                   
                   mlab = " ", # label
                   row=13, # location
                   atransf=exp,
                   efac = 0.5, # polygon size
                   col=colp, border=colp
  )
  
  text(-7.95, 13,
       tau_text(ma_bari),
       pos=4,
       cex=0.8)
  
  text(1.55, 13,
       PI_text(bari_pred),
       pos=4,
       cex = 0.85)
  
  abline(h=2) # horizontal line
  
  metafor::addpoly(x = sari_results$y,
                   ci.lb = sari_results$ymin,
                   ci.ub = sari_results$ymax,
                   
                   mlab = "Total [95% CrI]",
                   row=1,
                   efac = 1, # polygon size
                   atransf=exp)
  
  sari = dplyr::filter(d_sari_bari, treatment == "sarilumab")
  
  text(c(-4,-2.5),
       1,
       c(paste(sum(sari$trt_events), "/", sum(sari$trt_total)),
         paste(sum(sari$control_events), "/", sum(sari$control_total))),
       font = 2, cex = 0.75)
  
  # Prediction interval
  
  metafor::addpoly(x = sari_pred$y,
                   ci.lb = sari_pred$ymin,
                   ci.ub = sari_pred$ymax,
                   
                   pi.lb = sari_pred$ymin, # prediction interval
                   pi.ub = sari_pred$ymax, # prediction interval
                   annotate=FALSE,
                   
                   mlab = " ", # label
                   row=0, # location
                   atransf=exp,
                   efac = 0.5, # polygon size
                   col=colp, border=colp
  )
  
  
  
  text(-7.95, -0.15,
       tau_text(ma_sari),
       pos=4,
       cex=0.8)
  
  text(1.55, -0.15,
       PI_text(sari_pred),
       pos=4,
       cex = 0.85)
  
  # Headers
  text(mean(c(-4,-2.5)) - 0.1,
       23.5,
       "No of events / total",
       font=2,
       cex = 0.7)
  
  # Horizontal line
  segments(-4 - 0.6,
           23,
           -2.5 + 0.35,
           23)
  
  text(c(-4,-2.5), 22.5, c("Experimental","Control"), font = 2, cex = 0.7)
  
  text(c(-0.6,0.6), # X axis
       20.5, # Y axis
       cex = 0.7, # size
       c("Favors\nExperimental","Favors\nControl"),
       font = 1,
       pos=c(2,4), # Right + Left aligned
       offset=-1)
  
  text(-7.2,c(20.5, 11.5),
       c("Baricitinib", "Sarilumab"),
       font=4,
       cex = 0.9)
}

# Both together

pdf(file = here::here("output", "figures", 
                      'forest_toci_bari_sari.pdf'),
    height=6.5, width=14)

par(mfrow = c(1, 2), # 2 figures, side by side ("1 row, 2 columns")
    # Margins
    mar=c(6, # bottom
          4, # left
          3, # top
          3) # right
    ) 

forest_fun_toci()
forest_fun_bari_sari()

dev.off()


# Sensitivity analysis ----

d_bari_sens = d_logOR |>
  dplyr::filter(treatment == "baricitinib",
                study != "RECOVERY Bari (No Toci)") |> 
  arrange(desc(study))

bari_results_sens = ma_bari_sens |> brms::fixef(summary = F) |>
  ggdist::median_hdi()

bari_pred_sens = predictions(ma_bari_sens)

forest_sari_sens = function(){
  
  metafor::forest(d_bari_sens$yi, d_bari_sens$vi,
                  cex=0.75, # text size
                  ylim=c(12, 22.5), # Y length
                  rows=c(15:18),
                  alim=log(c(0.25, 4)),
                  xlim=c(-8, 3.5), 
                  at=log(c(0.25, 0.5, 1, 2, 4)),
                  slab=d_bari_sens$study,
                  ilab=cbind(paste(d_bari_sens$trt_events, "/", d_bari_sens$trt_total),
                             paste(d_bari_sens$control_events, "/", d_bari_sens$control_total)),
                  ilab.xpos=c(-4,-2.5),
                  header="Study",
                  atransf=exp,
                  order=d_bari_sens$treatment,
                  pch = 19) # Circles 
  
  abline(h=14) # horizontal line
  
  ### set font expansion factor (as in forest() above) and use a bold font
  par(cex=0.9)
  
  metafor::addpoly(x = bari_results_sens$y,
                   ci.lb = bari_results_sens$ymin,
                   ci.ub = bari_results_sens$ymax,
                   
                   mlab = "Total [95% CrI]",
                   row=13,
                   efac = 1, # polygon size
                   atransf=exp)
  
  text(c(-4,-2.5),
       13,
       c(paste(sum(d_bari_sens$trt_events), "/", sum(d_bari_sens$trt_total)),
         paste(sum(d_bari_sens$control_events), "/", sum(d_bari_sens$control_total))),
       font = 2, cex = 0.75)
  
  # Prediction interval
  colp = "#6b58a6"
  coll = "#a7a9ac"
  
  metafor::addpoly(x = bari_pred_sens$y,
                   ci.lb = bari_pred_sens$ymin,
                   ci.ub = bari_pred_sens$ymax,
                   
                   pi.lb = bari_pred_sens$ymin, # prediction interval
                   pi.ub = bari_pred_sens$ymax, # prediction interval
                   annotate=FALSE,
                   
                   mlab = " ", # label
                   row=12, # location
                   atransf=exp,
                   efac = 0.5, # polygon size
                   col=colp, border=colp
  )
  
  text(-7.95, 12,
       tau_text(ma_bari_sens),
       pos=4,
       cex=0.8)
  
  text(1.6, 12,
       PI_text(bari_pred_sens),
       pos=4,
       cex = 0.85)
  
  
  # Headers
  text(mean(c(-4,-2.5)) - 0.1,
       22.5,
       "No of events / total",
       font=2,
       cex = 0.7)
  
  # Horizontal line
  segments(-4 - 0.6,
           22,
           -2.5 + 0.35,
           22)
  
  text(c(-4,-2.5), 21.5, c("Experimental","Control"), font = 2, cex = 0.7)
  
  text(c(-0.6,0.6), # X axis
       19.5, # Y axis
       cex = 0.7, # size
       c("Favors\nExperimental","Favors\nControl"),
       font = 1,
       pos=c(2,4), # Right + Left aligned
       offset=-1)
  
  text(-7.2, 19.5,
       "Baricitinib",
       font=4,
       cex = 0.9)
}

pdf(file = here::here("output", "figures", "supplementary", 
                      'forest_sari_sens.pdf'),
    height=4.5, width=7)

forest_sari_sens()

dev.off()


