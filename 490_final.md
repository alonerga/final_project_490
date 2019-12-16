Introduction
============

This is an R markdown file showing analysis of research I conducted in
Dr. Andrew Steen’s biogeochemistry lab. This markdown can be used as a
reference for extracellular enzyme assay analyses, or in order to
reproduce the analysis we did.

Research background
-------------------

The goal of this project is to determine how extracellular enzymes are
reacting with substrates in the various hot springs in which the samples
were collected. These samples come from different geothermal features
across Panama, where there are unique biogeochemical interactions due to
the convergent margins, volcanic arcs, microbes, and carbon fluxes.
Enzymatic affects on carbon fluxes in these types of environments have
been studied very little. We use fluorogenic substrate proxies to
represent peptidases, polysaccharide hydrolases, sulfatases, and
phosphotases, to get an idea of what these enzymes are breaking down.
Most extracellular enzymes follow Michaelis-Menten kinetics, where the
initial hydrolysis rate of an enzyme increases with substrate
concentration, until it reaches an equillibrium hydrolysis rate.
Michaelis-Menten kinetics is shown in the equation and image below.

$$
\\large v\_o = \\frac{V\_{max} \\times \[S\]}{K\_m + \[S\]}
$$
<p align="center">
<img src = "https://upload.wikimedia.org/wikipedia/commons/8/83/Michaelis_Menten_curve_2.svg" width = "550" height = "275">
</p>
### **This analysis will provide us with *v*<sub>*o*</sub>, the rate of reaction of micromoles of substrate being broken down per liter per hour at each site**

Initial steps
=============

We will use a variety of packages and sets of packages to analyze this
data in R. Tidyverse is a set of packages meant for data science.
Lubridate helps ease the complications surrounding dates and times in R.
Devtools creates a directory structure for your packages. Lmstats
provides vector summaries of linear models.

**tidyverse**, **lubridate**, **devtools**, and **lmstats**: Install the
packages first, then load them in

    library(tidyverse)
    library(lubridate)
    library(devtools)
    library(lmstats)

Thermal degradation test
------------------------

A test of substrates and their stability at high temperatures. We want
to ensure that these substrates aren’t breaking apart on their own, that
way when we have live samples, any fluorometer readings we get can
represent enzymatic activity.  
Read in your data from wherever it is stored in your drive. You can use
glimpse to make sure it looks correct

    thermaldeg_data <- read_csv("data/2018_10_31_thermaldeg_markdown.csv")

    glimpse(thermaldeg_data)

\#\#\#\#**Note**  
Some of my code will not show in this document to keep it concise. For
example, I have ran some mutations using lubridate to ensure my time
measurements are placed on a similar scale (elapsed time) instead of the
exact time of day I took the measurements. The full code for this
analysis can be found on my github page.

Plot raw thermal degradation data
---------------------------------

Ggplot is a great package for data visualization. It is relatively
simple to use and creates professional looking plots. We want to show
the change in relative fluorescence units over time for our four
substrates and buffer. A linear model is a great way to plot this and we
can add a shaded region around each line for our standard error. The
code below allows us to show this.

    raw_thermalp <- ggplot(thermaldeg_data, aes(x = elapsed.time, 
                                               y = RFU, 
                                               colour = substrate, 
                                               shape = as.factor(replicate))) + 
      
      geom_point() +
      geom_smooth(method = "lm", se = TRUE) +
      expand_limits(ymin = 0) + 
      scale_shape_discrete(name = "replicate") +

      facet_grid(substrate~temperature, scales = "free_y")

    print(raw_thermalp)

![](490_final_files/figure-markdown_strict/unnamed-chunk-4-1.png)

Enzyme assay data
=================

After showing that our substrates are relatively stable at high
temperatures, we are ready to move forward with our enzyme analysis.
Read in your data from your live sample fluorometer measurements.

\#\#\#\#**Note**  
Similar to the thermal degradation data, I have made a few adjustments
to make this code more readily usable for the rest of my analysis. This
code can be found on my github.

Plot the raw data
-----------------

After you have loaded in your data and made any neccessary manipulations
to the time, names, or any other parts, you can now plot the raw data.
We will use ggplot again for this. Similar to the thermal degradation
plot, we want a linear model of the change in RFU over time. Using facet
wrap or facet grid is helpful to break apart your sites and substrates
and show them in an easily understood way.

    raw_panamap <- ggplot(panamadata, aes(x = elapsed.hr,
                                         y = RFU,
                                         color = filter.type)) +
      
      geom_point() +
      geom_smooth(aes(group = Cuvette),
                  method = "lm",
                  se = FALSE) +
      
      facet_wrap(Substrate~Site, scales = "free") +
      
      theme(text = element_text(size = 6))

    print(raw_panamap)

![](490_final_files/figure-markdown_strict/unnamed-chunk-7-1.png)

Begin analysis
==============

Here we create a function that will create a dataframe comparing RFU to
elapsed time. Lm\_stats provides us with information such as the slope,
error, and p value, which we will need later on.

    lm_fun <- function(df) {
      df <- as.data.frame(df)
      lm_stats(df, "elapsed.hr", "RFU")
    }

Calculate slope, standard error and p value
-------------------------------------------

Using the function we just made, we can use map() to extract a vector
summary of the slopes, standard errors, and p-values of our Panama data.
Nesting will create a list-column of of these data frames, making it
more tidy.

    slopes <- panamadata %>%
      nest() %>%
      mutate(lmstats_df = map(data, lm_fun)) %>%
      mutate(uncal.slope = map_dbl(lmstats_df, function(x) x[1, "slope"]),
             uncal.slope.se = map_dbl(lmstats_df, function(x) x[1, "slope.se"]),
             rsq = map_dbl(lmstats_df, function(x) x[1, "pval"]))

Make a column in the slopes data for which fluorophore goes with which
substrate. We used leucine-AMC, arginine-AMC, MUB-sulfatase,
MUB-b-cellobioside, and MUB-PO4. You should use the fluorophores which
are appropriate for your data. Then we merge our two dataframes that
have matching rows, by substrate.

    fluorophore_lookup <- data.frame(Substrate = unique(panamadata$Substrate), 
                                     fluorophore = c("AMC", "AMC", "MUB", "MUB", "MUB"))
    slopes <- slopes %>%
      left_join(fluorophore_lookup, by = "Substrate")

Calibration data
================

Our calibration data is how RFU changes with fluorophore concentration.
We will ultimately combine this data with our previous analysis to solve
for *v*<sub>*o*</sub>. Read in your data.

Plot calibration data
---------------------

Using ggplot in a similar fashion as the previous plots.

    plot_calib <- ggplot(calib, aes(x = conc.uM, 
                                 y = RFU, 
                                 colour = as.factor(temperature))) +
      
      geom_point() +
      geom_smooth(method = "lm", se = FALSE) + 
      scale_color_discrete(name = "temperature") +
      
      facet_wrap(~fluorophore)

    print(plot_calib) 

![](490_final_files/figure-markdown_strict/unnamed-chunk-12-1.png)

Create another function
-----------------------

This function is exactly like the last, just changing from time to
concentration of our substrates.

    calib_lm_fun <- function(df) {
      df <- as.data.frame(df)
      lm_stats(df, "conc.uM", "RFU")
    }

Calculate slopes of the calibration curve
-----------------------------------------

Again, a similar analysis as the RFU/time. We only need the type of
fluorophore and the slope, which is what select() provides us.

    calib_slopes <- calib %>%
      group_by(fluorophore) %>%
      nest() %>%
      mutate(lmstats_df = map(data, calib_lm_fun),
             calib.slope = map_dbl(lmstats_df, function(x) x[1, "slope"])) %>%
      select(fluorophore, calib.slope)

Combine datasets
================

Here we want to combine our two datasets, dividing our initial slopes by
the calibration slopes. We can use our RFU/time and RFU/concentration
data to calculate *v*<sub>*o*</sub>. If you used replicates with varying
filter weights, be sure to find the mean filter weight/measurement. Also
take into account the amount of sample fluid that was filtered through
each syringe, in my case 0.120 liters.

    mean.filter.weight = 10 * mean(panamadata$filter.weight)

    rates <- slopes %>%
      left_join(calib_slopes, by = "fluorophore") %>%
      mutate(v0.unnorm = uncal.slope / calib.slope,
             v0.unnorm.se = uncal.slope.se / calib.slope,
             v0.norm = uncal.slope * 10 / 0.120,
             v0.norm.se = uncal.slope.se * (mean.filter.weight / filter.weight) / 0.120) 

Calculate mean rates
--------------------

Here we want to summarize our rates by calculating the mean rate of
*v*<sub>*o*</sub>, based on where the sample was taken, the susbtrate
used, and whether it was a live sample or control sample.

    rates_summ <- rates %>%
      group_by(Site, Substrate, filter.type) %>%
      summarise(mean.v0 = mean(v0.norm, na.rm = TRUE),
                v0.sd = sd(v0.norm, na.rm  = TRUE))

Final plot
==========

We want the mean *v*<sub>*o*</sub> plotted based on site, showing both
the live and control samples. Our ranges for each point will take into
account standard deviation from the mean. Y lab can be used to change
the y axis label, and anything in the theme() section just alters the
appearance of the plot.

    p_rates <- ggplot(rates_summ, aes(x = Site, 
                                      y = mean.v0, 
                                      shape = filter.type)) + 
      
      geom_point() + 
      geom_pointrange(aes(ymin = mean.v0 - v0.sd,
                          ymax = mean.v0 + v0.sd)) +
      
      ylab(expression(paste(v[0], ", ", mu, "mol ", L^{-1}, " ", hr^{-1}))) +
      scale_shape_manual(values = c(1, 19), 
                         name = "Treatment") +
     
      facet_wrap(~Substrate, scales = "free_y") +
      
      theme_classic() +
      
      theme(
            strip.text = element_text(size = 10),
            legend.position = c(0.85, 0.10),
            text = element_text(size = 10),
            axis.text.x = element_text(angle = -45, hjust = 0))

    print(p_rates)         

![](490_final_files/figure-markdown_strict/unnamed-chunk-17-1.png)

Conclusion
==========

You should now have a plot showing your enzymatic activity for each
substrate based at each site. Refer to
<a href="https://github.com/alonerga/final_project_490.git" class="uri">https://github.com/alonerga/final_project_490.git</a>
for any questions about the code.

References
----------

Barry, P.H., Moor, J.M., Giovannelli, D. et al. Forearc carbon sink
reduces long-term volatile recycling into the mantle. Nature 568,
487–492 (2019)
<a href="doi:10.1038/s41586-019-1131-5" class="uri">doi:10.1038/s41586-019-1131-5</a>

Michaelis Menten Kinetics curve image, retrieved from
<a href="https://upload.wikimedia.org/wikipedia/commons/8/83/Michaelis_Menten_curve_2.svg" class="uri">https://upload.wikimedia.org/wikipedia/commons/8/83/Michaelis_Menten_curve_2.svg</a>

Schmidt, Jenna Marie, Microbial Extracellular Enzymes in Marine
Sediments: Methods Development and Potential Activities in the Baltic
Sea Deep Biosphere. Master’s Thesis, University of Tennessee (2016)
<a href="https://trace.tennessee.edu/utk_gradthes/4072" class="uri">https://trace.tennessee.edu/utk_gradthes/4072</a>
