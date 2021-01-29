# ------------------------------------------------------------------------------
# Description:  Script to derive field saturated hydraulic conductivity (kfs)
#               and sorptivity (S) from mini-disk infiltrometer data. Uses the
#               three curve fitting techniques discussed in Vandervaere 2000.
#               If negative values are found for either C1 or C2, the values
#               are dropped. The remaining values are averaged. The mean c1 and
#               c2 values are then divided by a1 and a2, repsectively, as 
#               proposed by Zhang 1997, to derive S and kfs.
# Date:         1/28/21
# Author:       Justin Johnson
# ------------------------------------------------------------------------------
# Install and load packages
pkgs <- c('tidyverse','broom', 'grid', 'gridExtra', 'stringr')

check <- sapply(pkgs,require,warn.conflicts = TRUE,character.only = TRUE)
if(any(!check)){
  pkgs.missing <- pkgs[!check]
  install.packages(pkgs.missing)
  check <- sapply(pkgs.missing,require,warn.conflicts = TRUE,character.only = TRUE)
}

#Inputs
mdifilename <- "raw_data/MDI_timeseries_Pinalenos_2020.csv"
outfilename <- "processed_data/MDI_Pinalenos_2020.csv"

#van Genutchen and soil moisture paramaters. Soil texture averaged sandy loam for all plots.
#Parameters taken from Carsel and Parrish 1988.
soil_texture <- "sandy loam"  #Average soil texture for each plot, hydrometer method
n <- 1.89                     #Carsel and Parrish 1988 for sandy loam
theta_i <- 0.025              #Estimate of initial soil moisture based on soil moisture probe samples.
theta_r <- 0.065
theta_s <- 0.41               #Carsel and Parrish 1988 for sandy loam
R <- 2.25                     #Radius of MDI (cm)
alpha <- 0.075                #Carsel and Parrish 1988 for sandy loam
h0 <- -1                      #Tension head of MDI (cm)

# ------------------------------------------------------------------------------
#Read in data
mdi <- read_csv(mdifilename, col_types = "Dfffnnnnn")

#Calculates cumulative infiltration and converts to cm from mL
mdi1 <- mdi %>%
  group_by(location_ID) %>%
  mutate(time0.5.sec0.5= time.sec^0.5, #Calculates square root of time
         I_time0.5.cm_sec0.5 = I.cm/time0.5.sec0.5, #Divides cumulative infiltration by square root of time
         geomean.time0.5.sec0.5=sqrt(time0.5.sec0.5*lead(time0.5.sec0.5)), #Calculates geometric mean of difference in time
         dIdtime0.5.cm_sec0.5=(lead(I.cm)-I.cm)/(lead(time0.5.sec0.5)-time0.5.sec0.5)) #Calculates derivative of infiltration by sqrt of time
  
# ------------------------------------------------------------------------------
# Cumulative infiltration calculations

# Runs non-linear regression on cumulative infiltration plots and extracts c1.
CI.c1 <- mdi1 %>%
  group_by(location_ID) %>%
  do(tidy(nls(I.cm~c1*time0.5.sec0.5+c2*time.sec,start=list(c1=10,c2=10), data=., control = list(maxiter = 500))))%>%
  filter(term=="c1") %>%
  mutate(c1.CI=estimate) %>%
  select(location_ID, c1.CI)

# Runs non-linear regression on cumulative infiltration plots and extracts c2.
CI.c2 <- mdi1 %>%
  group_by(location_ID) %>%
  do(tidy(nls(I.cm~c1*time0.5.sec0.5+c2*time.sec,start=list(c1=10,c2=10), data=., control = list(maxiter = 500))))%>%
  filter(term=="c2") %>%
  mutate(c2.CI=estimate) %>%
  select(location_ID, c2.CI)

# ------------------------------------------------------------------------------
# Cumulative linearization calculations.
# Visual examination of plots show a consisten breakpoint around 25 seconds,
# which appears to be related to the contact sand. Filtered out data points 
# before 25 seconds as suggested in Vandervaere et al. 2000.

# Runs linear regression on cumulative linearization plots and extracts c1 (intercept.)
CL.c1 <- mdi1 %>%
  group_by(location_ID) %>%
  filter(time.sec > 24) %>%
  do(tidy(lm(I_time0.5.cm_sec0.5~time0.5.sec0.5, data = .)))%>%
  filter(term=="(Intercept)") %>%
  mutate(c1.CL=estimate) %>%
  select(location_ID, c1.CL)

# Runs linear regression on cumulative linearization plots and extracts c2 (slope).
CL.c2 <- mdi1 %>%
  group_by(location_ID) %>%
  filter(time.sec > 24) %>%
  do(tidy(lm(I_time0.5.cm_sec0.5~time0.5.sec0.5, data = .)))%>%
  filter(term=="time0.5.sec0.5") %>%
  mutate(c2.CL=estimate) %>%
  select(location_ID, c2.CL)

# Runs linear regression on cumulative linearization plots and extracts r2
CL.r2 <- mdi1 %>%
  group_by(location_ID) %>%
  do(glance(lm(I_time0.5.cm_sec0.5~time0.5.sec0.5, data = .)))%>%
  mutate(r2.CL=r.squared) %>%
  select(location_ID, r2.CL)

# ------------------------------------------------------------------------------
# Differentiated linearization calculations

# Runs linear regression on differentiated linearization plots and extracts c1 (intercept.)
DL.c1 <- mdi1 %>%
  group_by(location_ID) %>%
  do(tidy(lm(dIdtime0.5.cm_sec0.5~geomean.time0.5.sec0.5, data = .)))%>%
  filter(term=="(Intercept)") %>%
  mutate(c1.DL=estimate) %>%
  select(location_ID, c1.DL)

# Runs linear regression on differentiated linearization plots and extracts c2 (1/2 slope).
DL.c2 <- mdi1 %>%
  group_by(location_ID) %>%
  do(tidy(lm(dIdtime0.5.cm_sec0.5~geomean.time0.5.sec0.5, data = .)))%>%
  filter(term=="geomean.time0.5.sec0.5") %>%
  mutate(c2.DL=estimate/2) %>%
  select(location_ID, c2.DL)

# Runs linear regression on differentiated linearization plots and extracts r2
DL.r2 <- mdi1 %>%
  group_by(location_ID) %>%
  do(glance(lm(dIdtime0.5.cm_sec0.5~geomean.time0.5.sec0.5, data = .)))%>%
  mutate(r2.DL=r.squared) %>%
  select(location_ID, r2.DL)

# ------------------------------------------------------------------------------
# Merges datasets and calculates sorptivity and Kfs
# Negative values of C1 and C2 are thrown out for each method. Remaining values are averaged.
# ------------------------------------------------------------------------------

# Completes merge and drops negative values. 
mdi.merge <- full_join(CI.c1, CI.c2, by="location_ID") %>%
  full_join(CL.c1, by="location_ID") %>%
  full_join(CL.c2, by="location_ID") %>%
  full_join(CL.r2, by="location_ID") %>%
  full_join(DL.c1, by="location_ID") %>%
  full_join(DL.c2, by="location_ID") %>%
  full_join(DL.r2, by="location_ID") %>%
  mutate(c1.CI = replace(c1.CI, c1.CI<0, NA),
         c2.CI = replace(c2.CI, c2.CI<0, NA),
         c1.CL = replace(c1.CL, c1.CL<0, NA),
         c2.CL = replace(c2.CL, c2.CL<0, NA),
         c1.DL = replace(c1.DL, c1.DL<0, NA),
         c2.DL = replace(c2.DL, c2.DL<0, NA))

# Calculates a1, a2, S, Kfs. Converts Kfs to mm/hr and sorptivity to mm/hr^0.5
mdi.calc <- mdi.merge %>%
  mutate(a1=((1.4*0.55^0.5)*((theta_s-theta_i)^0.25)*exp((3*(n-1.9)*alpha*h0)))/((alpha*R)^0.15),
         a2=ifelse(n>=1.9,
                   (11.65*((n^0.1)-1)*exp((2.92*(n-1.9)*alpha*h0)))/((alpha*R)^0.91),
                   (11.65*((n^0.1)-1)*exp((7.5*(n-1.9)*alpha*h0))))/((alpha*R)^0.91),
         c1=mean(c(c1.CI, c1.CL, c1.DL), na.rm=TRUE),
         c2=mean(c(c2.CI, c2.CL, c2.DL), na.rm=TRUE),
         S.cm_sec0.5 = c1/a1,
         S.mm_hr0.5 = S.cm_sec0.5*10*60, #60 sec^0.5 per hr^0.5; 10 mm per cm
         kfs.cm_sec = c2/a2,
         kfs.mm_hr = kfs.cm_sec*3600*10) #3600 sec per hr; 10 mm per sec

# ------------------------------------------------------------------------------

# Formats final dataset and adds plot, location, and burn severity data.
mdi.final <- mdi.calc %>%
  separate(location_ID, c("plot_ID", "x.m", "y.m"), sep = "_", remove = FALSE) %>%
  mutate(plot_ID = as.numeric(plot_ID),
         x.m = as.numeric(x.m),
         y.m = as.numeric(y.m),
         burn_severity = case_when(plot_ID < 4 ~ "high",
                                   plot_ID > 3 & plot_ID < 7 ~ "moderate",
                                   plot_ID > 6 & plot_ID < 10 ~ "low",
                                   plot_ID > 9 ~"unburned"),
         soil_texture_mdi = soil_texture,
         n = n,
         theta_i.cm3_cm3 = theta_i,
         theta_s.cm3_cm3 = theta_s,
         mdiradius.cm = R,
         alpha.percm = alpha,
         tension.cm = h0,
         hf.m = ((S.mm_hr0.5/((2*kfs.mm_hr*(theta_s-theta_r))^(1/2)))^2)/1000) %>%
  select(location_ID, plot_ID, x.m, y.m , burn_severity, S.mm_hr0.5, kfs.mm_hr, hf.m, everything())

# Assess the quality of the fit for each curve fitting technique.
mdi.assess <- mdi.final %>%
  mutate(CI.fit = case_when(plot_ID < 10 ~ "good",
                            location_ID == "10_0_10" ~ "good",
                            location_ID == "10_2_8" ~ "poor",
                            location_ID == "10_3_3" ~ "moderate",
                            location_ID == "10_8_8" ~ "moderate",
                            location_ID == "10_8_2" ~ "good",
                            location_ID == "10_10_0" ~ "good",
                            location_ID == "10_0_1" ~ "moderate",
                            location_ID == "10_0_3" ~ "good",
                            location_ID == "10_0_5" ~ "poor",
                            location_ID == "10_0_0" ~ "good",
                            location_ID == "10_1_0" ~ "good",
                            location_ID == "10_1_1" ~ "moderate",
                            location_ID == "10_3_0" ~ "good",
                            location_ID == "10_5_0" ~ "good",
                            location_ID == "10_10_10" ~ "good",
                            location_ID == "10_5_5" ~ "moderate",
                            location_ID == "11_8_8" ~ "good",
                            location_ID == "11_10_10" ~ "poor",
                            location_ID == "11_5_5" ~ "good",
                            location_ID == "11_3_3" ~ "good",
                            location_ID == "11_2_8" ~ "moderate",
                            location_ID == "11_0_10" ~ "moderate",
                            location_ID == "11_5_0" ~ "good",
                            location_ID == "11_10_0" ~ "good",
                            location_ID == "11_8_2" ~ "moderate",
                            location_ID == "11_0_0" ~ "moderate",
                            location_ID == "11_1_0" ~ "good",
                            location_ID == "11_3_0" ~ "good",
                            location_ID == "11_0_5" ~ "poor",
                            location_ID == "11_0_3" ~ "good",
                            location_ID == "11_0_1" ~ "good",
                            location_ID == "11_1_1" ~ "moderate",
                            location_ID == "12_2_8" ~ "moderate",
                            location_ID == "12_8_8" ~ "good",
                            location_ID == "12_8_2" ~ "good",
                            location_ID == "12_0_3" ~ "good",
                            location_ID == "12_0_5" ~ "good",
                            location_ID == "12_0_10" ~ "good",
                            location_ID == "12_0_0" ~ "good",
                            location_ID == "12_1_0" ~ "good",
                            location_ID == "12_3_0" ~ "moderate",
                            location_ID == "12_0_1" ~ "good",
                            location_ID == "12_1_1" ~ "good",
                            location_ID == "12_5_0" ~ "poor",
                            location_ID == "12_5_5" ~ "good",
                            location_ID == "12_10_0" ~ "good",
                            location_ID == "12_10_10" ~ "good",
                            location_ID == "12_3_3" ~ "good"),
         CL.fit= case_when(plot_ID < 10 ~ "good",
                           location_ID == "10_0_10" ~ "good",
                           location_ID == "10_2_8" ~ "poor",
                           location_ID == "10_3_3" ~ "moderate",
                           location_ID == "10_8_8" ~ "good",
                           location_ID == "10_8_2" ~ "good",
                           location_ID == "10_10_0" ~ "good",
                           location_ID == "10_0_1" ~ "moderate",
                           location_ID == "10_0_3" ~ "good",
                           location_ID == "10_0_5" ~ "poor",
                           location_ID == "10_0_0" ~ "good",
                           location_ID == "10_1_0" ~ "good",
                           location_ID == "10_1_1" ~ "good",
                           location_ID == "10_3_0" ~ "good",
                           location_ID == "10_5_0" ~ "poor",
                           location_ID == "10_10_10" ~ "poor",
                           location_ID == "10_5_5" ~ "moderate",
                           location_ID == "11_8_8" ~ "good",
                           location_ID == "11_10_10" ~ "poor",
                           location_ID == "11_5_5" ~ "good",
                           location_ID == "11_3_3" ~ "good",
                           location_ID == "11_2_8" ~ "moderate",
                           location_ID == "11_0_10" ~ "moderate",
                           location_ID == "11_5_0" ~ "good",
                           location_ID == "11_10_0" ~ "good",
                           location_ID == "11_8_2" ~ "poor",
                           location_ID == "11_0_0" ~ "good",
                           location_ID == "11_1_0" ~ "good",
                           location_ID == "11_3_0" ~ "good",
                           location_ID == "11_0_5" ~ "poor",
                           location_ID == "11_0_3" ~ "good",
                           location_ID == "11_0_1" ~ "poor",
                           location_ID == "11_1_1" ~ "moderate",
                           location_ID == "12_2_8" ~ "moderate",
                           location_ID == "12_8_8" ~ "poor",
                           location_ID == "12_8_2" ~ "good",
                           location_ID == "12_0_3" ~ "good",
                           location_ID == "12_0_5" ~ "good",
                           location_ID == "12_0_10" ~ "good",
                           location_ID == "12_0_0" ~ "good",
                           location_ID == "12_1_0" ~ "good",
                           location_ID == "12_3_0" ~ "moderate",
                           location_ID == "12_0_1" ~ "good",
                           location_ID == "12_1_1" ~ "good",
                           location_ID == "12_5_0" ~ "poor",
                           location_ID == "12_5_5" ~ "good",
                           location_ID == "12_10_0" ~ "poor",
                           location_ID == "12_10_10" ~ "poor",
                           location_ID == "12_3_3" ~ "good"),
         DL.fit= case_when(plot_ID < 10 ~ "good",
                           location_ID == "10_0_10" ~ "good",
                           location_ID == "10_2_8" ~ "moderate",
                           location_ID == "10_3_3" ~ "good",
                           location_ID == "10_8_8" ~ "good",
                           location_ID == "10_8_2" ~ "good",
                           location_ID == "10_10_0" ~ "good",
                           location_ID == "10_0_1" ~ "good",
                           location_ID == "10_0_3" ~ "good",
                           location_ID == "10_0_5" ~ "poor",
                           location_ID == "10_0_0" ~ "good",
                           location_ID == "10_1_0" ~ "good",
                           location_ID == "10_1_1" ~ "good",
                           location_ID == "10_3_0" ~ "good",
                           location_ID == "10_5_0" ~ "moderate",
                           location_ID == "10_10_10" ~ "poor",
                           location_ID == "10_5_5" ~ "moderate",
                           location_ID == "11_8_8" ~ "good",
                           location_ID == "11_10_10" ~ "poor",
                           location_ID == "11_5_5" ~ "good",
                           location_ID == "11_3_3" ~ "good",
                           location_ID == "11_2_8" ~ "good",
                           location_ID == "11_0_10" ~ "moderate",
                           location_ID == "11_5_0" ~ "good",
                           location_ID == "11_10_0" ~ "good",
                           location_ID == "11_8_2" ~ "moderate",
                           location_ID == "11_0_0" ~ "good",
                           location_ID == "11_1_0" ~ "good",
                           location_ID == "11_3_0" ~ "moderate",
                           location_ID == "11_0_5" ~ "poor",
                           location_ID == "11_0_3" ~ "good",
                           location_ID == "11_0_1" ~ "moderate",
                           location_ID == "11_1_1" ~ "moderate",
                           location_ID == "12_2_8" ~ "moderate",
                           location_ID == "12_8_8" ~ "poor",
                           location_ID == "12_8_2" ~ "good",
                           location_ID == "12_0_3" ~ "good",
                           location_ID == "12_0_5" ~ "good",
                           location_ID == "12_0_10" ~ "good",
                           location_ID == "12_0_0" ~ "good",
                           location_ID == "12_1_0" ~ "good",
                           location_ID == "12_3_0" ~ "moderate",
                           location_ID == "12_0_1" ~ "good",
                           location_ID == "12_1_1" ~ "good",
                           location_ID == "12_5_0" ~ "poor",
                           location_ID == "12_5_5" ~ "moderate",
                           location_ID == "12_10_0" ~ "moderate",
                           location_ID == "12_10_10" ~ "moderate",
                           location_ID == "12_3_3" ~ "good"))

#Write final dataset.
write_csv(mdi.assess, outfilename)

# ------------------------------------------------------------------------------
#Creates function to make a grob table that as a single grob instead of a grob for each table cell (makes processing very slow).
# grid.ftable <- function(d, padding = unit(4, "mm"), ...) {
# 
#   nc <- ncol(d)
#   nr <- nrow(d)
# 
#   ## character table with added row and column names
#   extended_matrix <- cbind(c(rownames(d)),
#                            rbind(colnames(d),
#                                  as.matrix(d)))
# 
#   ## string width and height
#   w <- apply(extended_matrix, 2, strwidth, "inch")
#   h <- apply(extended_matrix, 2, strheight, "inch")
# 
#   widths <- apply(w, 2, max)
#   heights <- apply(h, 1, max)
# 
#   padding <- convertUnit(padding, unitTo = "in", valueOnly = TRUE)
# 
#   x <- cumsum(widths + padding) - 0.5 * padding
#   y <- cumsum(heights + padding) - padding
# 
#   rg <- rectGrob(x = unit(x - widths/2, "in"),
#                  y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
#                  width = unit(widths + padding, "in"),
#                  height = unit(heights + padding, "in"))
# 
#   tg <- textGrob(c(t(extended_matrix)), x = unit(x - widths/2, "in"),
#                  y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
#                  just = "center")
# 
#   g <- gTree(children = gList(rg, tg), ...,
#              x = x, y = y, widths = widths, heights = heights)
# 
#   grid.draw(g)
#   invisible(g)
# }
# 
# #Subsets final data to print associated data with infiltration curves.
# mdi.print <- mdi.final %>%
#   select(location_ID, burn_severity, soil_texture_mdi, S.mm_hr0.5, kfs.mm_hr, hf.m, c1, c2, a1, a2)
# 
# #Writes figures of cumulative infiltration, cumulative linearization, and differenitated linearization plots for each MDI. Also includes
# #table with calculated values.
# z <- list()
# for(i in unique(mdi1$location_ID)){
#  CI <- ggplot(subset(mdi1,location_ID==i),aes(time.sec,I.cm))+
#     geom_point()+
#    geom_smooth(method="nls",
#                formula=y~c1*x^0.5+c2*x,
#                method.args = list(start=c(c1=10,c2=10)),
#                se=FALSE)+
#     labs(x="t (s)", y="I (cm)", title="Cumulative Infiltration")
# 
#  CL <-ggplot(subset(mdi1,location_ID==i),aes(time0.5.sec0.5,I_time0.5.cm_sec0.5))+
#    geom_point()+
#    geom_smooth(data = subset(mdi1,location_ID==i & time.sec > 24),method="lm")+
#    labs(x=expression("t"^{1/2}*" (s"^{1/2}*")"), y=expression("I/t"^{1/2}*" (cm·s"^{-1/2}*")"), title="Cumulative Linearization")
# 
#  DI <- ggplot(subset(mdi1,location_ID==i),aes(geomean.time0.5.sec0.5,dIdtime0.5.cm_sec0.5))+
#    geom_point()+
#    geom_smooth(method="lm")+
#    labs(x=expression("t"^{1/2}*" (s"^{1/2}*")"), y=expression("dI/dt"^{1/2}*" (cm·s"^{-1/2}*")"), title="Cumulative Differentiation")
# 
#  table <- grid.ftable(t(mdi.print[which(mdi.final$location_ID == i),]))
# 
#  z[[i]] <- arrangeGrob(CI,CL,DI,table,nrow=2,top=i)
# }
# 
# pdf("figures/exploratory_figures/MDI_plots_combined.pdf")
# for(p in z){
#   grid.arrange(p)
# }
# dev.off()

# ------------------------------------------------------------------------------
# Quality control kfs and S estimates.
drop_poor.CI <- mdi.assess %>%
  filter(CI.fit != "poor") %>%
  group_by(burn_severity) %>%
  summarize(meanC1.CI = mean(c1.CI, na.rm = TRUE),
            meanC2.CI = mean(c2.CI, na.rm = TRUE))

keep.CI <- mdi.assess %>%
  group_by(burn_severity) %>%
  summarize(meanC1.CI = mean(c1.CI, na.rm = TRUE),
            meanC2.CI = mean(c2.CI, na.rm = TRUE))

drop_poor_med.CI <- mdi.assess %>%
  filter(CI.fit == "good") %>%
  group_by(burn_severity) %>%
  summarize(meanC1.CI = mean(c1.CI, na.rm = TRUE),
            meanC2.CI = mean(c2.CI, na.rm = TRUE))
###
# Removing poor model fits results in dropping primarily low Kfs and S estimates.
# Although model fit isn't great, since error is systematic, removal of estimates
# biases sample. The trade off is measurment error vs sampling bias. Chose to keep 
# estimates.
###
