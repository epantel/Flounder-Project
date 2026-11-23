
## download data -------- 
setwd("/Users/davidtruhe/Desktop/Methods and Analyses of Fisheries Ecology")
flounder.det<-read.csv("Flounder_dets.csv")
flounder.tagspec<-read.csv("Flounder_Tag_Meta_Combined.csv")

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(lme4)   
library(readr)
library(tibble)    

flounder.tagspec<-flounder.tagspec %>%
  mutate(
    DateTime = mdy(DATE_TAGGED),
    Date = as.Date(DateTime),
    Full.ID = TRANSMITTER)

flounder.det2 <- full_join(
  
  flounder.det, flounder.tagspec, by = "Full.ID")



flounder.det2<- flounder.det2 %>%
  mutate(
    Est.Tag.Death = Date + EST_TAG_LIFE
  )



# r checking for weird detections --------

flounder.det2$DateTime <- ymd_hms(
  flounder.det2$Time..UTC.,
  tz = "UTC"
)
flounder.det2$det.Date <- as.Date(flounder.det2$Time..UTC.)

flounder.det2 <- flounder.det2 %>%
  group_by(Full.ID) %>%
  mutate(
    last.det = max(DateTime, na.rm = TRUE)
  )

post.death.summary <- flounder.det2 %>%
  filter(DateTime > Est.Tag.Death) %>%
  group_by(Full.ID) %>%
  summarise(
    n_post_death = n(),
    stations = paste(unique(Station), collapse = ", "),
    first_post_death = min(DateTime, na.rm = TRUE),
    last_post_death = max(DateTime, na.rm = TRUE),
    est_tag_death = first(Est.Tag.Death),
    time_after_death = last_post_death - est_tag_death,
    .groups = "drop"
  )

#check for detections on DateTime after est.tag.death
#Full.ID
#last date of detection on DateTime -> last.det
#number of detections after est tag death unique(DateTime)
#unique(Station) detected after est.tag.death 
#grace period is 3-4 months > last.det

#questionable fish
dead1<-subset(flounder.det2, Full.ID == "A69-1602-60632") #not dead

fake1<-subset(flounder.det2, Full.ID == "A69-1602-60636") #actually real

fake2<-subset(flounder.det2, Full.ID == "A69-1602-60671") #any detection after 11/2021 is probably fake 

fake3<-subset(flounder.det2, Full.ID == "A69-1602-60681") #anything after 4/2022 is probs fake

fake4<-subset(flounder.det2, Full.ID == "A69-1602-60693") #FAKE BITCH alll of 1 detection

#removing fake detections
flounder.det3 <- flounder.det2 %>%
  filter(
    !(Full.ID == "A69-1602-60671" & DateTime > ymd_hms("2021-11-19 23:59:59", tz = "UTC")),
    !(Full.ID == "A69-1602-60681" & DateTime > ymd_hms("2022-04-18 23:59:59", tz = "UTC")),
    Full.ID != "A69-1602-60693")

#simultaneous detections
#group by Full.ID
#detections (DateTime) within seconds/min of each other 
#Station on those detections 
#stations preceding proceeding the simultaneous ones to filter out fake 

simultaneous.det <- flounder.det2 %>%
  arrange(Full.ID, DateTime) %>%
  group_by(Full.ID) %>%
  mutate(
    next.time = lead(DateTime),
    next.station = lead(Station),
    time.diff.min = as.numeric(
      difftime(next.time, DateTime, units = "min")
    )
  ) %>%
  filter(
    time.diff.min <= 3,
    Station != next.station
  ) %>%
  ungroup()
simultaneous.det2 <- simultaneous.det %>%
  transmute(
    Full.ID,
    DateTime = DateTime,
    preceding_station = Station,
    following_station = next.station,
    time_diff_min = time.diff.min
  )
simultaneous.det3 <- simultaneous.det2 %>%
  filter(
    sub("_*", "", preceding_station) ==
      sub("_*", "", following_station)
  )
#checking to see if there are any obvious simultaneous detections based on detections not within the same waterbody
simultaneous.det3 <- simultaneous.det %>%
  filter(
    sub("_*", "", Station) ==
      sub("_*", "", next.station)
  )

#Receiver Serial without Station
#102350
#136874
#136875 could be MSS_GPC_LT_33
#129024 could be MSS_WG_Green_11? most likely not



#some tags that did not get detected
#A69-1602-60648
#A69-1602-60701
#A69-1602-60713
#A69-1602-60718
#A69-1604-9081


flounder.det.clean<-flounder.det3 %>%
  filter(
    Full.ID != "A69-1602-60648",
    Full.ID!="A69-1602-60701",
    Full.ID != "A69-1602-60713",
    Full.ID !="A69-1602-60718",
    Full.ID!= "A69-1604-9081"
  )



## r system detections ---------
#this code just checks for system transitions and how long it took the fish to be detected there. 
#i usually sort by time in descending order and give a brief look at the station changes just to see if some make sense
# long times inbetween stations are not uncommon for marine/esturarine fish that may travel offshore
#but long times in between inshore stations can be a red flag

system.transitions <- flounder.det.clean %>%
  arrange(Full.ID, DateTime) %>%
  group_by(Full.ID) %>%
  mutate(
    from_station = lag(Station),
    from_time = lag(DateTime),
    from_system = lag(sub("_.*", "", Station)),
    
    to_station = Station,
    to_time = DateTime,
    to_system = sub("_.*", "", Station),
    
    time_minutes = as.numeric(
      difftime(to_time, from_time, units = "mins")
    ),
    
    time_hours = as.numeric(
      difftime(to_time, from_time, units = "hours")
    )
  ) %>%
  ungroup() %>%
  filter(
    !is.na(from_station),
    from_system != to_system
  ) %>%
  select(
    Full.ID,
    from_system,
    from_station,
    from_time,
    to_system,
    to_station,
    to_time,
    time_minutes,
    time_hours
  ) %>%
  ungroup()

flounder.det.clean <- flounder.det.clean %>%
  filter(
    !(Full.ID == "A69-1602-60674" & DateTime > ymd_hms("2022-07-16 23:59:59", tz = "UTC")),
    !(Full.ID == "A69-1602-60681" & DateTime > ymd_hms("2022-08-09 23:59:59", tz = "UTC")))

#A69-1602-60674 7/17/2022 PR_Cumbest probs fake
#A69-1602-60704 8/10/2022 PR_Cumbest probs fake 

#Remove fish detected on fewer than 5 days ------------------------------
# Counts distinct calendar days with at least one detection (all receivers).
# Runs after the fake-detection trimming above so trimmed days are not counted.
min.det.days <- 5

fish.det.days <- flounder.det.clean %>%
  ungroup() %>%
  group_by(Full.ID) %>%
  summarise(n_det_days = n_distinct(as.Date(DateTime)), .groups = "drop")

short.fish <- fish.det.days %>% filter(n_det_days < min.det.days)
cat("Removing", nrow(short.fish), "fish detected on <", min.det.days, "days\n")
print(short.fish)

flounder.det.clean <- flounder.det.clean %>%
  filter(!Full.ID %in% short.fish$Full.ID)

movin.flounda<-unique(system.transitions$Full.ID)
flounda<-unique(flounder.det.clean$Full.ID)


movin.flounda2 <- data.frame(
  Full.ID = flounda[!flounda %in% movin.flounda]
)


# RI ------
denominator <- "detection_period"   # "detection_period", "study_period", "custom"
min_detections_per_day <- 1         # raise to 2+ to reduce false-detection influence

# command F for Station and switch to System if needed for the remainder of code

#remove all fish that move to barrier island if needed before running (didn't do this, not sure if I was supposed to)

# Clean and add date ---------------------------------------------------
flounder.ri <- flounder.det.clean %>%
  mutate(DateTime = as.POSIXct(DateTime, tz = "UTC"),
         date     = as_date(DateTime)) %>%
  filter(!is.na(DateTime), !is.na(Full.ID), !is.na(Station)) %>%
  distinct(Full.ID, Station, DateTime, .keep_all = TRUE)   # drop duplicates

# Daily detections ---------------------------------------------
daily <- flounder.ri %>%
  count(Full.ID, Station, date, name = "n_det") %>%
  filter(n_det >= min_detections_per_day)

# Monitoring period (denominator) per tag -----------------------------
windows <- switch(
  denominator,
  detection_period = flounder.ri %>%
    group_by(Full.ID) %>%
    summarise(start = min(date), end = max(date), .groups = "drop"),
  study_period = flounder.ri %>%
    distinct(Full.ID) %>%
    mutate(start = study_start, end = study_end),
  custom = tag_windows,
  stop("denominator must be 'detection_period', 'study_period' or 'custom'")
) %>%
  mutate(days_monitored = as.integer(end - start) + 1)

# Residency index: per tag x Station ----------------------------------
ri_Station <- daily %>%
  left_join(windows, by = "Full.ID") %>%
  filter(date >= start, date <= end) %>%                # keep days inside window
  group_by(Full.ID, Station, days_monitored) %>%
  summarise(days_detected = n_distinct(date), .groups = "drop") %>%
  mutate(RI = days_detected / days_monitored) %>%
  # fill in zeros for Station/tag combinations with no detections
  complete(Full.ID, Station, fill = list(days_detected = 0L, RI = 0)) %>%
  group_by(Full.ID) %>%
  mutate(days_monitored = max(days_monitored, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Full.ID, desc(RI))

# Residency index: per tag (all Stations combined) --------------------
ri_overall <- daily %>%
  left_join(windows, by = "Full.ID") %>%
  filter(date >= start, date <= end) %>%
  group_by(Full.ID, days_monitored) %>%
  summarise(days_detected = n_distinct(date), .groups = "drop") %>%
  mutate(RI_overall = days_detected / days_monitored)

# Residency index by month (tag x Station x month) --------------------
ri_monthly <- daily %>%
  mutate(month = floor_date(date, "month")) %>%
  group_by(Full.ID, Station, month) %>%
  summarise(days_detected = n_distinct(date), .groups = "drop") %>%
  mutate(days_in_month = days_in_month(month),
         RI_month = days_detected / days_in_month)

# Population summary by Station ---------------------------------------
ri_summary <- ri_Station %>%
  group_by(Station) %>%
  summarise(n_tags  = n_distinct(Full.ID),
            mean_RI = mean(RI),
            sd_RI   = sd(RI),
            min_RI  = min(RI),
            max_RI  = max(RI),
            .groups = "drop")

print(ri_Station)
print(ri_overall)
print(ri_summary)

# Plots ---------------------------------------------------------------
# (a) Heatmap of RI per tag and Station
p_heat <- ggplot(ri_Station, aes(Station, Full.ID, fill = RI)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", RI)), size = 3) +
  scale_fill_viridis_c(limits = c(0, 1), option = "C") +
  labs(x = "Station", y = "Tag ID", fill = "RI",
       title = "Residency index by individual and Station") +
  theme_minimal()

# (b) Boxplot of RI by Station
p_box <- ggplot(ri_Station, aes(Station, RI)) +
  geom_boxplot(outlier.shape = NA, fill = "grey85") +
  geom_jitter(width = 0.15, alpha = 0.7) +
  ylim(0, 1) +
  labs(x = "Station", y = "Residency index",
       title = "Residency index distribution by Station") +
  theme_minimal()

# (c) Monthly RI through time
p_month <- ggplot(ri_monthly, aes(month, RI_month, colour = Station, group = Station)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun = mean, geom = "point") +
  ylim(0, 1) +
  labs(x = NULL, y = "Mean monthly RI", colour = "Station",
       title = "Mean monthly residency index") +
  theme_minimal()

print(p_heat); print(p_box); print(p_month)

# ANOVA ---------------------------------------------------------------
ri_Station <- ri_Station %>%
  mutate(
    Station   = trimws(Station),
    System = case_when(
      is.na(Station) | Station == "" ~ NA_character_,
      grepl("LT", Station)           ~ "BI",
      TRUE                           ~ sub("_.*", "", Station)))

ri_Station_clean <- ri_Station |> filter(System %in% c("BB", "BSL", "DB", "MSS", "PR", "PRE"))

#Remove fihhh
ri_Station_clean <- ri_Station_clean %>%
  group_by(Full.ID) %>%
  filter(any(RI > 0)) %>%
  ungroup()

#Filter fihhh detected for less than 10 days 
ri_Station_clean <- ri_Station_clean %>%
  filter(days_monitored > 9)

ri_Station_all <- ri_Station %>%  #Keeps barrier island RI values
  filter(days_monitored > 9)

anova <- aov(RI ~ Station, data = ri_Station_all)
anova <- aov(RI ~ System, data = ri_Station_all)
anova <- aov(RI ~ System, data = ri_Station_clean)

summary(anova)

tukey <- TukeyHSD(anova, conf.level = 0.95)
print(tukey)


#plots 
aov_system <- aov(RI ~ System, data = ri_Station_clean)
summary(aov_system)
tukey <- TukeyHSD(aov_system, conf.level = 0.95)

# RI by System (boxplot + individual points + mean) ----
p_system <- ggplot(ri_Station_clean, aes(System, RI, fill = System)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  ylim(0, 1) +
  labs(x = "System", y = "Residency index",
       title = "Residency index by estuary system",
       subtitle = "Diamonds show group means") +
  theme_minimal() +
  theme(legend.position = "none")

# Tukey HSD pairwise differences ----
tukey_df <- as.data.frame(tukey$System) |>
  rownames_to_column("comparison") |>
  rename(p_adj = `p adj`) |>
  mutate(significant = p_adj < 0.05)

p_tukey <- ggplot(tukey_df, aes(diff, reorder(comparison, diff), colour = significant)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(xmin = lwr, xmax = upr)) +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "firebrick")) +
  labs(x = "Difference in mean RI (95% CI)", y = NULL,
       colour = "p < 0.05",
       title = "Tukey HSD pairwise comparisons") +
  theme_minimal()

print(p_system)
print(p_tukey)





# RI by system, spawning vs non-spawning season-----
spawn.months           <- c(11, 12, 1)  
min_detections_per_day <- 1                
min.days.monitored     <- 10               

#Remove barrier island receivers ------------------------------------
flounder.est <- flounder.det.clean %>%
  ungroup() %>%
  mutate(
    Station = trimws(Station),
    System  = if_else(grepl("LT", Station), "BI", sub("_.*", "", Station)),
    date    = as_date(DateTime)
  ) %>%
  filter(!is.na(DateTime), !is.na(Full.ID), !is.na(Station), Station != "",
         System != "BI") %>%
  distinct(Full.ID, Station, DateTime, .keep_all = TRUE)

# Drop fish with fewer than `min.det.days` detection days once barrier island
# receivers are gone (a fish may have had most of its days at the barrier islands).
min.det.days <- 5
est.short <- flounder.est %>%
  group_by(Full.ID) %>%
  summarise(n_det_days = n_distinct(date), .groups = "drop") %>%
  filter(n_det_days < min.det.days)
flounder.est <- flounder.est %>% filter(!Full.ID %in% est.short$Full.ID)

# 2. Assign each fish to its tagging system (DRAINAGE) -----------
fish.system <- flounder.est %>%
  group_by(Full.ID) %>%
  summarise(tag.system = toupper(trimws(first(na.omit(DRAINAGE)))), .groups = "drop") %>%
  # combine Biloxi and Back Bay into one drainage
  mutate(tag.system = if_else(tag.system %in% c("BILOXI", "BACK BAY"),
                              "BILOXI/BACK BAY", tag.system)) %>%
  filter(tag.system %in% c("BSL", "PASCAGOULA", "BILOXI/BACK BAY"))

flounder.est <- flounder.est %>% filter(Full.ID %in% fish.system$Full.ID)


home.map <- tribble(
  ~tag.system,        ~System,
  "BSL",              "BSL",
  "BILOXI/BACK BAY",  "BB",
  "BILOXI/BACK BAY",  "DB",     
  "PASCAGOULA",       "PR",
  "PASCAGOULA",       "PRE"
)

flounder.est %>%
  inner_join(fish.system, by = "Full.ID") %>%
  count(tag.system, System) %>%
  arrange(tag.system, desc(n)) %>%
  print(n = Inf)

# 3. Days monitored per fish per season (denominator) -------------------
windows <- flounder.est %>%
  group_by(Full.ID) %>%
  summarise(start = min(date), end = max(date), .groups = "drop") %>%
  inner_join(fish.system, by = "Full.ID")

monitored <- windows %>%
  rowwise() %>%
  mutate(date = list(seq(start, end, by = "day"))) %>%
  ungroup() %>%
  unnest(date) %>%
  mutate(season = if_else(month(date) %in% spawn.months, "Spawning", "Non-spawning"))

days.monitored <- monitored %>%
  count(Full.ID, tag.system, season, name = "days_monitored")

# Days detected in the fish's home drainage (numerator) -----------------
daily <- flounder.est %>%
  inner_join(fish.system, by = "Full.ID") %>%
  inner_join(home.map, by = c("tag.system", "System")) %>%   # home-drainage detections only
  count(Full.ID, tag.system, date, name = "n_det") %>%
  filter(n_det >= min_detections_per_day) %>%
  mutate(season = if_else(month(date) %in% spawn.months, "Spawning", "Non-spawning"))

days.detected <- daily %>%
  count(Full.ID, tag.system, season, name = "days_detected")

# Residency index per fish x season ----------------------------------
ri.season <- days.monitored %>%
  left_join(days.detected, by = c("Full.ID", "tag.system", "season")) %>%
  mutate(
    days_detected = replace_na(days_detected, 0L),
    RI = days_detected / days_monitored,
    season = factor(season, levels = c("Non-spawning", "Spawning"))
  ) %>%
  filter(days_monitored >= min.days.monitored)

# The 6 runs: system x season ----------------------------------------
ri.summary <- ri.season %>%
  group_by(tag.system, season) %>%
  summarise(
    n_fish    = n_distinct(Full.ID),
    mean_RI   = mean(RI),
    sd_RI     = sd(RI),
    median_RI = median(RI),
    min_RI    = min(RI),
    max_RI    = max(RI),
    .groups = "drop"
  )
print(ri.summary)

#  Plot ---------------------------------------------------------------
p_ri <- ggplot(ri.season, aes(season, RI, fill = season)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  facet_wrap(~ tag.system) +
  ylim(0, 1) +
  labs(x = NULL, y = "Residency index",
       title = "Residency index inside vs outside spawning season (Nov-Feb)",
       subtitle = "Barrier islands excluded; diamonds = means") +
  theme_minimal() +
  theme(legend.position = "none")
print(p_ri)

# Tests --------------------------------------------------------------
# (a) Within each system only (no between-system comparisons):
#     paired Wilcoxon signed-rank on fish that have both seasons.
ri.wide <- ri.season %>%
  select(Full.ID, tag.system, season, RI) %>%
  pivot_wider(names_from = season, values_from = RI) %>%
  drop_na()

spawn.tests <- ri.wide %>%
  group_by(tag.system) %>%
  summarise(
    n_fish        = n(),
    median_non    = median(`Non-spawning`),
    median_spawn  = median(Spawning),
    p_value       = if (n() >= 5)
      wilcox.test(Spawning, `Non-spawning`, paired = TRUE, exact = FALSE)$p.value
    else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(
    p_holm    = p.adjust(p_value, method = "holm"),
    direction = if_else(median_spawn < median_non,
                        "RI dropped in spawning season",
                        "RI increased in spawning season")
  )
print(spawn.tests)

m.ri <- glmer(
  cbind(days_detected, days_monitored - days_detected) ~ season +
    (1 | tag.system) + (1 | Full.ID),
  data = ri.season, family = binomial,
  control = glmerControl(optimizer = "bobyqa")
)
summary(m.ri)

# Three-season version: pre-spawn, spawn, post-spawn (summer ignored)
pre.months  <- c(8, 9, 10)   
post.months <- c(2)    

label3 <- function(m) {
  case_when(
    m %in% pre.months   ~ "Pre-spawn",
    m %in% spawn.months ~ "Spawn",
    m %in% post.months  ~ "Post-spawn",
    TRUE                ~ NA_character_
  )
}

days.monitored3 <- monitored %>%
  mutate(season3 = label3(month(date))) %>%
  filter(!is.na(season3)) %>%
  count(Full.ID, tag.system, season3, name = "days_monitored")

days.detected3 <- daily %>%
  mutate(season3 = label3(month(date))) %>%
  filter(!is.na(season3)) %>%
  count(Full.ID, tag.system, season3, name = "days_detected")

ri.season3 <- days.monitored3 %>%
  left_join(days.detected3, by = c("Full.ID", "tag.system", "season3")) %>%
  mutate(
    days_detected = replace_na(days_detected, 0L),
    RI = days_detected / days_monitored,
    season3 = factor(season3, levels = c("Pre-spawn", "Spawn", "Post-spawn"))
  ) %>%
  filter(days_monitored >= min.days.monitored)

ri.summary3 <- ri.season3 %>%
  group_by(tag.system, season3) %>%
  summarise(
    n_fish    = n_distinct(Full.ID),
    mean_RI   = mean(RI),
    sd_RI     = sd(RI),
    median_RI = median(RI),
    .groups = "drop"
  )
print(ri.summary3)

p_ri3 <- ggplot(ri.season3, aes(season3, RI, fill = season3)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white") +
  facet_wrap(~ tag.system) +
  ylim(0, 1) +
  labs(x = NULL, y = "Residency index",
       title = "Residency index: pre-spawn, spawn, post-spawn",
       subtitle = "Pre = Sep-Oct, Spawn = Nov-Feb, Post = Mar-Apr; summer excluded; diamonds = means") +
  theme_minimal() +
  theme(legend.position = "none")
print(p_ri3)

# Within-system tests only: Kruskal-Wallis, then Holm-adjusted pairwise Wilcoxon
ri.season3 %>%
  group_split(tag.system) %>%
  lapply(function(d) {
    list(
      system   = unique(d$tag.system),
      kruskal  = kruskal.test(RI ~ season3, data = d),
      pairwise = pairwise.wilcox.test(d$RI, d$season3,
                                      p.adjust.method = "holm", exact = FALSE)
    )
  }) %>%
  print()

# Mean RI tables (printed and saved to CSV)
#2
means.2season <- ri.season %>%
  group_by(tag.system, season) %>%
  summarise(n_fish = n_distinct(Full.ID),
            mean_RI = round(mean(RI), 3),
            sd_RI   = round(sd(RI), 3),
            se_RI   = round(sd(RI) / sqrt(n()), 3),
            .groups = "drop")
cat("\n--- Mean RI: spawning vs non-spawning ---\n")
print(means.2season)

#3
means.3season <- ri.season3 %>%
  group_by(tag.system, season3) %>%
  summarise(n_fish = n_distinct(Full.ID),
            mean_RI = round(mean(RI), 3),
            sd_RI   = round(sd(RI), 3),
            se_RI   = round(sd(RI) / sqrt(n()), 3),
            .groups = "drop")
cat("\n--- Mean RI: pre-spawn / spawn / post-spawn ---\n")
print(means.3season)

# Compact wide versions (one row per drainage)
means.2season %>%
  select(tag.system, season, mean_RI) %>%
  pivot_wider(names_from = season, values_from = mean_RI) %>%
  print()
means.3season %>%
  select(tag.system, season3, mean_RI) %>%
  pivot_wider(names_from = season3, values_from = mean_RI) %>%
  print()

#write.csv(means.2season, "mean_RI_2season.csv", row.names = FALSE)
#write.csv(means.3season, "mean_RI_3season.csv", row.names = FALSE)

