# Function to calculate SMD

library(ggpubr)
SMD <- function(work3) {
  if (nrow(work3) > 0) {
    SMD_ML <- rma.uni(yi = work3$SMD,
                      vi = work3$SMDv,
                      method = "REML",
                      test = "t",
                      data = work3)
    return(SMD_ML)
  } else {
    return(NULL)  # Return NULL if there are no studies
  }
}

# Define a function to calculate weighted R-squared
calculate_weighted_r_squared <- function(response, predictor, weights) {
  model <- lm(response ~ predictor, weight = weights)
  return(sprintf("%.2f", summary(model)$r.squared)) }
  
# get the file from df and make a new grouping
working <- df %>%
  mutate(GroupID = interaction(StudyId_I, ExperimentID_I, CohortId_I))

#get the n of cohorts with each outcome measure
group_characteristics <- working %>%
  group_by(GroupID, outcome_type) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = outcome_type, values_from = count, values_fill = 0) %>%
  ungroup()

group_characteristics$total <- rowSums(select(group_characteristics, -GroupID))
group_characteristics$GroupID <- as.character(group_characteristics$GroupID)

#put them in a grid
l <- ncol(group_characteristics) - 2
m <- l +1
outcomes <- group_characteristics[, -c(1, m)]
out_grid <- outcomes %>%
  mutate(across(everything(), ~ . > 0, .names = "{.col}_TF"))
out_grid <- out_grid[, -c(1:l)]

row_counts <- rowSums(out_grid, na.rm = TRUE)
result <- group_characteristics[row_counts >= 2, ]
result <- result[,c(1,7,8,6,3,2,5,9,4)]
colu <- as.data.frame(colnames(result))

for (j in 2:(nrow(colu)-1)){ #1
  for (k in (j+1):(nrow(colu))){ #2
    outcome1 <- colu[j, 1]
    outcome2 <- colu[k, 1]
    
    if (!outcome1 == outcome2){ #3
      target <- subset(result, result[[outcome1]] > 0 & result[[outcome2]] > 0)
      
    if (nrow(target) > 0){ #4
        res <- data.frame()
        
        for (i in 1:nrow(target)) { 
          work3 <- working %>% filter(GroupID %in% target[i, 1])
          work3 <- work3 %>% filter(outcome_type == outcome1)
          work3 <- work3 %>% filter(!is.na(SMDv))
          work3$SMD <- as.numeric(work3$SMD)
          work3$SMDv <- as.numeric(work3$SMDv)
          
          SMDx <- SMD(work3)
          
          res[i, 1] <- if (!is.null(SMDx)) SMDx$beta[[1]] else NA
          res[i, 2] <- if (!is.null(SMDx)) SMDx$se[[1]] else NA
          
          work4 <- working %>% filter(GroupID %in% target[i, 1])
          work4 <- work4 %>% filter(outcome_type == outcome2)
          work4 <- work4 %>% filter(!is.na(SMDv))
          work4$SMD <- as.numeric(work4$SMD)
          work4$SMDv <- as.numeric(work4$SMDv)
          
          SMDy <- SMD(work4)
          cn <- which(names(work4) == 'SortLabel')
          
          res[i, 3] <- if (!is.null(SMDy)) SMDy$beta[[1]] else NA
          res[i, 4] <- if (!is.null(SMDy)) SMDy$se[[1]] else NA
          res[i, 5] <- if (!is.null(SMDy)) SMDy$k[[1]] else NA
          res[i, 6] <- if (!is.null(SMDy)) work4[1,cn] else NA
        }
        
        resm <- subset(res, res$SortLabel == 'CvS')
        resi <- subset(res, res$SortLabel == 'TvC')
        labelx <- bquote(Delta ~ .(as.name(outcome1)))
        labely <- bquote(Delta ~ .(as.name(outcome2)))

        # Calculate weighted R-squared for the model and intervention
        if(nrow(resm)>2){
          if(nrow(resi)>2){
          model_r_squared <- calculate_weighted_r_squared(resm[, 3], resm[, 1], 1/sqrt((resm[, 4]^2) + (resm[, 2]^2)))
          intervention_r_squared <- calculate_weighted_r_squared(resi[, 3], resi[, 1], 1/sqrt((resi[, 4]^2) + (resi[, 2]^2)))
          labelr <- bquote(R^2~(model) == .(model_r_squared) * "; " * R^2~(intervention) == .(intervention_r_squared))
          
          p <- ggplot() +
            geom_point(aes(x = resi[, 1], y = resi[, 3], color = resi[,6], size = resi[, 5]), na.rm = TRUE) +  # Use resi[, 5] as the size
            geom_errorbar(aes(x = resi[, 1], y = resi[, 3], ymin = resi[, 3] - resi[, 4], ymax = resi[, 3] + resi[, 4], width = 0.2), na.rm = TRUE) +
            geom_errorbar(aes(x = resi[, 1], y = resi[, 3], xmin = resi[, 1] - resi[, 2], xmax = resi[, 1] + resi[, 2], width = 0.2), na.rm = TRUE) +
            geom_smooth(method = "lm", se = TRUE, formula = y ~ x, 
                        aes(x = resi[, 1], y = resi[, 3], weight = 1/sqrt((resi[, 4]^2) + (resi[,2]^2))), 
                        color = "black", linetype = "dashed") +
            geom_point(aes(x = resm[, 1], y = resm[, 3], color = resm[,6], size = resm[, 5]), na.rm = TRUE) +  # Use resm[, 5] as the size
            geom_errorbar(aes(x = resm[, 1], y = resm[, 3], ymin = resm[, 3] - resm[, 4], ymax = resm[, 3] + resm[, 4], width = 0.2), na.rm = TRUE) +
            geom_errorbar(aes(x = resm[, 1], y = resm[, 3], xmin = resm[, 1] - resm[, 2], xmax = resm[, 1] + resm[, 2], width = 0.2), na.rm = TRUE) +
            geom_smooth(method = "lm", se = TRUE, formula = y ~ x, 
                        aes(x = resm[, 1], y = resm[,3], weight = 1/sqrt((resm[, 4]^2) + (resm[,2]^2))), color = "black", linetype = "dashed") +
            labs(x = labelx, y = labely, title = labelr) +
            coord_cartesian(clip = "off") +
            scale_size_continuous(name = 'Number of comparisons', breaks = ~unique(round(pretty(.)))) +  
            scale_color_manual(name = 'Experiment type', 
                               values = c("CvS" = "red", "TvC" = "green"),
                               labels = c('Model', 'Intervention')) +  # Set colors for SortLabel values
            expand_limits(y = 0) +
            expand_limits(x = 0)
          
            ggsave(paste0("plot_", j, "_", k, ".png"), p, width = 8, height = 6)
          
          } else {
              model_r_squared <- calculate_weighted_r_squared(resm[, 3], resm[, 1], 1/sqrt((resm[, 4]^2) + (resm[, 2]^2)))
              labelr <- bquote(R^2~(model) == .(model_r_squared) * ".")
              
              p <- ggplot() +
                geom_point(aes(x = resi[, 1], y = resi[, 3], color = resi[,6], size = resi[, 5]), na.rm = TRUE) +  # Use resi[, 5] as the size
                geom_errorbar(aes(x = resi[, 1], y = resi[, 3], ymin = resi[, 3] - resi[, 4], ymax = resi[, 3] + resi[, 4], width = 0.2), na.rm = TRUE) +
                geom_errorbar(aes(x = resi[, 1], y = resi[, 3], xmin = resi[, 1] - resi[, 2], xmax = resi[, 1] + resi[, 2], width = 0.2), na.rm = TRUE) +
                geom_point(aes(x = resm[, 1], y = resm[, 3], color = resm[,6], size = resm[, 5]), na.rm = TRUE) +  # Use resm[, 5] as the size
                geom_errorbar(aes(x = resm[, 1], y = resm[, 3], ymin = resm[, 3] - resm[, 4], ymax = resm[, 3] + resm[, 4], width = 0.2), na.rm = TRUE) +
                geom_errorbar(aes(x = resm[, 1], y = resm[, 3], xmin = resm[, 1] - resm[, 2], xmax = resm[, 1] + resm[, 2], width = 0.2), na.rm = TRUE) +
                geom_smooth(method = "lm", se = TRUE, formula = y ~ x, 
                            aes(x = resm[, 1], y = resm[, 3], weight = 1/sqrt((resm[, 4]^2) + (resm[,2]^2))), color = "black", linetype = "dashed") +
                labs(x = labelx, y = labely, title = labelr) +
                coord_cartesian(clip = "off") +
                scale_size_continuous(name = 'Number of comparisons', breaks = ~unique(round(pretty(.)))) +  
                scale_color_manual(name = 'Experiment type', 
                                   values = c("CvS" = "red", "TvC" = "green"),
                                   labels = c('Model', 'Intervention')) +  # Set colors for SortLabel values
                expand_limits(y = 0) +
                expand_limits(x = 0)
              
                ggsave(paste0("plot_", j, "_", k, ".png"), p, width = 8, height = 6)
              }
          } else {
          if(nrow(resi)>2){
                intervention_r_squared <- calculate_weighted_r_squared(resi[, 3], resi[, 1], 1/sqrt((resi[, 4]^2) + (resi[, 2]^2)))
                labelr <- bquote(R^2~(intervention) == .(intervention_r_squared))
                
                p <- ggplot() +
                  geom_point(aes(x = resi[, 1], y = resi[, 3], color = resi[,6], size = resi[, 5]), na.rm = TRUE) +  # Use resi[, 5] as the size
                  geom_errorbar(aes(x = resi[, 1], y = resi[, 3], ymin = resi[, 3] - resi[, 4], ymax = resi[, 3] + resi[, 4], width = 0.2), na.rm = TRUE) +
                  geom_errorbar(aes(x = resi[, 1], y = resi[, 3], xmin = resi[, 1] - resi[, 2], xmax = resi[, 1] + resi[, 2], width = 0.2), na.rm = TRUE) +
                  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, 
                              aes(x = resi[, 1], y = resi[, 3], weight = 1/sqrt((resi[, 4]^2) + (resi[,2]^2))), 
                              color = "black", linetype = "dashed") +
                  geom_point(aes(x = resm[, 1], y = resm[, 3], color = resm[,6], size = resm[, 5]), na.rm = TRUE) +  # Use resm[, 5] as the size
                  geom_errorbar(aes(x = resm[, 1], y = resm[, 3], ymin = resm[, 3] - resm[, 4], ymax = resm[, 3] + resm[, 4], width = 0.2), na.rm = TRUE) +
                  geom_errorbar(aes(x = resm[, 1], y = resm[, 3], xmin = resm[, 1] - resm[, 2], xmax = resm[, 1] + resm[, 2], width = 0.2), na.rm = TRUE) +
                  labs(x = labelx, y = labely, title = labelr) +
                  coord_cartesian(clip = "off") +
                  scale_size_continuous(name = 'Number of comparisons', breaks = ~unique(round(pretty(.)))) +  
                  scale_color_manual(name = 'Experiment type', 
                                     values = c("CvS" = "red", "TvC" = "green"),
                                     labels = c('Model', 'Intervention')) +  # Set colors for SortLabel values
                  expand_limits(y = 0) +
                  expand_limits(x = 0)
                
                ggsave(paste0("plot_", j, "_", k, ".png"), p, width = 8, height = 6)
                
            } else {
                p <- ggplot() +
                    geom_point(aes(x = resi[, 1], y = resi[, 3], color = resi[,6], size = resi[, 5]), na.rm = TRUE) +  # Use resi[, 5] as the size
                    geom_errorbar(aes(x = resi[, 1], y = resi[, 3], ymin = resi[, 3] - resi[, 4], ymax = resi[, 3] + resi[, 4], width = 0.2), na.rm = TRUE) +
                    geom_errorbar(aes(x = resi[, 1], y = resi[, 3], xmin = resi[, 1] - resi[, 2], xmax = resi[, 1] + resi[, 2], width = 0.2), na.rm = TRUE) +
                    geom_point(aes(x = resm[, 1], y = resm[, 3], color = resm[,6], size = resm[, 5]), na.rm = TRUE) +  # Use resm[, 5] as the size
                    geom_errorbar(aes(x = resm[, 1], y = resm[, 3], ymin = resm[, 3] - resm[, 4], ymax = resm[, 3] + resm[, 4], width = 0.2), na.rm = TRUE) +
                    geom_errorbar(aes(x = resm[, 1], y = resm[, 3], xmin = resm[, 1] - resm[, 2], xmax = resm[, 1] + resm[, 2], width = 0.2), na.rm = TRUE) +
                    labs(x = labelx, y = labely) +
                    coord_cartesian(clip = "off") +
                    scale_size_continuous(name = 'Number of comparisons', breaks = ~unique(round(pretty(.)))) +  
                    scale_color_manual(name = 'Experiment type', 
                                       values = c("CvS" = "red", "TvC" = "green"),
                                       labels = c('Model', 'Intervention')) +  # Set colors for SortLabel values
                    expand_limits(y = 0) +
                    expand_limits(x = 0)
                  
                  ggsave(paste0("plot_", j, "_", k, ".png"), p, width = 8, height = 6)
      }
    }
  
          }}}}
    
#now take all neuro v all behavioural
  group_characteristics <- working %>%
  group_by(GroupID, out_cat) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = out_cat, values_from = count, values_fill = 0) %>%
  ungroup()

group_characteristics$total <- rowSums(select(group_characteristics, -GroupID))
group_characteristics$GroupID <- as.character(group_characteristics$GroupID)

#put them in a grid
l <- ncol(group_characteristics) - 2
m <- l +1
outcomes <- group_characteristics[, -c(1, m)]
out_grid <- outcomes %>%
  mutate(across(everything(), ~ . > 0, .names = "{.col}_TF"))
out_grid <- out_grid[, -c(1:l)]

row_counts <- rowSums(out_grid, na.rm = TRUE)
result <- group_characteristics[row_counts >= 2, ]
result <- result[,c(1,3,2)]
colu <- as.data.frame(colnames(result))

outcome1 <- 'Behavioural'
outcome2 <- 'Neurobiological'

target <- subset(result, result[[outcome1]] > 0 & result[[outcome2]] > 0)

res <- data.frame()

for (i in 1:nrow(target)) {
  work3 <- working %>% filter(GroupID %in% target[i, 1])
  work3 <- work3 %>% filter(out_cat == outcome1)
  work3 <- work3 %>% filter(!is.na(SMDv))
  work3$SMD <- as.numeric(work3$SMD)
  work3$SMDv <- as.numeric(work3$SMDv)
  
  SMDx <- SMD(work3)
  
  res[i, 1] <- if (!is.null(SMDx)) SMDx$beta[[1]] else NA
  res[i, 2] <- if (!is.null(SMDx)) SMDx$se[[1]] else NA
  
  work4 <- working %>% filter(GroupID %in% target[i, 1])
  work4 <- work4 %>% filter(out_cat == outcome2)
  work4 <- work4 %>% filter(!is.na(SMDv))
  work4$SMD <- as.numeric(work4$SMD)
  work4$SMDv <- as.numeric(work4$SMDv)
  
  SMDy <- SMD(work4)
  cn <- which(names(work4) == 'SortLabel')
  
  res[i, 3] <- if (!is.null(SMDy)) SMDy$beta[[1]] else NA
  res[i, 4] <- if (!is.null(SMDy)) SMDy$se[[1]] else NA
  res[i, 5] <- if (!is.null(SMDy)) SMDy$k[[1]] else NA
  res[i, 6] <- if (!is.null(SMDy)) work4[1,cn] else NA
}

resm <- subset(res, res$SortLabel == 'CvS')
resi <- subset(res, res$SortLabel == 'TvC')

# Calculate weighted R-squared for the model and intervention
model_r_squared <- calculate_weighted_r_squared(resm[, 3], resm[, 1], 1/sqrt((resm[, 4]^2) + (resm[, 2]^2)))
intervention_r_squared <- calculate_weighted_r_squared(resi[, 3], resi[, 1], 1/sqrt((resi[, 4]^2) + (resi[, 2]^2)))

# Create the label with superscript 2 using bquote
labelr <- bquote(R^2~(model) == .(model_r_squared) * "; " * R^2~(intervention) == .(intervention_r_squared))
labelx <- 'All Behavioural outcomes'
labely <- 'All Neurobiological outcomes'

p <- ggplot() +
  geom_point(aes(x = resi[, 1], y = resi[, 3], color = resi[,6], size = resi[, 5]), na.rm = TRUE) +  # Use resi[, 5] as the size
  geom_errorbar(aes(x = resi[, 1], y = resi[, 3], ymin = resi[, 3] - resi[, 4], ymax = resi[, 3] + resi[, 4], width = 0.2), na.rm = TRUE) +
  geom_errorbar(aes(x = resi[, 1], y = resi[, 3], xmin = resi[, 1] - resi[, 2], xmax = resi[, 1] + resi[, 2], width = 0.2), na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, 
              aes(x = resi[, 1], y = resi[, 3], weight = 1/sqrt((resi[, 4]^2) + (resi[,2]^2))), 
              color = "black", linetype = "dashed") +
  geom_point(aes(x = resm[, 1], y = resm[, 3], color = resm[,6], size = resm[, 5]), na.rm = TRUE) +  # Use resm[, 5] as the size
  geom_errorbar(aes(x = resm[, 1], y = resm[, 3], ymin = resm[, 3] - resm[, 4], ymax = resm[, 3] + resm[, 4], width = 0.2), na.rm = TRUE) +
  geom_errorbar(aes(x = resm[, 1], y = resm[, 3], xmin = resm[, 1] - resm[, 2], xmax = resm[, 1] + resm[, 2], width = 0.2), na.rm = TRUE) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ x, 
              aes(x = resm[, 1], y = resm[, 3], weight = 1/sqrt((resm[, 4]^2) + (resm[,2]^2))), color = "black", linetype = "dashed") +
  labs(x = labelx, y = labely, title = labelr) +
  coord_cartesian(clip = "off") +
  scale_size_continuous(name = 'Number of comparisons', breaks = ~unique(round(pretty(.)))) +  
  scale_color_manual(name = 'Experiment type', 
                     values = c("CvS" = "red", "TvC" = "green"),
                     labels = c('Model', 'Intervention')) +  # Set colors for SortLabel values
  expand_limits(y = 0) +
  expand_limits(x = 0)

  ggsave(paste0("BvN.png"), p, width = 8, height = 6)
