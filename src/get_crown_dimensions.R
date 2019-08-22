get_crown_dimensions <- function(field_data, validate_allometry = F){
  require(lme4)
  # field_data <- field_data %>% select(individualID, stemDiameter, height, taxonID,siteID, domainID, 
  #                                     maxCrownDiameter, ninetyCrownDiameter) 
  train_set <- field_data %>% filter(height > 0) %>%
    filter(maxCrownDiameter > 0) %>% filter(ninetyCrownDiameter > 0)
  train_set[c("maxCrownDiameter","ninetyCrownDiameter", "stemDiameter", "height")] <- 
    log(train_set[c("maxCrownDiameter","ninetyCrownDiameter", "stemDiameter", "height")])
  
  #evaluate the model on a  stratified train-test splitting. no need to do it every time
  if(validate_allometry){
    
    train_data <- train_set %>% 
      group_by_(.dots=c("taxonID", "siteID")) %>%
      sample_frac(0.8)
    test_data <- train_set[!(train_set$individualID %in% train_data$individualID),]
    
    lmmax <- lmer((maxCrownDiameter) ~ (height) + ( 1 | taxonID) + (1 | siteID) + (1 | domainID), 
                  data = train_data, REML = FALSE)
    lmm90 <- lmer((ninetyCrownDiameter) ~ (height) + ( 1 | taxonID) + (1 | siteID) + (1 | domainID), 
                  data = train_data, REML = FALSE)
    
    maxCr <- predict(lmmax, newdata = test_data, allow.new.levels=TRUE)
    Cr90 <- predict(lmm90, newdata = test_data, allow.new.levels=TRUE)
    plot(y = exp(maxCr), x = exp(test_data$maxCrownDiameter))
    #plot(exp(Cr90), exp(test_data$ninetyCrownDiameter))
    
  }
  
  # build the allometric functions
  lmmax <- lmer((maxCrownDiameter) ~ (height) + ( 1 | taxonID) + (1 | siteID) + (1 | domainID), data = train_set,
              REML = FALSE)
  lmm90 <- lmer((ninetyCrownDiameter) ~ (height) + ( 1 | taxonID) + (1 | siteID) + (1 | domainID), data = train_set,
                REML = FALSE)
  unknown_size <- field_data %>% filter(!individualID %in% train_set$individualID)
  unknown_size[c("maxCrownDiameter", "ninetyCrownDiameter","stemDiameter", "height")] <- 
    log(unknown_size[c("maxCrownDiameter","ninetyCrownDiameter", "stemDiameter", "height")])
  
  maxCr <- predict(lmmax, newdata = unknown_size, allow.new.levels=TRUE) %>% exp
  Cr90 <- predict(lmm90, newdata = unknown_size, allow.new.levels=TRUE) %>% exp

  measured <- train_set %>% mutate(stemDiameter = exp(stemDiameter), height = exp(height),
                                   maxCrownDiameter = exp(maxCrownDiameter), 
                                   ninetyCrownDiameter = exp(ninetyCrownDiameter), crown_size = "measured")
  final_data <- unknown_size %>% mutate(stemDiameter = exp(stemDiameter), height = exp(height),
                                        maxCrownDiameter = maxCr, ninetyCrownDiameter = Cr90, crown_size = "height") %>%
    rbind.data.frame(measured)
  return(final_data)
}
