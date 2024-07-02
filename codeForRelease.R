
library(janitor)
library(raster)
library(RColorBrewer)
library(tidyverse)
library(readxl)
library(lme4)
library(AICcmodavg)


#------------------------------------------------
#Data import and processing
#------------------------------------------------

rm(list=ls())

survey = read.csv('dataReleaseSurvey.csv')
station = read.csv('dataReleaseStation.csv')
site = read.csv('dataReleaseSite.csv')
behavior = read.csv('dataReleaseVideos.csv')

#Combining the site, station, and survey data into a single data frame for regression analysis
analysisData = site %>% 
  left_join(station, by='SiteNum') %>% 
  left_join(survey, by=c('SiteNum', 'StationNum')) %>% 
  filter(!is.na(SurveyDate)) %>% #Gets rid of stations B & C for Site 18-26 because there are no surveys
  mutate(StationNum = paste(SiteNum, StationNum, sep='_'))

# Filter out the surveys conducted at nest sites outside of the active nesting period.
tmp1 = analysisData %>% 
  filter(Treatment=='nest') %>%
  filter(SurveyDate >= NestInitiationDate & SurveyDate <= FledgeFailDate)

# Filter out the control sites
tmp2 = analysisData %>% 
  filter(Treatment=='control')

# Put all surveys at control sites back together with those at nest sites during active nesting
# Also center and scale the distance from the stations to the central tree
analysisData = rbind(tmp1, tmp2) %>% 
  arrange(StationNum) %>% 
  mutate(scaledDistance = scale(Distance, center=T, scale=T)) %>% 
  mutate(scaledDistance2 = scaledDistance^2)

# Store the mean and sd of the distances of stations from the nest tree for later interpretation
distSd = sd(analysisData$Distance)
distMean = mean(analysisData$Distance)

#------------------------------------------------
#Data for Table 1
#------------------------------------------------

analysisData %>% 
  group_by(Treatment, Year) %>% 
  summarise(sites = length(unique(SiteNum)),
            stations = length(unique(StationNum)),
            surveys = n(),
            occupiedSurveys = sum(occupied),
            presenceSurveys = sum(presence))

#------------------------------------------------
#Model fitting and data for Table 2
#------------------------------------------------

#Fit 4 logistic regression models using occupancy as the response variable
occMod1 = glmer(occupied ~ (1|SiteNum) + (1|StationNum), data=analysisData, family=binomial(link='logit'))
occMod2 = glmer(occupied ~ Treatment +
                  (1|SiteNum) + (1|StationNum), data=analysisData, family=binomial(link='logit'))
occMod3 = glmer(occupied ~ Treatment + scaledDistance + Treatment*scaledDistance +
                  (1|SiteNum) + (1|StationNum), data=analysisData, family=binomial(link='logit'))
occMod4 = glmer(occupied ~ Treatment + scaledDistance + scaledDistance2 + 
                  Treatment*scaledDistance + Treatment*scaledDistance2 +
                  (1|SiteNum) + (1|StationNum), data=analysisData, family=binomial(link='logit'))

aictab(list('occ~1' = occMod1, 'occ~treat' = occMod2, 'occ~treat*dist' = occMod3, 'occ~treat*dist^2' = occMod4))


#Fit 4 logistic regression models using presence as the response variable
presMod1 = glmer(presence ~ (1|StationNum), data=analysisData, family=binomial(link='logit'))
presMod2 = glmer(presence ~ Treatment +
                   (1|StationNum), data=analysisData, family=binomial(link='logit'))
presMod3 = glmer(presence ~ Treatment + scaledDistance + Treatment*scaledDistance +
                   (1|StationNum), data=analysisData, family=binomial(link='logit'))
presMod4 = glmer(presence ~ Treatment + scaledDistance + scaledDistance2 + 
                   Treatment*scaledDistance + Treatment*scaledDistance2 +
                   (1|StationNum), data=analysisData, family=binomial(link='logit'))

aictab(list('pres~1' = presMod1, 'pres~treat' = presMod2, 'pres~treat*dist' = presMod3, 'pres~treat*dist^2' = presMod4))



#------------------------------------------------
#Testing the fit of the top occupancy model
#------------------------------------------------

#Set number of data simulations from the fitted model
nSims = 500

#Identify the fitted model we want to simulate from
testMod = occMod1

#Create a vector for storing chi square values from simulated datasets
simResults = rep(NA, nSims)

#Simulate from the fitted model
simData = simulate(testMod, nsim=nSims)

#Refit the occupancy model using each of the 500 simulated response data sets
#Calculate a chi-square statistic from each fitted model and record the value
for(i in 1:nSims){
  newAnalysisData = analysisData %>%
    mutate(occupied = simData[,i])
  suppressWarnings(suppressMessages({newModel = update(testMod, .~., data=newAnalysisData)}))
  o = newAnalysisData$occupied
  e = predict(newModel, type='response')
  simResults[i] = sum((o-e)^2/e, na.rm=T)
  rm(o, e)
}

#Create a dataframe out of the results from the simulation
simResults = data.frame(simulation = 1:nSims, chiSquare = simResults)

#Calculate the chi-square statistic from the model fit on the original data
modChisq = sum((analysisData$occupied - predict(testMod, type='response'))^2/predict(testMod, type='response'))

#Plotting the distribution of chiSquare values
ggplot(simResults, aes(x=chiSquare))+
  geom_histogram(fill='white', color='black')+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_vline(xintercept=modChisq, color='red')

#What proportion of the chi-square values from the simulated data are greater than
#the chi-square value from the model fit to the original data?
#Note this number might be slightly different than reported in the manuscript due
#to randomization in the simulations.
sum(simResults[,2] > modChisq)/length(simResults[,2]) # p = 0.642



#------------------------------------------------
#Testing the fit of the top presence model
#------------------------------------------------

#Set number of data simulations from the fitted model
nSims = 500

#Identify the fitted model we want to simulate from
testMod = presMod2

#Create a vector for storing chi square values from simulated datasets
simResults = rep(NA, nSims)

#Simulate from the fitted model
simData = simulate(testMod, nsim=nSims)

#Refit the occupancy model using each of the 500 simulated response data sets
#Calculate a chi-square statistic from each fitted model and record the value
for(i in 1:nSims){
  newAnalysisData = analysisData %>%
    mutate(presence = simData[,i])
  suppressWarnings(suppressMessages({newModel = update(testMod, .~., data=newAnalysisData)}))
  o = newAnalysisData$presence
  e = predict(newModel, type='response')
  simResults[i] = sum((o-e)^2/e, na.rm=T)
  rm(o, e)
}

#Create a dataframe out of the results from the simulation
simResults = data.frame(simulation = 1:nSims, chiSquare = simResults)

#Calculate the chi-square statistic from the model fit on the original data
modChisq = sum((analysisData$presence - predict(testMod, type='response'))^2/predict(testMod, type='response'))

#Plotting the distribution of chiSquare values
ggplot(simResults, aes(x=chiSquare))+
  geom_histogram(fill='white', color='black')+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_vline(xintercept=modChisq, color='red')

#What proportion of the chi-square values from the simulated data are greater than
#the chi-square value from the model fit to the original data?
#Note this number might be slightly different than reported in the manuscript due
#to randomization in the simulations.
sum(simResults[,2] > modChisq)/length(simResults[,2]) # p = 0.624


#-----------------------------------------------------------------------------
#Figure 2 - Examining detections during surveys of nest sites when a murrelet
#was known to have arrived and/or departed based on camera survey data
#-----------------------------------------------------------------------------

#Selecting and filtering data needed from the analysis datataset and creating
#date/time stamps for when the selected surveys started and ended for comparison
#with the observations on video
nestSurveys = analysisData %>% 
  filter(Treatment=='nest') %>% 
  select(SiteNum, StationNum, SurveyDate, TotDetections, SigBehavior, SurveyStart, SurveyEnd, presence, occupied) %>% 
  mutate(startTmp = as.POSIXct(paste(SurveyDate, "00:00:00"), format='%Y-%m-%d %H:%M:%S'),
         endTmp = as.POSIXct(paste(SurveyDate, "00:00:00"), format='%Y-%m-%d %H:%M:%S')) %>% 
  mutate(startTmp = startTmp + SurveyStart*60) %>% 
  mutate(endTmp = endTmp + SurveyEnd*60) %>% 
  mutate(SurveyStart = startTmp,
         SurveyEnd = endTmp) %>% 
  select(-startTmp, -endTmp)

#Tabulating how many arrival/departure events occurred during each survey for
#any survey for which at least one occurred
tmp = nestSurveys %>% 
  left_join(behavior, by='SiteNum') %>% 
  filter(eventTime > SurveyStart & eventTime < SurveyEnd) %>% 
  group_by(SiteNum, StationNum, SurveyDate) %>% 
  summarise(events = n()) %>% 
  ungroup()

#Merging the known arrival/departure events back in with the survey data for nest sites
nestSurveys = nestSurveys %>% 
  left_join(tmp, by=c('SiteNum', 'StationNum', 'SurveyDate')) %>% 
  mutate(events = ifelse(is.na(events), 0, events))

#Filter to only surveys when a murrelet was known to have arrived or departed from the nest
#Also creating distance bin categories for the surveys
surveysWithActivity = nestSurveys %>% 
  filter(events > 0) %>% 
  left_join(station %>% mutate(StationNum = paste(SiteNum, StationNum, sep="_")),
            by=c("SiteNum", "StationNum")) %>% 
  mutate(distBin = ifelse(Distance < 67, "0-67", "67-133")) %>% 
  mutate(distBin = ifelse(Distance > 133, "133-200", distBin))

#Summarizing the results from these surveys with activity and formatting
#for plotting
tmp = surveysWithActivity %>% 
  group_by(distBin) %>% 
  summarise(total = n(),
            occ = sum(occupied),
            propOcc = mean(occupied),
            pres = sum(presence),
            propPres = mean(presence)) %>% 
  ungroup() %>% 
  select(distBin, propOcc, propPres) %>% 
  pivot_longer(cols=c(propOcc, propPres)) %>% 
  mutate(name = ifelse(name=='propOcc', 'Occupancy', 'Presence')) %>% 
  mutate(name = factor(name, levels=c('Presence', 'Occupancy'))) %>% 
  mutate(distBin = factor(distBin, levels=c('0-67', '67-133', '133-200')))

#Create the plot
ggplot(tmp, aes(x=distBin, y=value, fill=name))+
  theme_bw()+
  theme(panel.grid=element_blank())+
  geom_bar(stat='identity',
           position=position_dodge(width=c(0.9)))+
  scale_fill_brewer(palette='Set2')+
  ylab('Proportion of surveys recorded')+
  xlab('Distance to nest tree (m)')+
  theme(legend.title=element_blank(),
        legend.position='top')+
  geom_text(aes(label=c('n = 39', '', 'n = 41', '', 'n = 29', ''),
                y=c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9)),
            size=3)


#------------------------------------------------------
#Figure 3: Predicted probabilities from the top models
#------------------------------------------------------

#Creating a new dataset from which to generate predicted values from the fitted model
#This allows  us to plot the expected occupancy and presence based on our fitted models
#for a range of plausible distance-from-focal-tree values
#Note that we picked an arbitrary station and site for the simulations.

newData = expand.grid('scaledDistance' = c(0, seq(min(analysisData$scaledDistance), max(analysisData$scaledDistance), 0.1)),
                      'Treatment' = c('nest', 'control'),
                      'StationNum' = '21-4_B',
                      'SiteNum' = '21-4') %>% 
  mutate(scaledDistance2 = scaledDistance^2) %>% 
  mutate(Distance = scaledDistance*distSd + distMean)

#Simulating 500 predicted values from the top occupancy model
#Takes 1-2 minutes to run
bootResultsOcc1 = bootMer(occMod1, FUN=function(x){predict(x, newData)}, nsim=500)

#Simulationg 500 predicted values from the top presence model
#Takes 1-2 minutes to run
bootResultsPres2 = bootMer(presMod2, FUN=function(x){predict(x, newData)}, nsim=500)


#Summarizing the results from the occupancy model predictions and formatting for plot
bootResults = bootResultsOcc1
tmpOcc = newData %>% 
  mutate('predicted' = plogis(bootResults$t0)) %>% 
  mutate('sd' = apply(bootResults$t, MARGIN=2, FUN=function(x){sd(x)})) %>% 
  mutate('lcl' = plogis(bootResults$t0 - sd),
         'ucl' = plogis(bootResults$t0 + sd)) %>% 
  mutate(Treatment = ifelse(Treatment=='nest', 'Nest', 'Occupied control')) %>% 
  mutate(Treatment = factor(Treatment, levels=c('Occupied control', 'Nest'))) %>% 
  # mutate('predFalseNeg' = 1-predicted,
  #        'uclFalseNeg' = 1-lcl,
  #        'lclFalseNeg' = 1-ucl) %>% 
  filter(Treatment == 'Nest') %>% 
  mutate(Treatment = as.character(Treatment)) %>% 
  mutate(Treatment = 'Occupancy, all sites') %>% 
  mutate(Detection = 'Occupancy')



#Summarizing the results from the presence model predictions and formatting for plot
bootResults = bootResultsPres2
tmpPres = newData %>% 
  mutate('predicted' = plogis(bootResults$t0)) %>% 
  mutate('sd' = apply(bootResults$t, MARGIN=2, FUN=function(x){sd(x)})) %>% 
  mutate('lcl' = plogis(bootResults$t0 - sd),
         'ucl' = plogis(bootResults$t0 + sd)) %>% 
  mutate(Treatment = ifelse(Treatment=='nest', 'Nest', 'Occupied control')) %>% 
  mutate(Treatment = factor(Treatment, levels=c('Occupied control', 'Nest'))) %>% 
  # mutate('predFalseNeg' = 1-predicted,
  #        'uclFalseNeg' = 1-lcl,
  #        'lclFalseNeg' = 1-ucl) %>% 
  mutate(Treatment = as.character(Treatment)) %>% 
  # mutate(Treatment = ifelse(Treatment=='Occupied control', 'Control', Treatment)) %>% 
  mutate(Treatment = ifelse(Treatment=='Occupied control', 'Presence, control sites', 'Presence, nest sites')) %>% 
  mutate(Detection = 'Presence')


#Binding the occupancy and presence results together for plot
forPlot = rbind(tmpOcc, tmpPres) %>% 
  mutate(Treatment = factor(Treatment, levels=c('Presence, nest sites', 'Presence, control sites', 'Occupancy, all sites')))

#Generating the plot
ggplot(forPlot, aes(x=Distance, y=predicted, color=Treatment, fill=Treatment, linetype=Treatment))+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2, color=NA)+
  geom_line(linewidth=2, alpha=0.6)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  # theme(legend.position=c(0.8, 0.8))+
  ylab('Predicted probability')+
  xlab('Meters from focal tree')+
  ylim(0,1)+
  scale_color_brewer(palette='Set2')+
  scale_fill_brewer(palette='Set2')+
  scale_linetype_manual(values=c('dashed', 'dashed', 'solid'))+
  theme(legend.title=element_blank())


            