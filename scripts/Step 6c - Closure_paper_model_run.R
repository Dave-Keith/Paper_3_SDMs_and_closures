# THIS SCRIPT RUNS THE MODELS USING DATA UP TO 2019 rather than just to 2016. 

# Load er up...

rm(list=ls())
#direct.fun <- "Y:/Offshore/Assessment/"
direct.fun <- "D:/Github/Assessment_fns/"; dir.ftmp <- direct.fun
#direct.proj <- "Y:/Projects/GB_time_area_closure_SPERA/"
direct.proj <- "D:/Github/Paper_3_SDMs_and_closures/"; dir.tmp <- direct.proj

library(INLA)
library(boot)
library(fields)
library(PBSmapping)
require(tidyverse)
require(reshape2)
require(GGally)
require(boot)
require(bbmle)
require(MASS)
require(leaps)
require(COUNT)
require(cowplot)
require(viridis)
require(maptools)
library(maps)
require(gridExtra)
library(sp)
library(rgeos)
library(splancs)
library(marmap)
library(dplyr)
library(spatstat)
library(gridBase)
library(animation)
library(sf)
library(raster)
library(rgdal)
library(maptools)
library(mapdata)
#library(sdmTMB)
library(caret)
library(data.table)
library(units)
library(concaveman)
#library(SDMTools)
#devtools::install_github("pbs-assess/sdmTMB") # You may need to install the sdmTMB package from Sean Anderson
#inla.upgrade(testing = F)
# Set the number of threads INLA uses, this is done automatically and set the number of threads to number of logical procesors a computer has
# But this doesn't really seem to push the computer all that hard as a default...  Not clear if increasing this would be useful or not
# So the control option inside INLA for another way to tweak how many resources are dedicated to the model...
#inla.setOption( num.threads = 12) 
# Bring in the functions we'll need
source(paste(direct.fun,"Maps/pectinid_projector_sf.R",sep=""))
source(paste(direct.fun,"Maps/convert_coords.R",sep=""))
source(paste(direct.fun,"Maps/add_alpha_function.r",sep=""))
source(paste(direct.fun,"Maps/combo_shp.R",sep=""))
source(paste(direct.fun,"Maps/centre_of_gravity.R",sep=""))
# A function for calculating the cpo for a model...
fcpo <- function(m, id)
  -sum(log(m$cpo$cpo[id]), na.rm=TRUE) # Good to log this because of distributional reasons.
factor.2.number <- function(x) {as.numeric(levels(x))[x]}

# Here is the data we need, this comes from Step 3 INLA_mesh_for_gb_surveys_and_scallop_survey.R
load(paste0(direct.proj,"Data/INLA_mesh_input_data.RData"))
load(paste0(direct.proj,"Data/INLA_meshes.RData"))
load(paste0(direct.proj,"Data/Prediction_mesh.Rdata"))
#load(paste0(direct.proj,"Data/SST_and_Depth_covariates_and_boundary_for_prediction.RData"))
# This overwrites the above, need to tidy this all up...
load(paste0(direct.proj,"Data/Depth_SST_and_Sed_on_GB.RData"))
load(paste0(direct.proj,"Data/Survey_data_with_ALL_covars_2017_2020.RData"))
load(paste0(direct.proj,"Data/Depth_SST_and_Sed.RData"))


#load(paste0(direct.proj,"Results/INLA_st_3_output.RData"))
direct.proj <- dir.tmp 
direct.fun <- dir.ftmp 

# Now I need to stick my UTM coordinates into the data
dat.final$X <- loc.gf@coords[,1]
dat.final$Y <- loc.gf@coords[,2]

# Remove the nmfs fall data here as I won't run these additional diagonstics/validations on these data...
#dat.val <- dat.final %>% dplyr::filter(survey != "nmfs-fall") 
# At this point I really start the analysis, and what I want to do is split up each survey, I might as well loop this, or at least allow for a loop
# so we can move through each survey easily, this will get complicated with a number of models for each survey, especially when I split out the time...
# The  INLA parameters
sigma <- 0.5
s.alpha <- 0.05
r.alpha <- 0.05
# I need a loop for the number of "species" i'm looking at, Now all I care about are Cod and Yellowtail, everything 
# else is taken care of in Step 4.
species <- c("cod_PA","yt_PA") 
#species <- c("cod_PA_can","yt_PA_can") # These take the cod and YT data that are on the Canadian side of GB. 
num.species <- length(species)

# This predicts for all 6 of our final models, not including the intercept model because no need....

# Need to initialize some objects for later
st.mods <- c(5,3)
num.fields <- length(st.mods)
num.mods <- 1 # The number of models I'm gonna run validation on


dat.all <- dat.final %>% dplyr::select("SEDNUM","comldepth","lat","lon","survey","year","unique_set_ID","cod_PA","yt_PA","X","Y","sst_avg","years_3","years_5")

# Now I'm being lazy, but the new.dat.final doesn't have the sediment field in it, need to put those together

new.dat.final <- new.dat.final %>% dplyr::select("SEDNUM","comldepth","lat","lon","survey","year","unique_set_ID","cod_PA","yt_PA","X","Y","sst_avg","years_3","years_5")
# Need to save this new data file will everything in it...
# Order this by year here...
new.dat.final <- new.dat.final[order(new.dat.final$year),]
new.dat.final$cod_PA <- as.integer(new.dat.final$cod_PA)
new.dat.final$yt_PA <- as.integer(new.dat.final$yt_PA)
new.dat.final$SEDNUM <- as.integer(new.dat.final$SEDNUM)
#new.dat.final$comldepth <- -new.dat.final$comldepth
st_geometry(new.dat.final) <- NULL
dat.all.final <- rbind(dat.all,new.dat.final)
#saveRDS(dat.all.final,paste0(direct.proj,"/Data/survey_data_1970_2019.rds"))

# I want to set the data up here so we can get the folds identified now
# dat.cod <-  dat.all %>% dplyr::filter(survey == "RV")
# dat.yt <- dat.all %>% dplyr::filter(survey == "nmfs-spring")
# Number of folds, going for 5 fold cross validation
surveys <- unique(dat.all.final$survey)
num.surveys <- length(surveys)

# I need to redo the years 5 and year 3 since we've added 3 new years. For the Years 3 it's easy since it's just a new number
#dat.all$years_3[dat.all$year %in% c(2017:2019)] <- dat.all$years_3[dat.all$year %in% c(2017:2019)] + 1
# The year 5 is much more painful because it isn't even...
eras.5 <- c(seq(max(dat.all.final$year),min(dat.all.final$year), by = -5),min(dat.all.final$year))
n.eras.5 <- length(eras.5)
for(i in 1:(n.eras.5-1)) dat.all.final$years_5[dat.all.final$year <= eras.5[i] & dat.all.final$year > eras.5[i+1]] <- n.eras.5-i

eras.3 <- c(seq(max(dat.all.final$year),min(dat.all.final$year), by = -3),min(dat.all.final$year))
n.eras.3 <- length(eras.3)
for(i in 1:(n.eras.3-1)) dat.all.final$years_3[dat.all.final$year <= eras.3[i] & dat.all.final$year > eras.3[i+1]] <- n.eras.3-i



# Now we need to do the prediction fun...
############################# SECTION 6 ################################ SECTION 6  ###################################################################
############################# SECTION 6 ################################ SECTION 6  ###################################################################
############################# SECTION 6 ################################ SECTION 6  ###################################################################
############################# SECTION 6 ################################ SECTION 6  ###################################################################
# This used to be final section is MODEL PREDICTIONS, which is pretty simple, we take the best Cod/YT models
# for each survey and and we add a prediction stack to them which will predict at the vertex of
# every mesh node using the SST and depth information along with the random field.
# Problem is this blows up the A matrix and really grinds down the computation so
# we need to really trim our prediction box down to the real core part of GBa even more
# than our nice 'clp' poly...
# Here is a nicely trimmed up polygon we can use to clip down our clipped area...
#loc.gf.sf <- st_as_sf(loc.gf)
clp.poly <- st_as_sf(data.frame(X = c(508000,508000,900000,650000,600000,550000),
                                Y=c(4540000,4350000,4674000,4674000,4661000,4622000),ID=1),coords = c("X","Y"),crs= 32619)
# Now make this a polygon
clp.poly <- st_cast(st_combine(clp.poly),"POLYGON")
clp.pred <- st_intersection(clp,clp.poly)
# How does this compare to the GB grid size, these are about 44km2 so about 6.5 by 6.5..
# GB.grid <- read.csv("D:/NAS/Projects/GB_time_area_closure_SPERA/Results/GB_grid.csv")
# GB.grid <- as.PolySet(GB.grid,projection = "LL")
# GB.grid <- st_as_sf(maptools::PolySet2SpatialPolygons(GB.grid),crs = "4326")
# GB.grid <- st_transform(GB.grid,crs=32619)
# st_area(GB.grid[20,])/1e6 # around 44 km^2 so almost 7 km ^2.  I'd like 
# Now I can make the prediction grid, the n is the number of cells in the x and y direction
# 42 is 1200 boxes and a 6 x 6 resolution. 51 is 1800 points and a 5 x 5 resolution 
#63 is 2700 points and a 4 x 4 resolution, 84 is 4700 points and a 3 x 3 resolution,
# 127 is 11000 points and a 2 x 2 resolution, and 
# 200 gives boxes that are around 2 km by 2 km
mesh.grid <- st_make_grid(clp.pred,n = 84) # So we'd get almost 5 cells inside every 1 of the closure grids, that's decent...
st_area(mesh.grid[1,])/1e6
#ggplot(mesh.grid)+ geom_sf()
# Now I just want to predict in the middle of the mesh grid..
mesh.points <- st_centroid(mesh.grid)
# Both the mesh grid and the mesh points are very useful objects I want, so let's save them out
#save(mesh.grid,mesh.points,file = "D:/NAS/Projects/GB_time_area_closure_SPERA/Results/Prediction_mesh.RData")

# A plot of all this fun stuff
ggplot(clp) + geom_sf() + 
  #geom_sf(data=loc.gf.sf)  + 
  #geom_sf(data = clp.poly,fill=NA,colour='blue') + 
  geom_sf(data = clp.pred,fill=NA,colour='green',size=3) + 
  #geom_sf(data=mesh.grid)+ 
  #geom_sf(data=GB.grid) +
  geom_sf(data=mesh.points,size=1) + 
  coord_sf(datum=st_crs(32619))


# Now all I have to do is intersect the clip polying with the sst and depth fields and we got ourselves our covariates
pred.covars <- st_intersection(sst.gb,mesh.points)
names(pred.covars)[1] <- "sst_avg"
pred.covars <- st_intersection(depth.gb,pred.covars)
names(pred.covars)[1] <- 'comldepth'
sed.tmp <- sed.gb %>% dplyr::select(sednum,geometry)
sed.tmp$sednum <- as.character(sed.tmp$sednum)
sed.tmp$sednum[sed.tmp$sednum == 5] <- NA
pred.covars <- st_intersection(sed.tmp,pred.covars)
names(pred.covars)[1] <- 'SEDNUM'

pred.covars$X <- st_coordinates(pred.covars)[,1]
pred.covars$Y <- st_coordinates(pred.covars)[,2]
st_geometry(pred.covars) <- NULL
# If I have any positive depths toss those...
pred.covars$comldepth[pred.covars$comldepth >0] <- 0
# Need to initialize some objects for later
st.mods <- c(5,3)
num.fields <- length(st.mods)

surveys <- unique(dat.all.final$survey)
num.surveys <- length(surveys)
# We are going to run the 5 year model for cod and the 3 year model for yellowtail
# This will have everything ever.  Note the best model is always the depth-sst additive model, so that's
# at least handy...
mesh <- mesh.gf
#res.fd <- NULL
#mod.diag <- NULL
mod.output.pred <- NULL
mod.diagnostics.pred <- NULL
res.pred <- NULL
rand.field.pred <- NULL
pred.output.pred <- NULL

#pred.res.st.yr <- NULL

for(s in 1:num.species) 
{
  # Loop through each of the surveys
  for(i in 1:num.surveys) 
  {
    
    # Now lets get our input data sorted
    # Let's select the data for the particular survey of interest
    dat <- dat.all.final[dat.all.final$survey == surveys[i],]
    pred.c <- pred.covars
    # Rename the varialbe of interest to "response"
    if(species[s] == "cod_PA")
    {
      resp <- "cod_PA"
      #pred.c <- pred.c %>% dplyr::select(-SEDNUM)
    }
    if(species[s] == "yt_PA") resp <- "yt_PA"
    response <- which(names(dat) == resp)
    names(dat)[response] <- "response"
    # Next I need to add the era information to the pred data, which means I need
    # a full copy of it for every random field in the model, gross, yes!
    # if(species[s] == "cod_PA") 
    # {
    if(species[s] == 'yt_PA' && surveys[i] != 'nmfs-fall') 
    {
      eras.p <- unique(dat$years_3)
      num.eras.p <- length(eras.p)
      tmp <- NULL
      for(p in 1:num.eras.p) tmp[[p]] <- data.frame(pred.c,years_3 = eras.p[p],year=NA)
      covar.pred <- do.call('rbind',tmp)
      covar.pred$response <- NA
      dat <- dat %>% dplyr::select(comldepth,sst_avg,SEDNUM,X,Y,years_3,response,year)
      dat <- rbind(dat,covar.pred)
      dat$years_5 <- NA
      # Turn the Sediment numbers into a factor
      dat$SEDNUM[dat$SEDNUM %in% c(2,5,6,8,9)] <- NA
      dat$fSEDNUM <- as.factor(dat$SEDNUM)
      dat$years_3 <- as.factor(dat$years_3)
      eras <- as.numeric(dat$years_3)
    } else
    {
      eras.p <- unique(dat$years_5)
      num.eras.p <- length(eras.p)
      tmp <- NULL
      for(p in 1:num.eras.p) tmp[[p]] <- data.frame(pred.c,years_5 = eras.p[p],year=NA)
      covar.pred <- do.call('rbind',tmp)
      covar.pred$response <- NA
      dat <- dat %>% dplyr::select(SEDNUM,comldepth,sst_avg,X,Y,years_5,response,year)
      dat <- rbind(dat,covar.pred)
      dat$years_5 <- as.factor(dat$years_5)
      eras <- as.numeric(dat$years_5)
      dat$years_3 <- NA
      # Turn the Sediment numbers into a factor
      dat$SEDNUM[dat$SEDNUM %in% c(2,5,6,8,9)] <- NA
      dat$fSEDNUM <- as.factor(dat$SEDNUM)
    } #end the else 
    
    
    # Lets log transform depth and sst_avg, 
    dat$depth_log <- log(-dat$comldepth)
    dat$depth_cen <-  dat$depth_log - mean(dat$depth_log,na.rm=T) # Log transform should help with issues related to skew of the depth data.
    dat$sst_avg_cen <- as.vector(scale(dat$sst_avg))
    

    # Get the location of our data...
    loc <- cbind(dat$X,dat$Y)
    
    #We also need to decide Observation Likelihood
    fam <- "binomial"
    
    # The amount of data we have
    N = nrow(dat)
    # For both the beta and binomial families we'll need to determine the number of trials.
    Ntrials <- 1 # For each record there is only 1 trial.
    
    # For both our scenarios we are going to be using the logit model (note that this isn't scrictly necessary to write as the logit is the
    # 'canonical' link (to likely mis-use stats terminology) for the beta and binomial distributions.
    control.fam = list(control.link=list(model="logit"))
    # Now make the A matrix for the the model with the spatio-temporal random field we need a different A matrix.
    # First the 10 year era
    # Now the 3 and 5 year model
    #if(species[s] == "cod_PA") 
    
    #if(species[s] == "yt_PA") eras <- as.numeric(dat$years_3)
    
    era.names <- unique(eras)
    n.eras <- length(unique(eras))
    #if(species[s] == "cod_PA") 
    
    # Only do the commented bits if I'm interested in an AR1 model...
    # if(species[s] == 'yt_PA' && surveys[i] != 'nmfs-fall') 
    # {
      # A.era <- inla.spde.make.A(mesh, loc,group = eras,n.groups =n.eras) # 3 year models
    # } else {
    A.era <- inla.spde.make.A(mesh, loc,repl = eras)
    # } # 5 year models
    #if(species[s] == "yt_PA") A.era <- inla.spde.make.A(mesh, loc,group = eras,n.groups =n.eras)
    
    # While I have my range for my spatial priors I don't have my sigma or the probabilites for the prirors
    # The priors here can be pretty informative, the SPDE approximation improves if the corrleation length (i.e. roughly the range) of the process
    # is similar, but larger, than the maximum edge length ofthe mesh
    # so here we define our spde
    range <- range.gf
    spde <- inla.spde2.pcmatern(mesh,    
                                prior.sigma=c(sigma,s.alpha), # The probabiliy that the marginal standard deviation (first number) is larger than second number
                                prior.range=c(range,r.alpha)) # The Meidan range and the probability that the range is less than this...
    
    # and now we define the spatio-temporal random field.  
    # If I want to use group of the 3 year era so that I can make the temporal field be an AR1 process I need to do the bits that are commented out.
    # 
    #if(species[s] == "cod_PA") 
    #if(species[s] == 'yt_PA' && surveys[i] != 'nmfs-fall') 
    #{
    # w.index <- inla.spde.make.index(name = 'w',n.spde = spde$n.spde,n.group = n.eras)
    #} else {
    w.index <- inla.spde.make.index(name = 'w',n.spde = spde$n.spde,n.rep = n.eras)
    #}
    #if(species[s] == "yt_PA") w.index <- inla.spde.make.index(name = 'w',n.spde = spde$n.spde,n.group = n.eras)
    # Zuur never talks about this puppy I don't think, it is a penalised complexity prior but I'm not sure what for, Zuur only
    # discusses these in terms of the PCP's of the spatial field, this is a prior for precision, see inla.doc("pc.prec")
    # certainly isn't entirely clear to me!
    #pcprec <- list(prior='pc.prec', param=c(0.5, 0.01))
    
    
    # For all random walk models I use inla.group to group the data so that the minimum difference
    # between values is greater than the threshold for the random walk models, these groups are very
    # fine and really just bin things that are essentially the same/
    # First up we only need to do this with the variables that appear interesting from
    # the spatial model runs.  These likely differ by species so might need a YT and a Cod list
    # guess it depends if any of them come out from the modelling as significant?
    dat$depth_cen_g <- inla.group(dat$depth_cen,n=75)
    dat$sst_avg_cen_g <- inla.group(dat$sst_avg_cen,n=75)
    
    #  da <- dat #%>% dplyr::filter(year <2016)
    #  da$depth_cen_g <- inla.group(da$depth_cen,n=100)
    #  da$sst_avg_cen_g <- inla.group(da$sst_avg_cen,n=100)
    # # 
    #  ggplot(da) + geom_histogram(aes(y=comldepth)) + facet_wrap(~years_3)
    # ggplot(da) + geom_histogram(aes(y=sst_avg_cen_g)) + facet_wrap(~years_3)
    # ggplot(da) + geom_histogram(aes(y=factor.2.number(fSEDNUM))) + facet_wrap(~years_3)
    # First we need to make the model SEDNUM# First we need to make the model matrix, this needs to align with our formula below, this doesn't strictly need to be 
    # done unless we have a categorical covariate
    options(na.action='na.pass')# Need to do this so that the model matrix retains the NA's in it.
    # The nice thing here is that we can make this a complex as we want and just run submodels from within this model
    
    if(species[s] == "cod_PA")
    {
      X.matrix <- model.matrix(~ 0+ depth_cen_g +  sst_avg_cen_g , data = dat)
      
      # And then make a covariate matrix
      X <- data.frame(depth =        X.matrix[,1],
                      sst   =        X.matrix[,2])
    }
    
    if(species[s] == "yt_PA")
    {
      X.matrix <- model.matrix(~ 0+ depth_cen_g +  sst_avg_cen_g + fSEDNUM, data = dat)
      
      # Our matrix for the yt models...
      X <- data.frame(depth =        X.matrix[,1],
                      sst   =        X.matrix[,2],
                      fsed_3     =   X.matrix[,3],
                      fsed_4     =   X.matrix[,4])
    }
    
    # Next up we make the stack...
    # I need a stack, probably should eventually split this into an estimation, validation and prediction stack, but for now
    # will stick with the one estimation stack....
    # Make the stack for the spatial models without spatio-temporal correlation
    stk.e = inla.stack(tag="est",
                       data=list(y = dat$response[!is.na(dat$response)], link=1L),
                       effects=list(intercept = rep(1, length(dat$response[!is.na(dat$response)])), 
                                    X = X[!is.na(dat$response),],
                                    w = w.index),
                       A = list(1,1,A.era[!is.na(dat$response),]))
    
    stk.p = inla.stack(tag="pred",
                       data=list(y = dat$response[is.na(dat$response)], link=1L),
                       effects=list(intercept = rep(1, length(dat$response[is.na(dat$response)])), 
                                    X = X[is.na(dat$response),],
                                    w = w.index),
                       A = list(1,1,A.era[is.na(dat$response),]))
    
    
    #### join data stacks and extracts the data index
    stk <- inla.stack(stk.e, stk.p)
    e.id <- inla.stack.index(stk, 'est')$data
    p.id <- inla.stack.index(stk, 'pred')$data
    # Now let's make our formula, intercept with each observation set as a random variable, this makes 
    # As soon as you make a spatial model make your own intercept.  Here is an initial suite of models I like...
    intercept <- 1 # intercept
    # For the random walk models we need to set the priors for the random walk, Zuur recommmends rw2 as it seems to 
    # overfit less, to do this I need to bin the covariates using the inla.group function below, problem is 
    # the rw struggles with covariate values that are very close to each other (rw1 has same issue)
    # and he recommends these priors to make sure it doesn't get too funky
    U <- 0.5
    hyp.rw2 <- list(theta=list(prior = "pc.prec", param = c(U,0.05)))
    
    # The best cod and yellowtail models
    
    # This first big handles the yellowtail models
    # If I want to use group of the 3 year era so that I can make the temporal field be an AR1 processI need to make the last line
    # of the 3 year models this...         f(w,model=spde,group = w.group,control.group = list(model = 'iid')) # The w.repl is found inside w.index.X
    # For the moment without using the above this bit is all a little silly and unncessary, but I'm leaving it so I can see it later...
    if(species[s] == 'yt_PA' && surveys[i] != 'nmfs-fall')
    {
      model <- y ~ 0 + intercept + f(depth , model = "rw1", hyper = hyp.rw2)  + 
        f(sst , model = "rw1", hyper = hyp.rw2)  + 
        fsed_3 + fsed_4 +
        f(w,model=spde,replicate = w.repl) # The w.repl is found inside w.index.X
    } else { # This gets us the 5 year yellowtail model for the fall.
      model <- y ~ 0 + intercept + f(depth , model = "rw1", hyper = hyp.rw2)  + 
        f(sst , model = "rw1", hyper = hyp.rw2)  + 
        fsed_3 + fsed_4 +
        f(w,model=spde,replicate = w.repl) # The w.repl is found inside w.index.X
    } # end if(species[s] == "yt_PA" && st.mods[st] == 3)
    # 
    # ANd this is what we want for the cod models.#
    if(species[s] == 'cod_PA')
    {
      model <-  y ~ 0 + intercept + f(depth , model = "rw1", hyper = hyp.rw2)  +
        f(sst , model = "rw1", hyper = hyp.rw2)  +
        f(w,model=spde,replicate = w.repl) # The w.repl is found inside w.index.X
    }
    # Let's giver, make the spatial model.
    # Now we can loop this to run over each model and for each stack, 
    # I was going to run a stack loop, but that is amazingly slow, so instead I'm going to do the model selection with the
    # 10 year field, and using the best model from that I'll then run a one off script for the st.5
    # model and see which I prefer....
    
    #if(st == 2) {stk <- skt.5; st.mod <- "st.5"}
    run.name <- paste0(species[s]," ", surveys[i],"_survey")
    print(paste("Model run started at ",Sys.time()))
    r.out <- inla(model, family=fam, data = inla.stack.data(stk),
                  control.predictor=list(A=inla.stack.A(stk),link =1), # The link = 1 is the reason I was getting the transformed predictions, it predicts on the identity scale unless the link is specified (1 means the 1st likelihoods link, which is the only likelihood in this case)
                  control.inla=list(int.strategy='eb'), ## do not integrate over theta, makes the calculation quicker but not to be used for a final model run
                  #verbose=TRUE,
                  control.compute = list(dic=T,waic=T,cpo=T,openmp.strategy="huge"))  # The Openmp.strategy will use more cores, if set to 'huge' it uses everything...
    print(paste("Model run completed at ",Sys.time()))
    # The fitted model, residuals and the covariates, both the standardized and on the original scale.
    
    mo.out <- data.frame(fitted = r.out$summary.fitted.values[e.id,"mean"] , # The expected values can be found with this
                         resid = dat$response[!is.na(dat$response)] - r.out$summary.fitted.values[e.id,"mean"],
                         sd = r.out$summary.fitted.values[e.id,"sd"] ,
                         lci = r.out$summary.fitted.values[e.id,"0.025quant"] ,
                         uci = r.out$summary.fitted.values[e.id,"0.975quant"] ,
                         response = dat$response[!is.na(dat$response)],
                         dep = dat$depth_cen[!is.na(dat$response)],
                         sst = dat$sst_avg_cen[!is.na(dat$response)],
                         depth = dat$comldepth[!is.na(dat$response)],
                         sed = dat$SEDNUM[!is.na(dat$response)],
                         sst_avg = dat$sst_avg[!is.na(dat$response)],
                         years_3 = dat$years_3[!is.na(dat$response)],
                         years_5 = dat$years_5[!is.na(dat$response)],
                         X = dat$X[!is.na(dat$response)],
                         Y = dat$Y[!is.na(dat$response)],
                         year = dat$year[!is.na(dat$response)]
    )
    
    #  a couple of other variables to calculated
    mo.out$var.Y <- 1* mo.out$fitted * (1-mo.out$fitted) # Get the variance, for a Bernoulli it is n*p*(1-p), where n = 1 for a Bernoulli
    mo.out$resid.stan <- mo.out$resid / sqrt(mo.out$var.Y) # Now we can get Pearson residuals
    
    # Now get the predictions
    p.out <- data.frame(pred = r.out$summary.fitted.values[p.id,"mean"] , # The predictions  can be found with this, note these need to be transformed...
                        pred.sd = r.out$summary.fitted.values[p.id,"sd"] ,
                        pred.lci = r.out$summary.fitted.values[p.id,"0.025quant"] ,
                        pred.uci = r.out$summary.fitted.values[p.id,"0.975quant"] ,
                        dep = dat$depth_cen[is.na(dat$response)],
                        #sst = dat$sst_avg_cen[is.na(dat$response)],
                        depth = dat$comldepth[is.na(dat$response)],
                        sst_avg = dat$sst_avg[is.na(dat$response)],
                        sed = dat$SEDNUM[is.na(dat$response)],
                        years_3 = dat$years_3[is.na(dat$response)],
                        years_5 = dat$years_5[is.na(dat$response)],
                        X = dat$X[is.na(dat$response)],
                        Y = dat$Y[is.na(dat$response)],
                        year = dat$year[is.na(dat$response)]
    )
    
    # Now the model fits using dic and waic, results very similar.
    md.out <- data.frame(dic = r.out$dic$dic, 
                         dic.p.eff = r.out$dic$p.eff,
                         waic = r.out$waic$waic, 
                         waic.p.eff = r.out$waic$p.eff,
                         Dispersion = sum()) 
    md.out$Dispersion <- sum(mo.out$resid.stan^2)/ (N-md.out$waic.p.eff)
    
    
    res.pred[[run.name]] <- r.out
    #res[[mod.names[m]]]$model <- mod.names[m]
    mod.output.pred[[run.name]] <- mo.out
    mod.output.pred[[run.name]]$model <- run.name
    mod.output.pred[[run.name]]$species <- species[s]
    mod.output.pred[[run.name]]$survey <- surveys[i]
    #mod.output[[run.name]]$model.id <- mod.names[m]
    #mod.output[[run.name]]$st.era <- st.mods[st]
    pred.output.pred[[run.name]] <- p.out
    pred.output.pred[[run.name]]$model <- run.name
    pred.output.pred[[run.name]]$species <- species[s]
    pred.output.pred[[run.name]]$survey <- surveys[i]
    #mod.output[[run.name]]$model.id <- mod.names[m]
    #mod.output[[run.name]]$st.era <- st.mods[st]
    
    
    mod.diagnostics.pred[[run.name]] <- md.out
    mod.diagnostics.pred[[run.name]]$model <- run.name
    mod.diagnostics.pred[[run.name]]$species <- species[s]
    mod.diagnostics.pred[[run.name]]$survey <- surveys[i]
    #mod.diagnostics[[run.name]]$model.id <- mod.names[m]
    #mod.diagnostics[[run.name]]$st.era <-  st.mods[st]
    
    rand.field.pred[[run.name]] <- r.out$summary.random$w # THis will contain mutliple random fields for the spatio-temporal models.
    rand.field.pred[[run.name]]$model <- run.name
    rand.field.pred[[run.name]]$species <- species[s]
    rand.field.pred[[run.name]]$surveys <- surveys[i]
    #rand.field[[run.name]]$model.id <- mod.names[m]
    #rand.field[[run.name]]$st.era <- st.mods[st]
    
    
    # Save the image on the way through just in case the computer decides to shut down...
#    save.image(paste0(direct.proj,"Results/INLA_model_All_years_full_predicted_field_best_models_",species[s],"_",surveys[i],".RData"))
  }# end for(i in 1:n.surveys)
  # Save the image on the way through just in case the computer decides to shut down, worst case here I lose about half a day of modelling.
  # This thing is gonna be big...
  #save.image(paste0(direct.proj,"Results/INLA_output_full_predicted_field",species[s],"tmp.RData"))
} # end for(s in 1:n.species)

save.image(paste0(direct.proj,"Results/INLA_model_All_years_full_predicted_field_best_models.RData"))
load(paste0(direct.proj,"Results/INLA_model_All_years_full_predicted_field_best_models.RData"))

# Some simple post processing to make smaller files....

pred.dat <- pred.output.pred
# This should go into the Projects folder I think
save(pred.dat,rand.field.pred,file = "Y:/Projects/GB_time_area_closure_SPERA/Results/Paper_3_Closures_and_SDMs/Prediction_and_Rand_fields_all_models_including_2017_2019.RData")
# This sends to the data folder for the Paper 3 project.
save(pred.dat,rand.field.pred,file = paste0(direct.proj,"Data/Prediction_and_Rand_fields_all_models.RData"))
