require(brms)
require(tidyverse) 
require(parallel)
require(reshape2)

T1 = TRUE
filename = 'aseg'
db = read_csv(paste0('raw/db_',filename,'.csv'))
if (T1){
  db = subset(db, time == 'T1')
  mod = '1 + Clinical_Trajectory + Age + Gender + EstimatedTotalIntraCranialVol'
  ext = '_t1'
} else {
  mod = '1 + time * Clinical_Trajectory + Age + Gender + EstimatedTotalIntraCranialVol'
  ext = ''
}

iterations <- 4e4
chains <- 4
SCALE <- 1
ns <- iterations*chains/2

db['sub']  <- factor(db$sub)
db['Clinical_Trajectory'] <- factor(db$Clinical_Trajectory)
db['ROI'] <- factor(db$ROI)

# number of ROIs
print(paste0('Number of ROI: ',nlevels(db$ROI)))
print(paste0('Number of subs: ',nlevels(db$sub)))
print(paste0('Number of cores available: ', detectCores(all.tests = FALSE, logical = TRUE)))

print(getOption("mc.cores"))
options(mc.cores = parallel::detectCores())
print(getOption("mc.cores"))

modelForm = paste('Volume ~',mod,'+ (1 | gr(sub, dist= "student")) + (',mod,'| gr(ROI, dist="student"))')

priorRBA <- get_prior(formula = modelForm,data=db,family = 'student')

if (T1){
  priorRBA$prior[1] <- "student_t(3,0,10)"
  priorRBA$prior[6] <- "lkj(2)"
  priorRBA$prior[8:9] <- "gamma(3.325,0.1)"
  priorRBA$prior[11] <- "gamma(3.325,0.1)"
  print(priorRBA)
} else {
  priorRBA$prior[1] <- "student_t(3,0,10)"
  priorRBA$prior[10] <- "lkj(2)"
  priorRBA$prior[12:13] <- "gamma(3.325,0.1)"
  priorRBA$prior[15] <- "gamma(3.325,0.1)"
  print(priorRBA)
}

outdir <- paste0('data/')
if (!dir.exists(outdir)){dir.create(outdir,recursive = TRUE)}

# Generate the Stan code for our own reference
stan_code <- make_stancode(modelForm,
                           data=db,
                           chains = chains,
                           family = 'student',
                           prior = priorRBA)
cat(stan_code,file = paste0(outdir,"/stancode4.stan"),sep='\n')

# Following run the BML model
fm <- brm(modelForm,
          data=db,
          chains = chains,
          family = 'student',
          prior = priorRBA,
          inits=0, iter=iterations, 
          cores = 4,
          control = list(adapt_delta = 0.99, max_treedepth = 15))

save.image(file=paste0(outdir,filename,ext,'.RData'))
