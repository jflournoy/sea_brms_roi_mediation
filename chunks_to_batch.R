## @knitr load_data

library(dplyr)
library(tidyr)
library(brms)

model_obj_dir <- 'model_obj'
model_text_dir <- 'model_text'
if(!dir.exists(model_obj_dir)){
  message('Creating directory: ', file.path(model_obj_dir))
  dir.create(model_obj_dir)
}
if(!dir.exists(model_text_dir)){
  message('Creating directory: ', file.path(model_text_dir))
  dir.create(model_text_dir)
}

if(Sys.getenv('HOME') != '/users/jflournoy'){
  task_id <- 1
  cpus_per_task <- 4
  message('Not running on SLURM system')
  warmup <- 200
  iter <- 500
} else {
  task_id <- Sys.getenv('SLURM_ARRAY_TASK_ID')
  cpus_per_task <- Sys.getenv('SLURM_CPUS_PER_TASK')
  warmup <- 2500
  iter <- 7500
  message('Running on SLURM system')
}

seed <- 107 + task_id #random.org

rois_fearGTcalm <- readr::read_csv('data/Fear_GT_Calm04-29-2020.csv')
rois_fearGTcalm_l <- rois_fearGTcalm %>%
  extract(SubjectID, into = c('id', 'wave'), regex = '(\\d{4})_(\\d{2})') %>%
  gather(key = 'key', value = 'value', -id, -wave) %>%
  extract(key, into = c('model', 'number', 'roi'), regex = '([a-z._]+)(?:\\.|_)(\\w+)_X_(.*)', remove = FALSE) %>%
  mutate(roi = case_when(roi == 'FreeSurfer_a2009s' ~ 'All',
                         TRUE ~ roi),
         wave = as.numeric(wave)) %>%
  arrange(model, number, roi)

rois_fearGTcalm_l_disag <- rois_fearGTcalm_l %>%
  group_by(roi, model, number) %>%
  mutate(grand_mean = mean(value, na.rm = T)) %>%
  group_by(id, roi, model, number) %>%
  mutate(brain = value,
         wi_brain_mean = mean(brain, na.rm = T),
         wi_brain = brain - wi_brain_mean,
         bw_brain = wi_brain_mean - grand_mean) %>%
  select(-value)

psych <- readr::read_csv('data/stress_psych-midpoint5.csv') %>%
  mutate(idnum = as.character(idnum)) %>%
  select(idnum, time, TIMECENTER, 
         PHQ9_TOT, PHQ9_TOT_lead,
         GAD7_TOT, GAD7_TOT_lead,
         LSA_ANX, LSA_ANX_lead,
         matches('[WG]CEN_CHRONICSEV'),
         matches('[WG]CEN_EPISODICTOT')) %>%
  gather(key, X, 
         matches('[WG]CEN_CHRONICSEV'), 
         matches('[WG]CEN_EPISODICTOT')) %>%
  extract(key, into = c('disagg', 'X_name'), regex = '([WG]CEN)_(.*)') %>%
  mutate(disagg = case_when(disagg == 'WCEN' ~ 'wi_X',
                            disagg == 'GCEN' ~ 'bw_X',
                            TRUE ~ 'ERROR')) %>%
  spread(disagg, X) %>%
  gather(Y_name, Y, 
         PHQ9_TOT, PHQ9_TOT_lead, 
         GAD7_TOT, GAD7_TOT_lead,
         LSA_ANX, LSA_ANX_lead)

dat <- left_join(filter(rois_fearGTcalm_l_disag, roi == 'All'), psych, 
                 by = c('id' = 'idnum', 'wave' = 'time'))

variable_combos <- dplyr::bind_rows(list(
  expand.grid(model = 'chsev.fear', 
              number = unique(dat$number[dat$model == 'chsev.fear']),
              X_name = 'CHRONICSEV',
              Y_name = unique(dat$Y_name),
              stringsAsFactors = FALSE),
  expand.grid(model = 'eptot.fear', 
              number = unique(dat$number[dat$model == 'eptot.fear']),
              X_name = 'EPISODICTOT',
              Y_name = unique(dat$Y_name),
              stringsAsFactors = FALSE),
  expand.grid(model = 'gad_lead.fear', 
              number = unique(dat$number[dat$model == 'gad_lead.fear']),
              X_name = unique(dat$X_name),
              Y_name = 'GAD7_TOT_lead',
              stringsAsFactors = FALSE),
  expand.grid(model = 'socanx_lead.fear', 
              number = unique(dat$number[dat$model == 'socanx_lead.fear']),
              X_name = unique(dat$X_name),
              Y_name = 'LSA_ANX_lead',
              stringsAsFactors = FALSE)
))

if(task_id > dim(variable_combos)[1] | task_id < 1){
  stop(paste0('Task ID is out of bounds: ', task_id))
}

## @knitr estimate_model

varselect <- variable_combos[task_id,]
fname <- paste(varselect, collapse = '_')

medmod_dat <- filter(dat, 
                     model == varselect[['model']], 
                     number == varselect[['number']],
                     X_name == varselect[['X_name']],
                     Y_name == varselect[['Y_name']])


priors_medmods <- c(set_prior("normal(0, 20)", class = "b", resp = "brain"),
                    set_prior("normal(0, 20)", class = "b", resp = "Y"),
                    set_prior("normal(0, 20)", class = "Intercept", resp = "brain"),
                    set_prior("normal(0, 20)", class = "Intercept", resp = "Y"),
                    set_prior("weibull(2, 1)", class = "sigma", resp = "brain"),
                    set_prior("weibull(2, 1)", class = "sigma", resp = "Y"),
                    set_prior("lkj(1)", class = "cor"),
                    set_prior("lkj(1)", class = "rescor"),
                    set_prior("cauchy(0,2.5)", class = "sd", resp = "brain"),
                    set_prior("cauchy(0,2.5)", class = "sd", resp = "Y"))

ymod <- brms::bf(Y ~ 1 + TIMECENTER + 
                   bw_X + wi_X + 
                   bw_brain + wi_brain +
                   (1 + wi_X + wi_brain | ID | id))
mmod <- brms::bf(brain ~ 1 + TIMECENTER + 
                   bw_X + wi_X +
                   (1 + wi_X | ID | id))

message('Fitting model...')
medmod_fit_throwaway <- brms::brm(ymod + mmod, 
                                  data = medmod_dat,
                                  family = gaussian(),
                                  prior = priors_medmods, 
                                  chains = cpus_per_task, cores = cpus_per_task,
                                  seed = seed,
                                  warmup = warmup, iter = iter,
                                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                                  save_model = file.path(model_text_dir, fname),
                                  save_dso = TRUE,
                                  file = file.path(model_obj_dir, fname),
                                  silent = TRUE, open_progress = FALSE)
