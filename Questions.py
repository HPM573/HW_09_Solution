import scipy.stats as stat
import CalibrationClasses as Cls
import CalibrationSettings as CalibSets


# Problem 2:
print('Problem 2: A binomial distribution with q as the probability of success and N as the number of trials.')

# Problem 3:
weight = stat.binom.pmf(k=400, n=573, p=0.55)
print('Problem 3: The likelihood of the observed data:', weight)

# Problem 4:
# create a calibration object
calibration = Cls.Calibration()

# sample the posterior of the mortality probability
calibration.sample_posterior(n_samples=CalibSets.POST_N)

# estimate of annal mortality probability and the 95% credible interval
print('Problem 4: Estimate of annual mortality probability ({:.{prec}%} credible interval):'
      .format(1-CalibSets.ALPHA, prec=0),
      calibration.get_mortality_estimate_credible_interval(CalibSets.ALPHA))

# Problem 5:
# initialize and simulate a calibrated model
calibrated_model = Cls.CalibratedModel(csv_file_name='CalibrationResults.csv')
calibrated_model.simulate(num_of_simulated_cohorts=CalibSets.NUM_SIM_COHORTS,
                          cohort_size=CalibSets.SIM_POP_SIZE,
                          time_steps=CalibSets.TIME_STEPS)

# report mean and PI
print('Problem 5: Mean survival time and {:.{prec}%} projection interval:'.format(1-CalibSets.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(CalibSets.ALPHA))

# Problem 6:
CalibSets.OBS_N = 1146
CalibSets.OBS_ALIVE = 800

# create a calibration object
calibration = Cls.Calibration()
# sample the posterior of the mortality probability
calibration.sample_posterior(n_samples=CalibSets.POST_N)
# estimate of annal mortality probability and the 95% credible interval
print('Problem 6: Estimate of annual mortality probability ({:.{prec}%} credible interval):'
      .format(1-CalibSets.ALPHA, prec=0),
      calibration.get_mortality_estimate_credible_interval(CalibSets.ALPHA))

# initialize and simulate a calibrated model
calibrated_model = Cls.CalibratedModel(csv_file_name='CalibrationResults.csv')
calibrated_model.simulate(num_of_simulated_cohorts=CalibSets.NUM_SIM_COHORTS,
                          cohort_size=CalibSets.SIM_POP_SIZE,
                          time_steps=CalibSets.TIME_STEPS)

# report mean and PI
print('Problem 6: Mean survival time and {:.{prec}%} projection interval:'.format(1-CalibSets.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(CalibSets.ALPHA))
