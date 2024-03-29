import scipy.stats as stat

import CalibrationClasses as Cls
import CalibrationSettings as Sets

# Problem 2:
print('Problem 2: A binomial distribution with q as the probability of success and N as the number of trials.')

# Problem 3:
weight = stat.binom.pmf(k=400, n=573, p=0.55)
print('Problem 3: The likelihood of the observed data:', weight)

# Problem 4:
# create a calibration object
calibration = Cls.Calibration()

# sample the posterior of the mortality probability
calibration.sample_posterior(n_samples=Sets.PRIOR_N)

# Problem 5:
# initialize and simulate a calibrated model
calibrated_model = Cls.CalibratedModel(csv_file_name='CalibrationResults.csv')
calibrated_model.simulate(num_of_simulated_cohorts=Sets.NUM_SIM_COHORTS,
                          cohort_size=Sets.SIM_POP_SIZE,
                          time_steps=Sets.TIME_STEPS)

# report mean and PI
print('Problem 5a: Mean survival time and {:.{prec}%} projection interval:'.format(1 - Sets.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(Sets.ALPHA))
# estimate of annal mortality probability and the 95% credible interval
print('Problem 5b: Estimate of annual mortality probability ({:.{prec}%} credible interval):'
      .format(1 - Sets.ALPHA, prec=0),
      calibrated_model.get_mortality_estimate_credible_interval(Sets.ALPHA))

# Problem 6:
Sets.OBS_N = 1146
Sets.OBS_ALIVE = 800

# create a calibration object
calibration = Cls.Calibration()
# sample the posterior of the mortality probability
calibration.sample_posterior(n_samples=Sets.PRIOR_N)

# initialize and simulate a calibrated model
calibrated_model = Cls.CalibratedModel(csv_file_name='CalibrationResults.csv')
calibrated_model.simulate(num_of_simulated_cohorts=Sets.NUM_SIM_COHORTS,
                          cohort_size=Sets.SIM_POP_SIZE,
                          time_steps=Sets.TIME_STEPS)

# report mean and PI
print('Problem 6a: Mean survival time and {:.{prec}%} projection interval:'.format(1 - Sets.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(Sets.ALPHA))
# estimate of annal mortality probability and the 95% credible interval
print('Problem 6b: Estimate of annual mortality probability ({:.{prec}%} credible interval):'
      .format(1 - Sets.ALPHA, prec=0),
      calibrated_model.get_mortality_estimate_credible_interval(Sets.ALPHA))