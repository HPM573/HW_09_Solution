import scipy.stats as stat
import CalibrationClasses as Cls
import CalibrationSettings as CalibSets
import CalibrationClassesQ6 as ClsQ6
### Problem 2

print('Problem 2: The best distribution to represent the scenario is the binomial distribution, with a parameter q as the probability of five year survival.')



### Problem 3

weight = stat.binom.pmf(k=400, n=573, p=0.55, loc=0)

print('Problem 3: The likelihood of the observed data:', weight)


### Problem 4
# create a calibration object
calibration = Cls.Calibration()

# sample the posterior of the mortality probability
calibration.sample_posterior()

# estimate of annal mortality probability and the 95% credible interval
print('Problem 4: Estimate of annual mortality probability ({:.{prec}%} credible interval):'.format(1-CalibSets.ALPHA, prec=0),
      calibration.get_mortality_estimate_credible_interval(CalibSets.ALPHA, 4))  # return with 4 decimals

### Problem 5

# initialize and simulate a calibrated model
calibrated_model = Cls.CalibratedModel('CalibrationResults.csv')
calibrated_model.simulate(CalibSets.NUM_SIM_COHORTS, CalibSets.SIM_POP_SIZE, CalibSets.TIME_STEPS)

# report mean and PI
print('Problem 5: Mean survival time and {:.{prec}%} projection interval:',calibrated_model.get_mean_survival_time_proj_interval(CalibSets.ALPHA))

### Problem 6
calibration = ClsQ6.Calibration()
calibration.sample_posterior()

# report mean and PI from clinical trial
print('Problem 6: Estimate of annual mortality probability ({:.{prec}%} credible interval):'.format(1-CalibSets.ALPHA, prec=0),
      calibration.get_mortality_estimate_credible_interval(CalibSets.ALPHA, 4))  # return with 4 decimals

# initialize and simulate a calibrated model
calibrated_model = ClsQ6.CalibratedModel('CalibrationResults.csv')
calibrated_model.simulate(CalibSets.NUM_SIM_COHORTS, CalibSets.SIM_POP_SIZE, CalibSets.TIME_STEPS)

print('     Mean survival time and {:.{prec}%} projection interval:'.format(1 - CalibSets.ALPHA, prec=0),
      calibrated_model.get_mean_survival_time_proj_interval(CalibSets.ALPHA))