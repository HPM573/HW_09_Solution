from enum import Enum
import scipy.stats as stat
import numpy as np
import SimPy.InOutFunctions as InOutSupport
import SimPy.StatisticalClasses as Stat
import SurvivalModelClasses as SurvivalCls
import CalibrationSettings as CalibSets
import SimPy.FormatFunctions as FormatSupport


class CalibrationColIndex(Enum):
    """ indices of columns in the calibration results cvs file  """
    ID = 0          # cohort ID
    W = 1  # likelihood weight
    MORT_PROB = 2   # mortality probability

class Calibration:
    def __init__(self):
        np.random.seed(1) #specifying the seed of random # generator
        self.cohortIDs = range(CalibSets.POST_N)
        self.mortalitySamples = []
        self.mortalityResamples = []
        self.weights = []
        self.normalizedWeights = []
        self.csvRows = [['Cohort ID', 'Likelihood Weights' ,'Mortality Prob']]  # list containing the calibration results

    def sample_posterior(self):
        self.mortalitySamples = np.random.uniform(
            low = CalibSets.POST_L,
            high= CalibSets.POST_U,
            size = CalibSets.POST_N
        )

        multiCohort = SurvivalCls.MultiCohort(
            ids = self.cohortIDs,
            mortality_probs = self.mortalitySamples,
            pop_sizes = [CalibSets.SIM_POP_SIZE]*CalibSets.POST_N
        )

        # simulate the multi cohort
        multiCohort.simulate(CalibSets.TIME_STEPS)

        #calculate the likelihood of each simulated cohort
        for cohort_id in self.cohortIDs:

            #get the 5-year survival probability for this cohort
            proportion = multiCohort.get_cohort_prob_five_year_survival(cohort_id)

            weight = stat.binom.pmf(
                k=CalibSets.OBS_ALIVE,
                n=CalibSets.OBS_N,
                p=proportion,
                loc = 0)

            self.weights.append(weight)

        sum_weights = np.sum(self.weights)
        self.normalizedWeights = np.divide(self.weights,sum_weights)

        self.mortalityResamples = np.random.choice(
            a=self.mortalitySamples,
            size = CalibSets.NUM_SIM_COHORTS,
            replace = True,
            p=self.normalizedWeights)

        for i in range(0, len(self.mortalitySamples)):
            self.csvRows.append(
                [self.cohortIDs[i],self.normalizedWeights[i],self.mortalitySamples[i]])


        InOutSupport.write_csv('CalibrationResults.csv',self.csvRows)


    def get_mortality_resamples(self):
        return self.mortalityResamples


    def get_mortality_estimate_credible_interval(self, alpha, deci):

        sum_stat = Stat.SummaryStat('Posterior Samples',self.mortalityResamples)

        estimate = sum_stat.get_mean()
        credible_interval = sum_stat.get_PI(alpha)

        return FormatSupport.format_estimate_interval(estimate, credible_interval, deci)

    def get_effective_sample_size(self):
        return 1 / np.sum(np.power(self.normalizedWeights,2))


class CalibratedModel:
    def __init__(self, csv_file_name, drug_effectiveness_ratio=1):

        cols = InOutSupport.read_csv_cols(
            file_name=csv_file_name,
            n_cols=3,
            if_ignore_first_row=True,
            if_convert_float=True)

        self.cohortIDs = cols[CalibrationColIndex.ID.value].astype(int)
        self.weights = cols[CalibrationColIndex.W.value]
        self.mortalityProbs = cols[CalibrationColIndex.MORT_PROB.value] * drug_effectiveness_ratio
        self.multiCohorts = None

    def simulate(self, num_of_simulated_cohorts, cohort_size,time_steps, cohort_ids=None):

        #resample cohort IDs and mortality probabilities with replacement
        sampled_row_indices = np.random.choice(
            a=range(0,len(self.weights)),
            size=num_of_simulated_cohorts,
            replace=True,
            p=self.weights)

        resampled_ids = []
        resampled_prob = []
        for i in sampled_row_indices:
            resampled_ids.append(self.cohortIDs[i])
            resampled_prob.append(self.mortalityProbs[i])

        if cohort_ids is None:
            self.multiCohorts = SurvivalCls.MultiCohort(
                ids=resampled_ids,
                pop_sizes=[cohort_size]*num_of_simulated_cohorts,
                mortality_probs=resampled_prob)
        else:
            self.multiCohorts = SurvivalCls.MultiCohort(
                ids = cohort_ids,
                pop_sizes=[cohort_size]*num_of_simulated_cohorts,
                mortality_probs=resampled_prob)

        self.multiCohorts.simulate(time_steps)

    def get_mean_survival_time_proj_interval(self, alpha):
        """
        :param alpha: the significance level
        :param deci: decimal places
        :returns tuple in the form of (mean, [lower, upper]) of projection interval
        """
        mean = self.multiCohorts.sumStat_meanSurvivalTimes.get_mean()
        proj_interval = self.multiCohorts.sumStat_meanSurvivalTimes.get_PI(alpha=alpha)

        return mean, proj_interval