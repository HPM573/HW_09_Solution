from enum import Enum
import numpy as np
import SimPy.SamplePathClasses as PathCls
import SimPy.StatisticalClasses as Stat


class HealthStat(Enum):
    """ health status of patients  """
    ALIVE = 1
    DEAD = 0


class Patient:
    def __init__(self, id, mortality_prob):
        """ initiates a patient
        :param id: ID of the patient
        :param mortality_prob: probability of death during a time-step (must be in [0,1])
        """
        self.id = id
        self.rng = np.random.RandomState(seed=id)  # random number generator for this patient
        self.mortalityProb = mortality_prob
        self.healthState = HealthStat.ALIVE  # assuming all patients are alive at the beginning
        self.survivalTime = None   # won't be observed unless the patient dies

    def simulate(self, n_time_steps):
        """ simulate the patient over the specified simulation length """

        t = 0  # simulation current time

        # while the patient is alive and simulation length is not yet reached
        while self.healthState == HealthStat.ALIVE and t < n_time_steps:
            # determine if the patient will die during this time-step
            if self.rng.random_sample() < self.mortalityProb:
                # update the health state to death
                self.healthState = HealthStat.DEAD
                # record the survival time (assuming deaths occurs at the end of this period)
                self.survivalTime = t + 1

            # increment time
            t += 1


class Cohort:
    def __init__(self, id, pop_size, mortality_prob):
        """ create a cohort of patients
        :param id: cohort ID
        :param pop_size: population size of this cohort
        :param mortality_prob: probability of death for each patient in this cohort over a time-step (must be in [0,1])
        """
        self.id = id
        self.initialPopSize = pop_size  # initial population size
        self.patients = []  # list of patients
        self.cohortOutcomes = CohortOutcomes()  # outcomes of the this simulated cohort

        # populate the cohort
        for i in range(pop_size):
            # create a new patient (use id * pop_size + n as patient id)
            patient = Patient(id=id * pop_size + i, mortality_prob=mortality_prob)
            # add the patient to the cohort
            self.patients.append(patient)

    def simulate(self, n_time_steps):
        """ simulate the cohort of patients over the specified number of time-steps
        :param n_time_steps: number of time steps to simulate the cohort
        """
        # simulate all patients
        for patient in self.patients:
            # simulate
            patient.simulate(n_time_steps)

        # store outputs of this simulation
        self.cohortOutcomes.extract_outcomes(self.patients)


class CohortOutcomes:
    def __init__(self):

        self.survivalTimes = []    # survival times
        self.nSurvivedBeyond5yr = 0  # number of patients survived beyond 5 years
        self.meanSurvivalTime = None   # mean survival times
        self.propSurvivedBeyond5yar = 0
        self.nLivingPatients = None   # survival curve (sample path of number of alive patients over time)

    def extract_outcomes(self, simulated_patients):
        """ extracts outcomes of a simulated cohort
        :param simulated_patients: (list) of patients after being simulated """

        # record survival times
        for patient in simulated_patients:
            # if patient's survival time is observed, store it
            if patient.survivalTime is not None:
                self.survivalTimes.append(patient.survivalTime)

            # find if the patient survived beyond 5 years
            # (note that if patient's survival time is not observed, the patient has survived beyond
            # the simulation length)
            if patient.survivalTime is None or patient.survivalTime > 5:
                self.nSurvivedBeyond5yr += 1

        self.meanSurvivalTime = sum(self.survivalTimes)/len(self.survivalTimes)
        self.propSurvivedBeyond5yar = self.nSurvivedBeyond5yr/len(simulated_patients)

        # survival curve
        self.nLivingPatients = PathCls.PrevalencePathBatchUpdate(
           name='# of living patients',
           initial_size=len(simulated_patients),
           times_of_changes=self.survivalTimes,
           increments=[-1]*len(self.survivalTimes))


class MultiCohort:
    """ simulates multiple cohorts with different parameters """

    def __init__(self, ids, pop_sizes, mortality_probs):
        """
        :param ids: (list) of ids for cohorts to simulate
        :param pop_sizes: (list) of population sizes of cohorts to simulate
        :param mortality_probs: (list) of the mortality probabilities
        """
        self.ids = ids
        self.popSizes = pop_sizes
        self.mortalityProbs = mortality_probs
        self.multiCohortOutcomes = MultiCohortOutcomes()

    def simulate(self, n_time_steps):
        """ simulates all cohorts """

        for i in range(len(self.ids)):

            # create a cohort
            cohort = Cohort(id=self.ids[i], pop_size=self.popSizes[i], mortality_prob=self.mortalityProbs[i])

            # simulate the cohort
            cohort.simulate(n_time_steps=n_time_steps)

            # outcomes from simulating all cohorts
            self.multiCohortOutcomes.extract_outcomes(simulated_cohort=cohort)

        # calculate the summary statistics of from all cohorts
        self.multiCohortOutcomes.calculate_summary_stats()


class MultiCohortOutcomes:
    def __init__(self):

        self.survivalTimes = []  # two dimensional list of patient survival times from all simulated cohorts
        self.meanSurvivalTimes = []  # list of average patient survival time for all simulated cohorts
        self.probSurvivedBeyond5ry = []  # list of proportion of patients survived beyond 5yrs for all simulated cohorts
        self.survivalCurves = []  # list of survival curves from all simulated cohorts
        self.statMeanSurvivalTime = None  # summary statistics of mean survival time

    def extract_outcomes(self, simulated_cohort):
        """ extracts outcomes of a simulated cohort
        :param simulated_cohort: a cohort after being simulated"""

        # store all patient survival times from this cohort
        self.survivalTimes.append(simulated_cohort.cohortOutcomes.survivalTimes)

        # store proportion of patients survived beyond 5 yrs
        self.probSurvivedBeyond5ry.append(simulated_cohort.cohortOutcomes.propSurvivedBeyond5yar)

        # append the survival curve of this cohort
        self.survivalCurves.append(simulated_cohort.cohortOutcomes.nLivingPatients)

    def calculate_summary_stats(self):
        """
        calculate the summary statistics
        """

        # calculate average patient survival time for all simulated cohorts
        for obs_set in self.survivalTimes:
            self.meanSurvivalTimes.append(sum(obs_set)/len(obs_set))

        # summary statistics of mean survival time
        self.statMeanSurvivalTime = Stat.SummaryStat(name='Mean survival time',
                                                     data=self.meanSurvivalTimes)

    def get_cohort_CI_mean_survival(self, cohort_index, alpha):
        """
        :returns: the confidence interval of the mean survival time for a specified cohort
        :param cohort_index: integer over [0, 1, ...] corresponding to the 1st, 2nd, ... simulated cohort
        :param alpha: significance level
        """

        stat = Stat.SummaryStat(name='Summary statistics',
                                data=self.survivalTimes[cohort_index])

        return stat.get_t_CI(alpha=alpha)

    def get_cohort_PI_survival(self, cohort_index, alpha):
        """ :returns: the prediction interval of the survival time for a specified cohort
        :param cohort_index: integer over [0, 1, ...] corresponding to the 1st, 2ndm ... simulated cohort
        :param alpha: significance level
        """

        stat = Stat.SummaryStat(name='Summary statistics',
                                data=self.survivalTimes[cohort_index])

        return stat.get_PI(alpha=alpha)
