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
        self.meanSurvivalTime = None   # mean survival times
        self.nLivingPatients = None   # survival curve (sample path of number of alive patients over time)
        self.countSurvival_list = []
        self.percentageSurviveFive = None

    def extract_outcomes(self, simulated_patients):
        """ extracts outcomes of a simulated cohort
        :param simulated_patients: (list) of patients after being simulated """

        # record survival times
        for patient in simulated_patients:
            if not (patient.survivalTime is None):
                self.survivalTimes.append(patient.survivalTime)

        for k in self.survivalTimes:
            if k > 5:
                i=1
                self.countSurvival_list.append(i)
            else:
                i=0
                self.countSurvival_list.append(i)

        self.meanSurvivalTime = sum(self.survivalTimes)/len(self.survivalTimes)
        self.percentageSurviveFive = sum(self.countSurvival_list)/len(self.countSurvival_list)

        # survival curve
        self.nLivingPatients = PathCls.PrevalencePathBatchUpdate(
           name='# of living patients',
           initial_size=len(simulated_patients),
           times_of_changes=self.survivalTimes,
           increments=[-1]*len(self.survivalTimes))

    def get_survival_time(self):
        return self.survivalTimes

    def get_ave_survival_time(self):
        return self.meanSurvivalTime

    def get_five_year_survival_list(self):
        return self.countSurvival_list

    def get_percentage_survive_five(self):
        return self.percentageSurviveFive



class MultiCohort:
    def __init__(self, ids, pop_sizes, mortality_probs):
        self.ids = ids
        self.popSizes = pop_sizes
        self.mortalityProbs = mortality_probs

        self.survivalTimes = []
        self.meanSurvivalTimes = []
        self.sumStat_meanSurvivalTimes = None

        self.fiveYearSurvivals = []
        self.probFiveYearSurvival = []
        self.sumStat_probFiveYearSurvival = None

    def simulate(self, n_steps):

        for i in range(len(self.ids)):
            cohort = Cohort(self.ids[i], self.popSizes[i],self.mortalityProbs[i])
            cohort.simulate(n_steps)

            self.survivalTimes.append(cohort.cohortOutcomes.get_survival_time())
            self.meanSurvivalTimes.append(cohort.cohortOutcomes.get_ave_survival_time())
            self.fiveYearSurvivals.append(cohort.cohortOutcomes.get_five_year_survival_list())
            self.probFiveYearSurvival.append(cohort.cohortOutcomes.get_percentage_survive_five())

        self.sumStat_meanSurvivalTimes = Stat.SummaryStat('Mean Surviavl Time',self.meanSurvivalTimes)
        self.sumStat_probFiveYearSurvival = Stat.SummaryStat('Five Year Survival', self.probFiveYearSurvival)

    def get_cohort_mean_survival(self,cohort_index):
        return self.meanSurvivalTimes[cohort_index]

    def get_cohort_CI_mean_surviavl(self, cohort_index, alpha):
        st = Stat.SummaryStat('', self.survivalTimes[cohort_index])
        return st.get_t_CI(alpha)

    def get_cohort_prob_five_year_survival(self, cohort_index):
        return self.probFiveYearSurvival[cohort_index]

    def get_cohort_CI_prob_five_year_survival(self,cohort_index,alpha):
        st2 = Stat.SummaryStat('',self.fiveYearSurvivals[cohort_index])
        return st2.get_t_CI(alpha)

    def get_all_mean_survival(self):
        return self.meanSurvivalTimes

    def get_overall_mean_survival(self):
        return self.sumStat_meanSurvivalTimes.get_mean()

    def get_all_five_year_prob(self):
        """ :returns a list of five year survival probability for all simulated cohorts"""
        return self.probFiveYearSurvival

    def get_overall_five_year_prob(self):
        """ :returns the overall five year survival probability (the mean of the five year survival probabilities"""
        return self.sumStat_probFiveYearSurvival.get_mean()

    # Projection stats
    def get_cohort_PI_survival(self, cohort_index, alpha):
        """ :returns: the prediction interval of the survival time for a specified cohort
        :param cohort_index: integer over [0, 1, ...] corresponding to the 1st, 2ndm ... simulated cohort
        :param alpha: significance level
        """
        st = Stat.SummaryStat('', self.survivalTimes[cohort_index])
        return st.get_PI(alpha)

    def get_PI_mean_survival(self, alpha):
        """ :returns: the prediction interval of the mean survival time"""
        return self.sumStat_meanSurvivalTimes.get_PI(alpha)

    def get_cohort_PI_five_year_prob(self, cohort_index, alpha):
        st3 = Stat.SummaryStat('', self.probFiveYearSurvival[cohort_index])
        return st3.get_PI(alpha)

    def get_PI_five_year_prob(self,alpha):
        return self.sumStat_probFiveYearSurvival.get_PI(alpha)

