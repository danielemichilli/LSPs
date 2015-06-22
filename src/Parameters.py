#############################
#
# Parameters
#
# Written by Daniele Michilli
#
#############################


DM_MIN = 5.  #pc/cm3
SNR_MIN = 6.5
#DURATION_MAX = 0.049152  #s


F_MIN = 119.43  #MHz
F_MAX = F_MIN + 31.64  #MHz

#Grouping
ERR_FLOAT = 0.001  #Error on float numbers
STEPS_GROUP = 20  #Tollerance on the number of DM steps without a candidate
DURAT_GROUP = 0.03  #Tollerance in seconds on the temporal shift of two consecutive events

SIGMA_TOLL = 4
SIGMA_TOLL_IB = 2

RFI_percent = 3
PULS_LENGHT = 15


FILTERS = {
  'scattered': 0.09,
  'aligned': 0.045,
  'peak_central': 0.6,
  'duration_central': 1,
  'holes': 0.75,
  'variance': 0.06,
  'flat_duration': 0.9,
  'm': 23.684210526315788,
  'q': -60.86842105263156,
  'flat_SNR': 0.78,
  'flat_SNR_simmetric': 0.78,
  'DM_extremes': 0.9,
  'sigma_min': 0.95,
  'cum_scatter': 0.17,
  'sigma_std': 0.46,
  'sigma_scatter': 0.166,
  'sigma_scatter_max': 0.47,
  'sigma_std_largest': 0.55,
  'flat_fit0': 0.72,
  'flat_fit1': 0.68}


