#############################
#
# Parameters
#
# Written by Daniele Michilli
#
#############################


DM_MIN = 3.  #pc/cm3
SNR_MIN = 6.5
#DURATION_MAX = 0.049152  #s


F_MIN = 119.43  #MHz
F_MAX = F_MIN + 31.64  #MHz
RES = 491.52e-6  #s
DURATION = 3599.719  #s


#Grouping
ERR_FLOAT = 0.001  #Error on float numbers
STEPS_GROUP = 20  #Tollerance on the number of DM steps without a candidate
DURAT_GROUP = 0.03  #Tollerance in seconds on the temporal shift of two consecutive events

SIGMA_TOLL = 4
SIGMA_TOLL_IB = 2

RFI_percent = 2
PULS_LENGHT = 15

DS_OFFSET = 100000 #bin


#DM step changes
DM_STEP1 = 41.04   #40.49
DM_STEP2 = 133.14  #141.69


