# Aim: To be able to create a feature to ON & apply corrections & give corrected measurements
# For now applicable only to:  real-data TDM files which have SP3c files
# Author: ATK on 08/02/2022 - WORKING!
import os, json
from orbdetpy.utilities import importSP3wcorr
from orbdetpy.conversion import get_J2000_epoch_offset, get_UTC_string

# pythonic input: input TDM file (org. meas), input SP3c file
# pythonic output: corrected measurement
# later: amount of bias in (RA, Dec), LT corr. ON/OFF, Abbr. ON/OFF

curr_dir = os.getcwd()
test_data_dir = curr_dir
test_output_data_dir = curr_dir

final_testing = 1
if final_testing:
    input_tdm_filename = os.path.join(test_data_dir, "20200917_041957_001_03005A_UTA-ASTRIANet-02.tdm")
    input_sp3c_filename = os.path.join(test_data_dir, "COD21234.EPH_M")


turn_on_corr = 1
sp3corr_out = importSP3wcorr(tdmfilename=input_tdm_filename, sp3filename=input_sp3c_filename,
                             outfilepath=test_output_data_dir, corr_on=turn_on_corr)

output = list(sp3corr_out)
for line in output[:3]:
    print(get_UTC_string(line.time), line.values, line.true_state)


print("len of output: ", len(sp3corr_out))
print(get_UTC_string(sp3corr_out[0].time))
# Get GPS state from sp3c file
print(sp3corr_out[0].true_state)
# Get corrected RA/Dec for TDM measurements
print(sp3corr_out[0].values)
print(get_UTC_string(sp3corr_out[-1].time))
print(sp3corr_out[-1].true_state)

