from orbdetpy.multi_target import determine_orbit
from orbdetpy.simulation import simulate_measurements
from orbdetpy.plotting.estimation import plot as od_plot
import json


fit_data = determine_orbit(["data/input/ASTRIA_od.json"],
                           ["data/input/CombinedASTRIAData.json"],
                           output_file = "data/output/ASTRIA_output.json")


with open("data/output/00040A_OD_Results.json", 'w') as json_file:
    json.dump(fit_data["Estimation"][0], json_file, indent=4)

with open("data/output/18109A_OD_Results.json", 'w') as json_file:
    json.dump(fit_data["Estimation"][1], json_file, indent=4)

Unassociated = []
with open("data/output/Unassociated.json", 'w') as json_file:
	for measNum in fit_data["Unassociated"]:
		Unassociated.append(CombinedData[measNum])
	json.dump(Unassociated, json_file, indent=4)

MEO0_Associations = []
with open("data/output/00040A_Associations.json", 'w') as json_file:
	for measNum in fit_data["Associated"][0]:
		MEO0_Associations.append(CombinedData[measNum])
	json.dump(MEO0_Associations, json_file, indent=4)

MEO1_Associations = []
with open("data/output/18109A_Associations.json", 'w') as json_file:
	for measNum in fit_data["Associated"][1]:
		MEO1_Associations.append(CombinedData[measNum])
	json.dump(MEO1_Associations, json_file, indent=4)


od_plot("data/input/ASTRIA_od.json",
        "data/output/00040A_Associations.json",
        "data/output/00040A_OD_Results.json", interactive = True)

od_plot("data/input/MEO1_od.json",
        "data/output/18109A_Associations.json",
        "data/output/18109A_OD_Results.json", interactive = True)

