from orbdetpy.multi_target import determine_orbit
from orbdetpy.simulation import simulate_measurements
from orbdetpy.plotting.estimation import plot as od_plot
import json


simulate_measurements("data/input/MEO0_sim.json", output_file = "data/output/MEO0_sim_data.json")
simulate_measurements("data/input/MEO1_sim.json", output_file = "data/output/MEO1_sim_data.json")

with open("data/output/MEO0_sim_data.json") as json_file:
    MEO0 = json.load(json_file)
with open("data/output/MEO1_sim_data.json") as json_file:
    MEO1 = json.load(json_file)

for meas0 in MEO0:
    for meas1 in MEO1:
        if (meas0['Time'] == meas1['Time'] and "Station" in meas0 and
            "Station" in meas1 and meas0["Station"] == meas1["Station"]):
            MEO0.insert(MEO0.index(meas0)+1, meas1)
            MEO1.remove(meas1)

with open("data/output/Comb_sim_data.json", 'w') as json_file:
    json.dump(MEO0, json_file, indent=4)

fit_data = determine_orbit(["data/input/MEO0_od.json"],
                           ["data/output/Comb_sim_data.json"],
                           output_file = "data/output/od_output.json")

with open("data/output//MEO0_OD_Results.json", 'w') as json_file:
    json.dump(fit_data["Estimation"][0], json_file, indent=4)
'''
with open("data/output/MEO1_OD_Results.json", 'w') as json_file:
    json.dump(fit_data["Estimation"][1], json_file, indent=4)
'''
Unassociated = []
with open("data/output/Unassociated.json", 'w') as json_file:
	for measNum in fit_data["Unassociated"]:
		Unassociated.append(MEO0[measNum])
	json.dump(Unassociated, json_file, indent=4)

MEO0_Associations = []
with open("data/output/MEO0_Associations.json", 'w') as json_file:
	for measNum in fit_data["Associated"][0]:
		MEO0_Associations.append(MEO0[measNum])
	json.dump(MEO0_Associations, json_file, indent=4)
'''
MEO1_Associations = []
with open("data/output/MEO1_Associations.json", 'w') as json_file:
	for measNum in fit_data["Associated"][1]:
		MEO1_Associations.append(MEO0[measNum])
	json.dump(MEO1_Associations, json_file, indent=4)
'''

with open("data/output/Unassociated.json", 'w') as json_file:
    json.dump(fit_data["Unassociated"], json_file, indent=4)

od_plot("data/input/MEO0_od.json",
        "data/output/MEO0_sim_data.json",
        fit_data["Estimation"][0], interactive = True)
'''
od_plot("data/input/MEO1_od.json",
        "data/output/MEO1_sim_data.json",
        fit_data["Estimation"][1], interactive = True)
'''