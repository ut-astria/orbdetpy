from orbdetpy.multi_target import determine_orbit
from orbdetpy.simulation import simulate_measurements
from orbdetpy.plotting.estimation import plot as od_plot
import json


simulate_measurements("data/input/MEO0_sim.json", output_file = "data/output/MEO0_sim_data.json")
simulate_measurements("data/input/MEO1_sim.json", output_file = "data/output/MEO1_sim_data.json")
simulate_measurements("data/input/MEO2_sim.json", output_file = "data/output/MEO2_sim_data.json")

with open("data/output/MEO0_sim_data.json") as json_file:
    MEO0 = json.load(json_file)
with open("data/output/MEO1_sim_data.json") as json_file:
    MEO1 = json.load(json_file)
with open("data/output/MEO2_sim_data.json") as json_file:
    MEO2 = json.load(json_file)

MEO0.extend(MEO1)
MEO0.extend(MEO2)
CombinedData = sorted(MEO0, key = lambda x: (x['Time'], x['Station']))

with open("data/output/Comb_sim_data.json", 'w') as json_file:
    json.dump(CombinedData, json_file, indent=4)


fit_data = determine_orbit(["data/input/MEO0_od.json","data/input/MEO1_od.json"],
                           ["data/output/Comb_sim_data.json"],
                           output_file = "data/output/od_output.json")


with open("data/output/MEO0_OD_Results.json", 'w') as json_file:
    json.dump(fit_data["Estimation"][0], json_file, indent=4)

with open("data/output/MEO1_OD_Results.json", 'w') as json_file:
    json.dump(fit_data["Estimation"][1], json_file, indent=4)

Unassociated = []
with open("data/output/Unassociated.json", 'w') as json_file:
	for measNum in fit_data["Unassociated"]:
		Unassociated.append(CombinedData[measNum])
	json.dump(Unassociated, json_file, indent=4)

MEO0_Associations = []
with open("data/output/MEO0_Associations.json", 'w') as json_file:
	for measNum in fit_data["Associated"][0]:
		MEO0_Associations.append(CombinedData[measNum])
	json.dump(MEO0_Associations, json_file, indent=4)

MEO1_Associations = []
with open("data/output/MEO1_Associations.json", 'w') as json_file:
	for measNum in fit_data["Associated"][1]:
		MEO1_Associations.append(CombinedData[measNum])
	json.dump(MEO1_Associations, json_file, indent=4)


Unassociated = []
with open("data/output/Unassociated.json", 'w') as json_file:
	for measNum in fit_data["Unassociated"]:
		Unassociated.append(CombinedData[measNum])
	json.dump(Unassociated, json_file, indent=4)

od_plot("data/input/MEO0_od.json",
        "data/output/MEO0_Associations.json",
        "data/output/MEO0_OD_Results.json", interactive = True)

od_plot("data/input/MEO1_od.json",
        "data/output/MEO1_Associations.json",
        "data/output/MEO1_OD_Results.json", interactive = True)
