from orbdetpy.multi_target import determine_orbit
from orbdetpy.simulation import simulate_measurements
from orbdetpy.plotting.estimation import plot as od_plot
import json

'''
simulate_measurements("data/input/MEO0_sim.json", output_file = "data/output/MEO0_sim_data.json")

simulate_measurements("data/input/MEO1_sim.json", output_file = "data/output/MEO1_sim_data.json")

with open("data/output/MEO0_sim_data.json") as json_file:
    MEO0 = json.load(json_file)
with open("data/output/MEO1_sim_data.json") as json_file:
    MEO1 = json.load(json_file)

for meas0 in MEO0:
    for meas1 in MEO1:
        if meas0['Time'] == meas1['Time'] and "Station" in meas0 and "Station" in meas1 and meas0["Station"] == meas1["Station"]:
            MEO0.insert(MEO0.index(meas0)+1, meas1)
            MEO1.remove(meas1)
    
with open("data/output/Comb_sim_data.json", 'w') as json_file:
  json.dump(MEO0, json_file)
'''

determine_orbit(["data/input/MEO0_od.json",
    "data/input/MEO1_od.json"],
     ["data/output/Comb_sim_data.json"],
     output_file = "data/output/MEO_od_output.json")

od_plot("examples/data/input/MEO0_od.json",
 "data/output/MEO0_sim_data.json",
 "data/output/MEO_od_output.json", interactive = True)

od_plot("data/input/MEO1_od.json",
 "data/output/MEO1_sim_data.json",
 "data/output/MEO_od_output.json", interactive = True)
