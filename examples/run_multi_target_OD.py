from orbdetpy.multi_target import determine_orbit
from orbdetpy.simulation import simulate_measurements
from orbdetpy.plotting.estimation import plot as od_plot

#simulate_measurements("/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/input/MEO0_sim.json", output_file = "/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/MEO0_sim_data.json")

#simulate_measurements("/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/input/MEO1_sim.json", output_file = "/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/MEO1_sim_data.json")

determine_orbit(["/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/input/MEO1_od.json"], ["/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/MEO1_sim_data.json"], output_file = "/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/MEO_od_output.json")

#determine_orbit(["data/input/MEO0_od.json", "data/input/MEO1_od.json"], 
#                ["/data/output/MEO0_sim_data.json", "/data/output/MEO1_sim_data.json"],
#                output_file = "/data/output/MEO_od_output.json")

od_plot("/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/input/MEO1_od.json", "/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/MEO1_sim_data.json", "/h1/mreinhold/Desktop/Orbdetpy/localRepo/orbdetpy/examples/data/output/MEO_od_output.json", interactive = True)
