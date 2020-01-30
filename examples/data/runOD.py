from orbdetpy.simulation import simulate_measurements
from orbdetpy.estimation import determine_orbit
from orbdetpy.plotting.estimation import plot as od_plot

simulate_measurements("input/radar_sim_cfg.json", output_file = "output/sim_data.json")

determine_orbit("input/radar_od_cfg.json", "output/sim_data.json", output_file = "output/od_output.json")

od_plot("input/radar_od_cfg.json", "output/sim_data.json", "output/od_output.json", interactive = True)
