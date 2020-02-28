from orbdetpy.multi_target import determine_orbit
from orbdetpy.simulation import simulate_measurements
from orbdetpy.plotting.estimation import plot as od_plot

simulate_measurements("data/input/radar_sim_cfg.json",
                      output_file = "/tmp/sim_data.json")

determine_orbit(["data/input/radar_od_cfg.json"], ["/tmp/sim_data.json"],
                output_file = "/tmp/od_output.json")

od_plot("data/input/radar_od_cfg.json", "/tmp/sim_data.json", "/tmp/od_output.json",
        interactive = True)
