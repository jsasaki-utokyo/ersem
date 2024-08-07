import argparse
import json
import netCDF4 as nc
from numpy import ndarray, interp

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, ndarray):
            temp = [float(v) for v in obj.tolist()]
            return temp
        return json.JSONEncoder.default(self, obj)

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--data-path', type=str, required=True,
                    help='Path to output file from model run')
parser.add_argument('-m', '--model-run', type=str, required=True,
                    help='Type of model run that has been run',
                    choices=["gotm", "fabm0d"])
args, _ = parser.parse_known_args()
data_path = args.data_path
model_run = args.model_run

dir_mapping = {"fabm0d": "fabm0d-gotm-ersem", "gotm": "gotm-fabm-ersem"}
if model_run == "gotm":
    state_vars = [
        "N1_p", "N3_n", "N4_n", "N5_s", "O2_o", "O3_c", "O3_bioalk", "R1_c", "R1_n",
        "R1_p", "R2_c", "R3_c", "R4_c", "R4_n", "R4_p", "R6_c", "R6_n", "R6_p",
        "R6_s", "R8_c", "R8_n", "R8_p", "R8_s", "B1_c", "B1_n", "B1_p", "P1_c",
        "P1_n", "P1_p", "P1_Chl", "P1_s", "P2_c", "P2_n", "P2_p", "P2_Chl", "P3_c",
        "P3_n", "P3_p", "P3_Chl", "P4_c", "P4_n", "P4_p", "P4_Chl", "Z4_c", "Z5_c",
        "Z5_n", "Z5_p", "Z6_c", "Z6_n", "Z6_p", "L2_c", "Q1_c", "Q1_p", "Q1_n",
        "Q6_c", "Q6_p", "Q6_n", "Q6_s", "Q6_pen_depth_c", "Q6_pen_depth_n",
        "Q6_pen_depth_p", "Q6_pen_depth_s", "Q7_c", "Q7_p", "Q7_n", "Q7_pen_depth_c",
        "Q7_pen_depth_n", "Q7_pen_depth_p", "Q17_c", "Q17_p", "Q17_n", "bL2_c",
        "ben_col_D1m", "ben_col_D2m", "K1_p", "K3_n", "K4_n", "K5_s", "G2_o",
        "G2_o_deep", "G3_c", "ben_nit_G4n", "H1_c", "H2_c", "Y2_c", "Y3_c",
        "Y4_c"
    ]
    vars_of_interest = ["dates", "N1_p", "N3_n", "N5_s"]
    data_dict = {"expected": vars_of_interest, "expected_state": state_vars}

elif model_run == "fabm0d":
    vars_of_interest = [
        "dates", "light_parEIR", "temp", "salt", "N1_p", "N3_n", "B1_c", "P2_c"
    ]

    data_dict = {"expected": vars_of_interest}

for key, items in data_dict.items():
    data = nc.Dataset(data_path, 'r')
    expected_results = {}
    for v in items:
        try:
            if v == "dates":
                times = data.variables['time']
                dates = nc.num2date(times[:],
                                    units=times.units,
                                    calendar=times.calendar)
                dates = [str(d).split(" ")[0] for d in dates]
                expected_results[v] = dates
            elif key == "expected" and model_run == "gotm":
                depth = 0.0
                var = data.variables[v]
                zi = data.variables['zi'][:].squeeze()
                z = data.variables['z'][:].squeeze()
                var_time_series = []
                for i in range(var.shape[0]):
                    depth_offset = depth + zi[i, -1]
                    var_time_series.append(interp(depth_offset, z[i, :], var[i, :].squeeze()))
                expected_results[v] = var_time_series

            elif data.variables[v].ndim == 4:
                expected_results[v] = data.variables[v][:].squeeze() \
                        if model_run == "fabm0d" else expected_results[v][-1, :]
            elif data.variables[v].ndim == 3:
                expected_results[v] = data.variables[v][:].squeeze() \
                        if model_run == "fabm0d" else \
                        float(data.variables[v][:].squeeze()[-1])
            else:
                raise RuntimeError
        except Exception as e:
            print(f"Did not update {v} since {e}")

    with open(f'{dir_mapping[model_run]}/{key}.json', 'w') as f:
        json.dump(expected_results, f, cls=NumpyEncoder)
