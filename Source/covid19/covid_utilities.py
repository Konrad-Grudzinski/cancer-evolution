import numpy as np
import json
import os
import pandas as pd

# ================ common variables ===================


# base paths, for Windows and WSL

# where the downloaded data is located
DATA_PATH = r"D:\ncbi_dataset\ncbi_dataset\data"
DATA_PATH_WSL = "/mnt/d/ncbi_dataset/ncbi_dataset/data/"
# where meta data and variant information is stored
COVID_PATH = r"C:\Users\Konrad Grudzinski\OneDrive - University of Glasgow\Computing\4th Year\Individual Project\Source\covid19"
# where tools are located (if they are not in the PATH variable)
SOURCE_PATH_WSL = "/mnt/c/Users/Konrad Grudzinski/OneDrive - University of Glasgow/Computing/4th Year/Individual Project/Source/"

# variants available on Covariants.org on the 2nd March 2022
variant_names = [
    '20I (Alpha, V1)',
    '20H (Beta, V2)',
    '20J (Gamma, V3)',
    '21A (Delta)',
    '21I (Delta)',
    '21J (Delta)',
    '21K (Omicron)',
    '21L (Omicron)',
    '21B (Kappa)',
    '21D (Eta)',
    '21F (Iota)',
    '21G (Lambda)',
    '21H (Mu)',
    '20E (EU1)',
    '21C (Epsilon)',
    '20B/S:732A',
    '20A/S:126A',
    '20A/S:439K',
    'S:677H.Robin1',
    'S:677P.Pelican',
    '20A.EU2',
    '20A/S:98F',
    '20C/S:80Y',
    '20B/S:1122L',
    '20B/S:626S'
]
variant_colors = [
    # alpha
    "royalblue",
    # beta
    "violet",
    # gamma
    "teal",
    # delta
    "lightsalmon",
    "tomato",
    "salmon",
    # omicron
    "lightpink",
    "pink",
    # kappa
    "navajowhite",
    # eta
    "blanchedalmond",
    # iota
    "moccasin",
    # lambda
    "wheat",
    # mu
    "lemonchiffon",
    # EU1
    "gold",
    # epsilon
    "bisque",
    # other
    "lightyellow",
    "beige",
    "cornsilk",
    "khaki",
    "palegoldenrod",
    "peachpuff",
    "mistyrose",
    "paleturquoise",
    "lightcyan",
    "snow",
    "linen",
    "lavenderblush"
]
variant_color_mapping = {name:color for name,color in zip(variant_names, variant_colors)}


# ================== common objects =====================

# stores all figures to allow clearing them from memory
class FigureContainer:
    def __init__(self, plt):
        self.all_figures = []
        self.plt = plt
        
    def clear_figures(self):
        for figure in self.all_figures:
            self.plt.close(figure)
        self.all_figures = []
    
    def add(self, figure):
        self.all_figures.append(figure)


# ==================== common functions ============================

# reads the model probabilities
def read_probs(path):
    with open(path, "r") as f:
        return [float(p.split(",")[1]) for p in f.readlines()[1:]]


# finds which model is the most likely
def find_which_model(name):
    models_path = DATA_PATH + r"\ABC_results\{0}\finalpopulation\posterior\{0}-modelprobabilities.csv".format(name)
    probs = read_probs(models_path)
    return np.argmax(probs)


# parses a directory or filename into the components country, date, number of samples & read depth
def parse_name(name):
    _, country, date, n, rd = name.split("_")
    year, month = date.split("-")
    date = f"{year}-{int(month):02}"
    return country, date, n, rd


# loads a json file as a dictionary
def load_json(filename):
    with open(filename, "r") as f:
        dictionary = json.load(f)
    return dictionary


# writes a dictionary to a json file
def write_json(filename, dictionary):
    with open(filename, "w") as f:
        json_object = json.dumps(dictionary, indent = 4) 
        f.write(json_object)


# functions used to predict the virus spread for 1 and 2 subclones
def calculate_frequency_1(s, t1, growth_rate, t_end):
    dt = t_end - t1
    ex1 = growth_rate * s * dt
    ex2 = -growth_rate * t1
    b1 = np.exp(ex1)
    b2 = np.exp(ex2)
    b3 = b1 * b2
    return b3 / (b3 + 1 - b2)

def calculate_frequency_2(s1, t1, s2, t2, growth_rate, t_end):
    dt1 = t_end - t1
    dt2 = t_end - t2
    ex11 = growth_rate * s1 * dt1
    ex12 = growth_rate * s2 * dt2
    ex21 = -growth_rate * t1
    ex22 = -growth_rate * t2
    b1 = np.exp(ex11) * np.exp(ex21)
    b2 = np.exp(ex12) * np.exp(ex22)
    return b1 / (b1 + b2 + 1), b2 / (b1 + b2 + 1)


# get all model probabilities for one or more filtering values and with standard or special birth rate
def get_all_probs(filtering = 0.001, data = None, n_max = 2000, birth_rates = None):
    template_path = DATA_PATH + r"\ABC_results\{0}\finalpopulation\posterior\{0}-modelprobabilities.csv"
    if data is None:
        data = dict()
    if type(filtering) == type(0.01):
        filtering = [filtering]
    elif type(filtering) == type([]):
        filtering = filtering
    else:
        raise Exception()
    for f in filtering:
        if birth_rates:
            for b in birth_rates:
                properties = f"-filtered{f}_Nmax{n_max}_b-log({b})"
                l = len(properties)
                for name in list(map(lambda x:x[:-l], filter(lambda x:x.endswith(properties), os.listdir(DATA_PATH + r"\ABC_results")))):
                    probs = read_probs(template_path.format(name + properties))
                    country, date, _, _ =parse_name(name)
                    if country not in data:
                        data[country] = {date:{b:{f:probs}}}
                    elif date not in data[country]:
                        data[country][date] = {b:{f:probs}}
                    elif b not in data[country][date]:
                        data[country][date][b] =  {f:probs}
                    else:
                        data[country][date][b][f] =  probs
        else:
            properties = f"-filtered{f}_Nmax{n_max}"
            l = len(properties)
            for name in list(map(lambda x:x[:-l], filter(lambda x:x.endswith(properties), os.listdir(DATA_PATH + r"\ABC_results")))):
                probs = read_probs(template_path.format(name + properties))
                country, date, _, _ = parse_name(name)
                if country not in data:
                    data[country] = {date:{f:probs}}
                elif date not in data[country]:
                    data[country][date] = {f:probs}
                else:
                    data[country][date][f] = probs
    return data


# retrieves the inferred parameters
def get_inferred_parameters(name, n_clones):
    params_template_path = DATA_PATH + r"\ABC_results\{0}\finalpopulation\posterior\{0}-parameters-clone{1}.csv"
    path = params_template_path.format(name, n_clones)
    param_table = pd.read_csv(path)
    x = {c:(param_table[c].mean(), param_table[c].median()) for c in param_table.columns}
    del x["weight"]
    del x["cellularity"]
    return x


# calculate the frequency of each variant in a country for each month
def get_freqs(case_data, country):
    cases = case_data[country].copy()
    week_mapping = dict()
    total_sequence_counts = np.array(cases.pop("total_sequences"))
    for i, week in enumerate(cases.pop("week")):
        week_mapping[i] = week[:7] # extract year & month (yyyy-mm)
    variants = {month:{variant:0 for variant in cases.keys()} for month in week_mapping.values()}
    for i, count in enumerate(total_sequence_counts):
        month = week_mapping[i]
        variants[month]["count"] = variants[month].get("count", 0) + count
    for variant, counts in cases.items():
        for i, count in enumerate(counts):
            month = week_mapping[i]
            variants[month][variant] += count
    freqs = dict()
    for month, values in variants.items():
        total_count = values.pop("count")
        frequencies = {variant:var_count/total_count for variant, var_count in values.items()}
        freqs[month] = frequencies
    return freqs