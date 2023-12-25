import os
import pandas as pd
import re
import json
import operator
from Sample import Sample

class Bin:
    def __init__(self, feature_dir, hdr_dir, samples_path=None):

        fileConvention = r'D\d\d\d\d\d\d\d\dT\d\d\d\d\d\d'
        regex = re.compile(fileConvention)
        files = [(f.split('_')[0], f.split('_')[1]) for f in os.listdir(feature_dir) if regex.search(f)]

        self.samples_loaded = bool(samples_path)
        if self.samples_loaded:
            sample_file = open(samples_path, 'rb')
            samples = json.load(sample_file)["samples"]
            self.samples = [json.loads(json.dumps(s), object_hook=lambda d: Sample(**d)) for s in samples]

        else:
            self.samples = [Sample(f[0], feature_dir, hdr_dir, self, f[1]) for f in files]

        self.file_names = [s.name for s in self.samples]
        self.data = pd.DataFrame(columns=['mL_analyzed', 'max'] + [f'{i}μm' for i in range(0, 200)], index=self.file_names)
        self.fits = pd.DataFrame(columns=['a', 'k', 'R^2', 'max_ESD_diff', 'capture_percent', 'bead_run', 'humidity'],
                                 index=self.file_names)

    def add_data(self, file, datenum, data, mL_analyzed, maximum):

        formatted_data = {f'{i}μm': data[i] for i in range(0, 200)}
        formatted_data['mL_analyzed'] = mL_analyzed
        formatted_data['max'] = maximum
        formatted_data['datenum'] = datenum
        self.data.loc[file] = pd.Series(formatted_data)

    def add_fit(self, file, a, k, r_sqr, max_diff, capture_percent, bead_run, humidity):

        self.fits.loc[file] = pd.Series({'a': a, 'k': k, 'R^2': r_sqr, 'max_ESD_diff': max_diff,
                                         'capture_percent': capture_percent, 'bead_run': bead_run, 'humidity': humidity})

    def pick_start(self):
        files = self.file_names[:]
        months = [[] for m in range(12)]
        for m in range(1, 13):
            month = str(m).zfill(2)
            currFile = files[0]
            while currFile[5:7] == month:
                months[m - 1] += [currFile]
                files = files[1:]
        print(files)
        print()

    def plot_PSD(self, use_marker, save_graphs, start_fit):

        if save_graphs:
            os.mkdir(os.path.join(os.getcwd(), 'Graphs'))
        if not self.samples_loaded:
            for sample in self.samples:
                sample.create_histograms()
        for sample in self.samples:
            sample.plot_PSD(use_marker=use_marker, save_graph=save_graphs, start_fit=start_fit)
        print(f'Start fit: {start_fit}')

    def save_data(self, name, r_sqr=0.5, **kwargs):

        print(f'Saving Data')

        def flag(dataset, op, parameter, threshold, flag_name, priority, low_r_only):
            if low_r_only:
                dataset = dataset[self.fits['R^2'] < r_sqr]

            if type(op) != tuple:
                files = dataset[op(dataset[parameter], threshold)]
            else:
                files = dataset[op[0](dataset[parameter[0]], threshold[0])]
                for i in range(1, len(op)):
                    files = files[op[i](files[parameter[i]], threshold[i])]

            if flag_name == 'Beads':
                calculated_df = pd.DataFrame({'file': list(files.index), 'flag': [flag_name] * len(files), 'priority': priority})
                bead_runs = self.fits[self.fits['bead_run']]
                set_df = pd.DataFrame({'file': list(bead_runs.index), 'flag': [flag_name] * len(bead_runs), 'priority': priority})
                return pd.concat([calculated_df, set_df]).drop_duplicates('file')

            return pd.DataFrame({'file': list(files.index), 'flag': [flag_name] * len(files), 'priority': priority})

        flag_params = {
            'beads': (self.fits, operator.gt, 'a', 'Beads', 1),
            'bubbles': (self.data, operator.lt, 'max', 'Bubbles', 2),
            'incomplete': (self.data, (operator.lt, operator.lt), ('max', 'mL_analyzed'), 'Incomplete Run', 3),
            'missing_cells': (self.fits, operator.lt, 'capture_percent', 'Missing Cells', 4),
            'biomass': (self.data, operator.lt, 'max', 'Low Biomass', 5),
            'bloom': (self.fits, operator.gt, 'max_ESD_diff', 'Bloom', 6),
            'humidity': (self.fits, operator.gt, 'humidity', 'High Humidity', 7)
        }

        full_flags = pd.DataFrame({'file': [], 'flag': [], 'priority': 10000})
        r_limited_flags = ['biomass', 'bloom']
        esd_diff_flags = ['bloom']
        r_flag = flag(self.fits, operator.lt, 'R^2', r_sqr, 'Low R^2', 7, False)
        if len(r_flag['file']) > 0:
            full_flags = pd.concat([full_flags, r_flag], ignore_index=True)

        for key, value in kwargs.items():
            if key in esd_diff_flags:
                value = -value
            [dataset, op, parameter, flag_name, priority] = flag_params[key]
            key_flag = flag(dataset, op, parameter, value, flag_name, priority, key in r_limited_flags)
            if len(key_flag['file']) > 0:
                full_flags = pd.concat([full_flags, key_flag], ignore_index=True)

        flags = full_flags.sort_values('priority').drop_duplicates(subset=['file']).sort_values(by=['file'])
        flags = flags.drop('priority', axis=1)

        self.data.to_csv(f'{name}_data.csv')
        self.fits.to_csv(f'{name}_fits.csv')
        flags.to_csv(f'{name}_flags.csv')
