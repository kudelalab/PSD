import numpy as np
import os
import pandas as pd
from matplotlib import pyplot as plt
from math import floor, log10
from scipy.optimize import curve_fit
import datetime as dt
from Target import Target

class Sample:

    def __init__(self, name, feature_dir, roi_dir, overall_bin, ifcb):
        self.name = name
        self.ifcb = ifcb
        self.micron_factor = 1 / 2.7
        self.summary = {}
        self.data = []
        self.psd = {}
        self.coeffs = {}
        self.r_squared = {}
        self.max_diff = 0
        self.bin = overall_bin

        def getMATLabDate():
            datetime = dt.datetime.strptime(self.name, 'D%Y%m%dT%H%M%S')
            delta = (datetime - dt.datetime(1, 1, 1))
            datenum = delta.total_seconds() / (24 * 60 * 60) + 367
            return datenum

        self.datenum = getMATLabDate()

        print(f'Processing {ifcb}')

        self.features = pd.read_csv(f'{feature_dir}/{name}_{ifcb}_features.csv')
        self.metadata = self.read_metadata(roi_dir, name)
        self.targets = [Target(self, i, self.features) for i in range(len(self.features.Biovolume))]
        self.mL_analyzed = self.get_volume()
        self.capture_percent = self.get_capture_percent(roi_dir, name)
        self.humidity = float(self.metadata['humidity'][0])
        self.bead_run = self.metadata.get('runType', ['NORMAL'])[0] == 'BEADS'

        self.grouped_equiv_diameter = {}
        self.grouped_major_axis_length = {}
        self.grouped_minor_axis_length = {}

        self.grouped_equiv_diameter_export = []
        self.grouped_major_axis_length_export = []
        self.grouped_minor_axis_length_export = []

        self.grouped_equiv_diameter_json = []
        self.grouped_major_axis_length_json = []
        self.grouped_minor_axis_length_json = []

    def to_JSON(self):
        return {
            "name": self.name,
            "ifcb": self.ifcb,
            "micron_factor": self.micron_factor,
            "summary": self.summary,
            "data": self.data,
            "psd": self.psd,
            "coeffs": self.coeffs,
            "r_squared": self.r_squared,
            "max_diff": self.max_diff,
            "bin": {},
            "datenum": self.datenum,
            "features": {},
            "metadata": {},
            "targets": [t.json for t in self.targets],
            "mL_analyzed": self.mL_analyzed,
            "grouped_equiv_diameter": self.grouped_equiv_diameter_json,
            "grouped_major_axis": self.grouped_major_axis_length_json,
            "grouped_minor_axis": self.grouped_minor_axis_length_json,
        }

    def read_metadata(self, roi_dir, name):
        '''Extracts metadata from header file'''

        metadata = {}
        with open(f'{roi_dir}/{name[:9]}/{name}_{self.ifcb}.hdr', 'r') as hdr:
            raw_metadata = [line.strip() for line in hdr.readlines()]
            for line in raw_metadata:
                split = line.split(': ')
                if len(split) < 2:
                    continue
                try:
                    metadata[split[0]] = [float(v) for v in split[1].split(',')]
                except:
                    metadata[split[0]] = split[1].split(',')

        return metadata

    def get_volume(self):
        '''Determines the volume analyzed in the sample in mL'''

        flow_rate = 0.25
        looktime = self.metadata['runTime'][0] - self.metadata['inhibitTime'][0]

        return flow_rate * looktime / 60

    def get_capture_percent(self, roi_dir, name):
        '''Determines the ratio of triggers to images'''

        with open(f'{roi_dir}/{name[:9]}/{name}_{self.ifcb}.adc', 'r') as adc:
            trigger_count = 0
            for line in adc:
                trigger_count += 1

        return len(self.targets) / trigger_count

    def group(self, feature):
        '''Creates a dictionary representation for a histogram of the targets using some length feature'''

        groups = {num: [] for num in range(200)}
        json_groups = {num: [] for num in range(200)}
        for target in self.targets:
            diameter_group = floor(getattr(target, feature))
            if diameter_group > 199:
                diameter_group = 199
            groups[diameter_group] += [target]
            json_groups[diameter_group] += [target.json]

        setattr(self, f'grouped_{feature}', groups)
        setattr(self, f'grouped_{feature}_json', json_groups)
        return groups

    def export_specific_groups(self, feature):
        return [(len(targets) / self.mL_analyzed) * 1000 for targets in getattr(self, f'grouped_{feature}').values()]

    def export_groups(self):
        return self.data

    def create_histograms(self):
        self.group('equiv_diameter')
        self.group('major_axis_length')
        self.group('minor_axis_length')
        self.psd = self.histogram('equiv_diameter')

    def histogram(self, feature):
        '''Creates a dictionary representation for a histogram of biovolume using some length feature'''

        groups = getattr(self, f'grouped_{feature}')
        histogram = {num: 0 for num in range(200)}

        for diameter, targets in groups.items():
            if targets and self.mL_analyzed:
                histogram[diameter] = (len(targets) / (self.mL_analyzed)) * 1000

        return histogram

    def plot_PSD(self, use_marker, save_graph, start_fit):

        print(f'Graphing {self.name}')

        def power_curve(x, k, n):
            return k * (x ** n)

        def round_sig(x, sig=2):
            try:
                return round(x, sig - int(floor(log10(abs(x)))) - 1)
            except:
                return 0

        xdata, ydata = zip(*self.psd.items())
        maximum = max(ydata)
        max_diff = start_fit - ydata.index(maximum)

        if save_graph:
            fig, ax = plt.subplots()
            ax.set(xlabel="ESD [um]", ylabel="N'(D) [c/L⁻]")
            ax.set_ylim(bottom=-0.1 * maximum, top=1.1 * maximum)

        try:
            popt, pcov, infodict, mesg, ier = curve_fit(power_curve, xdata[start_fit:], ydata[start_fit:],
                                                        full_output=True, p0=[80000, -0.8])
            residuals = ydata[start_fit:] - power_curve(xdata[start_fit:], *popt)
            ss_res = np.sum(residuals ** 2)
            ss_tot = np.sum((ydata[start_fit:] - np.mean(ydata[start_fit:])) ** 2)
            r_sqr = 1 - (ss_res / ss_tot)
            self.bin.add_fit(self.name, round_sig(popt[0], 5), round_sig(popt[1], 5), r_sqr, max_diff, self.capture_percent, self.bead_run, self.humidity)
        except:
            popt = [0.0, 0.0]
            r_sqr = 0
            self.bin.add_fit(self.name, 0,0, r_sqr, max_diff, self.capture_percent, self.bead_run, self.humidity)

        self.bin.add_data(self.name, self.datenum, ydata, self.mL_analyzed, maximum)

        if use_marker:
            marker = 'o'
        else:
            marker = None

        if save_graph:
            psd_line = ax.plot(xdata[start_fit:], xdata[start_fit:],
                               color='#00afbf', marker=marker, linestyle='solid',
                               linewidth=1.25, markersize=4, label='PSD')
            if r_sqr > 0:
                curve_fit_line = ax.plot(xdata, power_curve(xdata[start_fit:], *popt), color='#516b6e', linestyle='dashed',
                                         label='Power Curve')
                ax.text(80, maximum * 0.75,
                        f'$y = ({round_sig(popt[0], 3)})x^{{{round_sig(popt[1], 3)}}}$, $R^{{2}} = {round_sig(r_sqr, 3)}$')

            ax.legend()
            ax.set_title(f'{self.name}')
            plt.savefig(os.path.join(os.getcwd(), 'Graphs', self.name))
            plt.close('all')