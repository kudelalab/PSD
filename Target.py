class Target:
    micron_factor = 1 / 2.7

    def __init__(self, sample, i, fea_file):
        self.sample = sample
        self.index = i
        self.biovolume = fea_file['Biovolume'][i] * pow(Target.micron_factor, 3)
        self.equiv_diameter = fea_file['EquivDiameter'][i] * Target.micron_factor
        self.major_axis_length = fea_file['MajorAxisLength'][i] * Target.micron_factor
        self.minor_axis_length = fea_file['MinorAxisLength'][i] * Target.micron_factor
        self.json = {
            "index": i,
            "biovolume": self.biovolume,
            "equiv_diameter": self.equiv_diameter,
            "major_axis_length": self.major_axis_length,
            "minor_axis_length": self.minor_axis_length
        }