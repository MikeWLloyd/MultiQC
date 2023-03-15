""" MultiQC module to parse log output from coverage metric computation """


import logging
import re
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Coverage Metrics",
            anchor="coverage_metrics",
            info="Amplicon coverage metric summary.",
        )

        # Parse logs
        # Find and load any Bowtie reports
        self.coverage_metrics_data = dict()
        for f in self.find_log_files("coverage_metrics"):
            self.parse_coverage_metrics_logs(f)

        # Filter to strip out ignored sample names
        self.coverage_metrics_data = self.ignore_samples(self.coverage_metrics_data)

        if len(self.coverage_metrics_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.coverage_metrics_data)))

        # Write parsed report data to a file
        self.write_data_file(self.coverage_metrics_data, "multiqc_coverage_metrics")

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.coverage_metrics_general_stats_table()

    def parse_coverage_metrics_logs(self, f):
        s_name = f["s_name"].replace("_amplicon_coverage_metrics", "")
        parsed_data = {}
        regexes = {
            "on_target_percent": r"on_target_percent\s+(\d+(?:\.\d+)?)",
            "coverage_uniformity": r"coverage_uniformity\s+(\d+(?:\.\d+)?)",
        }

        for l in f["f"].splitlines():
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    parsed_data[k] = float(match.group(1))

        if len(parsed_data) > 0:
            if s_name in self.coverage_metrics_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.coverage_metrics_data[s_name] = parsed_data

    def coverage_metrics_general_stats_table(self):
        """Take the parsed stats from the coverage_metrics log and add it to the
        basic stats table at the top of the report"""

        headers = OrderedDict()
        headers["on_target_percent"] = {
            "title": "On Target Percentage",
            "description": "The percentage of total aligned reads sequenced aligned to the intended target regions.",
            "min": 0,
            "scale": False,
        }
        headers["coverage_uniformity"] = {
            "title": "Coverage Uniformity",
            "description": "Coverage uniformity percent (CU%) is calculated as the ratio of 'number of target bases that have coverage at or above 20% (0.2μ)' and 'total number of target bases'",
            "min": 0,
            "scale": False,
        }
        self.general_stats_addcols(self.coverage_metrics_data, headers)
