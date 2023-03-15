""" MultiQC module to parse log output from Primerclip """


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
            name="Primerclip",
            anchor="primerclip",
            href="https://github.com/swiftbiosciences/primerclip",
            info="Primerclip™ is an alignment-based primer trimming tool designed to trim primer sequences for Swift Biosciences Accel-Amplicon™ panels.",
        )

        # Parse logs
        # Find and load any Bowtie reports
        self.primerclip_data = dict()
        for f in self.find_log_files("primerclip"):
            self.parse_primerclip_logs(f)

        # Filter to strip out ignored sample names
        self.primerclip_data = self.ignore_samples(self.primerclip_data)

        if len(self.primerclip_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.primerclip_data)))

        # Write parsed report data to a file
        self.write_data_file(self.primerclip_data, "multiqc_primerclip")

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.primerclip_general_stats_table()

    def parse_primerclip_logs(self, f):
        s_name = f["s_name"].replace("_primerclip_runstats", "")
        parsed_data = {}
        regexes = {
            "total_alignments": r"Total alignments processed:\s+(\d+)",
            "total_mapped_alignments": r"Total mapped alignments:\s+(\d+)",
            "trimmed_by_ge_one_base": r"^Alignments trimmed by >= 1 base:\s+(\d+)",
            "trimmed_to_zero": r"Alignments trimmed to zero aligned length:\s+(\d+.\d+)",
            "perc_trimmed_by_ge_one_base": r"% Alignments trimmed by >= 1 base:\s+(\d+.\d+)",
            "perc_mapped_after_trimming": r"% Alignments mapped after trimming:\s+(\d+.\d+)"
        }

        for l in f["f"].splitlines():
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    parsed_data[k] = float(match.group(1))

        if len(parsed_data) > 0:
            if s_name in self.primerclip_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.primerclip_data[s_name] = parsed_data

    def primerclip_general_stats_table(self):
        """Take the parsed stats from the primerclip log and add it to the
        basic stats table at the top of the report"""

        headers = OrderedDict()
        headers["total_alignments"] = {
            "title": "Total Alignments",
            "description": "The number total alignments in sample",
            "min": 0,
            "scale": False,
        }
        headers["total_mapped_alignments"] = {
            "title": "Mapped Alignments",
            "description": "The number of mapped alignments in the sample",
            "min": 0,
            "scale": False,
        }
        headers["trimmed_by_ge_one_base"] = {
            "title": "Trimmed Alignments",
            "description": "Alignments trimmed by >= 1 base",
            "min": 0,
            "scale": False,
        }
        headers["trimmed_to_zero"] = {
            "title": "Removed Alignments",
            "description": "Alignments trimmed to 0 length",
            "min": 0,
            "scale": False,
        }
        self.general_stats_addcols(self.primerclip_data, headers)
