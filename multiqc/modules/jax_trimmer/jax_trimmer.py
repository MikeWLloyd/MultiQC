""" MultiQC module to parse log output from the JAX custom FASTQ trimmer """


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
            name="JAX FASTQ Trimmer",
            anchor="jax_trimmer",
            href="https://bitbucket.org/jacksonlaboratory/pdx-nextflow-dsl2-conversion/src/main/bin/shared/filter_trim.py",
            info="a custom fastq trimmer.",
        )

        # Parse logs
        # Find and load any jax_trimmer reports
        self.jax_trimmer_data = dict()
        
        self.pe_bool = False

        for f in self.find_log_files("jax_trimmer"):
            self.parse_jax_trimmer_logs(f)

        # Filter to strip out ignored sample names
        self.jax_trimmer_data = self.ignore_samples(self.jax_trimmer_data)

        if len(self.jax_trimmer_data) == 0:
            raise UserWarning

        log.info("Found {} reports".format(len(self.jax_trimmer_data)))

        # Write parsed report data to a file
        self.write_data_file(self.jax_trimmer_data, "multiqc_jax_trimmer")

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.jax_trimmer_general_stats_table()

        # Alignment Rate Plot
        # self.jax_trimmer_stats_plot()
        ## NOTE: turn this on if you want a barplot of trim stats. Does not seem necessary. 

    def parse_jax_trimmer_logs(self, f):
        s_name = f["s_name"].rsplit("_", 1)[0] 

        parsed_data = {}
        regexes = {
            "perc_hq": r"Percentage of HQ reads\s+(\d+.\d+)%\s+(\d+.\d+)",
            "total_reads": r"Total number of reads\s+(\d+)\s+(\d+)",
            "total_hq_reads": r"Total number of HQ filtered reads\s+(\d+)\s+(\d+)",
            "reads_passing": r"Reads passing filter\s+(\d+)\s+(\d+)",
            "perc_passing": r"Percent reads passing filter\s+(\d+.\d+)%\s+(\d+.\d+)",
            "max_trim_len": r"Max Trimmed Length\s+((\d+(?:\.\d+)?))\s+((\d+(?:\.\d+)?))",
            "min_trim_len": r"Min Trimmed Length\s+((\d+(?:\.\d+)?))\s+((\d+(?:\.\d+)?))",
            "mean_trim_len": r"Mean Trimmed Length\s+((\d+(?:\.\d+)?))\s+((\d+(?:\.\d+)?))",
        }
        regexes_se = {
            "perc_hq": r"Percentage of HQ reads\s+(\d+.\d+)%",
            "total_reads": r"Total number of reads\s+(\d+)",
            "total_hq_reads": r"Total number of HQ filtered reads\s+(\d+)",
            "reads_passing": r"Reads passing filter\s+(\d+)",
            "perc_passing": r"Percent reads passing filter\s+(\d+.\d+)%",
            "max_trim_len": r"Max Trimmed Length\s+((\d+(?:\.\d+)?))",
            "min_trim_len": r"Min Trimmed Length\s+((\d+(?:\.\d+)?))",
            "mean_trim_len": r"Mean Trimmed Length\s+((\d+(?:\.\d+)?))",
        }

        self.pe_bool = False
        for l in f["f"].splitlines():
            if 'Read' in l:
                if 'Read 2' in l:
                    self.pe_bool = True
            else: 
                if self.pe_bool:
                    for k, r in regexes.items():
                        match = re.search(r, l)
                        if match:
                            parsed_data[k+"_1"] = float(match.group(1))
                            parsed_data[k+"_2"] = float(match.group(2))
                else:
                    for k, r in regexes_se.items():
                        match = re.search(r, l)
                        if match:
                            parsed_data[k+"_1"] = float(match.group(1))

        if len(parsed_data) > 0:
            if s_name in self.jax_trimmer_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.add_data_source(f, s_name)
            self.jax_trimmer_data[s_name] = parsed_data

    def jax_trimmer_general_stats_table(self):
        """Take the parsed stats from the Xenome log and add it to the
        basic stats table at the top of the report"""

        headers = OrderedDict()
        headers["total_reads_1"] = {
            "title": "Total Reads R1",
            "description": "The total number of reads in the sample",
            "min": 0,
            "scale": False,
        }
        headers["total_hq_reads_1"] = {
            "title": "Total HQ Reads R1",
            "description": "The total number of high quality reads",
            "min": 0,
            "scale": False,
        }
        headers["reads_passing_1"] = {
            "title": "Reads Passing Filter R1",
            "description": "The total number of reads passing filter",
            "min": 0,
            "scale": False,
        }
        headers["min_trim_len_1"] = {
            "title": "Min Trim Length R1",
            "description": "The minimum trim length of reads",
            "min": 0,
            "scale": False,
        }
        headers["mean_trim_len_1"] = {
            "title": "Mean Trim Length R1",
            "description": "The average trim length of reads",
            "min": 0,
            "scale": False,
        }
        headers["max_trim_len_1"] = {
            "title": "Max Trim Length R1",
            "description": "The max trim length of reads",
            "min": 0,
            "scale": False,
        }

        if self.pe_bool:
            headers["total_reads_2"] = {
                "title": "Total Reads R2",
                "description": "The total number of reads in the sample",
                "min": 0,
                "scale": False,
            }
            headers["total_hq_reads_2"] = {
                "title": "Total HQ Reads R2",
                "description": "The total number of high quality reads",
                "min": 0,
                "scale": False,
            }
            headers["reads_passing_2"] = {
                "title": "Reads Passing Filter R2",
                "description": "The total number of reads passing filter",
                "min": 0,
                "scale": False,
            }
            headers["min_trim_len_2"] = {
                "title": "Min Trim Length R2",
                "description": "The minimum trim length of reads",
                "min": 0,
                "scale": False,
            }
            headers["mean_trim_len_2"] = {
                "title": "Mean Trim Length R2",
                "description": "The average trim length of reads",
                "min": 0,
                "scale": False,
            }
            headers["max_trim_len_2"] = {
                "title": "Max Trim Length R2",
                "description": "The max trim length of reads",
                "min": 0,
                "scale": False,
            }

        self.general_stats_addcols(self.jax_trimmer_data, headers)

    def jax_trimmer_stats_plot(self):
        """Make the HighCharts HTML to plot the alignment rates"""

        # Specify the order of the different possible categories
        keys = OrderedDict()
        keys["total_reads_1"] = {"color": "#000034", "name": "Total Reads R1"}
        keys["total_hq_reads_1"] = {"color": "#328AE2", "name": "Total HQ Reads R1"}
        keys["reads_passing_1"] = {"color": "#7F9DA7", "name": "Reads Passing Filter R1"}

        if self.pe_bool:
            keys["total_reads_2"] = {"color": "#C1B8AA", "name": "Total Reads R2"}
            keys["total_hq_reads_2"] = {"color": "#DDD4C6", "name": "Total HQ Reads R2"}
            keys["reads_passing_2"] = {"color": "#db6d00", "name": "Reads Passing Filter R2"}

        # Config for the plot
        config = {
            "id": "jax_trimmer_stats",
            "title": "Jax Trimmer: Summary Statistics",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        self.add_section(
            description="This plot shows JAX trimmer summary read counts",
            helptext="""
            Trim read counts. 
            """,
            plot=bargraph.plot(self.jax_trimmer_data, keys, config),
        )
