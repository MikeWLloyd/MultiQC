import logging
import re
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module parses summary statistics generated by versions >= v2.1.0 where
    the command line option `--new-summary` has been specified.

    Note that running HISAT2 without this option (and older versions)
    gives log output identical to Bowtie2. These logs are indistinguishable
    and summary statistics will appear in MultiQC reports labelled as Bowtie2.
    See GitHub issues on the [HISAT2 repository](https://github.com/infphilo/hisat2/issues/48)
    and the [MultiQC repository](https://github.com/MultiQC/MultiQC/issues/221)
    for more information.

    HISAT2 does not report the input file names in the log, so MultiQC
    takes the filename as the sample. Note that if you specify
    `--summary-file` when running HISAT2 the same summary output
    appears both there and in the `stdout`. So if you save both with
    different names you may end up with duplicate samples in your
    MultiQC report.
    """

    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="HISAT2",
            anchor="hisat2",
            href="https://ccb.jhu.edu/software/hisat2/",
            info="Maps DNA or RNA reads against a genome or a population of genomes",
            doi=["10.1038/nmeth.3317", "10.1038/s41587-019-0201-4"],
        )

        # Find and load any HISAT2 reports
        self.hisat2_data: Dict = dict()
        for f in self.find_log_files("hisat2", filehandles=True):
            self.parse_hisat2_logs(f)

        # Filter to strip out ignored sample names
        self.hisat2_data = self.ignore_samples(self.hisat2_data)

        if len(self.hisat2_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.hisat2_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file
        self.write_data_file(self.hisat2_data, "multiqc_hisat2")

        # Basic Stats Table
        # Report table is immutable, so just updating it works
        self.hisat2_general_stats_table()

        # Alignment Rate Plot
        self.hisat2_alignment_plot()

    def parse_hisat2_logs(self, f):
        """
        Parse statistics generated by HISAT2 >= v2.1.0 that has been run
        with the --new-summary option. Older versions or logs from runs without
        that option are identical to that from bowtie2 and will be parsed
        by that module.
        """

        # Regexes
        regexes = {
            "unpaired_total": r"Total(?: unpaired)? reads: (\d+)",
            "unpaired_aligned_none": r"Aligned 0 times?: (\d+) \([\d\.]+%\)",
            "unpaired_aligned_one": r"Aligned 1 time: (\d+) \([\d\.]+%\)",
            "unpaired_aligned_multi": r"Aligned >1 times: (\d+) \([\d\.]+%\)",
            "paired_total": r"Total pairs: (\d+)",
            "paired_aligned_none": r"Aligned concordantly or discordantly 0 time: (\d+) \([\d\.]+%\)",
            "paired_aligned_one": r"Aligned concordantly 1 time: (\d+) \([\d\.]+%\)",
            "paired_aligned_multi": r"Aligned concordantly >1 times: (\d+) \([\d\.]+%\)",
            "paired_aligned_discord_one": r"Aligned discordantly 1 time: (\d+) \([\d\.]+%\)",
        }

        # Go through log file line by line
        s_name = f["s_name"]
        parsed_data = {}

        for line in f["f"]:
            # Attempt in vain to find original hisat2 command, logged by another program
            hscmd = re.search(r"hisat2 .+ -[1U] ([^\s,]+)", line)
            if hscmd:
                s_name = self.clean_s_name(hscmd.group(1), f)
                log.debug(f"Found a HISAT2 command, updating sample name to '{s_name}'")

            # Run through all regexes
            for k, r in regexes.items():
                match = re.search(r, line)
                if match:
                    parsed_data[k] = int(match.group(1))

            # Overall alignment rate
            overall = re.search(r"Overall alignment rate: ([\d\.]+)%", line)
            if overall:
                parsed_data["overall_alignment_rate"] = float(overall.group(1))

                # Save parsed data
                if s_name in self.hisat2_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, s_name)
                self.hisat2_data[s_name] = parsed_data

                # Reset in case we find more in this log file
                s_name = f["s_name"]
                parsed_data = {}

    def hisat2_general_stats_table(self):
        """Take the parsed stats from the HISAT2 report and add it to the
        basic stats table at the top of the report"""

        headers = {
            "overall_alignment_rate": {
                "title": "% Aligned",
                "description": "overall alignment rate",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "YlGn",
            }
        }
        self.general_stats_addcols(self.hisat2_data, headers)

    def hisat2_alignment_plot(self):
        # Split the data into SE and PE
        sedata = {}
        pedata = {}
        for s_name, data in self.hisat2_data.items():
            if "paired_total" in data:
                # Save half 'pairs' of mate counts
                m_keys = ["unpaired_total", "unpaired_aligned_none", "unpaired_aligned_one", "unpaired_aligned_multi"]
                for k in m_keys:
                    if k in data:
                        data[k] = float(data[k]) / 2.0
                pedata[s_name] = data
            else:
                sedata[s_name] = data

        # Two plots, don't mix SE with PE
        if len(sedata) > 0:
            sekeys = {
                "unpaired_aligned_one": {"color": "#20568f", "name": "SE mapped uniquely"},
                "unpaired_aligned_multi": {"color": "#f7a35c", "name": "SE multimapped"},
                "unpaired_aligned_none": {"color": "#981919", "name": "SE not aligned"},
            }
            pconfig = {
                "id": "hisat2_se_plot",
                "title": "HISAT2: SE Alignment Scores",
                "ylab": "# Reads",
                "cpswitch_counts_label": "Number of Reads",
            }
            self.add_section(plot=bargraph.plot(sedata, sekeys, pconfig))

        if len(pedata) > 0:
            pekeys = {
                "paired_aligned_one": {"color": "#20568f", "name": "PE mapped uniquely"},
                "paired_aligned_discord_one": {"color": "#5c94ca", "name": "PE mapped discordantly uniquely"},
                "unpaired_aligned_one": {"color": "#95ceff", "name": "PE one mate mapped uniquely"},
                "paired_aligned_multi": {"color": "#f7a35c", "name": "PE multimapped"},
                "unpaired_aligned_multi": {"color": "#ffeb75", "name": "PE one mate multimapped"},
                "unpaired_aligned_none": {"color": "#981919", "name": "PE neither mate aligned"},
            }
            pconfig = {
                "id": "hisat2_pe_plot",
                "title": "HISAT2: PE Alignment Scores",
                "ylab": "# Reads",
                "cpswitch_counts_label": "Number of Reads",
            }
            self.add_section(
                description="<em>Please note that single mate alignment counts are halved to tally with pair counts properly.</em>",
                plot=bargraph.plot(pedata, pekeys, pconfig),
            )
