#! usr/bin/python

"""Alright, now we're tying it all together. The user simply runs this script, and this should automatically
load the configuration file located in Desktop\Cluster Detection\config.txt, and perform cluster detection
followed by proximal gene analysis."""

import os
import sys
import analyze_clusters
import cluster_detect
import time


def main():

    desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
    working_path = os.path.join(desktop_path, 'Cluster Detection')
    config_file = os.path.join(working_path, 'config.txt')

    # establish vars
    nbin_size = -1  # number of nts per bin in the initial SOR read frequency histogram
    cbin_size = -1 # number of nts I want each cluster bin to be
    cluster_max_sep = -1  # furthest away two bins can be when considering inversion sites
    cluster_min_sep = -1  # closest two bins can be when considering inversion sites
    cluster_bin_cutoff_perc = -1  # eliminate data below this percentile in the cluster analysis
    ntpair_min_sep = -1
    ntpair_max_sep = -1
    max_genes = -1

    # load config file
    try:
        with open(config_file, 'r') as f:
            header = f.readline()
            print("Loaded:", header)
            f.readline(), f.readline(), f.readline()
            line = ''

            while '!' not in line:

                line = f.readline()

                try:
                    label = line.split('=')[0]
                    value = line.split('=')[1]

                    if label == 'nbin_size':
                        nbin_size = int(value)

                    elif label == 'cbin_size':
                        cbin_size = int(value)
                    elif label == 'cluster_max_sep':
                        cluster_max_sep = int(value)
                    elif label == 'cluster_min_sep':
                        cluster_min_sep = int(value)
                    elif label == 'cbin_cutoff':
                        cluster_bin_cutoff_perc = int(value)
                    elif label == 'ntpair_min_sep':
                        ntpair_min_sep = int(value)
                    elif label == 'ntpair_max_sep':
                        ntpair_max_sep = int(value)
                    elif label == 'max_genes':
                        max_genes = int(value)
                except IndexError:
                    pass

    except FileNotFoundError as e:
        sys.exit("Config file not found: {0}".format(e.errno))

    # do we have any bad vars?
    params = (nbin_size, cbin_size, cluster_min_sep, cluster_max_sep, cluster_bin_cutoff_perc, ntpair_min_sep,
              ntpair_max_sep, max_genes)

    print(params)

    for param in params:
        if param == -1:
            sys.exit("Warning! {0} not found in config file. Please ensure variable is set.".format(param))

    # okay! literally two lines here

    print("Detecting clusters....")
    cluster_detect.detect_inversion_clusters(nbin_size, cluster_bin_cutoff_perc, cbin_size, cluster_min_sep,
                                             cluster_max_sep, ntpair_min_sep, ntpair_max_sep)

    print("Analyzing clusters...")
    analyze_clusters.align_clusters_to_genes(max_genes)

    print("Done!")

    x = time.clock()
    print("Operation took", x, "seconds.")


if __name__ == "__main__":
    main()
