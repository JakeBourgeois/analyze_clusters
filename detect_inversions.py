"""
Insert useful information here...
"""

#IMPORTS

import seaborn as sns
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import csv


# loads SOR data as a frequency plot
def get_sor_freq_data(sor_file, dump_freqs=('n', 0)):
    pos_freq_dict = defaultdict(int)
    all_positions = list()

    with open(sor_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:

            # ignore TLEN values that are less than zero. Does this need to be within a read length?
            if int(row['TLEN']) != 0:
                pos_freq_dict[int(row['POS'])] += 1
                all_positions.append(int(row['POS']))

    if dump_freqs[0] == 'y':

        with open(dump_freqs[1], 'w') as d:
            writer = csv.writer(d, delimiter=',')
            writer.writerow(("GenomicPos", "Frequency"))

            for key in pos_freq_dict:
                writer.writerow((key, pos_freq_dict[key]))

    return pos_freq_dict, np.array(all_positions)


# makes a histogram array using numpy given a frequency dictionary and data bounds
def make_histogram_array(pos_freq_dict, datamin, datamax):
    a = list()
    for pos in pos_freq_dict:
        if (pos <= datamax) and (pos >= datamin):
            for i in range(0, pos_freq_dict[pos] - 1):
                a.append(pos)
    return np.array(a)


# filters an array based on data constraints
def filter_array(arr, dmin, dmax):
    f = list()
    for i in arr:
        if (i >= dmin) and (i <= dmax):
            f.append(i)
    return f


# draws an histogram of SOR data to get some threshold for a positive bin
def filter_sor_histogram(sor_array, ntbinsize=20000):
    # class HLineBuilder allows us to define a density cutoff in the initial screen.
    class HLineBuilder:
        def __init__(self, line, x_bin):
            self.line = line
            self.xs = (0, x_bin)
            self.ys = line.get_ydata()
            self.cid = line.figure.canvas.mpl_connect('button_press_event', self)

        def __call__(self, event):
            if event.inaxes != self.line.axes: return
            y1, y2 = event.ydata, event.ydata
            self.line.set_data(self.xs, (y1, y2))
            self.y_final = y1
            self.line.figure.canvas.draw()

    seq_size = sor_array.max() - sor_array.min()
    nbins = int(seq_size / ntbinsize)

    h_densities, den_bin_edges = np.histogram(sor_array, bins=nbins, density=True)

    # Plot the density histogram, providing visual representation of read densities
    fig, ax1 = plt.subplots()
    ax1.plot(h_densities)  # plot density histogram along axis.

    # Call a linebuilder class to allow the user to draw a cutoff visually.
    line, = ax1.plot([0], [0])  # empty line
    r = HLineBuilder(line, nbins)

    # Define plot parameters
    plt.title("Click to set read density cutoff")
    plt.xlabel("Bin")
    plt.ylabel("Bin read density")

    plt.show()

    # Our cutoff is equal to the y value of the last line drawn
    read_cutoff = r.y_final

    # add the left-sided bin edges to a list if they pass; these represent the left side of a potential cluster
    h_bin_left_pos_list = list()
    for i in range(0, len(h_densities)):
        if h_densities[i] > read_cutoff:
            h_bin_left_pos_list.append(float(den_bin_edges[i]))

    print("Screened", len(h_bin_left_pos_list), "potential clusters.")

    return h_bin_left_pos_list, read_cutoff


# generates all unique pairs of bins in a cluster
def generate_cluster_bin_pairs(cluster_pos_dict):
    print("Generating potential bin pairs...")
    cluster_bin_pairs = list()
    pass_pos_list = list()
    for pos in cluster_pos_dict:
        pass_pos_list.append(pos)

    for i in range(0, len(pass_pos_list) - 1):
        this_bin_pos = pass_pos_list[i]
        rest_bin_pos = pass_pos_list[i + 1:]
        for pos in rest_bin_pos:
            cluster_bin_pairs.append((this_bin_pos, pos))

    print("Bin pairs generated:", len(cluster_bin_pairs))

    return cluster_bin_pairs


# filters the initial cluster dictionary by read count percentile
def filter_by_read(cpos_dict, cperc):
    print("Eliminating bins that are below", cperc, "counts...")

    i = 0
    to_remove = list()
    for pos in cpos_dict:
        read_count = cpos_dict[pos]
        if read_count <= cperc:
            to_remove.append(pos)
            i += 1

    for pos in to_remove:
        cpos_dict.pop(pos)
    print("Eliminated", i, "bins.")

    return cpos_dict


# filters pairs based on location to each other
def filter_cluster_bin_pairs(cluster_pairs, dmin, dmax):
    print("Eliminating bin pairs that are further than", dmax, "nt or closer than", dmin, "nt from each other...")
    i = 0
    for bin_pair in cluster_pairs:
        bin1 = bin_pair[0]
        bin2 = bin_pair[1]
        if abs(bin2 - bin1) >= 6000.0 or abs(bin2 - bin1) <= 50.0:
            cluster_pairs.remove(bin_pair)
            i += 1
            # print("Removed:", bin_pair)
    print("Eliminated", i, "bins.")

    return cluster_pairs


# finds the maximally scoring pair of cluster bins
def find_max_pair(pairs, read_dict):
    print("Finding the bin pair that maximizes the reads...")
    bin_count_max, bin_max_pair = 0, (-1, -1)
    for bin_pair in pairs:
        bin1 = bin_pair[0]
        bin2 = bin_pair[1]
        read_count = read_dict[bin1] + read_dict[bin2]

        if read_count > bin_count_max:
            bin_count_max = read_count
            bin_max_pair = bin_pair

    print("Bin max pair:", bin_max_pair)
    return bin_max_pair


# finds the maximally scoring nucleotide
def find_best_nucleotide(pos_list, pos_freq_dict):

    best_nt = -1
    read_max = 0

    for pos in pos_list:
        if pos_freq_dict[pos] > read_max:
            best_nt = pos
            read_max = pos_freq_dict[pos]

    return best_nt, read_max


def detect_inversion_nt_positions(cluster_pos_start, cluster_pos_end, sor_pos_freqs, sor_pos_arr, cbin_size, cperc,
                                  cluster_min_sep, cluster_max_sep):

    # generate a histogram array of positions that lie within the cluster bounds
    cluster_pos_array = make_histogram_array(sor_pos_freqs, datamin=cluster_pos_start, datamax=cluster_pos_end)

    # make a histogram of this data based on cbin size
    cbins = 100
    cluster_histogram_counts, cluster_histogram_edges = np.histogram(cluster_pos_array, bins=cbins)

    # let's link the arrays with a dictionary of bin edges to read counts for easier manipulation later.
    # [{bin pos:read}]
    cluster_read_pos_dict = dict()
    print("Generating dictionary of cluster bin locations to cluster bin counts...")
    for i in range(0, len(cluster_histogram_counts)):
        cluster_read_pos_dict[cluster_histogram_edges[i]] = cluster_histogram_counts[i]

    # filter the bins here based off count percentile
    cluster_bin_cutoff_perc_val = np.percentile(cluster_histogram_counts, cperc)
    cluster_dict_read_screened = filter_by_read(cluster_read_pos_dict, cluster_bin_cutoff_perc_val)

    # we may have eliminated all of our bins!
    if len(cluster_dict_read_screened) == 0:
        print("All bins eliminated! False region identified. Skipping...\n")

    else:

        # generate cluster pairs
        cluster_bin_pairs = generate_cluster_bin_pairs(cluster_dict_read_screened)

        # filter cluster pairs based on proximity
        cluster_bin_pairs_screened = filter_cluster_bin_pairs(cluster_bin_pairs,
                                                              dmin=cluster_min_sep, dmax=cluster_max_sep)

        # now let's consider some cases...

        # 1. There are literally no bin pairs...likely nt spike
        if len(cluster_bin_pairs_screened) == 0:
            print("No bin pairs remaining! There is likely a single nucleotide signal.")
            arr1 = filter_array(sor_pos_arr, cluster_histogram_edges.min(), cluster_histogram_edges.max())
            sig_nt = find_best_nucleotide(arr1, sor_pos_freqs)
            return sig_nt, sig_nt

        # 2. We have bin pairs!
        else:

            # now let's find the best pair!
            max_cluster_pair = find_max_pair(cluster_bin_pairs_screened, cluster_dict_read_screened)

            # now let's find the best nucleotides in these pairs!
            inversion_nucleotides, inversion_nucleotides_read_scores = [-1, -1], [-1, -1]

            # generate arrays
            cluster_size = cluster_histogram_edges[2] - cluster_histogram_edges[1]
            c1_lb, c2_lb = max_cluster_pair[0] - 0.1 * cluster_size, max_cluster_pair[1] - 0.1 * cluster_size
            c1_ub, c2_ub = max_cluster_pair[0] + 1.1 * cluster_size, max_cluster_pair[1] + 1.1 * cluster_size

            arr1 = np.unique(filter_array(sor_pos_arr, dmin=c1_lb, dmax=c1_ub))
            arr2 = np.unique(filter_array(sor_pos_arr, dmin=c2_lb, dmax=c2_ub))

            # find the best nucleotides
            inversion_nucleotides[0], inversion_nucleotides_read_scores[0] = \
                find_best_nucleotide(arr1, pos_freq_dict=sor_pos_freqs)
            inversion_nucleotides[1], inversion_nucleotides_read_scores[1] = \
                find_best_nucleotide(arr2, pos_freq_dict=sor_pos_freqs)

            # now, is one of the nucleotides scoring nearly 100% of the data?
            per_dif = 100*abs((inversion_nucleotides_read_scores[1] - inversion_nucleotides_read_scores[0])) / \
                (inversion_nucleotides_read_scores[1] + inversion_nucleotides_read_scores[0])

            if per_dif > 95.0:

                print("Warning! One nucleotide exceeds 95% of read data.")
                dif = inversion_nucleotides_read_scores[1] - inversion_nucleotides_read_scores[0]
                if dif > 0.0:
                    return inversion_nucleotides[1], inversion_nucleotides[1]
                else:
                    return inversion_nucleotides[0], inversion_nucleotides[0]

            # if not, let's return the max pair
            else:
                return inversion_nucleotides[0], inversion_nucleotides[1]





