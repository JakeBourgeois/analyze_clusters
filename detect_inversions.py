"""
detect_inversions.py is meant to contain all the useful classes, methods, and functions necessary for
detection of inversions within an organism SOR file

Author: Jake Bourgeois
Email: jacob.bourgeois@tufts.edu
Affiliation: Tufts University, Camilli Lab
"""

#IMPORTS

import seaborn as sns
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
import csv
import sys


# class SOR holds the master SOR data
class SOR:

    def __init__(self, acc, sor_file, binsize=20000, ignore=[]):

        self.accession_num = acc                # accession number
        self.pos_freq_dict = defaultdict(int)   # contains pos:freq data
        self.pos_array = np.array([])           # contains all unique positions
        self.ignored_positions = ignore         # when loading the SOR file, ignore these positions
        self.data_sum = 0                       # sum of read counts in data

        self.load_sor(sor_file)                 # loads sor data into these attributes

        self.pos_min = self.pos_array.min()     # minimum position
        self.pos_max = self.pos_array.max()     # maximum position

        self.bin_size = binsize                 # how many nucleotides each bin should span
        self.final_bin_size = 0                 # what we ended up getting

        # useful output parameters
        self.clusters = list()                  # list of clusters filtered out of the initial screen
        self.signals = list()                   # list of all screened clusters
        self.true_clusters = list()             # signals that turn out to be a true pair of clusters
        self.spikes = list()                    # signals that turn out to be solitary spikes or otherwise
        self.clust_dist = list()                # list of distances between best nt
        self.clust_sum = list()                 # list of scores sums of best nt
        self.clust_proportion_perc = list()     # list of percent proportions of reads to cluster
        self.total_clust_prop_perc = list()     # list of percent proportions of reads to data total
        self.cluster_map_int = list()           # list of 1 and -1 reporting true or spikes; for output index mapping
        self.final_cbin_sizes = list()          # list of final cluster bin sizes
        self.read_cutoff = 0                    # ultimate density value used in thresholding

    # loads sor data from file into attributes
    def load_sor(self, sor_file):

        with open(sor_file, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    # ignore TLEN values that are less than zero.
                    if int(row['TLEN']) != 0:
                        if int(row['POS']) not in self.ignored_positions:
                            self.pos_freq_dict[int(row['POS'])] += 1
                            self.data_sum += 1

                except csv.Error as e:
                    print("Error occurred! Please ensure headers on file include TLEN and POS.")
                    sys.exit('file {}, line {}: {}'.format(sor_file, reader.line_num, e))
        self.pos_array = np.sort(np.array(list(self.pos_freq_dict)))
        return

    # returns a subset of the pos_freq_dict given a start and an end
    def subset(self, pos_start, pos_end):

        sub_dict = dict()
        for element in self.pos_array:
            if (element >= pos_start) and (element <= pos_end):
                sub_dict[element] = self.pos_freq_dict[element]
        return sub_dict

    # creates an interactive graph of histogram data useful for setting thresholds of cluster detection
    def make_interactive_graphical_threshold(self, save_path='n'):

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

        seq_size = self.pos_max - self.pos_min
        nbins = int(seq_size / self.bin_size)

        data = list()
        for pos in self.pos_array:
            for i in range(0, self.pos_freq_dict[pos]):
                data.append(pos)
        data = np.array(data)

        # make a frequency histogram

        h_densities, den_bin_edges = np.histogram(data, bins=nbins, density=True)
        self.final_bin_size = den_bin_edges[1] - den_bin_edges[0]

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
        self.read_cutoff = r.y_final
        plt.close()

        # If we have a save path, recreate the graph and save it
        if save_path != 'n':
            fig, ax1 = plt.subplots()
            ax1.plot(h_densities)  # plot density histogram along axis.
            plt.title(self.accession_num + " SOR Density Histogram; Bins=" + str(nbins))
            plt.xlabel("Bin Number")
            plt.ylabel("Bin read density")
            plt.axhline(self.read_cutoff, color='r')
            plt.savefig(save_path)
            plt.close()

        # add the left-sided bin edges to a list if they pass; these represent the left side of a potential cluster
        h_bin_left_pos_list = list()
        for i in range(0, len(h_densities)):
            if h_densities[i] >= self.read_cutoff:
                h_bin_left_pos_list.append(float(den_bin_edges[i]))

        for left_bin_edge in h_bin_left_pos_list:
            right_bin_edge = left_bin_edge + self.final_bin_size
            self.clusters.append((left_bin_edge, right_bin_edge))

        return


# class Cluster is given a position-frequency dictionary and can do the inversion nucleotide calculations
# I designed it to just do all the damn calculations when it is called. You can used class methods to
# draw graphs and stuff after that is all done anyways.
class Cluster:

    def __init__(self, pos_freq_dict, cbinsize=40, cperc=98,
                 clustersepmin=0, clustersepmax=10000, ntsepmin=0, ntsepmax=10000):

        self.pos_freq_dict = pos_freq_dict                  # position:frequency dictionary of this cluster
        self.pos_array = np.array(list(pos_freq_dict))      # unique position array of this cluster

        self.pos_min = self.pos_array.min()                 # minimum position
        self.pos_max = self.pos_array.max()                 # maximum position
        self.data_sum = 0                                   # sum of read counts in cluster

        # minimum and maximum frequency values
        self.freq_min = np.array(list(pos_freq_dict.values())).min()
        self.freq_max = np.array(list(pos_freq_dict.values())).max()
        self.pos_freq_max = list(pos_freq_dict.keys())[list(pos_freq_dict.values()).index(self.freq_max)]

        self.bin_size = cbinsize                            # how many nucleotides each bin should span
        self.c_sep_min = clustersepmin                      # limit to how close cluster pairs can be
        self.c_sep_max = clustersepmax                      # limit to how far cluster pairs can be
        self.n_sep_min = ntsepmin                           # limit to how close nucleotide pairs can be
        self.n_sep_max = ntsepmax                           # limit to how far nt pairs can be in the end
        self.count_percentile_threshold = cperc             # initial thresholding of counts for bins

        self.bin_size_tol = 5                               # nt size tolerance of cluster binning
        self.bins = 0                                       # what we eventually settled on for bins
        self.final_cbin_size = 0                            # what size we eventually got for the bins

        # NOTE! all clusters only has edges right now...
        self.all_cluster_bin_pairs = list()                 # list of all unique cluster bin pairs
        self.filtered_cluster_bin_pairs = list()            # lift of completely filtered bin pairs
        self.best_bin_pair = [(-1, -1), (-1, -1)]           # bin spans that best fulfill the conditions
        self.best_nt_pair = [(-1, 0), (-1, 0)]              # best scoring nucleotide pair with counts
        self.best_nt_pair_dist = 0                          # number of nt apart the pair is
        self.best_nt_pair_sum = 0                           # sum of scores of the best nt pair

        self.graph_nt_stream = 1000                         # amount of nt upstream and downstream when drawing

        self.filtered_cluster_bin_dictionary = dict()       # filtered cluster bin:frequency dict for indexing

        self.is_single_signal = 0                           # if one, may be a useless cluster (no buddy)
        self.per_dif_threshold = 95                         # if a single nt in a pair has 95% of the signal...
        self.signal = 0                                     # if a single hit, this is the nt causing the probs

        # =code to run on initialization=

        # make frequency histogram of data
        self.freq_histogram_counts, self.freq_histogram_edges, self.cluster_bin_dictionary = \
            self.make_freq_histogram()

        # filter out clusters by visualization

        # filter the dictionary by the cperc
        self.filter_by_read_count(self.count_percentile_threshold)

        # generate cluster bin pairs
        self.generate_cluster_bin_pairs()

        # filter the cluster bin pairs
        self.filter_bin_pairs()

        # if we are out of bin pairs, don't bother...
        if len(self.filtered_cluster_bin_pairs) == 0:
            self.is_single_signal = 1
            self.signal = (self.pos_freq_max, self.freq_max)
            pass

        else:

            # find the best cluster bin pair
            self.find_max_pair()

            # find the best nucleotides in there
            self.best_nt_pair[0] = self.find_best_nucleotide(self.sub_array(
                self.best_bin_pair[0][0], self.best_bin_pair[0][1]))
            self.best_nt_pair[1] = self.find_best_nucleotide(self.sub_array(
                self.best_bin_pair[1][0], self.best_bin_pair[1][1]))

            self.best_nt_pair_sum = self.best_nt_pair[0][1] + self.best_nt_pair[1][1]
            self.best_nt_pair_dist = abs(self.best_nt_pair[0][0] - self.best_nt_pair[1][0])

            # take a glance at the pair to see if we have a true inversion
            self.assess_nt_pair()

    # returns a frequency np histogram of a pos_freq_dict...so (10 10 10 10 20 20 30 30 30 30...etc.)
    def make_freq_histogram(self):

        data = list()
        for pos in self.pos_array:
            for i in range(0, self.pos_freq_dict[pos]):
                data.append(pos)
        data = np.array(data)
        bins = int((self.pos_max - self.pos_min) / self.bin_size)

        counts, edges = np.histogram(data, bins=bins)

        # check to make sure the bin size is right...sometimes, based on the pos array, it gets a bit small
        bin_size = edges[1] - edges[0]

        while abs(self.bin_size - bin_size) > self.bin_size_tol:

            if self.bin_size > bin_size:
                bins -= 1
            else:
                bins += 1

            counts, edges = np.histogram(data, bins=bins)
            bin_size = edges[1] - edges[0]

        # now that everything should be good, return the data array and the array of edges along with a
        # dictionary tying the two

        cbin_dict = dict()
        for i in range(0, len(counts)):
            cbin_dict[edges[i]] = counts[i]

        self.final_cbin_size = bin_size
        self.bins = bins

        for count in counts:
            self.data_sum += count

        return counts, edges, cbin_dict

    # returns an array subset based on data bounds
    def sub_array(self, pos_start, pos_end):

        a = list()
        for element in self.pos_array:
            if (element >= pos_start) and (element <= pos_end):
                a.append(element)
        return np.array(a)

    # filters the cluster bin dictionary based on cperc
    def filter_by_read_count(self, cperc):

        cluster_bin_cutoff_perc_val = np.percentile(self.freq_histogram_counts, cperc)

        for pos in self.cluster_bin_dictionary:
            read_count = self.cluster_bin_dictionary[pos]
            if read_count >= cluster_bin_cutoff_perc_val:
                self.filtered_cluster_bin_dictionary[pos] = self.cluster_bin_dictionary[pos]
        return

    # makes unique cluster bin pairs based on filtered cluster dictionary
    def generate_cluster_bin_pairs(self):

        pass_pos_list = list()
        for pos in self.filtered_cluster_bin_dictionary:
            pass_pos_list.append(pos)

        for i in range(0, len(pass_pos_list) - 1):
            this_bin_pos = pass_pos_list[i]
            rest_bin_pos = pass_pos_list[i + 1:]
            for pos in rest_bin_pos:
                self.all_cluster_bin_pairs.append((this_bin_pos, pos))
        return

    # filters bin pairs by class parameters
    def filter_bin_pairs(self):

        self.filtered_cluster_bin_pairs = self.all_cluster_bin_pairs
        for bin_pair in self.all_cluster_bin_pairs:
            bin1 = bin_pair[0]
            bin2 = bin_pair[1]
            if abs(bin2 - bin1) >= self.c_sep_max or abs(bin2 - bin1) <= self.c_sep_min:
                self.filtered_cluster_bin_pairs.remove(bin_pair)
        return

    # finds the maximally scoring pair of cluster bins
    def find_max_pair(self):
        bin_count_max, bin_max_pair = 0, (-1, -1)
        for bin_pair in self.filtered_cluster_bin_pairs:
            bin1 = bin_pair[0]
            bin2 = bin_pair[1]
            read_count = self.cluster_bin_dictionary[bin1] + self.cluster_bin_dictionary[bin2]

            if read_count > bin_count_max:
                bin_count_max = read_count
                bin_max_pair = bin_pair

        bin1_lb, bin1_ub = bin_max_pair[0], bin_max_pair[0]+ self.final_cbin_size
        bin2_lb, bin2_ub = bin_max_pair[1], bin_max_pair[1] + self.final_cbin_size
        self.best_bin_pair = ((bin1_lb, bin1_ub), (bin2_lb, bin2_ub))

        return

    # finds the best nucleotides in this cluster region
    def find_best_nucleotide(self, arr):

        best_nt = -1
        read_max = 0

        for pos in arr:
            if self.pos_freq_dict[pos] > read_max:
                best_nt = pos
                read_max = self.pos_freq_dict[pos]

        return best_nt, read_max

    # returns a frequency array over a position array using pos_freq_dict
    def make_freq_array(self, pos_array):
        data = list()
        for pos in pos_array:
            for i in range(0, self.pos_freq_dict[pos]):
                data.append(pos)
        return np.array(data)

    # looks at the nt pair and gives and idea of the legitness of the cluster based on class parameters.
    def assess_nt_pair(self):

        pos1 = self.best_nt_pair[0][0]
        pos2 = self.best_nt_pair[1][0]

        score1 = self.best_nt_pair[0][1]
        score2 = self.best_nt_pair[1][1]

        pos_dif = abs(pos1 - pos2)
        per_dif = 100 * (abs(score1 - score2) / (score1 + score2))

        # now, is one of the nucleotides scoring nearly 100% of the data?
        if per_dif > self.per_dif_threshold:
            self.is_single_signal = 1

        # or, are the nts waaay too close or far somehow despite our cluster distance thresholding?
        if (pos_dif >= self.n_sep_max) or (pos_dif <= self.n_sep_min):
            self.is_single_signal = 1

        # if the signal is up, find the offending nucleotide
        if self.is_single_signal == 1:

            if score2 > score1:
                self.signal = (pos2, score2)
            else:
                self.signal = (pos1, score1)
        return

    # uses matplotlib and sns to draw and save an illustration of the histogram data of the suggested inversion cluster
    def draw_inversion_site(self, save_path, show_fig='n'):

        # generate histogram data over this pos_array
        # first make a sub_array a little upstream and downstream of our positions
        pos1 = self.best_nt_pair[0][0]
        pos2 = self.best_nt_pair[1][0]

        start = pos1 - self.graph_nt_stream
        end = pos2 + self.graph_nt_stream
        arr = self.sub_array(start, end)

        # now make the frequency histogram
        dx = self.make_freq_array(arr)

        # use seabourn to draw our initial density histogram with a gaussian fit
        sns.distplot(dx, bins=30)

        # write vertical lines where we suspect the inversion pair to be
        plt.axvline(pos1, color='r')
        plt.axvline(pos2, color='r')

        # label our axes
        plt.xlabel('Nucleotide position')
        plt.ylabel('Read density')

        # draw arrows to annotate the vertical lines so you can see the exact position
        c_start = 'Cluster start: ' + str(pos1)
        c_end = 'Cluster end: ' + str(pos2)
        ymin, ymax = plt.ylim()
        plt.annotate(c_start, xy=(pos1, ymax), xycoords='data', xytext=(0.15, 0.95),
                     textcoords='figure fraction', arrowprops=dict(facecolor='black', shrink=0.05))

        plt.annotate(c_end, xy=(pos2, ymax), xycoords='data', xytext=(0.75, 0.95),
                     textcoords='figure fraction', arrowprops=dict(facecolor='black', shrink=0.05))

        # save our figure before showing
        plt.savefig(save_path)

        # show our figure if desired
        if show_fig == 'y':
            plt.show()

        # clear the figure
        plt.close()

        return


# append_to_csv takes a data tuple and appends to some csv file.
def append_to_csv(data, output):

    with open(output, 'a') as o:
        writer = csv.writer(o)
        writer.writerow(data)
    return

