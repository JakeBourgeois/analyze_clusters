from detect_inversions import *
import os

# define paths for file saving and loading
desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
working_path = os.path.join(desktop_path, "Cluster Detection")
results_path = os.path.join(working_path, 'Cluster Data')
sor_path = os.path.join(working_path, 'SOR Data')

result_file = os.path.join(working_path, 'cluster_pairs.csv')
freq_file = os.path.join(working_path, 'cluster_pos_freqs.csv')
accession_file = os.path.join(working_path, 'accession_list.txt')



# variables, make it so it loads by config file!
nbin_size = 5000  # number of nts per bin in the initial SOR read frequency histogram
cbin_size = 40  # number of nts I want each cluster bin to be
cluster_max_sep = 1000  # furthest away two bins can be when considering inversion sites
cluster_min_sep = 100  # closest two bins can be when considering inversion sites
cluster_bin_cutoff_perc = 98  # eliminate data below this percentile in the cluster analysis


accession_list = list()  # list of accession numbers to be processed
# load accession names from list
with open(accession_file, 'r') as a:
    try:
        for line in a:
            accession_list.append(line.split('\n')[0])
    except IOError as e:
        print("Error: Improper accession file.")
        quit(e.errno)

# if the accession list is empty, let the user know!
if len(accession_list) == 0:
    print("No accession numbers detected!")

# for each accession number...
for acc_num in accession_list:

    # define the input and output files
    sor_file = os.path.join(sor_path, acc_num+'.csv')

    acc_results_path = os.path.join(results_path, acc_num)
    if not os.path.exists(acc_results_path):
        os.makedirs(acc_results_path)

    graph_path = os.path.join(acc_results_path, "Cluster Graphs")
    if not os.path.exists(graph_path):
        os.makedirs(graph_path)

    analysis_file = os.path.join(acc_results_path, acc_num+' cluster analysis.csv')
    # if we have an analysis file here, delete it
    if os.path.exists(analysis_file):
        os.remove(analysis_file)

    # load SOR class, possible to make it so it loops through SOR files
    SOR_bug = SOR(acc_num, sor_file, binsize=nbin_size)

    # load up SOR thresholding
    SOR_bug.make_interactive_graphical_threshold(save_path=os.path.join(acc_results_path, acc_num + '_histogram'))

    print("Read density cutoff:", SOR_bug.read_cutoff)

    # for each cluster in the bug
    for cluster in SOR_bug.clusters:

        data_subset = SOR_bug.subset(cluster[0], cluster[1])

        # create a Cluster analysis class
        c = Cluster(data_subset, clustersepmin=100, ntsepmin=100)

        # if the signal is junk, print out some statement for now
        if c.is_single_signal == 1:
            print("Solitary signal found at:", c.signal)
            SOR_bug.spikes.append((c.signal[0], c.signal[0]))
            SOR_bug.signals.append((c.signal[0], c.signal[0]))
            SOR_bug.clust_dist.append(0)
            SOR_bug.clust_sum.append(c.signal[1])
            SOR_bug.clust_proportion_perc.append('{:.4}'.format(100*(c.signal[1] / c.data_sum)))
            SOR_bug.total_clust_prop_perc.append('{:.4}'.format(100*(c.signal[1] / SOR_bug.data_sum)))
            SOR_bug.cluster_map_int.append('N')
            SOR_bug.final_cbin_sizes.append(c.final_cbin_size)

        # otherwise, draw it here
        else:
            print("Cluster pair found at:", c.best_nt_pair[0], c.best_nt_pair[1])
            c_figname = os.path.join(graph_path, acc_num + '_cluster_' + str(c.best_nt_pair[0][0]))
            c.draw_inversion_site(c_figname, show_fig='n')
            SOR_bug.true_clusters.append(c.best_nt_pair)
            SOR_bug.signals.append((c.best_nt_pair[0][0], c.best_nt_pair[1][0]))
            SOR_bug.clust_dist.append(c.best_nt_pair_dist)
            SOR_bug.clust_sum.append(c.best_nt_pair_sum)
            SOR_bug.clust_proportion_perc.append('{:.4}'.format(100 * (c.best_nt_pair_sum / c.data_sum)))
            SOR_bug.total_clust_prop_perc.append('{:.4}'.format(100 * (c.best_nt_pair_sum/ SOR_bug.data_sum)))
            SOR_bug.cluster_map_int.append('Y')
            SOR_bug.final_cbin_sizes.append(c.final_cbin_size)

    # now lets dump the relevant data to disk

    # cluster stats and data
    num_signals = len(SOR_bug.clusters)                 # number of total signals detected
    all_signals = SOR_bug.signals                       # list of all signal positions
    local_props = SOR_bug.clust_proportion_perc         # list of percent read counts of inv. to cluster reads
    global_props = SOR_bug.total_clust_prop_perc        # list of percent read counts of inv. to all reads
    true_clusters = SOR_bug.true_clusters               # list of true clusters
    num_true_clusters = len(SOR_bug.true_clusters)      # number of true cluster pairs detected
    spikes = SOR_bug.spikes                             # list of spikes
    num_spikes = len(SOR_bug.spikes)                    # number of spikes detected
    is_inv_true = SOR_bug.cluster_map_int               # a list to tell us if we had a pair
    c_dist = SOR_bug.clust_dist                         # list of inversion pair nt span
    c_reads = SOR_bug.clust_sum                         # list of cluster read sums

    # run parameters not already defined at the top
    read_cutoff = SOR_bug.read_cutoff                   # read density cutoff
    final_ibin = SOR_bug.final_bin_size                 # the ultimate nt size of the bins
    final_cbins = SOR_bug.final_cbin_sizes              # final nt sizes of cluster bins list

    # let's write the cluster stats and data to a results file

    append_to_csv(['Accession Number:', acc_num], analysis_file)

    labels = ['Number of signals detected', 'Number of inversion pairs detected', 'Number of signal peaks detected']
    data = [num_signals, num_true_clusters, num_spikes]
    for i in range(0, len(labels)):
        d = (labels[i], data[i])
        append_to_csv(d, analysis_file)
    append_to_csv([''], analysis_file)

    header = ['Signal Start', 'Signal End', 'True Pair?', 'Inversion Length', 'Combined Read Count',
              'Percent Read to Cluster', 'Percent Read to All SORs']
    append_to_csv(header, analysis_file)
    for i in range(0, num_signals):
        data = [all_signals[i][0], all_signals[i][1], is_inv_true[i], c_dist[i], c_reads[i], local_props[i],
                global_props[i]]
        append_to_csv(data, analysis_file)
    append_to_csv([''], analysis_file)

    # now write run parameters
    append_to_csv(['RUN PARAMETERS:'], analysis_file)
    labels = ['Initial Density Cutoff', 'nt bin target for initial screen', 'nt bin achieved',
              'nt cluster bin target', 'nt cluster bins achieved', 'Minimum cluster bin distance',
              'Maximum cluster bin distance', 'Cluster bin count percentile cutoff', 'Minimum inversion size',
              'Maximum inversion size']
    data = [read_cutoff, nbin_size, final_ibin, 40, final_cbins, 100, 1000, 98, 100, 5000]
    for i in range(0, len(labels)):
        d = (labels[i], data[i])
        append_to_csv(d, analysis_file)












