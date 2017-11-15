from detect_inversions import *

# define paths for file saving and loading
desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
working_path = os.path.join(desktop_path, "Cluster Data")
results_path = os.path.join(working_path, 'Results')
sor_path = os.path.join(working_path, 'SOR Data')

sor_file = os.path.join(working_path, 'dehalogenans.csv')
result_file = os.path.join(working_path, 'cluster_pairs.csv')
freq_file = os.path.join(working_path, 'cluster_pos_freqs.csv')



# variables, make it so it loads by config file!
nbin_size = 5000  # number of nts per bin in the initial SOR read frequency histogram
cbin_size = 40  # number of nts I want each cluster bin to be
cluster_max_sep = 1000  # furthest away two bins can be when considering inversion sites
cluster_min_sep = 100  # closest two bins can be when considering inversion sites
cluster_bin_cutoff_perc = 98  # eliminate data below this percentile in the cluster analysis

# load SOR class, possible to make it so it loops through SOR files
SOR_bug = SOR(sor_file, binsize=nbin_size)

# load up SOR thresholding
SOR_bug.make_interactive_graphical_threshold(save_path=os.path.join(cluster_path, 'cdiff'))

# for each cluster in the bug
for cluster in SOR_bug.clusters:

    data_subset = SOR_bug.subset(cluster[0], cluster[1])
    data_spikes = list()

    # create a Cluster analysis class
    c = Cluster(data_subset, clustersepmin=100, ntsepmin=100)

    # if the signal is junk, print out some statement for now
    if c.is_single_signal == 1:
        print("Solitary signal found at:", c.signal)
        data_spikes.append(c.signal[0])

    # otherwise, draw it here
    else:
        print("Cluster pair found at:", c.best_nt_pair[0], c.best_nt_pair[1])
        c_figname = os.path.join(cluster_path, 'cluster_' + str(c.best_nt_pair[0][0]))
        c.draw_inversion_site(c_figname)






