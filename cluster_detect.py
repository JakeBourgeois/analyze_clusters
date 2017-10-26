from cluster_tools import *

desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
cluster_path = os.path.join(desktop_path, "DetectCluster")

sor_file = os.path.join(cluster_path, 'cdiff.csv')
result_file = os.path.join(cluster_path, 'cluster_pairs.csv')
freq_file = os.path.join(cluster_path, 'cluster_pos_freqs.csv')

#detect_clusters(sor_file, freq_file, result_file, fig_path=cluster_path, ctol=2000)

detect_inversion_clusters_interactive(sor_file, nbin_size=20000)
