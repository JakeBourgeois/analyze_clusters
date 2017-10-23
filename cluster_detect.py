from cluster_tools import *

desktop_path = os.path.join(os.path.expanduser('~'), 'Desktop')
cluster_path = os.path.join(desktop_path, "DetectCluster")

sor_file = os.path.join(cluster_path, 'data.csv')
result_file = os.path.join(cluster_path, 'checkme.csv')
freq_file = os.path.join(cluster_path, 'freqs.csv')

detect_clusters(sor_file, freq_file, result_file, ctol=100000)

