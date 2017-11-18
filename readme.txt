README for Jake's Clustering Tools.

1. Download git repo from https://github.com/JakeBourgeois/analyze_clusters to your Desktop

You can also type from the terminal:
$ git clone https://github.com/JakeBourgeois/analyze_clusters.git

2. Rename "analyze_clusters" folder to "Cluster Detection" (I know, stupid, but you gotta
do it)

3. Open terminal and navigate to "Cluster Detection"
ex. $ cd Desktop\"Cluster Detection"

4. Run the setup script
$ bash cluster_setup.sh

5. Place any SOR data into the "SOR Data" folder as its accession number.csv.
For example, for C.diff, name the SOR file "FN545816.csv"

6. Make sure the headers in the SOR file are there. The first row should read:
"colors,POS,CIGAR,TLEN" or at least have POS!

7. Run the script
$ python3 detect_inversion_clusters.py

8. That should be it!