import numpy as np


myTpmFile = "./data_analyses/01_Expression_level_screening/all_counts.csv"

Cutoff = 1.0

filteredOut = open("%s_%s_expression_filter.txt" %
                (myTpmFile.split(".csv")[0], Cutoff), "w")
filteredLogedOut = open("%s_%s_expression_filter.log2.txt" % (
    myTpmFile.split(".csv")[0], Cutoff), "w")
Cutoff = float(Cutoff)
grp_idx = {}
i = 0
for line in open(myTpmFile):
    line = line.strip()
    pa = line.split(",")
    if i == 0:
        filteredOut.write('\t'.join(pa) + "\n")
        filteredLogedOut.write('\t'.join(pa) + "\n")

        for c in pa[1:]:
            if "_" in c:

                grp = '_'.join(c.split("_")[ : -1])
            else:
                grp = c[ :2]
                
            if grp not in grp_idx:
                grp_idx[grp] = []
            grp_idx[grp].append(pa.index(c))
            
    else:
        count_grp = {}
        filterPa = []
        filterLogPa = []
        filterPa.append(pa[0])
        filterLogPa.append(pa[0])
        for key in grp_idx.keys():
            count_grp[key] = 0
            for each_idx in grp_idx[key]:
                if float(pa[each_idx]) >= Cutoff:
                    count_grp[key] += 1
        

        count = 0
        for key in count_grp.keys():
            #
            if count_grp[key] / float(len(grp_idx[key])) >= 0.5:
                # if count_grp[key] >= 1:
                count += 1
        

        if count >= 1:
            for c in pa[1:]:
                filterPa.append(c)
                filterLogPa.append(str(np.log2(float(c) + 1)))
            print("\t".join(filterPa),file=filteredOut)
            print("\t".join(filterLogPa),file=filteredLogedOut)
        
    i += 1
filteredOut.close()
filteredLogedOut.close()