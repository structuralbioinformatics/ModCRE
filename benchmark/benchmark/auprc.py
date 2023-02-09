import os,sys,re
import numpy as np

def interpolate(TPA, TPB, FPA, FPB, total_positives):

        TPA = float(TPA)
        TPB = float(TPB)

        FPA = float(FPA)
        FPB = float(FPB)

	if (TPA-TPB)==0:
		return []
	#	return [ (TPA/total_positives, TPA/(TPA+FPA)), (TPA/total_positives, TPA/(TPA+FPA))]

        skew = (FPB-FPA)/(TPB-TPA)

        points = []

        for x in xrange(int(TPB)-int(TPA)+1):
		if (TPA+x+FPA+skew*x)>0:
			points.append( ((TPA+x)/total_positives, (TPA+x)/(TPA+x+FPA+skew*x)) )

        return points


def get_aupr(scores, labels):


        paired_list = zip(scores,labels)
        paired_list.sort(key=lambda x: x[0], reverse=True)

        total_positives = len([x[1] for x in paired_list if x[1] is True])

        TP_dict = {}

        TPA = 0
        TPB = 0
        FPA = 0
        FPB = 0

        #points = [(0,1)]
        points = []
        for cutoff, label in paired_list:
                TP_dict.setdefault(cutoff,[0,0])
                if label is True:
                        TP_dict[cutoff][0]+=1
                else:
                        TP_dict[cutoff][1]+=1

        sorted_cutoffs = sorted(TP_dict.keys(),reverse=True)

	TPB = TP_dict[sorted_cutoffs[0]][0]
	FPB = TP_dict[sorted_cutoffs[0]][1]

	#FIRST POINT
	points.extend(interpolate(0,TPB,0,FPB, total_positives))

        for xcutoff in xrange(1,len(sorted_cutoffs)):
                TPA += TP_dict[sorted_cutoffs[xcutoff-1]][0]
                TPB = TPA + TP_dict[sorted_cutoffs[xcutoff]][0]
                FPA += TP_dict[sorted_cutoffs[xcutoff-1]][1]
                FPB = FPA+TP_dict[sorted_cutoffs[xcutoff]][1]
                p = interpolate(TPA, TPB, FPA, FPB, total_positives)
		points.extend(p)

        x,y = zip(*points)

        return np.trapz(x=x,y=y)

