import os, sys, re
import ConfigParser
import optparse
import subprocess
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import string
import argparse
from scipy import stats
import cPickle
from collections import Counter

#-------------#
# Options     #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.

    """

    parser = optparse.OptionParser("python enrichment_pvalue.py  -m number_of_models -n number_of_motifs -e enrichment --top top_selected [-o output]",
                                   epilog      = '@Oliva\'s lab 2018',
                                   description = "The program calculates a table of p-values as a function of the number of top selected motifs and the achieved enrichment ")

    parser.add_option("-m", action="store", default=None, type="int", dest="number_of_models", help="Number of models ", metavar="{integer}")
    parser.add_option("-o", action="store", default="output_pvalues.tab", type="string", dest="output_file", help="Output table (default = 'output_pvalues.tab')", metavar="{filename}")
    parser.add_option("-n", action="store", default=None, type="int", dest="number_of_motifs", help="Number of Motifs")
    parser.add_option("--top", action="store", default=None, type="int", dest="top_selected", help="Number of top selected motifs")
    parser.add_option("-e", "--enrichment", default=None, type="float", action="store", dest="enrichment", help="Ratio of enrichment of a motif")
    parser.add_option("-v","--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)")
   
    (options, args) = parser.parse_args()
#    if options.number_of_models is None or options.number_of_motifs is None  or options.top_selected is None  or options.enrichment is None :
    if options.number_of_models is None or options.number_of_motifs is None   :
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def log_p_value(m,e,n,x):
    y=None
    if x<n and e<=1:
       #fx = np.log(float(x)/float(n)) + (x-1)*np.log(float(n-1)/float(n))
       #fz = x * np.log(float(n-1)/float(n))
       #fx = log_combinatorial(x,1) - np.log(float(n)) + (x-1)* (np.log(float(n-1)) - np.log(float(n)))
       #fz = log_combinatorial(x,x) + (x)* (np.log(float(n-1)) - np.log(float(n)))
       fx = np.log(float(x)/float(n))
       fz = np.log(float(n-x)/float(n))
       y  = log_combinatorial(m,m*e)  +  m * ( e*fx  + (1-e)*fz )
    return y

def log_factorial(n):
    x=0
    if n>0:
       for i in xrange(1,int(n+1)): 
           x= x + np.log(float(i)) 
    return x

def log_combinatorial(n,m):
    x=0
    if n>m:
       x =  log_factorial(n) - log_factorial(m) - log_factorial(n-m)
    return x

def writetable(m,n,step,top_selected,enrichment,output_file):
    fo=open(output_file,"w")
    fo.write("#%14s\t%15s\t%15s\n"%("top","enrichment","p-value"))
    for p in xrange(1,101,1):
        e= (float)(p)/100 
        for x in xrange(1,n,step):
            log_p  = log_p_value(m,e,n,x)
            if log_p is not None:
               if log_p>-1000 :
                  if top_selected is None and enrichment is None:
                     fo.write("%15.5f\t%15.5f\t%15.5e\n"%(x,p,np.exp(log_p)))
                  else:
                     if top_selected is None: 
                         if abs(enrichment-p) <=5: fo.write("%15.5f\t%15.5f\t%15.5e\n"%(x,p,np.exp(log_p)))
                     if enrichment is None: 
                         if abs(top_selected-x)<=step: fo.write("%15.5f\t%15.5f\t%15.5e\n"%(x,p,np.exp(log_p)))
                     if enrichment is not None and top_selected is not None:
                         try:
                           if abs(enrichment-p) <=5 and abs(top_selected-x)<=step: fo.write("%15.5f\t%15.5f\t%15.5e\n"%(x,p,np.exp(log_p)))
                         except Exception as e:
                           sys.stderr.write("%s\n"%e)
        x=n-1
        log_p  = log_p_value(m,e,n,x)
        if log_p is not None:
                  if top_selected is None and enrichment is None:
                     fo.write("%15.5f\t%15.5f\t%15.5e\n"%(x,p,np.exp(log_p)))
                  else:
                     if top_selected is None: 
                         if abs(enrichment-p) <=5: fo.write("%15.5f\t%15.5f\t%15.5e\n"%(x,p,np.exp(log_p)))
                     if enrichment is None: 
                         if abs(top_selected-x)<=step: fo.write("%15.5f\t%15.5f\t%15.5e\n"%(x,p,np.exp(log_p)))
                     if enrichment is not None and top_selected is not None:
                         if abs(enrichment-p) <=5 and abs(top_selected-x)<=step: fo.write("%15.5f\t%15.5f\t%15.5e\n"%(x,p,np.exp(log_p)))
    fo.close()

def best_enrichment(m,n,x):
    max_LogpValue= -np.log(float(n))
    y = 1
    z = 1
    max_pvalue=0
    skip=False
    for p in xrange(100,0,-1):
        e= (float)(p)/100 
        log_p  = log_p_value(m,e,n,x)
        if np.exp(log_p)>max_pvalue: 
           z=p
           max_pvalue=np.exp(log_p)
        if log_p >= max_LogpValue and not skip:
           y=p
           skip=True
    return (y,z)

#-------------#
# Main        #
#-------------#

def main():

    # Arguments & Options #
    options = parse_options()
    m= options.number_of_models
    n= options.number_of_motifs
    top_selected=options.top_selected
    enrichment  =options.enrichment
    step=(int)(float(n)/100)
    if top_selected is not None: 
       y,z=best_enrichment(m,n,top_selected)
       sys.stdout.write("Minimum enrichment required: %f (maximum in %f)\n"%(y,z))
    writetable(m,n,step,top_selected,enrichment,options.output_file)

                    
if __name__ == "__main__":
    main() 
   
