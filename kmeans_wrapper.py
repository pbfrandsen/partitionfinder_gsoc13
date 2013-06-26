import sys
import time
import re
import os
# import util
import subprocess
import shlex
import csv

from kmeans import phyml_likelihood_parser
from kmeans import run_raxml
from kmeans import raxml_likelihood_parser
from kmeans import kmeans

from math import log
import numpy as np
from sklearn.preprocessing import scale
from sklearn.cluster import KMeans
from collections import defaultdict


def kmeans_wrapper(dictionary, max_ks = 25):
    '''This function performs kmeans on
    a specified number of different k's on a dictionary parsed
    from a *_phyml_lk.txt file. It then calculates
    the within sum of squares (wss) of the k's and returns
    the optimal number of k's and a dictionary of
    the k's and sites belonging to the k's as output
    determined by the first time the wss shows less than
    is less than 1/10 with the previous number of k's
    '''
    # Set variables for best num_ks and best dictionary
    # Set variable for most recent wss and new wss
    # Perform kmeans starting with a k of 2
    # When characters are returned, create list of sites and
    # calculate wss on each cluster
    pass

def ftest(likelihood_dict):
    '''This function takes a site likelihood dictionary for different
    clusters as input and returns the ratio of between-group variability
    over within-group variability
    '''
    pass

def multi_ftest(likelihood_dict):
    '''Takes a multi-dimensional dictionary of site likelihoods under
    different rate categories as input and returns the ratio of
    multivariate between-group variability over multivariate within-
    group variability using Wilks lambda
    '''
    pass

def ss(list_of_likelihoods):
    '''Input a list of likelihoods, and returns the sum of squares
    of the list
    '''
    sums_of_squares = 0
    # Need to log transform the data as you do during kmeans
    log_list_of_likelihoods = []
    for i in list_of_likelihoods:
        log_list_of_likelihoods.append(log(float(i)))
    # Calculate mean
    mean_likelihood = sum(log_list_of_likelihoods)/len(log_list_of_likelihoods)
    # Calculate and return sum of squares
    for i in log_list_of_likelihoods:
        sums_of_squares += (i - mean_likelihood)**2
    return sums_of_squares

def wss(likelihood_lists):
    '''Inputs a list that contains lists of likelihoods from different
    clusters, outputs the within sum of squares
    '''
    within_sum_of_squares = 0
    # Call ss on each list to get global wss
    for i in likelihood_lists:
        within_sum_of_squares += ss(i)
    return within_sum_of_squares

def make_likelihood_list(likelihood_dict, site_categories):
    '''Takes a scaled_array and site categories
    as input and returns a list of tuples with centroids 
    '''
    new_like_dict = {k: [likelihood_dict[i][0] for i in v] for k, v in site_categories.items()}
    rate_list = []
    for i in new_like_dict:
        rate_list.append(new_like_dict[i])
    return rate_list

if __name__ == "__main__":
    # list_of_likelihoods = [[1.421312, 5.12321, 8.3432, 3.4122], [3.43234, 4.1232, 1.532432, 2.53243]]
    # print wss(list_of_likelihoods)
    phylip_filename = sys.argv[1]
    outfile_command = "testing"
    run_raxml("-s " + str(phylip_filename) + " -m GTRGAMMA -n " + outfile_command + " -y -p 23456")
    outfile_command2 = "testing2"
    start = time.clock()
    run_raxml("-s " + str(phylip_filename) + " -m GTRGAMMA -n " + outfile_command2 + " -f g -p 23456 -z RAxML_parsimonyTree." + outfile_command)
    stop = time.clock()
    print "RAxML took " + str(stop-start) + " seconds!"
    likelihood_dict = raxml_likelihood_parser("RAxML_perSiteLLs." + outfile_command2)
    site_categories = kmeans(likelihood_dict, number_of_ks = 4)[1]
    likelihood_lists = make_likelihood_list(likelihood_dict, site_categories)
    print wss(likelihood_lists)

