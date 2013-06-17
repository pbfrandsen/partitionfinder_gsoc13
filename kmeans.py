import time
import sys
import re
import os
# import util
import subprocess
import shlex

from math import log
import numpy as np
from sklearn.preprocessing import scale
from sklearn.cluster import KMeans


_binary_name = 'phyml'
if sys.platform == 'win32':
    _binary_name += '.exe'

# class PhymlError(PhylogenyProgramError):
#     pass

def find_program():
    """Locate the binary ..."""
    # pth = util.program_path
    # pth = os.path.join(pth, _binary_name)
    # This is a temporary solution to call phyml without messing other PF files
    pth = _binary_name

    # log.debug("Checking for program %s in path %s", _binary_name, pth)
    if not os.path.exists(pth) or not os.path.isfile(pth):
        log.error("No such file: '%s'", pth)
        raise PhymlError
    # log.debug("Found program %s at '%s'", _binary_name, pth)
    return pth

_phyml_binary = None


def run_phyml(command):
    global _phyml_binary
    if _phyml_binary is None:
        _phyml_binary = find_program()

    #turn off any memory checking in PhyML - thanks Jess Thomas for pointing out this problem
    command = "%s --no_memory_check" % (command)
    print command

    # Add in the command file
    # log.debug("Running 'phyml %s'", command)
    command = "\"%s\" %s" % (_phyml_binary, command)
    print command

    # Note: We use shlex.split as it does a proper job of handling command
    # lines that are complex
    p = subprocess.Popen(
        shlex.split(command),
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)

    # Capture the output, we might put it into the errors
    stdout, stderr = p.communicate()
    print stdout
    print stderr
    # p.terminate()

    if p.returncode != 0:
        # log.error("Phyml did not execute successfully")
        # log.error("Phyml output follows, in case it's helpful for finding the problem")
        # log.error("%s", stdout)
        # log.error("%s", stderr)
        raise PhymlError

def rate_cat_likelhood_parser(phyml_lk_file):
	'''Takes the path of the phyml_lk file and returns
	a dictionary with sites as keys and a list of
	likelihoods under different rate categories as
	as the values
	'''
	# Open the file and assign it to a variable
	try:
		phyml_lk_file = open(phyml_lk_file, "r")
	except IOError:
		raise IOError("Could not find site likelihood file!")
	# Start the count of the line number
	count = 0
	# Begin site number count
	site = 1
	# Create dictionary for output
	all_rates_dict = {}
	# Read the file line by line
	for line in phyml_lk_file.readlines():
		count += 1
		if count > 7 and count < 1000007:
			# Split the data and objects into a list and remove 
			# unnecessary whitespace
			rate_list = " ".join(line.split()).split(" ")
			# This could also be accomplished with regex
			# rate_list = re.sub(' +', ' ', line).split(" ")
			# Raise error if rate variation is not modeled
			if len(rate_list) < 4:
				raise IOError("Rate variation was not modeled in "
					"PhyML analysis. Please choose > 1 rate category")
			# Remove first two and last numbers from list
			rate_list.pop(0)
			rate_list.pop(0)
			rate_list.pop(-1)
			# Convert all likelihoods to floats
			for i in range(len(rate_list)):
				rate_list[i] = float(rate_list[i])
			# Add site number and rate list to dictionary
			all_rates_dict[site] = rate_list
			site += 1
		# If there are more than 1,000,000 nts in a phyml_lk.txt 
		# file, the numbers merge with the site likelihoods
		elif count >= 1000007:
			# Split the data
			rate_list = " ".join(line.split()).split(" ")
			# Remove last and first number from list
			rate_list.pop(0)
			rate_list.pop(-1)
			# Convert all likelihoods to floats
			for i in range(len(rate_list)):
				rate_list[i] = float(rate_list[i])
			# Add site number and rate list to dictionary
			all_rates_dict[site] = rate_list
			site += 1
	# print all_rates_dict
	phyml_lk_file.close()
	return all_rates_dict

if __name__ == "__main__":
	phylip_filename = sys.argv[1]
	run_phyml("-i " + str(phylip_filename) + " -m GTR -o r --print_site_lnl -c 1")
	phyml_lk_file = str(phylip_filename) + "_phyml_lk.txt"
	likelihood_dictionary = rate_cat_likelhood_parser(phyml_lk_file)
	kmeans(likelihood_dictionary)
