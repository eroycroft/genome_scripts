#!/bin/bash

####################################################
#                                                  #
# Useful dependencies for biologists using Linux   #
# EMILY ROYCROFT                                   #
#                                                  #
####################################################


# convert text file line endings between CRLF and LF

sudo apt install dos2unix 

# command line trashcan utility 

sudo apt install trash-cli

# Python libraries essential for bioinformatics (implemented in Python 2)

sudo apt install python-biopython 

sudo apt-get install python-numpy 

# install Bio::SeqIO

cpan App::cpanminus
sudo cpanm Bio::SeqIO

#Alternatively;
#sudo cpanm --force Bio::SeqIO



####################################################
#                                                  #
# To look for new packages available through the   #
# ADVANCED PACKAGE TOOL (APT) use the code         #
# sudo apt-cache search <search term>              #
#                                                  #
####################################################

