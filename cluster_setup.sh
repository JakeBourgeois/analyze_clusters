#! usr/bin/bash

echo "Install script for clustering tools:"
echo "These tools require xcode command line tools, python3, and some dependencies."
echo "Installing...."

# install xcode tools
xcode-select --install

# install homebrew package manager
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"

# install python3
brew install python3

# install dependencies
pip3 install seaborn
pip3 install BioPython
pip3 install matplotlib
pip3 install numpy
pip3 install reportlab
pip3 install requests

echo "Setup completed! Now, once the SOR data is in the SOR folder, type:"
echo "$ python3 detect_inversion_clusters.py"
