=============================================================================

 Package: ChaLearn Connectomics Challenge Sample Code
 Source: http://connectomics.chalearn.org
 Authors: Javier Orlandi, Mehreen Saeed, Isabelle Guyon
 Date: February 2014
 Version: 1.0
 Contact: causality@chalearn.org
 License: GPL v3 see http://www.gnu.org/licenses/

=============================================================================

ALL INFORMATION, SOFTWARE, DOCUMENTATION, AND DATA ARE PROVIDED "AS-IS". CHALEARN, KAGGLE AND/OR OTHER ORGANIZERS DISCLAIM ANY EXPRESSED OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR ANY PARTICULAR PURPOSE, AND THE WARRANTY OF NON-INFRIGEMENT OF ANY THIRD PARTY'S INTELLECTUAL PROPERTY RIGHTS. IN NO EVENT SHALL THE ORGANIZERS BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF SOFTWARE, DOCUMENTS, MATERIALS, PUBLICATIONS, OR INFORMATION MADE AVAILABLE FOR THE CHALLENGE. 

=============================================================================

The goal of the challenge is to predict whether there is a (directed) connection between neuron i and neuron j in a network of 1000 neurons (self-connections permitted). We provide one hour of recording of the activity of all neurons as time series and the position of the neurons (arranged on a flat surface). The data, which are simulated, reproduce as faithfully as possible neural activity measured with calcium fluorescence imaging of neural cultures.

SAMPLE DATA
===========
The starter kit includes sample data files from 100 neuron networks, that are substitute for the larger 1000 neuron networks used for training, validation and testing in the challenge. They have base names mocktrain, mockvalid, and mocktest. The files, prefixed with fluorescence_, networkPositions_, and network_:

1) FLUORESCENCE = fluorescence_[basename].txt: comma separated files including time series of neural activity, each row representing a sample and each column a neuron. Signal is sampled at 20ms intervals.
2) NETWORK POSITIONS = networkPositions_[basename].txt: comma separated files, each row representing a neuron. First column is the X position and second column the Y position. Neurons span a 1mm2 square area.
3) NETWORK = network_[basename].txt: comma separated files representing the network architecture. Each row is a connection. The column structure is of the form I,J,W denoting a connection from I to J with weight W. Connections with positive weight (usually 1) are present. Pairs that are absent or have a weight -1 are inexistent or blocked in the simulations (which is the same thing as far as this challenge is concerned).

Note: For the "real" validation and test networks, the truth values of the connections (network_valid.txt and network_test.txt) are NOT provided during the challenge.

SAMPLE NETWORK RECONSTRUCTION CODE
==================================
The main entry point of the starter kit is:
> challengeMain;

The script takes a few minutes to run on the sample network of 100 neurons but it is very slow and requires a lot of memory for 1000 neurons. For that reason, we also provide a script that runs in a few minutes on the 1000 neuron networks, but includes only the correlation reconstruction method (a pretty good baseline though):
> challengeFastBaseline;

Both scripts perform the following functions:
1) Loads the fluorescence file as a matrix F, neurons in columns; each line is a time sample.
2) Performs various steps to compute scores for neuron i -> neuron j [correlation only for "challengeFastBaseline" and a choice of methods for "challengeMain", including the GTE algorithm, from Stetter, O., Battaglia, D., Soriano, J. & Geisel, T. Model-free reconstruction of excitatory neuronal connectivity from calcium imaging signals. PLoS Comput Biol 8, e1002653 (2012). 
The resulting "scores" matrix is a matrix N x N, N being the number of neurons, each entry (i, j) indicating the "confidence" that neuron i -> neuron j.
3) Writes the scores in Kaggle submission format as a 2-column csv file
NET_neuronI_neuronJ Strength indicating the "confidence" that neuron i -> neuron j.
4) Computes AUC performance [if the network architecture is provided, i.e. for training networks only].

CODE TO VISUALIZE DATA
======================
We provide code to visualize the time series, make a movie of the neural activity, and browse through the data:

> challengeVizualization;

CODE TO TRAIN A PREDICTIVE MODEL
================================
To combine several scores used as "features" we provide sample code training a simple linear model: 

> challengeTrain;

