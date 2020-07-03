ProRank
Version 1.0
=============

Cite: Nazar Zaki, Jose Berengueres, Dmitry Efimov (2012). "Detection Of Protein Complexes Using A Protein Ranking Algorithm". Proteins, 80(10): 2459-68, 2012. Wiley


Programs/packages you must have:
1.	Python (available at http://python.org/)
2.	Numpy (available at http://pypi.python.org/pypi/numpy)

From the commandline write the following commands:

$ python ProRank-0.1.py –i [Int_file] –s [Sim_file]

Int_file = the files contains the protein interactions
Sim_file = the files contains the similarity matrix of the proteins.


Other programs you will need to use to prepare the datasets: 
1.	Fasta34.exe or higher (Program and instructions to run it are available at http://faculty.virginia.edu/wrpearson/fasta/fasta_versions.html). Can be used to calculate the similarity between proteins in the network
2.	filterNadd_ppi.exe (Program and instructions to run it are available at http://www.comp.nus.edu.sg/~wongls/projects/complexprediction/CMC-26may09/). filterNadd_ppi.exe assigns weights to protein pairs in the interaction file.



SimMat
=============

$ perl SimMat.pl interaction.txt

SimMat.pl calculates the similarity between proteins in the network. The input to this program is the interactions.txt file. The program automatically retrieves the corresponding protein sequences from the file “all_seq.fasta” (available in the dataset folder). The output of this program is the similarity matrix between proteins in the network “similaritymatrix.txt”.

Programs you must have:
1.	Fasta34.exe or higher (Program and instructions to run it are available at http://faculty.virginia.edu/wrpearson/fasta/fasta_versions.html).