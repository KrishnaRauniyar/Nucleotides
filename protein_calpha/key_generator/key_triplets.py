import math
import pandas as pd
import requests
import sys
import csv
import multiprocessing
from joblib import Parallel, delayed

args= sys.argv[1:]
if len(args) < 3 :
    print("please specify all parameters")
    print("command : python key_triplets.py csv_file_path pdf_folder_path protein_path")
    quit(0)

CSV_FILE_PATH = args[0]
PDB_DIR_PATH = args[1]
PROTEIN_DIR_PATH = args[2]

class AminoAcidAnalyzer:
    def __init__(self, dtheta, dLen, numOfLabels):
        self.dtheta = dtheta
        self.dLen = dLen
        self.numOfLabels = numOfLabels
        # These are the different labels of the amino acids with their unique code (4 to 23)
        self.aminoAcidLabelWithCode = {
            "PHE": 4,
            "CYS": 5,
            "GLN": 6,
            "GLU": 7,
            "LEU": 8,
            "HIS": 9,
            "ILE": 10,
            "GLY": 11,
            "LYS": 12,
            "MET": 13,
            "ASP": 14,
            "PRO": 15,
            "SER": 16,
            "THR": 17,
            "TRP": 18,
            "TYR": 19,
            "VAL": 20,
            "ALA": 21,
            "ARG": 22,
            "ASN": 23
        }
        # Store animo acids code, x coordinate, y coordinate, z coordinate
        self.aminoAcidCode = {}
        self.xCoordinate = {}
        self.yCoordinate = {}
        self.zCoordinate = {}
        # This is to store the sequence number
        self.aminoSeqNum = {}
        # These three holds the label code of three amino acid, sorted them, and then store the sorted index 
        self.initAminoLabel = [0, 0, 0]
        self.sortedAminoLabel = [0, 0, 0]
        self.sortedAminoIndex = [0, 0, 0]

        # keys with its frequency
        self.keyFreq = {}
        # totak number of keys, and max distance list
        self.totalKeys = []
        self.maxDistList = []

    # This is reading the csv file and generating proteins list and separating chains and using below in the code function calcuTheteAndKey for filter in pdb file
        self.chain_dict = {}
        self.proteinList = []
    def readCSVProteinChain(self, csvFile):
        df = pd.read_csv(csvFile)
        self.chain_dict = dict([(p,c) for p,c in zip(df.protein, df.chain)])
        self.proteinList = df['protein'].tolist()
    
    # Code to download the data set from rcsb.org
    def downloadDataSet(self, download_dir = PDB_DIR_PATH):
        for file in self.proteinList:
            pdb_url = 'https://files.rcsb.org/download/'+str(file)+".pdb"  
            save_path = str(file)+'.pdb'
            print(save_path)
            response = requests.get(pdb_url)
            if response.status_code == 200:
                with open(download_dir+"/"+save_path, 'wb') as pdb_file:
                    pdb_file.write(response.content)
                print("PDB file for" ,file," downloaded and saved as ",save_path,".")
            else:
                print("Failed to download PDB file for", file ,". Status code:", response.status_code)

    def thetaClass(self, theta):
        if theta>=0 and theta<12.11:
            classT=1
        elif theta>=12.11 and theta<17.32:
            classT=2
        elif theta>=17.32 and theta<21.53:
            classT=3
        elif theta>=21.53 and theta<25.21:
            classT=4
        elif theta>=25.21 and theta<28.54:
            classT=5
        elif theta>=28.54 and theta<31.64:
            classT=6
        elif theta>=31.64 and theta<34.55:
            classT=7
        elif theta>=34.55 and theta<37.34:
            classT=8
        elif theta>=37.34 and theta<40.03:
            classT=9
        elif theta>=40.03 and theta<42.64:
            classT=10
        elif theta>=42.64 and theta<45.17:
            classT=11
        elif theta>=45.17 and theta<47.64:
            classT=12
        elif theta>=47.64 and theta<50.05:
            classT=13
        elif theta>=50.05 and theta<52.43:
            classT=14
        elif theta>=52.43 and theta<54.77:
            classT=15
        elif theta>=54.77 and theta<57.08:
            classT=16
        elif theta>=57.08 and theta<59.38:
            classT=17
        elif theta>=59.38 and theta<61.64:
            classT=18
        elif theta>=61.64 and theta<63.87:
            classT=19
        elif theta>=63.87 and theta<66.09:
            classT=20
        elif theta>=66.09 and theta<68.30:
            classT=21
        elif theta>=68.30 and theta<70.5:
            classT=22
        elif theta>=70.5 and theta<72.69:
            classT=23
        elif theta>=72.69 and theta<79.2:
            classT=24
        elif theta>=79.2 and theta<81.36:
            classT=25
        elif theta>=81.36 and theta<83.51:
            classT=26
        elif theta>=83.51 and theta<85.67:
            classT=27
        elif theta>=85.67 and theta<87.80:
            classT=28
        elif theta>=87.80 and theta<90.00:
            classT=29
        return classT

    def dist12Class(self, dist12):
        if (dist12<3.83):
            classL=1
        elif dist12>=3.83 and dist12<7.00:
            classL=2
        elif dist12>=7.00 and dist12<9.00:
            classL=3
        elif dist12>=9.00 and dist12<11.00:
            classL=4
        elif dist12>=11.00 and dist12<14.00:
            classL=5
        elif dist12>=14.00 and dist12<17.99:
            classL=6
        elif dist12>=17.99 and dist12<21.25:
            classL=7
        elif dist12>=21.25 and dist12<23.19:
            classL=8
        elif dist12>=23.19 and dist12<24.8:
            classL=9
        elif dist12>=24.8 and dist12<26.26:
            classL=10
        elif dist12>=26.26 and dist12<27.72:
            classL=11
        elif dist12>=27.72 and dist12<28.9:
            classL=12
        elif dist12>=28.9 and dist12<30.36:
            classL=13
        elif dist12>=30.36 and dist12<31.62:
            classL=14
        elif dist12>=31.62 and dist12<32.76:
            classL=15
        elif dist12>=32.76 and dist12<33.84:
            classL=16
        elif dist12>=33.84 and dist12<35.13:
            classL=17
        elif dist12>=35.13 and dist12<36.26:
            classL=18
        elif dist12>=36.26 and dist12<37.62:
            classL=19
        elif dist12>=37.62 and dist12<38.73:
            classL=20
        elif dist12>=38.73 and dist12<40.12:
            classL=21
        elif dist12>=40.12 and dist12<41.8:
            classL=22
        elif dist12>=41.8 and dist12<43.41:
            classL=23
        elif dist12>=43.41 and dist12<45.55:
            classL=24
        elif dist12>=45.55 and dist12<47.46:
            classL=25
        elif dist12>=47.46 and dist12<49.69:
            classL=26
        elif dist12>=49.69 and dist12<52.65:
            classL=27
        elif dist12>=52.65 and dist12<55.81:
            classL=28
        elif dist12>=55.81 and dist12<60.2:
            classL=29
        elif dist12>=60.2 and dist12<64.63:
            classL=30
        elif dist12>=64.63 and dist12<70.04:
            classL=31
        elif dist12>=70.04 and dist12<76.15:
            classL=32
        elif dist12>=76.15 and dist12<83.26:
            classL=33
        elif dist12>=83.26 and dist12<132.45:
            classL=34
        elif dist12>=132.45:
            classL=35
        return classL
    
    # Calculating the distance between the amino acids 
    def calDistance(self, l1_index, l2_index):

        x1 = self.xCoordinate[l1_index]
        x2 = self.xCoordinate[l2_index]
        y1 = self.yCoordinate[l1_index]
        y2 = self.yCoordinate[l2_index]
        z1 = self.zCoordinate[l1_index]
        z2 = self.zCoordinate[l2_index]
        return math.sqrt((x2 - x1) ** 2 + (y2 - y1) **2 + (z2 - z1) **2)

    def findTheIndex(self, l2_index, p1, q1, r1):
        if l2_index == p1:
            l1_index0 = q1
            l2_index1 = r1
        elif l2_index == q1:
            l1_index0 = p1
            l2_index1 = r1
        elif l2_index == r1:
            l1_index0 = p1
            l2_index1 = q1
        return l1_index0, l2_index1
    
    # This function is used to extract alpha carbons and then calculate theta and key
    def calcuTheteAndKey(self, fileName, chain, protein_path= PROTEIN_DIR_PATH):
        self.totalKeys = []
        self.maxDistList = []
        self.keyFreq = {}
        incrementVal=0
        self.newDictForLexicalCode = {}
        self.aminoAcidCode = {}
        self.aminoSeqNum = {}

        tripletsFile = open(protein_path+"/"+fileName+"_"+chain+".keys_theta29_dist35","w")
        keyFreqFile = open(protein_path+"/"+fileName+"_"+chain+".keys_Freq_theta29_dist35","w")
        # This is the extraction of CA atoms only from the pdb file
        incrementVal=0
        with open("./pdb_files/"+fileName+".pdb", "r") as pdbFile:
            for line in pdbFile:
                try:
                    # Do not take CA if the peptide bond is broken i.e after the TER
                    if ((line[0:6].rstrip()=="ENDMDL") or (line[0:6].rstrip()=='TER'and line[21].rstrip()==chain)):
                        break
                    if (line[0:6].rstrip()=="MODEL" and int(line[10:14].rstrip())>1):
                        break                       
                    if (line.startswith("ATOM") and (line[13:15].rstrip())=="CA" and (line[16]=='A'or line[16]==' ') and line[21:22].strip()==chain and line[17:20]!= "UNK"):
                        # Reading the lines in pdb file and then assigning residue (VAL) to its value (20)
                        self.aminoAcidCode[incrementVal]=int(self.aminoAcidLabelWithCode[line[17:20]])
                        # This is the sequence number of the amino acid stored (Residue seq number)
                        self.aminoSeqNum[incrementVal]=str(line[22:27])
                        self.xCoordinate[incrementVal]=(float(line[30:38]))
                        self.yCoordinate[incrementVal]=(float(line[38:46]))
                        self.zCoordinate[incrementVal]=(float(line[46:54]))
                        incrementVal+=1
                except Exception as e:
                    print("Their is an error in: ", line, pdbFile)
                    print(e)

        # This is the four rules that calculates the label, theta, and key (3 amino acids form a triplet)
        for i in range(0, len(self.aminoAcidCode) - 2):
            for j in range(i+1, len(self.aminoAcidCode) - 1):
                for k in range(j+1, len(self.aminoAcidCode)):
                    # This is a dictionary to keep the index and the labels
                    labelIndexToUse={}
                    # First, Second and Third label and Index
                    labelIndexToUse[self.aminoAcidCode[i]]=i
                    labelIndexToUse[self.aminoAcidCode[j]]=j
                    labelIndexToUse[self.aminoAcidCode[k]]=k
                    # First, Second and Third amino label list
                    self.initAminoLabel[0]=self.aminoAcidCode[i]
                    self.initAminoLabel[1]=self.aminoAcidCode[j]
                    self.initAminoLabel[2]=self.aminoAcidCode[k]
                    # Sorted labels from above list 
                    sortedAminoLabel=list(self.initAminoLabel)
                    # Reverse order from above sorted list
                    sortedAminoLabel.sort(reverse=True)

                    # The fourth case when l1=l2=l3
                    if (sortedAminoLabel[0] == sortedAminoLabel[1]) and (sortedAminoLabel[1]==sortedAminoLabel[2]):
                        distance1_2 = self.calDistance(i,j)
                        distance1_3 = self.calDistance(i,k)
                        distance2_3 = self.calDistance(j,k)
                        if distance1_2 >= (max(distance1_2,distance1_3,distance2_3)):
                            l1_index0=i
                            l2_index1=j
                            l3_index2=k
                        elif distance1_3 >= (max(distance1_2,distance1_3,distance2_3)):
                            l1_index0=i
                            l2_index1=k
                            l3_index2=j
                        else:
                            l1_index0=j
                            l2_index1=k
                            l3_index2=i

                    # Third condition when l1=l2>l3
                    elif(sortedAminoLabel[0]==sortedAminoLabel[1])and(sortedAminoLabel[1]!=sortedAminoLabel[2]):
                        l3_index2 = labelIndexToUse[sortedAminoLabel[2]]
                        indices = self.findTheIndex(l3_index2,i,j,k)
                        first = l3_index2
                        second = indices[0]
                        third  =indices[1]
                        distance1_3=self.calDistance(second,first)
                        distance2_3=self.calDistance(third,first)
                        if distance1_3>=distance2_3:
                            l1_index0=indices[0]
                            l2_index1=indices[1]	
                        else:
                            l1_index0=indices[1]
                            l2_index1=indices[0]

                    # Second condition when l1>l2=l3     
                    elif(sortedAminoLabel[0]!=sortedAminoLabel[1])and(sortedAminoLabel[1]==sortedAminoLabel[2]):
                        l1_index0=labelIndexToUse[sortedAminoLabel[0]]
                        indices = self.findTheIndex(l1_index0,i,j,k)
                        if self.calDistance(l1_index0,indices[0])>= self.calDistance(l1_index0,indices[1]):
                            l2_index1=indices[0]
                            l3_index2=indices[1]	
                        else:
                            l3_index2=indices[0]
                            l2_index1=indices[1]

                    # First condition when l1!=l2!=l3
                    elif(sortedAminoLabel[0]!=sortedAminoLabel[1])and(sortedAminoLabel[0]!=sortedAminoLabel[2])and(sortedAminoLabel[1]!=sortedAminoLabel[2]):
                        # Getting the index from the labelIndexToUse from sortedAminoLabel use
                        for index in range(0,3):
                            self.sortedAminoIndex[index]=labelIndexToUse[sortedAminoLabel[index]]
                        l1_index0=self.sortedAminoIndex[0]
                        l2_index1=self.sortedAminoIndex[1]
                        l3_index2=self.sortedAminoIndex[2]

                    distance01=self.calDistance(l1_index0,l2_index1)
                    # Calculating the mid distance
                    midDis01 = distance01/2
                    distance02=self.calDistance(l1_index0,l3_index2)
                    distance12=self.calDistance(l2_index1,l3_index2)
                    # Calculating the max distance (D)
                    maxDistance=max(distance01,distance02,distance12)
                    # Calculating the mid point 
                    m1 = (self.xCoordinate[l1_index0]+ self.xCoordinate[l2_index1])/2
                    m2 = (self.yCoordinate[l1_index0]+ self.yCoordinate[l2_index1])/2
                    m3 = (self.zCoordinate[l1_index0]+ self.zCoordinate[l2_index1])/2

                    # Calculating the d3 distance
                    d3 = math.sqrt((m1 - self.xCoordinate[l3_index2])**2+(m2 - self.yCoordinate[l3_index2])**2+(m3 - self.zCoordinate[l3_index2])**2)

                    # Calculating thetaAngle1
                    thetaAngle1 = 180*(math.acos((distance02**2-midDis01**2-d3**2)/(2*midDis01*d3)))/3.14

                    # Check in which category does the angle falls
                    if thetaAngle1<=90:
                        theta = thetaAngle1
                    else:
                        theta = abs(180-thetaAngle1)

                    # Calculating the bin values for theta and max distance
                    binTheta = self.thetaClass(theta)
                    binLength = self.dist12Class(maxDistance)
                    
                    # These are the residue of the protein, the amino three amino acids
                    aminoAcidR1 = list(self.aminoAcidLabelWithCode.keys())[list(self.aminoAcidLabelWithCode.values()).index(self.aminoAcidCode[l1_index0])]
                    aminoAcidR2 = list(self.aminoAcidLabelWithCode.keys())[list(self.aminoAcidLabelWithCode.values()).index(self.aminoAcidCode[l2_index1])]
                    aminoAcidR3 = list(self.aminoAcidLabelWithCode.keys())[list(self.aminoAcidLabelWithCode.values()).index(self.aminoAcidCode[l3_index2])]
                    # These are the sequence number of the three amino acids
                    seqNumber1 = list(self.aminoSeqNum.values())[l1_index0]
                    seqNumber2 = list(self.aminoSeqNum.values())[l2_index1]
                    seqNumber3 = list(self.aminoSeqNum.values())[l3_index2]
                    # These are the coordinates of the three amino acids
                    aminoAcidC10, aminoAcidC11, aminoAcidC12 = self.xCoordinate[l1_index0], self.yCoordinate[l1_index0], self.zCoordinate[l1_index0]
                    aminoAcidC20, aminoAcidC21, aminoAcidC22 = self.xCoordinate[l2_index1], self.yCoordinate[l2_index1], self.zCoordinate[l2_index1]
                    aminoAcidC30, aminoAcidC31, aminoAcidC32 = self.xCoordinate[l3_index2], self.yCoordinate[l3_index2], self.zCoordinate[l3_index2]

                    # Calculating the triplets key value
                    tripletKeys = dLen*dtheta*(numOfLabels**2)*(self.aminoAcidCode[l1_index0]-1)+dLen*dtheta*(numOfLabels)*(self.aminoAcidCode[l2_index1]-1)+dLen*dtheta*(self.aminoAcidCode[l3_index2]-1)+dtheta*(binLength-1)+(binTheta-1)
                    # print("Key Triplets", tripletKeys)

                    # Total number of keys and max distance list
                    self.totalKeys.append(tripletKeys)
                    self.maxDistList.append(maxDistance)

                    # Filtering out the distinct keys
                    if tripletKeys in self.keyFreq:
                        self.keyFreq[tripletKeys]+=1
                    else:
                        self.keyFreq[tripletKeys] = 1

                    # These are the info of all the triplets
                    tripletInfoAll = (str(tripletKeys)+"\t"+str(aminoAcidR1)+"\t"+str(seqNumber1)+"\t"+str(aminoAcidR2)+"\t"+str(seqNumber2)+"\t"+str(aminoAcidR3)+"\t"+str(seqNumber3)+"\t"+str(binTheta)+"\t"+str(theta)+"\t"+str(binLength)+"\t"+str(maxDistance)+"\t"+str(aminoAcidC10)+"\t"+str(aminoAcidC11)+"\t"+str(aminoAcidC12)+"\t"+str(aminoAcidC20)+"\t"+str(aminoAcidC21)+"\t"+str(aminoAcidC22)+"\t"+str(aminoAcidC30)+"\t"+str(aminoAcidC31)+"\t"+str(aminoAcidC32)+"\n")
                    tripletsFile.writelines(tripletInfoAll)

        # Storing the distinct keys in a file
        for values in self.keyFreq:
            keyFreqFile.writelines([str(values), '\t', str(self.keyFreq[values]), "\n"]) 

# Usage
if __name__ == "__main__":
    dtheta = 29
    dLen = 35
    numOfLabels = 20
    analyzer = AminoAcidAnalyzer(dtheta, dLen, numOfLabels)
    # Generate list of proteins from csv file containing list of proteins and chain
    analyzer.readCSVProteinChain(CSV_FILE_PATH)
    # analyzer.downloadDataSet()

    numOfCores = multiprocessing.cpu_count()
    print(numOfCores)
##### This is to write in a csv file
    csv_file_path = "proteinKeysDist.csv"
    header = ["Protein", "Chain", "#amino_acids", "#keys", "#keys_with_freq", "max_distance", "min_distance"]

    def calculate_data(fileName):
        analyzer.calcuTheteAndKey(fileName, analyzer.chain_dict[fileName].upper())
        aminoAcids = str(len(analyzer.aminoAcidCode))
        totalKeys = str(len(analyzer.totalKeys))
        keysWithFreq = str(len(analyzer.keyFreq))
        maxDistance = str(max(analyzer.maxDistList))
        minDistance = str(min(analyzer.maxDistList))
        data = [[fileName, analyzer.chain_dict[fileName].upper(), aminoAcids, totalKeys, keysWithFreq, maxDistance, minDistance]]

        with open(csv_file_path, "a", newline="") as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerows(data)

# Writing the header to the CSV file
with open(csv_file_path, "w", newline="") as csvfile:
    csv_writer = csv.writer(csvfile)
    csv_writer.writerow(header)

# Parallel processing of the protein list
results = Parallel(n_jobs=numOfCores, verbose=50)(delayed(calculate_data)(fileName) for fileName in analyzer.proteinList)