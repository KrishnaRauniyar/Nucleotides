import math
import csv
from utils.theta_utils import thetaClass
from utils.distance_utils import dist12Class
from joblib import Parallel, delayed
import argparse

class CrossKeys:
    def __init__(self):
        self.xCoordinate = {}
        self.yCoordinate = {}
        self.zCoordinate = {}
        self.totalKeys = []
        self.maxDistList = []
        self.keyFreq = {}
        self.newDictForLexicalCode = {}
        self.aminoSeqNum = {}
        self.initAminoLabel = [0, 0, 0]
        self.sortedAminoIndex = [0, 0, 0]
        self.drugAminoSeqNum = set()
        self.getAxisSeq = {}
        self.residueInfo = {}
        self.chainInfo = {}

    def thetaClass(self, theta):
        return thetaClass(theta)

    def dist12Class(self, dist12):
        return dist12Class(dist12)
    
    # Calculating the distance between the amino acids 
    def calDistance(self, l1_index, l2_index):
        x1 = self.xCoordinate[l1_index]
        x2 = self.xCoordinate[l2_index]
        y1 = self.yCoordinate[l1_index]
        y2 = self.yCoordinate[l2_index]
        z1 = self.zCoordinate[l1_index]
        z2 = self.zCoordinate[l2_index]
        return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

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

    # Function to extract sequence number from the third part of the residue string
    def extract_sequence_number(self, residue_string):
        parts = residue_string.split('_')
        return int(parts[-2])
    
    def extract_residue(self, residue_string):
        parts = residue_string.split('_')
        return parts[-3]
    
    def extract_chain(self, residue_string):
        parts = residue_string.split('_')
        return parts[-4]
    
        # Function to extract sequence number from the third part of the residue string
    def extract_drug_seq(self, residue_string):
        parts = residue_string.split('_')
        return f"{parts[-2]}_{parts[0]}_{parts[-4]}_{parts[-3]}"
    
    def getDrugSeqNum(self, file):
         with open(file, mode='r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip header row
            # Process each row in the CSV file
            for index, row in enumerate(reader):
                # Extract the residue sequence numbers
                drug_residue_num = self.extract_drug_seq(row[0])
                self.drugAminoSeqNum.add(drug_residue_num)

    # Read and process the CSV data
    def read_csv_data(self, filename, fileSeq):
         with open(filename, mode='r') as file:
            reader = csv.reader(file)
            next(reader)  # Skip header row
            # Process each row in the CSV file
            for index, row in enumerate(reader):
                if(self.extract_drug_seq(row[0]) == fileSeq):
                    # Extract atom names and sequence numbers
                    drug_atom = row[0].split('_')[-1]
                    drug_seq = int(row[1])
                    drug_coords = list(map(float, row[2].split('_')))

                    protein_atom = row[3].split('_')[-1]
                    protein_seq = int(row[4])
                    protein_coords = list(map(float, row[5].split('_')))

                    # Extract the residue sequence numbers
                    drug_residue_seq_num = self.extract_sequence_number(row[0])
                    protein_residue_seq_num = self.extract_sequence_number(row[3])

                    # Extract the residue 
                    drug_residue = self.extract_residue(row[0])
                    protein_residue = self.extract_residue(row[3])

                    # Extract the chain 
                    drug_chain = self.extract_chain(row[0])
                    protein_chain = self.extract_chain(row[3])

                    # Add entries to the dictionaries
                    self.getAxisSeq[f"{drug_atom}_D_cordinateds_seqnum"] = drug_coords + [drug_residue_seq_num, drug_seq, drug_residue, drug_chain] 
                    self.getAxisSeq[f"{protein_atom}_P_cordinateds_seqnum"] = protein_coords + [protein_residue_seq_num, protein_seq, protein_residue, protein_chain] 

    def create_separate_dictionaries(self):
        self.newDictForLexicalCode = {}
        self.xCoordinate = {}
        self.yCoordinate = {}
        self.zCoordinate = {}
        self.aminoSeqNum = {}
        self.residueInfo = {}
        self.chainInfo = {}
        increment = 0
        for key, value in self.getAxisSeq.items():
            atom = '_'.join(key.split('_')[:2])
            self.newDictForLexicalCode[atom] = value[-3]
            self.xCoordinate[increment] = value[0]
            self.yCoordinate[increment] = value[1]
            self.zCoordinate[increment] = value[2]
            self.aminoSeqNum[increment] = value[3]
            self.residueInfo[increment] = value[5]
            self.chainInfo[increment] = value[6]
            increment+=1
        ### delete this so that it does not store all the data when looped , delete after every loop for each calculation
        self.getAxisSeq = {}

    # Calculate keys based on your logic
    def calculatetheta(self, fileSeq, output_path):
        if len(self.newDictForLexicalCode) >= 3:
            self.keyFreq = {}
            ### For real script
            tripletsFile = open(f"{output_path}/{fileSeq}.keys_theta29_dist18", "w")
            keyFreqFile = open(f"{output_path}/{fileSeq}.keys_Freq_theta29_dist18", "w")
            for i in range(0, len(self.newDictForLexicalCode) - 2):
                for j in range(i + 1, len(self.newDictForLexicalCode) - 1):
                    for k in range(j + 1, len(self.newDictForLexicalCode)):
                        # Check if all three atoms are not from the same type (drug or protein)
                        if ("_D" in list(self.newDictForLexicalCode.keys())[i] and "_D" in list(self.newDictForLexicalCode.keys())[j] and "_D" in list(self.newDictForLexicalCode.keys())[k]) or ("_P" in list(self.newDictForLexicalCode.keys())[i] and "_P" in list(self.newDictForLexicalCode.keys())[j] and "_P" in list(self.newDictForLexicalCode.keys())[k]):
                            continue
                        # This is a dictionary to keep the index and the labels
                        labelIndexToUse = {}
                        # First, Second and Third label and Index
                        labelIndexToUse[list(self.newDictForLexicalCode.values())[i]] = i
                        labelIndexToUse[list(self.newDictForLexicalCode.values())[j]] = j
                        labelIndexToUse[list(self.newDictForLexicalCode.values())[k]] = k
                        # First, Second and Third amino label list
                        self.initAminoLabel[0] = list(self.newDictForLexicalCode.values())[i]
                        self.initAminoLabel[1] = list(self.newDictForLexicalCode.values())[j]
                        self.initAminoLabel[2] = list(self.newDictForLexicalCode.values())[k]
                        # Sorted labels from above list
                        sortedAminoLabel = list(self.initAminoLabel)
                        # Reverse order from above sorted list
                        sortedAminoLabel.sort(reverse=True)

                        # The fourth case when l1=l2=l3
                        if (sortedAminoLabel[0] == sortedAminoLabel[1]) and (sortedAminoLabel[1] == sortedAminoLabel[2]):
                            distance1_2 = self.calDistance(i, j)
                            distance1_3 = self.calDistance(i, k)
                            distance2_3 = self.calDistance(j, k)
                            if distance1_2 >= (max(distance1_2, distance1_3, distance2_3)):
                                l1_index0 = i
                                l2_index1 = j
                                l3_index2 = k
                            elif distance1_3 >= (max(distance1_2, distance1_3, distance2_3)):
                                l1_index0 = i
                                l2_index1 = k
                                l3_index2 = j
                            else:
                                l1_index0 = j
                                l2_index1 = k
                                l3_index2 = i

                        # Third condition when l1=l2>l3
                        elif (sortedAminoLabel[0] == sortedAminoLabel[1]) and (sortedAminoLabel[1] != sortedAminoLabel[2]):
                            l3_index2 = labelIndexToUse[sortedAminoLabel[2]]
                            indices = self.findTheIndex(l3_index2, i, j, k)
                            first = l3_index2
                            second = indices[0]
                            third = indices[1]
                            distance1_3 = self.calDistance(second, first)
                            distance2_3 = self.calDistance(third, first)
                            if distance1_3 >= distance2_3:
                                l1_index0 = indices[0]
                                l2_index1 = indices[1]
                            else:
                                l1_index0 = indices[1]
                                l2_index1 = indices[0]

                        # Second condition when l1>l2=l3
                        elif (sortedAminoLabel[0] != sortedAminoLabel[1]) and (sortedAminoLabel[1] == sortedAminoLabel[2]):
                            l1_index0 = labelIndexToUse[sortedAminoLabel[0]]
                            indices = self.findTheIndex(l1_index0, i, j, k)
                            if self.calDistance(l1_index0, indices[0]) >= self.calDistance(l1_index0, indices[1]):
                                l2_index1 = indices[0]
                                l3_index2 = indices[1]
                            else:
                                l3_index2 = indices[0]
                                l2_index1 = indices[1]

                        # First condition when l1!=l2!=l3
                        elif (sortedAminoLabel[0] != sortedAminoLabel[1]) and (sortedAminoLabel[0] != sortedAminoLabel[2]) and (
                                sortedAminoLabel[1] != sortedAminoLabel[2]):
                            # Getting the index from the labelIndexToUse from sortedAminoLabel use
                            for index in range(0, 3):
                                self.sortedAminoIndex[index] = labelIndexToUse[sortedAminoLabel[index]]
                            l1_index0 = self.sortedAminoIndex[0]
                            l2_index1 = self.sortedAminoIndex[1]
                            l3_index2 = self.sortedAminoIndex[2]

                        distance01 = self.calDistance(l1_index0, l2_index1)
                        # Calculating the mid distance
                        midDis01 = distance01 / 2
                        distance02 = self.calDistance(l1_index0, l3_index2)
                        distance12 = self.calDistance(l2_index1, l3_index2)
                        # Calculating the max distance (D)
                        maxDistance = max(distance01, distance02, distance12)
                        # Calculating the mid point 
                        m1 = (self.xCoordinate[l1_index0] + self.xCoordinate[l2_index1]) / 2
                        m2 = (self.yCoordinate[l1_index0] + self.yCoordinate[l2_index1]) / 2
                        m3 = (self.zCoordinate[l1_index0] + self.zCoordinate[l2_index1]) / 2

                        # Calculating the d3 distance
                        d3 = math.sqrt((m1 - self.xCoordinate[l3_index2])**2 + (m2 - self.yCoordinate[l3_index2])**2 + (m3 - self.zCoordinate[l3_index2])**2)

                        # Calculating thetaAngle1
                        thetaAngle1 = 180 * (math.acos((distance02**2 - midDis01**2 - d3**2) / (2 * midDis01 * d3))) / 3.14

                        # Check in which category does the angle falls
                        if thetaAngle1 <= 90:
                            theta = thetaAngle1
                        else:
                            theta = abs(180 - thetaAngle1)

                        # Calculating the bin values for theta and max distance
                        binTheta = self.thetaClass(theta)
                        binLength = self.dist12Class(maxDistance)
                        
                        # Additional code for extracting amino acid information
                        aminoAcidR1 = list(self.newDictForLexicalCode.keys())[l1_index0]
                        aminoAcidR2 = list(self.newDictForLexicalCode.keys())[l2_index1]
                        aminoAcidR3 = list(self.newDictForLexicalCode.keys())[l3_index2]

                        # These are the sequence number of the three amino acids
                        seqNumber1 = list(self.aminoSeqNum.values())[l1_index0]
                        seqNumber2 = list(self.aminoSeqNum.values())[l2_index1]
                        seqNumber3 = list(self.aminoSeqNum.values())[l3_index2]

                        # These are the residue type
                        residue1 = list(self.residueInfo.values())[l1_index0]
                        residue2 = list(self.residueInfo.values())[l2_index1]
                        residue3 = list(self.residueInfo.values())[l3_index2]

                        # These are the chain type
                        chain1 = list(self.chainInfo.values())[l1_index0]
                        chain2 = list(self.chainInfo.values())[l2_index1]
                        chain3 = list(self.chainInfo.values())[l3_index2]

                        # These are the coordinates of the three amino acids
                        aminoAcidC10, aminoAcidC11, aminoAcidC12 = self.xCoordinate[l1_index0], self.yCoordinate[l1_index0], self.zCoordinate[l1_index0]
                        aminoAcidC20, aminoAcidC21, aminoAcidC22 = self.xCoordinate[l2_index1], self.yCoordinate[l2_index1], self.zCoordinate[l2_index1]
                        aminoAcidC30, aminoAcidC31, aminoAcidC32 = self.xCoordinate[l3_index2], self.yCoordinate[l3_index2], self.zCoordinate[l3_index2]

                        # Calculating the triplets key value
                        tripletKeys = dLen * dtheta * (numOfLabels**2) * (self.newDictForLexicalCode[aminoAcidR1] - 1) + dLen * dtheta * (numOfLabels) * (self.newDictForLexicalCode[aminoAcidR2] - 1) + dLen * dtheta * (self.newDictForLexicalCode[aminoAcidR3] - 1) + dtheta * (binLength - 1) + (binTheta - 1)

                        # Total number of keys and max distance list
                        self.totalKeys.append(tripletKeys)
                        self.maxDistList.append(maxDistance)

                        # Filtering out the distinct keys
                        if tripletKeys in self.keyFreq:
                            self.keyFreq[tripletKeys] += 1
                        else:
                            self.keyFreq[tripletKeys] = 1

                        # These are the info of all the triplets
                        tripletInfoAll = (str(tripletKeys) + "\t" + str(aminoAcidR1) + "\t" + str(seqNumber1) + "\t" + str(residue1) + "\t" + str(chain1) + "\t" + str(aminoAcidR2) + "\t" + str(seqNumber2) + "\t" + str(residue2) + "\t" + str(chain2) + "\t" + str(aminoAcidR3) + "\t" + str(seqNumber3) + "\t" + str(residue3) + "\t" + str(chain3) + "\t"+ str(binTheta) + "\t" + str(theta) + "\t" + str(binLength) + "\t" + str(maxDistance) + "\t" + str(aminoAcidC10) + "\t" + str(aminoAcidC11) + "\t" + str(aminoAcidC12) + "\t" + str(aminoAcidC20) + "\t" + str(aminoAcidC21) + "\t" + str(aminoAcidC22) + "\t" + str(aminoAcidC30) + "\t" + str(aminoAcidC31) + "\t" + str(aminoAcidC32) + "\n")
                        tripletsFile.writelines(tripletInfoAll)

            # Storing the distinct keys in a file
            for values in self.keyFreq:
                keyFreqFile.writelines([str(values), '\t', str(self.keyFreq[values]), "\n"])
        else:
            print(f"Cannot make keys for {fileSeq} having only {len(self.newDictForLexicalCode)} atoms.")

def process_file_seq(cross_keys, fileSeq, csv_filename, output_path):
    cross_keys.read_csv_data(csv_filename, fileSeq)
    cross_keys.create_separate_dictionaries()
    cross_keys.calculatetheta(fileSeq, output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Drug-protein cross key calculation")
    parser.add_argument('-p', '--input_path', type=str, required=True, help="path to csv file")
    parser.add_argument('-o', '--output_path', type=str, required=True, help="path to output folder")
    args = parser.parse_args()
    dtheta = 29
    dLen = 18
    numOfLabels = 112
    cross_keys = CrossKeys()
    csv_filename  = args.input_path
    output_path = args.output_path
    cross_keys.getDrugSeqNum(file=csv_filename)

    # cross_keys.drugAminoSeqNum = ['2_2AC0_DG_F', '501_3VD2_DT_F']
    # for fileSeq in cross_keys.drugAminoSeqNum:
    #     cross_keys.read_csv_data(csv_filename, fileSeq)
    #     cross_keys.create_separate_dictionaries()
    #     cross_keys.calculatetheta(fileSeq, output_path)

    Parallel(n_jobs=-1, verbose=5)(delayed(process_file_seq)(cross_keys, fileSeq, csv_filename, output_path) for fileSeq in cross_keys.drugAminoSeqNum)

