#!/usr/bin/env python -tt

"""
This program prints the contacts for every atom.
The threshold should be specified in Angstroms, and can be a decimal.
"""

import sys
import math

#This method will read the PDB file and extract the atomic information.
def readPDBFile(fileName):
    
    #This matrix will hold all five pieces of information for each atom.
    valuesMatrix = []
    
    #Open the file
    f = open(filename, 'rU')
    
    #Walk through the file line by line
    for line in f:
        #Removes whitespace and creates an array containing each column entry
        splitLine = line.split()
        #Checks to make sure the line being read contains atomic information, and that the atom isn't a hydrogen
        if splitLine[0] == 'ATOM' and not splitLine[2].startswith('H'):
            #Puts the following information into a small array: atom type, residue number, and x, y, and z, coordinates, in that order.
            atomicInfo = [splitLine[2], splitLine[5], splitLine [6], splitLine [7], splitLine[8]]
            #Convert the residue number from a string to an integer
            try:
                atomicInfo[1] = int(atomicInfo[1])
            except ValueError:
                print "Hey, this residue number isn't a number at all!"
                print "You screwed up at: " + atomicInfo[1]
                sys.exit(1)
            #Convert coordinates from strings to floating point numbers
            try:
                atomicInfo[2] = float(atomicInfo[2])
                atomicInfo[3] = float(atomicInfo[3])
                atomicInfo[4] = float(atomicInfo[4])
            except ValueError:
                print "These coordinates are all screwed up!"
                print "Here they are:", atomicInfo[2], atomicInfo[3], atomicInfo[4]
            #Puts the information about this atom into the overall matrix
            valuesMatrix.append(atomicInfo)
        #If multiple models are contained within the same pdb file, consider only the first one.
        if splitLine[0] == 'ENDMDL':
            break
    
    #Pass the matrix back to the main method
    return valuesMatrix

def calc3dDistance(x1, y1, z1, x2, y2, z2):
    #Calculates the distance between the atoms with coordinates (x1, y1, z1) and (x2, y2, z2)
    return math.sqrt(math.pow(x1-x2, 2) + math.pow(y1-y2, 2) + math.pow(z1-z2, 2))

def calculateDistances(matrix):
    #Matrix containing each atom's distance to every other atom
    distanceInfo = []
    #Steps through each row in the matrix and calculates distances for that row
    for row in matrix:
        #Adds in atom type and residue number
        atomDistances = [row[0], row[1]]
        #Calculates distances between the current row and every other row (including the current row)
        for row2 in matrix:
            #Calculates the actual distance
            distance = calc3dDistance(row[2], row[3], row[4], row2[2], row2[3], row2[4])
            #Adds the distance info to this atom
            atomDistances.append(distance)
        #Adds all the distance information about this atom to the overall matrix
        distanceInfo.append(atomDistances)
    return distanceInfo

def findClosestContacts(matrix, threshold):
    #Contains the following info: residue number, atom type, and tuples containing (residue number, atom type, distance)
    contacts = []
    #Check for contacts atom by atom
    for row in matrix:
        #Add residue number and atom type to the start of the row
        contactsList = [row[1], row[0]]
        #Keeps track of the index being iterated over
        columnTracker = 0
        #Iterates over the distances to every other atom
        for element in row[2:]:
            #If the distance is less than the threshold distance, and the residues are more than four apart
            if element <= threshold and abs(matrix[columnTracker][1] - row[1]) > 4:
                #Add the tuple res#, atom type, distance
                contactsList.append((matrix[columnTracker][1], matrix[columnTracker][0], element))
            #Increment the index counter
            columnTracker += 1
        #Add info about this atom to the overall container
        contacts.append(contactsList)
    return contacts



def main():
    
    #if not enough or too many command line arguments are given, crash
    if len(sys.argv) != 3:
        print 'usage: python TestReadPDB.py file contactThreshold'
        sys.exit(1)

    #The first command line argument is the filename
    global filename
    filename = sys.argv[1]
    #The second command line argument is the threshold
    global threshold
    threshold = sys.argv[2]
    print("File name:", filename)
    try:
        threshold = float(threshold)
    except ValueError:
        print(sys.argv[2], "is not a number you fool")
        print("Please specify the threshold in Angstroms!")
        sys.exit(1)
    
    #atomicInfo contains atomic information
    atomicInfo = readPDBFile(filename)
    print("Done reading file...")
    
    #Use atomicInfo to calculate distances to other atoms
    distances = calculateDistances(atomicInfo)
    print("Done calculating all distances...")
    
    #Use distances to find the contacts
    matrix = findClosestContacts(distances, threshold)
    
    #Print it out
    print("Syntax: residue number, atom type, contact 1 res#, contact 1 type, contact 1 distance, etc")
    for row in matrix:
        print row
    
    print('All done!')
    sys.exit(0)


if __name__ == '__main__':
  main()