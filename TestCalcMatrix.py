#!/usr/bin/env python

"""
This program segment calculates a matrix showing the distance from every atom
to every other atom.
"""

import sys

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
        #Checks to make sure the line being read contains atomic information
        if splitLine[0] == 'ATOM':
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
    
    #Pass the matrix back to the main method
    return valuesMatrix

def calc3dDistance(x1, y1, z1, x2, y2, z2):
    import math
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

def main():
    
    #if not enough or too many command line arguments are given, crash
    if len(sys.argv) != 2:
        print 'usage: python TestReadPDB.py file'
        sys.exit(1)

    #The first command line argument is the filename
    global filename
    filename = sys.argv[1]
    print(filename)
    
    #AtomicInfo contains atomic information
    atomicInfo = readPDBFile(filename)
    
    #Use AtomicInfo to calculate distances to other atoms
    matrix = calculateDistances(atomicInfo)
    
    #Print it out
    print('atom type', 'residue number', 'distance to atom 1', 'to atom 2', 'etc')
    for row in matrix:
        print row
    
    print('All done!')
    sys.exit(0)


if __name__ == '__main__':
  main()