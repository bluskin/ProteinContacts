#!/usr/bin/env python

"""
This program is intended to read PDB files and extract the following information
about atoms: atom type, residue number, and x, y, and z coordinates.
Created by Dan Baum and Ben Luskin, 6/15/11
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



def main():
    
    #if not enough or too many command line arguments are given, crash
    if len(sys.argv) != 2:
        print 'usage: python TestReadPDB.py file'
        sys.exit(1)

    #The first command line argument is the filename
    global filename
    filename = sys.argv[1]
    print(filename)
    
    #Matrix contains atomic information
    matrix = readPDBFile(filename)
    
    #Print out matrix
    print("Atom type", "Residue number", "x", "y", "z")
    for element in matrix:
        print element
    
    print('All done!')
    sys.exit(0)


if __name__ == '__main__':
  main()