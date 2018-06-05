#!/usr/bin/env python -tt

"""
This program calculates contact orders.
It disregards all hydrogen atoms.
The command-line argument 'threshold' should be in Angstroms.
This is the BEN TEST PROGRAM
"""

import sys
import math
import csv

#This method will read the PDB file and extract the atomic information.
def readPDBFile(fileName):
    
    #This matrix will hold all five pieces of information for each atom.
    valuesMatrix = []
    
    #This matrix will hold information about helical content.
    helixMatrix = []
    
    #This matrix will hold information about sheet content.
    sheetMatrix = []
    
    #Open the file
    try:
        f = open(filename, 'rU')
    except IOError:
        print 'Whoops! Looks like \'' + filename +  '\' is either corrupt or nonexistent.'
        sys.exit(1)
        
    
    #Walk through the file line by line
    for line in f:
        #Removes whitespace and creates an array containing each column entry
        splitLine = line.split()
        
        #Checks to see if the line being read contains information about helix locations
        if splitLine[0] == 'HELIX':
            #This index contains the residue number where the helix starts
            start = splitLine[5]
            #This index contains the residue number where the helix ends
            end = splitLine[8]
            #Convert the residue numbers to integers and store them as a tuple: (start, end)
            try:
                helixMatrix.append((int(start), int(end)))
            except ValueError:
                print 'You got stuff here in the helix that ain\'t numbers, pal.'
                print 'Start residue:', start, 'End residue:', end
                sys.exit(1)
            
        #Checks to see if the line being read contains information about sheet locations
        if splitLine[0] == 'SHEET':
            #This index contains the residue number where the sheet starts
            start = splitLine[6]
            #This index contains the residue number where the sheet ends
            end = splitLine[9]
            #Convert the residue numbers to integers and store them as a tuple: (start, end)
            try:
                sheetMatrix.append((int(start), int(end)))
            except ValueError:
                print 'You got stuff here in the sheet that ain\'t numbers, pal.'
                print 'Start residue:', start, 'End residue:', end
                sys.exit(1)
                
        #Checks to see if the line being read contains atomic information, and that the atom isn't a hydrogen
        if splitLine[0] == 'ATOM' and not splitLine[2].startswith('H'):
            
            #If two entries run into each other, all the indexes will be off by one.
            if len(splitLine[2]) > 3:
                splitLine.insert(3, 'SPACER')
            
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
    
    #Pass the matrices back to the main method
    return valuesMatrix, helixMatrix, sheetMatrix

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
    #Maps residue numbers to their contacts
    contactsResiduesOnly = {}
    #Check for contacts atom by atom
    for row in matrix:
        #Add residue number and atom type to the start of the row
        contactsList = [row[1], row[0]]
        #Don't add anything to this list yet
        contactsListResiduesOnly = []
        #Keeps track of the index being iterated over
        columnTracker = 0
        #Iterates over the distances to every other atom
        for element in row[2:]:
            #If the distance is less than the threshold distance, and the residues are more than four apart
            if element <= threshold and abs(matrix[columnTracker][1] - row[1]) > 4:
                #Add the tuple res#, atom type, distance
                contactsList.append((matrix[columnTracker][1], matrix[columnTracker][0], element))
                #Add the res# to the list containing only res#s
                contactsListResiduesOnly.append(matrix[columnTracker][1])
            #Increment the index counter
            columnTracker += 1
        #Add info about this atom to the overall container
        contacts.append(contactsList)
        #Map the residue number to the list of residues it contacts
        if row[1] in contactsResiduesOnly:
            contactsResiduesOnly[row[1]].extend(contactsListResiduesOnly)
            contactsResiduesOnly[row[1]] = sorted(contactsResiduesOnly[row[1]])
        else:
            contactsResiduesOnly[row[1]] = contactsListResiduesOnly[:]
    return contacts, contactsResiduesOnly

def countContacts(dictOfContacts):
    #This dict will contain a mapping from residue numbers to the number of times that residue contacts every other residue
    contacts = {}
    #This line finds the total number of residues
    numberOfLines = sorted(dictOfContacts.keys())[-1]
    #For each residue, count the number of contacts
    for residue in dictOfContacts:
        #Has length equal to the number of residues, and each index corresponds to the number of contacts with that residue
        row = []
        #Since python uses zero-based counting, but PDB uses one-based counting, we need count(i+1)
        for i in range(0, numberOfLines):
            row.append(dictOfContacts[residue].count(i+1))
        contacts[residue] = row
    return contacts

#For the RCO formula, see Kurt, Mounce, Ellison, and Cavagnero, Biotechnology Progress, 24(3): 570-575 (2008)
def calcRCO(dictOfContacts):
    #Will map each residue to its residue-specific contact order
    dictOfRCOs = {}
    #L is the length of the protein
    global L
    L = len(dictOfContacts.keys())
    for residue in dictOfContacts:
        #listOfContacts contains info about how often this residue contacts every other residue
        listOfContacts = dictOfContacts[residue]
        #Ni is the total number of contacts this residue forms
        Ni = sum(listOfContacts)
        sigma = 0
        for i in range(0, len(listOfContacts)):
            #Sij is the residue separation between contacting residues
            Sij = abs(i + 1 - residue)
            #Nij is the number of contacts between residues i and j
            Nij = listOfContacts[i]
            sigma += Nij * Sij
        #This formula won't work when no contacts are formed, so define RCO to be 0 in that case
        try:
            dictOfRCOs[residue] = float(sigma) / (L * Ni)
        except ZeroDivisionError:
            dictOfRCOs[residue] = 0
    return dictOfRCOs
    
# Given a list of numbers, returns a list where
# all adjacent == elements have been reduced to a single element,
# so [1, 2, 2, 3] returns [1, 2, 3].
def remove_adjacent(nums):
  result = []
  for num in nums:
    if len(result) == 0 or num != result[-1]:
      result.append(num)
  return result

#For the RCB formula, see Kurt, Mounce, Ellison, and Cavagnero, Biotechnology Progress, 24(3): 570-575 (2008)
def calcRCB(residueContacts):
    #Maps each residue to its residue-specific contact breadth
    dictOfRCBs = {}
    #L is the length of the protein
    L = len(residueContacts.keys())
    #Calculate RCB for each residue
    for residue in residueContacts:
        #listOfResidues says which residues this one contacts
        listOfResidues = remove_adjacent(residueContacts[residue])
        #n is the total number of residues this one contacts
        n = len(listOfResidues)
        #Standard deviations aren't defined for less than 2 items - set this case to be 0
        if n == 0 or n == 1:
            dictOfRCBs[residue] = 0
            continue
        #Rbar is the average of all the residue numbers this one contacts
        Rbar = float(sum(listOfResidues)) / n
        sigma = 0.0
        #sigma is defined in the RCB calculation - see citation
        for contact in listOfResidues:
            sigma += math.pow(Rbar - contact, 2)
        dictOfRCBs[residue] = 1/float(L) * math.sqrt(sigma / (n-1))
    return dictOfRCBs

def calcAverages(RCO):
    NTerminusAverages = []
    middleAverages = []
    CTerminusAverages = []
    global tenPercent
    tenPercent = int(0.1 * len(RCO.keys()))
    listOfValues = []
    for residue in sorted(RCO.keys())[:tenPercent]:
        listOfValues.append(RCO[residue])
    NTerminusAverages.append(float(sum(listOfValues)) / len(listOfValues))
    listOfValues = []
    for residue in sorted(RCO.keys())[:(2 * tenPercent)]:
        listOfValues.append(RCO[residue])
    NTerminusAverages.append(float(sum(listOfValues)) / len(listOfValues))
    listOfValues = []
    for residue in sorted(RCO.keys())[:(3 * tenPercent)]:
        listOfValues.append(RCO[residue])
    NTerminusAverages.append(float(sum(listOfValues)) / len(listOfValues))
    listOfValues = []
    for residue in sorted(RCO.keys())[-tenPercent:]:
        listOfValues.append(RCO[residue])
    CTerminusAverages.append(float(sum(listOfValues)) / len(listOfValues))
    listOfValues = []
    for residue in sorted(RCO.keys())[(-2 * tenPercent):]:
        listOfValues.append(RCO[residue])
    CTerminusAverages.append(float(sum(listOfValues)) / len(listOfValues))
    listOfValues = []
    for residue in sorted(RCO.keys())[(-3 * tenPercent):]:
        listOfValues.append(RCO[residue])
    CTerminusAverages.append(float(sum(listOfValues)) / len(listOfValues))
    listOfValues = []
    for residue in sorted(RCO.keys())[tenPercent:-tenPercent]:
        listOfValues.append(RCO[residue])
    middleAverages.append(float(sum(listOfValues)) / len(listOfValues))
    listOfValues = []
    for residue in sorted(RCO.keys())[(2 * tenPercent):(-2 * tenPercent)]:
        listOfValues.append(RCO[residue])
    middleAverages.append(float(sum(listOfValues)) / len(listOfValues))
    listOfValues = []
    for residue in sorted(RCO.keys())[(3 * tenPercent):(-3 * tenPercent)]:
        listOfValues.append(RCO[residue])
    middleAverages.append(float(sum(listOfValues)) / len(listOfValues))
    return NTerminusAverages, middleAverages, CTerminusAverages

def calcSecondary(secondaryInfo):
    
    #The number of residues that belong to this type of secondary structure
    
    numResiduesN = 0
    numResiduesC = 0
    
    for tuple in secondaryInfo:
        #This part is for the N terminus.
        #If the secondary structure ends before the first ten percent, count the entire thing
        if tuple[1] <= tenPercent:
            #Add one because end - start doesn't count the very first residue (it gets subtracted off)
            #For example, if residues 1-10 were a helix, (10 - 1) would only count 9 residues, when actually there are 10
            numResiduesN += tuple[1] - tuple[0] + 1
        #If the ten percent mark occurs in the middle of the secondary structure element, only count up to the ten percent mark
        if tuple[1] > tenPercent and tuple[0] <= tenPercent:
            numResiduesN += tenPercent - tuple[0] + 1
        #This part is for the C terminus.
        #If the secondary structure starts after the last ten percent, count the entire thing
        if tuple[0] >= (L - tenPercent):
            numResiduesC += tuple[1] - tuple[0] + 1
        #If the secondary structure starts before the last ten percent, but crosses over, only count the part after the ten percent mark
        if tuple[0] < (L - tenPercent) and tuple[1] >= (L - tenPercent):
            numResiduesC += tuple[1] - (L - tenPercent) + 1
            
            
    percentSecondaryN = numResiduesN * 100 / float(tenPercent)
    percentSecondaryC = numResiduesC * 100 / float(tenPercent)
    
    return percentSecondaryN, percentSecondaryC
    
def shape(NTerminusAverages, middleAverages, CTerminusAverages):
    
    averagesmatrix = []
    averagesmatrix.append(NTerminusAverages[0])
    averagesmatrix.append(middleAverages[0])
    averagesmatrix.append(CTerminusAverages[0])
    
    shapematrix = [NTerminusAverages[1], middleAverages[1], CTerminusAverages[1]]
    
    shapesmatrix = [NTerminusAverages[2], middleAverages[2], CTerminusAverages[2]]
    
    shapes = [averagesmatrix, shapematrix, shapesmatrix]
    
    shapeCodes = []
    for row in shapes:
        if row[0] >= row[1] and row[2] >= row[1]:
            shapeCodes.append('Valley')
        if row[0] >= row[1] and row[2] < row[1]:
            shapeCodes.append('slopedown')
        if row[0] < row[1] and row[2] <= row[1]:
            shapeCodes.append('peak')
        if row[0] < row[1] and row[2] > row[1]:
            shapeCodes.append('slopeup')
            
            

    return shapeCodes



def writeFile(toWrite):
    
    newFileName = filename[:-4] + '-Output.csv'
    
    print 'Writing to file', newFileName
    
    wf=csv.writer(file(newFileName, 'wb'),dialect='excel')
    wf.writerows(toWrite)
    wf.writerow([])
    wf.writerow(['This program written by Dan Baum and Ben Luskin, Cavagnero Lab, University of Wisconsin-Madison, June 2011'])
    wf.writerow(['For equations used to calculate RCO and RCB, see Kurt, Mounce, Ellison, and Cavagnero, Biotechnology Progress, 24(3): 570-575 (2008)'])

def main():
    
    print sys.argv
    
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
    print "File name:", filename
    try:
        threshold = float(threshold)
    except ValueError:
        print "'" + sys.argv[2] +  "' is not a number."
        print "Please specify the threshold in Angstroms!"
        sys.exit(1)
    
    #atomicInfo contains atomic information
    print "Reading file..."
    atomicInfo, helixInfo, sheetInfo = readPDBFile(filename)
    print "Done reading file." 
    
    #Use atomicInfo to calculate distances to other atoms
    print "Calculating all distances..."
    distances = calculateDistances(atomicInfo)
    print "Done calculating distances."
    
    #Use distances to find the contacts
    print "Finding contacts..."
    atomContacts, residueContacts = findClosestContacts(distances, threshold)
    numContacts = countContacts(residueContacts)
    print "Done finding contacts."
    
    #Calculate residue-specific contact order and residue-specific contact breadth
    print "Calculating RCO and RCB..."
    RCO = calcRCO(numContacts)
    RCB = calcRCB(residueContacts)
    print "Done calculating RCO and RCB."
    
    #Calculate average RCO and RCB for 10%, 20%, and 30% termini.
    print "Calculating averages..."
    NTerminusRCO, middleRCO, CTerminusRCO = calcAverages(RCO)
    NTerminusRCB, middleRCB, CTerminusRCB = calcAverages(RCB)
    print "Done calculating averages."
    
    #Calculate percentage of each terminus that is a helix or sheet
    percentHelixN, percentHelixC = calcSecondary(helixInfo)
    percentSheetN, percentSheetC = calcSecondary(sheetInfo)
    
    #Calculate the shapes of each (peak, valley, etc)
    RCOshapes = shape(NTerminusRCO, middleRCO, CTerminusRCO)
    RCBshapes = shape(NTerminusRCB, middleRCB, CTerminusRCB)

    
    #Put everything into a matrix to print.
    NTerminusRCO.insert(0, 'N Terminus RCO')
    middleRCO.insert(0, 'Middle RCO')
    CTerminusRCO.insert(0, 'C Terminus RCO')
    NTerminusRCB.insert(0, 'N Terminus RCB')
    middleRCB.insert(0, 'Middle RCB')
    CTerminusRCB.insert(0, 'C Terminus RCB')
    RCOshapes.insert(0, '')
    RCBshapes.insert(0, '')
    printMatrix = [['% terminus','10%', '20%', '30%'], NTerminusRCO, middleRCO, CTerminusRCO, RCOshapes, NTerminusRCB, middleRCB, CTerminusRCB, RCBshapes]
    printMatrix.extend([[], ['# Residues', L], []])
    printMatrix.append(['With each terminus defined to be first or last 10%:'])
    printMatrix.append(['% Helix at N:', percentHelixN])
    printMatrix.append(['% Helix at C:', percentHelixC])
    printMatrix.append(['% Sheet at N:', percentSheetN])
    printMatrix.append(['% Sheet at C:', percentSheetC])
    printMatrix.append([])
    printMatrix.append([threshold, 'Your chosen contact threshold, in Angstroms'])
    
    #Write output to file.
    writeFile(printMatrix)
    
    print 'All done!'
    sys.exit(0)
    

if __name__ == '__main__':
  main()
