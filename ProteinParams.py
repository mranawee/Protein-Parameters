# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 20:12:13 2020

@author: mano_
"""

#!/usr/bin/env python3
# Name: Mano (mranawee)
# Group Members: Harrison Wismer (hwsismer), Chianna Li (chanli)

class ProteinParam :
    '''
    These tables will be used in methods for calculating:
     - molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
     - absorbance at 280 nm (aa2abs280)
     - pKa of positively charged Amino Acids (aa2chargePos)
     - pKa of negatively charged Amino acids (aa2chargeNeg)
     - constants aaNterm and aaCterm for pKa of the respective termini
     
     Example)
     Input(protein string): VLSPADKTNVKAAW
     
     Expected Output: Number of Amino Acids: 14
                      Molecular Weight: 1499.7
                      molar Extinction coefficient: 5500.00
                      mass Extinction coefficient: 3.67
                      Theoretical pI: 9.88
                      Amino acid composition:
                      A = 21.43%
                      C = 0.00%
                      D = 7.14%
                      E = 0.00%
                      F = 0.00%
                      G = 0.00%
                      H = 0.00%
                      I = 0.00%
                      K = 14.29%
                      L = 7.14%
                      M = 0.00%
                      N = 7.14%
                      P = 7.14%
                      Q = 0.00%
                      R = 0.00%
                      S = 7.14%
                      T = 7.14%
                      V = 14.29%
                      W = 7.14%
                      Y = 0.00%

    '''

 
    aa2mw = { #Molecular weight of each AA
             'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093,
            'C': 121.158,
             'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103,
            'I': 131.173,
             'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188,
            'Q': 146.145,
             'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201,
            'Y': 181.189
             }
 
    mwH2O = 18.015 
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}
    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34
 
    def __init__ (self, protein):
        '''
        This method initialized the protein attribute as well as the Amino Acid Composition dictionary that will be used for
        calculating the AA composition in the aaComposition method.
        '''
        self.protein = protein.upper() #Making the inputted protein string upper case
        
        self.aaCompDict = {key:0 for key in ProteinParam.aa2mw.keys()} #Taking keys from the aa2mw dictionary
        
        '''
        This for loop will iterate through the inputted protein string and add to the composition dictionary the number of each
        AA.  The main method will first call the aaComp method, which will just return the dictionary with the counts of each
        AA.  The main method will then do the appropriate calculations for the AA composition.
        '''
        for aa in self.protein:
            if aa in self.aaCompDict:
                self.aaCompDict[aa] +=1
    
    def aaCount (self):
        '''
        This method will count the number of amino acids from the inputted string. What's very important here is that
        this counts only valid amino acids from the dictionary already made and no spaces.  To make sure this happens,
        I added an if statement within the for loop that only allows counting if the AA is already found in the dictionary.
        '''
        count = 0
        for aa in self.protein:
            if aa in ProteinParam.aa2mw.keys(): #Count won't include invalid characters or spaces
                count += 1
        
        return (count)
 
    def pI (self): 
        '''
        This method calculates the theoretical pI by finding the pH that yields the closest net charge to 0.  This 
        for loop iterates over all pH values to find the charge closest to 0.  To find the best pH accurate to 2 decimal
        places, I made a range of 1401 followed by dividing the pH value by 100.  The value closest to 0 as a pH was the 
        best value.  Each pH that is iterated then calls the charge method.  
        '''
        bestCharge = 100000
        #for pH in range(14+1):
        for pH in range(1401):
            pH = pH / 100 #Makes the calculation more precise by increments of 0.01(Help by Dennis Mulligan)
            thisCharge = self._charge_(pH)
            if abs(thisCharge) < abs(bestCharge):
                bestCharge = thisCharge
                bestPH = pH
        
        return(bestPH)
           
    def aaComposition (self) :
        '''
        The dictionary is already made for calculating AA composition.  The for loop is already made in the init
        method for iteration and the proper calculations are done in the main method.  
        '''
        return (self.aaCompDict)
        
    def _charge_ (self,pH):
        '''
        This method is never used directly by the main method.  I did a for loop for positive and negative charges,
        followed by an if statement for each separately.  I then wrote equations for calculating the charges
        of the termini.
        '''
        pos = 0
        neg = 0
        for aa in self.protein:
            
            if aa in ProteinParam.aa2chargePos.keys():
                posTerm1 = float(10**ProteinParam.aa2chargePos.get(aa)) #numerator
                posTerm2 = float(posTerm1 + 10**pH) #denominator
                pos += float(posTerm1 / posTerm2)
            
            
            elif aa in ProteinParam.aa2chargeNeg.keys():
                negTerm1 = float(10**pH) #numerator
                negTerm2 = float(10**ProteinParam.aa2chargeNeg.get(aa) + 10**pH) #denominator
                neg += float(negTerm1 / negTerm2)
            
            else:
                netCharge = 0
               
        NTerminusTerm1 = float(10**ProteinParam.aaNterm)
        NTerminusTerm2 = float((10**ProteinParam.aaNterm) + 10**pH)
        NTerm = NTerminusTerm1 / NTerminusTerm2
        
                              
        CTerminusTerm1 = float((10**pH)) 
        CTerminusTerm2 = float(ProteinParam.aaCterm + 10**pH)
        CTerm = CTerminusTerm1 / CTerminusTerm2
        
        
        netCharge = float((pos + NTerm) - (neg + CTerm))
        
        return(netCharge)
        
    def molarExtinction (self):
        '''
        From the AA composition dictionary, I can use that to calculate the molar extinction coefficient of the 
        inputted protein.  From that dictionary, we will know the composition of Ys, Ws, and Cs.  The coefficients
        for each of those are already provided, and are calculated.  
        
        I first set this method up by taking the count of each of the three AAs in the protein.  
        
        In the next three lines, I multiplied the counts by each molar coefficient and added it all
        to get the total.
        '''
        Y = self.aaCompDict['Y'] #Count of each of the AAs
        W = self.aaCompDict['W']
        C = self.aaCompDict['C']
        
        molarY = Y * self.aa2abs280['Y']#Count multiplied by the molar coefficient of each
        molarW = W * self.aa2abs280['W']
        molarC = C * self.aa2abs280['C']
        
        totalMolar = molarY + molarW + molarC #Total molar extinction coefficient calculated
        return(totalMolar)

    def massExtinction (self):
        '''
        Mass extinction coefficient is calculated by dividing the molar extinction coefficient calculated by the 
        molecular weight. Molecular weight and molar extinction is called in this method.
        '''
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0
 
    def molecularWeight (self):
        '''
        This method used to calculate the molecular weight included a for loop that started by iterating through the 
        inputted string. If each amino acid inputted was found in the aa2mw dictionary, the value of the AA's mw was calculated
        and added to the next AA mw found.  The molecular weight of water was then subtracted.
        '''
        mw = 0
        for aa in self.protein:
            if aa in ProteinParam.aa2mw.keys():
                mw += ProteinParam.aa2mw.get(aa) - ProteinParam.mwH2O 
        newMW = sum([mw]) + ProteinParam.mwH2O
        return(newMW)
        

# Please do not modify any of the following. This will produce a standard output that can be parsed

import sys
def main():
 
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))

        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
 
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
 
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
 
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))

        print ("Amino acid composition:")
 
        myAAcomposition = myParamMaker.aaComposition()
 
        keys = list(myAAcomposition.keys())
 
        keys.sort()
 
        if myAAnumber == 0 : myAAnumber = 1 # handles the case where no AA are present

        for key in keys :

            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))
            
        inString = input('protein sequence?')

if __name__ == "__main__":
 
    main()