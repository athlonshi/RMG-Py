#!/usr/bin/python
# -*- coding: utf-8 -*-

################################################################################
#
#	RMG - Reaction Mechanism Generator
#
#	Copyright (c) 2002-2009 Prof. William H. Green (whgreen@mit.edu) and the
#	RMG Team (rmg_dev@mit.edu)
#
#	Permission is hereby granted, free of charge, to any person obtaining a
#	copy of this software and associated documentation files (the 'Software'),
#	to deal in the Software without restriction, including without limitation
#	the rights to use, copy, modify, merge, publish, distribute, sublicense,
#	and/or sell copies of the Software, and to permit persons to whom the
#	Software is furnished to do so, subject to the following conditions:
#
#	The above copyright notice and this permission notice shall be included in
#	all copies or substantial portions of the Software.
#
#	THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#	FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#	DEALINGS IN THE SOFTWARE.
#
################################################################################

"""
A module for generating chemkin-compatible mechanism file (chem.inp) and
thermodynamic file (therm.dat)
"""

import quantities as pq
import log as logging
import os
import math

import constants
import settings
import model
import structure
import data
import species
import reaction
import thermo
from time import localtime, strftime

class InvalidFormatException(Exception):
	"""
	An exception used when outputing Chemkin-compatible files if RMG input
        cannot satisfy Chemkin requirement. The `msg` parameter is used to
        specify what about the file caused the exception to be raised.
	"""

	def __init__(self, msg):
		self.msg = msg

	def __str__(self):
		return 'Invalid RMG input to Chemkin-compatible format requirement : ' + self.msg


def getChemkinName(speciesFormula, speciesID, radicalCount):
        """
        Return a chemkin-compatible species label, maximum label lenth is 10
        """
        if (radicalCount == 0):
            speciesFormula +='M'
        elif (radicalCount == 1):
            speciesFormula += 'J'
        elif (radicalCount == 2):
            speciesFormula += 'JJ'
        elif (radicalCount == 3):
            speciesFormula += 'JJJ'
            
        ChemkinName = speciesFormula + str(speciesID)

        if (len(ChemkinName) > 10):
            ChemkinName = 'SPC' + str(speciesID)
        return ChemkinName
    
def writeMech(fstr, reactionModel, reactionSystems):
	"""
	Write a chemkin-compatible mechanism file chem.inp at the location
        `fstr` based on the provided `reactionModel` and `reactionSystems`.
	"""
        path = os.path.join(fstr, 'chem.inp')
        try:
            chem_file = open(path,"w")
        except IOError, error:
            print "Issue occurred when creating mechanism file '%s' with error '%s'" %(path,error)
        else:
            chem_file.write("!This is an automatically generated mechanism using RMG\n")
            chem_file.write("!Created on\t"+strftime("%a, %d %b %Y %H:%M", localtime())+"\n")
#Create elements list
        chem_file.write("ELEMENTS\n")
        elementcount = 0
        element_list = ''
        for spec in reactionModel.core.species:
            for element in spec.getFormula():
                if not element.isdigit() and not element in element_list:
                    elementcount += 1
                    element_list += element + ' '
        chem_file.write(element_list)
        chem_file.write("\n!Total elements:"+str(elementcount)+"\nEND\n")

#Create species list
        chem_file.write("SPECIES\n")
        speciescount = 0
        for spec in reactionModel.core.species:
            ChemkinName = getChemkinName(spec.getFormula(),spec.id, spec.structure[0].getRadicalCount())
            chem_file.write(ChemkinName.upper() +"\t")
            speciescount += 1
            if (speciescount%5 == 0):chem_file.write("\n")

        if (speciescount%5 == 0):
            chem_file.write("!Total species:"+str(speciescount)+"\nEND\n")
        else:
            chem_file.write("\n!Total species:"+str(speciescount)+"\nEND\n")

#Create reactions list
        chem_file.write("REACTIONS KJOULES/MOLE MOLECULES\n")
        reactioncount = 0
        for rxn in reactionModel.core.reactions:
            reactioncount += 1
            string =''
            ChemkinReac = ''
            ChemkinProd = ''
            for reac in rxn.reactants:
                ChemkinReac = getChemkinName(reac.getFormula(),reac.id,reac.structure[0].getRadicalCount())
                string+= ChemkinReac + '+'
            if isinstance(rxn, reaction.PDepReaction):
                string = string[0:-1]+'(+M)<=>'
            else:
                string = string[0:-1] + '<=>'

            for prod in rxn.products:
                ChemkinProd = getChemkinName(prod.getFormula(),prod.id,prod.structure[0].getRadicalCount())
                string += ChemkinProd + '+'
            if isinstance(rxn, reaction.PDepReaction):
                string = string[0:-1]+'(+M)'
            else:
                string = string[0:-1]

            if isinstance(rxn, reaction.PDepReaction):
		string += "\t0.0\t0.0\t0.0 !Using Chebyshev coefficient, the Arrhenius expression is replaced"
                for T2 in range(rxn.kinetics.degreeT):
                    if (T2 == 0):
                        Pmin1 = rxn.kinetics.Pmin * pq.Pa
                        Pmin1.units = pq.atm
                        Pmax1 = rxn.kinetics.Pmax * pq.Pa
                        Pmax1.units = pq.atm
                        string += "\nTCHEB/%8.3f\t%8.3f/" %(rxn.kinetics.Tmin, rxn.kinetics.Tmax)
                        string += "\nPCHEB/%8.3f\t%8.3f/" %(Pmin1, Pmax1)
                        unitConvert = math.log10(100**((len(rxn.reactants)-1)*3)) #According to Cheb
                        string += "\nCHEB/%d\t%d\t%s/" % (rxn.kinetics.degreeT,rxn.kinetics.degreeP, \
                        str(rxn.kinetics.coeffs[0,0]+unitConvert)+"\t"+str(rxn.kinetics.coeffs[T2,1:rxn.kinetics.degreeP])[1:-1])
                    else:
                        string += "\nCHEB/%s/" % (str(rxn.kinetics.coeffs[T2,:])[1:-1])
            else:
                kcount = 0
                #Note chemkin use cm instead of m in RMG, so the unit of reaction rate has to be converted
                for k in rxn.kinetics:
                    kcount += 1
                    if kcount == 1:
                        string += "\t%12.3e\t%4.2f\t%12.3e\t" %(k.A*(100**((len(rxn.reactants)-1)*3)), k.n, -(k.E0+k.alpha*rxn.getEnthalpyOfReaction(298.15)))
                    else:
                        string +="\n!Alternative reaction constant set %d:" % (kcount-1) + "%12.3e\t%4.2f\t%12.3e\t" %(k.A*(100**((len(rxn.reactants)-1)*3)), k.n, -(k.E0+k.alpha*rxn.getEnthalpyOfReaction(298.15)))
            #Check duplicated reaction and flag it in chem.inp
            if (reactionModel.core.reactions.count(rxn)>1):
                string += "\nDUP"
            chem_file.write(string.upper()+"\n")
        chem_file.write("!Total reactions:"+str(reactioncount)+"\nEND\n")

        chem_file.close()
	# Print to log

	logging.info('')
	logging.info('chem.inp written to ' + fstr)

def writeTherm(fstr, reactionModel, reactionSystems):
        """
        Write a chemkin-compatible thermodynamical data file therm.dat at the
        location `fstr` based on the provided `reactionModel` and
        `reactionSystems`.
        """
        path = os.path.join(fstr, 'therm.dat')
        try:
            therm_file = open(path,"w")
        except IOError, error:
            print "Issue occurred when creating thermal file '%s' with error '%s'" %(path,error)
        else:
            therm_file.write("!This is an automatically generated thermal data file using RMG\n")
            therm_file.write("!Created on\t"+strftime("%a, %d %b %Y %H:%M", localtime())+"\n")


        therm_file.write("THERMO\n")
        therm_file.write("%10.3f%10.3f%10.3f" %(300.0, 1000.0, 5000.0))

        for spec in reactionModel.core.species:
            string = "\n"
            ChemkinName = getChemkinName(spec.getFormula(),spec.id, spec.structure[0].getRadicalCount())
            string +="%-16s" %(ChemkinName)+'        '  # 8 spaces
            elementcount = 0
            element_list = ''
            for atom in spec.structure[0].atoms():
                element = atom.atomType.element.symbol
                if not element in element_list:
                    element_list += element
                    elementcount += 1
                    string += "%-2s%-3d" % (element, spec.structure[0].countNumberOfAtomElement(atom,element))
            if (elementcount == 4):
                pass
            elif elementcount < 4:
                for i in range(4-elementcount):
                    string += '     '   #5 spaces
            else:
                raise InvalidFormatException("%d atoms found in species %s, but currently thermo data file allows for 4" % (elementcount, ChemkinName))

            string_thermo = spec.thermoData.toCHEMKIN()
            string += string_thermo
            therm_file.write(string.upper())

        therm_file.close()
        # Print to log

        logging.info('')
	logging.info('thermo.dat written to ' + fstr)