#!/usr/bin/python3

"""
// This is a helper QoMMMa Python file, which contains the functions used in execution of a Gaussian job. //
"""

# Global imports.
import os
import shutil

# Local imports.
from qomutil import qomend

def gauinp(imn, cwd, qmjob_prefix, qmkey, cha_mul, extra_basis, gau_head, extra_guess):
    """
    
    // Function which creates the input file required for a Gaussian job. //
    // This function is called in the function qmgaumain at each QoMMMa cycle. //
    
    Arguments
    ----------
    imn : integer
        If using the nudged elastic band, this is the image number.
    cwd : string
        The current working directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.
    qmkey : string
        Can either be None or the user input level of theory.
    cha_mul : list
        The charge and multiplicity given in the format of a list of two numbers (i.e., cha_mul = [0,1]).
    extra_basis : string
        User-specified basis set.
    gau_head : string
        Contains the details for assigning memory and number of processors for a QM job.
    extra_guess : string
        If the guess keyword is used, then it can be input.

    """
    
    # The Gaussian input file is created for the image.
    fd = open(('%s%d%s'%(qmjob_prefix, imn, '.in')), 'w') 

    # The variable gau_head is written to the input file, which contains details for assigning memory and processors.
    if gau_head.lower() != 'none':
        fd.write(gau_head.strip())  
        fd.write('\n')
    
    # The directory for the checkpoint file is written to the input file.
    cwd1 = cwd + '%s%d'%('/image', imn)
    chk = '%chk=' + cwd1 + ('%s%s%d%s'%('/', qmjob_prefix, imn, '.chk'))
    fd.write(chk)
    fd.write('\n')
    fd.write('\n')
    
    # Using a series of if and else statements, the qmkey (level of theory) is written to the input file.
    # In addition, the method of mixing the guess wavefunction is written along with some other important Gaussian parameters.
    if qmkey.lower() != 'none':
        if os.path.exists('%s%d%s'%(qmjob_prefix, imn, '.chk')):
            if extra_guess.lower() != 'none':
                jobkey = '# ' + qmkey + ' Nosymm Guess=(Read,' + extra_guess + ') Charge Force Prop=Grid'
            else:
                jobkey = '# ' + qmkey + ' Nosymm Guess=Read Charge Force Prop=Grid'
        else:
            if extra_guess.lower() != 'none':
                jobkey = '# ' + qmkey + 'Guess=' + extra_guess + ' Charge Nosymm Force Prop=Grid'
            else:
                jobkey = '# ' + qmkey + ' Nosymm Charge Force Prop=Grid'
    else:
        if os.path.exists('%s%d%s'%(qmjob_prefix, imn, '.chk')):
            jobkey = '# Guess=Read Nosymm Charge Force Prop=Grid'
        else:
            jobkey = '# Charge Nosymm Force Prop=Grid'
    fd.write(jobkey)
    fd.write('\n')
    fd.write('\n')
  
    # A placeholder title card is written to the input file.
    fd.write('TITLE')
    fd.write('\n')
    fd.write('\n')
    
    # The charge and multiplicity of the system is written to the input file.
    fd.write(str(cha_mul[0]).rjust(3))  
    fd.write(str(cha_mul[1]).rjust(5))  
    fd.write('\n')
 	  
    # The cartesian coordinates of the QM region are written to the input file.
    f = open(('%s%d%s'%('qmgeom', imn, '.xyz')), 'r') 
    for line in f:
        fd.write(line)
    f.close()
    fd.write('\n')
  
    # The point charge array is written to the input file.
    f = open(('%s%d%s'%('charges', imn, '.xyz')), 'r') 
    fresp = '%s%d%s'%('points', imn, 'gau.pts')
    numresp = len(open(fresp, 'r').readlines())
    for line in f:
        if line.strip() != '&pointch':
            if line.strip() != '&':
                line=line.split()
                fd.write(line[1].rjust(9))
                fd.write(line[2].rjust(13))
                fd.write(line[3].rjust(13))
                fd.write(line[0].rjust(10))
                fd.write('\n')
    f.close()
    fd.write('\n')
    
    # If a mixed bases set description is used, then it is written to the input file.
    if extra_basis.lower() != 'none':
        fd.write(extra_basis.lstrip())
        fd.write('\n')
    
    # The number of atoms with point charges as well as some numbers (??) are written to the input file.
    fd.write(str(numresp))
    fd.write(str(2).rjust(3))
    fd.write(str(20).rjust(4))
    fd.write(str(21).rjust(4))
    fd.write('\n')
    fd.write('\n')
    fd.close()
    shutil.copy(('%s%d%s'%('points', imn, 'gau.pts')), 'fort.20')

def gauout(fi, l, nqm, nlink, fortf):
    """
    
    // Function which extracts the energy, forces and Mulliken charges from a Gaussian output file. //
    
    Arguments
    ----------
    fi : string
        Gaussian job output file name.
    l : integer
        If using the nudged elastic band, this is the image number.
    nqm : integer
        Number of QM atoms.
    nlink : integer
        Number of link atoms.
    fortf : string
        Gaussian support file name containing the electric field.

    """
  
    # The output file containing QM energy and gradients is created.
    fd = open(('%s%d'%('ab_initio' ,l)), 'a')
    fd.write('Energy')
    fd.write('\n')
    
    # The output file containing QM Mulliken charges is created.
    fdm = open(('%s%d'%('mulliken', l)), 'a')
    
    # The Gaussian output file is opened.
    f = open(fi, 'r')
    
    # Total atoms are initialised.
    ntot = nqm + nlink
    
    # While loop iterates through the whole Gaussian output file.
    while 1:
        line = f.readline()
        if not line:break
        
        # The self energy of the charges is saved which represents the electrostatic interaction between point charges.
        if line[:29] == ' Self energy of the charges =':
            pointE = float(line.split()[6])
        
        # If the convergence criteria have not been met for the self consistent field (SCF) iterations, then the lines reporting this are skipped.
        elif line[:42] == ' >>>>>>>>>> Convergence criterion not met.':
            line = f.readline()
            line = f.readline()
        
        # Once the SCF is done, the energy is extracted, and the true energy is calculated by subtracting the point charge energy from the SCF energy.
        # This energy is written to the file ab_initio*.
        elif line[:9] == ' SCF Done':
            totE = float(line.split()[4])
            ene = totE - pointE 
            fd.write(str(ene))
            fd.write('\n')
        
        # The Mulliken charges are extracted from the Gaussian output file.
        # These Mulliken charges and the total charge are written to the file mulliken*.
        # Depending on the version of Gaussian, the Mulliken charges data is labelled different so two cases must be implemented.
        elif (line[:18] == ' Mulliken charges:') or (line[:37] == ' Mulliken charges and spin densities:'):
            if (line[:18] == ' Mulliken charges:'): temp = 4
            if (line[:37] == ' Mulliken charges and spin densities:'): temp = 5
            f.readline()
            for i in range(ntot):
                line = f.readline().split()
                fdm.write(line[2].rjust(10))
                fdm.write('\n')  
            line = f.readline().split()
            fdm.write('total charge')
            fdm.write('\n')
            fdm.write(line[temp].rjust(10)) 
            fdm.write('\n')
            fdm.close()

        # The gradients corresponding to each cartesian coordinate are extracted from the Gaussian output file.
        # This gradient is written to the file ab_initio*.
        elif line[:59] == ' Center     Atomic                   Forces (Hartrees/Bohr)': 
            fd.write('Gradient')
            fd.write('\n')
            f.readline()
            f.readline()
            for i in range(ntot):
                line = f.readline().split()
                fd.write(line[0].rjust(3))
                linex = float(line[2]) * -1
                fd.write(str(linex).rjust(20))
                liney = float(line[3]) * -1
                fd.write(str(liney).rjust(20))
                linez = float(line[4]) * -1
                fd.write(str(linez).rjust(20))
                fd.write('\n')  
            fd.close()
    f.close()
    
    # The point charge array is opened.
    ftc = open('%s%d%s'%('points_chg', l, '.pts'), 'r')
    
    # The output file containing gradients corressponding to the MM point charges is created.
    fg = open('%s%d%s'%('qmlatgrad', l, '.out'), 'a')
    
    # The electric field generated by Gaussian is opened.
    ff = open(fortf, 'r') 
    
    # Every point charge is added to the list cha.
    cha = []
    for line in ftc:
        tcha = float(line.split()[0])
        cha.append(tcha)
    ftc.close()
    
    # The number of point charges is written to qmlatgrad*.out.
    fg.write(' gradient output')
    fg.write('\n')
    fg.write(str(len(cha)))
    fg.write('\n')
    
    # Each gradient for MM point charges is calculated by multiplication of each point charge by the electric field.
    # This gradient is then written to qmlatgrad*.out.
    for i in range(len(cha)):
        ff.readline()
        linef = ff.readline().split()               
        for j in range(3):
            grad =- (float(cha[i]) * float(linef[j]))
            fg.write(str(grad).rjust(25))
        fg.write('\n')
    fg.close()

def qmgaumain(imn, cwd, usrdir, qmjob_prefix, qmgau_job, nqm, nlink, cln, qmkey, cha_mul, extra_basis, gau_head, extra_guess):
    """
    
    // Function which combines the functions in this file for execution of a Gaussian job. //
    // It will prepare the input file, run the Gaussian job, and then extract the results from the output. //
    
    Arguments
    ----------
    imn : integer
        If using the nudged elastic band, this is the image number.
    cwd : string
        The current working directory.
    usrdir : string
        The user directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.
    qmgau_job : string
        The path to the Gaussian program.
    nqm : integer
        Number of QM atoms.
    nlink : integer
        Number of link atoms.
    cln : integer
        The QoMMMa cycle number.
    qmkey : string
        Can either be None or the user input level of theory.
    cha_mul : list
        The charge and multiplicity given in the format of a list of two numbers (i.e., cha_mul = [0,1]).
    extra_basis : string
        User-specified basis set.
    gau_head : string
        Contains the details for assigning memory and number of processors for a QM job.
    extra_guess : string
        If the guess keyword is used, then it can be input.

    """
    
    # The files from the previous Gaussian run are rearranged and renamed with the prefix 'old_'.
    if os.path.exists(cwd + ('%s%d'%('/image', imn)) + ('%s%s%d%s'%('/', qmjob_prefix, imn, '.in'))):
        try:
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '.in')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '.in'))) 
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '.log')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '.log')))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '.in'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '.log'))
        except:
            qomend('ERROR, while rearranging previous cycle Gaussian files for image : ' + str(imn), cwd, usrdir)
   
    # The Gaussian input file is created.
    try:
        gauinp(imn, cwd, qmjob_prefix, qmkey, cha_mul, extra_basis, gau_head, extra_guess)	
    except:
        qomend('Error while creating Gaussian input for image : ' + str(imn), cwd, usrdir)
    
    # The Gaussian job is performed.    
    try:
        os.system(qmgau_job + ' ' + ('%s%d%s'%(qmjob_prefix, imn, '.in')))
    except:
        qomend('Error while running Gaussian job at cycle: ' + str(cln) + ' image :' + str(imn), cwd, usrdir)
    
    # The results from the Gaussian job are extracted and rearranged.
    try:
        fi = '%s%d%s'%(qmjob_prefix, imn, '.log')
        gauout(fi, imn, nqm, nlink, 'fort.21')	
        shutil.copy(('%s%d'%('ab_initio', imn)), cwd)
        shutil.copy(('%s%d'%('mulliken', imn)), cwd)
        shutil.copy(('%s%d%s'%('qmlatgrad', imn, '.out')), cwd) 
        os.remove('%s%d%s'%('qmlatgrad', imn, '.out'))
        os.remove('%s%d'%('ab_initio', imn))
        os.remove('%s%d'%('mulliken', imn))
    except:
        qomend('ERROR, while reading Energy or Gradient or Mulliken charge or electrostatic potential from QM output file of image : ' + str(imn), cwd, usrdir)    

def mgauinp(imn, cwd, qmjob_prefix, qmkey, cha_mul1, cha_mul2, extra_basis, gau_head, extra_guess):
    """
    
    // Function which creates the input files required for the Gaussian jobs for MECP calculations. //
    // This function is called in the function mecp_gaumain at each QoMMMa cycle. //
    
    Arguments
    ----------
    imn : integer
        If using the nudged elastic band, this is the image number.
    cwd : string
        The current working directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.
    qmkey : string
        Can either be None or the user input level of theory.
    cha_mul1 : list
        The charge and multiplicity for state A given in the format of a list of two numbers (i.e., cha_mul = [0,1]).
    cha_mul2 : list
        The charge and multiplicity for state B given in the format of a list of two numbers (i.e., cha_mul = [0,1]).
    extra_basis : string
        User-specified basis set.
    gau_head : string
        Contains the details for assigning memory and number of processors for a QM job.
    extra_guess : string
        If the guess keyword is used, then it can be input.

    """
    
    # The Gaussian input file is created for the image.
    # This operation is performed for both states A and B.
    fda = open(('%s%d%s'%(qmjob_prefix, imn, '_A.in')), 'w')  
    fdb = open(('%s%d%s'%(qmjob_prefix, imn, '_B.in')), 'w') 
    
    # The variable gau_head is written to the input file, which contains details for assigning memory and processors.
    # This operation is performed for both states A and B.
    if gau_head.lower() != 'none':
        fda.write(gau_head.strip())  
        fda.write('\n')
        fdb.write(gau_head.strip())  
        fdb.write('\n')
    
    # The directory for the checkpoint file is written to the input file.
    # This operation is performed for both states A and B.
    cwd1 = cwd + '%s%d'%('/image', imn)
    chk = '%chk=' + cwd1 + ('%s%s%d%s'%('/', qmjob_prefix, imn, '_A.chk'))
    fda.write(chk)
    chk = '%chk=' + cwd1 + ('%s%s%d%s'%('/', qmjob_prefix, imn, '_B.chk'))
    fdb.write(chk)
    fda.write('\n')
    fda.write('\n')
    fdb.write('\n')
    fdb.write('\n')
    
    # Using a series of if and else statements, the qmkey (level of theory) is written to the input file.
    # In addition, the method of mixing the guess wavefunction is written along with some other important Gaussian parameters.  
    # These operations are performed for both states A and B.
    if qmkey.lower() != 'none':
        if os.path.exists('%s%d%s'%(qmjob_prefix, imn, '_A.chk')):
            if extra_guess.lower() != 'none':
                jobkey = '# ' + qmkey + ' Nosymm Guess=(Read,' + extra_guess + ') Charge Force Prop=Grid'
            else:
                jobkey = '# ' + qmkey + ' Nosymm Guess=Read Charge Force Prop=Grid'
        else:
            jobkey = '# ' + qmkey + ' Charge Nosymm Force Prop=Grid'
    else:
        if os.path.exists('%s%d%s'%(qmjob_prefix, imn, '_A.chk')):
            jobkey = '# Guess=Read Nosymm Charge Force Prop=Grid'
        else:
            jobkey = '# Charge Nosymm Force Prop=Grid'
    fda.write(jobkey)
    fda.write('\n')
    fda.write('\n')
    if qmkey.lower() != 'none':
        if os.path.exists('%s%d%s'%(qmjob_prefix, imn, '_B.chk')):
            if extra_guess.lower() != 'none':
                jobkey = '# ' + qmkey + ' Nosymm Guess=(Read,' + extra_guess + ') Charge Force Prop=Grid'
            else:
                jobkey = '# ' + qmkey + ' Nosymm Guess=Read Charge Force Prop=Grid'
        else:
            jobkey = '# ' + qmkey + ' Charge Nosymm Force Prop=Grid'
    else:
        if os.path.exists('%s%d%s'%(qmjob_prefix, imn, '_B.chk')):
            jobkey = '# Guess=Read Nosymm Charge Force Prop=Grid'
        else:
            jobkey = '# Charge Nosymm Force Prop=Grid'
    fdb.write(jobkey)
    fdb.write('\n')
    fdb.write('\n')
    
    # A placeholder title card is written to the input file.
    # This operation is performed for both states A and B.
    fda.write('TITLE')
    fda.write('\n')
    fda.write('\n')
    fdb.write('TITLE')
    fdb.write('\n')
    fdb.write('\n')
    
    # The charge and multiplicity of the system is written to the input file.  
    # This operation is performed for both states A and B.
    fda.write(str(cha_mul1[0]).rjust(3))  
    fda.write(str(cha_mul1[1]).rjust(5))  
    fda.write('\n') 	  
    fdb.write(str(cha_mul2[0]).rjust(3))  
    fdb.write(str(cha_mul2[1]).rjust(5))  
    fdb.write('\n') 
    
	# The cartesian coordinates of the QM region are written to the input file.
    # This operation is performed for both states A and B.
    f = open(('%s%d%s'%('qmgeom', imn, '.xyz')), 'r') 
    for line in f:
        fda.write(line)
        fdb.write(line)
    f.close()
    fda.write('\n')
    fdb.write('\n')
  
    # The point charge array is written to the input file.
    # This operation is performed for both states A and B.
    f = open(('%s%d%s'%('charges', imn, '.xyz')), 'r') 
    ff1 = open('fort.20', 'w')
    ff2 = open('fort.30', 'w')
    ill = 0
    for line in f:
        if line.strip() != '&pointch':
            if line.strip() !='&':
                line = line.split()
                fda.write(line[1].rjust(13))
                ff1.write(line[1].rjust(20))
                fda.write(line[2].rjust(13))
                ff1.write(line[2].rjust(20))
                fda.write(line[3].rjust(13))
                ff1.write(line[3].rjust(20))
                fda.write(line[0].rjust(10))
                fda.write('\n')
                ff1.write('\n')
                fdb.write(line[1].rjust(13))
                ff2.write(line[1].rjust(20))
                fdb.write(line[2].rjust(13))
                ff2.write(line[2].rjust(20))
                fdb.write(line[3].rjust(13))
                ff2.write(line[3].rjust(20))
                fdb.write(line[0].rjust(10))
                fdb.write('\n')
                ff2.write('\n')
                ill += 1
    f.close()
    ff1.close()
    ff2.close()
    fda.write('\n')
    fdb.write('\n')
    
    # If a mixed bases set description is used, then it is written to the input file.
    # This operation is performed for both states A and B.
    if extra_basis.lower() != 'none':
        fda.write(extra_basis.lstrip())
        fda.write('\n')
        fdb.write(extra_basis.lstrip())
        fdb.write('\n')
    
    # The number of atoms with point charges as well as some numbers (??) are written to the input file.
    # These operations are performed for both states A and B.
    fda.write(str(ill))
    fda.write(str(2).rjust(3))
    fda.write(str(20).rjust(4))
    fda.write(str(21).rjust(4))
    fda.write('\n')
    fdb.write(str(ill))
    fdb.write(str(2).rjust(3))
    fdb.write(str(30).rjust(4))
    fdb.write(str(31).rjust(4))
    fdb.write('\n')
    fda.close()
    fdb.close()

def mecp_gaumain(imn, cwd, usrdir, qmjob_prefix, qmgau_job, nqm, nlink, cln, qmkey, cha_mul1, cha_mul2, extra_basis, gau_head, extra_guess):
    """
    
    // Function which combines the functions in this file for execution of Gaussian jobs for MECP calculations. //
    // It will prepare the input files, run the Gaussian jobs, and then extract the results from the outputs. //
    
    Arguments
    ----------
    imn : integer
        If using the nudged elastic band, this is the image number.
    cwd : string
        The current working directory.
    usrdir : string
        The user directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.
    qmgau_job : string
        The path to the Gaussian program.
    nqm : integer
        Number of QM atoms.
    nlink : integer
        Number of link atoms.
    cln : integer
        The QoMMMa cycle number.
    qmkey : string
        Can either be None or the user input level of theory.
    cha_mul1 : list
        The charge and multiplicity for state A given in the format of a list of two numbers (i.e., cha_mul = [0,1]).
    cha_mul2 : list
        The charge and multiplicity for state B given in the format of a list of two numbers (i.e., cha_mul = [0,1]).
    extra_basis : string
        User-specified basis set.
    gau_head : string
        Contains the details for assigning memory and number of processors for a QM job.
    extra_guess : string
        If the guess keyword is used, then it can be input.

    """
    
    # The files from the previous Gaussian runs are rearranged and renamed with the prefix 'old_'.
    # This operation is performed for both states A and B.
    if os.path.exists(cwd + ('%s%d'%('/image', imn)) + ('%s%s%d%s'%('/', qmjob_prefix, imn, '_A.in'))):
        try:
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_A.in')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_A.in'))) 
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_A.log')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_A.log')))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_A.in'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_A.log'))
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_B.in')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_B.in'))) 
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_B.log')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_B.log')))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_B.in'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_B.log'))
        except:
            qomend('error while rearrangeing Gaussian files created in previous run at cycle :' + str(cln) + ' for image: ' + str(imn), cwd, usrdir)
            
    # The Gaussian input files are created for both states A and B.        
    try:
        mgauinp(imn, cwd, qmjob_prefix, qmkey, cha_mul1, cha_mul2, extra_basis, gau_head, extra_guess)
    except:
        qomend('Error while trying to prepare Gaussian input files at cycle : ' + str(cln) + ' for image : ' + str(imn), cwd, usrdir)
        
    # The Gaussian job for state A is performed.       
    try:
        os.system(qmgau_job + ' ' + ('%s%d%s'%(qmjob_prefix, imn, '_A.in')))
    except:
        qomend('Error while running Gaussian job for state A at cycle : ' + str(cln) + ' for image : ' + str(imn), cwd, usrdir)
        
    # The Gaussian job for state B is performed.
    try:
        os.system(qmgau_job + ' ' + ('%s%d%s'%(qmjob_prefix, imn, '_B.in')))
    except:
        qomend('Error while running Gaussian job for state B at cycle : ' + str(cln) + ' for image : ' + str(imn), cwd, usrdir)
        
    # The results from the Gaussian jobs are extracted and rearranged.   
    # These operations are performed for both states A and B.
    try:
        fi = '%s%d%s'%(qmjob_prefix, imn, '_A.log')
        gauout(fi, imn, nqm, nlink, 'fort.21')	
        fi = '%s%d%s'%(qmjob_prefix, imn, '_B.log')
        gauout(fi, imn, nqm, nlink, 'fort.31')	
        shutil.copy(('%s%d'%('ab_initio', imn)), cwd)
        shutil.copy(('%s%d'%('mulliken', imn)), cwd)
        shutil.copy(('%s%d%s'%('qmlatgrad', imn, '.out')), cwd) 
        os.remove('%s%d%s'%('qmlatgrad', imn, '.out'))
        os.remove('%s%d'%('ab_initio', imn))
        os.remove('%s%d'%('mulliken', imn))
    except:
        qomend('Error while extracting results from Gaussian output files at cycle : ' + str(cln) + ' for image : ' + str(imn), cwd, usrdir)