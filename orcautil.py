#!/usr/bin/python3

"""
// This is a helper QoMMMa Python file, which contains the functions used in execution of an ORCA job. //
"""

# Global imports.
import os
import shutil

# Local imports.
from qomutil import qomend

def orcainp(imn, cwd, qmjob_prefix, qmkey, cha_mul, extra_basis, orca_head):
    """
    
    // Function which creates the input file required for an ORCA job. //
    // This function is called in the function qmorcamain at each QoMMMa cycle. //
    
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
    orca_head : string
        Contains the details for assigning memory and number of processors for a QM job.

    """
  
    # The ORCA input file is created for the image.
    fd = open(('%s%d%s'%(qmjob_prefix, imn, '.in')), 'w')
    
    # The variable orca_head is written to the input file, which contains details for assigning memory and processors.
    if orca_head.lower() != 'none':
        fd.write(orca_head.strip())  
        fd.write('\n')
        
    # Using a series of if and else statements, the qmkey (level of theory) is written to the input file.
    # If a .gbw wavefunction file is found, then this is used as a guess.
    cwd1 = cwd + '%s%d'%('/image', imn)
    if qmkey.lower() != 'none':
        if os.path.exists('%s%s%d%s'%('old_', qmjob_prefix, imn, '.gbw')):
            jobkey='!' + qmkey + ' moread engrad'
        else:
            jobkey='!' + qmkey + ' engrad'
    else:
        if os.path.exists('%s%s%d%s'%('old_', qmjob_prefix,imn,'.gbw')):
            jobkey = '! moread engrad'
        else:
            jobkey = '! engrad'
    fd.write(jobkey)
    fd.write('\n')
  
    # The point charge array is written to the input file.
    ptc = '%pointcharges "' + '%s%d%s'%('ptchg', imn, '.pc"')
    fd.write(ptc)
    fd.write('\n')
    
    # The guess wavefunction is written to the input file.
    if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'.gbw')):    
        wavefunction='%moinp '+'"'+'%s%s%d%s'%('old_',qmjob_prefix,imn,'.gbw'+'"')
        fd.write(wavefunction)
        fd.write('\n')
        
    # If a mixed bases set description is used, then it is written to the input file.     
    fd.write(extra_basis.strip())
    fd.write('\n')
    fd.write('\n')
    fd.write('* xyz ')
    
    # The charge and multiplicity of the system is written to the input file.
    fd.write(str(cha_mul[0]).rjust(3))  
    fd.write(str(cha_mul[1]).rjust(5))  
    fd.write('\n') 	 
    
    # The cartesian coordinates of the QM region are written to the input file.
    f = open(('%s%d%s'%('qmgeom', imn, '.xyz')), 'r') 
    for line in f:
        fd.write(line)
    f.close()
    fd.write('*')
    fd.write('\n')
    fd.write('\n')
    fd.write('\n')
    fd.write('\n')
    fd.write('\n')
    fd.close()
    
    # The point charge array is written to the file ptchg*.pc.
    fg = open(('%s%d%s'%('ptchg', imn, '.pc')), 'a')
    f = open(('%s%d%s'%('charges', imn, '.xyz')), 'r') 
    fresp = '%s%d%s'%('points', imn, 'gau.pts')
    numresp = len(open(fresp, 'r').readlines())
    fg.write(str(numresp))
    fg.write('\n')
    for line in f:
        if line.strip() != '&pointch':
            if line.strip() != '&':
                line = line.split()
                fg.write(line[0].rjust(10))
                fg.write(line[1].rjust(13))
                fg.write(line[2].rjust(13))
                fg.write(line[3].rjust(13))
                fg.write('\n')
    f.close()
    fg.close()

def orcaout(fi, l, nqm, nlink, pcg):
    """
    
    // Function which extracts the energy, forces and Mulliken charges from an ORCA output file. //
    
    Arguments
    ----------
    fi : string
        ORCA job output file name.
    l : integer
        If using the nudged elastic band, this is the image number.
    nqm : integer
        Number of QM atoms.
    nlink : integer
        Number of link atoms.
    pcg : string
        ORCA support file name containing the calculated MM gradients.

    """
    
    # The output file containing QM energy and gradients is created.
    fd = open(('%s%d'%('ab_initio', l)), 'a')
    fd.write('Energy')
    fd.write('\n')
    
    # The output file containing QM Mulliken charges is created.
    fdm = open(('%s%d'%('mulliken', l)), 'a')
    
    # The ORCA output file is opened.
    f = open(fi, 'r')
    
    # Total atoms are initialised.
    ntot = nqm + nlink
    
    # While loop iterates through the whole ORCA output file.
    while 1:
        line = f.readline()
        if not line:break
        
        # The energy is written to the file ab_initio*.
        if line[:25] == 'FINAL SINGLE POINT ENERGY':
            qmE = float(line.split()[4])
            fd.write(str(qmE))
            fd.write('\n')
            
        # The Mulliken charges are extracted from the ORCA output file.
        # These Mulliken charges and the total charge are written to the file mulliken*.
        # Depending on the version of ORCA, the Mulliken charges data is labelled different so two cases must be implemented.    
        elif (line[:23] == 'MULLIKEN ATOMIC CHARGES') or (line[:44] == 'MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS'):
            if (line[:23] == 'MULLIKEN ATOMIC CHARGES'): temp = 4
            if (line[:44] == 'MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS'): temp = 5
            f.readline()
            for i in range(ntot):
                line = f.readline().split()
                fdm.write(line[3].rjust(10))
                fdm.write('\n')  
            line = f.readline().split()
            fdm.write('total charge')
            fdm.write('\n')
            fdm.write(line[temp].rjust(10)) 
            fdm.write('\n')
            fdm.close()    
            
        # The gradients corresponding to each cartesian coordinate are extracted from the ORCA output file.
        # This gradient is written to the file ab_initio*.    
        elif line[:18] == 'CARTESIAN GRADIENT': 
            fd.write('Gradient')
            fd.write('\n')
            f.readline()
            f.readline()
            for i in range(ntot):
                line = f.readline().split()
                fd.write(line[0].rjust(3))
                linex = float(line[3])
                fd.write(str(linex).rjust(20))
                liney = float(line[4])
                fd.write(str(liney).rjust(20))
                linez = float(line[5])
                fd.write(str(linez).rjust(20))
                fd.write('\n')  
            fd.close()
    f.close()
    
    # The point charge array is opened.
    ftc=open('%s%d%s'%('points_chg',l,'.pts'),'r')
    
    # The output file containing gradients corressponding to the MM point charges is created.
    fg=open('%s%d%s'%('qmlatgrad',l,'.out'),'a')
    
    # The MM gradients generated by ORCA are opened.
    ff=open(pcg,'r') 
    
    # Every point charge is added to the list cha.
    cha=[]
    for line in ftc:
        tcha=float(line.split()[0])
        cha.append(tcha)     
    ftc.close()
    
    # The number of point charges is written to qmlatgrad*.out.
    fg.write(' gradient output')
    fg.write('\n')
    fg.write(str(len(cha)))
    fg.write('\n')
    
    # The MM gradients are written to qmlatgrad*.out.
    ff.readline()
    for line in ff:
        fg.write(line)
    f.close()
    fg.close()

def qmorcamain(imn, cwd, usrdir, qmjob_prefix, nqm, nlink, cln, qmkey, cha_mul, extra_basis, orca_head):
    """
    
    // Function which combines the functions in this file for execution of an ORCA job. //
    // It will prepare the input file, run the ORCA job, and then extract the results from the output. //
    
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
    orca_head : string
        Contains the details for assigning memory and number of processors for a QM job.

    """
    
    # The files from the previous ORCA run are rearranged and renamed with the prefix 'old_'.
    if os.path.exists(cwd + ('%s%d'%('/image', imn)) + ('%s%s%d%s'%('/', qmjob_prefix, imn, '.in'))):
     try:
       shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '.in')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '.in'))) 
       shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '.log')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '.log')))
       shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '.pcgrad')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '.pcgrad')))
       shutil.copy(('%s%d%s'%('ptchg', imn, '.pc')), ('%s%d%s'%('old_ptchg', imn, '.pc')))
       shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '.gbw')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '.gbw')))
       os.remove('%s%d%s'%(qmjob_prefix, imn, '.in'))
       os.remove('%s%d%s'%(qmjob_prefix, imn, '.log'))
       os.remove('%s%d%s'%(qmjob_prefix, imn, '.pcgrad'))
       os.remove('%s%d%s'%('ptchg', imn, '.pc'))
       os.remove('%s%d%s'%(qmjob_prefix, imn, '.gbw'))
     except:
       qomend('ERROR, while rearrangeing previous cycle ORCA files for image : ' + str(imn), cwd, usrdir)
    
    # The ORCA input file is created.
    try:
       orcainp(imn, cwd, qmjob_prefix, qmkey, cha_mul, extra_basis, orca_head)	
    except:
       qomend('Error while creating ORCA input for image : ' + str(imn), cwd, usrdir)
    
    # The ORCA job is performed.  
    try:
       os.system('orca ' + ('%s%d%s'%(qmjob_prefix, imn, '.in')) + ' > ' + ('%s%d%s'%(qmjob_prefix, imn, '.log')))
    except:
       qomend('Error while running ORCA job at cycle: ' + str(cln) + ' image :' + str(imn), cwd, usrdir)
    
    # The results from the ORCA job are extracted and rearranged.
    try:
       fi = '%s%d%s'%(qmjob_prefix, imn, '.log')
       pcg = '%s%d%s'%(qmjob_prefix, imn, '.pcgrad')
       orcaout(fi, imn, nqm, nlink, pcg)	
       shutil.copy(('%s%d'%('ab_initio', imn)), cwd)
       shutil.copy(('%s%d'%('mulliken', imn)), cwd)
       shutil.copy(('%s%d%s'%('qmlatgrad', imn, '.out')), cwd) 
       os.remove('%s%d%s'%('qmlatgrad', imn, '.out'))
       os.remove('%s%d'%('ab_initio', imn))
       os.remove('%s%d'%('mulliken', imn))
    except:
       qomend('ERROR, while reading Energy or Gradient or Mulliken charge or electrostatic potential from QM output file of image : ' + str(imn), cwd, usrdir)    

def morcainp(imn, cwd, qmjob_prefix, qmkey, cha_mul1, cha_mul2, extra_basis, orca_head):
    """
    
    // Function which creates the input files required for the ORCA jobs for MECP calculations. //
    // This function is called in the function mecp_orcamain at each QoMMMa cycle. //
    
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
    orca_head : string
        Contains the details for assigning memory and number of processors for a QM job.

    """
    
    # The ORCA input file is created for the image.
    # This operation is performed for both states A and B.
    fda = open(('%s%d%s'%(qmjob_prefix, imn, '_A.in')), 'w')  
    fdb = open(('%s%d%s'%(qmjob_prefix, imn, '_B.in')), 'w')
    
    # The variable orca_head is written to the input file, which contains details for assigning memory and processors.
    # This operation is performed for both states A and B.
    if orca_head.lower() != 'none':
        fda.write(orca_head.strip())  
        fda.write('\n')
        fdb.write(orca_head.strip())  
        fdb.write('\n')
     
    # Using a series of if and else statements, the qmkey (level of theory) is written to the input file.
    # If a .gbw wavefunction file is found, then this is used as a guess.
    # These operations are performed for both states A and B.    
    cwd1 = cwd + '%s%d'%('/image', imn)
    if qmkey.lower() != 'none':
        if os.path.exists('%s%s%d%s'%('old_', qmjob_prefix, imn, '_A.gbw')):
            jobkey = '!' + qmkey + ' moread engrad'
        else:
            jobkey = '!' + qmkey + ' engrad'
    else:
        if os.path.exists('%s%s%d%s'%('old_', qmjob_prefix, imn, '_A.gbw')):
            jobkey = '!' + qmkey + ' moread engrad'
        else:
            jobkey = '!' + qmkey + ' engrad'
        fda.write(jobkey)
        fda.write('\n')
    if qmkey.lower() != 'none':
        if os.path.exists('%s%s%d%s'%('old_', qmjob_prefix, imn, '_B.gbw')):
            jobkey = '!' + qmkey + ' moread engrad'
        else:
            jobkey = '!' + qmkey + ' engrad'
    else:
        if os.path.exists('%s%s%d%s'%('old_', qmjob_prefix, imn, '_B.gbw')):
            jobkey = '!' + qmkey + ' moread engrad'
        else:
            jobkey = '!' + qmkey + ' engrad'
    fdb.write(jobkey)
    fdb.write('\n')
  
    # The point charge array is written to the input file.
    # This operation is performed for both states A and B.
    ptca = '%pointcharges "' + '%s%d%s'%('ptchga', imn, '.pc"')
    ptcb = '%pointcharges "' + '%s%d%s'%('ptchgb', imn, '.pc"')
    fda.write(ptca)
    fda.write('\n')
    fdb.write(ptcb)
    fdb.write('\n')
    
    # The guess wavefunction is written to the input file.
    # This operation is performed for both states A and B.
    if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.gbw')):    
        wavefunctiona='%moinp '+'"'+'%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.gbw'+'"')
        fda.write(wavefunctiona)
        fda.write('\n')
    if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.gbw')):    
        wavefunctionb='%moinp '+'"'+'%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.gbw'+'"')
        fdb.write(wavefunctionb)
        fdb.write('\n')
        
    # If a mixed bases set description is used, then it is written to the input file.
    # This operation is performed for both states A and B.      
    fda.write(extra_basis.strip())
    fdb.write(extra_basis.strip())
    fda.write('\n')
    fdb.write('\n')
    fda.write('* xyz ')
    fdb.write('* xyz ')
  
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
    fda.write('*')
    fda.write('\n')
    fda.write('\n')
    fda.write('\n')
    fda.write('\n')
    fdb.write('*')
    fdb.write('\n')
    fdb.write('\n')
    fdb.write('\n')
    fdb.write('\n')
    fda.close()
    fdb.close()
  
    # The point charge arrays are written to the files ptchga*.pc and ptchgb*.pc.
    f=open(('%s%d%s'%('charges',imn,'.xyz')),'r') 
    ffa=open(('%s%d%s'%('ptchga',imn,'.pc')),'a')
    ffb=open(('%s%d%s'%('ptchgb',imn,'.pc')),'a')
    chx='%s%d%s'%('points',imn,'gau.pts')
    numchx=len(open(chx,'r').readlines())
    ffa.write(str(numchx))
    ffb.write(str(numchx))
    ffa.write('\n')
    ffb.write('\n')
    for line in f:
        if line.strip() != '&pointch':
            if line.strip() !='&':
                line=line.split()
                ffa.write(line[0].rjust(10))
                ffb.write(line[0].rjust(10))
                ffa.write(line[1].rjust(13))
                ffb.write(line[1].rjust(13))
                ffa.write(line[2].rjust(13))
                ffb.write(line[2].rjust(13))
                ffa.write(line[3].rjust(13))
                ffb.write(line[3].rjust(13))
                ffa.write('\n')
                ffb.write('\n')
    f.close()
    ffa.close()
    ffb.close()

def mecp_orcamain(imn, cwd, usrdir, qmjob_prefix, nqm, nlink, cln, qmkey, cha_mul1, cha_mul2, extra_basis, orca_head):
    """
    
    // Function which combines the functions in this file for execution of ORCA jobs for MECP calculations. //
    // It will prepare the input files, run the ORCA jobs, and then extract the results from the outputs. //
    
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
    orca_head : string
        Contains the details for assigning memory and number of processors for a QM job.

    """
    
    # The files from the previous ORCA runs are rearranged and renamed with the prefix 'old_'.
    # This operation is performed for both states A and B. 
    if os.path.exists(cwd + ('%s%d'%('/image', imn)) + ('%s%s%d%s'%('/', qmjob_prefix, imn, '_A.in'))):
        try:
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_A.in')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_A.in'))) 
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_A.log')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_A.log')))
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_A.pcgrad')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_A.pcgrad')))
            shutil.copy(('%s%d%s'%('ptchga', imn, '.pc')), ('%s%d%s'%('old_ptchga', imn, '.pc')))
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_A.gbw')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_A.gbw')))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_A.in'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_A.log'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_A.pcgrad'))
            os.remove('%s%d%s'%('ptchga', imn, '.pc'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_A.gbw'))
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_B.in')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_B.in'))) 
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_B.log')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_B.log')))
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_B.pcgrad')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_B.pcgrad')))
            shutil.copy(('%s%d%s'%('ptchgb', imn, '.pc')), ('%s%d%s'%('old_ptchgb', imn, '.pc')))
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '_B.gbw')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '_B.gbw')))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_B.in'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_B.log'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_B.pcgrad'))
            os.remove('%s%d%s'%('ptchgb', imn, '.pc'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '_B.gbw'))
        except:
            qomend('error while rearrangeing ORCA files created in previous run at cycle :' + str(cln) + ' for image: ' + str(imn), cwd, usrdir)
      
    # The ORCA input files are created for both states A and B.             
    try:
        morcainp(imn, cwd, qmjob_prefix, qmkey, cha_mul1, cha_mul2, extra_basis, orca_head)
    except:
        qomend('Error while trying to prepare ORCA input files at cycle : ' + str(cln) + ' for image : ' + str(imn), cwd, usrdir)
     
    # The ORCA job for state A is performed.        
    try:
        os.system('orca ' + ('%s%d%s'%(qmjob_prefix, imn, '_A.in')) + ' > ' + ('%s%d%s'%(qmjob_prefix, imn, '_A.log')))
    except:
        qomend('Error while running ORCA job for state A at cycle : ' + str(cln) + ' for image : ' + str(imn), cwd, usrdir)
    
    # The ORCA job for state B is performed.        
    try:
        os.system('orca ' + ('%s%d%s'%(qmjob_prefix, imn, '_B.in')) + ' > ' + ('%s%d%s'%(qmjob_prefix, imn, '_B.log')))
    except:
        qomend('Error while running ORCA job for state B at cycle : ' + str(cln) + ' for image : ' + str(imn), cwd, usrdir)
       
    # The results from the ORCA jobs are extracted and rearranged.   
    # These operations are performed for both states A and B.        
    try:
        fi = '%s%d%s'%(qmjob_prefix, imn, '_A.log')
        pcg = '%s%d%s'%(qmjob_prefix, imn, '_A.pcgrad')
        orcaout(fi, imn, nqm, nlink, pcg)	
        fi = '%s%d%s'%(qmjob_prefix, imn, '_B.log')
        pcg='%s%d%s'%(qmjob_prefix, imn, '_B.pcgrad')
        orcaout(fi, imn, nqm, nlink, pcg)	
        shutil.copy(('%s%d'%('ab_initio', imn)), cwd)
        shutil.copy(('%s%d'%('mulliken', imn)), cwd)
        shutil.copy(('%s%d%s'%('qmlatgrad', imn, '.out')), cwd) 
        os.remove('%s%d%s'%('qmlatgrad', imn, '.out'))
        os.remove('%s%d'%('ab_initio', imn))
        os.remove('%s%d'%('mulliken', imn))
    except:
        qomend('Error while extracting results from ORCA output files at cycle : ' + str(cln) + ' for image : ' + str(imn), cwd, usrdir)