#!/usr/bin/python3

"""
// This is a helper QoMMMa Python file, which contains the functions used in execution of an xTB job. //
"""

# Global imports.
import os
import shutil

# Local imports.
from qomutil import qomend

def xtbpc(imn, cwd, qmjob_prefix):
    """
    
    // Function which creates the point charge file required for an xTB job. //
    // This function is called in the function qmxtbmain at each QoMMMa cycle. //
    
    Arguments
    ----------
    imn : integer
        If using the nudged elastic band, this is the image number.
    cwd : string
        The current working directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.

    """
    
    # The point charge file is created for the image.
    fd = open(('%s%d%s'%(qmjob_prefix, imn, '.pc')), 'w')  
    
    # Get number of point charges.
    num_pc = sum(1 for _ in open(('%s%d%s'%('charges', imn, '.xyz')), 'r')) - 2
  
    # The point charge array is written to the .pc file.
    f = open(('%s%d%s'%('charges', imn, '.xyz')), 'r') 
    fd.write(str(num_pc))
    fd.write('\n')
    for line in f:
        if line.strip() != '&pointch':
            if line.strip() != '&':
                line=line.split()
                fd.write(line[0].rjust(10))
                fd.write(line[1].rjust(13))
                fd.write(line[2].rjust(13))
                fd.write(line[3].rjust(13)) 
                fd.write('\n')
    f.close()
    fd.close()
    
    # Copy electric field coordinates (for not now really used as can't calculate electric field with xTB).
    shutil.copy(('%s%d%s'%('points', imn, 'gau.pts')), 'fort.20')

def xtbout(fi, fchg, l, nqm, nlink, fortf):
    """
    
    // Function which extracts the energy, forces and Mulliken charges from xTB output files. //
    
    Arguments
    ----------
    fi : string
        xTB job output file name.
    fchg : string
        xTB Mulliken charges file.
    l : integer
        If using the nudged elastic band, this is the image number.
    nqm : integer
        Number of QM atoms.
    nlink : integer
        Number of link atoms.
    fortf : string
        xTB support file name containing the electric field.

    """
    
    # The output file containing QM energy and gradients is created.
    fd = open(('%s%d'%('ab_initio' ,l)), 'a')
    fd.write('Energy')
    fd.write('\n')

    # The output file containing QM Mulliken charges is created.
    fdm = open(('%s%d'%('mulliken', l)), 'a')
    
    # The xTB output files are opened.
    f = open(fi, 'r')  
    fchrg = open(fchg, 'r') 
    
    # Total atoms are initialised.
    ntot = nqm + nlink
    
    # Open xTB output file.
    linez = f.readlines()
    
    # Remove line breaks...
    lines = []
    for sub in linez:
        lines.append(sub.replace("\n",""))
        
    # Iterate through the xTB output file
    for counter,line in enumerate(lines):
        # The total energy is read from the xTB output file.
        # This energy is written to the file ab_initio*.
        if 'current total energy' in line:
            ene = float(lines[counter+2])
            fd.write(str(ene))
            fd.write('\n')
            
        # The gradients corresponding to each cartesian coordinate are extracted from the xTB output file.
        # This gradient is written to the file ab_initio*.
        if 'current gradient' in line: 
            fd.write('Gradient')
            fd.write('\n')
            temp_counter = counter + 2
            for j in range(ntot):
                fd.write(str(j+1).rjust(3))
                gradx = str(float(lines[temp_counter]))
                grady = str(float(lines[temp_counter+1]))
                gradz = str(float(lines[temp_counter+2]))
                fd.write('            ' + gradx+ '     ' + grady+ '     ' +gradz)
                fd.write('\n')  
                temp_counter+=3
            fd.close()
    f.close()    

    # Open Mulliken output file.
    linez = fchrg.readlines()
    
    # Remove line breaks...
    lines = []
    for sub in linez:
        lines.append(sub.replace("\n",""))
        
    # Iterate through Mulliken output file.  
    temp = 0
    for line in lines:
        fdm.write(str(line).rjust(10))
        fdm.write('\n')
    fdm.write('total charge=')
    fdm.write('\n')
    fdm.write(line[temp].rjust(10)) 
    fdm.write('\n')
    fdm.close()  
    fchrg.close()

    # The active point charge array is opened.
    ftc = open('%s%d%s'%('charges', l, 'mol.xyz'), 'r')
    
    # The total number of active point charges in the MM optimisation is initialised.
    # The location of each point charge is added to the list cha_lst.
    cha_lst = []
    cha_count = 0
    for counter,line in enumerate(ftc):
        split_line = line.split()
        if (len(split_line) > 1):
            is_active = int(split_line[4])
            if is_active == 1:
                cha_count += 1   
                cha_lst.append(counter)
    ftc.close()

    # The output file containing gradients corressponding to the MM point charges (calculated by xTB) is created.
    fg = open('%s%d%s'%('qmlatgrad', l, '.out'), 'a')

    # The number of point charges is written to qmlatgrad*.out.
    fg.write(' gradient output')
    fg.write('\n')
    fg.write(str(cha_count))
    fg.write('\n')
    
    # The xTB calculated point charge gradients is opened.
    fpc = open('pcgrad','r')
    linez = fpc.readlines()
    lines = []
    for sub in linez:
        lines.append(sub.replace("\n",""))
    fpc.close()

    # Each gradient for MM point charges is grepped from the xTB point charge gradient file.
    # This gradient is then written to qmlatgrad*.out.
    for item in cha_lst:
        temp_grad_lst = []
        temp_grad_lst = lines[item-1].split()
        gradx = str(float(temp_grad_lst[0]))
        grady = str(float(temp_grad_lst[1]))
        gradz = str(float(temp_grad_lst[2]))
        fg.write(gradx+ '     ' + grady+ '     ' +gradz)
        fg.write('\n')
    fg.close()

def qmxtbmain(imn, cwd, usrdir, qmjob_prefix, qmxtb_job, nqm, nlink, cln, cha_mul):
    """
    
    // Function which combines the functions in this file for execution of an xTB job. //
    // It will prepare the input file, run the xTB job, and then extract the results from the output. //
    
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
    qmxtb_job : string
        The path to the xTB program.
    nqm : integer
        Number of QM atoms.
    nlink : integer
        Number of link atoms.
    cln : integer
        The QoMMMa cycle number.
    cha_mul : list
        The charge and multiplicity given in the format of a list of two numbers (i.e., cha_mul = [0,1]).

    """
    
    # The files from the previous xTB run are rearranged and renamed with the prefix 'old_'.
    # Any files which are not required are simply removed.
    if os.path.exists(cwd + ('%s%d'%('/image', imn)) + ('%s%d%s'%(qmjob_prefix, imn, '.pc'))):
        try:
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '.pc')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '.pc'))) 
            shutil.copy(('%s%d%s'%(qmjob_prefix, imn, '.pcem')), ('%s%s%d%s'%('old_', qmjob_prefix, imn, '.pcem')))
            shutil.copy(('%s%d%s'%('xtb_geom', imn, '.xyx')), ('%s%s%d%s'%('old_xtb_geom', imn, '.xyz')))
            shutil.copy(('%s%d%s'%('xtb_geom', imn, '.engrad')), ('%s%s%d%s'%('old_xtb_geom', imn, '.engrad')))
            shutil.copy('charges', 'old_charges')
            os.remove('%s%d%s'%(qmjob_prefix, imn, '.pc'))
            os.remove('%s%d%s'%(qmjob_prefix, imn, '.pcem'))
            os.remove('%s%d%s'%('xtb_geom', imn, '.xyx'))
            os.remove('%s%d%s'%('xtb_geom', imn, '.engrad'))
            os.remove('charges')
        except:
            qomend('ERROR, while rearranging previous cycle xTB files for image : ' + str(imn), cwd, usrdir)
   
    # Create point charge input files.
    try:
        xtbpc(imn, cwd, qmjob_prefix)	
        fd = open(('%s%d%s'%(qmjob_prefix, imn, '.pcem')), 'w')
        fd.write('$embedding')
        fd.write('\n')
        fd.write('   input=')
        fd.write('%s%d%s'%(qmjob_prefix, imn, '.pc'))
        fd.write('\n')
        fd.write('$end')
        fd.close()
    except:
        qomend('Error while creating xTB point charge files for image : ' + str(imn), cwd, usrdir)
        
    # Create 'proper' xyz file (total atoms on line 1, and title on line 2).
    try:
        fd = open(('%s%d%s'%('xtb_geom', imn, '.xyz')), 'w')
        fd.write(str(nqm + nlink))
        fd.write('\n')
        fd.write('xyz file for xTB calculation')
        fd.write('\n')
        f = open(('%s%d%s'%('qmgeom', imn, '.xyz')), 'r') 
        for line in f:
            split_line = line.split()
            atom_indice = split_line[0]
            atom_name = ''.join([i for i in atom_indice if not i.isdigit()])
            fd.write(atom_name)
            fd.write(str(split_line[1]).rjust(20))
            fd.write(str(split_line[2]).rjust(20))
            fd.write(str(split_line[3]).rjust(20))
            fd.write('\n')
        f.close()
        fd.close()
    except:
        qomend('Error while creating xTB geometry file for image : ' + str(imn), cwd, usrdir)
        
    # The xTB job is performed.
    # Note that in xTB, the spin is specified by alpha electrons - beta electrons so (cha_mul[1] - 1) represents the spin.
    try:
        os.system(qmxtb_job + ' ' + ('%s%d%s'%('xtb_geom', imn, '.xyz')) + ' --grad --pop --gfn 2 --chrg ' + str(cha_mul[0]) + ' --uhf ' + str(cha_mul[1] - 1) + ' -I ' + ('%s%d%s'%(qmjob_prefix, imn, '.pcem')))
    except:
        qomend('Error while running xTB job at cycle: ' + str(cln) + ' image :' + str(imn), cwd, usrdir)
    
    # The results from the xTB job are extracted and rearranged.
    try:
        fi = '%s%d%s'%('xtb_geom', imn, '.engrad')
        fchg = 'charges'
        xtbout(fi, fchg, imn, nqm, nlink, 'fort.21')	
        shutil.copy(('%s%d'%('ab_initio', imn)), cwd)
        shutil.copy(('%s%d'%('mulliken', imn)), cwd)
        shutil.copy(('%s%d%s'%('qmlatgrad', imn, '.out')), cwd) 
        os.remove('%s%d%s'%('qmlatgrad', imn, '.out'))
        os.remove('%s%d'%('ab_initio', imn))
        os.remove('%s%d'%('mulliken', imn))
    except:
        qomend('ERROR, while reading Energy or Gradient or Mulliken charge or electrostatic potential from QM output file of image : ' + str(imn), cwd, usrdir)    