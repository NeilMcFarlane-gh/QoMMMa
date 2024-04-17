#!/usr/bin/python3

"""
// This is the main QoMMMa Python file, which collects the functions defined in other Python files for QoMMMa execution. //
"""

# Global imports.
import os
import re
import sys
import shutil
import platform
import math
import random
import glob
from time import asctime
from pathlib import Path

# Local imports.
import qomutil
import gsmutil


def fortinp(opt_dir):
    """
    
    // Function which creates two input files (fortinput and converg.data) required for execution of the QoMMMa Fortran code. //
    // These files (fortinput in particular) are the primary communications between Python and Fortran code in QoMMMa. //
    // All of these parameters in their default state are found in /bin/default_inp.py. //
    
    Arguments
    ----------
    opt_dir : string
        The directory to be optimised.

    """
  
    # The file fortinput is created.
    fd = open('fortinput','w')
  
    # Number of QM, MM, link, and total QM optimised atoms are written to fortinput.
    # Also perform a quick check that QM region atom numbers add up as they should.
    fd.write('natom     nqm      nlink   nopt')
    fd.write('\n')
    fd.write(str(natom))
    fd.write(str(('%s%d'%('   ', nqm))))
    fd.write(str(('%s%d'%('   ', nlink))))
    fd.write(str(('%s%d'%('   ', nopt))))
    fd.write('\n')
    if ((nlink + nqm) !=  nopt) and (hesoptmm == 0):
        qomutil.qomend('FATAL ERROR: There are no hessian optimised MM atoms and the number of QM optimised atoms is ' + str(nopt) + ' , but the number of QM and link atoms is ' + str(nlink + nqm), cwd, opt_dir)

    # Whether dispersion energy correction for the QM region is calculated is written to fortinput.
    fd.write('disp')
    fd.write('\n')
    fd.write(str((disp)))
    fd.write('\n') 

    # The number of cartesian coordinate constraints as well as whether the constraint is harmonic (1) or tannic (2) is written to fortinput.
    fd.write('ncon_cart     kcnstype')
    fd.write('\n')
    fd.write(str(ncon_cart))
    if kcnstype.lower() == 'harmonic':
        ikcnstype = 1
    elif kcnstype.lower() == 'tanh':
        ikcnstype = 2
    else:
        qomutil.qomlog('Error, unknown constrain type : ' + kcnstype + ' is requested', opt_dir)
    fd.write(str(('%s%d'%('   ', ikcnstype))))
    fd.write('\n') 
    
    # The number of primitive coordinate constraints is written to fortinput.
    fd.write('ncon_prim')
    fd.write('\n')
    fd.write(str(ncon_prim))
    fd.write('\n') 
    
    # The number of primitive coordinate displacements are written to fortinput.
    fd.write('disp_prim')
    fd.write('\n')
    fd.write(str(disp_prim))
    fd.write('\n')  
    
    # The number of additional primitive coordinates to be added to the primitive coordinate array is written to fortinput.
    # This can be useful when using a bonds, angles, and torsions definition of primitive internal coordinates.
    fd.write('add_prims')
    fd.write('\n')
    fd.write(str(add_prims))
    fd.write('\n')  
    
    # The nudged elastic band method used as well as the spring constant is written to fortinput.
    if nebtype.lower() == 'none':
        inebtyp = 0
    elif nebtype.lower() =='neb_bfgs': # NEB with BFGS optimisation.
        inebtyp = 1
    elif nebtype.lower() == 'cineb_bfgs': # Climbing image NEB with BFGS optimsation.
        inebtyp = 2
    elif nebtype.lower() == 'neb_cg': # NEB with congugate gradient optimisation.
        inebtyp = 3
    elif nebtype.lower() == 'none_cg': # Not NEB (??) but with congugate gradient optimisation.
        inebtyp = 4
    elif nebtype.lower() =='cineb_cg': # Climbing image NEB with congugate gradient optimisation.
        inebtyp = 5
    else:
        qomutil.qomend('FATAL ERROR: unknown nebtype : ' + nebtype + ' is requested.', cwd, opt_dir)
    fd.write('nimg    nebtype    kspring')
    fd.write('\n')
    fd.write(str(nimg))
    fd.write(str(('%s%d'%('   ', inebtyp))))
    fd.write(str(('%s%d'%('   ', kspring))))
    fd.write('\n')

    # The growing string method version is written to fortinput.
    if gsmtype.lower() == 'none':
        igsmtyp = 0
    elif gsmtype.lower() =='de_gsm': # Double-ended GSM.
        igsmtyp = 1
    elif gsmtype.lower() =='se_gsm': # Single-ended GSM.
        igsmtyp = 2
    else:
        qomutil.qomend('FATAL ERROR: unknown gsmtype : ' + gsmtype + ' is requested.', cwd, opt_dir)
    fd.write('gsmtype')
    fd.write('\n')
    fd.write(str(igsmtyp))
    fd.write('\n')
    
    # The phase that the growing string method is currently in is written to fortinput.
    if gsmphase.lower() == 'none':
        igsmphase = 0
    elif gsmphase.lower() =='growth':
        igsmphase = 1
    elif gsmphase.lower() =='opt':
        igsmphase = 2
    elif gsmphase.lower() =='reparam':
        igsmphase = 3
    else:
        qomutil.qomend('FATAL ERROR: unknown gsmphase : ' + gsmphase + ' is requested.', cwd, opt_dir)
    fd.write('gsmphase')
    fd.write('\n')
    fd.write(str(igsmphase))
    fd.write('\n')
 
    # The type of coordinates used is written to fortinput.
    if coordtype.lower() == 'cart': # Cartesian.
        icoordtyp = 0
    elif coordtype.lower() =='dlc': # Delocalised internal coordinates.
        icoordtyp = 1
    else:
        qomutil.qomend('FATAL ERROR: unknown coordinate type : ' + coordtype + ' is requested.', cwd, opt_dir)
    fd.write('coordtype')
    fd.write('\n')
    fd.write(str(icoordtyp))
    fd.write('\n') 
    
    # The type of primitive internal coordinates used to generate delocalised internal coordinates is written to fortinput.
    # Only relevant where DLC are used, but the format is still written to fortinput for continuity.
    if primtype.lower() == 'tc':
        iprimtyp = 0
    elif primtype.lower() =='full':
        iprimtyp = 1
    else:
        qomutil.qomend('FATAL ERROR: unknown primitive internal coordinate type : ' + primtype + ' is requested.', cwd, opt_dir)
    fd.write('primtype')
    fd.write('\n')
    fd.write(str(iprimtyp))
    fd.write('\n') 

    # If they have been read in already from a previous QoMMMa run, then write the primitive internal coordinate definitions to fortinput.
    # Only relevant where DLC are used, but the format is still written to fortinput for continuity.
    fd.write('Total prims followed by definitions')
    fd.write('\n')
    fd.write(str(prim_num))
    fd.write('\n')
    if prim_num > 0:
        for i in prim_defs:
            for j in i:
                fd.write(str(j).rjust(8))
            fd.write('\n')
   
    # If there are any additional primitive internal coordinates to be added to the definition, then write these to fortinput.
    # These may overlap with the primitive internal coordinates calcualted by the Fortran code as any duplicates will be removed automatically.
    # Only relevant where DLC are used, but the format is still written to fortinput for continuity.
    if add_prims > 0:
        fd.write('Additional primitive internal coordinates to be included in the definition')
        fd.write('\n')
        for i in prim_add_lst:
            for j in i:
                fd.write(str(j).rjust(8))
            fd.write('\n')

    # The atom indices and types relating to the QM atoms is written to fortinput using the function qmread defined in qomutil.py.
    if nqm > 0:
        fd.write('QM atoms list')
        fd.write('\n')
        fd.close()
        try:		
            qomutil.qmread(qm_lst, opt_dir)
        except:
            qomutil.qomend('FATAL ERROR: something went wrong reading the definition of QM atoms - check the format of your input file!', cwd, opt_dir)
        fd = open('fortinput', 'a')        
    else:
        qomutil.qomend('FATAL ERROR: QM atoms need to be defined for a QM/MM calculation.', cwd, opt_dir)
 
    # Link atom details are written to fortinput using the function link_write defined in qomutil.py.    
    if nlink > 0:
        fd.write('Link atom details')
        fd.write('\n')
        fd.close()
        try:
            qomutil.link_write(link_lst) 
        except:
            qomutil.qomend('FATAL ERROR: something went wrong reading the definition of link atoms - check the format of your input file!', cwd, opt_dir)
        fd = open('fortinput', 'a')
    else:
        qomutil.qomlog('NOTE: definition of link atoms is not required since nlink=0.', opt_dir) 
    
    # The hessian optimised MM atoms are written to fortinput using the function hesopt_write defined in qomutil.py.
    if (int(nopt) - int(nqm) - int(nlink)) > 0:	
        fd.write('Hessian optimized MM atoms list')	
        fd.write('\n')
        fd.close()
        try:
            qomutil.hesopt_write(hesoptmm, opt_dir)
        except:
            qomutil.qomend('FATAL ERROR: something went wrong reading the definition of hessian optimised MM atoms - check the format of your input file!', cwd, opt_dir)
        fd = open('fortinput', 'a')
    else:
        qomutil.qomlog('NOTE: no MM atoms are included in Hessian optimisation.', opt_dir)
       
    # Constraint details are written to fortinput using the function cons_write_cart (for Cartesians) or cons_write_prim (for primitive internal coordinates) defined in qomutil.py.
    if ncon_cart > 0:
        if (icoordtyp == 0):
            fd.write('Constrain details - cartesian')
            fd.write('\n')
            fd.close()
            try:
                qomutil.cons_write_cart(cart_constrain_lst, opt_dir)  
            except:
                qomutil.qomend('FATAL ERROR: something went wrong reading the definition of cartesian constraints - check the format of your input file!', cwd, opt_dir)
            fd = open('fortinput', 'a')
        else:
            qomutil.qomend('FATAL ERROR: you may not mix coordinate systems and constraint types.', cwd, opt_dir)
    elif ncon_prim > 0:
        if (icoordtyp == 1):
            fd.write('Constrain details - primitive internal coordinates')
            fd.write('\n')
            fd.close()
            try:
                qomutil.cons_write_prim(prim_constrain_lst, prim_coeff_lst, ncon_prim, opt_dir)  
            except:
                qomutil.qomend('FATAL ERROR: something went wrong reading the definition of primitive internal coordinate constraints - check the format of your input file!', cwd, opt_dir)
            fd = open('fortinput', 'a')
        else:
            qomutil.qomend('FATAL ERROR: you may not mix coordinate systems and constraint types.', cwd, opt_dir)
    else:
        qomutil.qomlog('NOTE: no constraints are enforced.', opt_dir)
        
    # Displacements in primitive internal coordinates are written to fortinput using the function disp_write_prim defined in qomutil.py.
    if disp_prim > 0:
        fd.write('Requested displacements in primitive internal coordinates')
        fd.write('\n')
        fd.close()
        try:
            qomutil.disp_write_prim(prim_displacement_lst, opt_dir)  
        except:
            qomutil.qomend('FATAL ERROR: something went wrong reading the displacements in primitive internal coordinates - check the format of your input file!', cwd, opt_dir)
        fd = open('fortinput', 'a')
    else:
        qomutil.qomlog('NOTE: no displacements in primitive internal coordinates are requested.', opt_dir)
         
    # The modified charges in the MM region are written to fortinput using the function newcha_write defined in qomutil.py.
    fd.write('Modified charge for MM atoms. First, number of such atoms, and then atom label and charge')
    fd.write('\n')
    if newcha_lst != 'None':
        try:
            fd.close()
            qomutil.newcha_write(newcha_lst)
            fd = open('fortinput', 'a')
        except:
            qomutil.qomend('FATAL ERROR: something went wrong reading the new charges on MM atoms - check the format of your input file!', cwd, opt_dir)
    else:
        fd.write(str(0))
        fd.write('\n')
        qomutil.qomlog('NOTE: no charges are modified on MM atoms.', opt_dir)

    # The inactive MM atoms are written to fortinput using the function inprange defined in qomutil.py.
    fd.write('Inactive atoms. First, number of such atoms then list')
    fd.write('\n')
    kk = ninact.split()
    fd.write(str(kk[0]))		       
    if int(kk[0]) != 0:
        fd.write('\n')
        kk = kk[1:]
        fd.close()
        try:
            qomutil.inprange(kk, 'Inactive', opt_dir)
        except:
            qomutil.qomend('FATAL ERROR: something went wrong reading the inactive MM atoms - check the format of your input file!', cwd, opt_dir)
    else: 
        qomutil.qomlog('NOTE: no atoms are set as inactive.', opt_dir)
        fd.close()
   
    # The file converg.data is created and all convergence criteria are written to converg.data as it is used in later Fortran calculations.  
    fc = open('converg.data', 'w')
    fc.write('change in energy')
    fc.write('\n')
    fc.write(str(dEmax))
    fc.write('\n')
    fc.write('change in gradient')
    fc.write('\n')
    fc.write(str(dGmax))
    fc.write('\n')
    fc.write('change in rms gradient')
    fc.write('\n')
    fc.write(str(dGrms))
    fc.write('\n')
    fc.write('change in displacement')
    fc.write('\n')
    fc.write(str(dXmax))
    fc.write('\n')
    fc.write('change in rms displacement')
    fc.write('\n')
    fc.write(str(dXrms))
    fc.write('\n')
    fc.write('Tolerence in Perpendicular Gradient for NEB')
    fc.write('\n')
    fc.write(str(dGper))
    fc.write('\n')
    fc.close()
    
def mmjob(mjob, opt_dir):
    """
    
    // Function which performs a series of operations for the first and second Tinker jobs. //
    // This function is executed separately from the main QoMMMa cycle as it does not fit exactly into the framework of a typical cycle. //
    
    Arguments
    ----------
    mjob : string
        If mjob is 'initial', then it executes the initial Tinker calculation, and if 'initial2', then it prepares for subsequent QoMMMa cycles.
    opt_dir : string
        The directory to be optimised.

    """
    
    # At present, QoMMMa operates only with Tinker, so the variable mmcode must be 'Tinker'.
    # This function runs either the initial MM job (mjob = 'initial') or prepares for the second MM job (mmjob = 'initial2') using results from the Fortran QoMMMa setup program.
    if mmcode == 'Tinker':
        if mjob == 'initial':

            # The default forcefield parameters are written to the files geom.key and tinker_header so that Tinker always has specified forcefield parameters.
            default_prm = qommmadir + '/lib/charmm27.prm'
            fd = open('geom.key','w')      
            fd.write('parameters ' + default_prm) 
            fd.write('\n')		
            fdt = open('tinker_header','w')
            fdt.write('parameters ' + default_prm)
            fdt.write('\n')   	
           
            # If a specific MM key is specified, then it is written to geom.key and tinker_header and thus overwrites the default.
            try:		
                fd.write(mmkey.lstrip())  
                fdt.write(mmkey.lstrip())
                qomutil.qomlog('NOTE: the default force field parameter file is charmm27.prm, this can be changed or modified through mmkey in user input.', opt_dir)           
            except: 
                qomutil.qomlog('NOTE: no mmkey is given by user, default force field parameter file, charmm27.prm will be used here', opt_dir) 
            fdt.close()
           
            # Only charges are obtained for the subsequent MM polarised QM calculation.
            fd.write('CHARGETERM ONLY')
            fd.write('\n')
            fd.write('EXTRATERM NONE')
            fd.write('\n')
            fd.close()
           
            # For every image (if using the nudged elastic band, that is), Tinker operations are carried out.
            for i in range(nimg): 
                # A directory is created for each image, and files are copied over.
                l = i + 1
                dirnam = '%s%d' % ('image', l)
                os.mkdir(dirnam)
                fn = '%s%d%s' % ('geom', l, '.key')  
                shutil.copyfile('geom.key', fn) 
                fin = '%s%d%s' % (inpgeom_prefix,l, '.xyz')
                fdn = '%s%d%s' % ('geom',l,'.xyz')
                try:
                    old = opt_dir + '/' + fin         
                    shutil.copy(old, fin) 
                    shutil.copy(fin, fdn)
                except:
                    qomutil.qomend('FATAL ERROR: input file geometry file : ' + fin + ' is not found in user directory. Check your input geometry.', cwd, opt_dir)
                try:
                    qomutil.geom_expl(fdn, l)
                except:
                    qomutil.qomend('FATAL ERROR: could not create geom_expl*.xyz. Check your input geometry.', cwd, opt_dir)   
                    
                # To obtain the charges for the upcoming QM job, the Tinker program analyze_grad is run for each image using the function tinkercharge defined in qomutil.py.
                try:
                    qomutil.tinkercharge(qommmadir) 
                except:
                    qomutil.qomend('FATAL ERROR: could not generate initial charges using Tinker program analyze_grad. Check your input geometry.', cwd, opt_dir) 
                    
        # Preparing MM files for all subsequent MM jobs.      
        elif mjob == 'initial2':   
           try:
               # The details in tinker_header and inactive_list are written to geom.key.
               fin1 = open('tinker_header', 'r')
               fin2 = open('inactive_list', 'r')
               fdg = open('geom.key', 'w')
               for line in fin1:
                   fdg.write(line)
               for line in fin2:
                   fdg.write(line)
               fin1.close()
               fin2.close()
               fdg.close()
              
               # For every image (if using the nudged elastic band, that is), geom.key is copied.
               l = 0
               for ir in range(nimg):
                   l = l + 1
                   shutil.copy('geom.key', ('%s%d%s'%('geom', l, '.key'))) 
               os.remove('inactive_list') 
           except:
               qomutil.qomend('FATAL ERROR: could not generate geom.key file after QoMMMa initial run. This error is most likely related to the definition of inactive atoms - check your input file!', cwd, opt_dir)
        else:
            qomutil.qomend('FATAL ERROR: Unknown MM job. This is handled within the source code so something must have been changed.', cwd, opt_dir)
    else:
        qomutil.qomend('FATAL ERROR: Unknown MM code ' + mmcode + '. At present, only Tinker can be used.', cwd, opt_dir)
        
def qminitial(opt_dir):
    """
    
    // Function which creates the initial inputs for a QM job. //
    // Allowed QM programs in QoMMMa are: Jaguar, Molpro, Gaussian, xTB, or ORCA. //
    
    Arguments
    ----------
    opt_dir : string
        The directory to be optimised.

    """
    
    # Jaguar job files are prepared by using functions from jagutil.py.
    if qmcode.lower() == 'jaguar':
        try:
            # Normal Jaguar job is prepared.
            if job.lower() != 'mecp': 
                # The file qmheader is made as it is used in operation of a Jaguar job.
                fds = open('qmheader', 'w')
                fds.write(qmjag_header)   
                if qmkey.strip() != 'None':
                    fds.write(qmkey.lstrip())		
                fds.write('&')
                fds.write('\n') 
                fds.write('&zmat')
                fds.close()
                
            # Minimum energy crossing point Jaguar job is prepared.
            elif job.lower() == 'mecp':
                from jagutil import mjag_initial
                mjag_initial(qmjag_header, qmkey, cha_mul1, cha_mul2)
            
            # Regardless of Jaguar job type, it is prepared the function atomic defined in jagutil.py.
            if atomic_section != 'None':
                fd = open('atomic_section','w')
                fd.write(atomic_section.lstrip())
                fd.close()
            else: 
                from jagutil import atomic
                atomic()
            
            # For every image (if using the nudged elastic band, that is), appropriate files are copied to run the Jaguar jobs.
            imn = 0
            for i in range(nimg):
                imn = imn + 1
                dst = cwd + ('%s%d'%('/image', imn))
            if os.path.exists(opt_dir + '%s%d'%('/orbs', imn)):
                shutil.copy(opt_dir + '%s%d'%('/orbs', imn), dst)
                qomutil.qomlog('initial orbital guess for Jaguar job is copied from user directory for image :' + str(imn), opt_dir)
            if job.lower()=='mecp':
                if os.path.exists(opt_dir + '%s%d'%('/orbsA', imn)):
                    shutil.copy(opt_dir + '%s%d'%('/orbsA', imn), dst)
                    qomutil.qomlog('initial orbital guess for Jaguar job is copied from user directory for first state of image :' + str(imn), opt_dir)
                if os.path.exists(opt_dir + '%s%d'%('/orbsB', imn)):
                    shutil.copy(opt_dir + '%s%d'%('/orbsB', imn), dst)
                    qomutil.qomlog('initial orbital guess for Jaguar job is copied from user directory for second state of image :' + str(imn), opt_dir)
        except:
            qomutil.qomend('FATAL ERROR: initial Jaguar files could not be created. This may be due to problems in input files or version of QM code. Check the files in /jobfiles/.', cwd, opt_dir)
                
    # Molpro job is prepared.       
    elif qmcode.lower() == 'molpro':
        try:
            # The file qmheader is made as it is used in operation of a Molpro job.
            fds = open('qmheader', 'w')
            fds.write(qmmol_header.lstrip())
            fds.close()

            # Normal Molpro job is prepared.
            if job.lower() != 'mecp':
                # The file qmmol_footer is made as it is used in operation of a Molpro job.
                fds = open('qmmol_footer', 'w')
                if qmkey.strip() != 'None':
                    fds.write(qmkey.lstrip())
                else:
                    qomutil.qomend('FATAL ERROR: for Molpro, you have to give at least one input in qmkey (i.e., HF).', cwd, opt_dir) 
                fds.write(mol_footer.lstrip())
                fds.close()
                
            # Minimum energy crossing point Molpro job is prepared.
            elif job.lower() == 'mecp':
                # The files qmmol_footer_A and qmmol_footer_B are made as they are used in operation of a Molpro job.
                fda = open('qmmol_footer_A', 'w')
                fdb = open('qmmol_footer_B', 'w')
                if qmkey.strip() != 'None':
                    fda.write(qmkey.lstrip())
                    fdb.write(qmkey.lstrip())
                else:
                    qomutil.qomend('FATAL ERROR: for Molpro, you have to give at least one input in qmkey (i.e., HF).', cwd, opt_dir) 
                fda.write(cha_mul1)
                fdb.write(cha_mul2)
                fda.write('\n')
                fdb.write('\n')
                fda.write(mol_footer.lstrip())
                fdb.write(mol_footer.lstrip())
                fda.close()
                fdb.close()
                
                # For every image (if using the nudged elastic band, that is), appropriate files are copied to run the Molpro jobs. 
                imn = 0 
                for i in range(nimg):
                    imn = imn + 1
                    dst = cwd + ('%s%d'%('/image', imn))
                if os.path.exists(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.intg')):
                    dst = cwd + ('%s%d'%('/image', imn))
                    shutil.copy(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.intg'), dst)
                    qomutil.qomlog('initial orbital guess file is copied from user directory for image :' + str(imn), opt_dir)
                if job.lower() == 'mecp':
                    if os.path.exists(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.intg')):
                        shutil.copy(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.intg'), dst)
                        qomutil.qomlog('initial orbital guess file is copied from user directory for first mecp state of image :' + str(imn), opt_dir)
                    if os.path.exists(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.intg')):
                        shutil.copy(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.intg'), dst)
                        qomutil.qomlog('initial orbital guess file is copied from user directory for second mecp state of image :' + str(imn), opt_dir)
        except:
            qomutil.qomend('FATAL ERROR: initial Molpro files could not be created. This may be due to problems in input files or version of QM code. Check the files in /jobfiles/.', cwd, opt_dir)
                    
    # Gaussian job is prepared.               
    elif qmcode.lower() == 'gaussian':
        try:
            if qmkey.strip() == 'None':
                qomutil.qomlog( 'NOTE: no extra QM option is given by user through qmkey, so Gaussian job will run at HF/sto-3G level', opt_dir)

            # For every image (if using the nudged elastic band, that is), appropriate files are copied to run the Gaussian jobs.     
            imn = 0 
            for i in range(nimg):
                imn = imn + 1
                dst = cwd + ('%s%d'%('/image', imn))
                if os.path.exists(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.chk')):
                    shutil.copy(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.chk'), dst)
                    qomutil.qomlog('initial Gaussian check file was copied from user directory for image :' + str(imn), opt_dir)
                if job.lower() ==' mecp':
                    if os.path.exists(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.chk')):
                        shutil.copy(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.chk'), dst)
                        qomutil.qomlog( 'initial Gaussian check file was taken from user directory for MECP state A of image :' + str(imn), opt_dir) 
                    if os.path.exists(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.chk')):
                        shutil.copy(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.chk'), dst)
                        qomutil.qomlog('initial Gaussian check file was taken from user directory for MECP state B of image :' + str(imn), opt_dir) 
        except:
            qomutil.qomend('FATAL ERROR: initial Gaussian files could not be created. This may be due to problems in input files or version of QM code. Check the files in /jobfiles/.', cwd, opt_dir)
                    
    # xTB job is prepared.               
    # NOTE: MECP and frequency calculations are not possible in the current implementation of xTB within QoMMMa. 
    elif qmcode.lower() == 'xtb':
        try:
            if qmkey.strip() == 'None':
                qomutil.qomlog('NOTE: No input for qmkey, but remember that no basis set or level of theory is required for xTB. No worries then!', opt_dir)
            else:
                qomutil.qomlog('NOTE: qmkey selection: ' + str(qmkey) + ' has been ignored since it is not required for xTB.', opt_dir)

            # For every image (if using the nudged elastic band, that is), appropriate files are copied to run the xTB jobs.  
            # For xTB, this is probably not really necessary, but nice to have.   
            imn = 0 
            for i in range(nimg):
                imn = imn + 1
                dst = cwd + ('%s%d'%('/image', imn))
                if os.path.exists(opt_dir + '/xtbrestart'):
                    shutil.copy(opt_dir + + '/xtbrestart', dst)
                    qomutil.qomlog('initial xTB restart file was copied from user directory for image :' + str(imn), opt_dir)
                if job.lower() ==' mecp':
                    pass
        except:
            qomutil.qomend('FATAL ERROR: initial xTB files could not be created. This may be due to problems in input files or version of QM code. Check the files in /jobfiles/.', cwd, opt_dir)
             
    # Orca job is prepared.                         
    elif qmcode.lower() == 'orca':
        try:
            if qmkey.strip() == 'None':
                qomutil.qomlog( 'NOTE: no extra QM option is given by user through qmkey, so ORCA job will run at HF/sto-3G level', opt_dir)

            # For every image (if using the nudged elastic band, that is), appropriate files are copied to run the Orca jobs. 
            imn = 0 
            for i in range(nimg):
                imn = imn + 1
                dst = cwd + ('%s%d'%('/image', imn))
                if os.path.exists(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.gbw')):
                    shutil.copy(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.gbw'), dst)
                    qomutil.qomlog('initial ORCA gbw file was copied from user directory for image :' + str(imn), opt_dir)
                if job.lower() == 'mecp':
                    if os.path.exists(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.gbw')):
                        shutil.copy(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.gbw'), dst)
                        qomutil.qomlog( 'initial ORCA gbw file was taken from user directory for MECP state A of image :' + str(imn), opt_dir) 
                    if os.path.exists(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.gbw')):
                        shutil.copy(opt_dir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.gbw'), dst)
                        qomutil.qomlog('initial ORCA gbw file was taken from user directory for MECP state B of image :' + str(imn), opt_dir) 
        except:
            qomutil.qomend('FATAL ERROR: initial ORCA files could not be created. This may be due to problems in input files or version of QM code. Check the files in /jobfiles/.', cwd, opt_dir)
    else:
        qomutil.qomend('FATAL ERROR: Unknown QM code: ' + qmcode + '. You can choose between Jaguar, Molpro, Gaussian, xTB, and ORCA. Check your input file.', cwd, opt_dir) 

def qmmm(cln, opt_dir):
    """
    
    // Function which performs a QoMMMa QM/MM cycle. The steps of a QoMMMa cycle are as follows... //
        1. Execute Tinker's analyze_grad to obtain MM gradient.
        2. Prepare input for QM program chosen, and subsequently run a QM job and extracts results.
        3. Using Fortran code (see Fortran code for more details), perform a QoMMMa optimisation.
        4a. Check convergence, and if achieved, move program files to user directory and end the QoMMMa program.
        4b. Check convergence, and if not achieved, the MM microiterative optimisation will follow, and create/rearrange files and return to 1.
        
    // This cycle continues until convergence is achieved or the maximum number of cycles has been reached. //
    
    Arguments
    ---------- 
    cln : integer
        The current QoMMMa cycle number.
    opt_dir : string
        The directory to be optimised.


    """
      
    # For every image (if using the nudged elastic band, that is), a QoMMMa cycle is performed. 
    imn = 0
    for irm in range(nimg):
        imn = imn + 1
        if os.path.exists(cwd + ('%s%d'%('/update_geom', imn))):

            # The Tinker program analyze_grad is used to obtain the MM gradient and it is copied to the appropriate location.
            try:
                qomutil.qomlog('About to run ' + mmcode +  ' job at cycle : ' + str(cln) +  '  for image : ' + str(imn), opt_dir)
                fch = os.popen(qommmadir + '/bin/analyze_grad', 'w')	
                fn = ('%s%d%s'%('geom', imn, '.xyz'))
                fch.write(fn + '\nE\n')
                fch.close()
                shutil.copy('mm_grad', ('%s%d'%('mm_grad', imn)))
                os.remove('mm_grad')
            except:
                qomutil.qomend('FATAL ERROR: problem with MM gradient extraction at step :' + str(cln), cwd, opt_dir)  
              
            # The QM job is performed using the program of the user's choice.
            os.chdir(('%s%d'%('image', imn)))
            qomutil.qomlog('About to run ' + qmcode +  ' job at cycle : ' + str(cln) +  '  for image : ' + str(imn), opt_dir)
              
            # Jaguar QM job.
            if qmcode.lower() == 'jaguar':
                try:
                    if job.lower() != 'mecp':
                        from jagutil import qmjagmain
                        qmjagmain(imn, cwd, opt_dir, qmjob_prefix, qmjag_job, nqm, nlink, cln)
                        qmc_job = qmjag_job    # Used later in frequency call.
                    elif job.lower() == 'mecp':
                        from jagutil import mecp_jagmain
                        mecp_jagmain(imn, cwd, opt_dir, qmjob_prefix, qmjag_job, nqm, nlink, cln)
                except:
                    qomutil.qomend('FATAL ERROR: problem with QM gradient extraction (Jaguar) at step :' + str(cln), cwd, opt_dir)
                      
            # Molpro QM job.       
            elif qmcode.lower() == 'molpro':
                try:
                    if job.lower() != 'mecp':
                        from molutil import qmmolmain
                        qmmolmain(imn, cwd, opt_dir, qmjob_prefix, qmmol_job, nqm, nlink, cln)
                        qmc_job = qmmol_job    # Used later in frequency call. 
                    elif job.lower() == 'mecp':
                        from molutil import mecp_molmain
                        mecp_molmain(imn, cwd, opt_dir, qmjob_prefix, qmmol_job, nqm, nlink, cln)
                except:
                    qomutil.qomend('FATAL ERROR: problem with QM gradient extraction (Molpro) at step :' + str(cln), cwd, opt_dir)
                    
            # Gaussian QM job.       
            elif qmcode.lower() == 'gaussian':
                try:
                    if job.lower() != 'mecp':
                        from gauutil import qmgaumain
                        qmgaumain(imn, cwd, opt_dir, qmjob_prefix, qmgau_job, nqm, nlink, cln, qmkey, cha_mul, extra_basis, gau_head, extra_guess)
                        qmc_job = qmgau_job    # Used later in frequency call. 
                    elif job.lower() == 'mecp':
                        from gauutil import mecp_gaumain
                        mecp_gaumain(imn, cwd, opt_dir, qmjob_prefix, qmgau_job, nqm, nlink, cln, qmkey, cha_mul1, cha_mul2, extra_basis, gau_head, extra_guess)
                except:
                    qomutil.qomend('FATAL ERROR: problem with QM gradient extraction (Gaussian) at step :' + str(cln), cwd, opt_dir)         

            # xTB QM job.
            # NOTE: MECP and frequency calculations are not possible in the current implementation of xTB within QoMMMa.
            elif qmcode.lower() == 'xtb':
                try:
                    if job.lower() != 'mecp':
                        from xtbutil import qmxtbmain
                        qmxtbmain(imn, cwd, opt_dir, qmjob_prefix, qmxtb_job, nqm, nlink, cln, cha_mul)   
                    elif job.lower() == 'mecp':
                        pass
                except:
                    qomutil.qomend('FATAL ERROR: problem with QM gradient extraction (xTB) at step :' + str(cln), cwd, opt_dir)    
                      
            # Orca QM job.         
            elif qmcode.lower() == 'orca':
                try:
                    if job.lower() != 'mecp':
                        from orcautil import qmorcamain
                        qmorcamain(imn, cwd, opt_dir, qmjob_prefix, nqm, nlink, cln, qmkey, cha_mul, extra_basis, orca_head)
                        qmc_job = qmorca_job    # Used later in frequency call. 
                    elif job.lower() == 'mecp':
                        from orcautil import mecp_orcamain
                        mecp_orcamain(imn, cwd, opt_dir, qmjob_prefix, nqm, nlink, cln, qmkey, cha_mul1, cha_mul2, extra_basis, orca_head)
                except:
                    qomutil.qomend('FATAL ERROR: problem with QM gradient extraction (ORCA) at step :' + str(cln), cwd, opt_dir)

            # From geom*.xyz, the file geom_expl*.xyz is created.
            os.chdir(cwd)
            fi = '%s%d%s'%('geom', imn, '.xyz')       
            try:
                qomutil.geom_expl(fi, imn) 	        
            except:
                qomutil.qomend('FATAL ERROR: could not create geom_expl*.xyz following MM job. Check your previous cycle geometry.', cwd, opt_dir)      
                os.remove(('%s%d'%('update_geom', imn)))
        else: pass
    
    # Starting QM/MM optimisation.
    # This is the point where the Fortran code is used, so it should be referenced to understand the optimisation process.
    qomutil.qomlog('About to run QM/MM optimization cycle :' + str(cln), opt_dir)
    try:		
        if job.lower() == 'mecp':						
            os.system(qommmadir + '/bin/qommma_mecp.x')
        else:
            os.system(qommmadir + '/bin/qommma8.x')  	
    except:
        qomutil.qomend('FATAL ERROR: could not run QoMMMa program at cycle :' + str(cln), cwd, opt_dir)
        
    # Update the report file with the results from the QoMMMa optimisation, and check for convergence either by reaching maxcycle or with a properly converged optimisation.
    try: 
        qomutil.qomreport(nimg, opt_dir, cwd, cln, qomout)
        converg_achieved = convergence(cln, opt_dir)
        if converg_achieved == True:
            return converg_achieved
    except:
        qomutil.qomend('FATAL ERROR: problem checking convergence after QoMMMa optimization at cycle :' + str(cln) + ' check QM, MM results or checkfile.', cwd, opt_dir)
        
    # If convergence has not been achieved, rearrange some files and prepare for microiteration.        
    try: 
        qomutil.rearrange(nimg, cwd, qmjob_prefix, cln)	 
    except:
        qomutil.qomend('FATAL ERROR: could not rearrange the files after QoMMMa optimization at cycle :' + str(cln) + '. Check the files in /jobfiles/.', cwd, opt_dir)
    
    # For every image (if using the nudged elastic band, that is), a microiteration is performed.
    im = 0
    for ir in range(nimg):				 
        try:
            im = im + 1
            if os.path.exists(cwd + ('%s%d'%('/update_geom', im))):
                dst = cwd + ('%s%d'%('/image', im))
                qomutil.qomlog('About to run Microoptimization of geometry :' + str(im), opt_dir)
                os.chdir(dst)

                # The Tinker program minimize is used to perform a minimisation of the MM region.
                try:				
                    fch = os.popen(qommmadir + '/bin/minimize > tmpmic.out', 'w')
                    fn = ('%s%d%s'%('microgeom', im, '.xyz'))
                    fch.write(fn + '\n' + str(acctink) + '\n')
                    fch.close()
                    shutil.copy('tmpmic.out', ('%s%d'%('micro_opt_log_', im)))
                    os.remove('tmpmic.out')  
                except:
                    qomutil.qomend('FATAL ERROR: Tinker minimization job at cycle :' + str(cln), cwd, opt_dir)
                 
                # Using microgeom.xyz_2 created by Tinker's minimize, geom_expl.xyz is created for further cycles.
                os.remove('gradcorrection')
                fi = '%s%d%s'%('microgeom', im, '.xyz_2')
                try:
                    qomutil.geom_expl(fi, im) 		  	 
                except:
                    qomutil.qomend('FATAL ERROR: could not create geom_expl*.xyz during microiterative step. Check your previous cycle geometry.', cwd, opt_dir)    
                shutil.copy(('%s%d%s'%('geom_expl', im, '.xyz')), cwd + ('%s%d%s'%('/geom_expl',im, '.xyz')))
                os.remove(('%s%d%s'%('microgeom', im, '.xyz_2')))
                os.remove('%s%d%s'%('geom_expl', im, '.xyz'))
                os.chdir(cwd)
            else:
                qomutil.qomlog('image Geometry  ' + str(im) + '  is fixed in this iteration', opt_dir)
        except:
             qomutil.qomend('FATAL ERROR: micro-iteration at QM/MM cycle :' + str(cln), cwd, opt_dir)
             
    # Using the Fortran code, the post microiteration manipulation of files is performed.     
    try:
        os.system(qommmadir + '/bin/postmicro8.x') 
    except:
        qomutil.qomend('FATAL ERROR: problem in QoMMMa post-processing job using /bin/postmicro8.x at cycle :' + str(cln), cwd, opt_dir)
    
    # Files are rearranged using the function postrearrange defined in qomutil.py.
    try:
        qomutil.postrearrange(nimg, cwd)
    except:
        qomutil.qomend('FATAL ERROR: problem rearranging the files after QoMMMa post-processing job at cycle :' + str(cln), cwd, opt_dir)  

def convergence(cln, opt_dir):
    """
    
    // Function which checks the convergence for a given QoMMMa run. //
    
    Arguments
    ----------
    cln : integer
        The current QoMMMa cycle number.    
    opt_dir : string
        The directory to be optimised.
    """

    # If the file 'convergence_ok' exists and hence convergence has been achieved, then the program is ended.
    if os.path.exists(cwd + '/convergence_ok'):
        converg_achieved = True

        # Copy geometry to the user directory.
        im = 0
        os.remove(cwd + '/convergence_ok')
        for ir in range(nimg):
            im = im + 1
            shutil.copy(cwd + ('%s%d%s'%('/geom', im, '.xyz')), opt_dir + ('%s%d%s'%('/geom', im, '.xyz')))

        # If requested, perform a frequency calculation.
        if job.lower() == 'freq':
            hlp = asctime()  
            qomutil.qomlog('About to start QoMMMa Frequency calculation at  ' + str(hlp), opt_dir) 
            from qomfreq import freqmain
            freqmain(qommmadir, qmkey, qmc_job, opt_dir, cha_mul, extra_basis, gau_head, prj_freq, qmjob_prefix, qmcode) 
            
        # QoMMMa optimisation complete!
        if gsmtype.lower() == 'none':
            return converg_achieved 
        elif (gsmtype.lower() == 'de_gsm' or 'se_gsm'):
            return converg_achieved   
        elif scantype.lower() == 'adiabatic': 
            return converg_achieved 

    # For the growing string method, there can be cases where a QoMMMa optimisation has run out of cycles and it has not converged.
    # In these cases, the calculation should continue even though it has not converged.
    elif (gsmtype.lower() == 'de_gsm' or 'se_gsm') and (cln == maxcycle):
        converg_achieved = True
        return converg_achieved

    else:
        converg_achieved = False
        return converg_achieved

def QoMMMa_opt(opt_dir):
    """
    
    // Function which executes a QoMMMa optimisation protocol for a given input geometry and qommma.in file. //
    
    Arguments
    ----------
    opt_dir : string
        The directory to be optimised.
        
    """
    
    # The time that the job starts is written to the log file.
    hlp = asctime()
    hlp1 = platform.uname()
    qomutil.qomlog('QoMMMa job starts at   ' + str(hlp) + '\n' + '  and running in node  ' + str(hlp1), opt_dir)
    
    # The directory 'qmmm_*' is created at TMPDIR. Typically, this will be a scratch directory on a node.
    # A random number between 1 and 200 is added to the end so that directories do not overlap.
    try:
        dirs = os.environ['TMPDIR'] + '/qmmm_'
        dirs = dirs + str(os.getpid()) + '_' + str(random.randint(1,200))
        if os.path.exists(dirs):
            shutil.rmtree(dirs)
        os.mkdir(dirs)
    except:
        qomutil.qomlog('FATAL ERROR: either TMPDIR is not properly defined or could not create temporary directory ' + dirs, opt_dir)
        sys.exit()
    
    # Directories are reset and the current working directory is set.
    os.chdir(dirs)
    global cwd
    cwd = os.getcwd()

    # The initial MM job is performed.
    mmjob('initial', opt_dir)
    
    # The primary file which communicates between Python and Fortran, fortinput, is generated.
    fortinp(opt_dir)	
    
    # The Fortran based QoMMMa "setup" program is run. This creates the building blocks for QM & MM jobs, and initializes output.
    try:
        os.system(qommmadir + '/bin/setup8.x')
    except:
        qomutil.qomend('FATAL ERROR: QoMMMa initial setup program /bin/setup8.x failed.', cwd, opt_dir) 

    # Initial QM input files are prepared.
    try:
        qminitial(opt_dir)	
    except:
        qomutil.qomend('FATAL ERROR: generation of initial QM input files failed.', cwd, opt_dir) 

    # Results from the QoMMMa setup program are moved.
    try:
        qomutil.setini(cwd, nimg, opt_dir) 	
    except:
        qomutil.qomend('FATAL ERROR: organisation of files before start of QoMMMa cycle failed.', cwd, opt_dir) 
       
    # New MM key files are creates using results from the QoMMMa setup program.
    try:
        mmjob('initial2', opt_dir)
    except:
        qomutil.qomend('FATAL ERROR: generation of MM key files failed.', cwd, opt_dir) 
       
    # Now begin the QM/MM cycles. Most of the work is done by the above function qmmm().
    cln = 0			
    for ic in range(maxcycle): 		
        cln = cln + 1
        # Run a QM/MM cycle and check the convergence.
        if gsmtype.lower() == 'none':
            converg_achieved = qmmm(cln, opt_dir)
            if converg_achieved == True:             
                qomutil.qomend('Congratulations,', cwd, opt_dir)
        elif gsmtype.lower() == 'se_gsm' or 'de_gsm':
            converg_achieved = qmmm(cln, opt_dir)
            if converg_achieved == True: 
                qomutil.qomsave(cwd, opt_dir)             
                return # for GSM, overall convergence is checked by a separate function.
        elif scantype.lower() == 'adiabatic':
            converg_achieved = qmmm(cln, opt_dir)
            if converg_achieved == True:  
                qomutil.qomsave(cwd, opt_dir)           
                return # for a scan, overall convergence is checked by a separate function.

        # For the case where the calculation does not converge in the maximum number of cycles, the program is stopped.
        if cln == maxcycle: 
            qomutil.qomend('Unfortunately, QoMMMa has reached maxcycle :' + str(cln) + '. Convergence has not been achieved so the program has been stopped.', cwd, opt_dir)

def QoMMMa_DE_gsm():
    """
    
    // Function which runs the DE_GSM workflow. //
    
    Arguments
    ----------

    """

    # The time that the GSM run starts is written to the GSM log file.
    hlp = asctime()
    hlp1 = platform.uname()
    gsmutil.gsmlog('GSM QoMMMa job starts at   ' + str(hlp) + '\n' + '  and all calculations are running in node  ' + str(hlp1), usrdir)

    # Details of GSM job.
    gsmutil.gsmlog('You have selected to use GSM type: ' + str(gsmtype), usrdir)
    if (total_nodes == 0):
        gsmutil.gsmend('FATAL ERROR: for a double-ended GSM calculation, you must provide the total number of nodes in qommma.in', usrdir)        
    gsmutil.gsmlog('With the GSM type: ' + str(gsmtype) + '  you have chosen to use ' + total_nodes + '  nodes in total.', usrdir)
    
    # The growth phase commences.
    gsmutil.gsmlog('        /// Growing the string... ///       ', usrdir)
    all_nodes = []
    is_grown = False
    while is_grown == False:

        # Initialising product and reactant node directories.
        nodeR_dir = usrdir + '/nodeR'
        if (os.path.exists(nodeR_dir)) is not True:
            gsmutil.gsmend('FATAL ERROR: for a double-ended GSM calculation, a reactant directory called ''nodeR'' must be in the working directory.', usrdir)
        nodeP_dir = usrdir + '/nodeP'
        if (os.path.exists(nodeP_dir)) is not True:
            gsmutil.gsmend('FATAL ERROR: for a double-ended GSM calculation,a product directory called ''nodeP'' must be in the working directory.', usrdir)

        # Initialising primitive internal coordinates for the reaction pathway so that they can be written to fortinput.
        # It is essential that the primitive internal coordinates remain the same across the pathway for coordinate continuity.
        prims_r = nodeR_dir + '/jobfiles/prim_list'
        if (os.path.exists(prims_r)) is not True:
            gsmutil.gsmend('FATAL ERROR: could not find the primitive definitions in /nodeR/jobfiles. The reactant node must first be optimised with delocalised internal coordinates.', usrdir)
        prims_p = nodeP_dir + '/jobfiles/prim_list'
        if (os.path.exists(prims_p)) is not True:
            gsmutil.gsmend('FATAL ERROR: could not find the primitive definitions in /nodeP/jobfiles. The product node must first be optimised with delocalised internal coordinates.', usrdir)
        global prim_num, prim_defs 
        prim_num, prim_defs = gsmutil.DE_initialise_prims(prims_r, prims_p, usrdir) 

        # If any additional primitive internal coordinates were requested in qommma.in, then these are added to prim_defs now.
        global additional_prims
        additional_prims=[]
        prim_num, prim_defs = gsmutil.additional_prims(additional_prims, prim_num, prim_defs, usrdir) 

        # Initially, the frontier nodes are the reactant and product nodes.
        frontierR_dir = nodeR_dir
        frontierP_dir = nodeP_dir
        all_nodes.append(nodeR_dir)
        all_nodes.append(nodeP_dir)

        # Now begin the loop of adding nodes. The current number of nodes is initially 2 due to reactant and product nodes.
        current_nodes = 2
        for (i,j) in zip(range(2, total_nodes, 1), range(total_nodes - 1, math.ceil(total_nodes / 2), -1)):
            # Creating new directories for new nodes.
            new_dirs = []
            frontier_dirs = []
            new_frontierR_dir = usrdir + '/node' + str(i)
            new_frontierP_dir = usrdir + '/node' + str(j)
            new_dirs.append(new_frontierR_dir)
            new_dirs.append(new_frontierP_dir)
            frontier_dirs.append(frontierR_dir)
            frontier_dirs.append(frontierP_dir)
            all_nodes.append(new_frontierR_dir)
            all_nodes.append(new_frontierP_dir)
            os.mkdir(new_frontierR_dir)
            os.mkdir(new_frontierP_dir)
            gsmutil.gsmlog('Node ' + str(i) + '   has been added to the string (reactant side)...', usrdir)
            gsmutil.gsmlog('Node ' + str(j) + '   has been added to the string (product side)...', usrdir)
            current_nodes += 2

            # To define the coordinates of the new nodes, obtain the difference in primitive internal coordinates between the reactant and product side frontier nodes.
            # New nodes are then added along this tangent.
            tangent, tangent_prims = gsmutil.DE_get_tangent(frontier_dirs, current_nodes, total_nodes, usrdir)

            # Now create two new qommma.in files for these new nodes and copy the geometries of the previous nodes over.
            gsmutil.DE_add_nodes(frontier_dirs, new_dirs, tangent, tangent_prims, usrdir, to_add=2)

            # Update the frontier directories to represent the new nodes.
            frontier_dirs = []
            frontierR_dir = new_frontierR_dir
            frontierP_dir = new_frontierP_dir
            frontier_dirs.append(frontierR_dir)
            frontier_dirs.append(frontierP_dir)

            # Now, the new nodes can be optimised subject to the conditions of a growth phase calculation.
            # To give the correct input parameters, read in the new qommma.in.
            for counter,dir in enumerate(frontier_dirs):
                inpf = dir + '/qommma.in'
                try:
                    exec(open(inpf).read(),globals())
                except:
                    gsmutil.gsmend('FATAL ERROR: could not open input file for node: ' + inpf + '  This may be due to error in user input file, or maybe it was not generated? Check for necessary inputs and its formats, see manual.', usrdir)
                gsmutil.gsmlog('Optimising node (growth-phase): ' + inpf + '........', usrdir) 
                if counter == 1:
                    QoMMMa_opt(frontierR_dir) # performing QoMMMa growth-phase optimisation.
                elif counter == 2:
                    QoMMMa_opt(frontierP_dir) # performing QoMMMa growth-phase optimisation.

            # Whenever new nodes are added, the string is reparameterised.
            # This is to say that it is ensured that the nodes are evenly spaced along the reaction path tangent.
            gsmutil.DE_reparam_g(all_nodes, current_nodes, total_nodes, usrdir)
            gsmutil.gsmlog('Reparameterising the string... All nodes should be evenly spaced', usrdir)

        # For an odd number of nodes, one final node from the reactant side must be added and then the growing terminated.
        if current_nodes < total_nodes:
            # Creating new directory for the central node.
            new_dirs = []
            frontier_dirs = []
            new_frontier_dir = usrdir + '/node' + str(math.ceil(total_nodes / 2))
            new_dirs.append(new_frontier_dir)
            frontier_dirs.append(frontierR_dir)
            frontier_dirs.append(frontierP_dir)
            all_nodes.append(new_frontier_dir)
            os.mkdir(new_frontier_dir)
            gsmutil.gsmlog('The final central node ' + str(math.ceil(total_nodes / 2)) + '   has been added to the string...', usrdir)
            current_nodes += 1

            # To define the coordinates of the central node, obtain the difference in primitive internal coordinates between the reactant and product side frontier nodes.
            # The node is then added along this tangent.
            tangent, tangent_prims = gsmutil.DE_get_tangent(frontier_dirs, current_nodes, total_nodes, usrdir)

            # Now create two new qommma.in files for these new nodes and copy the geometries of the previous nodes over.
            # The generation of new geometries (mathematically intensive) is handled by the Fortran code.
            gsmutil.DE_add_nodes(frontier_dirs, new_dirs, tangent, usrdir, to_add=1)

            # Update the frontier directory to represent the new node.
            old_frontier_dir = frontier_dir
            frontier_dir = new_frontier_dir

            # Now, the central node can be optimised subject to the conditions of a growth phase calculation.
            # To give correct input, read in the new qommma.in.
            inpf = frontier_dir + '/qommma.in'
            try:
                exec(open(inpf).read(), globals())
            except:
                gsmutil.gsmend('FATAL ERROR: could not open input file for central node: ' + inpf + '  This may be due to error in user input file, or maybe it was not generated? Check for necessary inputs and its formats, see manual', usrdir)
            gsmutil.gsmlog('Optimising node (growth-phase): ' + inpf + '........', usrdir)
            QoMMMa_opt(frontier_dir) # performing QoMMMa growth-phase optimisation.

            # Finally, the string is reparameterised one more time.
            gsmutil.DE_reparam_g(all_nodes, current_nodes, total_nodes, usrdir)
            gsmutil.gsmlog('Reparameterising the string... All nodes should be evenly spaced', usrdir)

        # Growth phase complete!
        if current_nodes == total_nodes:
            gsmutil.gsmlog('The growth phase is complete! There is now a total (including reactant and product) of: ' + str(current_nodes) + '   nodes...', usrdir)
            is_grown = True

    # Now move to the optimisation phase.
    gsmutil.gsmlog('        /// Optimising the string... ///       ', usrdir)
    is_optimised = False
    reparam_counter = 0
    while is_optimised == False:

        # First, generate a list of all the tangents for ease of working.
        # Tangents are defined between neighbouring nodes and the string is optimised with stricter convergence criteria than in the growth phase.
        tangent_list, tangent_prims_list = gsmutil.get_tangents_opt(all_nodes, usrdir, driving_coords)
        
        # Now, optimise each node with stricter convergence criteria than in the growth phase.
        # Create the qommma.in files for every node and then optimise each.
        gsmutil.gen_input_opt(all_nodes, tangent_list, tangent_prims_list, usrdir)
        for counter,dir in enumerate(all_nodes):
            if not ("nodeR" or "nodeP") in dir:
                inpf = dir + '/qommma.in'
                try:
                    exec(open(inpf).read(), globals())
                except:
                    gsmutil.gsmend('FATAL ERROR: Could not open input file for node: ' + inpf + '  This may be due to error in user input file, or maybe it was not generated? Check for necessary inputs and its formats, see manual', usrdir)
                gsmutil.gsmlog('Optimising node (optimisation-phase): ' + inpf + '........', usrdir)
                QoMMMa_opt(dir) # performing QoMMMa optimisation-phase optimisation.
         
        # Now, check for convergence of all nodes.
        # If all are converged, then the optimisation-phase is complete.
        # Otherwise, reparameterise the string and optimise all nodes again.
        is_converged = gsmutil.check_convergence(all_nodes, current_nodes, usrdir)

        if (all(is_converged) and (reparam_counter > 0)) is True: # we should at least reparameterise once due to nodes not being evenly spaced...
            gsmutil.gsmlog('The optimisation phase is complete!', usrdir)
            is_optimised = True
        else:
            reparam_counter += 1
            gsmutil.gsmlog('Performing reparameterisation number: ' + str(reparam_counter) + ', stand by...', usrdir)

            # First, generate the new qommma.in files.
            gsmutil.reparam_opt(all_nodes, current_nodes, usrdir)

            # Now, reparameterise every node.
            for counter,dir in enumerate(all_nodes):
                if not ("nodeR" or "nodeP") in dir:
                    inpf = dir + '/qommma.in'
                    try:
                        exec(open(inpf).read(), globals())
                    except:
                        gsmutil.gsmend('Could not open input file for node: ' + inpf + '  This may be due to error in user input file, or maybe it was not generated? Check for necessary inputs and its formats, see manual', usrdir)
                        
                    # The directory 'qmmm_*' is created at TMPDIR. Typically, this will be a scratch directory on a node.
                    try:
                        dirs = os.environ['TMPDIR'] + '/qmmm_'
                        dirs = dirs + str(os.getpid()) + '_' + str(counter)
                        if os.path.exists(dirs):
                            shutil.rmtree(dirs)
                        os.mkdir(dirs)
                    except:
                        gsmutil.gsmlog('Either Variable TMPDIR is not set or Could not be able to create temporary directory ' + dirs, usrdir)
                        sys.exit()
                    os.chdir(dirs)
                    
                    # Copy over the current geometries...
                    source = dir + '/jobfiles'
                    destination = dirs
                    shutil.copy(source + '/geom_expl1.xyz', destination + '/geom_expl1.xyz')
                    shutil.copy(source + '/geom1.xyz', destination + '/geom1.xyz')
                    
                    # Now, generate the new fortinput file.
                    try:
                        fortinp(dir)
                    except:
                        gsmutil.gsmend('Error while creating input file (fortinput or converg.data) for QoMMMa, check for necessary input and its format, see manual', usrdir)
                                        
                    # Reparameterise the given node.
                    try:
                        os.system(qommmadir + '/bin/gsm_reparam.x')
                    except:
                        gsmutil.gsmend('Could not reparameterise the string, as gsm_reparam.x failed.', dir, usrdir) 
                    
                    # Lastly, copy over the new geometry file and clean-up.
                    source = dirs
                    destination = dir
                    shutil.copy(source + '/geom1.xyz', destination + '/init_geom1.xyz')

        # We're all finished!
        gsmutil.gsmlog('The growing string method is complete!', usrdir)      

def QoMMMa_SE_gsm():
    """
    
    // Function which runs the SE_GSM workflow. //
    
    Arguments
    ----------

    """    

    # The time that the GSM run starts is written to the GSM log file.
    hlp = asctime()
    hlp1 = platform.uname()
    gsmutil.gsmlog('GSM QoMMMa job starts at   ' + str(hlp) + '\n' + '  and all calculations are running in node  ' + str(hlp1), usrdir)

    # Details of GSM job.
    gsmutil.gsmlog('You have selected to use GSM type: ' + str(gsmtype), usrdir)
    if (len(driving_coords) == 0):
        gsmutil.gsmend('FATAL ERROR: for a single-ended GSM calculation, you must provide at least 1 driving coordinate in qommma.in.', usrdir)
    else:
        gsmutil.gsmlog('With the GSM type: ' + str(gsmtype) + '   the driving coordinate(s) are as below...', usrdir)
        for driv in driving_coords:
            gsmutil.gsmlog(str(driv), usrdir)           

    # The growth phase commences.
    gsmutil.gsmlog('        /// Growing the string... ///       ', usrdir)
    all_nodes = []
    is_grown = False
    while is_grown == False:

        # Initialising reactant node directory.
        nodeR_dir = usrdir + '/nodeR'
        if (os.path.exists(nodeR_dir)) is not True:
            gsmutil.gsmend('FATAL ERROR: for a single-ended GSM calculation, a reactant directory called ''nodeR'' within the working directory must exist.', usrdir)
            
        # Initialising primitive internal coordinates for the reaction pathway so that they can be written to fortinput.
        # It is essential that the primitive internal coordinates remain the same across the pathway.
        prims_r = nodeR_dir + '/jobfiles/prim_list'
        if (os.path.exists(prims_r)) is not True:
            gsmutil.gsmend('FATAL ERROR: could not find the primitive definitions in /nodeR/jobfiles. The reactant node must first be optimised with delocalised internal coordinates.', usrdir)
        global prim_num, prim_defs
        prim_num, prim_defs = gsmutil.SE_initialise_prims(prims_r, usrdir) 
    
        # If any additional primitive internal coordinates were specified, then these are added to prim_defs now.
        global additional_prims
        additional_prims=[]
        prim_num, prim_defs = gsmutil.additional_prims(additional_prims, prim_num, prim_defs, usrdir) 

        # For the first cycle, the frontier node is the reactant node.
        frontier_dir = nodeR_dir
        all_nodes.append(nodeR_dir)
        
        # Now begin the loop of adding nodes.
        # When a stationary point (i.e., a barrier) is found, the growing is terminated.
        current_nodes = 1    
        SP_found = False
        while SP_found == False:

            # Creating directory for new node.
            new_frontier_dir = usrdir + '/node' + str(current_nodes + 1)
            current_nodes += 1
            all_nodes.append(new_frontier_dir)
            os.mkdir(new_frontier_dir)
            gsmutil.gsmlog('Node ' + str(current_nodes) + '   has been added to the string...', usrdir)

            # To define the coordinates of the new nodes, create tangent along the user-specified driving coordinates.
            # New nodes are then added along this tangent.
            tangent = gsmutil.SE_get_tangent(frontier_dir, driving_coords, usrdir)
            
            # Now create a new qommma.in file for the new node and copy the geometry of the previous node over.
            gsmutil.SE_add_node(frontier_dir, new_frontier_dir, tangent, driving_coords, usrdir)
            
            # Update the frontier directory to represent the new node.
            old_frontier_dir = frontier_dir
            frontier_dir = new_frontier_dir

            # Now, the new node can be optimised subject to the conditions of a growth phase calculation.
            # To give correct input, read in the new qommma.in.
            inpf = frontier_dir + '/qommma.in'
            try:
                exec(open(inpf).read(), globals())
            except:
                gsmutil.gsmend('FATAL ERROR: could not open input file for new node: ' + inpf + '  This may be due to error in user input file, or maybe it was not generated? Check for necessary inputs and its formats, see manual', usrdir)
            gsmutil.gsmlog('Optimising node (growth-phase): ' + inpf + '........', usrdir) 
            QoMMMa_opt(frontier_dir) # performing QoMMMa growth-phase optimisation.

            # Lastly, check if the current node is lower in energy than the previous node.
            # If so, then the stationary point (SP) has been found and the main part of the growth phase is over.
            SP_found = gsmutil.SE_check_delE(old_frontier_dir, frontier_dir, current_nodes)

            # If the number of nodes becomes greater than the maximum allowed number of nodes, then the calculation is ended.
            if ((current_nodes == max_nodes) and (SP_found == False)):
                gsmutil.gsmlog('Please note: the maximum number of allowed nodes is: ' + str(max_nodes) + ' , and this has been reached - adding final nodes...', usrdir)
                break

        # When the stationary point (SP) has been surpassed, 2 - 4 new nodes with large steps are created - hopefully finding the minimum on the other side.
        # The large steps are defined based on the total primitive internal coordinate change from the pre-SP nodes.
        # This is based on the assumption that the drop off in energy after the SP will be approximately at the same rate as pre-SP.
        tangent, total_new = gsmutil.SE_get_final_tangent(current_nodes, driving_coords, usrdir)
        gsmutil.gsmlog('Based on the magnitude of change in primitive internal coordinates, ' + str(total_new) + ' new nodes will be added...', usrdir)
        
        # Based on the magnitude of the tangent, new directories can be created. 
        # The criteria for these magnitudes are essentially based on trial and error, but can be easily modified.
        new_dirs = []
        for i in range(1,(total_new+1)):
            new_frontier_dir = usrdir + '/node' + str(current_nodes + i)
            new_dirs.append(new_frontier_dir)
            all_nodes.append(new_frontier_dir)
            os.mkdir(new_frontier_dir)
            gsmutil.gsmlog('Node ' + str(current_nodes + i) + '   has been added to the string...', usrdir)
        current_nodes += total_new  
        
        # Now create new qommma.in files for the new nodes and copy the geometries of the previous nodes over.
        # In this case, the starting geometry for all 2-4 nodes is the same and is the current frontier node.
        # The generation of the new geometries (mathematically intensive) is handled by the Fortran code.
        gsmutil.SE_add_final_nodes(frontier_dir, new_dirs, tangent, driving_coords, usrdir)
        
        # Copy the geometry from the previous frontier node to the new one.
        geom_s = frontier_dir + '/geom1.xyz'
        geom_d = new_dirs[0] + '/init_geom1.xyz'
        shutil.copy(geom_s, geom_d) 
        
        # Update the frontier directories to represent the new nodes.
        frontier_dirs = []
        frontier_dirs = new_dirs        
        # Lastly, the new nodes can be optimised subject to the conditions of a growth phase calculation.
        # To give correct input, read in the new qommma.in.
        for counter,dir in enumerate(frontier_dirs):
            inpf = dir + '/qommma.in'
            try:
                exec(open(inpf).read(), globals())
            except:
                gsmutil.gsmend('Could not open input file for node: ' + inpf + '  This may be due to error in user input file, or maybe it was not generated? Check for necessary inputs and its formats, see manual', usrdir)
            gsmutil.gsmlog('Optimising node (growth-phase): ' + inpf + '........', usrdir)
            QoMMMa_opt(dir) # performing QoMMMa growth-phase optimisation.
            
            # Now copy the geometry to the next new frontier node...
            # Unless it is the last node, of course, as this is the final frontier node and thus the newly obtained product node.
            if ((counter+1) != total_new):
                source = dir
                geom_s = source + '/geom1.xyz'
                destination = frontier_dirs[counter+1]
                geom_d = destination + '/init_geom1.xyz'  
                shutil.copy(geom_s, geom_d)  
    
        # Growth phase complete!
        gsmutil.gsmlog('The growth phase is complete! There is now a total (including reactant) of: ' + str(current_nodes) + '   nodes...', usrdir)
        is_grown = True 
    
    # Now move to the optimisation phase.
    gsmutil.gsmlog('        /// Optimising the string... ///       ', usrdir)
    is_optimised = False
    reparam_counter = 0
    while is_optimised == False:

        # First, generate a list of all the tangents for ease of working.
        # Tangents are defined between neighbouring nodes and the string is optimised with stricter convergence criteria than in the growth phase.
        tangent_list, tangent_prims_list = gsmutil.get_tangents_opt(all_nodes, usrdir, driving_coords)
        
        # Now, optimise each node with stricter convergence criteria than in the growth phase.
        # Create the qommma.in files for every node and then optimise each.
        gsmutil.gen_input_opt(all_nodes, tangent_list, tangent_prims_list, usrdir)
        for counter,dir in enumerate(all_nodes):
            if not ("nodeR" or "nodeP") in dir:
                inpf = dir + '/qommma.in'
                try:
                    exec(open(inpf).read(), globals())
                except:
                    gsmutil.gsmend('FATAL ERROR: Could not open input file for node: ' + inpf + '  This may be due to error in user input file, or maybe it was not generated? Check for necessary inputs and its formats, see manual', usrdir)
                gsmutil.gsmlog('Optimising node (optimisation-phase): ' + inpf + '........', usrdir)
                QoMMMa_opt(dir) # performing QoMMMa optimisation-phase optimisation.
         
        # Now, check for convergence of all nodes.
        # If all are converged, then the optimisation-phase is complete.
        # Otherwise, reparameterise the string and optimise all nodes again.
        is_converged = gsmutil.check_convergence(all_nodes, current_nodes, usrdir)

        if (all(is_converged) and (reparam_counter > 0)) is True: # we should at least reparameterise once due to nodes not being evenly spaced...
            gsmutil.gsmlog('The optimisation phase is complete!', usrdir)
            is_optimised = True
        else:
            reparam_counter += 1
            gsmutil.gsmlog('Performing reparameterisation number: ' + str(reparam_counter) + ', stand by...', usrdir)

            # First, generate the new qommma.in files.
            gsmutil.reparam_opt(all_nodes, current_nodes, usrdir)

            # Now, reparameterise every node.
            for counter,dir in enumerate(all_nodes):
                if not ("nodeR" or "nodeP") in dir:
                    # First read in variables from default.in and qommma.in.
                    inpf = dir + '/qommma.in'
                    try:
                        exec(open(inpf).read(), globals())
                    except:
                        gsmutil.gsmend('Could not open input file for node: ' + inpf + '  This may be due to error in user input file, or maybe it was not generated? Check for necessary inputs and its formats, see manual', usrdir)
                        
                    # The directory 'qmmm_*' is created at TMPDIR. Typically, this will be a scratch directory on a node.
                    try:
                        dirs = os.environ['TMPDIR'] + '/qmmm_'
                        dirs = dirs + str(os.getpid()) + '_' + str(counter)
                        if os.path.exists(dirs):
                            shutil.rmtree(dirs)
                        os.mkdir(dirs)
                    except:
                        gsmutil.gsmlog('Either Variable TMPDIR is not set or Could not be able to create temporary directory ' + dirs, usrdir)
                        sys.exit()
                    os.chdir(dirs)
                    
                    # Copy over the current geometries...
                    source = dir + '/jobfiles'
                    destination = dirs
                    shutil.copy(source + '/geom_expl1.xyz', destination + '/geom_expl1.xyz')
                    shutil.copy(source + '/geom1.xyz', destination + '/geom1.xyz')
                    
                    # Now, generate the new fortinput file.
                    try:
                        fortinp(dir)
                    except:
                        gsmutil.gsmend('Error while creating input file (fortinput or converg.data) for QoMMMa, check for necessary input and its format, see manual', usrdir)
                                        
                    # Reparameterise the given node.
                    try:
                        os.system(qommmadir + '/bin/gsm_reparam.x')
                    except:
                        gsmutil.gsmend('Could not reparameterise the string, as gsm_reparam.x failed.', dir, usrdir) 
                    
                    # Lastly, copy over the new geometry file and clean-up.
                    source = dirs
                    destination = dir
                    shutil.copy(source + '/geom1.xyz', destination + '/init_geom1.xyz')

    # We're all finished!
    gsmutil.gsmlog('The growing string method is complete!', usrdir)

    # Following the growing string method, the user has the option to perform a series of single point energy calculations on the pathway.
    if spe_pathway == 1:
        gsmutil.gsmlog('Now performing the SPE calculations along the pathway.', usrdir)
        gsmutil.gsmlog('The program used for this is ' + spe_qmcode + ', and the level of theory is ' + spe_qmkey, usrdir)
        QoMMMa_pathway_SPE(all_nodes)       

def QoMMMa_scan():
    """
    
    // Function which runs a series of constrained optimsations. //
    
    Arguments
    ----------

    """    

def QoMMMa_pathway_SPE(dirs):
    """
    
    // Function which peforms a series of single point energy calculations on a series of directories. //
    // Typically, this will be used following generation of a reaction pathway with the growing string method or a reaction coordinate scan. //
    // NOTE: for now, this only works with Gaussian, but will be extended to other QM codes in future. //
    
    Arguments
    ----------
    dirs : list
        Contains the paths for all the points on a reaction pathway for which the SPE calculations are to be performed.

    """ 

    # If any directories from previous SPE jobs are found, these are copied and retained.
    dir_spe = os.path.join(usrdir, 'SPE_jobs/')
    if (os.path.isdir(dir_spe)):
         if os.path.exists(usrdir + '/SPE_jobs_'):
                shutil.rmtree(usrdir + '/SPE_jobs_''/' + file + '_')
        shutil.copytree(usrdir + '/SPE_jobs', usrdir + '/SPE_jobs_')
        shutil.rmtree(usrdir + '/SPE_jobs')
        gsmutil.gsmlog('Note, previous SPE jobs directory was moved to SPE_jobs_; before submitting next QoMMMa job if you need this file rename it since it will be deleted', usrdir)
    os.mkdir(dir_spe)

    # Now, copy over the geometry, qommma.in and MM point charge files for all the points on the reaction path.
    for counter,dir in enumerate(dirs):
        shutil.copy(dir + '/jobfiles/image1/qmgeom1.xyz', dir_spe + 'node_' + str(counter+1) + '_qmgeom.xyz')
        shutil.copy(dir + '/qommma.in', dir_spe + 'node_' + str(counter+1) + '_qommma.in') # these should all be basically the same, anyway.
        shutil.copy(dir + '/jobfiles/image1/charges1.xyz', dir_spe + 'node_' + str(counter+1) + '_charges.pc')

    # Add each of the geometry, qommma.in and MM point charge files for each directory to a list.
    os.chdir(dir_spe)
    qommma_files = sorted(glob.glob(f'{os.getcwd()}/*.in'))
    geometry_files = sorted(glob.glob(f'{os.getcwd()}/*.xyz'))
    pointch_files = sorted(glob.glob(f'{os.getcwd()}/*.pc'))

    # For each geometry, qommma.in and MM point charge files, first read in the qommma.in parameters.
    for counter in range(len(qommma_files)):
        try:
            inpf = qommma_files[counter]
            exec(open(inpf).read(), globals())
        except:
            qomutil.qomlog('FATAL ERROR: Could not find or read input file: ' + inpf + '. This may be due to errors in the syntax of the input file, or it may simply be missing.', usrdir)
            sys.exit()

        # Initialise the geometry and point charge file.
        geom = geometry_files[counter]
        pointch = pointch_files[counter]
        qmjob_prefix = 'node_' + str(counter+1) + '_jobfile'

        # Now generate and run the Gaussian job.
        from gauutil import spe_gaumain
        imn=1
        cwd = usrdir
        spe_gaumain(imn, cwd, usrdir, qmjob_prefix, qmgau_job, nqm, nlink, spe_qmkey, cha_mul, extra_basis, gau_head, extra_guess, pointch, geom)

    # For continuity, manually rename the reactant node files to R.
    os.rename(dir_spe + 'node_1_qmgeom.xyz', dir_spe + 'node_R_qmgeom.xyz')
    os.rename(dir_spe + 'node_1_qommma.in', dir_spe + 'node_R_qommma.in')
    os.rename(dir_spe + 'node_1_charges.pc', dir_spe + 'node_R_charges.pc')
    os.rename(dir_spe + 'node_1_jobfile.in', dir_spe + 'node_R_jobfile.in')
    os.rename(dir_spe + 'node_1_jobfile.log', dir_spe + 'node_R_jobfile.log')

################
# MAIN PROGRAM #
################
if __name__ == "__main__":
    # The user directory is set.
    usrdir = os.getcwd()
    
    # If any directories from previous GSM runs are found, these are copied and retained.
    # The reactant and product side nodes are NOT copied and retained as these should be obtained and minimised prior to using GSM.
    for file in os.listdir(usrdir):
        dir = os.path.join(usrdir, file)
        if ((os.path.isdir(dir)) and ('node' in file) and ('nodeP' not in file) and ('nodeR' not in file) and ('_' not in file)):
            if os.path.exists(usrdir + '/' + file + '_'):
                shutil.rmtree(usrdir + '/' + file + '_')
            shutil.copytree(usrdir + '/' + file, usrdir + '/' + file + '_')
            shutil.rmtree(usrdir + '/' + file)
            gsmutil.gsmlog('Note, node directories from previous GSM QoMMMa runs have been moved to /nodei_; before submitting next QoMMMa job if you need this directory rename it since it will be deleted', usrdir)
     
    # Files from previous runs of QoMMMa are renamed and retained.
    if os.path.exists(usrdir + '/QoMMMa8.log'):
        shutil.copy(usrdir + '/QoMMMa8.log', usrdir + '/QoMMMa8.log_')
        os.remove(usrdir + '/QoMMMa8.log')
        qomutil.qomlog('Note, log file, QoMMMa8.log of your previous QoMMMa run was moved as QoMMMa8.log_; before submitting next QoMMMa job if you need this file rename it since it will be deleted', usrdir)

    # The directory for the program is set.
    try:
        qommmadir = os.environ['QOMMMA']
    except:
        qomutil.qomlog('FATAL ERROR: environment variable QOMMMA is not set. Set QOMMMA to the base directory of QoMMMa', usrdir)
        sys.exit()
    
    # Default parameters are read such that the parameters are saved as variables.
    try:
        exec(open(qommmadir + '/lib/default.in').read(), globals())
    except:
        qomutil.qomlog('FATAL ERROR: could not open default file: ' + qommmadir + '/lib/default.in', usrdir)
        sys.exit()

    # The user input file is the second command line argument, so it is set as the variable inpf.
    if len(sys.argv) < 2:
        qomutil.qomlog('Usage: ' + sys.argv[0] + ' <file.in>', usrdir)
        sys.exit()
    inpf = sys.argv[1]

    # The input file is read such that the parameters are saved as variables.
    # In this way, any parameters which are chosen by the user overwrite the default parameters.
    try:
        exec(open(inpf).read(), globals())
    except:
        qomutil.qomlog('FATAL ERROR: Could not find input file: ' + inpf + '. This may be due to error in user input file. Check for necessary inputs and its formats, see manual', usrdir)
        sys.exit()

    # Finally, the different options for a QoMMMa run are:
    #  1. double-ended growing string method (DE-GSM)
    #  2. single-ended growing string method (SE-GSM)
    #  3. adiabatic mapping reaction coordinate scan (NOTE: not yet implemented)
    #  4. normal QM/MM optimisation
    if gsmtype.lower() == 'de_gsm':
        QoMMMa_DE_gsm()
    elif gsmtype.lower() == 'se_gsm':
        QoMMMa_SE_gsm()
    elif scantype.lower() == 'adiabatic':
        QoMMMa_scan()
    else:  
        QoMMMa_opt(usrdir)
