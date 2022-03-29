#!/usr/bin/python3

"""
// This is the main QoMMMa Python file, which collects the functions defined in other Python files for QoMMMa execution. //
"""

# Global imports.
import os
import sys
import shutil
import qomutil
import platform
from time import asctime

def fortinp():
    """
    
    // Function which creates two input files (fortinput and converg.data) required for execution of the QoMMMa Fortran code. //
    // This function takes no arguments as they have been read in through default and user input files. //
    
    Arguments
    ----------
    NONE

    """
  
    # The file fortinput is created.
    fd = open('fortinput','w')
  
    # Numbers relating to QM, MM, link, and total QM optimised atoms are written to fortinput.
    fd.write('natom     nqm      nlink   nopt')
    fd.write('\n')
    fd.write(str(natom))
    fd.write(str(('%s%d'%('   ', nqm))))
    fd.write(str(('%s%d'%('   ', nlink))))
    fd.write(str(('%s%d'%('   ', nopt))))
    fd.write('\n')
    
    # Whether dispersion is to be calculated (1) or not (0) is written to fortinput.
    fd.write('disp')   # dispersion keyword
    fd.write('\n')
    fd.write(str((disp)))   # dispersion = 1 or 0 (on/off)
    fd.write('\n') 

    # The number of constraints, is used, as well as whether the constraint is harmonic (1) or tannic (2) is written to fortinput.
    fd.write('ncon     kcnstype')
    fd.write('\n')
    fd.write(str(ncon))
    if kcnstype.lower() == 'harmonic':
        ikcnstype = 1
    elif kcnstype.lower() == 'tanh':
        ikcnstype = 2
    else:
        qomutil.qomlog('Error, unknown constrain type : ' + kcnstype + ' is requested', usrdir)
    fd.write(str(('%s%d'%('   ', ikcnstype))))
    fd.write('\n') 
      
    # If used, the nudged elastic band method used as well as the spring constant is written to fortinput.
    if nebtype.lower() == 'none':
        inebtyp = 0
    elif nebtype.lower() =='neb_bfgs':
        inebtyp = 1
    elif nebtype.lower() == 'cineb_bfgs':
        inebtyp = 2
    elif nebtype.lower() == 'neb_cg':
        inebtyp = 3
    elif nebtype.lower() == 'none_cg':
        inebtyp = 4
    elif nebtype.lower() =='cineb_cg':
        inebtyp = 5
    else:
        qomutil.qomlog('Error, unknown nebtype : ' + nebtype + ' is requested', usrdir)
    fd.write('nimg    nebtype    kspring')
    fd.write('\n')
    fd.write(str(nimg))
    fd.write(str(('%s%d'%('   ', inebtyp))))
    fd.write(str(('%s%d'%('   ', kspring))))
    fd.write('\n')
      
    # If used, the growing string method is written to fortinput.
    if gsmtype.lower() == 'none':
        igsmtyp = 0
    elif gsmtype.lower() =='de_gsm':
        igsmtyp = 1
    elif gsmtype.lower() =='se_gsm':
        igsmtyp = 2
    else:
        qomutil.qomlog('Error, unknown gsmtype : ' + gsmtype + ' is requested', usrdir)
    qomutil.qomlog('The gsm type : ' + gsmtype + ' is requested', usrdir)
    fd.write('gsmtype')
    fd.write('\n')
    fd.write(str(igsmtyp))
    fd.write('\n')

    # The type of coordinates used is written to fortinput.
    if coordtype.lower() == 'cart':
        icoordtyp = 0
    elif coordtype.lower() =='dlc':
        icoordtyp = 1
    else:
        qomutil.qomlog('Error, unknown coordtype : ' + coordtype + ' is requested', usrdir)
    qomutil.qomlog('The coordinate type : ' + coordtype + ' is requested', usrdir)
    fd.write('coordtype')
    fd.write('\n')
    fd.write(str(icoordtyp))
    fd.write('\n') 
    
    # If DLC are used, then the type of primitive internal coordinates used to generate them is written to fortinput.
    if primtype.lower() == 'tc':
        iprimtyp = 0
    elif primtype.lower() =='full':
        iprimtyp = 1
    else:
        qomutil.qomlog('Error, unknown primtype : ' + primtype + ' is requested', usrdir)
    qomutil.qomlog('The coordinate type : ' + primtype + ' is requested', usrdir)
    fd.write('primtype')
    fd.write('\n')
    fd.write(str(iprimtyp))
    fd.write('\n') 
    
    # The atom indices and types relating to the QM atoms is written to fortinput using the function qmread defined in qomutil.py.
    # If the old format of defining qm atoms is used, then it will still work, but this will be removed in future versions.
    if nqm > 0:
        fd.write('QM atoms list')
        fd.write('\n')
        if qm_lst != 'None':
            fd.close()		
            qomutil.qmread(qm_lst, usrdir)
            fd = open('fortinput', 'a')
        else:
            fd.write((qm.lstrip()))		
            qomutil.qomlog('''Note, to provide QM atom list try to use 'qm_lst' option (see manual), instead of 'qm'. Since 'qm' option will be removed in the next version of QoMMMa.''', usrdir)
    else:
        qomutil.qomend('ERROR: QM atoms need to be defined for a QM/MM calculation.', cwd, usrdir)
 
    # Link atom details are written to fortinput using the function link_write defined in qomutil.py.    
    # If the old format of defining link atoms is used, then it will still work, but this will be removed in future versions.
    if nlink > 0:
        fd.write('Link atom details')
        fd.write('\n')
        if link_lst != 'None':
            fd.close()
            qomutil.link_write(link_lst) 
            fd = open('fortinput', 'a')
        else:
            fd.write((link.lstrip()))
            qomutil.qomlog('''Note, to provide link atom details try to use 'link_lst' option (see manual), instead of 'link'. Since 'link' option will be removed in the next version of QoMMMa.''', usrdir)
    else:
        qomutil.qomlog('No Link atom details is taken, since nlink=0.', usrdir) 
    
    # If any, the hessian optimised MM atoms are written to fortinput using the function hesopt_write defined in qomutil.py.
    if (int(nopt) - int(nqm) - int(nlink)) > 0:	
        fd.write('Hessian optimized MM atoms list')	
        fd.write('\n')
        fd.close()
        qomutil.hesopt_write(hesoptmm, usrdir)
        fd = open('fortinput', 'a')
    else:
        qomutil.qomlog('No MM atom is included in Hessian optimization',usrdir)
    
    # If used, constraint details are written to fortinput using the function cons_write defined in qomutil.py.
    # If the old format of defining constraints is used, then it will still work, but this will be removed in future versions.
    if ncon > 0:
        fd.write('Constrain details')
        fd.write('\n')
        if constrain_lst != 'None':
            fd.close()
            qomutil.cons_write(constrain_lst, usrdir)  
            fd = open('fortinput', 'a')
        else:
            fd.write((constrain.lstrip()))     
            qomutil.qomlog('''Note, to provide constrain details try to use 'constrain_lst' option (see manual), instead of 'constrain'. Since 'constrain' option will be removed in the next version of QoMMMa.''', usrdir)
    else:
        qomutil.qomlog('No constrain is taken into account, since ncon=0.', usrdir)
    
    # The charge dispersion in the MM region required by the inclusion of link atoms is written to fortinput using the function newcha_write defined in qomutil.py.
    # If the old format of defining modified charges is used, then it will still work, but this will be removed in future versions.
    fd.write('Modified charge for MM atoms. First, number of such atoms, and then atom label and charge')
    fd.write('\n')
    if newcha_lst != 'None':
         fd.close()
         qomutil.newcha_write(newcha_lst)
         fd = open('fortinput', 'a')
    else:
         fd.write(newcha.lstrip()) 
         if int(newcha.lstrip()[0]) > 0:          
             qomutil.qomlog('''Note, to provide new charge details try to use 'newcha_lst' option (see manual), instead of 'newcha'. Since 'newcha' option will be removed in the next version of QoMMMa.''', usrdir)
    
    # If any, the inactive MM atoms are written to fortinput using the function inprange defined in qomutil.py.
    fd.write('Inactive atoms. First, number of such atoms then list')
    fd.write('\n')
    kk = ninact.split()
    fd.write(str(kk[0]))		       
    if int(kk[0]) != 0:
       fd.write('\n')
       kk = kk[1:]
       fd.close()
       qomutil.inprange(kk, 'Inactive', usrdir)
    else: 
       qomutil.qomlog('No atom is set as inactive', usrdir)
       fd.close()
   
    # The file converg.data is created and set to filename variable 'fc'.  
    # All convergence criteria are written to converg.data as it is used in later Fortran calculations.
    
    #########################################
    # GSM CONVERGENCE CRITERIA WILL GO HERE #
    #########################################
    
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

def mmjob(mjob):
    """
    
    // Function which performs a series of operations for the first and second Tinker jobs. //
    // These jobs are performed separately from the main QoMMMa cycle as they do not fit exactly into the framework of a typical cycle. //
    
    Arguments
    ----------
    mjob : string
        If mjob is 'initial', then it executes the initial Tinker calculation, and if 'initial2', then it prepares for subsequent QoMMMa cycles.

    """
    
    # At present, QoMMMa operates only with Tinker, so the variable mmcode must be 'Tinker'.
    # This function runs either the initial MM job (mjob = 'initial') or prepares for the second MM job (mmjob = 'initial2') using results from the Fortran QoMMMa setup program.
    # Running these independently of the main QM/MM cycle is important for correct setup.
    if mmcode == 'Tinker':
        # The initial MM job is setup and subsequently carried out.
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
                qomutil.qomlog( 'Note, default force field parameter file is charmm27.prm, this could be changed or modified through mmkey in user input', usrdir)           
            except: 
                qomutil.qomlog('Note, no mmkey is given by user, default force field parameter file, charmm27.prm will be used here', usrdir) 
            fdt.close()
           
            # Important details for the Tinker job are written to geom.key.
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
                     old = usrdir + '/' + fin         
                     shutil.copy(old, fin) 
                     shutil.copy(fin, fdn)
                except:
                    qomutil.qomend('Error, input file : ' + fin + ' is not found in user directory', cwd, usrdir)
                try:
                    qomutil.geom_expl(fdn, l)
                except:
                    qomutil.qomend('ERROR, while creating geom_expl*.xyz for image: ' + str(l), cwd, usrdir)   
                    
                # To obtain the charges for the upcoming QM job, the Tinker program analyze_grad is run for each image using the function tinkercharge defined in qomutil.py.
                try:
                    qomutil.tinkercharge(qommmadir) 
                except:
                    qomutil.qomend('Error while getting initial charges using Tinker program analyze_grad', cwd, usrdir) 
                    
        # The second MM job is setup for actual QM/MM cycles.      
        elif mjob == 'initial2':   
           try:
               # The details in tinker_header and inactive_list are written to geom.key
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
               qomutil.qomend('Error: while preparing geom.key file after QoMMMa initial run. This step is used to include inactive atom list (including QM atoms) in tinker key file', cwd, usrdir)
        else:
            qomutil.qomend('Unknown MM job: ' + job, cwd, usrdir)
    else:
        qomutil.qomend('Unknown MM code: ' + mmcode, cwd, usrdir)


def qminitial():
    """
    
    // Function which creates the initial inputs for a QM job. //
    // This function takes no arguments as they have been read in through default and user input files. //
    
    Arguments
    ----------
    NONE

    """
    
    # Jaguar job is prepared by using functions from jagutil.py.
    # Different operations are performed if doing a minimum energy crossing point calculation.
    if qmcode.lower() == 'jaguar':
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
        if os.path.exists(usrdir + '%s%d'%('/orbs', imn)):
            shutil.copy(usrdir + '%s%d'%('/orbs', imn), dst)
            qomutil.qomlog('initial orbital guess for Jaguar job is copied from user directory for image :' + str(imn), usrdir)
        if job.lower()=='mecp':
            if os.path.exists(usrdir + '%s%d'%('/orbsA', imn)):
                shutil.copy(usrdir + '%s%d'%('/orbsA', imn), dst)
                qomutil.qomlog('initial orbital guess for Jaguar job is copied from user directory for first state of image :' + str(imn), usrdir)
            if os.path.exists(usrdir + '%s%d'%('/orbsB', imn)):
                shutil.copy(usrdir + '%s%d'%('/orbsB', imn), dst)
                qomutil.qomlog('initial orbital guess for Jaguar job is copied from user directory for second state of image :' + str(imn), usrdir)
                
    # Molpro job is prepared.
    # Different operations are performed if doing a minimum energy crossing point calculation.              
    elif qmcode.lower() == 'molpro':
        # The file qmheader is made as it is used in operation of a Molpro job.
        fds = open('qmheader', 'w')
        fds.write(qmmol_header.lstrip())
        fds.close()
        
        # Minimum energy crossing point Molpro job is prepared.
        if job.lower() == 'mecp':
            # The files qmmol_footer_A and qmmol_footer_B are made as they are used in operation of a Molpro job.
            fda = open('qmmol_footer_A', 'w')
            fdb = open('qmmol_footer_B', 'w')
            if qmkey.strip() != 'None':
                fda.write(qmkey.lstrip())
                fdb.write(qmkey.lstrip())
            else:
                qomutil.qomend('''For Molpro, you have to give atleast one input in 'qmkey' for instance 'hf' is must and may be basis set option''', cwd, usrdir) 
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
                if os.path.exists(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.intg')):
                    shutil.copy(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.intg'), dst)
                    qomutil.qomlog('initial orbital guess file is copied from user directory for first mecp state of image :' + str(imn), usrdir)
                if os.path.exists(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.intg')):
                    shutil.copy(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.intg'), dst)
                    qomutil.qomlog('initial orbital guess file is copied from user directory for second mecp state of image :' + str(imn), usrdir)
                    
        # Normal Molpro job is prepared.
        else:
            # The file qmmol_footer is made as it is used in operation of a Molpro job.
            fds = open('qmmol_footer', 'w')
            if qmkey.strip() != 'None':
                fds.write(qmkey.lstrip())
            else:
                qomutil.qomend('''For Molpro, you have to give atleast one input in 'qmkey' for instance 'hf' is must''', cwd, usrdir) 
            fds.write(mol_footer.lstrip())
            fds.close()
            
            # For every image (if using the nudged elastic band, that is), appropriate files are copied to run the Molpro jobs. 
            imn = 0 
            for i in range(nimg):
                imn = imn + 1
                if os.path.exists(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.intg')):
                    dst = cwd + ('%s%d'%('/image', imn))
                    shutil.copy(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.intg'), dst)
                    qomutil.qomlog('initial orbital guess file is copied from user directory for image :' + str(imn), usrdir)
                    
    # Gaussian job is prepared.               
    # Different operations are performed if doing a minimum energy crossing point calculation.   
    elif qmcode.lower() == 'gaussian':
        if qmkey.strip() == 'None':
            qomutil.qomlog( 'Note, no extra QM option is given by user through qmkey, so Gaussian job will run at HF/sto-3G level', usrdir)
        # For every image (if using the nudged elastic band, that is), appropriate files are copied to run the Gaussian jobs.     
        imn = 0 
        for i in range(nimg):
            imn = imn + 1
            dst = cwd + ('%s%d'%('/image', imn))
            if os.path.exists(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.chk')):
                shutil.copy(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.chk'), dst)
                qomutil.qomlog('initial Gaussian check file was copied from user directory for image :' + str(imn), usrdir)
            if job.lower() ==' mecp':
                if os.path.exists(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.chk')):
                    shutil.copy(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.chk'), dst)
                    qomutil.qomlog( 'initial Gaussian check file was taken from user directory for MECP state A of image :' + str(imn), usrdir) 
                if os.path.exists(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.chk')):
                    shutil.copy(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.chk'), dst)
                    qomutil.qomlog('initial Gaussian check file was taken from user directory for MECP state B of image :' + str(imn), usrdir) 
             
    # Orca job is prepared.               
    # Different operations are performed if doing a minimum energy crossing point calculation.             
    elif qmcode.lower() == 'orca':
        if qmkey.strip() == 'None':
            qomutil.qomlog( 'Note, no extra QM option is given by user through qmkey, so ORCA job will run at HF/sto-3G level', usrdir)
        # For every image (if using the nudged elastic band, that is), appropriate files are copied to run the Orca jobs. 
        imn = 0 
        for i in range(nimg):
            imn = imn + 1
            dst = cwd + ('%s%d'%('/image', imn))
            if os.path.exists(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.gbw')):
                shutil.copy(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '.gbw'), dst)
                qomutil.qomlog('initial ORCA gbw file was copied from user directory for image :' + str(imn), usrdir)
            if job.lower() == 'mecp':
                if os.path.exists(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.gbw')):
                    shutil.copy(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_A.gbw'), dst)
                    qomutil.qomlog( 'initial ORCA gbw file was taken from user directory for MECP state A of image :' + str(imn), usrdir) 
                if os.path.exists(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.gbw')):
                    shutil.copy(usrdir + '%s%s%d%s'%('/', qmjob_prefix, imn, '_B.gbw'), dst)
                    qomutil.qomlog('initial ORCA gbw file was taken from user directory for MECP state B of image :' + str(imn), usrdir) 
    else:
        qomutil.qomend('Unknown QM code: ' + qmcode, cwd, usrdir) 


def qmmm():
    """
    
    // Function which performs a QoMMMa cycle. The steps of a QoMMMa cycle are as follows... //
        1. Execute Tinker's analyze_grad to obtain MM gradient.
        2. Prepare input for QM program chosen, and subsequently run a QM job and extracts results.
        3. Using Fortran code, perform a QoMMMa optimisation.
        4a. Check convergence, and if achieved, move program files to user directory and end the QoMMMa program.
        4b. Check convergence, and if not achieved, start MM microiterative optimisation will follow, and create/rearrange files and return to 1.
        5. If the option of a frequency calculation is selected, start the frequency calculation.
        
    // This cycle continues until convergence is achieved or the maximum number of cycles has been reached. //
    // This function takes no arguments as they have been read in through default and user input files. //
    
    Arguments
    ----------
    NONE

    """
      
    # For every image (if using the nudged elastic band, that is), a QoMMMa cycle is performed. 
    imn = 0
    for irm in range(nimg):
        imn = imn + 1
        if os.path.exists(cwd + ('%s%d'%('/update_geom', imn))):
            # The Tinker program analyze_grad is used to obtain the MM gradient and it is copied to the appropriate location.
            try:
                fch = os.popen(qommmadir + '/bin/analyze_grad', 'w')	
                fn = ('%s%d%s'%('geom', imn, '.xyz'))
                fch.write(fn + '\nE\n')
                fch.close()
                shutil.copy('mm_grad', ('%s%d'%('mm_grad', imn)))
                os.remove('mm_grad')
            except:
                qomutil.qomend('Problem with MM gradient extraction at step :' + str(cln), cwd, usrdir)  
              
            # The QM job is performed using the program of the user's choice.
            os.chdir(('%s%d'%('image', imn)))
            qomutil.qomlog('About to run ' + qmcode +  ' job at cycle : ' + str(cln) +  '  for image : ' + str(imn), usrdir)
              
            # Jaguar QM job.
            if qmcode.lower() == 'jaguar':
                if job.lower() != 'mecp':
                    from jagutil import qmjagmain
                    qmjagmain(imn, cwd, usrdir, qmjob_prefix, qmjag_job, nqm, nlink, cln)
                    qmc_job = qmjag_job    # Used later in frequency call.
                elif job.lower() == 'mecp':
                    from jagutil import mecp_jagmain
                    mecp_jagmain(imn, cwd, usrdir, qmjob_prefix, qmjag_job, nqm, nlink, cln)
                      
            # Molpro QM job.       
            elif qmcode.lower() == 'molpro':
                if job.lower() != 'mecp':
                    from molutil import qmmolmain
                    qmmolmain(imn, cwd, usrdir, qmjob_prefix, qmmol_job, nqm, nlink, cln)
                    qmc_job = qmmol_job    # Used later in frequency call. 
                elif job.lower() == 'mecp':
                    from molutil import mecp_molmain
                    mecp_molmain(imn, cwd, usrdir, qmjob_prefix, qmmol_job, nqm, nlink, cln)
                    
            # Gaussian QM job.       
            elif qmcode.lower() == 'gaussian':
                if job.lower() != 'mecp':
                    from gauutil import qmgaumain
                    qmgaumain(imn, cwd, usrdir, qmjob_prefix, qmgau_job, nqm, nlink, cln, qmkey, cha_mul, extra_basis, gau_head, extra_guess)
                    qmc_job = qmgau_job    # Used later in frequency call. 
                elif job.lower() == 'mecp':
                    from gauutil import mecp_gaumain
                    mecp_gaumain(imn, cwd, usrdir, qmjob_prefix, qmgau_job, nqm, nlink, cln, qmkey, cha_mul1, cha_mul2, extra_basis, gau_head, extra_guess)
                      
            # Orca QM job.         
            elif qmcode.lower() == 'orca':
                if job.lower() != 'mecp':
                    from orcautil import qmorcamain
                    qmorcamain(imn, cwd, usrdir, qmjob_prefix, nqm, nlink, cln, qmkey, cha_mul, extra_basis, orca_head)
                    qmc_job = qmorca_job    # Used later in frequency call. 
                elif job.lower() == 'mecp':
                    from orcautil import mecp_orcamain
                    mecp_orcamain(imn, cwd, usrdir, qmjob_prefix, nqm, nlink, cln, qmkey, cha_mul1, cha_mul2, extra_basis, orca_head)
          
            # From geom*.xyz, the file geom_expl*.xyz is created.
            os.chdir(cwd)
            fi = '%s%d%s'%('geom', imn, '.xyz')       
            try:
                qomutil.geom_expl(fi, imn) 	        
            except:
                qomutil.qomend('ERROR, while creating geom_expl*.xyz for image: ' + str(l), cwd, usrdir)    
                os.remove(('%s%d'%('update_geom', imn)))
        else: pass
    
    # Starting QM/MM optimisation.
    # This is the point where the Fortran code is used, so it should be referenced to understand the optimisation process.
    qomutil.qomlog('About to run QM/MM optimization cycle :' + str(cln), usrdir)
    try:		
        if job.lower() == 'mecp':						
            os.system(qommmadir + '/bin/qommma_mecp.x')
        else:
            os.system(qommmadir + '/bin/qommma8.x')  		
    except:
        qomutil.qomend('Error while running QoMMMa program at cycle :' + str(cln), cwd, usrdir)
        
    # Updates the report file with the results from the QoMMMa optimisation, and if convergence has been achieved, then a dummy file 'convergence_ok' is made.
    try: 
        qomutil.qomreport(nimg, usrdir, cwd, cln, qomout)
    except:
        qomutil.qomend('Error while looking for convergence after QoMMMa optimization at cycle :' + str(cln) + ' check QM, MM results or checkfile', cwd, usrdir)
        
    # If the file 'convergence_ok' exists and hence convergence has been achieved, then the program is ended.
    # If a frequency calculation is desired, then it is performed prior to ending the program.
    if os.path.exists(cwd + '/convergence_ok'):
        # For every image (if using the nudged elastic band, that is), a frequency calculation is performed.
        im = 0
        os.remove(cwd + '/convergence_ok')
        for ir in range(nimg):
            im = im + 1
            shutil.copy(cwd + ('%s%d%s'%('/geom', im, '.xyz')), usrdir + ('%s%d%s'%('/geom', im, '.xyz')))
        # Frequency calculation starts.
        if job.lower() == 'freq':
            hlp = asctime()  
            qomutil.qomlog('About to start QoMMMa Frequency calculation at  ' + str(hlp), usrdir) 
            from qomfreq import freqmain
            freqmain(qommmadir, qmkey, qmc_job, usrdir, cha_mul, extra_basis, gau_head, prj_freq, qmjob_prefix, qmcode) 
            
        # QoMMMa optimisation complete!
        qomutil.qomend('Congratulations,', cwd, usrdir)
            
    try:   			# rearranging files 
        qomutil.rearrange(nimg, cwd, qmjob_prefix, cln)	 
    except:
        qomutil.qomend('Error while rearranging the files after QoMMMa optimization at cycle :' + str(cln), cwd, usrdir)
    
    # For every image (if using the nudged elastic band, that is), a microiteration is performed.
    im = 0
    for ir in range(nimg):				 
        try:
            im = im + 1
            if os.path.exists(cwd + ('%s%d'%('/update_geom', im))):
                dst = cwd + ('%s%d'%('/image', im))
                qomutil.qomlog('About to run Microoptimization of geometry :' + str(im), usrdir)
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
                    qomutil.qomend('Error while running Tinker minimization job at cyle :' + str(cln), cwd, usrdir)
                 
                # Using microgeom.xyz_2 created by Tinker's minimize, geom_expl.xyz is created.
                os.remove('gradcorrection')
                fi = '%s%d%s'%('microgeom', im, '.xyz_2')
                try:
                    qomutil.geom_expl(fi, im) 		  	 
                except:
                    qomutil.qomend('ERROR, while creating geom_expl*.xyz for image: ' + str(l), cwd, usrdir)    
                shutil.copy(('%s%d%s'%('geom_expl', im, '.xyz')), cwd + ('%s%d%s'%('/geom_expl',im, '.xyz')))
                os.remove(('%s%d%s'%('microgeom', im, '.xyz_2')))
                os.remove('%s%d%s'%('geom_expl', im, '.xyz'))
                os.chdir(cwd)
            else:
                qomutil.qomlog('image Geometry  ' + str(im) + '  is fixed in this iteration', usrdir)
        except:
             qomutil.qomend('Error while doing micro-iteration at QM/MM cycle :' + str(cln), cwd, usrdir)
             
    # Using the Fortran code, the post microiteration manipulation of files is performed.          
    try:
        os.system(qommmadir + '/bin/postmicro8.x') 
    except:
        qomutil.qomend('Problem in QoMMMa post-processing job at cycle :' + str(cln), cwd, usrdir)
    
    # Files are rearranged using the function postrearrange defined in qomutil.py.
    try:
        qomutil.postrearrange(nimg, cwd)
    except:
        qomutil.qomend('Error, while rearranging the files after QoMMMa post-processing  job at cycle :' + str(cln), cwd, usrdir)  
            

################
# MAIN PROGRAM #
################
if __name__ == "__main__":
    # User directory is set.
    usrdir = os.getcwd()
    
    # Files from previous runs of QoMMMa are renamed and retained.
    if os.path.exists(usrdir + '/QoMMMa8.log'):
        shutil.copy(usrdir + '/QoMMMa8.log', usrdir + '/QoMMMa8.log_')
        os.remove(usrdir + '/QoMMMa8.log')
        qomutil.qomlog('Note, log file, QoMMMa8.log of your previous QoMMMa run was moved as QoMMMa8.log_; before submitting next QoMMMa job if you need this file rename it since it will be deleted', usrdir)
    
    # The time that the job starts is written to the log file.
    hlp = asctime()
    hlp1 = platform.uname()
    qomutil.qomlog('QoMMMa job starts at   ' + str(hlp) + '\n' + '  and running in node  ' + str(hlp1), usrdir)
    
    # The base directory is set.
    try:
        qommmadir = os.environ['QOMMMA']
    except:
        qomutil.qomlog('Environment variable QOMMMA is not set. Set QOMMMA to the base directory of QoMMMa', usrdir)
        sys.exit()
    
    # Default parameters are read such that the parameters are saved as variables.
    try:
        exec(open(qommmadir + '/lib/default.in').read())
    except:
        qomutil.qomlog('Could not be able to open default file: ' + qommmadir + '/lib/default.in', usrdir)
        sys.exit()
    
    # The user input file is the second command line argument, so it is set as the variable inpf.
    if len(sys.argv) < 2:
        qomutil.qomlog('Usage: ' + sys.argv[0] + ' <file.in>', usrdir)
        sys.exit()
    inpf = sys.argv[1]
    
    # The input file is read such that the parameters are saved as variables.
    try:
        exec(open(inpf).read())
    except:
        qomutil.qomlog('Could not be able to open input file: ' + inpf + '  This may be due to error in user input file. Check for necessary inputs and its formats, see manual', usrdir)
        sys.exit()
    
    # The directory 'qmmm_*' is created at TMPDIR. Typically, this will be a scratch directory on a node.
    try:
        dirs = os.environ['TMPDIR'] + '/qmmm_'
        dirs = dirs + str(os.getpid())
        os.mkdir(dirs)
    except:
        qomutil.qomlog('Either Variable TMPDIR is not set or Could not be able to create temporary directory ' + dirs, usrdir)
        sys.exit()
    
    # Directories are reset and the current working directory is set.
    os.chdir(dirs)
    cwd = os.getcwd()
    
    # The initial MM job is performed.
    mmjob('initial')
    
    # The primary file which communicates between Python and Fortran, fortinput, is generated.
    try:
        fortinp()		
    except:
        qomutil.qomend('Error while creating input file (fortinput or converg.data) for QoMMMa, check for necessary input and its format, see manual', cwd, usrdir)
    
    # The Fortran based QoMMMa "setup" program is run. This creates building blocks for QM & MM jobs, and initializes output.
    try:
        os.system(qommmadir + '/bin/setup8.x')
    except:
        qomutil.qomend('QoMMMa initial setup program fails', cwd, usrdir) 
    
    # Initial QM input files are prepared.
    qminitial()	
    	 	
    # Results from the QoMMMa setup program are moved.
    qomutil.setini(cwd, nimg, usrdir) 	
    
    # New MM key files are creates using results from the QoMMMa setup program.
    mmjob('initial2')
    
    # Main QM/MM run. The variable 'cln' is the QM/MM cycle number.
    # Performing simply a frequency calculation.
    if maxcycle == 0 and job.lower() == 'freq':
        if len(prj_freq) == 0:
            qomutil.qomlog('Frequency job is started with the geometry given in inpgeom_prefix*.xyz file, here qmcode used is Gaussian', usrdir) 
            from qomfreq import freqmain
            qmc_job = qmgau_job
            freqmain(qommmadir, qmkey, qmc_job, usrdir, cha_mul, extra_basis, gau_head, prj_freq, qmjob_prefix, qmcode) 
        else:
            qomutil.qomend('Error, Projected Frequency calculation is requested; set maxcycle=1 or request full optimization and frequency calculation', cwd, usrdir)
        
    # Performing normal QM/MM job. This also allows for one QM/MM cycle to be performed with a subsequent frequency calculation.
    else:
        cln = 0			
        for ic in range(maxcycle): 		
            cln = cln + 1
            qmmm()
        if maxcycle == 1 and job.lower() == 'freq':
            qomutil.qomlog('Frequency job is started after first QoMMMa run, Frequency job will be performed with the geometry given in ' + inpgeom_prefix + '*.xyz file, here qmcode used is Gaussian', usrdir) 
            from qomfreq import freqmain
            qmc_job = qmgau_job
            freqmain(qommmadir, qmkey, qmc_job, usrdir, cha_mul, extra_basis, gau_head, prj_freq, qmjob_prefix, qmcode) 
        # For the case where the calculation does not converge in the maximum number of cycles, the program is stopped.
        elif cln == maxcycle: 
            qomutil.qomend('QoMMMa run reaches maxcycle :' + str(cln), cwd, usrdir)