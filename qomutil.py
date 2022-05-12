#!/usr/bin/python3

"""
// This is a helper QoMMMa Python file, which contains utility functions used in QoMMMa execution. //
"""

# Global imports.
import os
import sys
import linecache
import shutil
from time import asctime

def qomlog(s, usrdir):
    """
    
    // Function which writes the progress or error message of the given QoMMMa calculation to the log file QoMMMa8.log. //
    
    Arguments
    ----------
    s : string
        The message to be written to the file QoMMMa8.log.
    usrdir : string
        The user directory.

    """
    
    # Simply opens QoMMMa8.log and adds the string, s.
    fd = open(usrdir + '/QoMMMa8.log', 'a') 
    fd.write(s)
    fd.write('\n')
    fd.write('\n')
    fd.close()

def qomend(s, cwd, usrdir):
    """
    
    // Function which ends the given QoMMMa calculation and writes the error to the log file QoMMMa8.log. //
    
    Arguments
    ----------
    s : string
        The message to be written to the file QoMMMa8.log.
    cwd : string
        The current working directory.
    usrdir : string
        The user directory.

    """
 
    # The function asctime from the module time is used to record when the program is ended.
    hlp = asctime()
    qomlog(s + '\n' + '  QoMMMa job ends at ' + str(hlp), usrdir)
    
    # The jobfiles from the failed QoMMMa run are retained for diagnostic purposes.
    # If the jobfiles directory is found from a previous QoMMMa run, then it is copied and renamed to jobfiles_old.
    if os.path.exists(usrdir + '/jobfiles'):
        if os.path.exists(usrdir + '/jobfiles_old'):
            shutil.rmtree(usrdir + '/jobfiles_old')
        shutil.copytree(usrdir + '/jobfiles', usrdir + '/jobfiles_old') 
        shutil.rmtree(usrdir + '/jobfiles')
        qomlog('Note, jobfiles directory is renamed as jobfiles_old if you want to keep this directory rename it with different name before submitting next QoMMMa job', usrdir)
    shutil.copytree(cwd, usrdir + '/jobfiles')
    shutil.rmtree(cwd) 
    
    # Ending QoMMMa job.
    sys.exit()
    
def qomsave(cwd, usrdir):
    # The jobfiles from the QoMMMa run are retained for diagnostic purposes.
    # If the jobfiles directory is found from a previous QoMMMa run, then it is copied and renamed to jobfiles_old.
    if os.path.exists(usrdir + '/jobfiles'):
        if os.path.exists(usrdir + '/jobfiles_old'):
            shutil.rmtree(usrdir + '/jobfiles_old')
        shutil.copytree(usrdir + '/jobfiles', usrdir + '/jobfiles_old') 
        shutil.rmtree(usrdir + '/jobfiles')
        qomlog('Note, jobfiles directory is renamed as jobfiles_old if you want to keep this directory rename it with different name before submitting next QoMMMa job', usrdir)
    shutil.copytree(cwd, usrdir + '/jobfiles')
    if not os.path.exists(usrdir + '/geom1.xyz'):
        shutil.copy(usrdir + '/jobfiles/geom1.xyz', usrdir + '/geom1.xyz')
    shutil.rmtree(cwd) 
    

def addatmnum(imn):
    """
    
    // Function which adds atoms numbers to the QM atoms. //
    // This function is called within the functions rearrange() and setini() found in this file. //
    
    Arguments
    ----------
    imn : integer
        If using the nudged elastic band, this is the image number.

    """
  
    # The QM geometry returned from the QM calculation is opened.
    f = open(('%s%d%s'%('nqmgeom', imn, '.xyz')), 'r')
  
    # The new QM geometry with the correct atom numbers is created.
    fd = open(('%s%d%s'%('qmgeom', imn, '.xyz')), 'w')
    ll = 0
    
    # The atom numbers are added to the file qmgeom*.xyz.
    for line in f:					
        ll = ll+1
        kline = line.split()
        atom = ('%s%d'%(kline[0], ll)).ljust(6)
        res = line[2:]
        fd.write(str(atom))
        fd.write(res)
    f.close()
    fd.close()

def rearrange(nimg, cwd, qmjob_prefix, cln):
    """
    
    // Function which  rearranges the files generated by QoMMMa and QM output to use in subsequent micro-optimisation and QM jobs. //
    // This function is called in the function qmmm found in qommmma.py at each QoMMMa cycle. //
    
    Arguments
    ----------
    nimg : integer
        If using the nudged elastic band, this is the total number of images.
    cwd : string
        The current working directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.
    cln : integer
        The QoMMMa cycle number.        
    """
    
    # For every image (if using the nudged elastic band, that is), appropriate files are rearranged.    
    im = 0
    for ir in range(nimg):
        im = im + 1
        if os.path.exists(cwd + ('%s%d'%('/update_geom', im))):
            dst = cwd + ('%s%d'%('/image', im))
                        
            # Checkfiles and QM region geometries are copied, moved, and deleted.
            shutil.copy(('%s%d'%('CheckFile', im)), dst + (('%s%d'%('/old_CheckFile', im))))
            shutil.copy(('%s%d'%('nCheckFile', im)), ('%s%d'%('CheckFile', im)))
            shutil.copy((dst + ('%s%d%s'%('/qmgeom', im, '.xyz'))), (dst + ('%s%d%s'%('/old_qmgeom', im, '.xyz'))))
            os.remove('%s%d'%('nCheckFile', im))
            os.remove(dst + ('%s%d%s'%('/qmgeom', im, '.xyz'))) 
            addatmnum(im)
            shutil.copy(('%s%d%s'%('qmgeom', im, '.xyz')), dst + (('%s%d%s'%('/qmgeom', im, '.xyz'))))
            os.remove('%s%d%s'%('qmgeom', im, '.xyz'))
            
            # For the next microiteration, the details from tinker_header and inactive_microiter are written to the new file microgeom.key.
            fm = open(('%s%d%s'%('microgeom', im, '.key')), 'w')
            ft = open('tinker_header', 'r')
            fi = open(('%s%d'%('inactive_microiter', im)), 'r')
            for line in ft:
                fm.write(line)
            ft.close()
            for line in fi:
                fm.write(line)
            fi.close()
            fm.close()
            
            # Appropriate files are copied, moved, and deleted.
            shutil.copy(('%s%d%s'%('geom', im, '.xyz')), dst + ('%s%d%s'%('/microgeom', im, '.xyz')))  
            shutil.copy(('%s%d'%('gradcorrection', im)), dst + '/gradcorrection')  
            shutil.copy(('%s%d%s'%('microgeom', im, '.key')), dst + ('%s%d%s'%('/microgeom', im, '.key')))
            os.remove('%s%d%s'%('geom', im, '.xyz'))
            os.remove('%s%d'%('gradcorrection', im))
            os.remove('%s%d%s'%('microgeom', im, '.key')) 
        else:
            shutil.copy(('%s%d'%('nCheckFile', im)), ('%s%d'%('CheckFile', im)))
            os.remove('%s%d'%('nCheckFile', im))

def postrearrange(nimg, cwd):
    """
    
    // Function which rearranges the files after micro-optimisation. //
    // This function is called in the function qmmm found in qommmma.py at each QoMMMa cycle. //
    
    Arguments
    ----------
    nimg : integer
        If using the nudged elastic band, this is the total number of images.
    cwd : string
        The current working directory.
        
    """
    
    # For every image (if using the nudged elastic band, that is), appropriate files are rearranged.    
    im = 0
    for ir in range(nimg):
        im = im + 1
        if os.path.exists(cwd + ('%s%d'%('/update_geom', im))):
            dst = cwd + ('%s%d'%('/image', im))
            
            # Appropriate files are copied, moved, and deleted.
            shutil.copy(('%s%d%s'%('ngeom', im, '.xyz')), ('%s%d%s'%('geom', im, '.xyz')))  
            shutil.copy(('%s%d%s'%('charges', im, '.xyz')),dst + ('%s%d%s'%('/charges', im, '.xyz')))  
            shutil.copy(('%s%d%s'%('charges', im, 'mol.xyz')), dst + ('%s%d%s'%('/charges', im, 'mol.xyz')))  
            shutil.copy(('%s%d%s'%('points', im, '.pts')),dst + ('%s%d%s'%('/points', im, '.pts')))  
            shutil.copy(('%s%d%s'%('points_chg', im, '.pts')), dst + ('%s%d%s'%('/points_chg', im, '.pts')))  
            shutil.copy(('%s%d%s'%('points', im, 'gau.pts')), dst + ('%s%d%s'%('/points', im, 'gau.pts')))  
            os.remove('%s%d%s'%('ngeom', im, '.xyz'))
            os.remove('%s%d%s'%('charges', im, '.xyz'))
            os.remove('%s%d%s'%('charges', im, 'mol.xyz'))
            os.remove('%s%d%s'%('points', im, '.pts'))
            os.remove('%s%d%s'%('points_chg', im, '.pts'))
            os.remove('%s%d%s'%('points', im, 'gau.pts'))

def qomreport(nimg, usrdir, cwd, cln, qomout):
    """
    
    // Function which adds to the new report, checks the convergence and, if the convergence is achieved, create a dummy file with name 'convergence_ok'. //
    // This function is called in the function qmmm found in qommmma.py at each QoMMMa cycle. //
    
    Arguments
    ----------
    nimg : integer
        If using the nudged elastic band, this is the total number of images.
    usrdir : string
        The user directory.
    cwd : string
        The current working directory.
    cln : integer
        The QoMMMa cycle number.  
    qomout : string
        If 'deep' is used, then QM results, geometries, and gradients will be saved at each cycle. 
        If 'None' (as is default) is used, then only results from the previous two (or one) cycles will be saved.
        
    """
  
    # For every image (if using the nudged elastic band, that is), appropriate details are added to the report.  
    im = 0
    result = 'FALSE'
    for i in range(nimg):
        im = im + 1
        
        # The report, and what details are to be added to the report is opened.
        fd = open(usrdir + ('%s%d'%('/report', im)), 'a')
        fi = open(('%s%d'%('add_to_report', im)), 'r')
        
        # If the report contains the message from the Fortran code which indicates that the optimisation has converged, then this is written to QoMMMa8.log.
        for line in fi:
            fd.write(line)
            if nimg == 1:
                if line.strip() == 'Congratulations!!':
                    result = 'TRUE'
                    qomlog('Optimization converged', usrdir) 
        fd.close()
        fi.close()
        os.remove(('%s%d'%('add_to_report', im)))
        
        # The convergence details which are found in the file add_to_update are added to the file update.
        fd = open(cwd + ('%s%d'%('/image', im)) + ('%s%d'%('/update', im)), 'a')
        fi = open(('%s%d'%('add_to_update', im)), 'r')
        for line in fi:
           fd.write(line)
        fd.close()
        fi.close()
        os.remove(('%s%d'%('add_to_update',im)))
        
        # If qomout is 'deep', then QM results, geometries, and gradients will be saved at each cycle. 
        # If qomout is 'None' (as is default), then only results from the previous two (or one) cycles will be saved.
        if os.path.exists(cwd + ('%s%d'%('/update_geom', im))):
           dst = cwd + ('%s%d'%('/image', im))
           if qomout.lower() == 'deep':
             tcln = cln - 1 
             shutil.copy(('%s%d'%('ab_initio', im)), dst + '/old_ab_initio_' + str(tcln))
             shutil.copy(('%s%d'%('mulliken', im)), dst + '/old_mulliken_' + str(tcln))
             shutil.copy(('%s%d%s'%('qmlatgrad', im, '.out')), dst + (('%s%d%s'%('/old_qmlatgrad', im, '.out_'))) + str(tcln))
             shutil.copy(('%s%d%s'%('geom', im, '.xyz')), dst + (('%s%d%s'%('/old_geom', im, '.xyz_'))) + str(tcln)) 
             shutil.copy(('%s%d'%('gradients', im)), dst + '/old_gradients_' + str(tcln))
             shutil.copy(('%s%d'%('CheckFile', im)), dst + (('%s%d%s%d'%('/CheckFile', im, '_', tcln))))
           else:
             shutil.copy(('%s%d'%('ab_initio', im)), dst + '/old_ab_initio')
             shutil.copy(('%s%d'%('mulliken', im)), dst + '/old_mulliken')
             shutil.copy(('%s%d%s'%('qmlatgrad', im, '.out')), dst + (('%s%d%s'%('/old_qmlatgrad', im, '.out'))))
             shutil.copy(('%s%d%s'%('geom', im, '.xyz')), dst + (('%s%d%s'%('/old_geom', im, '.xyz')))) 
           os.remove('%s%d%s'%('qmlatgrad', im, '.out'))
           os.remove('%s%d'%('ab_initio', im))
           os.remove('%s%d'%('mulliken', im))
    
    # If the report contains the message from the Fortran code which indicates that the optimisation has converged, then this is written to QoMMMa8.log. 
    # This is fort the case of the nudged elastic band.
    if nimg != 1:
        fm = open(usrdir + '/master_report', 'a')
        fi = open('add_to_master_report', 'r')
        for line in fi:
            fm.write(line)
            if line.strip() == 'Congratulations!!':		
                result='TRUE'
                qomlog('Optimization converged', usrdir)
            if line[:30] == 'Geomtries of all images fixed ':
                result = 'TRUE'
                qomlog('All images are fiexd', usrdir)   
        fi.close()
        fm.close()
        os.remove('add_to_master_report')
    
    # If the string result is 'TRUE', then optimisation has converged and the dummy file convergence_ok is created.
    if result=='TRUE':
        fc=open(cwd+'/convergence_ok','w')
        fc.close() 
 
def setini(cwd, nimg, usrdir):
    """
    
    // Function which creates the file update_geom, and moves the initial files created by the Fortran QoMMMa setup program to the image directories. //
    // This function is called in qommmma.py only once. //
    
    Arguments
    ----------
    cwd : string
        The current working directory.
    nimg : integer
        If using the nudged elastic band, this is the total number of images.        
    usrdir : string
        The user directory.

    """
      
    # For every image (if using the nudged elastic band, that is), the initial setup is performed.  
    l = 0  		
    for i in range(nimg):     
        l = l + 1
        fu = open(('%s%d'%('update_geom', l)), 'w')
        dst = cwd + ('%s%d'%('/image', l))
        
        # Atom number are added to the QM atoms.
        addatmnum(l)
        shutil.copy(('%s%d%s'%('qmgeom', l, '.xyz')), dst + ('%s%d%s' % ('/qmgeom', l, '.xyz')))
        
        # If an report file from a previous QoMMMa run is found, then it is renamed and retained.
        if os.path.exists(usrdir + ('%s%d'%('/report', l))):
            shutil.copy((usrdir + ('%s%d'%('/report', l))), (usrdir + ('%s%d%s'%('/report', l, '_'))))
            qomlog('Note, Previous report file in user directory was moved to report*_, before submitting next job if you want rename it or else it will be deleted', usrdir)
            
        # Files generates by the Fortran based setup are copied, moved, and deleted.
        shutil.copy(('%s%d'%('report_header', l)), usrdir + ('%s%d' % ('/report', l)))
        shutil.copy(('%s%d%s'%('charges', l, '.xyz')), dst)
        shutil.copy(('%s%d%s'%('charges', l, 'mol.xyz')), dst)
        shutil.copy(('%s%d%s'%('points', l, '.pts')), dst)
        shutil.copy(('%s%d%s'%('points_chg', l, '.pts')), dst)
        shutil.copy(('%s%d%s'%('points', l, 'gau.pts')), dst)
        os.remove('%s%d%s'%('qmgeom', l, '.xyz'))
        os.remove('%s%d'%('report_header', l))
        os.remove('%s%d%s'%('charges', l, '.xyz'))
        os.remove('%s%d%s'%('charges', l, 'mol.xyz'))
        os.remove('%s%d%s'%('points', l, '.pts'))    
        os.remove('%s%d%s'%('points_chg', l, '.pts'))    
        os.remove('%s%d%s'%('points', l, 'gau.pts'))    
    fu.close()

def tinkercharge(qommmadir):
    """
    
    // Function which runs Tinker's analyze_grad, and extracts the intial charges from the output file. //
    // This function is called in the function mmjob found in qommmma.py. //
    
    Arguments
    ----------
    qommmadir : string
        The base QoMMMa directory defined in the file qommma.

    """
    
    # The Tinker program analyze_grad is called.
    fch = os.popen(qommmadir + '/bin/analyze_grad > tmpch.out', 'w')
    fch.write('geom1.xyz\nP\nALL\n')
    fch.close()
    
    # The file produced by analyze_grad - tmpch.out - is opened.
    fil = open('tmpch.out', 'r')
    
    # A new file containing the point charges is created.
    fdch = open('DefaultCharges', 'w')
    
    # While loop iterates through tmpch.out and extracts the relevant information on atomic partial charges.
    while 1:
        line = fil.readline()
        if not line:break
        if line[:33] == ' Atomic Partial Charge Parameters':
            fil.readline()
            fil.readline()
            fil.readline()
            im = 0
            while 1:
                line = fil.readline()
                im = im + 1
                if not line:break
                lnu = int(line.split()[1])
                cha = float(line.split()[2])
                if lnu == im:
                    fdch.write(str(im).rjust(6))
                    fdch.write(str(cha).rjust(11))
                    fdch.write('\n')
                else:
                    kk = 0
                    kk = lnu - im
                    for i in range(kk + 1):
                        fdch.write(str(im).rjust(6))
                        if im == lnu: 
                            fdch.write(str(cha).rjust(11))
                            fdch.write('\n')
                        else:
                            fdch.write(('0.0000').rjust(11))
                            im = im + 1
                            fdch.write('\n') 		          
    fil.close()
    fdch.close()

def geom_expl(f, l):
    """
    
    // Function which creates a coordinate file (geom_expl*.xyz) to use as a coordinate file in QoMMMa optimisation and setup program. //
    // This function is called in the functions mmjob and qmmm found in qommmma.py. //
    
    Arguments
    ----------
    f : string
        The file name for the input geometry file.
    l : integer
        If using the nudged elastic band, this is the image number.

    """
    
    # A coordinate file is opened and formatted to be used in QoMMMa optimsation.
    fi = open(f, 'r')         		
    fd = open(('%s%d%s'%('geom_expl', l, '.xyz')), 'w')         		  
    line = fi.readline()
    line = line.split()
    fd.write(line[0].rjust(6))
    fd.write('\n') 
    n = int(line[0])
    fi.seek(0, 1)
    for il in range(n):
        line = fi.readline()
        k = len(line.split()[6:])
        fd.write(str(k).rjust(8))
        fd.write('\n')
        fd.write(line)
    fi.close()
    fd.close()

def qmread(kk, usrdir):
    """
    
    // Function which generates the QM atom list. //
    // This function is called in the function fortinp found in qommmma.py. //
    
    Arguments
    ----------
    kk : list/range of integers
        Represents the QM atoms which have been selected.
    usrdir : string
        The user directory.

    """
  
    # Files are opened, their names initialised, and the atom labels are initialised.
    fds = open('fortinput', 'a')
    fi = 'geom1.xyz'
    PTab = ['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os','Ir','Pt','Au','Hg','Tl','bp','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg']
    qomlog('Following QM atoms list is generated from qmread module in qomutil based on user input given in qm_lst key. If this list is not the one as you except check qm_lst input or break QoMMMa calculation and give your full QM atom list using qm key. Note, to use this key, if your atom name consist of two characters then use small case for second character in tinker xyz (initial geometry) file. For example; Fe, Cl, He (for iron, chlorine, helium) is correct but FE, CL,HE is not correct.', usrdir)
    nqmlst = []

    # The atom number are appended to the list nqmlst.
    # The if block handles a range as input [1-5], and the else block handles a simple list as input [1,2,3,4,5].
    for i in range(len(kk)):
        if '-' in str(kk[i]):
            a, b = kk[i].split('-')
            int_a, int_b = int(a), int(b)
            nqmlst.extend(range(int_a, (int_b + 1)))
        else:
            nqmlst.append(int(kk[i]))

    # If the atoms numbers are not in numerical order, then they are reordered such that they are.
    for i in range(len(nqmlst)):
        j = i
        for k in range(len(nqmlst) - 1 - i):
            j += 1
            if nqmlst[i] > nqmlst[j]:
                tt = nqmlst[i]
                nqmlst[i] = nqmlst[j]
                nqmlst[j] = tt

    # Each QM atom number along with its atom type is written to fortinput.
    for i in range(len(nqmlst)):
        iline = nqmlst[i]
        line = (linecache.getline(fi, int(iline) + 1)).split()
        aaname = line[1]
        try:
            aname = str.upper(aaname[0]) + aaname[1]
        except:
            aname = str.upper(aaname[0])
        if aname in PTab:
            fds.write(str(iline).rjust(6))
            fds.write(aname.rjust(10))
            qomlog(str(iline).rjust(6) + str(aname).rjust(10), usrdir)
        elif aname[0] in PTab:
            fds.write(str(iline).rjust(6))
            fds.write(str.upper(aname[0]).rjust(10))
            qomlog(str(iline).rjust(6) + str(aname[0]).rjust(10), usrdir)
        else:
            fds.write(str(iline).rjust(6))
            fds.write(aname.rjust(10))
            qomlog('Note, Unknown atom name : ' + aname + ' at line number : ' + str(iline) + ' is taken', usrdir)    
        fds.write('\n')
    qomlog('Total number of QM atoms in the above list : ' + str(len(nqmlst)), usrdir)
    fds.close()

def hesopt_write(hesopt, usrdir):
    """
    
    // Function which generates the list of hessian optimised MM atoms. //
    // This function is called in the function fortinp found in qommmma.py. //
    
    Arguments
    ----------
    hesopt : list/range of integers
        Represents the MM atoms which have been selected for hessian optimisation.
    usrdir : string
        The user directory.

    """
    
    # The file fortinput is opened.
    fd = open('fortinput', 'a')
    
    # The hessian optimised MM atoms are written to fortinput.
    # The if block handles a range as input [1-5], and the else block handles a simple list as input [1,2,3,4,5].
    itot = 0
    for i in range(len(hesopt)):
        if '-' in str(hesopt[i]):
            ii = hesopt[i]
            it = 0
            for element in ii:
                it = it + 1
                if element == '-':
                    break
            inl = int(ii[it:]) - int(ii[0:it-1])
            for k in range(inl + 1):
                ival = int(ii[0:it - 1]) + k
                fd.write(str(ival).rjust(8))
                itot += 1
                fd.write('\n')
        else:
            itot += 1  
            fd.write(str(hesopt[i]).rjust(8))
            fd.write('\n')
    qomlog('Total number of Hessian optimized MM atoms is : ' + str(itot), usrdir)
    fd.close()

def link_write(link):
    """
    
    // Function which generates the details for link atoms. //
    // This function is called in the function fortinp found in qommmma.py. //
    
    Arguments
    ----------
    link : tuple
        Represents the details for each link atom.

    """
    
    # The file fortinput is opened.
    fd = open('fortinput','a')
    
    # If the link atom numbers are not in numerical order, then they are reordered such that they are.
    # This is important for a frequency calculation.
    if len(link) > 1:
        for i in range(len(link)):
            j = i
            for k in range(len(link) - 1 - i):
                j += 1
                if link[i][1] > link[j][1]:
                    tt = link[i]
                    link[i] = link[j]
                    link[j] = tt
                    
    # Link atom details are written to fortinput.
    for i in range(len(link)):
        for j in range(len(link[i])):
            fd.write(str(link[i][j]).rjust(8))
        fd.write('\n')  
    fd.close()
    
def cons_write_prim(cons, usrdir):
    """
    
    // Function which generates the constraint details for primitive internal coordinates. //
    // This function is called in the function fortinp found in qommmma.py. //
    
    Arguments
    ----------
    cons : tuple
        Represents the details for each constraint.
    usrdir : string
        The user directory.
        
    """
    
    # The file fortinput is opened.
    fd = open('fortinput', 'a')

    for i in range(len(cons)):
        # The constraint details are written to fortinput.
        fd.write(str(cons[i][0]).rjust(8)) # The desired change in primitive coordinates.
        fd.write('\n')
        for j in range(len(cons[i][1])):
            fd.write(str(cons[i][1][j]).rjust(8))
        fd.write('\n')
    fd.close()
    
def cons_write_cart(cons, usrdir):
    """
    
    // Function which generates the constraint details for cartesian coordinates. //
    // This function is called in the function fortinp found in qommmma.py. //
    
    Arguments
    ----------
    cons : tuple
        Represents the details for each constraint.
    usrdir : string
        The user directory.
        
    """
    
    # The file fortinput is opened.
    fd = open('fortinput', 'a')
    
    # Depending on the constraint type used, the integer ctyp is set.
    for i in range(len(cons)):
        if cons[i][0].upper() == 'RAB':
            ctyp = 1
        elif cons[i][0].upper() == 'RAB-RCD':
            ctyp = 2
        elif cons[i][0].upper() == 'RAB+RCD':
            ctyp = 3
        elif cons[i][0].upper() == 'RAB+RCD-REF':
            ctyp = 4
        else:
            qomlog('ERROR, Unknown constraint type :' + cons[i][0], usrdir)
      
        # The constraint type along with its details are written to fortinput.
        fd.write(str(ctyp).rjust(3))  
        fd.write(str(cons[i][1]).rjust(8))
        fd.write(str(cons[i][2]).rjust(8))
        fd.write('\n')
        for j in range(len(cons[i][3])):
            fd.write(str(cons[i][3][j]).rjust(8))
        fd.write('\n')
    fd.close()

def newcha_write(newcha):
    """
    
    // Function which generates a list of MM atoms with their newly assigned charge. //
    // This function is called in the function fortinp found in qommmma.py. //
    
    Arguments
    ----------
    newcha : tuple
        Represents the details for each new charge for MM atoms.

    """
    
    # The file fortinput is opened.
    fd = open('fortinput', 'a')
    
    # The modified MM charges are written to fortinput.
    fd.write(str(len(newcha)).rjust(6))
    fd.write('\n')
    for i in range(len(newcha)):
        fd.write(str(newcha[i][0]).rjust(8))
        fd.write(str(newcha[i][1]).rjust(8))
        fd.write('\n')
    fd.close()

def inprange(kk, job, usrdir):
    """
    
    // Function which writes atom labels to the file fortinput. //
    // This function is called in the function fortinp found in qommmma.py. //
    
    Arguments
    ----------
    kk : list/range of integers
        Represents the atoms which have been selected.
    job : string
        The purpose of this input.
    usrdir : string
        The user directory.    

    """
    
    # The file fortinput is opened.
    fd = open('fortinput', 'a')
    
    # The atom labels for the given jobtype are written to fortinput.
    # The if block handles a range as input [1-5], and the else block handles a simple list as input [1,2,3,4,5].
    itot = 0
    for i in range(len(kk)):
     if '-' in kk[i]:
       ii = kk[i]
       it = 0
       for element in ii:
          it = it + 1
          if element == '-':
            break
       inl = (int(ii[it:]) - int(ii[0:it - 1])) + 1
       for k in range(inl):
           iline = int(ii[0:it - 1]) + k
           itot = itot + 1
           fd.write(str(iline))
           fd.write('\n')
     else:
       itot = itot + 1
       iline = kk[i]
       fd.write(str(iline))
       fd.write('\n')
    qomlog('Number of  ' + job + ' atoms is ' + str(itot), usrdir)
    fd.close()