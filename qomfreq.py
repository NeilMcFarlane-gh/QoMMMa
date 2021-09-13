#!/usr/bin/python3

"""
// This is a helper QoMMMa Python file, which contains the functions used in execution of a QoMMMa frequency calculation. //
"""

# Global imports.
import os
import shutil
    
# Local imports.
from qomutil import qomend, qomlog   
    
def initial(nimg, cwd, qommmadir):
    """
    
    // Function which sets up the intial arrangements for a frequency calculation. //
    
    Arguments
    ----------
    nimg : integer
        If using the nudged elastic band, this is the total number of images.
    cwd : string
        The current working directory.
    qommmadir : string
        The base QoMMMa directory defined in the file qommma.

    """
    
    im = 0
    for i in range(nimg):
        im = im + 1
        
        # The geometry is opened.
        fi = '%s%d%s'%('geom', im, '.xyz')
        
        # The function geom_exp1 found in qomutil.py is imported to get geom_exp1*.xyz.
        from qomutil import geom_expl
        geom_expl(fi, im)
        
        # A directory /freq* is created.
        dst = cwd + ('%s%d'%('/image', im)) + ('%s%d'%('/freq', im))
        os.mkdir(dst)
        
        # Relevant files are copied to the newly created /freq*
        shutil.copy(cwd + ('%s%d%s'%('/geom', im, '.xyz')), dst + ('%s%d%s'%('/geom', im, '.xyz')))  
        shutil.copy('geom1.key', dst + ('/mgrad.key'))
    
    # The Fortran based QoMMMa setup program is called.
    os.system(qommmadir + '/bin/setup8.x') 
    
    im = 0
    for i in range(nimg):
        im = im + 1
        
        # The function addatmnum found in qomutil.py is imported to add atom numbers.
        from qomutil import addatmnum
        addatmnum(im)
        
        # Files created by the QoMMMa setup program are copied to /freq* and subsequently removed.
        dst = cwd + ('%s%d'%('/image', im)) + ('%s%d'%('/freq', im))
        shutil.copy(('%s%d%s'%('qmgeom', im, '.xyz')), dst + ('%s%d%s'%('/qmgeom', im, '.xyz')))
        shutil.copy(('%s%d%s'%('nqmgeom', im, '.xyz')), dst + ('%s%d%s'%('/nqmgeom', im, '.xyz')))
        shutil.copy(('%s%d%s'%('charges', im, '.xyz')), dst)
        os.remove('%s%d'%('report_header', im))
        os.remove('%s%d%s'%('charges', im, '.xyz'))
        os.remove('%s%d%s'%('qmgeom', im, '.xyz'))
        os.remove('%s%d%s'%('nqmgeom', im, '.xyz'))

def mminpco(im, qmlst, delq):
    """
    
    // Function which generates the 2*3*nqm number of Tinker coordinate files to calculate the MM gradients used for a frequency calculation. //
    // This function is called in freqmain for each image, if using the nudged elastic band. //
    
    Arguments
    ----------
    imn : integer
        If using the nudged elastic band, this is the image number.
    qmlst : list
        The list of QM atoms.
    delq : integer
        Real constant representing the finite difference used when generating coordinates.
    
    """
    
    # The number of QM atoms is initialised.
    nqm = len(qmlst)
    
    # The geometry file is opened and the total number of atoms is initialised.
    fi = open(('%s%d%s'%('geom', im, '.xyz')), 'r')
    n = int(fi.readline().split()[0])
    
    # The tuples alabel, aconect, and acoord are initialised and completed with the relevant data found in geom*.xyz.
    ii = 0
    alabel = {}
    aconect = {}
    acoord = {}
    for i in range(n):
        ii = ii + 1
        line=fi.readline().split()
        alabel[ii] = line[1]
        acoord[((ii-1)*3)+1] = line[2]
        acoord[((ii-1)*3)+2] = line[3]
        acoord[((ii-1)*3)+3] = line[4]
        aconect[ii] = line[5:]
    fi.close()
    
    # Some counters are initialised.
    k = 0
    sn = 0
    irun = 0
    
    # While loop iterates until the counter irun has reached the total number of QM atoms.
    while 1:
        k = k + 1
        if k in qmlst:
            
            # The tuples ancoord and apcoord are initialised.
            ancoord = {}
            apcoord = {}
            
            # The counter irun is used as the end condition.
            irun = irun + 1
            
            # For-loop iterates over the total number of atoms, n.
            ii = 0
            for i in range(n):
                ii = ii + 1 
                if k == ii:
                    ik = 0
                    for il in range(3):
                        ik = ik + 1
                        sn = sn + 1 
                        
                        # The coordinates from acoord are coped to the tuples ancoord and apcoord.
                        ancoord[((ii - 1) * 3) + 1] = float(acoord[((ii - 1) * 3) + 1])
                        apcoord[((ii - 1) * 3) + 1] = float(acoord[((ii - 1) * 3) + 1])
                        ancoord[((ii - 1) * 3) + 2] = float(acoord[((ii - 1) * 3) + 2])
                        apcoord[((ii - 1) * 3) + 2] = float(acoord[((ii - 1) * 3) + 2])
                        ancoord[((ii - 1) * 3) + 3] = float(acoord[((ii - 1) * 3) + 3])
                        apcoord[((ii - 1) * 3) + 3] = float(acoord[((ii - 1) * 3) + 3])
                        
                        # Using the finite difference, new coordinates are generated.
                        if ik == 1:
                            ancoord[((ii - 1) * 3) + 1] = ancoord[((ii - 1) * 3) + 1] - delq
                            apcoord[((ii - 1) * 3) + 1] = apcoord[((ii - 1) * 3) + 1] + delq
                        if ik == 2:
                            ancoord[((ii - 1) * 3) + 2] = ancoord[((ii - 1) * 3) + 2] - delq
                            apcoord[((ii - 1) * 3) + 2] = apcoord[((ii - 1) * 3) + 2] + delq
                        if ik == 3:
                            ancoord[((ii - 1) * 3) + 3] = ancoord[((ii - 1) * 3) + 3] - delq
                            apcoord[((ii - 1) * 3) + 3] = apcoord[((ii - 1) * 3) + 3] + delq
                        irem = n - ii
                        iin = ii
                        for ir in range(irem):
                            iin = iin + 1
                            ancoord[((iin - 1) * 3) + 1] = acoord[((iin - 1) * 3) + 1]
                            apcoord[((iin - 1) * 3) + 1] = acoord[((iin - 1) * 3) + 1]
                            ancoord[((iin - 1) * 3) + 2] = acoord[((iin - 1) * 3) + 2]
                            apcoord[((iin - 1) * 3) + 2] = acoord[((iin - 1) * 3) + 2]
                            ancoord[((iin - 1) * 3) + 3] = acoord[((iin - 1) * 3) + 3]
                            apcoord[((iin - 1) * 3) + 3] = acoord[((iin - 1) * 3) + 3]
                        f = '%s%d%s%d%s'%('mmgradco', im, '_qmn', sn, '.xyz')
                        fn = open(f, 'w')
                        f1 = '%s%d%s%d%s'%('mmgradco', im, '_qmp', sn, '.xyz')
                        fp = open(f1, 'w')
                        iir = 0
                        fn.write(str(n).rjust(6))
                        fn.write('\n')
                        fp.write(str(n).rjust(6))
                        fp.write('\n')
                        for ip in range(n):
                            iir = iir + 1
                            fn.write(str(iir).rjust(6))
                            fn.write(str(alabel[iir]).rjust(4))
                            fn.write(str(ancoord[((iir - 1) * 3) + 1]).rjust(12))
                            fn.write(str(ancoord[((iir - 1) * 3) + 2]).rjust(12))
                            fn.write(str(ancoord[((iir - 1) * 3) + 3]).rjust(12))
                            for item in aconect[iir]: 
                                fn.write(item.rjust(6))
                            fn.write('\n')
                            fp.write(str(iir).rjust(6))
                            fp.write(str(alabel[iir]).rjust(4))
                            fp.write(str(apcoord[((iir - 1) * 3) + 1]).rjust(12))
                            fp.write(str(apcoord[((iir - 1) * 3) + 2]).rjust(12))
                            fp.write(str(apcoord[((iir - 1) * 3) + 3]).rjust(12))
                            for item in aconect[iir]: 
                                fp.write(item.rjust(6))
                            fp.write('\n')
                        fn.close()
                        fp.close()
                else: 
                    ancoord[((ii - 1) * 3) + 1] = acoord[((ii - 1) * 3) + 1]
                    apcoord[((ii - 1) * 3) + 1] = acoord[((ii - 1) * 3) + 1]
                    ancoord[((ii - 1) * 3) + 2] = acoord[((ii - 1) * 3) + 2]
                    apcoord[((ii - 1) * 3) + 2] = acoord[((ii - 1) * 3) + 2]
                    ancoord[((ii - 1) * 3) + 3] = acoord[((ii - 1) * 3) + 3]
                    apcoord[((ii - 1) * 3) + 3] = acoord[((ii - 1) * 3) + 3]
        # End condition.
        if irun == nqm: break

def analgrad(im, qmlst, qommmadir):
    """
    
    // Function which runs Tinker's analyze_grad for the structures generated by mminpco. //
    
    Arguments
    ----------
    im : integer
        If using the nudged elastic band, this is the image number.
    qmlst : list
        The list of QM atoms.
    qommmadir : string
        The base QoMMMa directory defined in the file qommma.
    
    """
    
    # The number of QM atoms is initialised.
    nqm = len(qmlst)
    
    # For-loop which iterates over all structures.
    ii = 0
    for i in range(nqm * 3):
        ii = ii + 1
        
        # Performing an MM gradient evaluation for mmgradco*_qmp.xyz.
        fip = '%s%d%s%d%s'%('mmgradco', im, '_qmp', ii, '.xyz')
        shutil.copy(fip, 'mgrad.xyz')
        fch = os.popen(qommmadir + '/bin/analyze_grad', 'w')
        fn = 'mgrad.xyz'
        fch.write(fn + '\nE\n')
        fch.close()
        
        # The evaluated gradient is copied and saved as mmgrad*_qmp.
        shutil.copy('mm_grad', ('%s%d%s%d'%('mmgrad', im, '_qmp', ii)))
        os.remove('mgrad.xyz')
        os.remove('mm_grad') 
        os.remove(fip)
        
        # Performing an MM gradient evaluation for mmgradco*_qmn.xyz.        
        fin = '%s%d%s%d%s'%('mmgradco', im, '_qmn', ii, '.xyz')
        shutil.copy(fin, 'mgrad.xyz')
        fch = os.popen(qommmadir + '/bin/analyze_grad', 'w')
        fn = 'mgrad.xyz'
        fch.write(fn + '\nE\n')
        fch.close()
        
        # The evaluated gradient is copied and saved as mmgrad*_qmm.
        shutil.copy('mm_grad', ('%s%d%s%d'%('mmgrad', im, '_qmn', ii)))
        os.remove('mgrad.xyz')
        os.remove('mm_grad') 
        os.remove(fin)

def mmhes(im, qmlst, n, delq):
    """
    
    // Function which calculates the MM hessian from the gradients calculated by Tinker's analyze_grad. //
    // This function is called in freqmain for each image, if using the nudged elastic band. //
    
    Arguments
    ----------
    im : integer
        If using the nudged elastic band, this is the image number.
    qmlst : list
        The list of QM atoms.
    n : integer
        The total number of atoms.
    delq : integer
        Real constant representing the finite difference used when generating coordinates.
    
    """
    
    # Conversions between kcal/mol and a.u., and bohr and angstroms are initialised.
    kcal_au = 627.5095
    bo_ang = 0.529177249
    
    # The number of QM atoms is initialised.
    nqm = len(qmlst)
    
    # The dimensions of the MM hessian is initialised.
    mhes=[[{} for i in range (nqm*3)] for i in range(nqm*3)]
    
    # For-loop iterates over all gradients calculated by analyze_grad.
    ii = 0 
    for i in range(nqm * 3):
        ii = ii + 1
        
        # The files generated for the suffix _qmp and _qmn are both opened.
        fp = '%s%d%s%d'%('mmgrad', im, '_qmp', ii)
        fn = '%s%d%s%d'%('mmgrad', im, '_qmn', ii) 
        fip = open(fp, 'r')
        fin = open(fn, 'r')   
        fip.readline()
        fin.readline()
        
        # Some counters are initialised.
        iql=0
        mp=0
        
        # The relevant terms are added to their positions in the MM hessian using a for-loop.
        # Units are converted at this stage.
        for il in range(n):
            iql = iql + 1
            linep = fip.readline().split()
            linen = fin.readline().split()
            if iql in qmlst:
                mp = mp + 1
                mhes[ii - 1][mp - 1] = ((float(linep[2]) - float(linen[2])) / kcal_au / (2 * delq) * bo_ang * bo_ang)
                mp = mp + 1
                mhes[ii - 1][mp - 1] = ((float(linep[3]) - float(linen[3])) / kcal_au / (2 * delq) * bo_ang * bo_ang)
                mp = mp + 1
                mhes[ii - 1][mp - 1] = ((float(linep[4]) - float(linen[4])) / kcal_au / (2 * delq) * bo_ang * bo_ang)
        fip.close()
        fin.close()
        os.remove(fp)
        os.remove(fn)
    
    # The symmetry of the MM hessian is ensured by calculating the average of the symmetric positions.
    for i in range(nqm * 3):
        for j in range(nqm * 3):
            mhes[i][j] = ((mhes[i][j]) + (mhes[j][i])) / 2
            mhes[j][i] = mhes[i][j]
     
    # The file containing the MM hessian is created and the data found in mhes is formatted appropriately.
    fh = open(('%s%d'%('mmhes', im)), 'w')
    indx = 0
    for i in range(nqm * 3):
        indx = indx + 1
        hlp = 0
        fh.write(str(indx).rjust(4))
        inum = 0         
        for j in range(i + 1):
            inum = inum + 1
            hlp = hlp + 1
            fh.write(str(mhes[i][j]).rjust(20))
            if (hlp == 4) and (inum != (i + 1)):
                fh.write('\n')
                fh.write(str(indx).rjust(4))         
                hlp=0
        fh.write('\n')
    fh.close()

def freqinp(im, qmlst, cwd):
    """
    
    // Function which creates the input file (hessian_input) required for a frequency calculation. //
    
    Arguments
    ----------
    im : integer
        If using the nudged elastic band, this is the image number.
    qmlst : list
        The list of QM atoms.
    cwd : string
        The current working directory.
    
    """
    
    # The number of QM atoms is initialised.
    nqm = len(qmlst)
    
    # The file hessian_input is created.
    fh = open('hessian_input', 'w')
    
    # The total number of QM atoms is written to the input file.
    fh.write('Number of QM atoms')
    fh.write('\n')
    fh.write(str(nqm))
    fh.write('\n')
    
    # The cartesian coordinates of the QM region are written to the input file.
    fh.write('Coordinates (in angstrom) of QM atoms')
    fh.write('\n')
    f = open('%s%d%s'%('nqmgeom', im, '.xyz'), 'r')
    for i in range(nqm):
        line = f.readline()
        fh.write(line)
    f.close()
    
    # The QM and MM hessian are combined by simply adding them together and written to the input file as the total hessian.
    fh.write('Total Hessian')
    fh.write('\n')
    fq = open('%s%d'%('qmhes', im,), 'r')
    fm = open('%s%d'%('mmhes', im,), 'r')
    while 1:
        lineq = fq.readline().split()
        if not lineq:break
        linem = fm.readline().split()
        if not linem:break
        fh.write(lineq[0].rjust(3))
        for i in range((len(lineq)) - 1): 
            res = float(lineq[i + 1]) + float(linem[i + 1])
            fh.write(str(res).rjust(20))
        fh.write('\n')
    fm.close()
    fq.close()
    
    # If the gradient file is found, then it is written to the input file.
    try:
        f = open(cwd + ('%s%d'%('/gradients', im)), 'r')
        fh.write('Total Gradient')
        fh.write('\n')
        f.readline()
        ln = 0
        chk = 0 
        while 1:
            line = f.readline()
            if not line:break   
            ln = ln + 1
            if ln in qmlst:
                chk = chk + 1
                fh.write(line)
            if chk == nqm: break
        f.close()
    except:pass
    fh.close()

def hescorrect(qhes, ntotqm, nqm, cwd, im):
    """
    
    // Function which corrects the QM hessian if any link atoms are inlcuded in the QoMMMa frequency calculation. //
    
    Arguments
    ----------
    qhes : array
        The QM hessian array.
    ntotqm : integer
        Number of QM atoms including link atoms.
    nqm : integer
        Number of QM atoms.
    cwd : string
        The current working directory.
    im : integer
        If using the nudged elastic band, this is the image number.

    """
    
    # The number of link atoms is initialised by subtracting the user input QM atoms from the total QM atoms.
    nlink = ntotqm - nqm 
    
    # The link atom details are opened.
    ff = open(cwd + '/link_details')
    
    # The list qm_lk contains details of which QM atom is connected to which link atom.
    qm_lk = []
    
    # The list alpha contains the alpha value corresponding to each link atom.
    alpha = []
    
    # By reading the file link_details, the QM atom and the alpha value are read and added to appropriate lists.
    alpha.append(0) 		
    for i in range(nlink):
        line = ff.readline().split()
        qm_lk.append(int(line[0]))
        alpha.append(float(line[1]))     
    ff.close() 
    
    # The list qmllab contains the label of Hessian corresponding to QM atom connected with link atom.	
    qmllab = []	
    qlh = [[] for i in range(1) for j in range(nlink + 1)]
    k = 0	       
    for item in qm_lk:
        hlp = (item * 3) - 3	
        k = k + 1
        for i in range(3):
            qmllab.append(hlp + i + 1)	
            qlh[k].append(hlp + i + 1)  # This is useful to find alpha value 
    
    # The hessian matrix is corrected with consideration to link atoms included in the calculation.
    i = 0
    for ir in range(nqm * 3):
        i = i + 1
        j = 0
        for jr in range(ir + 1):
            j = j + 1 
            if i in qmllab:
                ip = qmllab.index(i) + 1
                if j in qmllab:
                    jp = qmllab.index(j) + 1
                    if i == j:
                        rr = 0
                        for r in range(nlink):
                            rr = rr + 1
                            if i in qlh[rr]:
                                const = alpha[rr]
                        hlp = (1 - const) * (1 - const) * (float(qhes[nqm * 3 + ip][nqm * 3 + jp])) 
                        qhes[i][j] = float(qhes[i][j]) + (2 * (1 - const) * float(qhes[nqm * 3 + ip][j])) + hlp
                    else:
                        rr = 0
                        for r in range(nlink):
                            rr = rr + 1
                            if i in qlh[rr]:
                                consti = alpha[rr]
                            if j in qlh[rr]:
                                constj = alpha[rr]
                        hlp = (1 - consti) * (1 - constj) * (float(qhes[nqm * 3 + ip][nqm * 3 + jp]))  
                        qhes[i][j] = float(qhes[i][j]) + (1 - consti) * (float(qhes[nqm * 3 + ip][j])) + (1 - constj) * (float(qhes[nqm * 3 + jp][i])) + hlp
                else:
                    rr = 0
                    for r in range(nlink):
                        rr = rr + 1
                        if i in qlh[rr]:
                            const = alpha[rr]
                    qhes[i][j] = float(qhes[i][j]) + ((1 - const) * float(qhes[nqm * 3 + ip][j]))
            elif j in qmllab:
                rr = 0
                for r in range(nlink):
                    rr = rr + 1
                    if j in qlh[rr]:
                        const = alpha[rr]
                jp = qmllab.index(j) + 1
                qhes[i][j] = float(qhes[i][j]) + ((1 - const) * float(qhes[nqm * 3 + jp][i]))
                
    # The hessian is appropriately formatted and written to an output file. 
    fh = open(('%s%d'%('qmhes', im)), 'w')
    indx = 0
    for i in range(nqm * 3):
        indx = indx + 1
        hlp = 0
        fh.write(str(indx).rjust(4))
        inum = 0         
        for j in range(i + 1):
            inum = inum + 1
            hlp = hlp + 1
            fh.write(str(qhes[i + 1][j + 1]).rjust(20))
            if (hlp == 4) and (inum != (i + 1)):
                fh.write('\n')
                fh.write(str(indx).rjust(4))         
                hlp = 0
        fh.write('\n')
    fh.close()

def qmjaginp(im, qmkey, cwd, qmjob_prefix):
    """
    
    // Function which creates the input file required for a Jaguar frequency calculation. //
    // This function is called in the function freqmain for each image. //
    
    Arguments
    ----------
    im : integer
        If using the nudged elastic band, this is the image number.
    qmkey : string
        Can either be None or the user input level of theory.
    cwd : string
        The current working directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.

    """
    
    # The Jaguar input file is created for the image.    
    fj = open(('%s%d%s'%(qmjob_prefix, im, '.in')), 'w')
    
    # The details found in qmkey are written to the input file.
    fj.write('&gen')
    fj.write('\n')
    if qmkey.strip() != 'None':
        fj.write(qmkey.lstrip())
    fj.write('ifreq=1')
    fj.write('\n')
    fj.write('&')
    fj.write('\n')
    fj.write('&zmat')
    fj.write('\n')
    
    # The cartesian coordinates of the QM region are written to the input file.       
    f = open(('%s%d%s'%('qmgeom', im, '.xyz')), 'r') 
    for line in f:
        fj.write(line)
    f.close()
    fj.write('&')
    fj.write('\n')
    
    # The details in the file atomic_section (generated by the function atomic() in jagutil.py) are written to the input file.    
    f = open((cwd + '/atomic_section'), 'r')
    for line in f:
        fj.write(line)
    f.close()
    fj.write('\n')
    
    # The point charge array is written to the input file.    
    f = open(('%s%d%s'%('charges', im, '.xyz')), 'r') 
    for line in f:
        fj.write(line)
    f.close()
    fj.close() 

def jaghesread(im, ntotqm, nlink, cwd, qmjob_prefix):
    """
    
    // Function which extracts the QM hessian from a Jaguar output file. //
    // This function is called in the function freqmain for each image. //
    
    Arguments
    ----------
    im : integer
        If using the nudged elastic band, this is the image number.
    ntotqm : integer
        Number of QM atoms including link atoms.
    nlink : integer
        Number of link atoms.
    cwd : string
        The current working directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.

    """
    
    # The Jaguar output file is opened for the image.    
    f = open('%s%d%s'%(qmjob_prefix, im, '.01.in'), 'r')
    
    # The QM hessian matrix is initialised.
    qhes=[[] for i in range(ntotqm * 3 + 1)]
    
    # While loop iterates through the output file.
    while 1:
        line = f.readline()
        if not line:break
        
        # The hessian matrix is extracted from the Jaguar output file.
        if line[:5] == '&hess':
            while 1:  
                line = f.readline().split()
                if not line:break
                if line == '&':break 
                if len(line) >= 2:
                    lline = len(line)
                    place = (int(line[0]))
                    for k in range(lline - 1):
                        qhes[place].append(line[k + 1])
    f.close()
    
    # The QM hessian matrix is appropriately formatted.
    for i in range(ntotqm * 3):
        qhes[i + 1].insert(0, [])
      
    # The number of QM atoms is initialised.
    nqm = ntotqm - nlink
    
    # If link atoms are included in the frequency calculation, then the hessian must be corrected using the function hescorrect found in this file.
    # Within the function hescorrect, the hessian is appropriately formatted and written to an output file.
    if nlink != 0:
        hescorrect(qhes, ntotqm, nqm, cwd, im)
    
    # If only QM atoms are included in the frequency calculation, then the hessian is appropriately formatted and written to an output file.
    else:
        fh = open(('%s%d'%('qmhes', im)), 'w')
        indx = 0
        for i in range(nqm * 3):
            indx = indx + 1
            hlp = 0
            fh.write(str(indx).rjust(4))
            inum = 0         
            for j in range(i + 1):
                inum = inum + 1
                hlp = hlp + 1
                fh.write(str(qhes[i + 1][j + 1]).rjust(20))
                if (hlp == 4) and (inum != (i+1)):
                    fh.write('\n')
                    fh.write(str(indx).rjust(4))         
                    hlp = 0
            fh.write('\n')
        fh.close()

def qmmolinp(im, qmkey, cwd, ntotqm, qmjob_prefix):
    """
    
    // Function which creates the input file required for a Molpro frequency calculation. //
    // This function is called in the function freqmain for each image. //
    
    Arguments
    ----------
    im : integer
        If using the nudged elastic band, this is the image number.
    qmkey : string
        Can either be None or the user input level of theory.
    cwd : string
        The current working directory.
    ntotqm : integer
        Number of QM atoms including link atoms.
    qmjob_prefix : string
        A prefix used to label the QM job input file.

    """
    
    # The Molpro input file is created for the image.
    fd = open(('%s%d%s'%(qmjob_prefix, im, '.in')), 'w')
    
    # The file qmheader is used to write header details to the input file.
    f = open((cwd + '/qmheader'), 'r')
    for line in f:
        fd.write(line)
    f.close()
    
    # The total number of atoms is written to the input file.
    fd.write(str(ntotqm))
    fd.write('\n')
    fd.write('\n')
    
    # The cartesian coordinates of the QM region are written to the input file.    
    f = open(('%s%d%s'%('qmgeom', im, '.xyz')), 'r') 
    for line in f:
        fd.write(line)
    f.close()
    fd.write('}')
    fd.write('\n')
    
    # The point charges and lattice gradient file names are written to the input file.    
    latinp = '%s%d%s'%('lat', im, '.dat')
    latout = '%s%d%s'%('qmlatgrad', im, '.out')
    latline = 'lattice' + ' ' + '''\'''' + latinp + '''\'''' + ' ' + '''\'''' + latout + '''\''''
    fd.write(latline)
    fd.write('\n')
    fd.write(qmkey.lstrip())
    fd.write('{Frequency;Print,Hessian}')
    fd.close()

def molhesread(im, ntotqm, nlink, cwd, qmjob_prefix):
    """
    
    // Function which extracts the QM hessian from a Molpro output file. //
    // This function is called in the function freqmain for each image. //
    
    Arguments
    ----------
    im : integer
        If using the nudged elastic band, this is the image number.
    ntotqm : integer
        Number of QM atoms including link atoms.
    nlink : integer
        Number of link atoms.
    cwd : string
        The current working directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.

    """
    
    # The number of QM atoms is initialised.
    nqm = ntotqm - nlink
    
    # The Molpro output file is opened for the image.  
    f = open('%s%d%s'%(qmjob_prefix, im, '.out'), 'r')
    
    # The QM hessian matrix is initialised.
    qhes=[[] for i in range(ntotqm * 3 + 1)]
    
    # While loop iterates through the output file.    
    it = ntotqm * 3 
    while 1:
        line = f.readline()
        if not line:break
        
        # The hessian matrix is extracted from the Molpro output file.        
        if line[:61] == ' Force Constants (Second Derivatives of the Energy) in [a.u.]':
            f.readline()
            ll = 0
            irun = 1
            while 1:
                line = f.readline()
                ll = ll + 1
                if len(line.strip()) == 0:break
                nline = line.split()
                lline = len(nline)
                for k in range(lline - 1):
                    if nline[0][1] == 'X':
                        ind = 1
                    if nline[0][1] == 'Y':
                        ind = 2 
                    if nline[0][1] == 'Z':
                        ind = 3 
                    place = (int(nline[0][2:]) - 1) * 3 + ind 
                    qhes[place].append(nline[k + 1])
                if ll == it:
                    f.readline()
                    ll = 5 * irun
                    irun = irun + 1
    f.close()
    
    # The QM hessian matrix is appropriately formatted.    
    for i in range(ntotqm * 3):
        qhes[i + 1].insert(0, [])
        
    # If link atoms are included in the frequency calculation, then the hessian must be corrected using the function hescorrect found in this file.
    # Within the function hescorrect, the hessian is appropriately formatted and written to an output file.     
    if nlink != 0:
        hescorrect(qhes, ntotqm, nqm, cwd, im)
        
    # If only QM atoms are included in the frequency calculation, then the hessian is appropriately formatted and written to an output file.        
    else:
        fh = open(('%s%d'%('qmhes', im)), 'w')
        indx = 0
        for i in range(nqm * 3):
            indx = indx + 1
            hlp = 0
            fh.write(str(indx).rjust(4))
            inum = 0         
            for j in range(i + 1):
                inum = inum + 1
                hlp = hlp + 1
                fh.write(str(qhes[i + 1][j + 1]).rjust(20))
                if (hlp == 4) and (inum != (i + 1)):
                    fh.write('\n')
                    fh.write(str(indx).rjust(4))         
                    hlp = 0
            fh.write('\n')
        fh.close()

def qmgauinp(im, qmkey, cwd, cha_mul, extra_basis, gau_head, qmjob_prefix):
    """
    
    // Function which extracts the QM hessian from a Gaussian output file. //
    // This function is called in the function freqmain for each image. //
    
    Arguments
    ----------
    im : integer
        If using the nudged elastic band, this is the image number.
    qmkey : string
        Can either be None or the user input level of theory.
    cwd : string
        The current working directory.
    cha_mul : list
        The charge and multiplicity given in the format of a list of two numbers (i.e., cha_mul = [0,1]).
    extra_basis : string
        User-specified basis set.
    gau_head : string
        Contains the details for assigning memory and number of processors for a QM job.
    qmjob_prefix : string
        A prefix used to label the QM job input file.

    """
    
    # The Gaussian input file is created for the image.    
    fd = open(('%s%d%s'%(qmjob_prefix, im, '.in')), 'w')
    
    # The variable gau_head is written to the input file, which contains details for assigning memory and processors.    
    if gau_head.lower() != 'none':
        fd.write(gau_head.strip())  
        fd.write('\n')
        
    # The directory for the checkpoint file is written to the input file.
    cwd1 = cwd + '%s%d'%('/image', im) + '%s%d'%('/freq', im)
    chk = '%chk=' + cwd1 + ('%s%s%d%s'%('/', qmjob_prefix, im, '.chk'))
    fd.write(chk)
    fd.write('\n')
    fd.write('\n')
    
    # Using a series of if and else statements, the qmkey (level of theory) is written to the input file.
    # In addition, the method of mixing the guess wavefunction is written along with some other important Gaussian parameters.    
    if qmkey.lower() != 'none':
        if os.path.exists(cwd + '%s%d%s%s%d%s'%('/image', im, '/', qmjob_prefix, im, '.chk')):
            shutil.copy(cwd + '%s%d%s%s%d%s'%('/image', im, '/', qmjob_prefix, im, '.chk'), cwd1)
            jobkey = '# ' + qmkey + ' Nosymm Guess=Read Charge Freq'
        else:
            jobkey = '# ' + qmkey + ' Charge Nosymm Freq'
    else:
        if os.path.exists(cwd + '%s%d%s%s%d%s'%('/image', im, '/', qmjob_prefix, im, '.chk')):
            shutil.copy(cwd + '%s%d%s%s%d%s'%('/image', im, '/', qmjob_prefix, im, '.chk'), cwd1)
            jobkey = '# Guess=Read Nosymm Charge Freq'
        else:
            jobkey = '# Charge Nosymm Freq'
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
    f = open(('%s%d%s'%('qmgeom', im, '.xyz')), 'r') 
    for line in f:
        fd.write(line)
    f.close()
    fd.write('\n')
    
    # The point charge array is written to the input file.    
    f = open(('%s%d%s'%('charges', im, '.xyz')), 'r') 
    for line in f:
        if line.strip() != '&pointch':
            if line.strip() != '&':
                line=line.split()
                fd.write(line[1].rjust(13))
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
    fd.write('\n')
    fd.close()

def gauhesread(im, ntotqm, nlink, cwd, qmjob_prefix):
    """
    
    // Function which extracts the QM hessian from a Gaussian output file. //
    // This function is called in the function freqmain for each image. //
    
    Arguments
    ----------
    im : integer
        If using the nudged elastic band, this is the image number.
    ntotqm : integer
        Number of QM atoms including link atoms.
    nlink : integer
        Number of link atoms.
    cwd : string
        The current working directory.
    qmjob_prefix : string
        A prefix used to label the QM job input file.

    """
    
    # The Gaussian output file is opened for the image.      
    f = open('%s%d%s'%(qmjob_prefix, im, '.log'), 'r')
    
    # The QM hessian matrix is initialised.
    qhes=[[] for i in range(ntotqm * 3 + 1)]
    
    # While loop iterates through the output file.      
    while 1:
        line = f.readline()
        if not line:break
        
        # The hessian matrix is extracted from the Gaussian output file.    
        if line[:23] == ' Test job not archived.':
            line = f.readline()
            resl = line.strip()
            while 1:  
                line = f.readline()
                if not line:break
                resl = resl + line.strip()
    f.close()
    for i in range(len(resl)):
        if resl[i:(i + 5)] == 'NImag':
            res = resl[(i + 6):]
            break
    for i in range(len(res)):
        if res[i:(i + 1)] == '\\':
            resl = res[(i + 2):]
            break   
    res = resl.split(',')
    ip = 0
    for i in range(ntotqm * 3):
        for j in range(i + 1):
            if '\\' in res[ip]:
                hlp = res[ip]
                for k in range(len(res[ip])):
                    if hlp[k:(k + 1)] == '\\':
                        res[ip] = hlp[:(k - 1)]
            qhes[i + 1].append(str(res[ip]))  
            ip = ip + 1
            
    # The QM hessian matrix is appropriately formatted.                
    for i in range(ntotqm * 3):
        qhes[i + 1].insert(0, [])
        
    # The number of QM atoms is initialised.        
    nqm = ntotqm - nlink
    
    # If link atoms are included in the frequency calculation, then the hessian must be corrected using the function hescorrect found in this file.
    # Within the function hescorrect, the hessian is appropriately formatted and written to an output file.     
    if nlink != 0:
        hescorrect(qhes, ntotqm, nqm, cwd, im)
        
    # If only QM atoms are included in the frequency calculation, then the hessian is appropriately formatted and written to an output file.     
    else:
        fh = open(('%s%d'%('qmhes', im)), 'w')
        indx = 0
        for i in range(nqm * 3):
            indx = indx + 1
            hlp = 0
            fh.write(str(indx).rjust(4))
            inum = 0         
            for j in range(i + 1):
                inum = inum + 1
                hlp = hlp + 1
                fh.write(str(qhes[i + 1][j + 1]).rjust(20))
                if (hlp == 4) and (inum != (i + 1)):
                    fh.write('\n')
                    fh.write(str(indx).rjust(4))         
                    hlp = 0
            fh.write('\n')
        fh.close()

def freqmain(qommmadir, qmkey, qmc_job, usrdir, cha_mul, extra_basis, gau_head, prj_freq, qmjob_prefix, qmcode):
    """
    
    // Function which combines the functions in this file for execution of a frequency calculation. //
    // It will prepare the input file, run the frequency calculation, and then extract the results from the output. //
    // This function is called in the main QoMMMa program file qommma.py if the user requests a frequency calculation. //
    
    Arguments
    ----------
    qommmadir : string
        The base QoMMMa directory defined in the file qommma.
    qmkey : string
        Can either be None or the user input level of theory.
    qmc_job : string
        The path to the QM program of the user's choice.
    usrdir : string
        The user directory.
    cha_mul : list
        The charge and multiplicity given in the format of a list of two numbers (i.e., cha_mul = [0,1]) - used only for Gaussian jobs.
    extra_basis : string
        User-specified basis set - used only for Gaussian jobs.  
    gau_head : string
        Contains the details for assigning memory and number of processors for a Gaussian job.
    prj_freq : array
        Used to requency a projected frequency calculation for a given image.
    qmjob_prefix : string
        A prefix used to label the QM job input file.
    qmcode : string
        Represents the user selected QM package.

    """    

    # Initialising the finite difference in the MM hessian.
    delq = 0.0001	
    
    # The current working directory is initialised.
    cwd = os.getcwd()
    
    # The file fortinput is opened and is used to extract relevant information about the system.
    f = open(cwd + '/fortinput', 'r')
    f.readline()
    line = f.readline().split()
    
    # Initialising total number of atoms.
    n = int(line[0])
    
    # Initialising the number of QM atoms.
    nqm = int(line[1])
    
    # Initialising the number of link atoms.
    nlink = int(line[2])
    
    # Initialising the total number of QM atoms, including link atoms.
    ntotqm = nqm + nlink   
    
    # Skipping some lines...
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    f.readline()
    
    # Initialising the number of images, if using the nudged elastic band, that is.
    nimg = int(f.readline().split()[0])
    f.readline()
    
    # The list of QM atom labels is initialised.
    qmlst = []
    for i in range(nqm):
        line = f.readline().split()
        qmlst.append(int(line[0]))
    f.close()

    if nlink > 1:
        # The list qm_lk contains details of which QM atom is connected to which link atom.
        qm_lk = []
        
        # By reading the file link_details, the QM atom connected to a link atom is read and added to the list qm_lk.
        ff = open(cwd + '/link_details')
        for i in range(nlink):
            line = ff.readline().split()
            qm_lk.append(int(line[0]))
        ff.close() 
        for i in range(nlink - 1):
            j = i
            for k in range(nlink - 1 - i):  
                j += 1  
                if qm_lk[i] > qm_lk[j]:
                    qomend('Error, link atoms details should be in increasing order with respect to QM atom number, In Link key give first details for link atom connected with QM atom ' + str(qm_lk[j]) + '  then for QM atom  ' + str(qm_lk[i]), cwd, usrdir)
    
    # The initial arrangements for a frequency calculation are performed.
    try:
        initial(nimg, cwd, qommmadir)
    except:
        qomend('Error: Initial setup program for QoMMMa frequency job failed', cwd, usrdir) 
   
    # The frequency calculation is performed on every image using a for-loop.
    im = 0
    for i in range(nimg):
        im = im + 1 
        qomlog('About to start MM Hessian calculation for QoMMMa Frequency job for image ' + str(im), usrdir)
        dst = cwd + ('%s%d'%('/image', im)) + ('%s%d'%('/freq', im))
        os.chdir(dst)
        
        # Tinker coordinate files required for the frequency calculation are created.
        try:
            mminpco(im, qmlst, delq)
        except:
            qomend('Error, while creating input geometry for mm gradient job for image : ' + str(im), cwd, usrdir) 
         
        # For the generated Tinker coordinate files, the gradient is obtained.
        try:
            analgrad(im, qmlst, qommmadir)
        except:
            qomend('Error, while running mm gradient job for image : ' + str(im), cwd, usrdir)
        
        # From the calculated gradients, the MM hessians are created.
        try:
            mmhes(im, qmlst, n, delq)
        except:
            qomend('Error, while calculating mm hessian for image : ' + str(im), cwd, usrdir)
     
    
        # Preparing input files for QM frequency calculation, which uses different functions depending on the QM package.
        qomlog('About to start QM frequency calculation for image : ' + str(im), usrdir)    
        try: 
            if qmcode.lower() == 'jaguar':
                qmjaginp(im, qmkey, cwd, qmjob_prefix)
            elif qmcode.lower() == 'gaussian': 
                qmgauinp(im, qmkey, cwd, cha_mul, extra_basis, gau_head, qmjob_prefix)
            elif qmcode.lower() == 'molpro':
                # The function latdatin needs to be used to prepare for a Molpro job.
                from molutil import latdatin 
                latdatin(im) 
                qmmolinp(im,qmkey,cwd,ntotqm,qmjob_prefix)
        except:
            qomend('Error, while creating input file for QM job for image : ' + str(im), cwd, usrdir)
            
        # Performing QM frequency calculation.
        try:
            os.system(qmc_job + ' ' + ('%s%d%s'%(qmjob_prefix, im, '.in')))
        except:
            qomend('Error while running QM frequency job for image :' + str(im), cwd, usrdir)
            
        # The QM hessian matrix is extracted from the output file, which uses different functions depending on the QM package.    
        try:
            if qmcode.lower() == 'jaguar':
                jaghesread(im, ntotqm, nlink, cwd, qmjob_prefix)
            elif qmcode.lower() == 'gaussian':
                gauhesread(im, ntotqm, nlink, cwd, qmjob_prefix)
            else:
                molhesread(im, ntotqm, nlink, cwd, qmjob_prefix)
        except:
            qomend('Error while reading QM hessian from QM output for image :' + str(im), cwd, usrdir)
        
        # The input file required for execution of the Fortran code is created.
        try:
            freqinp(im, qmlst, cwd)
        except:
            qomend('Error while preparing hessian_input for image :' + str(im), cwd, usrdir)
        
        # Depending on the user input, either a projected or normal frequency calculation is performed.
        # This is the point where the Fortran code is used, so it should be referenced to understand the frequency calculation.
        if im in prj_freq:
            qomlog('Note, For image ' + str(im) + 'projected frequency calculation is requested', usrdir)
            try:
                os.system(qommmadir + '/bin/prj_freq.x')     
            except:
                qomend('Error while running QoMMMa projected frequency program for image :' + str(im), cwd, usrdir)
        else:
            try:
                 os.system(qommmadir + '/bin/freq.x')     
            except:
                qomend('Error while running QoMMMa Frequency program for image :' + str(im), cwd, usrdir)
                
        # The results from the frequency calculation are rearranged.         
        try:
            shutil.copy('freqs.xyz', ('%s%d%s'%('freqs', im, '.xyz')))  
            shutil.copy('freqs_results', ('%s%d'%('freqs_results', im)))  
            os.remove('freqs.xyz')
            os.remove('freqs_results')
            shutil.copy(('%s%d%s'%('freqs', im, '.xyz')), usrdir + ('%s%d%s'%('/freqs', im, '.xyz')))
            shutil.copy(('%s%d'%('freqs_results', im)), usrdir + ('%s%d'%('/freqs_results', im)))
            os.chdir(cwd)
        except:
            qomend('Error in QoMMMa Frequency program for image :' + str(im), cwd, usrdir)
            
    # Frequency calculation complete!
    qomend('Congratulations, Frequency job for all images completed', cwd, usrdir)
