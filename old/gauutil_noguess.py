#!/usr/bin/python
import os
import sys
import shutil
from qomutil import qomend


"""This python program file is used to generate Gauusian input, to run the
Gaussian job and to extract the results from Gaussian output. It consists of the following routines. Here qomend routine is imported from qomutil to print the error message and to stop the program if error occurs."""

def gauinp(imn,cwd,qmjob_prefix,qmkey,cha_mul,extra_basis,gau_head):
  """This routine is to prepare the input file for Gaussian job. It takes arguments: imn (integer, the image number); cwd (string, current working directory), qmjob_prefix (string, prefix for QM job input file), qmkey (string, either None or user given level of theory), cha_mul is the charge and multiplicity given as list of two numbers (default is cha_mul=[0,1]), extra_basis is user specified basis set as a string (default='None') and gau_head is used to assign memory and number of processors for Gaussian jobs, given as a string (default='None'). This routine is called from qmgaumain routine for each cycle of QoMMMa."""
  fd=open(('%s%d%s'%(qmjob_prefix,imn,'.in')),'w') 
  if gau_head.lower()!='none':
      fd.write(gau_head.strip())  
      fd.write('\n')
  cwd1=cwd+'%s%d'%('/image',imn)
  chk='%chk='+cwd1+('%s%s%d%s'%('/',qmjob_prefix,imn,'.chk'))
  fd.write(chk)
  fd.write('\n')
  fd.write('\n')
  if qmkey.lower()!='none':
     if os.path.exists('%s%d%s'%(qmjob_prefix,imn,'.chk')):
          jobkey='# '+qmkey+' Nosymm Guess=Read Charge Force Prop=Grid'
     else:
          jobkey='# '+qmkey+' Charge Nosymm Force Prop=Grid'
  else:
     if os.path.exists('%s%d%s'%(qmjob_prefix,imn,'.chk')):
          jobkey='# Guess=Read Nosymm Charge Force Prop=Grid'
     else:
          jobkey='# Charge Nosymm Force Prop=Grid'
  fd.write(jobkey)
  fd.write('\n')
  fd.write('\n')
  fd.write('TITLE')
  fd.write('\n')
  fd.write('\n')
  fd.write(str(cha_mul[0]).rjust(3))  
  fd.write(str(cha_mul[1]).rjust(5))  
  fd.write('\n') 	  
  f=open(('%s%d%s'%('qmgeom',imn,'.xyz')),'r') 
  for line in f:
     fd.write(line)
  f.close()
  fd.write('\n')
  f=open(('%s%d%s'%('charges',imn,'.xyz')),'r') 
  fresp='%s%d%s'%('points',imn,'gau.pts')
  numresp=len(open(fresp,'r').readlines())
  for line in f:
    if line.strip() != '&pointch':
      if line.strip() !='&':
         line=line.split()
         fd.write(line[1].rjust(13))
         fd.write(line[2].rjust(13))
         fd.write(line[3].rjust(13))
         fd.write(line[0].rjust(10))
         fd.write('\n')
  f.close()
  fd.write('\n')
  if extra_basis.lower()!='none':
     fd.write(extra_basis.lstrip())
     fd.write('\n')
  fd.write(str(numresp))
  fd.write(str(2).rjust(3))
  fd.write(str(20).rjust(4))
  fd.write(str(21).rjust(4))
  fd.write('\n')
  fd.write('\n')
  fd.close()
  shutil.copy(('%s%d%s'%('points',imn,'gau.pts')),'fort.20')

def gauout(fi,l,nqm,nlink,fortf):
  '''This program is to extract energy, forces and Mulliken charges from Gaussian output. The arguments are: fi (string, Gaussian job output file name); l (integer, image number); nqm (integer, number of QM atoms); nlink (integer, number of link atoms). Gradients are calculated by multiplying forces with -1. Energy and gradients are stored in file ab_initio* and Mulliken charges is stored in file mulliken*. The gradients corresponding to MM point charges is calculated from electric field given in Gaussian support file (here file name is given as string 'fortf'), and are stored in a file 'qmlatgrad*.out'.'''
  fd=open(('%s%d'%('ab_initio',l)),'a')		# here 'a' is used, useful for mecp
  fd.write('Energy')
  fd.write('\n')
  fdm=open(('%s%d'%('mulliken',l)),'a')
  f=open(fi,'r')
  ntot=nqm+nlink
  while 1:
      line=f.readline()
      if not line:break
      if line[:29]==' Self energy of the charges =':
          pointE=float(line.split()[6])
      elif line[:42]==' >>>>>>>>>> Convergence criterion not met.':
          line=f.readline()
          line=f.readline()
      elif line[:9]==' SCF Done':
        totE=float(line.split()[4])
        ene=totE-pointE 
        fd.write(str(ene))
        fd.write('\n')
      elif line[:18]==' Mulliken charges:':
         f.readline()
         for i in range(ntot):
            line=f.readline().split()
            fdm.write(line[2].rjust(10))
            fdm.write('\n')  
         line=f.readline().split()
         fdm.write('total charge')
         fdm.write('\n')
         fdm.write(line[5].rjust(10)) 
         fdm.write('\n')
         fdm.close()
      elif line[:37]==' Mulliken charges and spin densities:':
         f.readline()
         for i in range(ntot):
            line=f.readline().split()
            fdm.write(line[2].rjust(10))
            fdm.write('\n')
         line=f.readline().split()
         fdm.write('total charge')
         fdm.write('\n')
         fdm.write(line[5].rjust(10))
         fdm.write('\n')
         fdm.close()
      elif line[:59]==' Center     Atomic                   Forces (Hartrees/Bohr)': 
        fd.write('Gradient')
        fd.write('\n')
        f.readline()
        f.readline()
        for i in range(ntot):
	   line=f.readline().split()
           fd.write(line[0].rjust(3))
           linex=float(line[2])*-1
           fd.write(str(linex).rjust(20))
           liney=float(line[3])*-1
           fd.write(str(liney).rjust(20))
           linez=float(line[4])*-1
           fd.write(str(linez).rjust(20))
           fd.write('\n')  
        fd.close()
  f.close()
  ftc=open('%s%d%s'%('points_chg',l,'.pts'),'r')
  fg=open('%s%d%s'%('qmlatgrad',l,'.out'),'a')
  ff=open(fortf,'r') 
  cha=[]
  for line in ftc:
     tcha=float(line.split()[0])
     cha.append(tcha)
  ftc.close()
  fg.write(' gradient output')
  fg.write('\n')
  fg.write(str(len(cha)))
  fg.write('\n')
  for i in range(len(cha)):
    ff.readline()
    linef=ff.readline().split()               
    for j in range(3):
       grad=-(float(cha[i])*float(linef[j]))
       fg.write(str(grad).rjust(25))
    fg.write('\n')
  fg.close()

def qmgaumain(imn,cwd,usrdir,qmjob_prefix,qmgau_job,nqm,nlink,cln,qmkey,cha_mul,extra_basis,gau_head):
    '''This part of program is to prepare the input file, and then to run the Gaussian and to extract the results from output. It uses two subprograms (gauinp, gauout). It takes arguments: imn (integer, the image number); cwd, (string, current working directory), usrdir (string, user directory), qmjob_prefix (string, prefix for QM job input file), qmgau_job (is string, representing Gaussian path), nqm (integer, number of QM atoms); nlink (integer, number of link atoms), cln (integer, QoMMMa cycle number), qmkey (string, any extra keys specified by user for Gaussian job, default='None', i.e.'HF/STO-3G'), cha_mul is the charge and multiplicity given as list of two numbers (default is cha_mul=[0,1]), extra_basis is user specified basis set as a string (default='None'), and gau_head is used to assign memory and number of processors for Gaussian jobs, given as a string (default='None'). This main routine is called from QoMMMa main program at each cycle. Before preparing new input files this routine will rename all Gaussian files created in previous cycle.''' 
    if os.path.exists(cwd+('%s%d'%('/image',imn))+('%s%s%d%s'%('/',qmjob_prefix,imn,'.in'))):
     try:			# rearrangeing previous cycle Gaussian files
       shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'.in')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'.in'))) 
       shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'.log')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'.log')))
       os.remove('%s%d%s'%(qmjob_prefix,imn,'.in'))
       os.remove('%s%d%s'%(qmjob_prefix,imn,'.log'))
     except:
       qomend('ERROR, while rearrangeing previous cycle Gaussian files for image : '+str(imn),cwd,usrdir)
    try: 				# preparing input
       gauinp(imn,cwd,qmjob_prefix,qmkey,cha_mul,extra_basis,gau_head)	
    except:
       qomend('Error while creating Gaussian input for image : '+str(imn),cwd,usrdir)
    try:			               		# Gaussian job
       os.system(qmgau_job+' '+('%s%d%s'%(qmjob_prefix,imn,'.in')))
    except:
       qomend('Error while running Gaussian job at cycle: '+str(cln)+' image :'+str(imn),cwd,usrdir)
    try:				# extracting QM results
       fi='%s%d%s'%(qmjob_prefix,imn,'.log')
       gauout(fi,imn,nqm,nlink,'fort.21')	
       shutil.copy(('%s%d'%('ab_initio',imn)),cwd)
       shutil.copy(('%s%d'%('mulliken',imn)),cwd)
       shutil.copy(('%s%d%s'%('qmlatgrad',imn,'.out')),cwd) 
       os.remove('%s%d%s'%('qmlatgrad',imn,'.out'))
       os.remove('%s%d'%('ab_initio',imn))
       os.remove('%s%d'%('mulliken',imn))
    except:
       qomend('ERROR, while reading Energy or Gradient or Mulliken charge or electrostatic potential from QM output file of image : '+str(imn),cwd,usrdir)    

def mgauinp(imn,cwd,qmjob_prefix,qmkey,cha_mul1,cha_mul2,extra_basis,gau_head):
  """This routine is to prepare the input file for Gaussian jobin MECP QoMMMa job. It will create two input files for two states A and B. It takes arguments: imn (integer, the image number); cwd (string, current working directory), qmjob_prefix (string, prefix for QM job input file), qmkey (string, either None or user given level of theory), cha_mul1 and cha_mul2 are the charge and multiplicity of states A and B, given as list, extra_basis is user specified basis set as a string (default='None') and gau_head is used to assign memory and number of processors for Gaussian jobs, given as a string (default='None'). This routine is called from mecp_gaumain routine for each cycle of QoMMMa."""
  fda=open(('%s%d%s'%(qmjob_prefix,imn,'_A.in')),'w')  
  fdb=open(('%s%d%s'%(qmjob_prefix,imn,'_B.in')),'w')  
  if gau_head.lower()!='none':
      fda.write(gau_head.strip())  
      fda.write('\n')
      fdb.write(gau_head.strip())  
      fdb.write('\n')
  cwd1=cwd+'%s%d'%('/image',imn)
  chk='%chk='+cwd1+('%s%s%d%s'%('/',qmjob_prefix,imn,'_A.chk'))
  fda.write(chk)
  chk='%chk='+cwd1+('%s%s%d%s'%('/',qmjob_prefix,imn,'_B.chk'))
  fdb.write(chk)
  fda.write('\n')
  fda.write('\n')
  fdb.write('\n')
  fdb.write('\n')
  if qmkey.lower()!='none':
     if os.path.exists('%s%d%s'%(qmjob_prefix,imn,'_A.chk')):
          jobkey='# '+qmkey+' Nosymm Guess=Read Charge Force Prop=Grid'
     else:
          jobkey='# '+qmkey+' Charge Nosymm Force Prop=Grid'
  else:
     if os.path.exists('%s%d%s'%(qmjob_prefix,imn,'_A.chk')):
          jobkey='# Guess=Read Nosymm Charge Force Prop=Grid'
     else:
          jobkey='# Charge Nosymm Force Prop=Grid'
  fda.write(jobkey)
  fda.write('\n')
  fda.write('\n')
  if qmkey.lower()!='none':
     if os.path.exists('%s%d%s'%(qmjob_prefix,imn,'_B.chk')):
          jobkey='# '+qmkey+' Nosymm Guess=Read Charge Force Prop=Grid'
     else:
          jobkey='# '+qmkey+' Charge Nosymm Force Prop=Grid'
  else:
     if os.path.exists('%s%d%s'%(qmjob_prefix,imn,'_B.chk')):
          jobkey='# Guess=Read Nosymm Charge Force Prop=Grid'
     else:
          jobkey='# Charge Nosymm Force Prop=Grid'
  fdb.write(jobkey)
  fdb.write('\n')
  fdb.write('\n')
  fda.write('TITLE')
  fda.write('\n')
  fda.write('\n')
  fdb.write('TITLE')
  fdb.write('\n')
  fdb.write('\n')
  fda.write(str(cha_mul1[0]).rjust(3))  
  fda.write(str(cha_mul1[1]).rjust(5))  
  fda.write('\n') 	  
  fdb.write(str(cha_mul2[0]).rjust(3))  
  fdb.write(str(cha_mul2[1]).rjust(5))  
  fdb.write('\n') 	  
  f=open(('%s%d%s'%('qmgeom',imn,'.xyz')),'r') 
  for line in f:
     fda.write(line)
     fdb.write(line)
  f.close()
  fda.write('\n')
  fdb.write('\n')
  f=open(('%s%d%s'%('charges',imn,'.xyz')),'r') 
  ff1=open('fort.20','w')
  ff2=open('fort.30','w')
  ill=0
  for line in f:
    if line.strip() != '&pointch':
      if line.strip() !='&':
         line=line.split()
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
         ill+=1
  f.close()
  ff1.close()
  ff2.close()
  fda.write('\n')
  fdb.write('\n')
  if extra_basis.lower()!='none':
     fda.write(extra_basis.lstrip())
     fda.write('\n')
     fdb.write(extra_basis.lstrip())
     fdb.write('\n')
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



def mecp_gaumain(imn,cwd,usrdir,qmjob_prefix,qmgau_job,nqm,nlink,cln,qmkey,cha_mul1,cha_mul2,extra_basis,gau_head):
  '''This part of program is to prepare the input file, and then to run the Gaussian and to extract the results from output. It uses two subprograms (mgauinp, gauout). It takes arguments: imn (integer, the image number); cwd, (string, current working directory), usrdir (string, user directory), qmjob_prefix (string, prefix for QM job input file), qmgau_job (is string, representing Gaussian path), nqm (integer, number of QM atoms); nlink (integer, number of link atoms), cln (integer, QoMMMa cycle number), qmkey (string, any extra keys specified by user for Gaussian job, default='None', i.e.'HF/STO-3G'), cha_mul1 and cha_mul2 are the charge and multiplicity of two states in MECP, given as list of two numbers, extra_basis is user specified basis set as a string (default='None'), and gau_head is used to assign memory and number of processors for Gaussian jobs, given as a string (default='None'). This main routine is called from QoMMMa main program at each cycle if MECP is requested. Before preparing new input files this routine will rename all Gaussian files created in previous cycle.''' 
  if os.path.exists(cwd+('%s%d'%('/image',imn))+('%s%s%d%s'%('/',qmjob_prefix,imn,'_A.in'))):
   try:			#rearrangeing previous cycle files
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_A.in')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.in'))) 
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_A.log')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.log')))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_A.in'))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_A.log'))
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_B.in')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.in'))) 
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_B.log')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.log')))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_B.in'))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_B.log'))
   except:
    qomend('error while rearrangeing Gaussian files created in previous run at cycle :'+str(cln)+' for image: '+str(imn),cwd,usrdir)
  try:			# preparing input file for both states 
    mgauinp(imn,cwd,qmjob_prefix,qmkey,cha_mul1,cha_mul2,extra_basis,gau_head)
  except:
    qomend('Error while trying to prepare Gaussian input files at cycle : '+str(cln)+' for image : '+str(imn),cwd,usrdir)
  try:			               		# Gaussian job
    os.system(qmgau_job+' '+('%s%d%s'%(qmjob_prefix,imn,'_A.in')))
  except:
    qomend('Error while running Gaussian job for state A at cycle : '+str(cln)+' for image : '+str(imn),cwd,usrdir)
  try:			               		# Gaussian job
    os.system(qmgau_job+' '+('%s%d%s'%(qmjob_prefix,imn,'_B.in')))
  except:
    qomend('Error while running Gaussian job for state B at cycle : '+str(cln)+' for image : '+str(imn),cwd,usrdir)
  try:				# extracting QM results
    fi='%s%d%s'%(qmjob_prefix,imn,'_A.log')
    gauout(fi,imn,nqm,nlink,'fort.21')	
    fi='%s%d%s'%(qmjob_prefix,imn,'_B.log')
    gauout(fi,imn,nqm,nlink,'fort.31')	
    shutil.copy(('%s%d'%('ab_initio',imn)),cwd)
    shutil.copy(('%s%d'%('mulliken',imn)),cwd)
    shutil.copy(('%s%d%s'%('qmlatgrad',imn,'.out')),cwd) 
    os.remove('%s%d%s'%('qmlatgrad',imn,'.out'))
    os.remove('%s%d'%('ab_initio',imn))
    os.remove('%s%d'%('mulliken',imn))
  except:
    qomend('Error while extracting results from Gaussian output files at cycle : '+str(cln)+' for image : '+str(imn),cwd,usrdir)

