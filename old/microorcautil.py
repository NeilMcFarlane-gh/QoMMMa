#!/usr/bin/python
import os
import sys
import shutil
from qomutil import qomend


"""This python program file is used to generate ORCA input, to run the
ORCA job and to extract the results from ORCA output. It consists of the following routines. Here qomend routine is imported from qomutil to print the error message and to stop the program if error occurs."""


def microorcainp(imn,cwd,qmjob_prefix,qmkey,cha_mul,extra_basis,orca_head):
  """This routine is to prepare the input file for the microoptimization job for ORCA. It takes arguments: imn (integer, the image number); cwd (string, current working directory), qmjob_prefix (string, prefix for QM job input file), qmkey (string, either None or user given level of theory), cha_mul is the charge and multiplicity given as list of two numbers (default is cha_mul=[0,1]), extra_basis is user specified basis set as a string (default='None') and orca_head is used to assign memory and number of processors for ORCA jobs, given as a string (default='None'). This routine is called from qmorcamain routine for each cycle of QoMMMa."""
  fd=open(('%s%s%d%s'%('micro',qmjob_prefix,imn,'.in')),'w') 
  if orca_head.lower()!='none':
      fd.write(orca_head.strip())  
      fd.write('\n')
  cwd1=cwd+'%s%d'%('/image',imn)
  if qmkey.lower()!='none':
     if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'.gbw')):
          jobkey='!'+qmkey+' moread looseopt xyzfile'
     else:
          jobkey='!'+qmkey+' looseopt xyzfile'
  else:
     if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'.gbw')):
          jobkey='! moread looseopt xyzfile'
     else:
          jobkey='! looseopt xyzfile'
  fd.write(jobkey)
  fd.write('\n')
  ptc='%pointcharges "'+'%s%d%s'%('ptchg',imn,'.pc"')
  fd.write(ptc)
  fd.write('\n')
  if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'.gbw')):    
      wavefunction='%moinp '+'"'+'%s%s%d%s'%('old_',qmjob_prefix,imn,'.gbw'+'"')
      fd.write(wavefunction)
      fd.write('\n')
  fd.write(extra_basis.strip())
  fd.write('\n')
  fd.write('\n')
  fd.write('* xyz ')
  fd.write(str(cha_mul[0]).rjust(3))  
  fd.write(str(cha_mul[1]).rjust(5))  
  fd.write('\n') 	  
  f=open(('%s%d%s'%('qmgeom',imn,'.xyz')),'r') 
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
  fg=open(('%s%d%s'%('ptchg',imn,'.pc')),'a')
  f=open(('%s%d%s'%('charges',imn,'.xyz')),'r') 
  fresp='%s%d%s'%('points',imn,'gau.pts')
  numresp=len(open(fresp,'r').readlines())
  fg.write(str(numresp))
  fg.write('\n')
  for line in f:
    if line.strip() != '&pointch':
      if line.strip() !='&':
         line=line.split()
         fg.write(line[0].rjust(10))
         fg.write(line[1].rjust(13))
         fg.write(line[2].rjust(13))
         fg.write(line[3].rjust(13))
         fg.write('\n')
  f.close()
  fg.close()




def orcainp(imn,cwd,qmjob_prefix,qmkey,cha_mul,extra_basis,orca_head):
  """This routine is to prepare the input file for ORCA job. It takes arguments: imn (integer, the image number); cwd (string, current working directory), qmjob_prefix (string, prefix for QM job input file), qmkey (string, either None or user given level of theory), cha_mul is the charge and multiplicity given as list of two numbers (default is cha_mul=[0,1]), extra_basis is user specified basis set as a string (default='None') and orca_head is used to assign memory and number of processors for ORCA jobs, given as a string (default='None'). This routine is called from qmorcamain routine for each cycle of QoMMMa."""
  fd=open(('%s%d%s'%(qmjob_prefix,imn,'.in')),'w') 
  if orca_head.lower()!='none':
      fd.write(orca_head.strip())  
      fd.write('\n')
  cwd1=cwd+'%s%d'%('/image',imn)
  if qmkey.lower()!='none':
     if os.path.exists('%s%s%d%s'%('micro',qmjob_prefix,imn,'.gbw')):
          jobkey='!'+qmkey+' moread engrad'
     else:
          jobkey='!'+qmkey+' engrad'
  else:
     if os.path.exists('%s%s%d%s'%('micro',qmjob_prefix,imn,'.gbw')):
          jobkey='! moread engrad'
     else:
          jobkey='! engrad'
  fd.write(jobkey)
  fd.write('\n')
  ptc='%pointcharges "'+'%s%d%s'%('ptchg',imn,'.pc"')
  fd.write(ptc)
  fd.write('\n')
  if os.path.exists('%s%s%d%s'%('micro',qmjob_prefix,imn,'.gbw')):    
      wavefunction='%moinp '+'"'+'%s%s%d%s'%('micro',qmjob_prefix,imn,'.gbw'+'"')
      fd.write(wavefunction)
      fd.write('\n')
  fd.write('\n')
  fd.write('* xyz ')
  fd.write(str(cha_mul[0]).rjust(3))  
  fd.write(str(cha_mul[1]).rjust(5))  
  fd.write('\n') 	  
  f=open(('%s%s%d%s'%('micro',qmjob_prefix,imn,'.xyz')),'r') 
  f.readline()
  f.readline()
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

def orcaout(fi,l,nqm,nlink,pcg):
  '''This program is to extract energy, gradient and Mulliken charges from ORCA output. The arguments are: fi (string, ORCA job output file name); l (integer, image number); nqm (integer, number of QM atoms); nlink (integer, number of link atoms), pcg (MM gradients calculated by ORCA). Energy and gradients are stored in file ab_initio* and Mulliken charges is stored in file mulliken*. The gradients corresponding to MM point charges are written by ORCA in the file '.pcgrad' and stored in a file 'qmlatgrad*.out'.'''
  fd=open(('%s%d'%('ab_initio',l)),'a')		# here 'a' is used, useful for mecp
  fd.write('Energy')
  fd.write('\n')
  fdm=open(('%s%d'%('mulliken',l)),'a')
  f=open(fi,'r')
  ntot=nqm+nlink
  while 1:
      line=f.readline()
      if not line:break
      if line[:25]=='FINAL SINGLE POINT ENERGY':
          qmE=float(line.split()[4])
          fd.write(str(qmE))
          fd.write('\n')
      elif line[:23]=='MULLIKEN ATOMIC CHARGES':
         f.readline()
         for i in range(ntot):
            line=f.readline().split()
            fdm.write(line[3].rjust(10))
            fdm.write('\n')  
         line=f.readline().split()
         fdm.write('total charge')
         fdm.write('\n')
         fdm.write(line[4].rjust(10)) 
         fdm.write('\n')
         fdm.close()
      elif line[:44]=='MULLIKEN ATOMIC CHARGES AND SPIN POPULATIONS':
         f.readline()
         for i in range(ntot):
            line=f.readline().split()
            fdm.write(line[3].rjust(10))
            fdm.write('\n')
         line=f.readline().split()
         fdm.write('total charge')
         fdm.write('\n')
         fdm.write(line[5].rjust(10))
         fdm.write('\n')
         fdm.close()
      elif line[:18]=='CARTESIAN GRADIENT': 
        fd.write('Gradient')
        fd.write('\n')
        f.readline()
        f.readline()
        for i in range(ntot):
	   line=f.readline().split()
           fd.write(line[0].rjust(3))
           linex=float(line[3])
           fd.write(str(linex).rjust(20))
           liney=float(line[4])
           fd.write(str(liney).rjust(20))
           linez=float(line[5])
           fd.write(str(linez).rjust(20))
           fd.write('\n')  
        fd.close()
  f.close()
  ftc=open('%s%d%s'%('points_chg',l,'.pts'),'r')
  fg=open('%s%d%s'%('qmlatgrad',l,'.out'),'a')
  ff=open(pcg,'r') 
  cha=[]
  for line in ftc:
     tcha=float(line.split()[0])
     cha.append(tcha)
  ftc.close()
  fg.write(' gradient output')
  fg.write('\n')
  fg.write(str(len(cha)))
  fg.write('\n')
  ff.readline()
  for line in ff:
     fg.write(line)
  f.close()
  fg.close()

def qmorcamain(imn,cwd,usrdir,qmjob_prefix,nqm,nlink,cln,qmkey,cha_mul,extra_basis,orca_head):
    '''This part of program is to prepare the input file, and then to run ORCA and to extract the results from output. It uses two subprograms (orcainp, orcaout). It takes arguments: imn (integer, the image number); cwd, (string, current working directory), usrdir (string, user directory), qmjob_prefix (string, prefix for QM job input file), nqm (integer, number of QM atoms); nlink (integer, number of link atoms), cln (integer, QoMMMa cycle number), qmkey (string, any extra keys specified by user for ORCA job, default='None', i.e.'HF/STO-3G'), cha_mul is the charge and multiplicity given as list of two numbers (default is cha_mul=[0,1]), extra_basis is user specified basis set as a string (default='None'), and orca_head is used to assign memory and number of processors for ORCA jobs, given as a string (default='None'). This main routine is called from QoMMMa main program at each cycle. Before preparing new input files this routine will rename all ORCA files created in previous cycle.''' 
    if os.path.exists(cwd+('%s%d'%('/image',imn))+('%s%s%d%s'%('/',qmjob_prefix,imn,'.in'))):
     try:			# rearrangeing previous cycle ORCA files
       shutil.copy(('%s%s%d%s'%('micro',qmjob_prefix,imn,'.in')),('%s%s%d%s'%('old_micro',qmjob_prefix,imn,'.in'))) 
       shutil.copy(('%s%s%d%s'%('micro',qmjob_prefix,imn,'.log')),('%s%s%d%s'%('old_micro',qmjob_prefix,imn,'.log')))
       shutil.copy(('%s%s%d%s'%('micro',qmjob_prefix,imn,'.gbw')),('%s%s%d%s'%('old_micro',qmjob_prefix,imn,'.gbw')))
       shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'.in')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'.in'))) 
       shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'.log')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'.log')))
       shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'.pcgrad')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'.pcgrad')))
       shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'.gbw')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'.gbw')))
       shutil.copy(('%s%d%s'%('ptchg',imn,'.pc')),('%s%d%s'%('old_ptchg',imn,'.pc')))
       os.remove('%s%s%d%s'%('micro',qmjob_prefix,imn,'.in'))
       os.remove('%s%s%d%s'%('micro',qmjob_prefix,imn,'.log'))
       os.remove('%s%s%d%s'%('micro',qmjob_prefix,imn,'.gbw'))
       os.remove('%s%d%s'%(qmjob_prefix,imn,'.in'))
       os.remove('%s%d%s'%(qmjob_prefix,imn,'.log'))
       os.remove('%s%d%s'%(qmjob_prefix,imn,'.pcgrad'))
       os.remove('%s%d%s'%(qmjob_prefix,imn,'.gbw'))
       os.remove('%s%d%s'%('ptchg',imn,'.pc'))
     except:
       qomend('ERROR, while rearrangeing previous cycle ORCA files for image : '+str(imn),cwd,usrdir)
    try: 				# preparing input for microoptimization
       microorcainp(imn,cwd,qmjob_prefix,qmkey,cha_mul,extra_basis,orca_head)	
    except:
       qomend('Error while creating ORCA microoptimization input for image : '+str(imn),cwd,usrdir)
    try:			               		# ORCA microoptimization job
       os.system('orca '+('%s%s%d%s'%('micro',qmjob_prefix,imn,'.in'))+' > '+('%s%s%d%s'%('micro',qmjob_prefix,imn,'.log')))
    except:
       qomend('Error while running ORCA microoptimization job at cycle: '+str(cln)+' image :'+str(imn),cwd,usrdir)
    try: 				# preparing input for energy and gradient calculation
       orcainp(imn,cwd,qmjob_prefix,qmkey,cha_mul,extra_basis,orca_head)	
    except:
       qomend('Error while creating ORCA engrad input for image : '+str(imn),cwd,usrdir)
    try:			               		# ORCA engrad job
       os.system('orca '+('%s%d%s'%(qmjob_prefix,imn,'.in'))+' > '+('%s%d%s'%(qmjob_prefix,imn,'.log')))
    except:
       qomend('Error while running ORCA engrad job at cycle: '+str(cln)+' image :'+str(imn),cwd,usrdir)
    try:				# extracting QM results
       fi='%s%d%s'%(qmjob_prefix,imn,'.log')
       pcg='%s%d%s'%(qmjob_prefix,imn,'.pcgrad')
       orcaout(fi,imn,nqm,nlink,pcg)	
       shutil.copy(('%s%d'%('ab_initio',imn)),cwd)
       shutil.copy(('%s%d'%('mulliken',imn)),cwd)
       shutil.copy(('%s%d%s'%('qmlatgrad',imn,'.out')),cwd) 
       os.remove('%s%d%s'%('qmlatgrad',imn,'.out'))
       os.remove('%s%d'%('ab_initio',imn))
       os.remove('%s%d'%('mulliken',imn))
    except:
       qomend('ERROR, while reading Energy or Gradient or Mulliken charge or electrostatic potential from QM output file of image : '+str(imn),cwd,usrdir)    

def morcainp(imn,cwd,qmjob_prefix,qmkey,cha_mul1,cha_mul2,extra_basis,orca_head):
  """This routine is to prepare the input file for ORCA job in MECP QoMMMa job. It will create two input files for two states A and B. It takes arguments: imn (integer, the image number); cwd (string, current working directory), qmjob_prefix (string, prefix for QM job input file), qmkey (string, either None or user given level of theory), cha_mul1 and cha_mul2 are the charge and multiplicity of states A and B, given as list, extra_basis is user specified basis set as a string (default='None') and orca_head is used to assign memory and number of processors for ORCA jobs, given as a string (default='None'). This routine is called from mecp_orcamain routine for each cycle of QoMMMa."""
  fda=open(('%s%d%s'%(qmjob_prefix,imn,'_A.in')),'w')  
  fdb=open(('%s%d%s'%(qmjob_prefix,imn,'_B.in')),'w')  
  if orca_head.lower()!='none':
      fda.write(orca_head.strip())  
      fda.write('\n')
      fdb.write(orca_head.strip())  
      fdb.write('\n')
  cwd1=cwd+'%s%d'%('/image',imn)
  if qmkey.lower()!='none':
     if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.gbw')):
          jobkey='!'+qmkey+' moread engrad'
     else:
          jobkey='!'+qmkey+' engrad'
  else:
     if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.gbw')):
          jobkey='!'+qmkey+' moread engrad'
     else:
          jobkey='!'+qmkey+' engrad'
  fda.write(jobkey)
  fda.write('\n')
  if qmkey.lower()!='none':
     if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.gbw')):
          jobkey='!'+qmkey+' moread engrad'
     else:
          jobkey='!'+qmkey+' engrad'
  else:
     if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.gbw')):
          jobkey='!'+qmkey+' moread engrad'
     else:
          jobkey='!'+qmkey+' engrad'
  fdb.write(jobkey)
  fdb.write('\n')
  ptca='%pointcharges "'+'%s%d%s'%('ptchga',imn,'.pc"')
  ptcb='%pointcharges "'+'%s%d%s'%('ptchgb',imn,'.pc"')
  fda.write(ptca)
  fda.write('\n')
  fdb.write(ptcb)
  fdb.write('\n')
  if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.gbw')):    
      wavefunctiona='%moinp '+'"'+'%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.gbw'+'"')
      fda.write(wavefunctiona)
      fda.write('\n')
  if os.path.exists('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.gbw')):    
      wavefunctionb='%moinp '+'"'+'%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.gbw'+'"')
      fdb.write(wavefunctionb)
      fdb.write('\n')
  fda.write('\n')
  fdb.write('\n')
  fda.write('* xyz ')
  fdb.write('* xyz ')
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

def mecp_orcamain(imn,cwd,usrdir,qmjob_prefix,nqm,nlink,cln,qmkey,cha_mul1,cha_mul2,extra_basis,orca_head):
  '''This part of program is to prepare the input file, and then to run ORCA and to extract the results from output. It uses two subprograms (morcainp, orcaout). It takes arguments: imn (integer, the image number); cwd, (string, current working directory), usrdir (string, user directory), qmjob_prefix (string, prefix for QM job input file), nqm (integer, number of QM atoms); nlink (integer, number of link atoms), cln (integer, QoMMMa cycle number), qmkey (string, any extra keys specified by user for ORCA job, default='None', i.e.'HF/STO-3G'), cha_mul1 and cha_mul2 are the charge and multiplicity of two states in MECP, given as list of two numbers, extra_basis is user specified basis set as a string (default='None'), and gau_head is used to assign memory and number of processors for ORCA jobs, given as a string (default='None'). This main routine is called from QoMMMa main program at each cycle if MECP is requested. Before preparing new input files this routine will rename all ORCA files created in previous cycle.''' 
  if os.path.exists(cwd+('%s%d'%('/image',imn))+('%s%s%d%s'%('/',qmjob_prefix,imn,'_A.in'))):
   try:			#rearrangeing previous cycle files
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_A.in')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.in'))) 
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_A.log')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.log')))
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_A.pcgrad')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.pcgrad')))
    shutil.copy(('%s%d%s'%('ptchga',imn,'.pc')),('%s%d%s'%('old_ptchga',imn,'.pc')))
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_A.gbw')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_A.gbw')))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_A.in'))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_A.log'))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_A.pcgrad'))
    os.remove('%s%d%s'%('ptchga',imn,'.pc'))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_A.gbw'))
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_B.in')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.in'))) 
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_B.log')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.log')))
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_B.pcgrad')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.pcgrad')))
    shutil.copy(('%s%d%s'%('ptchgb',imn,'.pc')),('%s%d%s'%('old_ptchgb',imn,'.pc')))
    shutil.copy(('%s%d%s'%(qmjob_prefix,imn,'_B.gbw')),('%s%s%d%s'%('old_',qmjob_prefix,imn,'_B.gbw')))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_B.in'))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_B.log'))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_B.pcgrad'))
    os.remove('%s%d%s'%('ptchgb',imn,'.pc'))
    os.remove('%s%d%s'%(qmjob_prefix,imn,'_B.gbw'))
   except:
    qomend('error while rearrangeing ORCA files created in previous run at cycle :'+str(cln)+' for image: '+str(imn),cwd,usrdir)
  try:			# preparing input file for both states 
    morcainp(imn,cwd,qmjob_prefix,qmkey,cha_mul1,cha_mul2,extra_basis,orca_head)
  except:
    qomend('Error while trying to prepare ORCA input files at cycle : '+str(cln)+' for image : '+str(imn),cwd,usrdir)
  try:			               		# ORCA job
    os.system('orca '+('%s%d%s'%(qmjob_prefix,imn,'_A.in'))+' > '+('%s%d%s'%(qmjob_prefix,imn,'_A.log')))
  except:
    qomend('Error while running ORCA job for state A at cycle : '+str(cln)+' for image : '+str(imn),cwd,usrdir)
  try:			               		# ORCA job
    os.system('orca '+('%s%d%s'%(qmjob_prefix,imn,'_B.in'))+' > '+('%s%d%s'%(qmjob_prefix,imn,'_B.log')))
  except:
    qomend('Error while running ORCA job for state B at cycle : '+str(cln)+' for image : '+str(imn),cwd,usrdir)
  try:				# extracting QM results
    fi='%s%d%s'%(qmjob_prefix,imn,'_A.log')
    pcg='%s%d%s'%(qmjob_prefix,imn,'_A.pcgrad')
    orcaout(fi,imn,nqm,nlink,pcg)	
    fi='%s%d%s'%(qmjob_prefix,imn,'_B.log')
    pcg='%s%d%s'%(qmjob_prefix,imn,'_B.pcgrad')
    orcaout(fi,imn,nqm,nlink,pcg)	
    shutil.copy(('%s%d'%('ab_initio',imn)),cwd)
    shutil.copy(('%s%d'%('mulliken',imn)),cwd)
    shutil.copy(('%s%d%s'%('qmlatgrad',imn,'.out')),cwd) 
    os.remove('%s%d%s'%('qmlatgrad',imn,'.out'))
    os.remove('%s%d'%('ab_initio',imn))
    os.remove('%s%d'%('mulliken',imn))
  except:
    qomend('Error while extracting results from ORCA output files at cycle : '+str(cln)+' for image : '+str(imn),cwd,usrdir)

