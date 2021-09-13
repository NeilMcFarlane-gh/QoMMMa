#!/usr/bin/python
import os
import sys
import shutil
import qomutil

"""This is QoMMMa program version 8 writern in python. This program is to couple and run QM, MM and QoMMMa ForTran codes. It consist of following modules"""

def fortinp():
  '''This program is to create two input files; 'fortinput' this file contains all information to setup qmmm calculation, used in read_initfile.f90, and 'converg.data' file containing convergence criteria for QoMMMa optimization used in read_converg.f90. The input for this part is already readed from default and user input files. Input should be preovided through specific keywords.'''
  fd=open('fortinput','w')
  fd.write('natom     nqm      nlink   nopt')
  fd.write('\n')
  fd.write(str(natom))
  fd.write(str(('%s%d'%('   ',nqm))))
  fd.write(str(('%s%d'%('   ',nlink))))
  fd.write(str(('%s%d'%('   ',nopt))))
  fd.write('\n')
  fd.write('disp')   # dispersion keyword
  fd.write('\n')
  fd.write(str((disp)))   # dispersion = 1 or 0 (on/off)
  fd.write('\n')          
  fd.write('ncon     kcnstype')
  fd.write('\n')
  fd.write(str(ncon))
  if kcnstype.lower()=='harmonic':
    ikcnstype=1
  elif kcnstype.lower()=='tanh':
    ikcnstype=2
  else:
     qomutil.qomlog('Error, unknown constrain type : '+kcnstype+' is requested',usrdir)
  fd.write(str(('%s%d'%('   ',ikcnstype))))
  fd.write('\n') 
  if nebtype.lower()=='none':
     inebtyp=0
  elif nebtype.lower()=='neb_bfgs':
     inebtyp=1
  elif nebtype.lower()=='cineb_bfgs':
     inebtyp=2
  elif nebtype.lower()=='neb_cg':
     inebtyp=3
  elif nebtype.lower()=='none_cg':
     inebtyp=4
  elif nebtype.lower()=='cineb_cg':
     inebtyp=5
  else:
     qomutil.qomlog('Error, unknown nebtype : '+nebtype+' is requested',usrdir)
  fd.write('nimg    nebtype    kspring')
  fd.write('\n')
  fd.write(str(nimg))
  fd.write(str(('%s%d'%('   ',inebtyp))))
  fd.write(str(('%s%d'%('   ',kspring))))
  fd.write('\n')
  fd.write('QM atoms list')
  fd.write('\n')
  if qm_lst!='None':			# qm atoms list
    fd.close()			
    qomutil.qmread(qm_lst,usrdir)
    fd=open('fortinput','a')
  else:
    fd.write((qm.lstrip()))		
    qomutil.qomlog('''Note, to provide QM atom list try to use 'qm_lst' option (see manual), instead of 'qm'. Since 'qm' option will be removed in the next version of QoMMMa.''',usrdir)
  if nlink > 0:
    fd.write('Link atom details')
    fd.write('\n')
    if link_lst!='None':
      fd.close()
      qomutil.link_write(link_lst) 
      fd=open('fortinput','a')
    else:
      fd.write((link.lstrip()))		# link atom(s) details (if any)
      qomutil.qomlog('''Note, to provide link atom details try to use 'link_lst' option (see manual), instead of 'link'. Since 'link' option will be removed in the next version of QoMMMa.''',usrdir)
  else:
    qomutil.qomlog('No Link atom details is taken, since nlink=0',usrdir) 
  if (int(nopt)-int(nqm)-int(nlink))>0:	
    fd.write('Hessian optimized MM atoms list')	
    fd.write('\n')			# Hessian optmized MM atoms (if any)
    fd.close()
    qomutil.hesopt_write(hesoptmm,usrdir)
    fd=open('fortinput','a')
  else:
    qomutil.qomlog('No MM atom is included in Hessian optimization',usrdir)
  if ncon > 0:				# for constrain (if any)
    fd.write('Constrain details')
    fd.write('\n')
    if constrain_lst!='None':
      fd.close()
      qomutil.cons_write(constrain_lst,usrdir)  
      fd=open('fortinput','a')
    else:
       fd.write((constrain.lstrip()))     
       qomutil.qomlog('''Note, to provide constrain details try to use 'constrain_lst' option (see manual), instead of 'constrain'. Since 'constrain' option will be removed in the next version of QoMMMa.''',usrdir)
  else:
    qomutil.qomlog('No constrain is taken into account, since ncon=0',usrdir)
  fd.write('Modified charge for MM atoms. First, number of such atoms, and then atom label and charge')
  fd.write('\n')			# for modified charges (default=0)
  if newcha_lst!='None':	# this kind of statements need to be removed/changed	
     fd.close()				# when 'newcha'	key is removed from QoMMMa
     qomutil.newcha_write(newcha_lst)
     fd=open('fortinput','a')
  else:
     fd.write(newcha.lstrip()) 
     if int(newcha.lstrip()[0])>0:          
       qomutil.qomlog('''Note, to provide new charge details try to use 'newcha_lst' option (see manual), instead of 'newcha'. Since 'newcha' option will be removed in the next version of QoMMMa.''',usrdir)
  fd.write('Inactive atoms. First, number of such atoms then list')
  fd.write('\n')
  kk=ninact.split()                    # for inactive atoms list (default=0)
  fd.write(str(kk[0]))		       
  if int(kk[0])!=0:
   fd.write('\n')
   kk=kk[1:]
   fd.close()
   qomutil.inprange(kk,'Inactive',usrdir)
  else: 
   qomutil.qomlog('No atom is set as inactive',usrdir)
   fd.close()
  fc=open('converg.data','w')
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
    '''This program is to create basic key files for Tinker job, creating image directories, checking number of input geometries, creating geom_expl*.xyz and to run Tinker analyze program to get initial charge. And  creating new geom*.key files from mm keys and inactive_list created by QoMMMa setup program.'''
    if mmcode=='Tinker':
       if mjob=='initial':
          default_prm=qommmadir+'/lib/charmm27.prm'
          fd=open('geom.key','w')      
          fd.write('parameters '+default_prm) 
          fd.write('\n')		
          fdt=open('tinker_header','w')
          fdt.write('parameters '+default_prm)
          fdt.write('\n')   	
          try:		
             fd.write(mmkey.lstrip())  
             fdt.write(mmkey.lstrip())
             qomutil.qomlog( 'Note, default force field parameter file is charmm27.prm, this could be changed or modified through mmkey in user input',usrdir)           
          except: 
             qomutil.qomlog('Note, no mmkey is given by user, default force field parameter file, charmm27.prm will be used here',usrdir) 
          fdt.close()
          fd.write('CHARGETERM ONLY')
          fd.write('\n')
          fd.write('EXTRATERM NONE')
          fd.write('\n')
          fd.close()
          for i in range(nimg): 
             l=i+1
             dirnam='%s%d' % ('image',l)
             os.mkdir(dirnam)
             fn= '%s%d%s' % ('geom',l,'.key')  
             shutil.copyfile('geom.key',fn) 
             fin='%s%d%s' % (inpgeom_prefix,l,'.xyz')
             fdn='%s%d%s' % ('geom',l,'.xyz')
             try:				# copying initial geometry file
                old=usrdir+'/'+fin         
                shutil.copy(old,fin) 
                shutil.copy(fin,fdn)
             except:
                qomutil.qomend('Error, input file : '+fin+' is not found in user directory',cwd,usrdir)
             try:
	        qomutil.geom_expl(fdn,l)
 	     except:
                qomutil.qomend('ERROR, while creating geom_expl*.xyz for image: '+str(l),cwd,usrdir)    
          try:			# Getting initial charge
              qomutil.tinkercharge(qommmadir) 
          except:
              qomutil.qomend('Error while getting initial charges using Tinker program analyze_grad',cwd,usrdir) 
       elif mjob=='initial2':
          try:
             fin1=open('tinker_header','r')
             fin2=open('inactive_list','r')
             fdg=open('geom.key','w')
             for line in fin1:
                   fdg.write(line)
             for line in fin2:
                   fdg.write(line)
             fin1.close()
             fin2.close()
             fdg.close()
             l=0
             for ir in range(nimg):
                 l=l+1
                 shutil.copy('geom.key',('%s%d%s'%('geom',l,'.key'))) 
             os.remove('inactive_list') 
          except:
             qomutil.qomend('Error: while preparing geom.key file after QoMMMa initial run. This step is used to include inactive atom list (including QM atoms) in tinker key file',cwd,usrdir)
       else:
          qomutil.qomend('Unknown MM job: '+job,cwd,usrdir)
    else:
          qomutil.qomend('Unknown MM code: '+mmcode,cwd,usrdir)

def qminitial():
   ''' This program creates the initial inputs for QM job, details were already readed from default and user input files.'''
   if qmcode.lower()=='jaguar':
     if job.lower()!='mecp': 
       fds=open('qmheader','w')
       fds.write(qmjag_header)   
       if qmkey.strip()!='None':
           fds.write(qmkey.lstrip())		
       fds.write('&')
       fds.write('\n') 
       fds.write('&zmat')
       fds.close()
     elif job.lower()=='mecp':
       from jagutil import mjag_initial
       mjag_initial(qmjag_header,qmkey,cha_mul1,cha_mul2)
     if atomic_section!='None':
          fd=open('atomic_section','w')
          fd.write(atomic_section.lstrip())
          fd.close()
     else: 
         from jagutil import atomic
         atomic()
     imn=0
     for i in range(nimg):
       imn=imn+1
       dst=cwd+('%s%d'%('/image',imn))
       if os.path.exists(usrdir+'%s%d'%('/orbs',imn)):
          shutil.copy(usrdir+'%s%d'%('/orbs',imn),dst)
          qomutil.qomlog('initial orbital guess for Jaguar job is copied from user directory for image :'+str(imn),usrdir)
       if job.lower()=='mecp':
         if os.path.exists(usrdir+'%s%d'%('/orbsA',imn)):
           shutil.copy(usrdir+'%s%d'%('/orbsA',imn),dst)
           qomutil.qomlog('initial orbital guess for Jaguar job is copied from user directory for first state of image :'+str(imn),usrdir)
         if os.path.exists(usrdir+'%s%d'%('/orbsB',imn)):
           shutil.copy(usrdir+'%s%d'%('/orbsB',imn),dst)
           qomutil.qomlog('initial orbital guess for Jaguar job is copied from user directory for second state of image :'+str(imn),usrdir)
   elif qmcode.lower()=='molpro':
       fds=open('qmheader','w')
       fds.write(qmmol_header.lstrip())
       fds.close()
       if job.lower()=='mecp':
         fda=open('qmmol_footer_A','w')
         fdb=open('qmmol_footer_B','w')
         if qmkey.strip()!='None':
           fda.write(qmkey.lstrip())
           fdb.write(qmkey.lstrip())
         else:
           qomutil.qomend('''For Molpro, you have to give atleast one input in 'qmkey' for instance 'hf' is must and may be basis set option''',cwd,usrdir) 
         fda.write(cha_mul1)
         fdb.write(cha_mul2)
         fda.write('\n')
         fdb.write('\n')
         fda.write(mol_footer.lstrip())
         fdb.write(mol_footer.lstrip())
         fda.close()
         fdb.close()
         imn=0 
         for i in range(nimg):
          imn=imn+1
          dst=cwd+('%s%d'%('/image',imn))
          if os.path.exists(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'_A.intg')):
             shutil.copy(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'_A.intg'),dst)
             qomutil.qomlog('initial orbital guess file is copied from user directory for first mecp state of image :'+str(imn),usrdir)
          if os.path.exists(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'_B.intg')):
             shutil.copy(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'_B.intg'),dst)
             qomutil.qomlog('initial orbital guess file is copied from user directory for second mecp state of image :'+str(imn),usrdir)
       else:
         fds=open('qmmol_footer','w')
         if qmkey.strip()!='None':
           fds.write(qmkey.lstrip())
         else:
           qomutil.qomend('''For Molpro, you have to give atleast one input in 'qmkey' for instance 'hf' is must''',cwd,usrdir) 
         fds.write(mol_footer.lstrip())
         fds.close()
         imn=0 
         for i in range(nimg):
          imn=imn+1
          if os.path.exists(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'.intg')):
             dst=cwd+('%s%d'%('/image',imn))
             shutil.copy(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'.intg'),dst)
             qomutil.qomlog('initial orbital guess file is copied from user directory for image :'+str(imn),usrdir)
   elif qmcode.lower()=='gaussian':
      if qmkey.strip()=='None':
         qomutil.qomlog( 'Note, no extra QM option is given by user through qmkey, so Gaussian job will run at HF/sto-3G level',usrdir)
      imn=0 
      for i in range(nimg):
          imn=imn+1
          dst=cwd+('%s%d'%('/image',imn))
          if os.path.exists(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'.chk')):
             shutil.copy(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'.chk'),dst)
             qomutil.qomlog('initial Gaussian check file was copied from user directory for image :'+str(imn),usrdir)
          if job.lower()=='mecp':
           if os.path.exists(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'_A.chk')):
             shutil.copy(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'_A.chk'),dst)
             qomutil.qomlog( 'initial Gaussian check file was taken from user directory for MECP state A of image :'+str(imn),usrdir) 
           if os.path.exists(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'_B.chk')):
             shutil.copy(usrdir+'%s%s%d%s'%('/',qmjob_prefix,imn,'_B.chk'),dst)
             qomutil.qomlog('initial Gaussian check file was taken from user directory for MECP state B of image :'+str(imn),usrdir) 
   else:
       qomutil.qomend('Unknown QM code: '+qmcode,cwd,usrdir) 

def qmmm():
      '''This is the main QoMMMa program part, first it will execute MM program to get MM gradient, then it will prepare the input for QM program, runs the QM program and extracts the results from the QM run. Then QoMMMa optimization will begin. The convegence will be checked. Once convergence is achived it will move files to user directory and stops the QoMMMa program [or it will start the frequency calculation if the job='Freq' (default job='opt') using qomfreq program. If convergence is not achived the MM optimization will follow, then a postmicro job will create/rearrange the files to use for next QoMMMa cycle. This will continue until the convergence is achived or maxcycle.''' 
      imn=0
      for irm in range(nimg):
        imn=imn+1
        if os.path.exists(cwd+('%s%d'%('/update_geom',imn))):
          try:    			# MM gradient run
              fch=os.popen(qommmadir+'/bin/analyze_grad','w')	
              fn=('%s%d%s'%('geom',imn,'.xyz'))
              fch.write(fn+'\nE\n')
              fch.close()
              shutil.copy('mm_grad',('%s%d'%('mm_grad',imn)))
              os.remove('mm_grad')
          except:
	      qomutil.qomend('Problem with MM gradient extraction at step :'+str(cln),cwd,usrdir)         
          os.chdir(('%s%d'%('image',imn)))
          qomutil.qomlog('About to run '+qmcode+ ' job at cycle : '+str(cln)+ '  for image : '+str(imn),usrdir)
          if qmcode.lower()=='jaguar':
           if job.lower()!='mecp':				# qmjob
             from jagutil import qmjagmain
             qmjagmain(imn,cwd,usrdir,qmjob_prefix,qmjag_job,nqm,nlink,cln)
             qmc_job=qmjag_job			# used in frequency call
           elif job.lower()=='mecp':
	     from jagutil import mecp_jagmain
             mecp_jagmain(imn,cwd,usrdir,qmjob_prefix,qmjag_job,nqm,nlink,cln)
          elif qmcode.lower()=='molpro':
           if job.lower()!='mecp':				# qmjob
             from molutil import qmmolmain
             qmmolmain(imn,cwd,usrdir,qmjob_prefix,qmmol_job,nqm,nlink,cln)
             qmc_job=qmmol_job 
           elif job.lower()=='mecp':				# qmjob
             from molutil import mecp_molmain
             mecp_molmain(imn,cwd,usrdir,qmjob_prefix,qmmol_job,nqm,nlink,cln)
          elif qmcode.lower()=='gaussian':
           if job.lower()!='mecp':				# qmjob
             from gauutil import qmgaumain
             qmgaumain(imn,cwd,usrdir,qmjob_prefix,qmgau_job,nqm,nlink,cln,qmkey,cha_mul,extra_basis,gau_head)
             qmc_job=qmgau_job 
           elif job.lower()=='mecp':
             from gauutil import mecp_gaumain
	     mecp_gaumain(imn,cwd,usrdir,qmjob_prefix,qmgau_job,nqm,nlink,cln,qmkey,cha_mul1,cha_mul2,extra_basis,gau_head)
          os.chdir(cwd)
          fi='%s%d%s'%('geom',imn,'.xyz')               
          try:				# from geom*.xyz creating geom_expl*.xyz
              qomutil.geom_expl(fi,imn) 	        
	  except:
              qomutil.qomend('ERROR, while creating geom_expl*.xyz for image: '+str(l),cwd,usrdir)    
          os.remove(('%s%d'%('update_geom',imn)))
        else: pass
      qomutil.qomlog('About to run QM/MM optimization cyle :'+str(cln),usrdir)
      try:		# QoMMMa optimization
        if job.lower()=='mecp':						
          os.system(qommmadir+'/bin/qommma_mecp.x')
        else:
          os.system(qommmadir+'/bin/qommma8.x')  		
          sys.exit
      except:
          qomutil.qomend('Error while running QoMMMa program at cycle :'+str(cln),cwd,usrdir)
      try: 
          qomutil.qomreport(nimg,usrdir,cwd,cln,qomout)
      except:
          qomutil.qomend('Error while looking for convergence after QoMMMa optimization at cycle :'+str(cln)+' check QM, MM results or checkfile',cwd,usrdir)
      if os.path.exists(cwd+'/convergence_ok'):
          im=0
          os.remove(cwd+'/convergence_ok')
          for ir in range(nimg):
            im=im+1
            shutil.copy(cwd+('%s%d%s'%('/geom',im,'.xyz')),usrdir+('%s%d%s'%('/geom',im,'.xyz')))
          if job.lower()=='freq':
               hlp=asctime()  
               qomutil.qomlog('About to start QoMMMa Frequency calculation at  '+str(hlp),usrdir) 
               from qomfreq import freqmain
               freqmain(qommmadir,qmkey,qmc_job,usrdir,cha_mul,extra_basis,gau_head,prj_freq,qmjob_prefix,qmcode) 
          qomutil.qomend('Congratulations,',cwd,usrdir)
      try:   			# rearranging files 
          qomutil.rearrange(nimg,cwd,qmjob_prefix,cln)	 
      except:
          qomutil.qomend('Error while rearranging the files after QoMMMa optimization at cycle :'+str(cln),cwd,usrdir)
      im=0				# micro-optimization
      for ir in range(nimg):				 
          try:
             im=im+1
             if os.path.exists(cwd+('%s%d'%('/update_geom',im))):
	       dst=cwd+('%s%d'%('/image',im))
               qomutil.qomlog('About to run Microoptimization of geometry :'+str(im),usrdir) 
               os.chdir(dst)
               try:					
		  fch=os.popen(qommmadir+'/bin/minimize > tmpmic.out','w')
                  fn=('%s%d%s'%('microgeom',im,'.xyz'))
                  fch.write(fn+'\n'+str(acctink)+'\n')
                  fch.close()
                  shutil.copy('tmpmic.out',('%s%d'%('micro_opt_log_',im)))
                  os.remove('tmpmic.out')  
	       except:
                     qomutil.qomend('Error while running Tinker minimization job at cyle :'+str(cln),cwd,usrdir)
	       os.remove('gradcorrection')
               fi='%s%d%s'%('microgeom',im,'.xyz_2')     
               try:	# from microgeom*.xyz_2 creating geom_expl*.xyz
                  qomutil.geom_expl(fi,im) 		  	 
	       except:
                  qomutil.qomend('ERROR, while creating geom_expl*.xyz for image: '+str(l),cwd,usrdir)    
	       shutil.copy(('%s%d%s'%('geom_expl',im,'.xyz')),cwd+('%s%d%s'%('/geom_expl',im,'.xyz')))
               os.remove(('%s%d%s'%('microgeom',im,'.xyz_2')))
               os.remove('%s%d%s'%('geom_expl',im,'.xyz'))
               os.chdir(cwd)
             else:
                 qomutil.qomlog('image Geometry  '+str(im)+'  is fixed in this iteration',usrdir)
          except:
               qomutil.qomend('Error while doing mico-iteration at QM/MM cycle :'+str(cln),cwd,usrdir)
      try:
           os.system(qommmadir+'/bin/postmicro8.x')   
      except:
           qomutil.qomend('Problem in QoMMMa post-processing job at cyle :'+str(cln),cwd,usrdir)
      try:
            qomutil.postrearrange(nimg,cwd)
      except:
           qomutil.qomend('Error, while rearranging the files after QoMMMa post-processing  job at cyle :'+str(cln),cwd,usrdir) 

#########################################################################################################

usrdir=os.getcwd()

if os.path.exists(usrdir+'/QoMMMa8.log'):
  shutil.copy(usrdir+'/QoMMMa8.log',usrdir+'/QoMMMa8.log_')
  os.remove(usrdir+'/QoMMMa8.log')
  qomutil.qomlog('Note, log file, QoMMMa8.log of your previous QoMMMa run was moved as QoMMMa8.log_; before submitting next QoMMMa job if you need this file rename it since it will be deleted',usrdir)

from time import asctime
hlp=asctime()
hlp1= os.uname()
qomutil.qomlog('QoMMMa job starts at   '+str(hlp)+'\n'+'  and running in node  '+str(hlp1),usrdir)

# seting base directory
try:
    qommmadir=os.environ['QOMMMA']
except:
    qomutil.qomlog('Environment variable QOMMMA is not set. Set QOMMMA to the base directory of QoMMMa',usrdir)
    sys.exit()

# Reading defaults

try:
    execfile(qommmadir+'/lib/default.in')
except:
    qomutil.qomlog('Could not be able to open default file: '+qommmadir+'/lib/default.in',usrdir)
    sys.exit()

# check cmd line for input file

if len(sys.argv)<2:
    qomutil.qomlog('Usage: '+sys.argv[0]+' <file.in>',usrdir)
    sys.exit()
inpf=sys.argv[1]

# Reading input file

try:
    execfile(inpf)
except:
    qomutil.qomlog('Could not be able to open input file: '+inpf+'  This may be due to error in user input file. Check for necessary inputs and its formats, see manual',usrdir)
    sys.exit()

# creating 'qmmm_*' directory at TMPDIR

try:
    dirs=os.environ['TMPDIR']+'/qmmm_'
    dirs=dirs+str(os.getpid())
    os.mkdir(dirs)
except:
    qomutil.qomlog('Either Variable TMPDIR is not set or Could not be able to create temporary directory '+dirs,usrdir)
    sys.exit()

os.chdir(dirs)
cwd=os.getcwd()
mmjob('initial')     	# initial MM job

try:		# to create an input file for QoMMMa ForTran programs
    fortinp()		
except:
    qomutil.qomend('Error while creating input file (fortinput or converg.data) for QoMMMa, check for necessary input and its format, see manual',cwd,usrdir)

# QoMMMa "setup" program. This creates building blocks for QM & MM jobs, and initializes output.
try:
  os.system(qommmadir+'/bin/setup8.x')
except:
  qomutil.qomend('QoMMMa initial setup program fails',cwd,usrdir) 

qminitial()		 		# preparing QM initial input  
qomutil.setini(cwd,nimg,usrdir) 	# moving results of QoMMMa setup program
mmjob('initial2')	 		# preparing new MM key files by using QoMMMa setup output

# Main QM/MM run: 'cln' is the qmmm cycle number

if maxcycle==0 and job.lower()=='freq':
  if len(prj_freq)==0:
    qomutil.qomlog('Frequency job is started with the geometry given in inpgeom_prefix*.xyz file, here qmcode used is Gaussian',usrdir) 
    from qomfreq import freqmain
    qmc_job=qmgau_job
    freqmain(qommmadir,qmkey,qmc_job,usrdir,cha_mul,extra_basis,gau_head,prj_freq,qmjob_prefix,qmcode) 
  else:
    qomutil.qomend('Error, Projected Frequency calculation is requested; set maxcycle=1 or request full optimization and frequency calculation',cwd,usrdir)
else:
 cln=0			
 for ic in range(maxcycle): 		
   cln=cln+1
   qmmm()
 if maxcycle==1 and job.lower()=='freq':
    qomutil.qomlog('Frequency job is started after first QoMMMa run, Frequency job will be performed with the geometry given in '+inpgeom_prefix+'*.xyz file, here qmcode used is Gaussian',usrdir) 
    from qomfreq import freqmain
    qmc_job=qmgau_job
    freqmain(qommmadir,qmkey,qmc_job,usrdir,cha_mul,extra_basis,gau_head,prj_freq,qmjob_prefix,qmcode) 
 elif cln == maxcycle: 
    qomutil.qomend('QoMMMa run reaches maxcycle :'+str(cln),cwd,usrdir)
