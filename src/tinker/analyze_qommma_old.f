c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  program analyze  --  energy partitioning and analysis  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "analyze" computes and displays the total potential; options
c     are provided to partition the energy by atom or by potential
c     function type; parameters used in computing interactions can
c     also be displayed by atom; output of large energy interactions
c     and of electrostatic and inertial properties is available
c
c
      program analyze
      implicit none
      include 'sizes.i'
      include 'atoms.i'
      include 'files.i'
      include 'inform.i'
      include 'iounit.i'
      integer i,j,ixyz
      integer frame
      integer freeunit
      integer trimtext
      integer list(20)
      logical doenergy,doatom
      logical dolarge,dodetail
      logical doprops,doparam
      logical exist
      logical active(maxatm)
      character*1 letter
      character*60 xyzfile
      character*80 record
      character*80 string
c
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     get the desired types of analysis to be performed
c
      call nextarg (string,exist)
      if (.not. exist) then
         write (iout,10)
   10    format (/,' The TINKER Analysis Facility can Provide :',
     &           //,' Total Potential Energy and its Components [E]',
     &           /,' Energy Breakdown over each of the Atoms [A]',
     &           /,' List of the Large Individual Interactions [L]',
     &           /,' Details for All Individual Interactions [D]',
     &           /,' Electrostatic, Inertial & Virial Properties [M]',
     &           /,' Force Field Parameters for Interactions [P]')
   20    continue
         write (iout,30)
   30    format (/,' Enter the Desired Analysis Types',
     &              ' [E,A,L,D,M,P] :  ',$)
         read (input,40,err=20)  string
   40    format (a80)
      end if
c
c     set option control flags based desired analysis types
c
      doenergy = .false.
      doatom = .false.
      dolarge = .false.
      dodetail = .false.
      doprops = .false.
      doparam = .false.
      call upcase (string)
      do i = 1, trimtext(string)
         letter = string(i:i)
         if (letter .eq. 'E')  doenergy = .true.
         if (letter .eq. 'A')  doatom = .true.
         if (letter .eq. 'L')  dolarge = .true.
         if (letter .eq. 'D')  dodetail = .true.
         if (letter .eq. 'M')  doprops = .true.
         if (letter .eq. 'P')  doparam = .true.
      end do
c
c     get the list of atoms for which output is desired
c
      if (doatom .or. doparam) then
         do i = 1, 20
            list(i) = 0
         end do
         if (exist) then
            do i = 1, 20
               call nextarg (string,exist)
               if (.not. exist)  goto 50
               read (string,*,err=50,end=50)  list(i)
            end do
   50       continue
         else
            write (iout,60)
   60       format (/,' List Atoms for which Output is Desired',
     &                 ' [ALL] :  '/,'    >  ',$)
            read (input,70)  record
   70       format (a80)
            read (record,*,err=80,end=80)  (list(i),i=1,20)
   80       continue
         end if
         do i = 1, n
            active(i) = .true.
         end do
         i = 1
         dowhile (list(i) .ne. 0)
            if (i .eq. 1) then
               do j = 1, n
                  active(j) = .false.
               end do
            end if
            if (list(i) .gt. 0) then
               active(list(i)) = .true.
               i = i + 1
            else
               do j = abs(list(i)), abs(list(i+1))
                  active(j) = .true.
               end do
               i = i + 2
            end if
         end do
      end if
c
c     setup to write out the large individual energy terms
c
      if (dolarge) then
         verbose = .true.
      else
         verbose = .false.
      end if
c
c     setup to write out all of the individual energy terms
c
      if (dodetail) then
         doenergy = .true.
         debug = .true.
      else
         debug = .false.
      end if
c
c     reopen the coordinates file and read the first structure
c
      frame = 0
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz')
      call version (xyzfile,'old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     perform analysis for each successive coordinate structure
c
      dowhile (.not. abort)
         frame = frame + 1
         if (frame .gt. 1) then
            write (iout,90)  frame
   90       format (/,' Analysis for Archive Structure :',8x,i8)
         end if
c
c     make the call to compute the potential energy
c
         if (doenergy .or. doatom .or. dolarge)  call enrgyze
c
c     energy partitioning by potential energy components
c
         if (doenergy) then
            if (digits .ge. 8) then
               call analyz8
            else if (digits .ge. 6) then
               call analyz6
            else
               call analyz4
            end if
         end if
c
c     get various electrostatic and inertial properties
c
         if (doprops)  call propyze
c
c     energy partitioning over the individual atoms
c
         if (doatom)  call atomyze (active)
c
c     list parameters used for molecular mechanics potentials
c
         if (doparam)  call paramyze (active)
c
c     attempt to read next structure from the coordinate file
c
         call readxyz (ixyz)
      end do
c
c     perform any final tasks before program exit
c
      close (unit=ixyz)
      if (dodetail)  debug = .false.
      call final
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine enrgyze  --  compute & report energy analysis  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "energyze" is an auxiliary routine for the analyze program
c     that performs the energy analysis and prints the total and
c     intermolecular energies
c
c
      subroutine enrgyze
      implicit none
      integer freeunit
      integer i,j,igrd
      include 'sizes.i'
      include 'atoms.i'
      include 'atmtyp.i'
      include 'inform.i'
      include 'inter.i'
      include 'iounit.i'
      include 'molcul.i'
      real*8 energy
      real*8 derivs(3,maxatm)
c
c     this is my hack - to print the gradient on each atom!
c
      call gradient (energy,derivs)
      igrd = freeunit ()
      open (unit=igrd,file="mm_grad")
      write (igrd,'(A,5X,F25.17)') "Energy:",energy
      do i=1,n
          write (igrd,'(I5,A4,3F18.10)') i,name(i),derivs(1,i),
     &                 derivs(2,i),derivs(3,i)
      end do
      close(igrd)

c
c     perform the energy analysis by atom and component
c
      call analysis (energy)
c
c     print out the total potential energy of the system
c
      if (digits .ge. 8) then
         if (abs(energy) .lt. 1.0d10) then
            write (iout,10)  energy
   10       format (/,' Total Potential Energy :',8x,f20.8,
     &                    ' Kcal/mole')
         else
            write (iout,20)  energy
   20       format (/,' Total Potential Energy :',8x,d20.8,
     &                    ' Kcal/mole')
         end if
      else if (digits .ge. 6) then
         if (abs(energy) .lt. 1.0d10) then
            write (iout,30)  energy
   30       format (/,' Total Potential Energy :',8x,f18.6,
     &                    ' Kcal/mole')
         else
            write (iout,40)  energy
   40       format (/,' Total Potential Energy :',8x,d18.6,
     &                    ' Kcal/mole')
         end if
      else
         if (abs(energy) .lt. 1.0d10) then
            write (iout,50)  energy
   50       format (/,' Total Potential Energy :',8x,f16.4,
     &                    ' Kcal/mole')
         else
            write (iout,60)  energy
   60       format (/,' Total Potential Energy :',8x,d16.4,
     &                    ' Kcal/mole')
         end if
      end if
c
c     intermolecular energy for systems with multiple molecules
c
      if (nmol.gt.1 .and. nmol.lt.n) then
         if (digits .ge. 8) then
            if (abs(einter) .lt. 1.0d10) then
               write (iout,70)  einter
   70          format (/,' Intermolecular Energy :',9x,f20.8,
     &                       ' Kcal/mole')
            else
               write (iout,80)  einter
   80          format (/,' Intermolecular Energy :',9x,d20.8,
     &                       ' Kcal/mole')
            end if
         else if (digits .ge. 6) then
            if (abs(einter) .lt. 1.0d10) then
               write (iout,90)  einter
   90          format (/,' Intermolecular Energy :',9x,f18.6,
     &                       ' Kcal/mole')
            else
               write (iout,100)  einter
  100          format (/,' Intermolecular Energy :',9x,d18.6,
     &                       ' Kcal/mole')
            end if
         else
            if (abs(einter) .lt. 1.0d10) then
               write (iout,110)  einter
  110          format (/,' Intermolecular Energy :',9x,f16.4,
     &                       ' Kcal/mole')
            else
               write (iout,120)  einter
  120          format (/,' Intermolecular Energy :',9x,d16.4,
     &                       ' Kcal/mole')
            end if
         end if
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine analyz4  --  low precision energy components  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "analyz4" prints the energy to 4 decimal places and number
c     of interactions for each component of the potential energy
c
c
      subroutine analyz4
      include 'action.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'iounit.i'
      include 'potent.i'
      character*8 label
c
c
c     write out each energy component to the desired precision
c
      write (iout,10)
   10 format (/,' Energy Component Breakdown :',11x,'Kcal/mole',
     &          6x,'Interactions'/)
      if (use_bond .and. neb.ne.0) then
         write (iout,20)  eb,neb
   20    format (' Bond Stretching',17x,f16.4,i15)
      end if
      if (use_angle .and. nea.ne.0) then
         write (iout,30)  ea,nea
   30    format (' Angle Bending',19x,f16.4,i15)
      end if
      if (use_strbnd .and. neba.ne.0) then
         write (iout,40)  eba,neba
   40    format (' Stretch-Bend',20x,f16.4,i15)
      end if
      if (use_urey .and. neub.ne.0) then
         write (iout,50)  eub,neub
   50    format (' Urey-Bradley',20x,f16.4,i15)
      end if
      if (use_angang .and. neaa.ne.0) then
         write (iout,60)  eaa,neaa
   60    format (' Angle-Angle',21x,f16.4,i15)
      end if
      if (use_opbend .and. neopb.ne.0) then
         write (iout,70)  eopb,neopb
   70    format (' Out-of-Plane Bend',15x,f16.4,i15)
      end if
      if (use_opdist .and. neopd.ne.0) then
         write (iout,80)  eopd,neopd
   80    format (' Out-of-Plane Distance',11x,f16.4,i15)
      end if
      if (use_improp .and. neid.ne.0) then
         write (iout,90)  eid,neid
   90    format (' Improper Dihedral',15x,f16.4,i15)
      end if
      if (use_imptor .and. neit.ne.0) then
         write (iout,100)  eit,neit
  100    format (' Improper Torsion',16x,f16.4,i15)
      end if
      if (use_tors .and. net.ne.0) then
         write (iout,110)  et,net
  110    format (' Torsional Angle',17x,f16.4,i15)
      end if
      if (use_strtor .and. nebt.ne.0) then
         write (iout,120)  ebt,nebt
  120    format (' Stretch-Torsion',17x,f16.4,i15)
      end if
      if (use_tortor .and. nett.ne.0) then
         write (iout,130)  ett,nett
  130    format (' Torsion-Torsion',17x,f16.4,i15)
      end if
      if (use_vdw .and. nev.ne.0) then
         if (abs(ev) .lt. 1.0d10) then
            write (iout,140)  ev,nev
  140       format (' Van der Waals',19x,f16.4,i15)
         else
            write (iout,150)  ev,nev
  150       format (' Van der Waals',19x,d16.4,i15)
         end if
      end if
      if (use_charge .and. nec.ne.0) then
         label = '        '
         if (use_ewald)  label = ' (Ewald)'
         if (abs(ec) .lt. 1.0d10) then
            write (iout,160)  label,ec,nec
  160       format (' Charge-Charge',a8,11x,f16.4,i15)
         else
            write (iout,170)  label,ec,nec
  170       format (' Charge-Charge',a8,11x,d16.4,i15)
         end if
      end if
      if (use_chgdpl .and. necd.ne.0) then
         if (abs(ecd) .lt. 1.0d10) then
            write (iout,180)  ecd,necd
  180       format (' Charge-Dipole',19x,f16.4,i15)
         else
            write (iout,190)  ecd,necd
  190       format (' Charge-Dipole',19x,d16.4,i15)
         end if
      end if
      if (use_dipole .and. ned.ne.0) then
         if (abs(ed) .lt. 1.0d10) then
            write (iout,200)  ed,ned
  200       format (' Dipole-Dipole',19x,f16.4,i15)
         else
            write (iout,210)  ed,ned
  210       format (' Dipole-Dipole',19x,d16.4,i15)
         end if
      end if
      if (use_mpole .and. nem.ne.0) then
         label = '        '
         if (use_ewald)  label = ' (Ewald)'
         if (abs(em) .lt. 1.0d10) then
            write (iout,220)  label,em,nem
  220       format (' Atomic Multipoles',a8,7x,f16.4,i15)
         else
            write (iout,230)  label,em,nem
  230       format (' Atomic Multipoles',a8,7x,d16.4,i15)
         end if
      end if
      if (use_polar .and. nep.ne.0) then
         label = '        '
         if (use_ewald)  label = ' (Ewald)'
         if (abs(ep) .lt. 1.0d10) then
            write (iout,240)  label,ep,nep
  240       format (' Polarization',a8,12x,f16.4,i15)
         else
            write (iout,250)  label,ep,nep
  250       format (' Polarization',a8,12x,d16.4,i15)
         end if
      end if
      if (use_rxnfld .and. ner.ne.0) then
         write (iout,260)  er,ner
  260    format (' Reaction Field',18x,f16.4,i15)
      end if
      if (use_solv .and. nes.ne.0) then
         write (iout,270)  es,nes
  270    format (' Continuum Solvation',13x,f16.4,i15)
      end if
      if (use_metal .and. nelf.ne.0) then
         write (iout,280)  elf,nelf
  280    format (' Metal Ligand Field',14x,f16.4,i15)
      end if
      if (use_geom .and. neg.ne.0) then
         write (iout,290)  eg,neg
  290    format (' Geometric Restraints',12x,f16.4,i15)
      end if
      if (use_extra .and. nex.ne.0) then
         write (iout,300)  ex,nex
  300    format (' Extra Energy Terms',14x,f16.4,i15)
      end if
      return
      end
c
c
c     ##################################################################
c     ##                                                              ##
c     ##  subroutine analyz6  --  medium precision energy components  ##
c     ##                                                              ##
c     ##################################################################
c
c
c     "analyz6" prints the energy to 6 decimal places and number
c     of interactions for each component of the potential energy
c
c
      subroutine analyz6
      include 'action.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'iounit.i'
      include 'potent.i'
      character*8 label
c
c
c     write out each energy component to the desired precision
c
      write (iout,10)
   10 format (/,' Energy Component Breakdown :',13x,'Kcal/mole',
     &          6x,'Interactions'/)
      if (use_bond .and. neb.ne.0) then
         write (iout,20)  eb,neb
   20    format (' Bond Stretching',17x,f18.6,i15)
      end if
      if (use_angle .and. nea.ne.0) then
         write (iout,30)  ea,nea
   30    format (' Angle Bending',19x,f18.6,i15)
      end if
      if (use_strbnd .and. neba.ne.0) then
         write (iout,40)  eba,neba
   40    format (' Stretch-Bend',20x,f18.6,i15)
      end if
      if (use_urey .and. neub.ne.0) then
         write (iout,50)  eub,neub
   50    format (' Urey-Bradley',20x,f18.6,i15)
      end if
      if (use_angang .and. neaa.ne.0) then
         write (iout,60)  eaa,neaa
   60    format (' Angle-Angle',21x,f18.6,i15)
      end if
      if (use_opbend .and. neopb.ne.0) then
         write (iout,70)  eopb,neopb
   70    format (' Out-of-Plane Bend',15x,f18.6,i15)
      end if
      if (use_opdist .and. neopd.ne.0) then
         write (iout,80)  eopd,neopd
   80    format (' Out-of-Plane Distance',11x,f18.6,i15)
      end if
      if (use_improp .and. neid.ne.0) then
         write (iout,90)  eid,neid
   90    format (' Improper Dihedral',15x,f18.6,i15)
      end if
      if (use_imptor .and. neit.ne.0) then
         write (iout,100)  eit,neit
  100    format (' Improper Torsion',16x,f18.6,i15)
      end if
      if (use_tors .and. net.ne.0) then
         write (iout,110)  et,net
  110    format (' Torsional Angle',17x,f18.6,i15)
      end if
      if (use_strtor .and. nebt.ne.0) then
         write (iout,120)  ebt,nebt
  120    format (' Stretch-Torsion',17x,f18.6,i15)
      end if
      if (use_tortor .and. nett.ne.0) then
         write (iout,130)  ett,nett
  130    format (' Torsion-Torsion',17x,f18.6,i15)
      end if
      if (use_vdw .and. nev.ne.0) then
         if (abs(ev) .lt. 1.0d10) then
            write (iout,140)  ev,nev
  140       format (' Van der Waals',19x,f18.6,i15)
         else
            write (iout,150)  ev,nev
  150       format (' Van der Waals',19x,d18.6,i15)
         end if
      end if
      if (use_charge .and. nec.ne.0) then
         label = '        '
         if (use_ewald)  label = ' (Ewald)'
         if (abs(ec) .lt. 1.0d10) then
            write (iout,160)  label,ec,nec
  160       format (' Charge-Charge',a8,11x,f18.6,i15)
         else
            write (iout,170)  label,ec,nec
  170       format (' Charge-Charge',a8,11x,d18.6,i15)
         end if
      end if
      if (use_chgdpl .and. necd.ne.0) then
         if (abs(ecd) .lt. 1.0d10) then
            write (iout,180)  ecd,necd
  180       format (' Charge-Dipole',19x,f18.6,i15)
         else
            write (iout,190)  ecd,necd
  190       format (' Charge-Dipole',19x,d18.6,i15)
         end if
      end if
      if (use_dipole .and. ned.ne.0) then
         if (abs(ed) .lt. 1.0d10) then
            write (iout,200)  ed,ned
  200       format (' Dipole-Dipole',19x,f18.6,i15)
         else
            write (iout,210)  ed,ned
  210       format (' Dipole-Dipole',19x,d18.6,i15)
         end if
      end if
      if (use_mpole .and. nem.ne.0) then
         label = '        '
         if (use_ewald)  label = ' (Ewald)'
         if (abs(em) .lt. 1.0d10) then
            write (iout,220)  label,em,nem
  220       format (' Atomic Multipoles',a8,7x,f18.6,i15)
         else
            write (iout,230)  label,em,nem
  230       format (' Atomic Multipoles',a8,7x,d18.6,i15)
         end if
      end if
      if (use_polar .and. nep.ne.0) then
         label = '        '
         if (use_ewald)  label = ' (Ewald)'
         if (abs(ep) .lt. 1.0d10) then
            write (iout,240)  label,ep,nep
  240       format (' Polarization',a8,12x,f18.6,i15)
         else
            write (iout,250)  label,ep,nep
  250       format (' Polarization',a8,12x,d18.6,i15)
         end if
      end if
      if (use_rxnfld .and. ner.ne.0) then
         write (iout,260)  er,ner
  260    format (' Reaction Field',18x,f18.6,i15)
      end if
      if (use_solv .and. nes.ne.0) then
         write (iout,270)  es,nes
  270    format (' Continuum Solvation',13x,f18.6,i15)
      end if
      if (use_metal .and. nelf.ne.0) then
         write (iout,280)  elf,nelf
  280    format (' Metal Ligand Field',14x,f18.6,i15)
      end if
      if (use_geom .and. neg.ne.0) then
         write (iout,290)  eg,neg
  290    format (' Geometric Restraints',12x,f18.6,i15)
      end if
      if (use_extra .and. nex.ne.0) then
         write (iout,300)  ex,nex
  300    format (' Extra Energy Terms',14x,f18.6,i15)
      end if
      return
      end
c
c
c     ################################################################
c     ##                                                            ##
c     ##  subroutine analyz8  --  high precision energy components  ##
c     ##                                                            ##
c     ################################################################
c
c
c     "analyz8" prints the energy to 8 decimal places and number
c     of interactions for each component of the potential energy
c
c
      subroutine analyz8
      include 'action.i'
      include 'cutoff.i'
      include 'energi.i'
      include 'iounit.i'
      include 'potent.i'
      character*8 label
c
c
c     write out each energy component to the desired precision
c
      write (iout,10)
   10 format (/,' Energy Component Breakdown :',15x,'Kcal/mole',
     &          6x,'Interactions'/)
      if (use_bond .and. neb.ne.0) then
         write (iout,20)  eb,neb
   20    format (' Bond Stretching',17x,f20.8,i15)
      end if
      if (use_angle .and. nea.ne.0) then
         write (iout,30)  ea,nea
   30    format (' Angle Bending',19x,f20.8,i15)
      end if
      if (use_strbnd .and. neba.ne.0) then
         write (iout,40)  eba,neba
   40    format (' Stretch-Bend',20x,f20.8,i15)
      end if
      if (use_urey .and. neub.ne.0) then
         write (iout,50)  eub,neub
   50    format (' Urey-Bradley',20x,f20.8,i15)
      end if
      if (use_angang .and. neaa.ne.0) then
         write (iout,60)  eaa,neaa
   60    format (' Angle-Angle',21x,f20.8,i15)
      end if
      if (use_opbend .and. neopb.ne.0) then
         write (iout,70)  eopb,neopb
   70    format (' Out-of-Plane Bend',15x,f20.8,i15)
      end if
      if (use_opdist .and. neopd.ne.0) then
         write (iout,80)  eopd,neopd
   80    format (' Out-of-Plane Distance',11x,f20.8,i15)
      end if
      if (use_improp .and. neid.ne.0) then
         write (iout,90)  eid,neid
   90    format (' Improper Dihedral',15x,f20.8,i15)
      end if
      if (use_imptor .and. neit.ne.0) then
         write (iout,100)  eit,neit
  100    format (' Improper Torsion',16x,f20.8,i15)
      end if
      if (use_tors .and. net.ne.0) then
         write (iout,110)  et,net
  110    format (' Torsional Angle',17x,f20.8,i15)
      end if
      if (use_strtor .and. nebt.ne.0) then
         write (iout,120)  ebt,nebt
  120    format (' Stretch-Torsion',17x,f20.8,i15)
      end if
      if (use_tortor .and. nett.ne.0) then
         write (iout,130)  ett,nett
  130    format (' Torsion-Torsion',17x,f20.8,i15)
      end if
      if (use_vdw .and. nev.ne.0) then
         if (abs(ev) .lt. 1.0d10) then
            write (iout,140)  ev,nev
  140       format (' Van der Waals',19x,f20.8,i15)
         else
            write (iout,150)  ev,nev
  150       format (' Van der Waals',19x,d20.8,i15)
         end if
      end if
      if (use_charge .and. nec.ne.0) then
         label = '        '
         if (use_ewald)  label = ' (Ewald)'
         if (abs(ec) .lt. 1.0d10) then
            write (iout,160)  label,ec,nec
  160       format (' Charge-Charge',a8,11x,f20.8,i15)
         else
            write (iout,170)  label,ec,nec
  170       format (' Charge-Charge',a8,11x,d20.8,i15)
         end if
      end if
      if (use_chgdpl .and. necd.ne.0) then
         if (abs(ecd) .lt. 1.0d10) then
            write (iout,180)  ecd,necd
  180       format (' Charge-Dipole',19x,f20.8,i15)
         else
            write (iout,190)  ecd,necd
  190       format (' Charge-Dipole',19x,d20.8,i15)
         end if
      end if
      if (use_dipole .and. ned.ne.0) then
         if (abs(ed) .lt. 1.0d10) then
            write (iout,200)  ed,ned
  200       format (' Dipole-Dipole',19x,f20.8,i15)
         else
            write (iout,210)  ed,ned
  210       format (' Dipole-Dipole',19x,d20.8,i15)
         end if
      end if
      if (use_mpole .and. nem.ne.0) then
         label = '        '
         if (use_ewald)  label = ' (Ewald)'
         if (abs(em) .lt. 1.0d10) then
            write (iout,220)  label,em,nem
  220       format (' Atomic Multipoles',a8,7x,f20.8,i15)
         else
            write (iout,230)  label,em,nem
  230       format (' Atomic Multipoles',a8,7x,d20.8,i15)
         end if
      end if
      if (use_polar .and. nep.ne.0) then
         label = '        '
         if (use_ewald)  label = ' (Ewald)'
         if (abs(ep) .lt. 1.0d10) then
            write (iout,240)  label,ep,nep
  240       format (' Polarization',a8,12x,f20.8,i15)
         else
            write (iout,250)  label,ep,nep
  250       format (' Polarization',a8,12x,d20.8,i15)
         end if
      end if
      if (use_rxnfld .and. ner.ne.0) then
         write (iout,260)  er,ner
  260    format (' Reaction Field',18x,f20.8,i15)
      end if
      if (use_solv .and. nes.ne.0) then
         write (iout,270)  es,nes
  270    format (' Continuum Solvation',13x,f20.8,i15)
      end if
      if (use_metal .and. nelf.ne.0) then
         write (iout,280)  elf,nelf
  280    format (' Metal Ligand Field',14x,f20.8,i15)
      end if
      if (use_geom .and. neg.ne.0) then
         write (iout,290)  eg,neg
  290    format (' Geometric Restraints',12x,f20.8,i15)
      end if
      if (use_extra .and. nex.ne.0) then
         write (iout,300)  ex,nex
  300    format (' Extra Energy Terms',14x,f20.8,i15)
      end if
      return
      end
c
c
c     #################################################################
c     ##                                                             ##
c     ##  subroutine propyze  --  electrostatic & inertial analysis  ##
c     ##                                                             ##
c     #################################################################
c
c
c     "propyze" finds and prints the total charge, dipole moment
c     components, radius of gyration and moments of inertia
c
c
      subroutine propyze
      include 'sizes.i'
      include 'chgpot.i'
      include 'iounit.i'
      include 'moment.i'
      include 'virial.i'
      real*8 rg,energy
      real*8 derivs(3,maxatm)
c
c
c     get the total charge, dipole and quadrupole moments
c
      call moments
      write (iout,10)  netchg
   10 format (/,' Total Electric Charge :',13x,f12.5,' Electrons')
      write (iout,20)  netdpl,xdpl,ydpl,zdpl
   20 format (/,' Dipole Moment Magnitude :',11x,f12.3,' Debyes',
     &        //,' Dipole X,Y,Z-Components :',11x,3f12.3)
      write (iout,30)  xxqdp,xyqdp,xzqdp,yxqdp,yyqdp,
     &                 yzqdp,zxqdp,zyqdp,zzqdp
   30 format (/,' Quadrupole Moment Tensor :',10x,3f12.3,
     &        /,6x,'(Buckinghams)',18x,3f12.3,
     &        /,37x,3f12.3)
      write (iout,40)  netqdp(1),netqdp(2),netqdp(3)
   40 format (/,' Principal Axes Quadrupole :',9x,3f12.3)
      if (dielec .ne. 1.0d0) then
         write (iout,50)  dielec
   50    format (/,' Dielectric Constant :',15x,f12.3)
         write (iout,60)  netchg/sqrt(dielec)
   60    format (' Effective Total Charge :',12x,f12.5,' Electrons')
         write (iout,70)  netdpl/sqrt(dielec)
   70    format (' Effective Dipole Moment :',11x,f12.3,' Debyes')
      end if
c
c     get the radius of gyration and moments of inertia
c
      call gyrate (rg)
      write (iout,80)  rg
   80 format (/,' Radius of Gyration :',16x,f12.3,' Angstroms')
      call inertia (1)
c
c     get the internal virial tensor via gradient calculation
c
      call gradient (energy,derivs)
      write (iout,90)  (vir(1,i),vir(2,i),vir(3,i),i=1,3)
   90 format (/,' Internal Virial Tensor :',12x,3f12.3,
     &        /,37x,3f12.3,/,37x,3f12.3)
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine atomyze  --  individual atom energy analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "atomyze" prints the potential energy components broken
c     down by atom and to a choice of precision
c
c
      subroutine atomyze (active)
      include 'sizes.i'
      include 'analyz.i'
      include 'atoms.i'
      include 'inform.i'
      include 'iounit.i'
      integer i
      logical active(maxatm)
c
c
c     energy partitioning over the individual atoms
c
      write (iout,10)
   10 format (/,' Potential Energy Breakdown over Atoms :')
      if (digits .ge. 8) then
         write (iout,20)
   20    format (/,'  Atom',9x,'EB',14x,'EA',14x,'EBA',13x,'EUB',
     &           /,15x,'EAA',13x,'EOPB',12x,'EOPD',12x,'EID',
     &           /,15x,'EIT',13x,'ET',14x,'EBT',13x,'ETT',
     &           /,15x,'EV',14x,'EC',14x,'ECD',13x,'ED',
     &           /,15x,'EM',14x,'EP',14x,'ER',14x,'ES',
     &           /,15x,'ELF',13x,'EG',14x,'EX')
         do i = 1, n
            if (active(i)) then
               write (iout,30)  i,aeb(i),aea(i),aeba(i),aeub(i),
     &                          aeaa(i),aeopb(i),aeopd(i),aeid(i),
     &                          aeit(i),aet(i),aebt(i),aett(i),
     &                          aev(i),aec(i),aecd(i),aed(i),aem(i),
     &                          aep(i),aer(i),aes(i),aelf(i),aeg(i),
     &                          aex(i)
   30          format (/,i6,4f16.8,/,6x,4f16.8,/,6x,4f16.8,
     &                 /,6x,4f16.8,/,6x,4f16.8,/,6x,3f16.8)
            end if
         end do
      else if (digits .ge. 6) then
         write (iout,40)
   40    format (/,'  Atom',8x,'EB',12x,'EA',12x,'EBA',11x,'EUB',
     &              11x,'EAA',
     &           /,14x,'EOPB',10x,'EOPD',10x,'EID',11x,'EIT',
     &              11x,'ET',
     &           /,14x,'EBT',11x,'ETT',11x,'EV',12x,'EC',12x,'ECD',
     &           /,14x,'ED',12x,'EM',12x,'EP',12x,'ER',12x,'ES',
     &           /,14x,'ELF',11x,'EG',12x,'EX')
         do i = 1, n
            if (active(i)) then
               write (iout,50)  i,aeb(i),aea(i),aeba(i),aeub(i),
     &                          aeaa(i),aeopb(i),aeopd(i),aeid(i),
     &                          aeit(i),aet(i),aebt(i),aett(i),
     &                          aev(i),aec(i),aecd(i),aed(i),aem(i),
     &                          aep(i),aer(i),aes(i),aelf(i),aeg(i),
     &                          aex(i)
   50          format (/,i6,5f14.6,/,6x,5f14.6,/,6x,5f14.6,
     &                 /,6x,5f14.6,/,6x,3f14.6)
            end if
         end do
      else
         write (iout,60)
   60    format (/,'  Atom',8x,'EB',10x,'EA',10x,'EBA',9x,'EUB',
     &              9x,'EAA',9x,'EOPB',
     &           /,14x,'EOPD',8x,'EID',9x,'EIT',9x,'ET',10x,'EBT',
     &              9x,'ETT',
     &           /,14x,'EV',10x,'EC',10x,'ECD',9x,'ED',10x,'EM',
     &              10x,'EP',
     &           /,14x,'ER',10x,'ES',10x,'ELF',9x,'EG',10x,'EX')
         do i = 1, n
            if (active(i)) then
               write (iout,70)  i,aeb(i),aea(i),aeba(i),aeub(i),
     &                          aeaa(i),aeopb(i),aeopd(i),aeid(i),
     &                          aeit(i),aet(i),aebt(i),aett(i),
     &                          aev(i),aec(i),aecd(i),aed(i),aem(i),
     &                          aep(i),aer(i),aes(i),aelf(i),aeg(i),
     &                          aex(i)
   70          format (/,i6,6f12.4,/,6x,6f12.4,/,6x,6f12.4,
     &                 /,6x,5f12.4)
            end if
         end do
      end if
      return
      end
c
c
c     ###############################################################
c     ##                                                           ##
c     ##  subroutine paramyze  --  force field parameter analysis  ##
c     ##                                                           ##
c     ###############################################################
c
c
c     "paramyze" prints the force field parameters used in the
c     computation of each of the potential energy terms
c
c
      subroutine paramyze (active)
      implicit none
      include 'sizes.i'
      include 'angang.i'
      include 'angle.i'
      include 'atmtyp.i'
      include 'atoms.i'
      include 'bond.i'
      include 'charge.i'
      include 'dipole.i'
      include 'improp.i'
      include 'imptor.i'
      include 'iounit.i'
      include 'kvdws.i'
      include 'mpole.i'
      include 'opbend.i'
      include 'opdist.i'
      include 'polar.i'
      include 'polgrp.i'
      include 'potent.i'
      include 'solute.i'
      include 'strbnd.i'
      include 'strtor.i'
      include 'tors.i'
      include 'units.i'
      include 'urey.i'
      include 'vdw.i'
      integer i,j,k
      integer ia,ib,ic,id
      integer izaxe,ixaxe,iyaxe
      integer fold(6)
      real*8 bla,blc
      real*8 ampli(6)
      real*8 phase(6)
      logical header
      logical active(maxatm)
c
c
c     list parameters used for molecular mechanics atom types
c
      header = .true.
      do i = 1, n
         if (active(i)) then
            if (header) then
               header = .false.
               write (iout,10)
   10          format (/,' Atom Type Definition Parameters :',
     &                 //,3x,'Atom',4x,'Symbol',3x,'Type',
     &                    2x,'Class',2x,'Atomic',4x,'Mass',
     &                    3x,'Valence',2x,'Description',/)
            end if
            write (iout,20)  i,name(i),type(i),class(i),atomic(i),
     &                       mass(i),valence(i),story(i)
   20       format (i6,7x,a3,3i7,f11.3,i6,5x,a20)
         end if
      end do
c
c     list parameters used for van der Waals interactions
c
      if (use_vdw) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,30)
   30             format (/,' Van der Waals Parameters :',
     &                    //,10x,'Atom Number',5x,'Radius',
     &                       3x,'Epsilon',3x,'Rad 1-4',
     &                       3x,'Eps 1-4',3x,'Reduction',/)
               end if
               j = class(i)
               if (rad(j).eq.rad4(j) .and. eps(j).eq.eps4(j)) then
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,40)  k,i,rad(j),eps(j)
   40                format (i6,3x,i6,7x,2f10.4)
                  else
                     write (iout,50)  k,i,rad(j),eps(j),reduct(j)
   50                format (i6,3x,i6,7x,2f10.4,20x,f10.4)
                  end if
               else
                  if (reduct(j) .eq. 0.0d0) then
                     write (iout,60)  k,i,rad(j),eps(j),rad4(j),
     &                                eps4(j)
   60                format (i6,3x,i6,7x,4f10.4)
                  else
                     write (iout,70)  k,i,rad(j),eps(j),rad4(j),
     &                                eps4(j),reduct(j)
   70                format (i6,3x,i6,7x,5f10.4)
                  end if
               end if
            end if
         end do
      end if
c
c     list parameters used for bond stretching interactions
c
      if (use_bond) then
         header = .true.
         do i = 1, nbond
            ia = ibnd(1,i)
            ib = ibnd(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,80)
   80             format (/,' Bond Stretching Parameters :',
     &                    //,10x,'Atom Numbers',25x,'KS',
     &                       5x,'Length',/)
               end if
               write (iout,90)  i,ia,ib,bk(i),bl(i)
   90          format (i6,3x,2i6,13x,f16.4,f10.4)
            end if
         end do
      end if
c
c     list parameters used for angle bending interactions
c
      if (use_angle) then
         header = .true.
         do i = 1, nangle
            ia = iang(1,i)
            ib = iang(2,i)
            ic = iang(3,i)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,100)
  100             format (/,' Angle Bending Parameters :',
     &                    //,13x,'Atom Numbers',22x,'KB',
     &                       6x,'Angle',/)
               end if
               if (angtyp(i) .eq. 'HARMONIC') then
                  write (iout,110)  i,ia,ib,ic,ak(i),anat(i)
  110             format (i6,3x,3i6,7x,f16.4,f10.4)
               else if (angtyp(i) .eq. 'IN-PLANE') then
                  write (iout,120)  i,ia,ib,ic,ak(i),anat(i)
  120             format (i6,3x,3i6,7x,f16.4,f10.4,2x,'In-Plane')
               end if
            end if
         end do
      end if
c
c     list parameters used for stretch-bend interactions
c
      if (use_strbnd) then
         header = .true.
         do i = 1, nstrbnd
            k = isb(1,i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            if (active(ia) .or. active(ib) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,130)
  130             format (/,' Stretch-Bend Parameters :',
     &                    //,13x,'Atom Numbers',11x,'KSB',
     &                       6x,'Angle',3x,'Length1',3x,'Length2',/)
               end if
               bla = 0.0d0
               blc = 0.0d0
               if (isb(2,i) .ne. 0)  bla = bl(isb(2,i))
               if (isb(3,i) .ne. 0)  blc = bl(isb(3,i))
               write (iout,140)  i,ia,ib,ic,ksb(i),anat(k),bla,blc
  140          format (i6,3x,3i6,f13.4,3f10.4)
            end if
         end do
      end if
c
c     list parameters used for Urey-Bradley interactions
c
      if (use_urey) then
         header = .true.
         do i = 1, nurey
            ia = iury(1,i)
            ib = iury(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,150)
  150             format (/,' Urey-Bradley Parameters :',
     &                    //,10x,'Atom Numbers',24x,'KUB',
     &                       4x,'Distance',/)
               end if
               write (iout,160)  i,ia,ib,uk(i),ul(i)
  160          format (i6,3x,2i6,13x,f16.4,f10.4)
            end if
         end do
      end if
c
c     list parameters used for out-of-plane bending interactions
c
      if (use_opbend) then
         header = .true.
         do i = 1, nopbend
            k = iopb(i)
            ia = iang(1,k)
            ib = iang(2,k)
            ic = iang(3,k)
            id = iang(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,170)
  170             format (/,' Out-of-Plane Bending Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KOPB',/)
               end if
               write (iout,180)  i,id,ib,ia,ic,kopb(i)
  180          format (i6,3x,4i6,9x,f10.4)
            end if
         end do
      end if
c
c     list parameters used for out-of-plane distance interactions
c
      if (use_opdist) then
         header = .true.
         do i = 1, nopdist
            ia = iopd(1,i)
            ib = iopd(2,i)
            ic = iopd(3,i)
            id = iopd(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,190)
  190             format (/,' Out-of-Plane Distance Parameters :',
     &                    //,17x,'Atom Numbers',18x,'KOPD',/)
               end if
               write (iout,200)  i,ia,ib,ic,id,kopd(i)
  200          format (i6,3x,4i6,8x,f10.4)
            end if
         end do
      end if
c
c     list parameters used for improper dihedral interactions
c
      if (use_improp) then
         header = .true.
         do i = 1, niprop
            ia = iiprop(1,i)
            ib = iiprop(2,i)
            ic = iiprop(3,i)
            id = iiprop(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,210)
  210             format (/,' Improper Dihedral Parameters :',
     &                    //,17x,'Atom Numbers',19x,'KID',
     &                       4x,'Dihedral',/)
               end if
               write (iout,220)  i,ia,ib,ic,id,kprop(i),vprop(i)
  220          format (i6,3x,4i6,9x,2f10.4)
            end if
         end do
      end if
c
c     list parameters used for improper torsion interactions
c
      if (use_imptor) then
         header = .true.
         do i = 1, nitors
            ia = iitors(1,i)
            ib = iitors(2,i)
            ic = iitors(3,i)
            id = iitors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,230)
  230             format (/,' Improper Torsion Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (itors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = itors1(1,i)
                  phase(j) = itors1(2,i)
               end if
               if (itors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = itors2(1,i)
                  phase(j) = itors2(2,i)
               end if
               if (itors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = itors3(1,i)
                  phase(j) = itors3(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,240)  i,ia,ib,ic,id
  240             format (i6,3x,4i6)
               else if (j .eq. 1) then
                  write (iout,250)  i,ia,ib,ic,id,
     &                              ampli(1),phase(1),fold(1)
  250             format (i6,3x,4i6,10x,f10.3,f8.1,i4)
               else if (j .eq. 2) then
                  write (iout,260)  i,ia,ib,ic,id,(ampli(k),
     &                              phase(k),fold(k),k=1,j)
  260             format (i6,3x,4i6,2x,2(f10.3,f6.1,i4))
               else
                  write (iout,270)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  270             format (i6,3x,4i6,4x,3(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     list parameters used for torsional interactions
c
      if (use_tors) then
         header = .true.
         do i = 1, ntors
            ia = itors(1,i)
            ib = itors(2,i)
            ic = itors(3,i)
            id = itors(4,i)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,280)
  280             format (/,' Torsional Angle Parameters :',
     &                    //,17x,'Atom Numbers',11x,
     &                       'Amplitude, Phase and Periodicity',/)
               end if
               j = 0
               if (tors1(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = tors1(1,i)
                  phase(j) = tors1(2,i)
               end if
               if (tors2(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = tors2(1,i)
                  phase(j) = tors2(2,i)
               end if
               if (tors3(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = tors3(1,i)
                  phase(j) = tors3(2,i)
               end if
               if (tors4(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 4
                  ampli(j) = tors4(1,i)
                  phase(j) = tors4(2,i)
               end if
               if (tors5(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 5
                  ampli(j) = tors5(1,i)
                  phase(j) = tors5(2,i)
               end if
               if (tors6(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 6
                  ampli(j) = tors6(1,i)
                  phase(j) = tors6(2,i)
               end if
               if (j .eq. 0) then
                  write (iout,290)  i,ia,ib,ic,id
  290             format (i6,3x,4i6)
               else
                  write (iout,300)  i,ia,ib,ic,id,(ampli(k),
     &                              nint(phase(k)),fold(k),k=1,j)
  300             format (i6,3x,4i6,4x,6(f8.3,i4,'/',i1))
               end if
            end if
         end do
      end if
c
c     list parameters used for stretch-torsion interactions
c
      if (use_strtor) then
         header = .true.
         do i = 1, nstrtor
            k = ist(1,i)
            ia = itors(1,k)
            ib = itors(2,k)
            ic = itors(3,k)
            id = itors(4,k)
            if (active(ia) .or. active(ib) .or.
     &             active(ic) .or. active(id)) then
               if (header) then
                  header = .false.
                  write (iout,310)
  310             format (/,' Stretch-Torsion Parameters :',
     &                    //,17x,'Atom Numbers',10x,'Length',
     &                       5x,'Torsion Terms',/)
               end if
               j = 0
               if (kst(1,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 1
                  ampli(j) = kst(1,i)
                  phase(j) = tors1(2,k)
               end if
               if (kst(2,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 2
                  ampli(j) = kst(2,i)
                  phase(j) = tors2(2,k)
               end if
               if (kst(3,i) .ne. 0.0d0) then
                  j = j + 1
                  fold(j) = 3
                  ampli(j) = kst(3,i)
                  phase(j) = tors3(2,k)
               end if
               write (iout,320)  i,ia,ib,ic,id,bl(ist(2,i)),
     &                           (ampli(k),nint(phase(k)),
     &                           fold(k),k=1,j)
  320          format (i6,3x,4i6,2x,f10.4,1x,3(f8.3,i4,'/',i1))
            end if
         end do
      end if
c
c     list parameters used for atomic partial charges
c
      if (use_charge .or. use_chgdpl) then
         header = .true.
         do i = 1, nion
            ia = iion(i)
            ib = jion(i)
            ic = kion(i)
            if (active(ia) .or. active(ic)) then
               if (header) then
                  header = .false.
                  write (iout,330)
  330             format (/,' Atomic Partial Charge Parameters :',
     &                    /,45x,'Neighbor',3x,'Cutoff',
     &                    /,10x,'Atom Number',13x,'Charge',
     &                       7x,'Site',6x,'Site',/)
               end if
               if (ia.eq.ib .and. ia.eq.ic) then
                  write (iout,340)  i,ia,pchg(i)
  340             format (i6,3x,i6,15x,f10.4)
               else
                  write (iout,350)  i,ia,pchg(i),ib,ic
  350             format (i6,3x,i6,15x,f10.4,5x,i6,4x,i6)
               end if
            end if
         end do
      end if
c
c     list parameters used for bond dipole moments
c
      if (use_dipole .or. use_chgdpl) then
         header = .true.
         do i = 1, ndipole
            ia = idpl(1,i)
            ib = idpl(2,i)
            if (active(ia) .or. active(ib)) then
               if (header) then
                  header = .false.
                  write (iout,360)
  360             format (/,' Bond Dipole Moment Parameters :',
     &                    //,10x,'Atom Numbers',22x,'Dipole',
     &                       3x,'Position',/)
               end if
               write (iout,370)  i,ia,ib,bdpl(i),sdpl(i)
  370          format (i6,3x,2i6,13x,f16.4,f10.4)
            end if
         end do
      end if
c
c     list parameters used for atomic multipole moments
c
      if (use_mpole) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,380)
  380             format (/,' Atomic Multipole Parameters :',
     &                    //,12x,'Atom',4x,'Coordinate Frame',
     &                       ' Definition',7x,'Multipole Moments',/)
               end if
               izaxe = zaxis(i)
               ixaxe = xaxis(i)
               iyaxe = yaxis(i)
               if (izaxe .gt. n)  izaxe = 0
               if (ixaxe .gt. n)  ixaxe = 0
               if (iyaxe .lt. 0)  iyaxe = -iyaxe
               do j = 2, 4
                  pole(j,i) = pole(j,i) / bohr
               end do
               do j = 5, 13
                  pole(j,i) = 3.0d0 * pole(j,i) / bohr**2
               end do
               if (iyaxe .eq. 0) then
                  write (iout,390)  i,ia,izaxe,ixaxe,polaxe(i),
     &                              (pole(j,i),j=1,5),pole(8,i),
     &                              pole(9,i),(pole(j,i),j=11,13)
  390             format (i6,3x,i6,5x,2i6,6x,a8,2x,f9.5,/,48x,3f9.5,
     &                    /,48x,f9.5,/,48x,2f9.5,/,48x,3f9.5)
               else
                  write (iout,400)  i,ia,izaxe,ixaxe,iyaxe,polaxe(i),
     &                              (pole(j,i),j=1,5),pole(8,i),
     &                              pole(9,i),(pole(j,i),j=11,13)
  400             format (i6,3x,i6,2x,3i6,3x,a8,2x,f9.5,/,48x,3f9.5,
     &                    /,48x,f9.5,/,48x,2f9.5,/,48x,3f9.5)
               end if
            end if
         end do
      end if
c
c     list parameters used for dipole polarizability
c
      if (use_polar) then
         header = .true.
         do i = 1, npole
            ia = ipole(i)
            if (active(ia)) then
               if (header) then
                  header = .false.
                  write (iout,410)
  410             format (/,' Dipole Polarizability Parameters :',
     &                    //,10x,'Atom Number',9x,'Alpha',8x,
     &                       'Polarization Group',/)
               end if
               write (iout,420)  i,ia,polarity(i),
     &                           (ip11(j,ia),j=1,np11(ia))
  420          format (i6,3x,i6,10x,f10.4,5x,20i6)
            end if
         end do
      end if
c
c     list parameters used for empirical solvation
c
      if (use_solv) then
         header = .true.
         k = 0
         do i = 1, n
            if (active(i)) then
               k = k + 1
               if (header) then
                  header = .false.
                  write (iout,430)
  430             format (/,' Empirical Solvation Parameters :',
     &                    //,10x,'Atom Number',13x,'Radius',
     &                       3x,'ASP Value',/)
               end if
               write (iout,440)  k,i,rsolv(i),asolv(i)
  440          format (i6,3x,i6,15x,2f10.4)
            end if
         end do
      end if
      return
      end
