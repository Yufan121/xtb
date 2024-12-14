! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

module xtb_readparam
   use xtb_xtb_data ! 
   use xtb_xtb_gfn0
   use xtb_xtb_gfn1
   use xtb_xtb_gfn2
   use xtb_xtb_gfn2_per_atom
   use xtb_paramset
   use xtb_disp_dftd3param, only : copy_c6, reference_c6
   use xtb_disp_dftd4, only : newD4Model, p_refq_gfn2xtb, p_refq_goedecker
   use xtb_type_molecule, only : TMolecule

contains

subroutine read2Param &
      (env, iunit,iunitPerAtom ,globpar,xtbData,initialize, mol) ! Yufan added mol
   use xtb_mctc_accuracy, only : wp

   use xtb_readin, only : getline => strip_line
   use xtb_type_environment, only : TEnvironment
   use xtb_type_param, only : TxTBParameter
   use xtb_param_paulingen, only : paulingEN
   use xtb_param_atomicrad, only : atomicRad
   use xtb_mctc_param, only: chemical_hardness

   implicit none

   logical, parameter :: debug = .false.

   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: iunit, iunitPerAtom 
   type(TxTBParameter), intent(inout) :: globpar
   type(TxTBData), intent(out) :: xtbData
   type(TMolecule), intent(in) :: mol ! Yufan
   logical, intent(in) :: initialize

   character(len=1), parameter :: equal = '='
   character(len=1), parameter :: space = ' '
   character(len=1), parameter :: flag = '$'
   character(len=*), parameter :: flag_end = flag//'end'
   integer, parameter :: p_str_length = 48
   integer, parameter :: p_arg_length = 24

   integer, parameter :: max_elem = 118
   integer :: nShell(max_elem)
   integer :: principalQuantumNumber(10,max_elem)
   integer :: angShell(10,max_elem)
   real(wp) :: shellPoly(1:4,max_elem)
   real(wp) :: selfEnergy(10,max_elem)
   real(wp) :: slaterExponent(10,max_elem)
   real(wp) :: thirdOrderAtom(max_elem)
   real(wp) :: atomicHardness(max_elem)
   real(wp) :: shellHardness(1:4,max_elem)
   real(wp) :: electronegativity(max_elem)
   real(wp) :: repAlpha(max_elem)
   real(wp) :: repZeff(max_elem)
   real(wp) :: halogenBond(max_elem)
   real(wp) :: dipKernel(max_elem)
   real(wp) :: quadKernel(max_elem)
   real(wp) :: eeqkcn(max_elem)
   real(wp) :: chargeWidth(max_elem)
   real(wp) :: kqat2(max_elem)
   real(wp) :: kqat(0:2,max_elem)
   real(wp) :: kpair(max_elem,max_elem)
   real(wp) :: kcnat(0:2,max_elem)
   real(wp) :: eeqEN(max_elem)
   real(wp) :: kExpLight, kExp
   character(len=30) :: timestp(max_elem)
   type(dftd_parameter) :: disp

   character(len=:), allocatable :: line

   integer :: mShell, iSh, jSh
   integer :: level
   integer :: err
   logical :: newFormat

   disp = dftd_parameter(s6=1.0_wp, s8=0.0_wp, a1=0.0_wp, a2=0.0_wp, s9=0.0_wp)
   globpar = TxTBParameter()

   principalQuantumNumber=0
   angShell  =0
   nShell  =0
   selfEnergy=0.0_wp
   slaterExponent=0.0_wp
   shellPoly =0.0_wp
   halogenBond = 0.0_wp
   atomicHardness = chemical_hardness(:max_elem)
   shellHardness  =0.0_wp
   thirdOrderAtom  =0.0_wp
   kqat2 =0.0_wp
   eeqkcn=0.0_wp
   eeqen =0.0_wp
   kcnat =0.0_wp

   electronegativity = paulingEN(:max_elem)
   repAlpha = 0.0_wp
   repZeff = 0.0_wp

   dipKernel = 0.0_wp ! read values are scaled by 0.01
   quadKernel = 0.0_wp !  "     "     "    "    "   "

   kpair =1.0_wp


   !> read all to xtb_xtb_gfn2
   kExpLight = 0.0_wp
   level = -1
   newFormat = .false.
   call getline(iunit,line,err) ! iunit is filename
   if (debug) print'(">",a)',line
   readgroups: do
      if (index(line,flag).eq.1) then 
         select case(line(2:))
         case('info')
            newFormat = .true.
            call read_info
         case('globpar')
            call read_globpar
         case('pairpar')
            call read_pairpar
         case default
            if (index(line,'Z').eq.2) then 
               call read_elempar !!!! read_elempar is called here !!!!!
            else
               call getline(iunit,line,err)
               if (debug) print'(">",a)',line
            endif
         end select
      else
         call getline(iunit,line,err)
         if (debug) print'(">",a)',line
      endif
      if (err.ne.0) exit readgroups
      !if (index(line,flag_end).gt.0) exit readgroups
   enddo readgroups

   if (.not.newFormat) then
      call env%error("Old format parameter file is not supported anymore")
   end if

   call setpair(level, kpair) 



   mShell = maxval(nShell)
   xtbData%level = level ! set to level
   xtbData%nShell = nShell ! 1-dim array
   xtbData%ipeashift = globpar%ipeashift * 0.1_wp

   ! Repulsion
   call init(xtbData%repulsion, kExp, kExpLight, 1.0_wp, globpar%renscale, &
      & repAlpha, repZeff, electronegativity)    ! Yufan, element parameter repAlpha repZeff

   ! Coulomb
   xtbData%coulomb%gExp = globpar%alphaj
   xtbData%coulomb%chemicalHardness = atomicHardness(:max_elem)  ! Yufan: element parameter, atomicHardness
   allocate(xtbData%coulomb%shellHardness(mShell, max_elem))
   call setGFN1ShellHardness(xtbData%coulomb%shellHardness, nShell, angShell, & ! Yufan: element parameter, shellHardness
      & atomicHardness, shellHardness)
   xtbData%coulomb%thirdOrderAtom = thirdOrderAtom(:max_elem)  ! Yufan: element parameter
   xtbData%coulomb%electronegativity = eeqEN(:max_elem)
   xtbData%coulomb%kCN = eeqkCN(:max_elem)      ! Yufan: element parameter kCN
   xtbData%coulomb%chargeWidth = chargeWidth(:max_elem)

   ! Dispersion
   xtbData%dispersion%dpar = disp
   xtbData%dispersion%g_a = 3.0_wp
   xtbData%dispersion%g_c = 2.0_wp
   xtbData%dispersion%wf  = 6.0_wp

   ! Hamiltonian
   xtbData%hamiltonian%angShell = angShell(:mShell, :)

   do iSh = 0, 3
      do jSh = 0, 3
         xtbData%hamiltonian%kScale(jSh, iSh) = 0.5_wp * (globpar%kShell(iSh) &
            & + globpar%kShell(jSh))
      end do
   end do
   if (globpar%ksp > 0.0_wp) then
      xtbData%hamiltonian%kScale(0,1) = globpar%ksp
      xtbData%hamiltonian%kScale(1,0) = globpar%ksp
   end if
   if (globpar%ksd > 0.0_wp) then
      xtbData%hamiltonian%kScale(0,2) = globpar%ksd
      xtbData%hamiltonian%kScale(2,0) = globpar%ksd
   end if
   if (globpar%kpd > 0.0_wp) then
      xtbData%hamiltonian%kScale(1,2) = globpar%kpd
      xtbData%hamiltonian%kScale(2,1) = globpar%kpd
   end if
   xtbData%hamiltonian%kDiff = globpar%kDiff

   ! enScale, not in gfn2
   do iSh = 0, 3
      do jSh = 0, 3
         xtbData%hamiltonian%enScale(jSh, iSh) = 0.005_wp * (globpar%enShell(iSh) &
            & + globpar%enShell(jSh))
      end do
   end do
   xtbData%hamiltonian%enScale4 = globpar%enscale4

   xtbData%hamiltonian%electronegativity = electronegativity(:)
   xtbData%hamiltonian%atomicRad = atomicRad(:)
   xtbData%hamiltonian%shellPoly = shellPoly(:, :) ! Yufan: element parameter
   xtbData%hamiltonian%pairParam = kpair(:, :)
   xtbData%hamiltonian%selfEnergy = selfEnergy(:mShell, :)  ! Yufan: element parameter, selfEnergy
   xtbData%hamiltonian%slaterExponent = slaterExponent(:mShell, :)   ! Yufan: element parameter, slaterExponent
   xtbData%hamiltonian%principalQuantumNumber = principalQuantumNumber(:mShell, :)

   allocate(xtbData%hamiltonian%valenceShell(mShell, max_elem))
   call generateValenceShellData(xtbData%hamiltonian%valenceShell, &
      & xtbData%nShell, xtbData%hamiltonian%angShell)



   select case(level) ! GFN 0 or 1 or 2
   case(0)
      ! Hamiltonian
      xtbData%hamiltonian%wExp = 1.0_wp

      allocate(xtbData%hamiltonian%kCN(mShell, max_elem))
      call angToShellData(xtbData%hamiltonian%kCN, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, kcnat)

      allocate(xtbData%hamiltonian%kQShell(mShell, max_elem))
      call angToShellData(xtbData%hamiltonian%kQShell, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, kqat)
      xtbData%hamiltonian%kQAtom = kqat2(:)

      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN2ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN0NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%srb)
      xtbData%srb%shift = globpar%srbshift
      xtbData%srb%prefactor = globpar%srbpre
      xtbData%srb%steepness = globpar%srbexp
      xtbData%srb%enScale = globpar%srbken

      ! Dispersion
      call newD4Model(xtbData%dispersion%dispm, xtbData%dispersion%g_a, &
         & xtbData%dispersion%g_c, p_refq_goedecker)

   case(1)
      ! Halogen
      allocate(xtbData%halogen)
      call init(xtbData%halogen, globpar%xbrad, globpar%xbdamp, halogenBond)

      ! Hamiltonian
      xtbData%hamiltonian%wExp = 0.0_wp

      allocate(xtbData%hamiltonian%kCN(mShell, max_elem))
      call setGFN1kCN(xtbData%hamiltonian%kCN, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, xtbData%hamiltonian%selfEnergy, &
         & globpar%cnshell)

      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN1ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN1NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

      ! Dispersion
      if (.not.allocated(reference_c6)) call copy_c6(reference_c6)

   case(2) !!!!! GFN2
      ! Coulomb
      if (any(globpar%gam3shell > 0.0_wp)) then
         allocate(xtbData%Coulomb%thirdOrderShell(mShell, max_elem))
         call setGFN2ThirdOrderShell(xtbData%Coulomb%thirdOrderShell, &
            & xtbData%nShell, xtbData%hamiltonian%angShell, thirdOrderAtom, &
            & globpar%gam3shell)
      end if

      ! Multipole
      allocate(xtbData%multipole)
      call init(xtbData%multipole, globpar%aesshift, globpar%aesexp, &
         & globpar%aesrmax, globpar%aesdmp3, globpar%aesdmp5, &
         & dipKernel, quadKernel) ! Yufan: element parameter dipKernel, quadKernel

      ! Hamiltonian
      xtbData%hamiltonian%wExp = 0.5_wp

      allocate(xtbData%hamiltonian%kCN(mShell, max_elem))
      call angToShellData(xtbData%hamiltonian%kCN, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, kcnat)                      ! Yufan: element parameter kCN, kCN holds all elements' data
         ! this kcn is processed using "kcn" from param file

      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN2ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN2NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

      ! Dispersion
      call newD4Model(xtbData%dispersion%dispm, xtbData%dispersion%g_a, &
         & xtbData%dispersion%g_c, p_refq_gfn2xtb)



      !!!!!!!!! Yufan: lines loop structure for PerAtom !!!!!!!!!
      call getline(iunitPerAtom,line,err) ! iunit is filename
      allocate(xtbData%perAtomXtbData)
      !!! Allocate fields of xtbData%perAtomXtbData
      allocate(xtbData%perAtomXtbData%hamiltonian%valenceShell(mShell, max_elem))
      if (any(globpar%gam3shell > 0.0_wp)) then
         allocate(xtbData%perAtomXtbData%Coulomb%thirdOrderShell(mShell, max_elem))
      end if
      allocate(xtbData%perAtomXtbData%multipole)
      allocate(xtbData%perAtomXtbData%hamiltonian%kCN(mShell, max_elem))
      allocate(xtbData%perAtomXtbData%hamiltonian%referenceOcc(mShell, max_elem))
      allocate(xtbData%perAtomXtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      !!! copy fields of xtbData to xtbData%perAtomXtbData
      ! xtbData%perAtomXtbData%level = xtbData%level ! Should be read in loop below
      ! xtbData%perAtomXtbData%name = xtbData%name ! Should be read in loop below
      ! xtbData%perAtomXtbData%srb = xtbData%srb ! not for level 2
      ! xtbData%perAtomXtbData%halogen = xtbData%halogen ! not for level 2
      xtbData%perAtomXtbData%coulomb = xtbData%coulomb 
      xtbData%perAtomXtbData%doi = xtbData%doi
      xtbData%perAtomXtbData%hamiltonian = xtbData%hamiltonian
      xtbData%perAtomXtbData%dispersion = xtbData%dispersion
      xtbData%perAtomXtbData%ipeashift = xtbData%ipeashift
      xtbData%perAtomXtbData%multipole = xtbData%multipole
      xtbData%perAtomXtbData%nshell = xtbData%nshell
      xtbData%perAtomXtbData%repulsion = xtbData%repulsion

      !!! Allocate xtb_gfn2_per_atom 
      ! Allocate thirdOrderAtom with size mol%n
      allocate(thirdOrderAtomPerAtom(mol%n))
      ! Allocate repAlpha with size mol%n
      allocate(repAlphaPerAtom(mol%n))
      ! Allocate dipKernel with size mol%n
      allocate(dipKernelPerAtom(mol%n))
      ! Allocate quadKernel with size mol%n
      allocate(quadKernelPerAtom(mol%n)) 
      ! Allocate Effective nuclear charge with size mol%n
      allocate(repZeffPerAtom(mol%n)) 
      ! Allocate nshell 
      allocate(nShellPerAtom(mol%n))
      ! Allocate ElemId
      allocate(ElemIdPerAtom(mol%n))
      ! Allocate ENPerAtom
      allocate(electronegativityPerAtom(mol%n))
      ! Allocate atomicHardnessPerAtom
      allocate(atomicHardnessPerAtom(mol%n))
      ! Allocate shellPoly with shape (4, mol%n)
      allocate(shellPolyPerAtom(4, mol%n))
      ! Allocate selfEnergy with shape (3, mol%n)
      allocate(selfEnergyPerAtom(3, mol%n))
      ! Allocate slaterExponent with shape (3, mol%n)
      allocate(slaterExponentPerAtom(3, mol%n))
      ! Allocate kCN with shape (4, mol%n)
      allocate(kCNPerAtom(4, mol%n))
      ! Allocate Shell Hardness with shape (3, mol%n)
      allocate(shellHardnessPerAtom(3, mol%n))
      ! Allocate angshell with shape (3, mol)
      allocate(angshellPerAtom(3, mol%n))
      ! Allocate kcnat with shape (3, mol%n)
      allocate(kcnatPerAtom(3, mol%n))
      ! Allocate principalQuantumNumberPerAtom (3, mol%n)
      allocate(principalQuantumNumberPerAtom(3, mol%n))


      !!! loop structure, load file prarm to xtb_xtb_gfn2_perAtom (overwrite original file)
      if (debug) print'(">",a)',line
      readgroupsperatom: do
         if (index(line,flag).eq.1) then 
            select case(line(2:))
            case('info')
               call read_info_per_atom
            ! case('globpar')
            !    call read_globpar
            ! case('pairpar')
            !    call read_pairpar
            case default
               if (index(line,'Z').eq.2) then
                  call read_elempar_per_atom !!!! read_elempar is called here !!!!!
               else
                  call getline(iunitPerAtom,line,err)
                  if (debug) print'(">",a)',line
               endif
            end select
         else
            call getline(iunitPerAtom,line,err)
            if (debug) print'(">",a)',line
         endif
         if (err.ne.0) exit readgroupsperatom
         !if (index(line,flag_end).gt.0) exit readgroupsperatom
      enddo readgroupsperatom
      xtbData%perAtomXtbData%level = level
      if (.not.newFormat) then
         call env%error("Old format parameter file is not supported anymore")
      end if
      !!! EndYufan loop structure

      ! Verify that nAtom and ElemIdPerAtom match
      call checkElemIdPerAtomMatch(env, ElemIdPerAtom, mol)

      !!! Assign some values from xtb_xtb_gfn2_perAtom to perAtomXtbData
      xtbData%perAtomXtbData%nShell = nShellPerAtom ! 1-dim array
      
      xtbData%perAtomXtbData%hamiltonian%electronegativity = electronegativityPerAtom(:)
      ! xtbData%perAtomXtbData%hamiltonian%atomicRad = atomicRad(:)
      xtbData%perAtomXtbData%hamiltonian%shellPoly = shellPolyPerAtom(:, :) ! Yufan: element parameter
      ! xtbData%perAtomXtbData%hamiltonian%pairParam = kpair(:, :)
      xtbData%perAtomXtbData%hamiltonian%selfEnergy = selfEnergyPerAtom(:mShell, :)  ! Yufan: element parameter, selfEnergy
      xtbData%perAtomXtbData%hamiltonian%slaterExponent = slaterExponentPerAtom(:mShell, :)   ! Yufan: element parameter, slaterExponent
      xtbData%perAtomXtbData%hamiltonian%principalQuantumNumber = principalQuantumNumberPerAtom(:mShell, :)

      xtbData%perAtomXtbData%hamiltonian%angShell = angShellPerAtom(:mShell, :)


      !!! Reallocate perAtomXtbData and assign value to it

      ! Repulsion
      ! Yufan, element parameter repAlpha repZeff electronegativity (temporarily use the fixed), does not change irrelavant
      ! Change self%electronegativity, self%alpha, self%zeff and other global pars
      call init(xtbData%perAtomXtbData%repulsion, kExp, kExpLight, 1.0_wp, globpar%renscale, &
         & repAlphaPerAtom, repZeffPerAtom, electronegativityPerAtom)    
      
      xtbData%perAtomXtbData%coulomb%chemicalHardness = atomicHardnessPerAtom(:mol%n)  ! Yufan: element parameter, atomicHardness


      ! Done, shell hardness tweak
      deallocate(xtbData%perAtomXtbData%coulomb%shellHardness)
      allocate(xtbData%perAtomXtbData%coulomb%shellHardness(mShell, mol%n))
      call setGFN1ShellHardnessPerAtom(xtbData%perAtomXtbData%coulomb%shellHardness, nShellPerAtom, angShellPerAtom, & ! Yufan: element parameter, shellHardness
         & atomicHardnessPerAtom, atomicHardness, shellHardnessPerAtom, shellHardness, ElemIdPerAtom)
      
      xtbData%perAtomXtbData%coulomb%thirdOrderAtom = thirdOrderAtomPerAtom(:mol%n)  ! Yufan: element parameter
      ! xtbData%coulomb%electronegativity = eeqEN(:max_elem) ! neglect, no data for gfn2
      ! xtbData%coulomb%kCN = eeqkCN(:max_elem)      ! neglect, no data for gfn2

      deallocate(xtbData%perAtomXtbData%hamiltonian%valenceShell)       ! May not be needed anymore
      allocate(xtbData%perAtomXtbData%hamiltonian%valenceShell(mShell, mol%n))
      call generateValenceShellData(xtbData%perAtomXtbData%hamiltonian%valenceShell, &
         & xtbData%perAtomXtbData%nShell, xtbData%perAtomXtbData%hamiltonian%angShell)


      ! Coulomb
      if (any(globpar%gam3shell > 0.0_wp)) then
         deallocate(xtbData%perAtomXtbData%coulomb%thirdOrderShell)
         allocate(xtbData%perAtomXtbData%coulomb%thirdOrderShell(mShell, mol%n))
         call setGFN2ThirdOrderShellPerAtom(xtbData%perAtomXtbData%Coulomb%thirdOrderShell, & ! 
            & xtbData%perAtomXtbData%nShell, xtbData%perAtomXtbData%hamiltonian%angShell, thirdOrderAtomPerAtom, &
            & globpar%gam3shell, ElemIdPerAtom)
      end if

      ! Multipole
      ! deallocate(xtbData%perAtomXtbData)
      ! allocate(xtbData%perAtomXtbData%multipole)  
      call init(xtbData%perAtomXtbData%multipole, globpar%aesshift, globpar%aesexp, &
         & globpar%aesrmax, globpar%aesdmp3, globpar%aesdmp5, &
         & dipKernelPerAtom, quadKernelPerAtom, ElemIdPerAtom)

      ! Hamiltonian
      xtbData%perAtomXtbData%hamiltonian%wExp = 0.5_wp

      deallocate(xtbData%perAtomXtbData%hamiltonian%kCN)
      allocate(xtbData%perAtomXtbData%hamiltonian%kCN(mShell, mol%n))
      call angToShellData(xtbData%perAtomXtbData%hamiltonian%kCN, xtbData%perAtomXtbData%nShell, &
         & xtbData%perAtomXtbData%hamiltonian%angShell, kcnatPerAtom)  ! no need for new subroutine

      ! setGFN2ReferenceOcc, just retrive using index, no element_specific params
      deallocate(xtbData%perAtomXtbData%hamiltonian%referenceOcc)
      allocate(xtbData%perAtomXtbData%hamiltonian%referenceOcc(mShell, mol%n))
      call setGFN2ReferenceOccPerAtom(xtbData%perAtomXtbData%hamiltonian, xtbData%perAtomXtbData%nShell, ElemIdPerAtom)

      deallocate(xtbData%perAtomXtbData%hamiltonian%numberOfPrimitives)
      allocate(xtbData%perAtomXtbData%hamiltonian%numberOfPrimitives(mShell, mol%n))
      call setGFN2NumberOfPrimitivesPerAtom(xtbData%perAtomXtbData%hamiltonian, xtbData%perAtomXtbData%nShell)

      ! Dispersion
      call newD4Model(xtbData%perAtomXtbData%dispersion%dispm, xtbData%dispersion%g_a, &
         & xtbData%dispersion%g_c, p_refq_gfn2xtb)

      
   end select

contains

! subroutine read_per_atom_param(env, iunit) ! Yufan
!    use xtb_mctc_strings
!    use xtb_readin
!    implicit none
!    type(TEnvironment), intent(inout) :: env
!    integer, intent(in) :: iunit
!    character(len=:), allocatable :: key, val, line
!    integer :: iz, ie, err, natom, nelement
!    logical :: in_section

!    in_section = .false.

!    do
!       call getline(iunit, line, err)
!       if (err.ne.0) exit
!       if (index(line, '$') == 1) then
!             if (index(line, 'info') > 0) then
!                call read_info_atom_par(env, iunit)
!             else if (index(line, 'Z=') == 2) then
!                call read_atom_param(env, iunit, line)
!             end if
!       end if
!    end do
! end subroutine read_per_atom_param

! subroutine read_atom_param(env, iunit, line) ! Yufan
!     use xtb_mctc_strings
!     use xtb_readin
!     implicit none
!     type(TEnvironment), intent(inout) :: env
!     integer, intent(in) :: iunit
!     character(len=*), intent(in) :: line
!     character(len=:), allocatable :: key, val
!     integer :: iz, ie, err

!     if (getValue(env, line(4:5), iz)) then
!         timestp(iz) = line(7:len_trim(line))
!         do
!             call getline(iunit, line, err)
!             if (err.ne.0) exit
!             if (index(line, '$end') > 0) exit
!             ie = index(line, '=')
!             if (line == '') cycle
!             if (ie == 0) cycle
!             key = lowercase(trim(line(:ie-1)))
!             val = trim(adjustl(line(ie+1:)))
!             call gfn_atom_param(key, val, iz)
!         end do
!     else
!         call getline(iunit, line, err)
!     end if
! end subroutine read_atom_param

! subroutine gfn_atom_param(key, val, iz) ! Yufan
!     use xtb_mctc_strings
!     use xtb_readin
!     implicit none
!     character(len=*), intent(in) :: key, val
!     integer, intent(in) :: iz
!     integer :: narg
!     character(len=p_str_length), dimension(p_arg_length) :: argv
!     integer :: i, ii
!     integer :: idum
!     real(wp) :: ddum

!     select case(key)
!     case default
!         call env%warning('Unknown key "'//key//'" for "$Z"')
!     case('ele_id')
!         if (getValue(env, val, idum)) ele_id(iz) = idum
!     case('ao')
!         if (mod(len(val), 2) == 0) then
!             nShell(iz) = len(val) / 2
!             do i = 1, nShell(iz)
!                 ii = 2 * i - 1
!                 if (getValue(env, val(ii:ii), idum)) then
!                     principalQuantumNumber(i, iz) = idum
!                     select case(val(ii+1:ii+1))
!                     case('s'); angShell(i, iz) = 0
!                     case('p'); angShell(i, iz) = 1
!                     case('d'); angShell(i, iz) = 2
!                     case('f'); angShell(i, iz) = 3
!                     case('g'); angShell(i, iz) = 4
!                     case('S'); angShell(i, iz) = 0
!                     end select
!                 end if
!             end do
!         end if
!     case('lev')
!         call parse(val, space, argv, narg)
!         if (narg == nShell(iz)) then
!             do i = 1, nShell(iz)
!                 if (getValue(env, trim(argv(i)), ddum)) selfEnergy(i, iz) = ddum
!             end do
!         end if
!     case('exp')
!         call parse(val, space, argv, narg)
!         if (narg == nShell(iz)) then
!             do i = 1, nShell(iz)
!                 if (getValue(env, trim(argv(i)), ddum)) slaterExponent(i, iz) = ddum  ! slaterExponent
!             end do
!         end if
!     case('gam')
!         if (getValue(env, val, ddum)) atomicHardness(iz) = ddum  ! seems only in gfn1
!     case('gam3')
!         if (getValue(env, val, ddum)) thirdOrderAtom(iz) = ddum * 0.1_wp ! Third order Hubbard derivatives
!     case('kcns')
!         if (getValue(env, val, ddum)) kcnat(0, iz) = ddum * 0.1_wp ! coordination number dependence
!     case('kcnp')
!         if (getValue(env, val, ddum)) kcnat(1, iz) = ddum * 0.1_wp
!     case('dpol')
!         if (getValue(env, val, ddum)) dipKernel(iz) = ddum * 0.01_wp ! multipole
!     case('qpol')
!         if (getValue(env, val, ddum)) quadKernel(iz) = ddum * 0.01_wp ! multipole
!     case('repa')
!         if (getValue(env, val, ddum)) repAlpha(iz) = ddum ! repulsion
!     case('repb')
!         if (getValue(env, val, ddum)) repZeff(iz) = ddum ! repulsion
!     case('polys')
!         if (getValue(env, val, ddum)) shellPoly(1, iz) = ddum ! shell polynomial for hamiltonian
!     case('polyp')
!         if (getValue(env, val, ddum)) shellPoly(2, iz) = ddum
!     case('polyd')
!         if (getValue(env, val, ddum)) shellPoly(3, iz) = ddum
!     case('polyf')
!         if (getValue(env, val, ddum)) shellPoly(4, iz) = ddum
!     case('lpars')
!         if (getValue(env, val, ddum)) shellHardness(1, iz) = ddum * 0.1_wp ! shell poly for electrostatics
!     case('lparp')
!         if (getValue(env, val, ddum)) shellHardness(2, iz) = ddum * 0.1_wp
!     case('lpard')
!         if (getValue(env, val, ddum)) shellHardness(3, iz) = ddum * 0.1_wp
!     case('lparf')
!         if (getValue(env, val, ddum)) shellHardness(4, iz) = ddum * 0.1_wp
!     end select
! end subroutine gfn_atom_param

! subroutine read_info_atom_par(env, iunit) ! Yufan
!     use xtb_readin, only : getValue
!     implicit none
!     type(TEnvironment), intent(inout) :: env
!     integer, intent(in) :: iunit
!     character(len=:), allocatable :: key, val, line
!     integer :: ie, err, idum

!     do
!         call getline(iunit, line, err)
!         if (err.ne.0) exit
!         if (index(line, '$end') > 0) exit
!         ie = index(line, ' ')
!         if (line == '') cycle
!         if (ie == 0) cycle
!         key = trim(line(:ie-1))
!         val = trim(adjustl(line(ie+1:)))
!         select case(key)
!         case('level')
!             if (getValue(env, val, idum)) level = idum
!         case('name')
!             xtbData%name = val
!         case('natom')
!             if (getValue(env, val, idum)) natom = idum
!         case('nelement')
!             if (getValue(env, val, idum)) nelement = idum
!         case default
!             call env%warning('Unknown key "'//key//'" for "$info"')
!         end select
!     end do
! end subroutine read_info_atom_par


subroutine read_info
   use xtb_readin, only : getValue
   character(len=:), allocatable :: key, val
   integer :: ie, idum ! Yufan: these are only intermediate values
   do
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) exit
      if (index(line,flag).gt.0) exit

      ie = index(line,space)
      if (line.eq.'') cycle ! skip empty lines
      if (ie.eq.0) cycle

      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))

      select case(key)
      case default
         call env%warning("Unknown key '"//key//"' for '"//flag//"info'")
      case('level')
         if (getValue(env,val,idum)) level = idum
      case('name')
         xtbData%name = val
      case('doi')
         xtbData%doi = val
      end select
   enddo
end subroutine read_info

subroutine read_info_per_atom
   use xtb_readin, only : getValue
   character(len=:), allocatable :: key, val
   integer :: ie, idum ! Yufan: these are only intermediate values
   do
      call getline(iunitPerAtom,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) exit
      if (index(line,flag).gt.0) exit

      ie = index(line,space)
      if (line.eq.'') cycle ! skip empty lines
      if (ie.eq.0) cycle

      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))

      select case(key)
      case default
         call env%warning("Unknown key '"//key//"' for '"//flag//"info'")
      case('level')
         if (getValue(env,val,idum)) level = idum ! Yufan: this level is only intermediate, so use the same one for perAtom
      case('name')
         xtbData%perAtomXtbData%name = val
      case('doi')
         xtbData%perAtomXtbData%doi = val
      ! case('ele_id')
      !    xtbData%perAtomXtbData%ele_id = val
      case('natom')
         if (getValue(env,val,idum)) xtbData%perAtomXtbData%nAtom = idum
      case('nelement')
         if (getValue(env,val,idum)) xtbData%perAtomXtbData%nElement = idum

      end select
   enddo
end subroutine read_info_per_atom


subroutine read_globpar
   implicit none
   character(len=:), allocatable :: key, val
   integer :: ie
   do
      call getline(iunit,line,err) ! where is the file
      if (debug) print'("->",a)',line
      if (err.ne.0) exit
      if (index(line,flag).gt.0) exit

      ie = index(line,space)
      if (line.eq.'') cycle ! skip empty lines
      if (ie.eq.0) cycle

      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:))) ! Yufan: already got data in str format

      call gfn_globpar(key,val,globpar)

   enddo
end subroutine read_globpar

subroutine gfn_globpar(key,val,globpar) ! Yufan : this is only for converting to real and set to globpar 
   use xtb_readin, only : getValue
   implicit none
   character(len=*), intent(in) :: key, val
   type(TxTBParameter), intent(inout) :: globpar
   real(wp) :: ddum
   select case(key)
   case default
      call env%warning("Unknown key '"//key//"' for '"//flag//"globpar'")
   case('kexp'); if (getValue(env,val,ddum)) kExp = ddum
   case('kexplight'); if (getValue(env,val,ddum)) kExpLight = ddum
   case('ks'); if (getValue(env,val,ddum)) globpar%kshell(0) = ddum ! read from env
   case('kp'); if (getValue(env,val,ddum)) globpar%kshell(1) = ddum
   case('kd'); if (getValue(env,val,ddum)) globpar%kshell(2) = ddum
   case('kf'); if (getValue(env,val,ddum)) globpar%kshell(3) = ddum
   case('ksp'); if (getValue(env,val,ddum)) globpar%ksp = ddum
   case('ksd'); if (getValue(env,val,ddum)) globpar%ksd = ddum
   case('kpd'); if (getValue(env,val,ddum)) globpar%kpd = ddum
   case('kdiff'); if (getValue(env,val,ddum)) globpar%kdiff = ddum
   case('kdiffa'); if (getValue(env,val,ddum)) globpar%kdiffa = ddum
   case('kdiffb'); if (getValue(env,val,ddum)) globpar%kdiffb = ddum
   case('ens'); if (getValue(env,val,ddum)) globpar%enshell(0) = ddum
   case('enp'); if (getValue(env,val,ddum)) globpar%enshell(1) = ddum
   case('end'); if (getValue(env,val,ddum)) globpar%enshell(2) = ddum
   case('enf'); if (getValue(env,val,ddum)) globpar%enshell(3) = ddum
   case('enscale'); if (getValue(env,val,ddum)) globpar%enshell = ddum
   case('enscale4'); if (getValue(env,val,ddum)) globpar%enscale4 = ddum
   case('renscale'); if (getValue(env,val,ddum)) globpar%renscale = ddum
   case('cns'); if (getValue(env,val,ddum)) globpar%cnshell(:, 0) = ddum
   case('cnp'); if (getValue(env,val,ddum)) globpar%cnshell(:, 1) = ddum
   case('cnd'); if (getValue(env,val,ddum)) globpar%cnshell(:, 2) = ddum
   case('cnf'); if (getValue(env,val,ddum)) globpar%cnshell(:, 3) = ddum
   case('cnd1'); if (getValue(env,val,ddum)) globpar%cnshell(1, 2) = ddum
   case('cnd2'); if (getValue(env,val,ddum)) globpar%cnshell(2, 2) = ddum
   case('gam3s'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 0) = ddum
   case('gam3p'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 1) = ddum
   case('gam3d'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 2) = ddum
   case('gam3f'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 3) = ddum
   case('gam3d1'); if (getValue(env,val,ddum)) globpar%gam3shell(1, 2) = ddum
   case('gam3d2'); if (getValue(env,val,ddum)) globpar%gam3shell(2, 2) = ddum
   case('srbshift'); if (getValue(env,val,ddum)) globpar%srbshift = ddum
   case('srbpre'); if (getValue(env,val,ddum)) globpar%srbpre = ddum
   case('srbexp'); if (getValue(env,val,ddum)) globpar%srbexp = ddum
   case('srbken'); if (getValue(env,val,ddum)) globpar%srbken = ddum
   case('wllscal'); if (getValue(env,val,ddum)) globpar%wllscal = ddum
   case('ipeashift'); if (getValue(env,val,ddum)) globpar%ipeashift = ddum
   case('gscal'); if (getValue(env,val,ddum)) globpar%gscal = ddum
   case('zcnf'); if (getValue(env,val,ddum)) globpar%zcnf = ddum
   case('tscal'); if (getValue(env,val,ddum)) globpar%tscal = ddum
   case('kcn'); if (getValue(env,val,ddum)) globpar%kcn = ddum
   case('fpol'); if (getValue(env,val,ddum)) globpar%fpol = ddum
   case('ken'); if (getValue(env,val,ddum)) globpar%ken = ddum
   case('lshift'); if (getValue(env,val,ddum)) globpar%lshift = ddum
   case('lshifta'); if (getValue(env,val,ddum)) globpar%lshifta = ddum
   case('split'); if (getValue(env,val,ddum)) globpar%split = ddum
   case('zqf'); if (getValue(env,val,ddum)) globpar%zqf = ddum
   case('alphaj'); if (getValue(env,val,ddum)) globpar%alphaj = ddum
   case('kexpo'); if (getValue(env,val,ddum)) globpar%kexpo = ddum
   case('dispa'); if (getValue(env,val,ddum)) globpar%dispa = ddum
   case('dispb'); if (getValue(env,val,ddum)) globpar%dispb = ddum
   case('dispc'); if (getValue(env,val,ddum)) globpar%dispc = ddum
   case('dispatm'); if (getValue(env,val,ddum)) globpar%dispatm = ddum
   case('a1'); if (getValue(env,val,ddum)) disp%a1 = ddum
   case('a2'); if (getValue(env,val,ddum)) disp%a2 = ddum
   case('s6'); if (getValue(env,val,ddum)) disp%s6 = ddum
   case('s8'); if (getValue(env,val,ddum)) disp%s8 = ddum
   case('s9'); if (getValue(env,val,ddum)) disp%s9 = ddum
   case('xbdamp'); if (getValue(env,val,ddum)) globpar%xbdamp = ddum
   case('xbrad'); if (getValue(env,val,ddum)) globpar%xbrad = ddum
   case('aesdmp3'); if (getValue(env,val,ddum)) globpar%aesdmp3 = ddum
   case('aesdmp5'); if (getValue(env,val,ddum)) globpar%aesdmp5 = ddum
   case('aesshift'); if (getValue(env,val,ddum)) globpar%aesshift = ddum
   case('aesexp'); if (getValue(env,val,ddum)) globpar%aesexp = ddum
   case('aesrmax'); if (getValue(env,val,ddum)) globpar%aesrmax = ddum
   end select
end subroutine gfn_globpar

subroutine read_pairpar
   use xtb_mctc_strings
   use xtb_readin, only : getValue
   implicit none
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   integer  :: iAt,jAt
   real(wp) :: ddum
   do
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) exit
      if (index(line,flag).gt.0) exit

      call parse(line,space,argv,narg)
      if (narg .ne. 3) cycle
      if (getValue(env,trim(argv(1)),iAt) .and. &
         &getValue(env,trim(argv(2)),jAt) .and. &
         &getValue(env,trim(argv(3)),ddum)) then
            kpair(iAt,jAt) = ddum
            kpair(jAt,iAt) = ddum
      endif
   enddo
end subroutine read_pairpar

subroutine read_elempar
   use xtb_mctc_strings
   use xtb_readin
   implicit none
   character(len=:), allocatable :: key, val  
   integer :: iz, ie
   if (getValue(env,line(4:5),iz)) then
      timestp(iz) = line(7:len_trim(line))
      do
         call getline(iunit,line,err)
         if (debug) print'("->",a)',line
         if (err.ne.0) exit
         if (index(line,flag).gt.0) exit

         ie = index(line,equal)
         if (line.eq.'') cycle ! skip empty lines
         if (ie.eq.0) cycle

         key = lowercase(trim(line(:ie-1)))
         val = trim(adjustl(line(ie+1:)))

         call gfn_elempar(key,val,iz) ! call gfn_elempar to get each parameter

      enddo
   else
      call getline(iunit,line,err)
   endif
end subroutine read_elempar



subroutine gfn_elempar(key,val,iz)  ! iz is the element number
   use xtb_mctc_strings
   use xtb_readin
   implicit none
   character(len=*), intent(in) :: key, val
   integer, intent(in) :: iz
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   integer :: i, ii
   integer :: idum
   real(wp) :: ddum

   select case(key)
   case default
      call env%warning("Unknown key '"//key//"' for '"//flag//"Z'")
   case('ao')
      !print'(a,":",a)',key,val
      if (mod(len(val),2).eq.0) then
         nShell(iz) = len(val)/2
         do i = 1, nShell(iz)
            ii = 2*i-1
            !print*,i,ii,val(ii:ii),val(ii+1:ii+1)
            if (getValue(env,val(ii:ii),idum)) then
               principalQuantumNumber(i,iz) = idum
               select case(val(ii+1:ii+1))
               case('s'); angShell(i,iz) = 0
               case('p'); angShell(i,iz) = 1
               case('d'); angShell(i,iz) = 2
               case('f'); angShell(i,iz) = 3
               case('g'); angShell(i,iz) = 4
               case('S'); angShell(i,iz) = 0
               end select
            endif
         enddo
      endif
   case('lev')
      call parse(val,space,argv,narg)
      if (narg .eq. nShell(iz)) then
         do i = 1, nShell(iz)
            if (getValue(env,trim(argv(i)),ddum)) selfEnergy(i,iz) = ddum
         enddo
      endif
   case('exp')
      call parse(val,space,argv,narg)
      if (narg .eq. nShell(iz)) then
         do i = 1, nShell(iz)
            if (getValue(env,trim(argv(i)),ddum)) slaterExponent(i,iz) = ddum
         enddo
      endif
   case('en');  if (getValue(env,val,ddum)) electronegativity(iz)    = ddum
   case('gam'); if (getValue(env,val,ddum)) atomicHardness(iz)   = ddum
   case('xi');  if (getValue(env,val,ddum)) eeqEN(iz) = ddum
   case('alpg'); if (getValue(env,val,ddum)) chargeWidth(iz)   = ddum
   case('gam3');  if (getValue(env,val,ddum)) thirdOrderAtom(iz)    = ddum * 0.1_wp
   case('kappa'); if (getValue(env,val,ddum)) eeqkCN(iz)  = ddum
   case('cxb');   if (getValue(env,val,ddum)) halogenBond(iz) = ddum * 0.1_wp
   case('kqat2'); if (getValue(env,val,ddum)) kqat2(iz)   = ddum
   case('dpol');  if (getValue(env,val,ddum)) dipKernel(iz)   = ddum * 0.01_wp
   case('qpol');  if (getValue(env,val,ddum)) quadKernel(iz)   = ddum * 0.01_wp
   case('repa');  if (getValue(env,val,ddum)) repAlpha(iz)   = ddum
   case('repb');  if (getValue(env,val,ddum)) repZeff(iz)   = ddum
   case('polys'); if (getValue(env,val,ddum)) shellPoly(1,iz) = ddum !!! this line means that the shell poly is not supported?
   case('polyp'); if (getValue(env,val,ddum)) shellPoly(2,iz) = ddum
   case('polyd'); if (getValue(env,val,ddum)) shellPoly(3,iz) = ddum
   case('polyf'); if (getValue(env,val,ddum)) shellPoly(4,iz) = ddum
   case('lpars'); if (getValue(env,val,ddum)) shellHardness(1,iz)  = ddum * 0.1_wp
   case('lparp'); if (getValue(env,val,ddum)) shellHardness(2,iz)  = ddum * 0.1_wp
   case('lpard'); if (getValue(env,val,ddum)) shellHardness(3,iz)  = ddum * 0.1_wp
   case('lparf'); if (getValue(env,val,ddum)) shellHardness(4,iz)  = ddum * 0.1_wp
   case('kcns');  if (getValue(env,val,ddum)) kcnat(0,iz) = ddum * 0.1_wp
   case('kcnp');  if (getValue(env,val,ddum)) kcnat(1,iz) = ddum * 0.1_wp
   case('kcnd');  if (getValue(env,val,ddum)) kcnat(2,iz) = ddum * 0.1_wp
   case('kqs');   if (getValue(env,val,ddum)) kqat(0,iz)  = ddum
   case('kqp');   if (getValue(env,val,ddum)) kqat(1,iz)  = ddum
   case('kqd');   if (getValue(env,val,ddum)) kqat(2,iz)  = ddum
   end select
end subroutine gfn_elempar

subroutine read_elempar_per_atom
   use xtb_mctc_strings
   use xtb_readin
   implicit none
   character(len=:), allocatable :: key, val  
   integer :: iz, ie
   if (getValue(env,line(4:5),iz)) then
      timestp(iz) = line(7:len_trim(line))
      do
         call getline(iunitPerAtom,line,err)
         if (debug) print'("->",a)',line
         if (err.ne.0) exit 
         if (index(line,flag).gt.0) exit

         ie = index(line,equal)
         if (line.eq.'') cycle ! skip empty lines
         if (ie.eq.0) cycle

         key = lowercase(trim(line(:ie-1)))
         val = trim(adjustl(line(ie+1:)))

         call gfn_elempar_per_atom(key,val,iz) ! call gfn_elempar to get each parameter

      enddo
   else
      call getline(iunitPerAtom,line,err)
   endif
end subroutine read_elempar_per_atom


subroutine gfn_elempar_per_atom(key,val,iz)  ! iz is the atomic id
   use xtb_mctc_strings
   use xtb_readin
   implicit none
   character(len=*), intent(in) :: key, val
   integer, intent(in) :: iz
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   integer :: i, ii
   integer :: idum
   real(wp) :: ddum
   ! character(len=32) :: repZeff_str ! remove after use

   select case(key)
   case default
      call env%warning("Unknown key '"//key//"' for '"//flag//"Z'")
   case('ao')                                                         ! Yufan: changed
      !print'(a,":",a)',key,val
      if (mod(len(val),2).eq.0) then
         nShellPerAtom(iz) = len(val)/2                                      ! Yufan: changed
         do i = 1, nShellPerAtom(iz)
            ii = 2*i-1
            !print*,i,ii,val(ii:ii),val(ii+1:ii+1)
            if (getValue(env,val(ii:ii),idum)) then
               principalQuantumNumberPerAtom(i,iz) = idum
               select case(val(ii+1:ii+1))
               case('s'); angShellPerAtom(i,iz) = 0
               case('p'); angShellPerAtom(i,iz) = 1
               case('d'); angShellPerAtom(i,iz) = 2
               case('f'); angShellPerAtom(i,iz) = 3
               case('g'); angShellPerAtom(i,iz) = 4
               case('S'); angShellPerAtom(i,iz) = 0
               end select
            endif
         enddo
      endif
   case('lev')                                                       ! Yufan: changed
      call parse(val,space,argv,narg)
      if (narg .eq. nShellPerAtom(iz)) then
         do i = 1, nShellPerAtom(iz)
            if (getValue(env,trim(argv(i)),ddum)) selfEnergyPerAtom(i,iz) = ddum
         enddo
      endif
   case('exp')                                                       ! Yufan: changed
      call parse(val,space,argv,narg)
      if (narg .eq. nShellPerAtom(iz)) then
         do i = 1, nShellPerAtom(iz)
            if (getValue(env,trim(argv(i)),ddum)) slaterExponentPerAtom(i,iz) = ddum
         enddo
      endif
   case('en');  if (getValue(env,val,ddum)) electronegativityPerAtom(iz)    = ddum ! not in gfn2?  ! Yufan: changed
   case('gam'); if (getValue(env,val,ddum)) atomicHardnessPerAtom(iz)   = ddum ! Yufan: changed
   ! case('xi');  if (getValue(env,val,ddum)) eeqEN(iz) = ddum ! nof in gfn2
   ! case('alpg'); if (getValue(env,val,ddum)) chargeWidth(iz)   = ddum    !  nof in gfn2
   case('gam3');  if (getValue(env,val,ddum)) thirdOrderAtomPerAtom(iz)    = ddum * 0.1_wp   ! Yufan: changed
   ! case('kappa'); if (getValue(env,val,ddum)) eeqkCN(iz)  = ddum
   ! case('cxb');   if (getValue(env,val,ddum)) halogenBond(iz) = ddum * 0.1_wp
   ! case('kqat2'); if (getValue(env,val,ddum)) kqat2(iz)   = ddum      ! not in gfn2
   case('dpol');  if (getValue(env,val,ddum)) dipKernelPerAtom(iz)   = ddum * 0.01_wp ! Yufan: changed
   case('qpol');  if (getValue(env,val,ddum)) quadKernelPerAtom(iz)   = ddum * 0.01_wp ! Yufan: changed
   case('repa');  if (getValue(env,val,ddum)) repAlphaPerAtom(iz)   = ddum ! Yufan: changed
   case('repb');  if (getValue(env,val,ddum)) repZeffPerAtom(iz)   = ddum ! Yufan: changed
   ! case('repb'); if (getValue(env, val, ddum)) then
   !    print *, "Before modification: repZeff(", iz, ") = ", repZeff(iz)
   !    repZeff(iz) = ddum
   !    print *, "After modification: repZeff(", iz, ") = ", repZeff(iz)
   !    call env%warning("val and repZeff(iz) '" // trim(val) // "' for '" // trim(adjustl(repZeff_str)) // "Z'")
   ! end if
   case('polys'); if (getValue(env,val,ddum)) shellPolyPerAtom(1,iz) = ddum ! Yufan: changed
   case('polyp'); if (getValue(env,val,ddum)) shellPolyPerAtom(2,iz) = ddum 
   case('polyd'); if (getValue(env,val,ddum)) shellPolyPerAtom(3,iz) = ddum
   case('polyf'); if (getValue(env,val,ddum)) shellPolyPerAtom(4,iz) = ddum
   case('lpars'); if (getValue(env,val,ddum)) shellHardnessPerAtom(1,iz)  = ddum * 0.1_wp ! Yufan: changed
   case('lparp'); if (getValue(env,val,ddum)) shellHardnessPerAtom(2,iz)  = ddum * 0.1_wp
   case('lpard'); if (getValue(env,val,ddum)) shellHardnessPerAtom(3,iz)  = ddum * 0.1_wp
   case('lparf'); if (getValue(env,val,ddum)) shellHardnessPerAtom(4,iz)  = ddum * 0.1_wp
   case('kcns');  if (getValue(env,val,ddum)) kcnatPerAtom(1,iz) = ddum * 0.1_wp !! Yufan: changed, chaos in existing code...
   case('kcnp');  if (getValue(env,val,ddum)) kcnatPerAtom(2,iz) = ddum * 0.1_wp  !! important, since using allocatable, index changed to start from 1
   case('kcnd');  if (getValue(env,val,ddum)) kcnatPerAtom(3,iz) = ddum * 0.1_wp 
   ! case('kqs');   if (getValue(env,val,ddum)) kqat(0,iz)  = ddum  ! not in gfn2
   ! case('kqp');   if (getValue(env,val,ddum)) kqat(1,iz)  = ddum  ! not in gfn2
   ! case('kqd');   if (getValue(env,val,ddum)) kqat(2,iz)  = ddum  ! not in gfn2
   case('ele_id');  if (getValue(env,val,idum)) ElemIdPerAtom(iz) = idum 
   end select
end subroutine gfn_elempar_per_atom

subroutine checkElemIdPerAtomMatch(env, ElemIdPerAtom, mol)
  implicit none
   type(TEnvironment), intent(inout) :: env
  integer, intent(in) :: ElemIdPerAtom(:)
  type(TMolecule), intent(in) :: mol
  integer :: i

  ! Check if the lengths match
  if (size(ElemIdPerAtom) /= size(mol%at)) then
    call env%error('Error: Length of ElemIdPerAtom does not match length of mol%at')
    return
  end if

  ! Check if each element matches
  do i = 1, size(ElemIdPerAtom)
    if (ElemIdPerAtom(i) /= mol%at(i)) then
      call env%error("Error: Element mismatch at index ?")
      return
    end if
  end do

  print *, 'Success: ElemIdPerAtom matches mol%at'
end subroutine checkElemIdPerAtomMatch


end subroutine read2Param







! Yufan: below is old version

subroutine readParam &
      (env, iunit,globpar,xtbData,initialize)
   use xtb_mctc_accuracy, only : wp

   use xtb_readin, only : getline => strip_line
   use xtb_type_environment, only : TEnvironment
   use xtb_type_param, only : TxTBParameter
   use xtb_param_paulingen, only : paulingEN
   use xtb_param_atomicrad, only : atomicRad
   use xtb_mctc_param, only: chemical_hardness

   implicit none

   logical, parameter :: debug = .false.

   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: iunit
   type(TxTBParameter), intent(inout) :: globpar
   type(TxTBData), intent(out) :: xtbData
   logical, intent(in) :: initialize

   character(len=1), parameter :: equal = '='
   character(len=1), parameter :: space = ' '
   character(len=1), parameter :: flag = '$'
   character(len=*), parameter :: flag_end = flag//'end'
   integer, parameter :: p_str_length = 48
   integer, parameter :: p_arg_length = 24

   integer, parameter :: max_elem = 118
   integer :: nShell(max_elem)
   integer :: principalQuantumNumber(10,max_elem)
   integer :: angShell(10,max_elem)
   real(wp) :: shellPoly(1:4,max_elem)
   real(wp) :: selfEnergy(10,max_elem)
   real(wp) :: slaterExponent(10,max_elem)
   real(wp) :: thirdOrderAtom(max_elem)
   real(wp) :: atomicHardness(max_elem)
   real(wp) :: shellHardness(1:4,max_elem)
   real(wp) :: electronegativity(max_elem)
   real(wp) :: repAlpha(max_elem)
   real(wp) :: repZeff(max_elem)
   real(wp) :: halogenBond(max_elem)
   real(wp) :: dipKernel(max_elem)
   real(wp) :: quadKernel(max_elem)
   real(wp) :: eeqkcn(max_elem)
   real(wp) :: chargeWidth(max_elem)
   real(wp) :: kqat2(max_elem)
   real(wp) :: kqat(0:2,max_elem)
   real(wp) :: kpair(max_elem,max_elem)
   real(wp) :: kcnat(0:2,max_elem)
   real(wp) :: eeqEN(max_elem)
   real(wp) :: kExpLight, kExp
   character(len=30) :: timestp(max_elem)
   type(dftd_parameter) :: disp

   character(len=:), allocatable :: line

   integer :: mShell, iSh, jSh
   integer :: level
   integer :: err
   logical :: newFormat

   disp = dftd_parameter(s6=1.0_wp, s8=0.0_wp, a1=0.0_wp, a2=0.0_wp, s9=0.0_wp)
   globpar = TxTBParameter()

   principalQuantumNumber=0
   angShell  =0
   nShell  =0
   selfEnergy=0.0_wp
   slaterExponent=0.0_wp
   shellPoly =0.0_wp
   halogenBond = 0.0_wp
   atomicHardness = chemical_hardness(:max_elem)
   shellHardness  =0.0_wp
   thirdOrderAtom  =0.0_wp
   kqat2 =0.0_wp
   eeqkcn=0.0_wp
   eeqen =0.0_wp
   kcnat =0.0_wp

   electronegativity = paulingEN(:max_elem)
   repAlpha = 0.0_wp
   repZeff = 0.0_wp

   dipKernel = 0.0_wp ! read values are scaled by 0.01
   quadKernel = 0.0_wp !  "     "     "    "    "   "

   kpair =1.0_wp

   kExpLight = 0.0_wp
   level = -1
   newFormat = .false.
   call getline(iunit,line,err)
   if (debug) print'(">",a)',line
   readgroups: do
      if (index(line,flag).eq.1) then
         select case(line(2:))
         case('info')
            newFormat = .true.
            call read_info
         case('globpar')
            call read_globpar
         case('pairpar')
            call read_pairpar
         case default
            if (index(line,'Z').eq.2) then
               call read_elempar
            else
               call getline(iunit,line,err)
               if (debug) print'(">",a)',line
            endif
         end select
      else
         call getline(iunit,line,err)
         if (debug) print'(">",a)',line
      endif
      if (err.ne.0) exit readgroups
      !if (index(line,flag_end).gt.0) exit readgroups
   enddo readgroups

   if (.not.newFormat) then
      call env%error("Old format parameter file is not supported anymore")
   end if

   call setpair(level, kpair)

   mShell = maxval(nShell)
   xtbData%level = level
   xtbData%nShell = nShell
   xtbData%ipeashift = globpar%ipeashift * 0.1_wp

   ! Repulsion
   call init(xtbData%repulsion, kExp, kExpLight, 1.0_wp, globpar%renscale, &
      & repAlpha, repZeff, electronegativity)

   ! Coulomb
   xtbData%coulomb%gExp = globpar%alphaj
   xtbData%coulomb%chemicalHardness = atomicHardness(:max_elem)
   allocate(xtbData%coulomb%shellHardness(mShell, max_elem))
   call setGFN1ShellHardness(xtbData%coulomb%shellHardness, nShell, angShell, &
      & atomicHardness, shellHardness)
   xtbData%coulomb%thirdOrderAtom = thirdOrderAtom(:max_elem)
   xtbData%coulomb%electronegativity = eeqEN(:max_elem)
   xtbData%coulomb%kCN = eeqkCN(:max_elem)
   xtbData%coulomb%chargeWidth = chargeWidth(:max_elem)

   ! Dispersion
   xtbData%dispersion%dpar = disp
   xtbData%dispersion%g_a = 3.0_wp
   xtbData%dispersion%g_c = 2.0_wp
   xtbData%dispersion%wf  = 6.0_wp

   ! Hamiltonian
   xtbData%hamiltonian%angShell = angShell(:mShell, :)

   do iSh = 0, 3
      do jSh = 0, 3
         xtbData%hamiltonian%kScale(jSh, iSh) = 0.5_wp * (globpar%kShell(iSh) &
            & + globpar%kShell(jSh))
      end do
   end do
   if (globpar%ksp > 0.0_wp) then
      xtbData%hamiltonian%kScale(0,1) = globpar%ksp
      xtbData%hamiltonian%kScale(1,0) = globpar%ksp
   end if
   if (globpar%ksd > 0.0_wp) then
      xtbData%hamiltonian%kScale(0,2) = globpar%ksd
      xtbData%hamiltonian%kScale(2,0) = globpar%ksd
   end if
   if (globpar%kpd > 0.0_wp) then
      xtbData%hamiltonian%kScale(1,2) = globpar%kpd
      xtbData%hamiltonian%kScale(2,1) = globpar%kpd
   end if
   xtbData%hamiltonian%kDiff = globpar%kDiff

   do iSh = 0, 3
      do jSh = 0, 3
         xtbData%hamiltonian%enScale(jSh, iSh) = 0.005_wp * (globpar%enShell(iSh) &
            & + globpar%enShell(jSh))
      end do
   end do
   xtbData%hamiltonian%enScale4 = globpar%enscale4

   xtbData%hamiltonian%electronegativity = electronegativity(:)
   xtbData%hamiltonian%atomicRad = atomicRad(:)
   xtbData%hamiltonian%shellPoly = shellPoly(:, :)
   xtbData%hamiltonian%pairParam = kpair(:, :)
   xtbData%hamiltonian%selfEnergy = selfEnergy(:mShell, :)
   xtbData%hamiltonian%slaterExponent = slaterExponent(:mShell, :)
   xtbData%hamiltonian%principalQuantumNumber = principalQuantumNumber(:mShell, :)

   allocate(xtbData%hamiltonian%valenceShell(mShell, max_elem))
   call generateValenceShellData(xtbData%hamiltonian%valenceShell, &
      & xtbData%nShell, xtbData%hamiltonian%angShell)

   select case(level)
   case(0)
      ! Hamiltonian
      xtbData%hamiltonian%wExp = 1.0_wp

      allocate(xtbData%hamiltonian%kCN(mShell, max_elem))
      call angToShellData(xtbData%hamiltonian%kCN, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, kcnat)

      allocate(xtbData%hamiltonian%kQShell(mShell, max_elem))
      call angToShellData(xtbData%hamiltonian%kQShell, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, kqat)
      xtbData%hamiltonian%kQAtom = kqat2(:)

      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN2ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN0NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%srb)
      xtbData%srb%shift = globpar%srbshift
      xtbData%srb%prefactor = globpar%srbpre
      xtbData%srb%steepness = globpar%srbexp
      xtbData%srb%enScale = globpar%srbken

      ! Dispersion
      call newD4Model(xtbData%dispersion%dispm, xtbData%dispersion%g_a, &
         & xtbData%dispersion%g_c, p_refq_goedecker)

   case(1)
      ! Halogen
      allocate(xtbData%halogen)
      call init(xtbData%halogen, globpar%xbrad, globpar%xbdamp, halogenBond)

      ! Hamiltonian
      xtbData%hamiltonian%wExp = 0.0_wp

      allocate(xtbData%hamiltonian%kCN(mShell, max_elem))
      call setGFN1kCN(xtbData%hamiltonian%kCN, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, xtbData%hamiltonian%selfEnergy, &
         & globpar%cnshell)

      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN1ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN1NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

      ! Dispersion
      if (.not.allocated(reference_c6)) call copy_c6(reference_c6)

   case(2)
      ! Coulomb
      if (any(globpar%gam3shell > 0.0_wp)) then
         allocate(xtbData%Coulomb%thirdOrderShell(mShell, max_elem))
         call setGFN2ThirdOrderShell(xtbData%Coulomb%thirdOrderShell, &
            & xtbData%nShell, xtbData%hamiltonian%angShell, thirdOrderAtom, &
            & globpar%gam3shell)
      end if

      ! Multipole
      allocate(xtbData%multipole)
      call init(xtbData%multipole, globpar%aesshift, globpar%aesexp, &
         & globpar%aesrmax, globpar%aesdmp3, globpar%aesdmp5, &
         & dipKernel, quadKernel)

      ! Hamiltonian
      xtbData%hamiltonian%wExp = 0.5_wp

      allocate(xtbData%hamiltonian%kCN(mShell, max_elem))
      call angToShellData(xtbData%hamiltonian%kCN, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, kcnat)

      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN2ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN2NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

      ! Dispersion
      call newD4Model(xtbData%dispersion%dispm, xtbData%dispersion%g_a, &
         & xtbData%dispersion%g_c, p_refq_gfn2xtb)

   end select

contains

subroutine read_info
   use xtb_readin, only : getValue
   character(len=:), allocatable :: key, val
   integer :: ie, idum
   do
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) exit
      if (index(line,flag).gt.0) exit

      ie = index(line,space)
      if (line.eq.'') cycle ! skip empty lines
      if (ie.eq.0) cycle

      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))

      select case(key)
      case default
         call env%warning("Unknown key '"//key//"' for '"//flag//"info'")
      case('level')
         if (getValue(env,val,idum)) level = idum
      case('name')
         xtbData%name = val
      case('doi')
         xtbData%doi = val
      end select
   enddo
end subroutine read_info

subroutine read_globpar
   implicit none
   character(len=:), allocatable :: key, val
   integer :: ie
   do
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) exit
      if (index(line,flag).gt.0) exit

      ie = index(line,space)
      if (line.eq.'') cycle ! skip empty lines
      if (ie.eq.0) cycle

      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))

      call gfn_globpar(key,val,globpar)

   enddo
end subroutine read_globpar

subroutine gfn_globpar(key,val,globpar)
   use xtb_readin, only : getValue
   implicit none
   character(len=*), intent(in) :: key, val
   type(TxTBParameter), intent(inout) :: globpar
   real(wp) :: ddum
   select case(key)
   case default
      call env%warning("Unknown key '"//key//"' for '"//flag//"globpar'")
   case('kexp'); if (getValue(env,val,ddum)) kExp = ddum
   case('kexplight'); if (getValue(env,val,ddum)) kExpLight = ddum
   case('ks'); if (getValue(env,val,ddum)) globpar%kshell(0) = ddum
   case('kp'); if (getValue(env,val,ddum)) globpar%kshell(1) = ddum
   case('kd'); if (getValue(env,val,ddum)) globpar%kshell(2) = ddum
   case('kf'); if (getValue(env,val,ddum)) globpar%kshell(3) = ddum
   case('ksp'); if (getValue(env,val,ddum)) globpar%ksp = ddum
   case('ksd'); if (getValue(env,val,ddum)) globpar%ksd = ddum
   case('kpd'); if (getValue(env,val,ddum)) globpar%kpd = ddum
   case('kdiff'); if (getValue(env,val,ddum)) globpar%kdiff = ddum
   case('kdiffa'); if (getValue(env,val,ddum)) globpar%kdiffa = ddum
   case('kdiffb'); if (getValue(env,val,ddum)) globpar%kdiffb = ddum
   case('ens'); if (getValue(env,val,ddum)) globpar%enshell(0) = ddum
   case('enp'); if (getValue(env,val,ddum)) globpar%enshell(1) = ddum
   case('end'); if (getValue(env,val,ddum)) globpar%enshell(2) = ddum
   case('enf'); if (getValue(env,val,ddum)) globpar%enshell(3) = ddum
   case('enscale'); if (getValue(env,val,ddum)) globpar%enshell = ddum
   case('enscale4'); if (getValue(env,val,ddum)) globpar%enscale4 = ddum
   case('renscale'); if (getValue(env,val,ddum)) globpar%renscale = ddum
   case('cns'); if (getValue(env,val,ddum)) globpar%cnshell(:, 0) = ddum
   case('cnp'); if (getValue(env,val,ddum)) globpar%cnshell(:, 1) = ddum
   case('cnd'); if (getValue(env,val,ddum)) globpar%cnshell(:, 2) = ddum
   case('cnf'); if (getValue(env,val,ddum)) globpar%cnshell(:, 3) = ddum
   case('cnd1'); if (getValue(env,val,ddum)) globpar%cnshell(1, 2) = ddum
   case('cnd2'); if (getValue(env,val,ddum)) globpar%cnshell(2, 2) = ddum
   case('gam3s'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 0) = ddum
   case('gam3p'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 1) = ddum
   case('gam3d'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 2) = ddum
   case('gam3f'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 3) = ddum
   case('gam3d1'); if (getValue(env,val,ddum)) globpar%gam3shell(1, 2) = ddum
   case('gam3d2'); if (getValue(env,val,ddum)) globpar%gam3shell(2, 2) = ddum
   case('srbshift'); if (getValue(env,val,ddum)) globpar%srbshift = ddum
   case('srbpre'); if (getValue(env,val,ddum)) globpar%srbpre = ddum
   case('srbexp'); if (getValue(env,val,ddum)) globpar%srbexp = ddum
   case('srbken'); if (getValue(env,val,ddum)) globpar%srbken = ddum
   case('wllscal'); if (getValue(env,val,ddum)) globpar%wllscal = ddum
   case('ipeashift'); if (getValue(env,val,ddum)) globpar%ipeashift = ddum
   case('gscal'); if (getValue(env,val,ddum)) globpar%gscal = ddum
   case('zcnf'); if (getValue(env,val,ddum)) globpar%zcnf = ddum
   case('tscal'); if (getValue(env,val,ddum)) globpar%tscal = ddum
   case('kcn'); if (getValue(env,val,ddum)) globpar%kcn = ddum
   case('fpol'); if (getValue(env,val,ddum)) globpar%fpol = ddum
   case('ken'); if (getValue(env,val,ddum)) globpar%ken = ddum
   case('lshift'); if (getValue(env,val,ddum)) globpar%lshift = ddum
   case('lshifta'); if (getValue(env,val,ddum)) globpar%lshifta = ddum
   case('split'); if (getValue(env,val,ddum)) globpar%split = ddum
   case('zqf'); if (getValue(env,val,ddum)) globpar%zqf = ddum
   case('alphaj'); if (getValue(env,val,ddum)) globpar%alphaj = ddum
   case('kexpo'); if (getValue(env,val,ddum)) globpar%kexpo = ddum
   case('dispa'); if (getValue(env,val,ddum)) globpar%dispa = ddum
   case('dispb'); if (getValue(env,val,ddum)) globpar%dispb = ddum
   case('dispc'); if (getValue(env,val,ddum)) globpar%dispc = ddum
   case('dispatm'); if (getValue(env,val,ddum)) globpar%dispatm = ddum
   case('a1'); if (getValue(env,val,ddum)) disp%a1 = ddum
   case('a2'); if (getValue(env,val,ddum)) disp%a2 = ddum
   case('s6'); if (getValue(env,val,ddum)) disp%s6 = ddum
   case('s8'); if (getValue(env,val,ddum)) disp%s8 = ddum
   case('s9'); if (getValue(env,val,ddum)) disp%s9 = ddum
   case('xbdamp'); if (getValue(env,val,ddum)) globpar%xbdamp = ddum
   case('xbrad'); if (getValue(env,val,ddum)) globpar%xbrad = ddum
   case('aesdmp3'); if (getValue(env,val,ddum)) globpar%aesdmp3 = ddum
   case('aesdmp5'); if (getValue(env,val,ddum)) globpar%aesdmp5 = ddum
   case('aesshift'); if (getValue(env,val,ddum)) globpar%aesshift = ddum
   case('aesexp'); if (getValue(env,val,ddum)) globpar%aesexp = ddum
   case('aesrmax'); if (getValue(env,val,ddum)) globpar%aesrmax = ddum
   end select
end subroutine gfn_globpar

subroutine read_pairpar
   use xtb_mctc_strings
   use xtb_readin, only : getValue
   implicit none
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   integer  :: iAt,jAt
   real(wp) :: ddum
   do
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) exit
      if (index(line,flag).gt.0) exit

      call parse(line,space,argv,narg)
      if (narg .ne. 3) cycle
      if (getValue(env,trim(argv(1)),iAt) .and. &
         &getValue(env,trim(argv(2)),jAt) .and. &
         &getValue(env,trim(argv(3)),ddum)) then
            kpair(iAt,jAt) = ddum
            kpair(jAt,iAt) = ddum
      endif
   enddo
end subroutine read_pairpar

subroutine read_elempar
   use xtb_mctc_strings
   use xtb_readin
   implicit none
   character(len=:), allocatable :: key, val
   integer :: iz, ie
   if (getValue(env,line(4:5),iz)) then
      timestp(iz) = line(7:len_trim(line))
      do
         call getline(iunit,line,err)
         if (debug) print'("->",a)',line
         if (err.ne.0) exit
         if (index(line,flag).gt.0) exit

         ie = index(line,equal)
         if (line.eq.'') cycle ! skip empty lines
         if (ie.eq.0) cycle

         key = lowercase(trim(line(:ie-1)))
         val = trim(adjustl(line(ie+1:)))

         call gfn_elempar(key,val,iz)

      enddo
   else
      call getline(iunit,line,err)
   endif
end subroutine read_elempar

subroutine gfn_elempar(key,val,iz)
   use xtb_mctc_strings
   use xtb_readin
   implicit none
   character(len=*), intent(in) :: key, val
   integer, intent(in) :: iz
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   integer :: i, ii
   integer :: idum
   real(wp) :: ddum
   select case(key)
   case default
      call env%warning("Unknown key '"//key//"' for '"//flag//"Z'")
   case('ao')
      !print'(a,":",a)',key,val
      if (mod(len(val),2).eq.0) then
         nShell(iz) = len(val)/2
         do i = 1, nShell(iz)
            ii = 2*i-1
            !print*,i,ii,val(ii:ii),val(ii+1:ii+1)
            if (getValue(env,val(ii:ii),idum)) then
               principalQuantumNumber(i,iz) = idum
               select case(val(ii+1:ii+1))
               case('s'); angShell(i,iz) = 0
               case('p'); angShell(i,iz) = 1
               case('d'); angShell(i,iz) = 2
               case('f'); angShell(i,iz) = 3
               case('g'); angShell(i,iz) = 4
               case('S'); angShell(i,iz) = 0
               end select
            endif
         enddo
      endif
   case('lev')
      call parse(val,space,argv,narg)
      if (narg .eq. nShell(iz)) then
         do i = 1, nShell(iz)
            if (getValue(env,trim(argv(i)),ddum)) selfEnergy(i,iz) = ddum
         enddo
      endif
   case('exp')
      call parse(val,space,argv,narg)
      if (narg .eq. nShell(iz)) then
         do i = 1, nShell(iz)
            if (getValue(env,trim(argv(i)),ddum)) slaterExponent(i,iz) = ddum
         enddo
      endif
   case('en');  if (getValue(env,val,ddum)) electronegativity(iz)    = ddum
   case('gam'); if (getValue(env,val,ddum)) atomicHardness(iz)   = ddum
   case('xi');  if (getValue(env,val,ddum)) eeqEN(iz) = ddum
   case('alpg'); if (getValue(env,val,ddum)) chargeWidth(iz)   = ddum
   case('gam3');  if (getValue(env,val,ddum)) thirdOrderAtom(iz)    = ddum * 0.1_wp
   case('kappa'); if (getValue(env,val,ddum)) eeqkCN(iz)  = ddum
   case('cxb');   if (getValue(env,val,ddum)) halogenBond(iz) = ddum * 0.1_wp
   case('kqat2'); if (getValue(env,val,ddum)) kqat2(iz)   = ddum
   case('dpol');  if (getValue(env,val,ddum)) dipKernel(iz)   = ddum * 0.01_wp
   case('qpol');  if (getValue(env,val,ddum)) quadKernel(iz)   = ddum * 0.01_wp
   case('repa');  if (getValue(env,val,ddum)) repAlpha(iz)   = ddum
   case('repb');  if (getValue(env,val,ddum)) repZeff(iz)   = ddum
   case('polys'); if (getValue(env,val,ddum)) shellPoly(1,iz) = ddum
   case('polyp'); if (getValue(env,val,ddum)) shellPoly(2,iz) = ddum
   case('polyd'); if (getValue(env,val,ddum)) shellPoly(3,iz) = ddum
   case('polyf'); if (getValue(env,val,ddum)) shellPoly(4,iz) = ddum
   case('lpars'); if (getValue(env,val,ddum)) shellHardness(1,iz)  = ddum * 0.1_wp
   case('lparp'); if (getValue(env,val,ddum)) shellHardness(2,iz)  = ddum * 0.1_wp
   case('lpard'); if (getValue(env,val,ddum)) shellHardness(3,iz)  = ddum * 0.1_wp
   case('lparf'); if (getValue(env,val,ddum)) shellHardness(4,iz)  = ddum * 0.1_wp
   case('kcns');  if (getValue(env,val,ddum)) kcnat(0,iz) = ddum * 0.1_wp
   case('kcnp');  if (getValue(env,val,ddum)) kcnat(1,iz) = ddum * 0.1_wp
   case('kcnd');  if (getValue(env,val,ddum)) kcnat(2,iz) = ddum * 0.1_wp
   case('kqs');   if (getValue(env,val,ddum)) kqat(0,iz)  = ddum  
   case('kqp');   if (getValue(env,val,ddum)) kqat(1,iz)  = ddum
   case('kqd');   if (getValue(env,val,ddum)) kqat(2,iz)  = ddum
   end select
end subroutine gfn_elempar

end subroutine readParam



end module xtb_readparam


      logical function maingroup(i)
      integer i
      logical main_group(107)
      data main_group /&
     &  2*.true.,&                              ! H  - He
     &  8*.true.,&                              ! Li - Ne
     &  8*.true.,&                              ! Na - Ar
     &  2*.true., 9*.false., 7*.true.,&         ! K  - Kr
     &  2*.true., 9*.false., 7*.true.,&         ! Rb - Xe
     &  2*.true.,23*.false., 7*.true.,&         ! Cs - Rn
     & 21*.true.                               / ! Fr - Tv

      maingroup = main_group(i)

      end function maingroup
