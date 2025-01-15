! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb. If not, see <https://www.gnu.org/licenses/>.

module xtb_xtb_gfn2_per_atom
  use xtb_mctc_accuracy, only: wp
  ! use xtb_xtb_gfn2, only: gfn2Kinds
  use xtb_xtb_data
  implicit none
  public :: setGFN1ShellHardnessPerAtom, setGFN2ThirdOrderShellPerAtom, setGFN2ReferenceOccPerAtom, setGFN2NumberOfPrimitivesPerAtom


  integer, parameter :: maxElem = 86

  integer, parameter :: gfn2Kinds(118) = [&
  &  1,                                                 1, &! H-He
  &  1, 1,                               1, 1, 1, 1, 1, 1, &! Li-Ne
  &  1, 1,                               1, 1, 1, 1, 1, 1, &! Na-Ar
  &  1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, &! K-Kr
  &  1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, &! Rb-Xe
  &  1, 1, &! Cs/Ba
  &        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &!La-Lu
  &        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, &! Lu-Rn
  &  1, 1, &
  &        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &!Fr-
  &        2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1 ]! -Og

  !> Reference occupation of the atom
  real(wp), parameter :: referenceOcc(0:2, 1:maxElem) = reshape([&
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp,  1.0_wp, 3.0_wp, 0.0_wp, &
      & 1.5_wp, 3.5_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp, &
      & 2.0_wp, 6.0_wp, 0.0_wp,  1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  1.5_wp, 2.5_wp, 0.0_wp,  1.5_wp, 3.5_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp,  1.0_wp, 1.0_wp, 3.0_wp,  1.0_wp, 1.0_wp, 4.0_wp, &
      & 1.0_wp, 1.0_wp, 5.0_wp,  1.0_wp, 1.0_wp, 6.0_wp,  1.0_wp, 1.0_wp, 7.0_wp, &
      & 1.0_wp, 1.0_wp, 8.0_wp,  1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  1.5_wp, 2.5_wp, 0.0_wp,  1.5_wp, 3.5_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp,  1.0_wp, 1.0_wp, 3.0_wp,  1.0_wp, 1.0_wp, 4.0_wp, &
      & 1.0_wp, 1.0_wp, 5.0_wp,  1.0_wp, 1.0_wp, 6.0_wp,  1.0_wp, 1.0_wp, 7.0_wp, &
      & 1.0_wp, 1.0_wp, 8.0_wp,  1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 1.0_wp, 0.0_wp,  2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp, &
      & 2.0_wp, 4.0_wp, 0.0_wp,  2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 1.0_wp,  1.0_wp, 1.0_wp, 2.0_wp, &
      & 1.0_wp, 1.0_wp, 3.0_wp,  1.0_wp, 1.0_wp, 4.0_wp,  1.0_wp, 1.0_wp, 5.0_wp, &
      & 1.0_wp, 1.0_wp, 6.0_wp,  1.0_wp, 1.0_wp, 7.0_wp,  1.0_wp, 1.0_wp, 8.0_wp, &
      & 1.0_wp, 0.0_wp,10.0_wp,  2.0_wp, 0.0_wp, 0.0_wp,  2.0_wp, 1.0_wp, 0.0_wp, &
      & 2.0_wp, 2.0_wp, 0.0_wp,  2.0_wp, 3.0_wp, 0.0_wp,  2.0_wp, 4.0_wp, 0.0_wp, &
      & 2.0_wp, 5.0_wp, 0.0_wp,  2.0_wp, 6.0_wp, 0.0_wp], shape(referenceOcc))


  ! Third order Hubbard derivatives
  real(wp), allocatable :: thirdOrderAtomPerAtom(:)

  ! Exponents of repulsion term
  real(wp), allocatable :: repAlphaPerAtom(:)

  ! Dipole exchange-correlation kernel
  real(wp), allocatable :: dipKernelPerAtom(:)

  ! Quadrupole exchange-correlation kernel
  real(wp), allocatable :: quadKernelPerAtom(:)

  !> Effective nuclear charge
  real(wp), allocatable :: repZeffPerAtom(:)
  
  ! Number of shells (for kcn)
  integer, allocatable :: nShellPerAtom(:)

  ! Element Species
  integer, allocatable :: ElemIdPerAtom(:)

  ! electronegativityPerAtom
  real(wp), allocatable :: electronegativityPerAtom(:)

  ! valanceCNPerAtom
  real(wp), allocatable :: valanceCNPerAtom(:)

  ! multiRadPerAtom
  real(wp), allocatable :: multiRadPerAtom(:)


  ! atomicHardnessPerAtom
  real(wp), allocatable :: atomicHardnessPerAtom(:)

  ! Shell polynomials to scale Hamiltonian elements, (4, natom)
  real(wp), allocatable :: shellPolyPerAtom(:, :)

  !> Scaling factors for shell electrostatics (3, natom)
  real(wp), allocatable :: shellHardnessPerAtom(:, :)

  ! Atomic level information, (3, natom)
  real(wp), allocatable :: selfEnergyPerAtom(:, :)

  ! Exponent of the Slater function, (3, natom)
  real(wp), allocatable :: slaterExponentPerAtom(:, :)

  ! Angular momentum of each shell, (3, natom)
  integer, allocatable :: angShellPerAtom(:, :)

  ! Coordination number dependence of the atomic levels, (4, natom)
  real(wp), allocatable :: kCNPerAtom(:, :) ! all zero at (4, :)

  !> kcnat, (3, natom)
  real(wp), allocatable :: kcnatPerAtom(:, :)

  ! Principal quantum number of each shell, (3, natom)
  integer, allocatable :: principalQuantumNumberPerAtom(:, :)



  contains

  subroutine setGFN1ShellHardnessPerAtom(shellHardness, nShell, angShell, atomicHardness, atomicHardnessPerElem, &
      & angHardness, angHardnessPerElem, ElemIdPerAtom)

  ! Yyufan: modification 1. change index end to eElem 


   real(wp), intent(out) :: shellHardness(:, :)

   integer, intent(in) :: nShell(:)

   integer, intent(in) :: angShell(:, :)

   real(wp), intent(in) :: atomicHardness(:), atomicHardnessPerElem(:) ! access by elem_id

   real(wp), intent(in) :: angHardness(0:, :), angHardnessPerElem(0:, :) ! access by elem_id

   integer, intent(in) :: ElemIdPerAtom(:)

   integer :: nAtom, iAt, iSh, lAng

   nAtom = min(size(shellHardness, dim=2), size(nShell), size(angShell, dim=2))

   shellHardness(:, :) = 0.0_wp
   do iAt = 1, nAtom    ! natom
      do iSh = 1, nShell(iAt)   ! nShell
         lAng = angShell(iSh, iAt)    
         shellHardness(iSh, iAt) = (atomicHardness(iAt) + atomicHardnessPerElem(ElemIdPerAtom(iAt))) * &
            (1.0_wp + angHardness(lAng, iAt) + angHardnessPerElem(lAng, ElemIdPerAtom(iAt)))
      end do
   end do

  end subroutine setGFN1ShellHardnessPerAtom


  subroutine setGFN2ThirdOrderShellPerAtom(thirdOrderShell, nShell, angShell, &
      & thirdOrderAtom, gam3Shell, ElemIdPerAtom)

  ! Yyufan: modification 1. change index end to eElem 2. add a mapping from atom id to kind

   real(wp), intent(out) :: thirdOrderShell(:, :)

   integer, intent(in) :: nShell(:)

    integer, intent(in) :: ElemIdPerAtom(:)

   integer, intent(in) :: angShell(:, :)

   real(wp), intent(in) :: thirdOrderAtom(:)

   real(wp), intent(in) :: gam3Shell(:, 0:)

   integer :: nElem, iZp, iSh, lAng, iKind

   nElem = min(size(thirdOrderShell, dim=2), size(nShell), size(angShell, dim=2), &
      & size(thirdOrderAtom))

   thirdOrderShell(:, :) = 0.0_wp
   do iZp = 1, nElem
      iKind = gfn2Kinds(ElemIdPerAtom(iZp)) ! atomic number
      do iSh = 1, nShell(iZp)
         lAng = angShell(iSh, iZp)
         thirdOrderShell(iSh, iZp) = thirdOrderAtom(iZp) * gam3Shell(iKind, lAng)
      end do
   end do

  end subroutine setGFN2ThirdOrderShellPerAtom



  subroutine setGFN2ReferenceOccPerAtom(self, nShell, ElemIdPerAtom)

    !> Data instance
    type(THamiltonianData), intent(inout) :: self

    !> Number of shells
    integer, intent(in) :: nShell(:)
    integer, intent(in) :: ElemIdPerAtom(:)


    integer :: lAng, iZp, iSh, nAtom
    logical :: valShell(0:3)

    nAtom = size(nShell)

    self%referenceOcc(:, :) = 0.0_wp
    do iZp = 1, nAtom
        do iSh = 1, nShell(iZp)
          lAng = self%angShell(iSh, iZp)
          if (self%valenceShell(iSh, iZp) /= 0) then
              self%referenceOcc(iSh, iZp) = referenceOcc(lAng, ElemIdPerAtom(iZp))
          end if
        end do
    end do

  end subroutine setGFN2ReferenceOccPerAtom

  subroutine setGFN2NumberOfPrimitivesPerAtom(self, nShell)

    !> Data instance
    type(THamiltonianData), intent(inout) :: self

    !> Number of shells
    integer, intent(in) :: nShell(:)

    integer :: nPrim, iZp, iSh, nAtom

    nAtom = size(nShell)


    do iZp = 1, nAtom
        do iSh = 1, nShell(iZp)
          nPrim = 0
          if (iZp <= 2) then
              select case(self%angShell(iSh, iZp))
              case(0)
                nPrim = 3
              case(1)
                nPrim = 4
              end select
          else
              select case(self%angShell(iSh, iZp))
              case(0)
                if (self%principalQuantumNumber(iSh, iZp) > 5) then
                    nPrim = 6
                else
                    nPrim = 4
                end if
              case(1)
                if (self%principalQuantumNumber(iSh, iZp) > 5) then
                    nPrim = 6
                else
                    nPrim = 4
                end if
              case(2)
                nPrim = 3
              case(3)
                nPrim = 4
              end select
          end if
          self%numberOfPrimitives(iSh, iZp) = nPrim
        end do
    end do

  end subroutine setGFN2NumberOfPrimitivesPerAtom


end module xtb_xtb_gfn2_per_atom


