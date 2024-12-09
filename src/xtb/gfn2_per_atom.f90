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
  implicit none
  public :: thirdOrderAtomPerAtom

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


  

end module xtb_xtb_gfn2_per_atom


