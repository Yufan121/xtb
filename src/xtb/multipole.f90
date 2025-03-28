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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> TODO
module xtb_xtb_multipole
   use xtb_mctc_accuracy, only : wp
   use xtb_xtb_data, only : TMultipoleData
   implicit none
   private

   public :: TxTBMultipole, init


   type, extends(TMultipoleData) :: TxTBMultipole

      real(wp), allocatable :: gab3(:, :)

      real(wp), allocatable :: gab5(:, :)

      type(TMultipoleData), allocatable :: multipoledataPerAtom

   end type


   interface init
      module procedure :: initMultipole
      module procedure :: initMultipolePerAtom
   end interface init


contains


subroutine initMultipole(self, input)

   type(TxTBMultipole), intent(out) :: self

   type(TMultipoleData), intent(in) :: input

   self%TMultipoleData = input


end subroutine initMultipole

subroutine initMultipolePerAtom(self, input, inputPerAtom)

   type(TxTBMultipole), intent(out) :: self

   type(TMultipoleData), intent(in) :: input, inputPerAtom

   self%TMultipoleData = input

   allocate(self%multipoledataPerAtom)
   self%multipoledataPerAtom = inputPerAtom

end subroutine initMultipolePerAtom



end module xtb_xtb_multipole
