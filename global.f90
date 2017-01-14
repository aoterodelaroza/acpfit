! Copyright (c) 2017 Alberto Otero de la Roza <aoterodelaroza@gmail.com>
!
! acpfit is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at
! your option) any later version.
!
! acpfit is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see
! <http://www.gnu.org/licenses/>.
module global
  implicit none

  public

  ! file prefix
  character(len=:), allocatable :: fileroot !< file prefix

  ! global input data
  character(len=:), allocatable :: datapath !< path to the data files
  integer :: natoms !< number of atoms
  character*2, allocatable :: atom(:) !< atomic symbols
  integer, allocatable :: lmax(:) !< maximum angular momentum for each atom
  integer :: maxlmax !< maximum lmax

  ! list of exponents and n-values
  integer :: nexp ! number of exponents
  integer, allocatable :: nval(:) ! n-value associated to the exponent
  real*8, allocatable :: eexp(:) ! exponent value

  ! fitting set
  integer :: nfit
  integer :: nset
  character*255, allocatable :: iset_label(:)
  integer, allocatable :: iset_ini(:), iset_n(:), iset_step(:)

  ! named files
  character(len=:), allocatable :: emptyfile !< name of the empty file
  character(len=:), allocatable :: reffile !< name of the reference energy file
  character(len=:), allocatable :: wfile !< name of the weight file
  character(len=:), allocatable :: namesfile !< name of the names file
  character*255, allocatable :: efile(:,:,:) !< name of the energy terms files
  character*255, allocatable :: efilei(:) !< name of the energy terms files in manual input
  character*255, allocatable :: subfile(:) !< name of the subtract energy file
  integer :: nsubfiles
  integer :: nefilesi
  

  ! labels for the angular momentum channels
  character*1, parameter :: lname(6) = (/"l","s","p","d","f","g"/)

contains
  
  !> Initialize all global variables
  subroutine global_init()

    natoms = 0
    datapath = ""
    emptyfile = "empty.dat"
    reffile = "ref.dat"
    wfile = "w.dat"
    namesfile = "names.dat"
    nexp = 0
    nsubfiles = 0
    nefilesi = 0
    nfit = 0
    nset = 0
    if (allocated(atom)) deallocate(atom)
    if (allocated(lmax)) deallocate(lmax)
    if (allocated(nval)) deallocate(nval)
    if (allocated(eexp)) deallocate(eexp)
    if (allocated(subfile)) deallocate(subfile)
    if (allocated(efilei)) deallocate(efilei)
    if (allocated(iset_label)) deallocate(iset_label)
    if (allocated(iset_ini)) deallocate(iset_ini)
    if (allocated(iset_n)) deallocate(iset_n)
    if (allocated(iset_step)) deallocate(iset_step)
    allocate(atom(1))
    allocate(lmax(1))
    allocate(nval(1))
    allocate(eexp(1))
    allocate(subfile(1))
    allocate(efilei(1))
    allocate(iset_label(1))
    allocate(iset_ini(1))
    allocate(iset_n(1))
    allocate(iset_step(1))

  end subroutine global_init

  !> Check the sanity of the global variables
  subroutine global_check(natoms_ang)
    use tools_io, only: ferror, faterr
    integer :: natoms_ang

    ! check and output data path
    if (len_trim(datapath) < 1) &
       call ferror("acpfit","no path to the data; use DATAPATH",faterr)

    ! check atoms
    if (natoms <= 0) call ferror("acpfit","no atoms in input; use ATOM",faterr)

    ! check angular momentum lmax
    if (natoms_ang == 0) & 
       call ferror("acpfit","no maximum angular momentum; use LMAX",faterr)
    if (natoms /= natoms_ang) & 
       call ferror("acpfit","number of atoms in ATOM and LMAX inconsistent",faterr)
    maxlmax = maxval(lmax)

    ! check exponents
    if (nexp == 0) &
       call ferror("acpfit","no exponent information; use EXP",faterr)
    if (any(nval < 0) .or. any(nval > 2)) &
       call ferror("acpfit","n values can only be 0, 1, or 2",faterr)

    ! check number of molecules
    if (nfit == 0) &
       call ferror("acpfit","no number of molecules in fitting set; use NFIT",faterr)

    

  end subroutine global_check

end module global
