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
  use types, only: column
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
  integer :: nexp !< number of exponents
  integer, allocatable :: nval(:) !< n-value associated to the exponent
  real*8, allocatable :: eexp(:) !< exponent value

  ! fitting set, and subsets
  integer :: nfit !< number of systems in the fitting set
  integer :: nset !< number of subsets
  character*255, allocatable :: iset_label(:) !< labels for the subsets
  integer, allocatable :: iset_ini(:), iset_n(:), iset_step(:) !< initial, number, and step for the subsets

  ! named files
  character(len=:), allocatable :: emptyfile !< name of the empty file
  character(len=:), allocatable :: reffile !< name of the reference energy file
  character(len=:), allocatable :: wfile !< name of the weight file
  character(len=:), allocatable :: namesfile !< name of the names file
  character*255, allocatable :: efile(:,:,:) !< name of the energy terms files
  character*255, allocatable :: subfile(:) !< name of the subtract energy file
  integer :: nsubfiles

  ! labels for the angular momentum channels
  character*1, parameter :: lname(6) = (/"l","s","p","d","f","g"/) !< names of the angmom channels

  ! column information
  integer :: ncols
  type(column), allocatable :: col(:)

  ! data for the fits
  real*8, allocatable :: yempty(:), yref(:), ydisp(:), w(:) !< empty, ref, disp, and weight data
  real*8, allocatable :: x(:,:) !< matrix containing the BSIP term energies
  real*8, allocatable :: ytarget(:) !< target of the fit
  character*128, allocatable :: names(:) !< names of the systems in the fitting set

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
    nfit = 0
    nset = 0
    if (allocated(atom)) deallocate(atom)
    if (allocated(lmax)) deallocate(lmax)
    if (allocated(nval)) deallocate(nval)
    if (allocated(eexp)) deallocate(eexp)
    if (allocated(subfile)) deallocate(subfile)
    if (allocated(iset_label)) deallocate(iset_label)
    if (allocated(iset_ini)) deallocate(iset_ini)
    if (allocated(iset_n)) deallocate(iset_n)
    if (allocated(iset_step)) deallocate(iset_step)
    allocate(atom(1))
    allocate(lmax(1))
    allocate(nval(1))
    allocate(eexp(1))
    allocate(subfile(1))
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

  !> Write information about the output data to the output
  subroutine global_printinfo()
    use tools_io, only: string, ioj_left, uout

    integer :: i, j, k

    write (uout,'("+ Summary of input data")')
    write (uout,'("Atomic informaton: ")')
    write (uout,'("# Id At lmax")')
    do i = 1, natoms
       write (uout,'(2X,99(A,X))') string(i,2,ioj_left), string(atom(i),3), string(lname(lmax(i)),2)
    end do
    write (uout,'("List of exponents: ")')
    write (uout,'("# Id n      exponent ")')
    do i = 1, nexp
       write (uout,'(2X,99(A,X))') string(i,2,ioj_left), string(nval(i),3), string(eexp(i),'f',12,8,4)
    end do
    write (uout,'("List of evaluation sets: ")')
    write (uout,'("#Id  ini  num  step name")')
    do i = 1, nset
       write (uout,'(2X,99(A,X))') string(i,2,ioj_left), string(iset_ini(i),5), string(iset_n(i),4), &
          string(iset_step(i),3), string(iset_label(i))
    end do
    write (uout,'("Size of the fitting set: ",A)') string(nfit)
    write (uout,'("Data path: ",A)') string(datapath)
    write (uout,'("Names file: ",A)') string(namesfile)
    write (uout,'("Weight file: ",A)') string(wfile)
    write (uout,'("Empty file: ",A)') string(emptyfile)
    write (uout,'("Reference file: ",A)') string(reffile)
    if (nsubfiles > 0) then
       write (uout,'("List of subtraction files: ")')
       do i = 1, nsubfiles
          write (uout,'(2X,A,": ",A)') string(i), string(subfile(i))
       end do
    end if
    write (uout,'("List of energy term files: ")')
    write (uout,'("#At ang iexp n   exponent  filename")')
    do i = 1, natoms
       do j = 1, lmax(i)
          do k = 1, nexp
             write (uout,'(2X,99(A,X))') string(atom(i),2), string(lname(j),2), string(k,4), &
                string(nval(k),2), string(eexp(k),'f',10,4,4), string(efile(i,j,k))
          end do
       end do
    end do
    write (uout,*)

  end subroutine global_printinfo

  !> Fill the information in the column variable and ncols based on
  !> the input data.
  subroutine global_fillcol()

    integer :: i, j, k, n

    ! number of columns
    ncols = 0
    do i = 1, natoms
       do j = 1, lmax(i)
          do k = 1, nexp
             ncols = ncols + 1
          end do
       end do
    end do

    ! allocate the column variable
    if (allocated(col)) deallocate(col)
    allocate(col(ncols))

    ! assign the column information
    n = 0
    do i = 1, natoms
       do j = 1, lmax(i)
          do k = 1, nexp
             n = n + 1
             col(n)%atom = atom(i)
             col(n)%iatom = i
             col(n)%l = j
             col(n)%iexp = k
             col(n)%n = nval(k)
             col(n)%eexp = eexp(k)
          end do
       end do
    end do

  end subroutine global_fillcol

end module global
