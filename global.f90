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
  integer :: maxnamelen !< maximum length of the names

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
    use tools_io, only: ferror, faterr, string
    integer, intent(in) :: natoms_ang

    integer :: i

    ! check and output data path
    if (len_trim(datapath) < 1) &
       call ferror("global_check","no path to the data; use DATAPATH",faterr)

    ! check atoms
    if (natoms <= 0) call ferror("global_check","no atoms in input; use ATOM",faterr)

    ! check angular momentum lmax
    if (natoms_ang == 0) & 
       call ferror("global_check","no maximum angular momentum; use LMAX",faterr)
    if (natoms /= natoms_ang) & 
       call ferror("global_check","number of atoms in ATOM and LMAX inconsistent",faterr)
    maxlmax = maxval(lmax)

    ! check exponents
    if (nexp == 0) &
       call ferror("global_check","no exponent information; use EXP",faterr)
    if (any(nval < 0) .or. any(nval > 2)) &
       call ferror("global_check","n values can only be 0, 1, or 2",faterr)

    ! check number of molecules
    if (nfit == 0) &
       call ferror("global_check","no number of molecules in fitting set; use NFIT",faterr)

    ! check that the sets are sane
    do i = 1, nset
       if (iset_ini(i) < 1 .or. iset_ini(i) > nfit) &
          call ferror("global_check","erroneous initial system in set " // string(iset_label(i)),faterr)
       if (iset_ini(i)+iset_n(i)*iset_step(i)-1 > nfit) &
          call ferror("global_check","erroneous final system in set " // string(iset_label(i)),faterr)
    end do

  end subroutine global_check

  !> Write information about the output data to the output
  subroutine global_printinfo()
    use tools_io, only: string, ioj_left, uout

    integer :: i, j, k

    write (uout,'("+ Summary of input data")')
    write (uout,'("Number of atoms: ",A)') string(natoms)
    write (uout,'("List of atoms:")')
    write (uout,'("# Id At lmax")')
    do i = 1, natoms
       write (uout,'(2X,99(A,X))') string(i,2,ioj_left), string(atom(i),3), string(lname(lmax(i)),2)
    end do
    write (uout,'("List of exponents: ")')
    write (uout,'("# Id n      exponent ")')
    do i = 1, nexp
       write (uout,'(2X,99(A,X))') string(i,2,ioj_left), string(nval(i),3), string(eexp(i),'f',12,8,4)
    end do
    write (uout,'("Size of the fitting set: ",A)') string(nfit)
    write (uout,'("List of evaluation sets: ")')
    write (uout,'("#Id  ini  num  step name")')
    do i = 1, nset
       write (uout,'(2X,99(A,X))') string(i,2,ioj_left), string(iset_ini(i),5), string(iset_n(i),4), &
          string(iset_step(i),3), string(iset_label(i))
    end do
    write (uout,'("Number of columns: ",A)') string(ncols)
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

  ! Print the results of an evaluation
  subroutine global_printeval(label,y,stat,ofile)
    use tools_io, only: ferror, faterr, uout, string, ioj_left, ioj_center, ioj_right, &
       fopen_write, fclose
    use types, only: stats
    
    character*(*), intent(in) :: label !< label 
    type(stats), intent(in), optional :: stat !< statistics
    real*8, intent(in), optional :: y(:) !< list of individual results 
    character*(*), intent(in), optional :: ofile !< output file
    
    integer :: i, lu
    logical :: doclose

    lu = uout
    doclose = .false.
    if (present(ofile)) then
       if (len_trim(ofile) > 0) then
          lu = fopen_write(ofile)
          doclose = .true.
          write (uout,'("+ Evaluation (",A,") written to file: ",A/)') string(label), string(ofile)
       end if
    end if

    write (lu,'("# Evaluation: ",A)') string(label)

    if (present(stat)) then
       write (lu,'("# Statistics: ")')
       write (lu,'("#   norm =    ",F12.6)') stat%norm
       write (lu,'("#   maxcoef = ",F12.6)') stat%maxcoef
       write (lu,'("#   wrms =    ",F14.8)') stat%wrms
       do i = 1, nset
          write (lu,'("#",3X,A," rms = ",A," mae = ",A)') &
             string(iset_label(i),10,ioj_left), string(stat%rms(i),'f',14,8), &
             string(stat%mae(i),'f',14,8)
       end do
    end if

    if (present(y)) then
       if (size(y,1) /= nfit) &
          call ferror("global_printeval","inconsistent size of y",faterr)

       write (lu,'("# Id                     Name                     weight        yscf                ytotal                yref                 diff")')
       do i = 1, nfit
          write (lu,'(99(A,X))') string(i,6,ioj_left), string(names(i),maxnamelen,ioj_center), &
             string(w(i),'f',5,1,ioj_right), string(y(i),'f',20,10,8), string(y(i)+ydisp(i),'f',20,10,8), &
             string(yref(i),'f',20,10,8), string(y(i)+ydisp(i)-yref(i),'f',20,10,8)
       end do
       write (lu,*)
    end if

    if (doclose) call fclose(lu)

  end subroutine global_printeval

end module global
