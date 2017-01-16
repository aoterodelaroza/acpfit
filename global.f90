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

  public :: global_init
  private :: global_check
  public :: global_printinfo
  public :: global_fillcol
  public :: global_printeval
  public :: global_printacp
  public :: global_input
  public :: whichatom
  public :: whichl
  public :: whichcol

  ! file prefix
  character(len=:), allocatable :: fileroot !< file prefix

  ! global input data
  character(len=:), allocatable :: datapath !< path to the data files
  integer :: natoms !< number of atoms
  character*2, allocatable :: atom(:) !< atomic symbols
  integer, allocatable :: lmax(:) !< maximum angular momentum for each atom
  integer :: maxlmax !< maximum lmax
  real*8 :: coef0 !< coefficients with which the ACP terms were evaluated

  ! list of exponents and n-values
  integer :: nexp !< number of exponents
  integer, allocatable :: nval(:) !< n-value associated to the exponent
  real*8, allocatable :: eexp(:) !< exponent value

  ! fitting set, and subsets
  integer :: nfit !< number of systems in the fitting set
  integer :: nfitw !< number of systems in the fitting set with non-zero wegiht
  integer :: nset !< number of subsets
  character*255, allocatable :: iset_label(:) !< labels for the subsets
  integer, allocatable :: iset_ini(:), iset_n(:), iset_step(:) !< initial, number, and step for the subsets

  ! named files
  character(len=:), allocatable :: emptyfile !< name of the empty file
  character(len=:), allocatable :: reffile !< name of the reference energy file
  character(len=:), allocatable :: wfile !< name of the weight file
  character(len=:), allocatable :: namesfile !< name of the names file
  character(len=:), allocatable :: outempty !< name of the output file
  character(len=:), allocatable :: outeval !< name of the evaluation output file
  character(len=:), allocatable :: outacp !< name of the acp output file
  character(len=:), allocatable :: inacp !< name of the acp input file
  character*255, allocatable :: efile(:,:,:) !< name of the energy terms files
  character*255, allocatable :: subfile(:) !< name of the subtract energy file
  integer :: nsubfiles

  ! labels for the angular momentum channels
  character*1, parameter :: lname(6) = (/"l","s","p","d","f","g"/) !< names of the angmom channels

  ! column information
  integer :: ncols
  type(column), allocatable :: col(:)

  ! data for the fits
  logical, allocatable :: wmask(:) !< mask of non-zero weights
  real*8, allocatable :: yempty(:), yref(:), ydisp(:), w(:) !< empty, ref, disp, and weight data
  real*8, allocatable :: x(:,:) !< matrix containing the BSIP term energies
  real*8, allocatable :: xw(:,:) !< weighted x
  real*8, allocatable :: ytarget(:) !< target of the fit
  real*8, allocatable :: ywtarget(:) !< weighted target y
  character*128, allocatable :: names(:) !< names of the systems in the fitting set
  integer :: maxnamelen !< maximum length of the names

  ! working space for lapack
  integer, allocatable :: jpvt(:) !< work array for lapack
  real*8, allocatable :: work(:) !< work array for lapack

  ! run modes
  integer, parameter :: imode_fit = 1
  integer, parameter :: imode_fit_manual = 2
  integer, parameter :: imode_eval = 3
  integer, parameter :: imode_eval_file = 4
  integer, parameter :: imode_test = 5
  integer, parameter :: imode_no = 6

contains
  
  !> Initialize all global variables
  subroutine global_init()

    natoms = 0
    datapath = ""
    emptyfile = "empty.dat"
    reffile = "ref.dat"
    wfile = "w.dat"
    namesfile = "names.dat"
    outempty = ""
    outeval = ""
    outacp = ""
    nexp = 0
    nsubfiles = 0
    nfit = 0
    nset = 0
    coef0 = 1d-3
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

    integer :: i

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
    write (uout,'("Coefficent ACP term evaluation: ",A)') string(coef0,'f',12,6)
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
       write (lu,'("#   norm =    ",A)') string(stat%norm,'f',12,6,ioj_left)
       write (lu,'("#   maxcoef = ",A)') string(stat%maxcoef,'f',12,6,ioj_left)
       write (lu,'("#   wrms =    ",A)') string(stat%wrms,'f',14,8,ioj_left)
       do i = 1, nset
          write (lu,'("#",3X,A," rms = ",A," mae = ",A)') &
             string(iset_label(i),10,ioj_left), string(stat%rms(i),'f',14,8), &
             string(stat%mae(i),'f',14,8)
       end do
    end if

    if (present(y)) then
       if (size(y,1) /= nfit) &
          call ferror("global_printeval","inconsistent size of y",faterr)

       write (lu,'("# Id                     Name                     weight       yempty               yscf                 ytotal               yref                 diff")')
       do i = 1, nfit
          write (lu,'(99(A,X))') string(i,6,ioj_left), string(names(i),maxnamelen,ioj_center), &
             string(w(i),'f',5,1,ioj_right), string(yempty(i),'f',20,10,8), &
             string(yempty(i)+y(i),'f',20,10,8), string(yempty(i)+y(i)+ydisp(i),'f',20,10,8), &
             string(yref(i),'f',20,10,8), string(yempty(i)+y(i)+ydisp(i)-yref(i),'f',20,10,8)
       end do
       write (lu,*)
    end if

    if (doclose) call fclose(lu)

  end subroutine global_printeval

  !> Print an ACP in Gaussian form
  subroutine global_printacp(label,ndim,idx,coef,ofile)
    use tools_io, only: uout, string, fopen_write, fclose
    character*(*), intent(in) :: label !< label 
    integer, intent(in) :: ndim !< number of terms (total)
    integer, intent(in) :: idx(ndim) !< columns for the terms
    real*8, intent(in) :: coef(ndim) !< coefficients
    character*(*), intent(in), optional :: ofile !< output file

    integer :: catom, cang, i, id, j, n, lu, nl(natoms)
    logical :: doclose

    ! open the output lu
    lu = uout
    doclose = .false.
    if (present(ofile)) then
       if (len_trim(ofile) > 0) then
          lu = fopen_write(ofile)
          doclose = .true.
          write (uout,'("+ ACP (",A,") written to file: ",A/)') string(label), string(ofile)
       end if
    end if

    if (lu == uout) then
       write (lu,'("# ACP (",A,")")') string(label)
    end if

    ! maximum l for each atom
    nl = 0
    do i = 1, ndim
       nl(col(idx(i))%iatom) = max(nl(col(idx(i))%iatom),col(idx(i))%l-1)
    end do

    catom = -1
    cang = 0
    do i = 1, ndim
       id = idx(i)
       if (col(id)%iatom /= catom) then
          catom = col(id)%iatom
          cang = 0
          write (lu,'("-",A," 0")') string(atom(catom))
          write (lu,'(A,X,A," 0")') string(atom(catom)), string(nl(catom))
       end if
       if (col(id)%l /= cang) then
          do j = cang+1, col(id)%l-1
             write (lu,'(A)') lname(j)
             write (lu,'("0")') 
          end do
          cang = col(id)%l
          write (lu,'(A)') lname(cang)
          n = 1
          do j = i+1, ndim
             if (col(idx(j))%iatom == catom .and. col(idx(j))%l == cang) then
                n = n + 1
             else
                exit
             end if
          end do
          write (lu,'(A)') string(n)
       end if
       write (lu,'(3(A,X))') string(col(id)%n), string(col(id)%eexp,'f',20,12),&
          string(coef(i)*coef0,'f',20,12,6)
    end do

    if (doclose) then
       call fclose(lu)
    else
       write (uout,*)
    end if

  end subroutine global_printacp

  ! read and parse input
  subroutine global_input(nefilesi,efilei,imode,ifit_n,fit_maxnorm,fit_maxcoef)
    use tools_io, only: uin, getline, lgetword, equal, getword, faterr, ferror,&
       isinteger, isreal
    use types, only: realloc
    integer, intent(out) :: nefilesi ! number of files
    character*255, allocatable, intent(out) :: efilei(:) ! name of files
    integer, intent(out) :: ifit_n 
    integer, intent(out) :: imode
    real*8, intent(out) :: fit_maxnorm
    real*8, intent(out) :: fit_maxcoef

    character(len=:), allocatable :: word, line, subline, aux
    integer :: natoms_ang, lp, idum, i, idx, lp2
    real*8 :: rdum
    logical :: ok

    ! initialize
    natoms_ang = 0
    nefilesi = 0
    if (allocated(efilei)) deallocate(efilei)
    allocate(efilei(1))
    fit_maxnorm = huge(1d0)
    fit_maxcoef = huge(1d0)
    imode = imode_no

    ! parse the input
    do while (getline(uin,line))
       lp=1
       word = lgetword(line,lp)
       subline = line(lp:)
       if (equal(word,'atom').or.equal(word,'atoms')) then
          ! ATOM|ATOMS at1.s at2.s ...
          do while (.true.)
             word = getword(line,lp)
             if (len_trim(word) < 1) exit
             natoms = natoms + 1
             if (natoms > ubound(atom,1)) call realloc(atom,2*natoms)
             atom(natoms) = word
          end do
          call realloc(atom,natoms)
       elseif (equal(word,'lmax')) then
          ! LMAX at.s maxl.i
          do while (.true.)
             word = lgetword(line,lp)
             if (len_trim(word) == 0) exit
             idum = 0
             do i = 1, ubound(lname,1)
                if (word(1:1) == lname(i)) then
                   idum = i
                   exit
                end if
             end do
             if (idum == 0) &
                call ferror("acpfit","wrong angular momentum label in LMAX",faterr)

             natoms_ang = natoms_ang + 1
             if (natoms_ang > ubound(lmax,1)) call realloc(lmax,2*natoms_ang)
             lmax(natoms_ang) = idum
          end do
          call realloc(lmax,natoms_ang)
       elseif (equal(word,'datapath')) then
          ! DATAPATH path.s
          datapath = trim(adjustl(subline))
          idx = len_trim(datapath)
          if (datapath(idx:idx) /= "/") then
             aux = trim(datapath) // "/"
             datapath = aux
          end if
       elseif (equal(word,'exp').or.equal(word,'exponent').or.equal(word,'exponents')) then
          ! EXP|EXPONENT|EXPONENTS n.i exp1.r exp2.r ...
          ok = isinteger(idum,line,lp)
          if (.not.ok) &
             call ferror("acpfit","wrong syntax in EXP",faterr)

          do while (.true.)
             ok = isreal(rdum,line,lp)
             if (.not.ok) exit
             nexp = nexp + 1
             if (nexp > ubound(nval,1)) then
                call realloc(nval,2*nexp)
                call realloc(eexp,2*nexp)
             end if
             nval(nexp) = idum
             eexp(nexp) = rdum
          end do
          call realloc(nval,nexp)
          call realloc(eexp,nexp)
       elseif (equal(word,'coef0')) then
          ! COEF0 coef0.r
          ok = isreal(coef0,line,lp)
          if (.not.ok) &
             call ferror("acpfit","wrong COEF0 syntax",faterr)
       elseif (equal(word,'nfit')) then
          ! NFIT nfit.i
          ok = isinteger(nfit,line,lp)
          if (.not.ok) &
             call ferror("acpfit","wrong NFIT syntax",faterr)
       elseif (equal(word,'set')) then
          ! SET label.s ini.i end.i step.i
          nset = nset + 1
          if (nset > ubound(iset_ini,1)) then
             call realloc(iset_label,2*nset)
             call realloc(iset_ini,2*nset)
             call realloc(iset_n,2*nset)
             call realloc(iset_step,2*nset)
          end if
          iset_label(nset) = getword(line,lp)
          ok = isinteger(iset_ini(nset),line,lp)
          ok = ok .and. isinteger(iset_n(nset),line,lp)
          if (.not.ok) &
             call ferror("acpfit","wrong SET syntax",faterr)
          ok = isinteger(iset_step(nset),line,lp)
          if (.not.ok) iset_step(nset) = 1
       elseif (equal(word,'file')) then
          ! FILE ...
          word = lgetword(line,lp)
          if (equal(word,'empty')) then
             ! FILE EMPTY emptyfile.s
             emptyfile = trim(line(lp:))
          elseif (equal(word,'ref')) then
             ! FILE REF reffile.s 
             reffile = trim(line(lp:))
          elseif (equal(word,'w')) then
             ! FILE W wfile.s
             wfile = trim(line(lp:))
          elseif (equal(word,'names')) then
             ! FILE NAMES namesfile.s
             namesfile = trim(line(lp:))
          elseif (equal(word,'sub')) then
             ! FILE SUB subfile.s
             nsubfiles = nsubfiles + 1
             if (nsubfiles > ubound(subfile,1)) &
                call realloc(subfile,2*nsubfiles)
             subfile(nsubfiles) = trim(line(lp:))
          elseif (equal(word,'eterm')) then
             ! FILE ETERM at.s l.i n.i exp.r efile.s
             nefilesi = nefilesi + 1
             if (nefilesi > ubound(efilei,1)) &
                call realloc(efilei,2*nefilesi)
             efilei(nefilesi) = trim(line(lp:))
          else
             call ferror("acpfit","unknown FILE keyword: " // word,faterr)
          end if
       elseif (equal(word,'output')) then
          ! OUTPUT ...
          word = lgetword(line,lp)
          if (equal(word,'empty')) then
             ! OUTPUT EMPTY file.s
             outempty = trim(line(lp:))
          elseif (equal(word,'eval')) then
             ! OUTPUT EMPTY file.s
             outeval = trim(line(lp:))
          elseif (equal(word,'acp')) then
             ! OUTPUT ACP file.s
             outacp = trim(line(lp:))
          else
             call ferror("acpfit","unknown OUTPUT keyword: " // word,faterr)
          end if
       elseif (equal(word,'run')) then
          ! RUN ...
          word = lgetword(line,lp)
          if (equal(word,'fit')) then
             ! RUN FIT ...
             imode = imode_fit
             lp2 = lp
             word = lgetword(line,lp)
             ifit_n = 0
             if (equal(word,"inf")) then
                ! RUN FIT INF ...
                ifit_n = -1
             elseif (len_trim(word) < 1) then
                imode = imode_fit_manual
                ifit_n = 0
             else
                lp = lp2
                ok = isinteger(idum,line,lp)
                if (ok) then
                   if (idum > 0) then
                      ifit_n = idum
                   end if
                end if
                if (ifit_n == 0) &
                   call ferror("acpfit","wrong RUN FIT syntax",faterr)
                ! RUN ... [MAXNORM norm.r] [MAXCOEF coef.r]
                do while(.true.)
                   word = lgetword(line,lp)
                   if (equal(word,"maxnorm")) then
                      ok = isreal(fit_maxnorm,line,lp)
                      if (.not.ok) &
                         call ferror("acpfit","wrong RUN FIT MAXNORM syntax",faterr)
                   elseif (equal(word,"maxcoef")) then
                      ok = isreal(fit_maxcoef,line,lp)
                      if (.not.ok) &
                         call ferror("acpfit","wrong RUN FIT MAXCOEF syntax",faterr)
                   elseif (len_trim(word) > 0) then
                      call ferror("acpfit","unknown RUN FIT keyword: " // word,faterr)
                   else
                      exit
                   end if
                end do
             end if
          elseif (equal(word,'eval')) then
             word = getword(line,lp)
             if (len_trim(word) == 0) then
                ! RUN EVAL ...
                imode = imode_eval
             else
                ! RUN EVAL file.acp
                imode = imode_eval_file
                inacp = word
             endif
          elseif (equal(word,'test')) then
             imode = imode_test
          elseif (len_trim(word) > 0) then
             call ferror("acpfit","unknown RUN keyword: " // word,faterr)
          else
             ! RUN
             imode = imode_no
          endif
          exit
       elseif (len_trim(word) > 0) then
          call ferror("acpfit","unknown keyword: " // word,faterr)
       end if
    end do

    ! check the input data for sanity
    call global_check(natoms_ang)

  end subroutine global_input

  ! Identify an atom by its name. Zero if the atom is not found.
  function whichatom(at)
    use tools_io, only: lower, equal
    character*(*), intent(in) :: at
    integer :: whichatom

    integer :: i
    character(len=:), allocatable :: at0, at1

    at0 = lower(trim(adjustl(at)))
    if (at0(1:1) == "-") then
       at1 = at0(2:)
    else
       at1 = at0
    endif

    whichatom = 0
    do i = 1, natoms
       if (equal(lower(trim(atom(i))),at1)) then
          whichatom = i
          exit
       end if
    end do

  end function whichatom

  ! Identify an angular momentum channel by its name. Zero if it
  ! is not found.
  function whichl(ll)
    use tools_io, only: lower
    character*1, intent(in) :: ll
    integer :: whichl

    integer :: i

    whichl = 0
    do i = 1, maxlmax
       if (lower(ll) == lname(i)) then
          whichl = i
          exit
       end if
    end do

  end function whichl

  ! Determine the column that corresponds to the given atom (iat)
  ! angular momentum (il), n value (n), and exponent (e), or 0 if it
  ! does not exist.
  function whichcol(iat,l,n,e)
    integer, intent(in) :: iat, l, n
    real*8, intent(in) :: e
    integer :: whichcol

    integer :: i

    whichcol = 0
    do i = 1, ncols
       if (col(i)%iatom == iat .and. col(i)%l == l .and.&
          col(i)%n == n .and. abs(col(i)%eexp - e) < 1d-10) then
          whichcol = i
          exit
       end if
    end do

  end function whichcol

end module global
