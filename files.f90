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
module files
  implicit none

  private

  public :: makefilenames
  public :: readfiles
  private :: gety

contains

  !> Construct the default file names containing the data for the run,
  !> then parse the user's input for manually changing any of these
  !> names.  Requires a complete set of input data (natoms, atom, lmax, etc).
  subroutine makefilenames(nefilesi,efilei)
    use tools_io, only: lgetword, isinteger, isreal, ferror, faterr, lower, equal, string
    use global, only: natoms, atom, lmax, maxlmax, lname, nexp, nval, eexp, efile
    integer, intent(in) :: nefilesi 
    character*(*), intent(in) :: efilei(:) 
    
    integer :: i, j, k, lp, iatom, il, ii, iexp, idum, idum2
    character(len=:), allocatable :: line, word, word2
    logical :: ok
    real*8 :: rdum

    ! build the default energy file names
    if (allocated(efile)) deallocate(efile)
    allocate(efile(natoms,maxlmax,nexp))
    do i = 1, natoms
       do j = 1, lmax(i)
          do k = 1, nexp
             efile(i,j,k) = lower(trim(atom(i))) // "_" // lname(j) // "_" // string(k) // ".dat"
          end do
       end do
    end do

    ! parse user's input
    do i = 1, nefilesi
       lp = 1
       line = efilei(i)
       word = lgetword(line,lp)
       word2 = lgetword(line,lp)
       ok = isinteger(idum2,line,lp)
       ok = ok .and. isreal(rdum,line,lp)
       if (.not.ok) &
          call ferror("acpfit","incorrect syntax in FILE ETERM",faterr)

       ! atom - iatom
       iatom = 0
       do ii = 1, natoms
          if (equal(word,lower(atom(ii)))) then
             iatom = ii
             exit
          end if
       end do
       if (iatom == 0) &
          call ferror("acpfit","unknown atom in FILE ETERM",faterr)

       ! the angular momentum channel - il
       il = 0
       do ii = 1, lmax(iatom)
          if (word2(1:1) == lname(ii)) then
             il = ii
             exit
          end if
       end do
       if (iatom == 0) &
          call ferror("acpfit","unknown ang. mom. label in FILE ETERM",faterr)

       ! the exponent
       iexp = 0
       do ii = 1, nexp
          if (idum2 == nval(ii) .and. abs(rdum-eexp(ii)) < 1d-10) then
             iexp = ii
             exit
          end if
       end do
       if (iexp == 0) &
          call ferror("acpfit","unknown n/exponent in FILE ETERM",faterr)

       ! the file
       efile(iatom,idum,iexp) = trim(line(lp:))
    end do

  end subroutine makefilenames

  subroutine readfiles()
    use global, only: datapath, emptyfile, yempty, yref, reffile, w, wfile, &
       namesfile, names, nfit, ydisp, nsubfiles, subfile, natoms, lmax, nexp, &
       x, efile, ncols, ytarget, maxnamelen, wmask, nfitw, xw, xwork, ywtarget
    use tools_io, only: ferror, faterr, fopen_read, getline, fclose, string

    character(len=:), allocatable :: file, line
    logical :: ok
    integer :: lu, i, j, k, n, istat
    real*8 :: xaux(nfit)

    ! w
    if (allocated(w)) deallocate(w)
    allocate(w(nfit),stat=istat)
    call checkstat(istat,"w")
    call gety(trim(datapath) // trim(wfile),"weights",nfit,w)

    ! empty
    if (allocated(yempty)) deallocate(yempty)
    allocate(yempty(nfit),stat=istat)
    call checkstat(istat,"empty")
    call gety(trim(datapath) // trim(emptyfile),"empty",nfit,yempty)

    ! ref
    if (allocated(yref)) deallocate(yref)
    allocate(yref(nfit),stat=istat)
    call checkstat(istat,"ref")
    call gety(trim(datapath) // trim(reffile),"reference",nfit,yref)

    ! names
    if (allocated(names)) deallocate(names)
    allocate(names(nfit),stat=istat)
    call checkstat(istat,"names")
    file = trim(datapath) // trim(namesfile)
    inquire(file=file,exist=ok)
    if (.not.ok) &
       call ferror("readfiles","File not found (names): " // trim(file),faterr)
    lu = fopen_read(file)
    maxnamelen = 0
    do i = 1, nfit
       ok = getline(lu,line,.true.)
       names(i) = adjustl(line)
       maxnamelen = max(len_trim(names(i)),maxnamelen)
    end do
    call fclose(lu)

    ! disp files
    if (allocated(ydisp)) deallocate(ydisp)
    allocate(ydisp(nfit),stat=istat)
    call checkstat(istat,"ydisp")
    ydisp = 0d0
    do i = 1, nsubfiles
       call gety(trim(datapath) // trim(subfile(i)),"subfile-"//string(i),nfit,xaux)
       ydisp = ydisp + xaux
    end do

    ! make the target
    if (allocated(ytarget)) deallocate(ytarget)
    allocate(ytarget(nfit),stat=istat)
    call checkstat(istat,"ytarget")
    ytarget = yref - ydisp - yempty

    ! energy terms 
    if (allocated(x)) deallocate(x)
    allocate(x(nfit,ncols),stat=istat)
    call checkstat(istat,"x")
    n = 0
    do i = 1, natoms
       do j = 1, lmax(i)
          do k = 1, nexp
             n = n + 1
             call gety(trim(datapath) // trim(efile(i,j,k)),&
                "eterm-"//string(i)//"-"//string(j)//"-"//string(k),&
                nfit,xaux)
             x(:,n) = xaux - yempty
          end do
       end do
    end do
    
    ! build the mask of non-zero weights
    if (allocated(wmask)) deallocate(wmask)
    allocate(wmask(nfit),stat=istat)
    call checkstat(istat,"wmask")
    wmask = .false.
    do i = 1, nfit
       if (abs(w(i)) > 1d-10) then
          wmask(i) = .true.
       end if
    end do
    nfitw = count(wmask)

    ! build the weighted variables
    if (allocated(xw)) deallocate(xw)
    if (allocated(xwork)) deallocate(xwork)
    if (allocated(ywtarget)) deallocate(ywtarget)
    allocate(xw(nfitw,ncols),stat=istat)
    call checkstat(istat,"xw")
    allocate(xwork(nfitw,ncols),stat=istat)
    call checkstat(istat,"xwork")
    allocate(ywtarget(nfitw),stat=istat)
    call checkstat(istat,"ywtarget")
    n = 0
    do i = 1, nfit
       if (abs(w(i)) > 1d-10) then
          n = n + 1
          ywtarget(n) = ytarget(i) * sqrt(w(i))
          xw(n,:) = x(i,:) * sqrt(w(i))
       endif
    end do
    xwork = 0d0

  contains
    subroutine checkstat(istat,label)
      integer, intent(in) :: istat
      character*(*), intent(in) :: label

      if (istat /= 0) &
         call ferror("readfiles","could not allocate memory for "//trim(label),faterr)

    end subroutine checkstat
  end subroutine readfiles

  ! Read n real numbers from the file (identified by label) and return
  ! them in y. Each number is assumed to be in its own line.
  subroutine gety(file,label,n,y)
    use tools_io, only: ferror, faterr, fopen_read, fclose
    character*(*), intent(in) :: file
    character*(*), intent(in) :: label
    integer, intent(in) :: n
    real*8, intent(out) :: y(n)
    
    logical :: ok
    integer :: lu, i

    inquire(file=file,exist=ok)
    if (.not.ok) &
       call ferror("gety","File not found (" // trim(label) // "): " // trim(file),faterr)
    lu = fopen_read(file)
    do i = 1, n
       read(lu,*,err=999) y(i)
    end do
    call fclose(lu)

    return
999 continue
    call ferror("gety","Error reading file ("//trim(label)//"): "//trim(file),faterr)

  end subroutine gety

end module files
