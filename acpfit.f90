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
program acpfit
  use tools_io, only: ioinit, stdargs, help_me, uin, uout, getline, tictac,&
     string, nwarns, ncomms, start_clock, print_clock, lgetword, equal,&
     ferror, faterr, getword, isinteger, ioj_left
  use global, only: fileroot, natoms, atom, realloc, datapath, lmax, lmax0,&
     lname
  implicit none

  character(len=:), allocatable :: optv, subline, word, line, aux
  integer :: lp, i, idx, idum
  integer :: natoms_ang
  logical :: ok

  ! open input and output files
  call ioinit()
  call stdargs(optv,fileroot)

  ! help message and banner
  write (uout,'("* ACPFIT: fit atom-centered potentials using least squares")')
  if (index(optv,'h') > 0) then
     call help_me()
     goto 999
  end if
  
  ! timing information
  call tictac('ACPFIT started')
  call start_clock()
  write (uout,*)

  ! initialize
  natoms = 0
  natoms_ang = 0
  datapath = ""

  ! parse the input
  do while (getline(uin,line))
     lp=1
     word = lgetword(line,lp)
     subline = line(lp:)
     if (equal(word,'atom').or.equal(word,'atoms')) then
        ! ATOM|ATOMS at1.s at2.s ...
        natoms = 0
        if (allocated(atom)) deallocate(atom)
        allocate(atom(1))
        do while (.true.)
           word = getword(line,lp)
           if (len_trim(word) < 1) exit
           natoms = natoms + 1
           if (natoms > ubound(atom,1)) call realloc(atom,2*natoms)
           atom(natoms) = word
        end do
        call realloc(atom,natoms)
     elseif (equal(word,'lmax')) then
        ! ANGULAR at.s maxl.i
        natoms_ang = 0
        if (allocated(lmax)) deallocate(lmax)
        allocate(lmax(1))
        do while (.true.)
           ok = isinteger(idum,line,lp)
           if (.not.ok) exit
           natoms_ang = natoms_ang + 1
           if (natoms_ang > ubound(lmax,1)) call realloc(lmax,2*natoms_ang)
           lmax(natoms_ang) = idum
        end do
        call realloc(atom,natoms)
     elseif (equal(word,'datapath')) then
        ! DATAPATH path.s
        datapath = trim(adjustl(subline))
        idx = len_trim(datapath)
        if (datapath(idx:idx) /= "/") then
           aux = trim(datapath) // "/"
           datapath = aux
        end if
     end if
  end do

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
  if (any(lmax < 0) .or. any(lmax > lmax0)) &
     call ferror("acpfit","Invalid value in LMAX",faterr)

  ! some output
  write (uout,'("+ Summary of input data")')
  write (uout,'("  Data path: ",A)') string(datapath)
  write (uout,'("  Atomic informaton: ")')
  write (uout,'("# Id At lmax")')
  do i = 1, natoms
     write (uout,'(2X,99(A,X))') string(i,2,ioj_left), string(atom(i),3), string(lname(lmax(i)),2)
  end do
  write (uout,*)

  write (uout,'("ACPFIT ended succesfully (",A," WARNINGS, ",A," COMMENTS)")')&
     string(nwarns), string(ncomms)
  call print_clock()
  call tictac('ACPFIT finished')

999 continue

end program acpfit

