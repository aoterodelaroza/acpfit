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
  use tools_io
  use tools
  use global
  implicit none

  character(len=:), allocatable :: optv, subline, word, word2, line, aux
  integer :: lp, i, ii, j, k, idx, idum, idum2, iatom, iexp, il
  integer :: natoms_ang
  logical :: ok
  real*8 :: rdum

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
  natoms_ang = 0
  call global_init()

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
        word = getword(line,lp)
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
        end if
     elseif (len_trim(word) > 0) then
        call ferror("acpfit","unknown keyword: " // word,faterr)
     end if
  end do

  ! check the input data for sanity
  call global_check(natoms_ang)

  ! build the default energy file names
  allocate(efile(natoms,maxlmax,nexp))
  do i = 1, natoms
     do j = 1, lmax(i)
        do k = 1, nexp
           efile(i,j,k) = lower(trim(atom(i))) // "_" // lname(j) // "_" // string(k) // ".dat"
        end do
     end do
  end do

  ! parse the efile lines, if any available
  if (nefilesi > 0) then
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
  end if

  ! add the "all" set
  nset = nset + 1
  call realloc(iset_label,nset)
  call realloc(iset_ini,nset)
  call realloc(iset_n,nset)
  call realloc(iset_step,nset)
  iset_label(nset) = "all"
  iset_ini(nset) = 1
  iset_n(nset) = nfit
  iset_step(nset) = 1

  ! some output
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

  ! do stuff

  write (uout,'("ACPFIT ended succesfully (",A," WARNINGS, ",A," COMMENTS)")')&
     string(nwarns), string(ncomms)
  call print_clock()
  call tictac('ACPFIT finished')

999 continue

end program acpfit

