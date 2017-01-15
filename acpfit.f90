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
  use calc, only: calc_stats, runfit_inf
  use files, only: makefilenames, readfiles
  use global
  use types, only: realloc, stats
  use tools_io
  implicit none

  character(len=:), allocatable :: optv, subline, word, word2, line, aux
  integer :: lp, lp2, i, ii, j, k, idx, idum, idum2, iatom, iexp, il
  integer :: natoms_ang
  logical :: ok
  real*8 :: rdum
  type(stats) :: statempty
  character(len=:), allocatable :: outempty, outeval, outacp

  integer :: imode
  integer, parameter :: imode_fit = 1
  integer, parameter :: imode_fit_manual = 2
  integer, parameter :: imode_eval = 3
  integer, parameter :: imode_no = 4
  integer :: ifit_n 
  real*8 :: fit_maxnorm
  real*8 :: fit_maxcoef
  real*8, allocatable :: ydum(:)

  ! energy terms files in manual input
  integer :: nefilesi ! number of files
  character*255, allocatable :: efilei(:) ! name of files

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
  nefilesi = 0
  call global_init()
  if (allocated(efilei)) deallocate(efilei)
  allocate(efilei(1))
  outempty = ""
  outeval = ""
  outacp = ""
  fit_maxnorm = -1d0
  fit_maxcoef = -1d0
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
           else
              lp = lp2
              ok = isinteger(idum,line,lp)
              if (ok) then
                 ! RUN FIT N
                 ifit_n = idum
              else
                 ! RUN FIT (manual op.)
                 imode = imode_fit_manual
              end if
           end if

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
        elseif (equal(word,'eval')) then
           ! RUN EVAL ...
           imode = imode_eval
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

  ! build the file names and process user input re specific energy term file names
  call makefilenames(nefilesi,efilei)

  ! fill the column information
  call global_fillcol()

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

  ! write the input data information
  call global_printinfo()

  ! read the information from the external files
  call readfiles()

  ! calculate the statistics for the empty
  allocate(ydum(nfit))
  ydum = 0d0
  call calc_stats(ydum,statempty)
  call global_printeval("empty",ydum,statempty,outempty)

  ! run the calculation
  if (imode == imode_fit) then
     if (ifit_n < 0) then
        ! all terms in the ACP
        call runfit_inf(outeval,outacp)
     else
     end if
  elseif (imode == imode_fit_manual) then
  elseif (imode == imode_eval) then
  end if

  write (uout,'("ACPFIT ended succesfully (",A," WARNINGS, ",A," COMMENTS)")')&
     string(nwarns), string(ncomms)
  call print_clock()
  call tictac('ACPFIT finished')

999 continue

end program acpfit

