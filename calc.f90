module calc
  implicit none

  private
  public :: calc_stats
  public :: runfit_inf
  public :: runfit_scanatom
  public :: runfit_manual
  private :: choose
  private :: comb

  integer, parameter :: ilong = selected_int_kind(32)

contains
  
  ! calculate the statistics for a given evaluation
  subroutine calc_stats(y,stat,coef)
    use types, only: stats
    use global, only: nfit, maxnamelen, w, ydisp, yref, names, nset, ytarget, &
       iset_ini, iset_step, iset_n, coef0

    real*8, intent(in) :: y(:) ! results of the fit (no dispersion)
    type(stats), intent(inout) :: stat ! statistics of the fit
    real*8, intent(in), optional :: coef(:) ! coefficients

    integer :: i, j, id
    real*8 :: dy(nfit)

    ! write down some constants
    if (present(coef)) then
       stat%norm = sqrt(sum(coef**2)) * coef0
       stat%maxcoef = maxval(abs(coef)) * coef0
    else
       stat%norm = 0d0
       stat%maxcoef = 0d0
    end if
    stat%nset = nset

    ! calculate the wrms
    if (abs(sum(w)) < 1d-10) then
       stat%wrms = -1d0
    else
       stat%wrms = sqrt(sum(w * (y-ytarget)**2) / sum(w))
    end if

    ! allocate rms and mae
    if (allocated(stat%rms)) deallocate(stat%rms)
    if (allocated(stat%mae)) deallocate(stat%mae)
    allocate(stat%rms(nset))
    allocate(stat%mae(nset))

    ! calculate rms and mae for all sets
    dy = y - ytarget
    do i = 1, nset
       id = iset_ini(i) - iset_step(i)
       stat%rms(i) = 0d0
       stat%mae(i) = 0d0
       do j = 1, iset_n(i)
          id = id + iset_step(i)
          stat%rms(i) = stat%rms(i) + (dy(id))**2
          stat%mae(i) = stat%mae(i) + abs(dy(id))
       end do
       stat%rms(i) = sqrt(stat%rms(i) / iset_n(i))
       stat%mae(i) = stat%mae(i) / iset_n(i)
    end do

  end subroutine calc_stats

  ! Fit the ACP with all possible terms. Good for testing.
  subroutine runfit_inf(outeval,outacp)
    use global, only: ncols, nfitw, nfit, global_printeval, x, yempty, global_printacp
    use types, only: stats
    use tools_io, only: uout, faterr, ferror, string
    
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: outacp

    integer :: idx(ncols), i
    real*8 :: coef(nfitw), y(nfit)
    type(stats) :: stat

    ! header
    write (uout,'("+ Generating the ACP that uses all available terms, ACP(Inf)")')

    if (ncols > nfitw) &
       call ferror("runfit_inf","Too many columns ("//string(ncols)//") for current data ("//string(nfitw)//")",faterr)

    ! all possible terms
    idx = 0
    do i = 1, ncols
       idx(i) = i
    end do

    ! run least squares
    call lsqr(ncols,idx,coef)
    y = matmul(x(:,idx),coef)

    ! print stats and evaluation
    call calc_stats(y,stat,coef)
    call global_printeval("final",y,stat,outeval)
    
    ! print resulting acp
    call global_printacp("final",ncols,idx,coef,outacp)

  end subroutine runfit_inf

  ! Fit an ACP using an iterative procedure over the atoms. For each
  ! atom, select the combination of nuse terms that minimizes the error
  ! subject to the maximum norm (maxnorm) and/or maximum absolute coefficient
  ! (maxcoef) constraints.
  subroutine runfit_scanatom(nuse,maxnorm,maxcoef,outeval,outacp)
    use global, only: ncols, nfitw, nfit, global_printeval, x, yempty, global_printacp, natoms,&
       atom, maxlmax, lmax, nexp, w, ywtarget, xw, coef0, col, lname
    use types, only: stats
    use tools_io, only: uout, string, ferror, faterr, string, ioj_left
    
    integer, intent(in) :: nuse
    real*8, intent(in) :: maxnorm, maxcoef
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: outacp

    character(len=:), allocatable :: str, aux, aux2
    integer :: i, j, k, ncyc, n1, n2, nanuse, nanz, n, nsame, id
    real*8 :: coef(nfitw), y(nfitw), y0(nfit), wrms, minwrms, norm, minnorm
    real*8 :: acoef, minacoef
    type(stats) :: stat
    integer, allocatable :: iatidx(:,:), natidx(:)
    integer :: co(nuse), idx(natoms*nuse), idx0(natoms*nuse)
    integer :: idxnz(natoms*nuse), idxsave(nuse)

    ! header
    write (uout,'("+ Generating ACP using atom iterations, with ",A," terms per atom")') string(nuse)
    if (maxnorm < huge(1d0)) &
       write (uout,'("  Maximum norm: ",A)') string(maxnorm,'f',12,6)
    if (maxcoef < huge(1d0)) &
       write (uout,'("  Maximum abs(coef): ",A)') string(maxcoef,'f',12,6)
    write (uout,*)

    ! some dimensions
    nanuse = natoms*nuse
    if (nanuse > nfitw) &
       call ferror("runfit_inf","Too many columns ("//string(ncols)//") for current data ("//string(nfitw)//")",faterr)

    ! prepare the sets of atomic indices
    allocate(iatidx(maxlmax*nexp,natoms),natidx(natoms))
    natidx = 0
    iatidx = 0
    n1 = 0
    do i = 1, natoms
       n2 = 0
       do j = 1, lmax(i)
          do k = 1, nexp
             n1 = n1 + 1
             n2 = n2 + 1
             iatidx(n2,i) = n1
             natidx(i) = natidx(i) + 1
          end do
       end do
    end do

    ! run the main loop (over cycles)
    ncyc = 0
    idx0 = 0
    nanz = 0
    nsame = 0
    main: do while(.true.)
       ncyc = ncyc + 1
       write (uout,'("# Cycle ",A)') string(ncyc)

       ! run the loop over atoms
       do i = 1, natoms
          write (uout,'("# Atom ",A," (",A,")")') string(i), string(atom(i))
          idx = idx0
          idxsave = idx0(nuse*(i-1)+1:nuse*i)
          nanz = min(nanz + nuse,nanuse)
          minwrms = huge(1d0)
          minnorm = huge(1d0)
          minacoef = huge(1d0)
          
          ! Run over all the combinations. To save memory and allow
          ! parallelization, use Buckles' algorithm.
          !$omp parallel do private(coef,y,wrms,norm,acoef,co) firstprivate(idx)
          do j = 1, binom(natidx(i),nuse)
             ! get the indices for this combination
             call comb(natidx(i),nuse,j,co)

             ! put them together with the rest of the atoms
             idx(nuse*(i-1)+1:nuse*i) = iatidx(co,i)
             
             ! run least squares and calculate wrms
             call lsqr(nanz,idx(1:nanz),coef,y)
             wrms = sqrt(sum((y-ywtarget)**2) / sum(w))
             norm = sqrt(sum(coef(1:nanz)**2)) * coef0
             acoef = maxval(abs(coef(1:nanz))) * coef0

             ! apply the discard criteria; save minimum wrms
             !$omp critical (save)
             if (wrms < minwrms .and. norm < maxnorm .and. acoef < maxcoef) then
                idx0 = idx
                minwrms = wrms
                minnorm = norm
                minacoef = acoef
             end if
             !$omp end critical (save)
          end do
          !$omp end parallel do

          ! check we had at least one combination
          if (minwrms == huge(1d0)) &
             call ferror("runfit_scanatom","Could not find any combination matching the input criteria",faterr)

          ! increase nsame if the new indices are the same as the old ones
          if (all(idxsave == idx0(nuse*(i-1)+1:nuse*i))) then
             nsame = nsame + 1
          else
             nsame = 1
          end if

          ! some output
          do j = 1, natoms
             str = ""
             do k = 1, nuse
                id = idx0(nuse*(j-1)+k)
                if (id > 0) then
                   aux =  string(id,4,ioj_left) // "(" // string(col(id)%atom) // &
                      "," // string(lname(col(id)%l)) // "," // string(col(id)%n) // &
                      "," // trim(string(col(id)%eexp,'f',10,2,ioj_left)) // ")"
                   aux2 = str // "  " // string(aux,20,ioj_left)
                   str = aux2
                end if
             end do
             if (len_trim(str) > 0) &
                write (uout,'(2X,A,2X)') str
          end do
          write (uout,'("  wrms    = ",A)') string(minwrms,'f',14,8)
          write (uout,'("  norm    = ",A)') string(minnorm,'f',14,8)
          write (uout,'("  maxcoef = ",A)') string(minacoef,'f',14,8)
          write (uout,'("  nsame   = ",A)') string(nsame)

          ! exit if we're done
          if (nsame == natoms) exit main
       end do
    end do main
    write (uout,*)

    ! final results
    call lsqr(nanuse,idx0,coef)
    y0 = matmul(x(:,idx0),coef(1:nanuse))

    ! print stats and evaluation
    call calc_stats(y0,stat,coef(1:nanuse))
    call global_printeval("final",y0,stat,outeval)
    
    ! print resulting acp
    call global_printacp("final",nanuse,idx0,coef(1:nanuse),outacp)

  end subroutine runfit_scanatom

  subroutine runfit_manual(outeval,outacp)
    character*(*), intent(in) :: outeval
    character*(*), intent(in) :: outacp

    write (*,*) "bleh!"
    stop 1

  end subroutine runfit_manual

  ! Run least squares using ndim columns given by idx. Returns the
  ! coefficients in coef(1:ndim).
  subroutine lsqr(ndim,idx,coef,yout)
    use global, only: natoms, ncols, nfitw, xw, ywtarget
    integer, intent(in) :: ndim
    integer, intent(in) :: idx(ndim)
    real*8, intent(out) :: coef(nfitw)
    real*8, intent(out), optional :: yout(nfitw)

    ! work data for lapack
    integer :: jpvt(natoms*ncols) 
    real*8 :: work(natoms*ncols + 3 * nfitw + 1), xwork(nfitw,ndim)
    integer :: lwork, rank, info

    ! initialize lapack
    jpvt = 0
    lwork = natoms*ncols + 3 * nfitw + 1

    ! build the work slice
    ! xwork(:,1:ndim) = xw(:,idx)
    xwork = xw(:,idx)
    coef = ywtarget

    ! run the least squares 
    call dgelsy(nfitw,ndim,1,xwork,nfitw,coef,nfitw,jpvt,1d-10,rank,work,lwork,info)
    if (present(yout)) then
       xwork = xw(:,idx)
       yout = matmul(xwork(:,1:ndim),coef(1:ndim))
    end if

  end subroutine lsqr

  !> Calculate n choose k using the multiplicative formula.
  function choose(n,k)
    integer*8 :: choose
    integer, intent(in) :: n, k
 
    integer, parameter :: qp = selected_real_kind(32)
    integer :: i, s
    real(kind=qp) :: x
 
    if (k > n .or. k < 0) then
       choose = 0
       return
    end if

    s = min(k,n-k)
    x = 1.
    do i = 1, s
       x = (x * (n + 1 - i)) / i
    end do
    choose = nint(x,8)

  end function choose
  
  !> Given a combination of k in n elements (idx), calculate the next
  !> combination in the lexicographical order.
  subroutine nextcomb(idx,n,k)
    integer, intent(in) :: n, k
    integer, intent(inout) :: idx(k)

    integer :: i, j
    
    do i = k, 1, -1
       if (idx(i) < n-k+i)  then
          idx(i) = idx(i) + 1
          do j = i+1, k
             idx(j) = idx(j-1)+1
          end do
          return
       end if
    end do

  end subroutine nextcomb

  ! From TOMS/515
  ! Generates a vector from a lexicographical index That is, let C1,
  ! C2, ... Cm be the set of combinations of n items taken p at a
  ! time arranged in lexographical order. Given an integer i, this
  ! routine finds Ci.
  ! B.P. Buckles and M. Lybanon, ACM TOMS 3 (1977) 180-182
  pure subroutine comb(N, P, L, C)                                       
    ! THIS SUBROUTINE FINDS THE COMBINATION SET OF N THINGS
    ! TAKEN P AT A TIME FOR A GIVEN LEXICOGRAPHICAL INDEX.
    ! N - NUMBER OF THINGS IN THE SET
    ! P - NUMBER OF THINGS IN EACH COMBINATION
    ! L - LEXICOGRAPHICAL INDEX OF COMBINATION SOUGHT
    ! C - OUTPUT ARRAY CONTAINING THE COMBINATION SET
    ! THE FOLLOWING RELATIONSHIPS MUST EXIST AMONG THE INPUT
    ! VARIABLES.  L MUST BE GREATER THAN OR EQUAL TO 1 AND LESS
    ! THAN OR EQUAL TO THE MAXIMUM LEXICOGRAPHICAL INDEX.
    ! P MUST BE LESS THAN OR EQUAL TO N AND GREATER THAN ZERO.
    integer, intent(in) :: n, p, l
    integer, intent(out) :: c(p)

    integer :: K, R, P1, I

    ! SPECIAL CASE CODE IF P = 1
    IF(P.EQ.1)THEN
       C(1) = L
       RETURN
    ENDIF

    ! INITIALIZE LOWER BOUND INDEX AT ZERO
    K = 0
    ! LOOP TO SELECT ELEMENTS IN ASCENDING ORDER
    P1 = P - 1
    C(1) = 0
    DO I=1,P1
       ! SET LOWER BOUND ELEMENT NUMBER FOR NEXT ELEMENT VALUE
       IF (I.NE.1) C(I) = C(I-1)
       ! LOOP TO CHECK VALIDITY OF EACH ELEMENT VALUE
10     C(I) = C(I) + 1
       R = BINOM(N-C(I),P-I)
       K = K + R
       IF (K.LT.L) GO TO 10
       K = K - R
    end DO
    C(P) = C(P1) + L - K

  end subroutine comb

  pure function binom(M, N)
    ! ACM ALGORITHM 160 TRANSLATED TO FORTRAN.  CALCULATES THE
    ! NUMBER OF COMBINATIONS OF M THINGS TAKEN N AT A TIME.
    integer, intent(in) :: m, n
    integer :: binom

    INTEGER :: P, I, N1, R

    N1 = N
    P = M - N1
    IF (N1.GE.P) GO TO 10
    P = N1
    N1 = M - P
10  R = N1 + 1
    IF (P.EQ.0) R = 1
    IF (P.LT.2) GO TO 30
    DO I=2,P
       R = (R*(N1+I))/I
    end DO
30  BINOM = R
  end FUNCTION BINOM

end module calc
