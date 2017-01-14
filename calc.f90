module calc
  implicit none

  private
  public :: calc_stats
  public :: runfit_inf

contains
  
  subroutine calc_stats(y,norm,maxcoef,stat)
    use types, only: stats
    use global, only: nfit, maxnamelen, w, ydisp, yref, names, nset, ytarget, &
       iset_ini, iset_step, iset_n

    real*8, intent(in) :: y(:) ! results of the fit (no dispersion)
    real*8, intent(in) :: norm, maxcoef ! norm and maxcoef for this evaluation
    type(stats), intent(inout) :: stat ! statistics of the fit

    integer :: i, j, id
    real*8 :: dy(nfit)

    ! write down some constants
    stat%norm = norm
    stat%maxcoef = maxcoef
    stat%nset = nset

    ! calculate the wrms
    stat%wrms = sqrt(sum(w * (y-ytarget)**2) / sum(w))
    
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

  subroutine runfit_inf()

    write (*,*) "bleh!"
    stop 1

    ! integer :: i, j, k, l, n1, rank, lwork, info, ihere, nall
    ! integer, allocatable :: icur(:)
    ! real*8, allocatable :: xsub(:,:), xsubw(:,:), dy(:), dyw(:), coef(:), work(:)
    ! integer, allocatable :: jpvt(:)
    ! real*8 :: norm, wrms, rms, mae
    
    ! ! count terms in the least squares
    ! allocate(icur(size(idx)))
    ! n1 = count(idx /= 0)
    ! icur(1:n1) = pack(idx,idx/=0)

    ! ! work space for lapack
    ! allocate(jpvt(n1))
    ! jpvt = 0
    ! lwork = n1 + 3 * n + 1
    ! allocate(work(lwork))

    ! ! allocate and populate the arrays
    ! nall = size(x0,1)
    ! allocate(xsubw(size(xw,1),n1),xsub(nall,n1),dy(nall),dyw(size(yw)),coef(size(yw)))

    ! ! least squares
    ! if (n1 > 0) then
    !    xsubw = xw(:,icur(1:n1))
    !    coef = yw
    !    call dgelsy(n,n1,1,xsubw,n,coef,n,jpvt,1d-10,rank,work,lwork,info)
    !    norm = sqrt(sum(coef(1:n1)**2) * 1d-6)
    ! else
    !    coef = 0d0
    !    norm = 0d0
    ! endif
    
    ! ! calculate the wrms
    ! xsubw(:,1:n1) = xw(:,icur(1:n1))
    ! dyw = yw
    ! call dgemv('n',n,n1,-1d0,xsubw,n,coef(1:n1),1,1d0,dyw,1)
    ! if (abs(sum(w)) > 1d-10) then
    !    wrms = sqrt(sum(dyw**2) / sum(w))
    ! else
    !    wrms = 0d0
    ! endif

    ! ! calculate the residues for the whole set
    ! xsub(:,1:n1) = x0(:,icur(1:n1))
    ! dy = y0
    ! call dgemv('n',nall,n1,-1d0,xsub,nall,coef(1:n1),1,1d0,dy,1)

    ! write (*,'("  Statistics of the results: ")')
    ! write (*,'("    norm = ",F12.6)') norm
    ! write (*,'("    wrms(all) = ",F14.8)') wrms
    ! do i = 1, size(iset)
    !    rms = sqrt(sum(dy(iset(i)%idx)**2) / size(iset(i)%idx))
    !    mae = sum(abs(dy(iset(i)%idx))) / size(iset(i)%idx)
    !    write (*,'(4X,A8," rms = ",F14.8," mae = ",F14.8)') &
    !       iset(i)%name, rms, mae
    ! end do
    ! write (*,*)

    ! if (doeval) then
    !    write (*,'("  Evaluation: ")')
    !    write (*,'("  | id | Name | w | yref | yempty | ycalc |")')
    !    do i = 1, size(dy)
    !       write (*,'("    | ",I4," | ",A," | ",F5.1," | ",3(F14.8," | "))') &
    !          i, trim(name(i)), w(i), yref(i), yempty(i), yempty(i)+y0(i)-dy(i)
    !    end do
    !    write (*,*)
    ! end if

    ! if (dodcp) then
    !    write (*,'("  DCP in Gaussian form: ")')
    !    l = 0
    !    do i = 1, natom
    !       write (*,'(A,X,"0")') trim(atomlist(i))
    !       write (*,'(A,X,"3 0")') trim(atomlist(i))
    !       do j = 1, nchan
    !          write (*,'(A)') trim(channame(j))
    !          write (*,'(I2)') count((icur >= l+1) .and. (icur <= l+nterms))
    !          do k = 1, nterms
    !             l = l + 1
    !             if (any(icur == l)) then
    !                ihere = maxloc(icur,1,icur == l)
    !                write (*,'("2 ",F14.10,X,F14.10)') explist(k), coef(ihere) * 0.001d0
    !             end if
    !          end do
    !       end do
    !    end do
    !    write (*,*)
    ! end if

    ! deallocate(xsub,xsubw,dy,dyw,coef,work,jpvt,icur)


  end subroutine runfit_inf

end module calc
