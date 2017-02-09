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
module types
  implicit none
  private
  public :: realloc
  public :: column
  public :: stats

  ! overloaded functions
  interface realloc
     module procedure realloc1l
     module procedure realloc1r
     module procedure realloc2r
     module procedure realloc3r
     module procedure realloc4r
     module procedure realloc5r
     module procedure realloc1i2
     module procedure realloc1i
     module procedure realloc2i
     module procedure realloc1c
     module procedure realloc1cmplx4
     module procedure realloc2cmplx4
     module procedure realloc4cmplx4
     module procedure realloc1cmplx8
  end interface

  ! column information
  type column
     character*2 :: atom ! atom
     integer :: iatom ! atom index
     integer :: l ! angular momentum
     integer :: n ! n value
     integer :: iexp ! exponent index
     real*8 :: eexp ! exponent value
  end type column

  ! statistics for an evaluation
  type stats
     real*8 :: norm ! norm of the coefficients
     real*8 :: maxcoef ! maximum absolute value of the coefficients
     real*8 :: wrms ! wrms for the fit
     integer :: nset ! number of subsets
     real*8, allocatable :: rms(:) ! rms for each of the subsets
     real*8, allocatable :: mae(:) ! mae for each of the subsets
  end type stats

contains

  !> Adapt the size of an allocatable 1D logical array
  subroutine realloc1l(a,nnew)
    use tools_io, only: ferror, faterr

    logical, intent(inout), allocatable :: a(:) !< Input array, logical, 1D
    integer, intent(in) :: nnew !< new dimension
    
    logical, allocatable :: temp(:)
    integer :: l1, u1
    
    if (.not.allocated(a)) &
       call ferror('realloc1l','array not allocated',faterr)
    l1 = lbound(a,1)
    u1 = ubound(a,1)
    if (u1 == nnew) return
    allocate(temp(l1:nnew))
    
    temp(l1:min(nnew,u1)) = a(l1:min(nnew,u1))
    call move_alloc(temp,a)

  end subroutine realloc1l

  !> Adapt the size of an allocatable 1D real*8 array
  subroutine realloc1r(a,nnew)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    real*8, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1r','array not allocated',faterr)
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1r

  !> Adapt the size of an allocatable 2D real*8 array
  subroutine realloc2r(a,n1,n2)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    real*8, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) &
       call ferror('realloc2r','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    if (nold(1) == n1 .and. nold(2) == n2) return
    allocate(temp(n1,n2))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc3r(a,n1,n2,n3)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3 !< new dimension
    
    real*8, allocatable :: temp(:,:,:)
    integer :: nold(3)
    
    if (.not.allocated(a)) &
       call ferror('realloc3r','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    if (nold(1) == n1 .and. nold(2) == n2 .and. nold(3) == n3) return
    allocate(temp(n1,n2,n3))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)))
    call move_alloc(temp,a)

  end subroutine realloc3r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc4r(a,n1,n2,n3,n4)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:,:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3, n4 !< new dimension
    
    real*8, allocatable :: temp(:,:,:,:)
    integer :: nold(4)
    
    if (.not.allocated(a)) &
       call ferror('realloc4r','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    nold(4) = size(a,4)
    if (nold(1) == n1 .and. nold(2) == n2 .and. nold(3) == n3 .and.&
        nold(4) == n4) return
    allocate(temp(n1,n2,n3,n4))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)))
    call move_alloc(temp,a)

  end subroutine realloc4r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc5r(a,n1,n2,n3,n4,n5)
    use tools_io, only: ferror, faterr

    real*8, intent(inout), allocatable :: a(:,:,:,:,:) !< Input array, real*8, 3D
    integer, intent(in) :: n1, n2, n3, n4, n5 !< new dimension
    
    real*8, allocatable :: temp(:,:,:,:,:)
    integer :: nold(5)
    
    if (.not.allocated(a)) &
       call ferror('realloc5r','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    nold(3) = size(a,3)
    nold(4) = size(a,4)
    nold(5) = size(a,5)
    if (nold(1) == n1 .and. nold(2) == n2 .and. nold(3) == n3 .and.&
        nold(4) == n4 .and. nold(5) == n5) return
    allocate(temp(n1,n2,n3,n4,n5))
    
    temp = 0d0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)),1:min(n5,nold(5))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)),1:min(n5,nold(5)))
    call move_alloc(temp,a)

  end subroutine realloc5r

  !> Adapt the size of an allocatable 1D integer*2 array
  subroutine realloc1i2(a,nnew)
    use tools_io, only: ferror, faterr

    integer*2, intent(inout), allocatable :: a(:) !< Input array, integer, 1D
    integer, intent(in) :: nnew !< New dimension
    
    integer*2, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1i2','array not allocated',faterr)
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1i2

  !> Adapt the size of an allocatable 1D integer array
  subroutine realloc1i(a,nnew)
    use tools_io, only: ferror, faterr

    integer, intent(inout), allocatable :: a(:) !< Input array, integer, 1D
    integer, intent(in) :: nnew !< New dimension
    
    integer, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1i','array not allocated',faterr)
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1i

  !> Adapt the size of an allocatable 1D real*8 array
  subroutine realloc2i(a,n1,n2)
    use tools_io, only: ferror, faterr

    integer, intent(inout), allocatable :: a(:,:) !< Input array, integer, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    integer, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) &
       call ferror('realloc2i','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    if (nold(1) == n1 .and. nold(2) == n2) return
    allocate(temp(n1,n2))
    
    temp = 0
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2i

  !> Adapt the size of an allocatable 1D character array
  subroutine realloc1c(a,nnew)
    use tools_io, only: ferror, faterr

    character*(*), intent(inout), allocatable :: a(:) !< Input array, character, 1D
    integer, intent(in) :: nnew !< New dimension
    
    character*(len(a)), allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1c','array not allocated',faterr)
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1c

  !> Adapt the size of an allocatable 1D complex*8 array
  subroutine realloc1cmplx4(a,nnew)
    use tools_io, only: ferror, faterr

    complex*8, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    complex*8, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1cmplx4','array not allocated',faterr)
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1cmplx4

  !> Adapt the size of an allocatable 2D complex*8 array
  subroutine realloc2cmplx4(a,n1,n2)
    use tools_io, only: ferror, faterr

    complex*8, intent(inout), allocatable :: a(:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2 !< new dimension
    
    complex*8, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) &
       call ferror('realloc2cmplx4','array not allocated',faterr)
    nold(1) = size(a,1)
    nold(2) = size(a,2)
    if (nold(1) == n1 .and. nold(2) == n2) return
    allocate(temp(n1,n2))
    
    temp = cmplx(0,8)
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2cmplx4

  !> Adapt the size of an allocatable 2D complex*8 array
  subroutine realloc4cmplx4(a,n1,n2,n3,n4)
    use tools_io, only: ferror, faterr

    complex*8, intent(inout), allocatable :: a(:,:,:,:) !< Input array, real*8, 2D
    integer, intent(in) :: n1, n2, n3, n4 !< new dimension
    
    complex*8, allocatable :: temp(:,:,:,:)
    integer :: nold(4), i
    
    if (.not.allocated(a)) &
       call ferror('realloc4cmplx4','array not allocated',faterr)
    do i = 1, 4
       nold(i) = size(a,i)
    end do
    if (nold(1) == n1 .and. nold(2) == n2 .and. nold(3) == n3 .and.&
        nold(4) == n4) return
    allocate(temp(n1,n2,n3,n4))
    
    temp = cmplx(0,8)
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)),1:min(n4,nold(4)))
    call move_alloc(temp,a)

  end subroutine realloc4cmplx4

  !> Adapt the size of an allocatable 1D complex*16 array
  subroutine realloc1cmplx8(a,nnew)
    use tools_io, only: ferror, faterr

    complex*16, intent(inout), allocatable :: a(:) !< Input array, real*8, 1D
    integer, intent(in) :: nnew !< new dimension
    
    complex*16, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call ferror('realloc1cmplx8','array not allocated',faterr)
    nold = size(a)
    if (nold == nnew) return
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1cmplx8

end module types
