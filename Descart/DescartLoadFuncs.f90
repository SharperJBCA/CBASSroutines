! type test
!    integer value
! end type test


module mymod
  implicit none
  private
  public :: loadOffsetsFromFilename
  public :: ds_get_lun
  public :: loadOffsetsFromFileNum
  public :: countOffsetsFromFile
contains

  subroutine countOffsetsFromFile(filename, count, length)
    implicit none
     
    character(*), intent(in) :: filename
    integer :: unit
    integer, intent(out):: count,length

    unit=ds_get_lun()
    open(unit=unit,file=filename,form='unformatted',action='read', status='old')
    read(unit) length
    read(unit) count
    close(unit)
  end subroutine countOffsetsFromFile
  

  subroutine loadOffsetsFromFilename(filename, values) !result(values)
    implicit none
     
    !type(ds_offsets) :: offsets
    character(*) :: filename
    integer :: unit
    real*8, dimension(:), intent(inout) :: values

    unit=ds_get_lun()
    open(unit=unit,file=filename,form='unformatted',action='read', status='old')
    call loadOffsetsFromFileNum(unit, values)
    close(unit)
  end subroutine loadOffsetsFromFilename

  subroutine loadOffsetsFromFileNum(unit, values) !result(values)
    implicit none

    integer unit 
    integer :: length, na,n_az
    real*8, dimension(:), intent(inout) :: values
    logical :: azimuth_flag
    real*8,allocatable,dimension(:) :: amplitudes  !size n_az
    real(8),allocatable,dimension(:) :: values_read
    integer :: i

    read(unit) length
    read(unit) na
    allocate(values_read(na))
    read(unit) values_read
    read(unit) azimuth_flag
    if (azimuth_flag) then
       read(unit) n_az
       allocate(amplitudes(n_az))
       read(unit) amplitudes
    endif

    do i=1, size(values_read)
       values(i) = values_read(i)
    enddo
    deallocate(values_read)
  end subroutine loadOffsetsFromFileNum


  function ds_get_lun() result(f)
    implicit none
    !Get an unused logical unit number for a file.
    integer :: f
    logical :: inUse
    f=50
    do
       inquire(UNIT=f, opened=inUse)
       if (.not. inUse) exit
       f=f+1
    enddo
  end function ds_get_lun

end module mymod

subroutine wrapper(filename,values, n)
  use mymod
  implicit none
  
  integer :: n
  character(*) :: filename
  real*8, intent(inout) :: values(n)
!f2py intent(in) filename
!f2py intent(in) n

  call loadOffsetsFromFilename(filename, values)
end subroutine wrapper

subroutine countvals(filename,count, length)
  use mymod
  implicit none
  
  character(*) :: filename
  integer, intent(out) :: count, length
!f2py intent(in) filename
!f2py intent(out) count ,length

  call countOffsetsFromFile(filename , count, length)
end subroutine countvals
