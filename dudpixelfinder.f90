program dudpixelfinder

  !! Find dud pixels in a series of perkin elmer frames by finding pixels 
  !! whose values do not change

  use pnm_class
  use iso_varying_string
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use string_functions, only: join, real
  use file_functions, only: stderr, stdout, open, freeunit
  use sort_functions
  use variable_array

  implicit none

  ! Note that we rely on recl being in bytes, if it is in words this will
  ! need to be changed to myrecl = 1
  integer, parameter :: myrecl = 4

  integer :: nx, ny, error, unit, i, j, irec

  character(len=500) :: fname, scalefile

  real, allocatable :: realdata(:,:), badpixeldata(:,:)
  integer, allocatable :: badpixelcount(:,:)

  type (pnm_object) :: pnm

  real :: tolerance
  real(kind=4) :: realtmp

  logical :: verbose, quiet, have_tolerance, firstimethroughloop = .TRUE.
  logical :: fixed_nx, fixed_ny

  type (varying_string) :: myoptions(5)

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'verbose'
  myoptions(3) = 'tolerance'
  myoptions(4) = 'nx'
  myoptions(5) = 'ny'

  ! This call parses the command line arguments for command line options
  call get_options(myoptions, error)

  ! Check we weren't passed duff options -- spit the dummy if we were
  if (error > 0) then
     write(stderr,*) 'ERROR! Unknown options: ',join(bad_options()," ")
     call usage
     STOP
  end if

  ! Check if we just want to print the usage
  if (option_exists('help')) then
     call usage
     STOP
  end if

  if (num_args() < 1) then
     write(stderr,*) 'ERROR! Must supply perkin elmer data file(s) as command-line arguments'
     call usage
     STOP
  end if

  verbose = option_exists('verbose')

  if (option_exists('tolerance')) then
     ! Make sure we have a value
     if (.NOT. has_value('tolerance')) then
        write(stderr,*) 'Option tolerance must specify a value!'
        call usage
        stop
     end if
     tolerance = get_value('tolerance')
     have_tolerance = .TRUE.
  else
     tolerance = 0.01
     have_tolerance = .FALSE.
  end if

  ! See if we have specified a width
  if (option_exists('nx')) then
     ! Make sure we have a value
     if (.NOT. has_value('nx')) then
        write(stderr,*) 'Option nx must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     nx = get_value('nx')
     fixed_nx = .TRUE.
  else
     fixed_nx = .FALSE.
  end if

  ! See if we have specified a height
  if (option_exists('ny')) then
     ! Make sure we have a value
     if (.NOT. has_value('ny')) then
        write(stderr,*) 'Option ny must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     ny = get_value('ny')
     fixed_ny = .TRUE.
  else
     fixed_ny = .FALSE.
  end if

  if (fixed_nx) then
     if (.not. fixed_ny) then
        ny = nx
     end if
  else
     if (fixed_ny) then
        nx = ny
     else
        nx = 2048
        ny = 2048
     end if
  end if

  ! Allocate the space required for the data
  allocate(realdata(nx,ny),badpixeldata(nx,ny),badpixelcount(nx,ny))

  do while (have_args())

     ! The mar files are on the command line
     fname = next_arg()

     if (verbose) write(stdout,'(A,A,A)') 'Opening ',trim(fname),' ....'
     ! write(stdout,'(A,A)') trim(fname)

     ! The data is a tiff file, but we rely on the fact that the first 8 bytes
     ! are the standard tiff header 'II' followed by the number 42. This is 
     ! relying on little endian byte ordering. All bets are off it is big
     ! endian or this is run on a big endian machine.

     unit = freeunit()
     open(unit=unit, file=trim(fname), form='UNFORMATTED', &
         action='READ', status='OLD', access='DIRECT', recl=myrecl) 
     
     realdata = 0

     irec=3

     do i = 1, ny
        do j = 1, nx
           read(1, rec=irec) realtmp
           realdata(j,i) = realtmp
           irec=irec+1
        end do
     end do

     close(unit)


     if (firstimethroughloop) then
        firstimethroughloop = .FALSE.
        ! Initialise background to a very large number
        badpixelcount = 0
        badpixeldata = realdata
        cycle
     end if

     where (abs(realdata - badpixeldata)/realdata < tolerance) badpixelcount = badpixelcount + 1

     print *,minval(badpixelcount),maxval(badpixelcount)

  end do

  pnm = badpixelcount
  call write(pnm,'badpixels.pgm')

  deallocate(realdata, badpixeldata, badpixelcount)

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Usage: marbackground --help marfile(s)'
    write(stderr,*)
    write(stderr,*) '  --help     - print this message.'
    write(stderr,*) '  --verbose  - print more information'
    write(stderr,*)

  end subroutine usage

  real function get_scale(file) result(scale)

    ! Interface variables
    character(len=*) :: file
    
    ! Local variables
    character(len=1000) :: buffer
    logical, save :: opened_ok = .FALSE.
    logical :: finished
    integer, save :: unit

    if (.not. opened_ok) unit = open(trim(file),opened_ok,status='old')
    
    scale = 1.0

    if (opened_ok) then
       do
          ! Wrap this in a do loop in case we come to the end
          ! of the scale file, then we will rewind the file and
          ! traverse the loop again
          read(unit,'(F)',end=1) scale
          exit
          ! We only get here if we are at the end of the file
1         rewind(unit)
       end do
    end if

  end function get_scale

end program dudpixelfinder
