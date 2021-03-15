program raw2real

  use cmdline_arguments
  use string_functions, only: join, real
  use file_functions, only: stderr, freeunit
  use pnm_class
  use precision, only: i8_kind, rd_kind 
  use binary_io, only: read, open, close, binary_filehandle
  use iso_varying_string
  use mar_class, only: new, mar_object, write
  use statistics, only: mean
  use mask_functions
!  use rannum_module, only: rseed, rannum

  implicit none

  integer, parameter :: dynamic_range = 999999 !16353

  ! Note that we rely on recl being in bytes, if it is in words this will
  ! need to be changed to myrecl = 1
  integer, parameter :: myrecl = 4

  integer :: i, j, k, error, nx, ny, unit, tmp, tmp2
  integer :: filenum, radius, numpixels, neg, irec, minimum

  real(kind=4) :: realtmp
  character(len=4) :: buf

  character(len=2000) :: fname, outname

  logical :: verbose = .FALSE., autoscale_min, autoscale_max, fixed_min, fixed_max
  logical :: fixed_nx, fixed_ny, have_dark, have_badpixel, add_noise, data_is_real

  type (pnm_object) :: pnm

  real :: maximum, darkscale, lambda, distance, phistart, phiinc, noise

  integer, allocatable :: rawdata(:)
  logical, allocatable :: goodpixels(:,:)

  integer(i8_kind), dimension(:,:), allocatable :: data, newdata
  real(rd_kind), dimension(:,:), allocatable :: realdata, dark
  logical, dimension(:,:), allocatable :: overexposed
  
  type (binary_filehandle) :: fh

  type (mar_object) :: mar345

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  type (varying_string) :: myoptions(1)

  myoptions(1) = 'help'

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
     write(stderr,*) 'ERROR! Must supply raw file(s) as command-line arguments'
     call usage
     STOP
  end if
     
  do while (have_args())

     ! The files are on the command line
     fname = next_arg()

     ! The data is a tiff file, but we rely on the fact that the first 8 bytes
     ! are the standard tiff header 'II' followed by the number 42. This is 
     ! relying on little endian byte ordering. All bets are off it is big
     ! endian or this is run on a big endian machine.

     unit = freeunit()
     open(unit=unit, file=trim(fname), form='UNFORMATTED', &
          action='READ', status='OLD', access='DIRECT', recl=myrecl) 
     
     data = 0
     irec=1
     do i = 1, ny
        do j = 1, nx
           read(1, rec=irec) realtmp
           realdata(j,i) = realtmp
           irec=irec+1
        end do
     end do

     close(unit)
     
  end do

  deallocate(data)

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Convert Perkin Elmer TIFF file to mar format'
    write(stderr,*)
    write(stderr,*) 'Usage: pe2mar [OPTIONS] pefile(s)'
    write(stderr,*)
    write(stderr,*) ' OPTIONS:'
    write(stderr,*) '  --help           - print this message'
    write(stderr,*) '  --verbose        - verbose output'
    write(stderr,*) '  --int            - raw pe file stored as 32-bit integers'
    write(stderr,*) '  --minimum=val    - minimum to add to the data'
    write(stderr,*) '  --maximum=val    - maximum to scale the data to'
    write(stderr,*) '  --nx=val         - number of pixels in horizontal (def = 2048)'
    write(stderr,*) '  --ny=val         - number of pixels in vertical (def = 2048)'
    write(stderr,*) '  --dark=file      - dark current file to subtract'
    write(stderr,*) '  --badpixel=file  - mask file containing bad pixels to be ignored'
    write(stderr,*) '  --darkscale=val  - divide dark current by this value before subtracting'
    write(stderr,*) '  --wavelength=val - add this to mar file meta data'
    write(stderr,*) '  --distance=val   - add this to mar file meta data'
    write(stderr,*) '  --phistart=val   - add this to mar file meta data'
    write(stderr,*) '  --phiinc=val     - add this to mar file meta data'
    write(stderr,*)

  end subroutine usage

end program pe2mar
