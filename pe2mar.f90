program pe2mar

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
  type (varying_string) :: myoptions(16)

  myoptions(1) = 'help'
  myoptions(2) = 'verbose'
  myoptions(3) = 'minimum'
  myoptions(4) = 'nx'
  myoptions(5) = 'ny'
  myoptions(6) = 'maximum'
  myoptions(7) = 'dark'
  myoptions(8) = 'darkscale'
  myoptions(9) = 'wavelength'
  myoptions(10) = 'distance'
  myoptions(11) = 'phistart'
  myoptions(12) = 'phiinc'
  myoptions(13) = 'badpixel'
  myoptions(14) = 'noise'
  myoptions(15) = 'int'
  myoptions(16) = 'maskfile'

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

  verbose = option_exists('verbose')

  data_is_real = .not. option_exists('int')

  autoscale_min = .FALSE.
  autoscale_max = .FALSE.
  fixed_min = .FALSE.
  fixed_max = .FALSE.

  maximum=2**16 - 1
  ! See if we have specified a maximum we wish to scale the data to
  if (option_exists('maximum')) then
     ! Make sure we have a value
     if (has_value('maximum')) then
        ! Get the maximum value we will scale to
        maximum = get_value('maximum')
        fixed_max = .TRUE.
     else
        autoscale_max = .TRUE.
     end if
  end if

  minimum = 0
  ! See if we have specified a minimum we wish to scale the data to
  if (option_exists('minimum')) then
     ! Make sure we have a value
     if (has_value('minimum')) then
        ! Get the maximum value we will scale to
        minimum = get_value('minimum')
        fixed_min = .TRUE.
     else
        autoscale_min = .TRUE.
     end if
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
     if (fixed_ny) then
     else
        ny = nx
     end if
  else
     if (.not. fixed_ny) then
        nx = 2048
        ny = 2048
     else
        nx = ny
     end if
  end if

  ! Allocate the space required for the data
  allocate(realdata(nx,ny),data(nx,ny),rawdata(nx*ny),overexposed(nx,ny))

  allocate(goodpixels(nx,ny))

  ! See if we have specified a badpixel image
  if (option_exists('badpixel')) then
     ! Make sure we have a value
     if (.NOT. has_value('badpixel')) then
        write(stderr,*) 'Option badpixel must specify a pgm file to read'
        call usage
        stop
     end if
     fname = get_value('badpixel')
     have_badpixel = .TRUE.

     call read(pnm,trim(fname))
     goodpixels = (as_array_2d(pnm) /= 0)
     print '(A,I0,A,I0,A)','Masking ',(nx*ny - count(goodpixels)),' out of ',nx*ny,' pixels'

  else
     goodpixels = .TRUE.
     have_badpixel = .FALSE.
  end if

  ! See if we have specified a maskfile
  if (option_exists('maskfile')) then
     ! Make sure we have a value
     if (.NOT. has_value('maskfile')) then
        write(stderr,*) 'Option maskfile must specify a file to read'
        call usage
        stop
     end if
     fname = get_value('maskfile')
     have_badpixel = .TRUE.

     call read_mask(trim(fname),goodpixels)
     ! Mask file is true for bad pixels, so invert
     goodpixels = .not. goodpixels
     print '(A,I0,A,I0,A)','Masking ',(nx*ny - count(goodpixels)),' out of ',nx*ny,' pixels'

  else
     goodpixels = .TRUE.
     have_badpixel = .FALSE.
  end if

  ! See if we have specified a dark current image
  if (option_exists('dark')) then
     ! Make sure we have a value
     if (.NOT. has_value('dark')) then
        write(stderr,*) 'Option dark must specify a file name!'
        call usage
        stop
     end if
     fname = ''
     fname = get_value('dark')
     have_dark = .TRUE.

     ! Allocate the space required for the data
     allocate(dark(nx,ny))

     dark = 0

     if (data_is_real) then

        unit = freeunit()
        open(unit=unit, file=trim(fname), form='UNFORMATTED', &
             action='READ', status='OLD', access='DIRECT', recl=myrecl) 
        
        irec = 3
        
        do i = 1, ny
           do j = 1, nx
              read(1, rec=irec) realtmp
              dark(j,i) = realtmp
              irec = irec + 1
           end do
        end do
        
        close(unit)

     else

        call open(fh, trim(fname))
        
        ! Read in header info
        call read(fh, tmp,  4, littleendian=.true.)
        call read(fh, tmp2,  4, littleendian=.true.)

        do i = 1, ny
           do j = 1, nx
              call read(fh, tmp,  2, littleendian=.true.)
              call read(fh, tmp2,  2, littleendian=.true.)
              dark(j,i) = real(ishft(tmp2,16) + tmp)
           end do
        end do

        call close(fh)

     end if

     if (option_exists('darkscale')) then
        if (.NOT. has_value('darkscale')) then
           write(stderr,*) 'Option darkscale must specify a value!'
           call usage
           stop
        end if
        ! Get the maximum value we will scale to
        darkscale = get_value('darkscale')
        dark = dark/darkscale
     end if

  else
     have_dark = .FALSE.
  end if

  if (option_exists('wavelength')) then
     if (.NOT. has_value('wavelength')) then
        write(stderr,*) 'Option wavelength must specify a value!'
        call usage
        stop
     end if
     lambda = get_value('wavelength')
  else
     lambda = 99.999
  end if

  if (option_exists('distance')) then
     if (.NOT. has_value('distance')) then
        write(stderr,*) 'Option distance must specify a value!'
        call usage
        stop
     end if
     distance = get_value('distance')
  else
     distance = 9999.99
  end if

  if (option_exists('phistart')) then
     if (.NOT. has_value('phistart')) then
        write(stderr,*) 'Option phistart must specify a value!'
        call usage
        stop
     end if
     phistart = get_value('phistart')
  else
     phistart = 9999.99
  end if

  if (option_exists('phiinc')) then
     if (.NOT. has_value('phiinc')) then
        write(stderr,*) 'Option phiinc must specify a value!'
        call usage
        stop
     end if
     phiinc = get_value('phiinc')
  else
     phiinc = 9999.99
  end if

!!$  ! See if we have specified we want to add some random noise to the background
!!$  if (option_exists('noise')) then
!!$     ! Make sure we have a value
!!$     if (.NOT. has_value('noise')) then
!!$        write(stderr,*) 'Option noise must have a value!'
!!$        call usage
!!$        stop
!!$     end if
!!$     ! Get the maximum value we will scale to
!!$     noise = get_value('noise')
!!$     add_noise = .TRUE.
!!$     call rseed(123,456)
!!$  end if

  filenum = 0
     
  do while (have_args())

     filenum = filenum + 1

     ! The files are on the command line
     fname = next_arg()

     ! The data is a tiff file, but we rely on the fact that the first 8 bytes
     ! are the standard tiff header 'II' followed by the number 42. This is 
     ! relying on little endian byte ordering. All bets are off it is big
     ! endian or this is run on a big endian machine.

     if (data_is_real) then

        unit = freeunit()
        open(unit=unit, file=trim(fname), form='UNFORMATTED', &
             action='READ', status='OLD', access='DIRECT', recl=myrecl) 
     
        data = 0
        irec=3
        do i = 1, ny
           do j = 1, nx
              read(1, rec=irec) realtmp
              realdata(j,i) = realtmp
              irec=irec+1
           end do
        end do

        close(unit)

     else

        call open(fh, trim(fname))
        
        ! Read in header info
        call read(fh, tmp,  4, littleendian=.true.)
        call read(fh, tmp2,  4, littleendian=.true.)

        do i = 1, ny
           do j = 1, nx
              call read(fh, tmp,  2, littleendian=.true.)
              call read(fh, tmp2,  2, littleendian=.true.)
              ! data(j,i) = ishft(tmp2,16) + tmp
              ! realdata(j,i) = real(tmp)
              realdata(j,i) = real(ishft(tmp2,16) + tmp)
           end do
           ! print *,i
        end do

        call close(fh)

     end if

     close(unit)

     if (have_dark) realdata = realdata - dark

     overexposed = (realdata >= dynamic_range)
     neg = count(realdata < 0.)

     print '(I0,A)',count(overexposed),' overexposed pixels'
     print '(I0,A)',neg,' negative pixels'
     print *,'Minimum value: ', minval(realdata)
     print *,'Maximum value: ', maxval(realdata)

     if (autoscale_min .or. fixed_min) then
        if (autoscale_min) minimum = minval(realdata,mask=goodpixels)
        realdata = realdata - minimum
        print *,'Removing background of ',minimum
     end if
     if (autoscale_max .or. fixed_max) then
        if (autoscale_max) maximum = maxval(realdata,mask=goodpixels)
        print *,'Scaling to a maximum of ',maximum
        realdata = realdata * (real(2**16 - 1) / maximum)
     end if

     where (realdata < 0) realdata = 0

     if (have_badpixel) where (.not. goodpixels) realdata = 0

     where (overexposed) realdata = 999999

     data = nint(realdata)

     call new(mar345, int(data), pixlength=200., pixheight=200., wavelength=lambda, &
          distance=distance, phistart=phistart+phiinc*(filenum-1), phiend=phistart+phiinc*(filenum),status=error)

     ! Use the last input filename as a template for the output filename.
     ! Search for a ".***" suffix .. if we find one than 
     ! replace it with ".mar****"
     i = index(fname,".",back=.true.)
     if (i == 0) then
        ! Couldn't find a matching suffix, so look for the end of the
        ! filename string in the character variable
        i = verify(fname," ",back=.TRUE.) + 1
     end if
     outname = fname
     write(outname(i:),'(A,I0)') ".mar",nx
     
     if (verbose) print *,'Wrote ',trim(outname)

     ! Write the pgm file
     call write(mar345, trim(outname))

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

!!$  real function rand_number
!!$
!!$    real :: x(1)
!!$
!!$    call rannum(x, 1)
!!$
!!$    rand_number = x(1)
!!$
!!$  end function rand_number

end program pe2mar
