program darkcleaner

  ! This program is designed to clean up dark current files from the
  ! perkin elmer detector

  use cmdline_arguments
  use string_functions, only: join, real
  use file_functions, only: stderr, freeunit
  use pnm_class
  use precision, only: i8_kind, rs_kind, rd_kind, i4_kind
  use binary_io, only: read, open, close, binary_filehandle
  use iso_varying_string
  use mar_class
  use statistics, only: median
  use sort_functions, only: sort
  use rannum_module, only: rseed, rannum

  implicit none

  integer, parameter :: dynamic_range = 999999 !16353
  integer, parameter :: subwidth = 128 ! This is the width of the individual 
                                       ! elements which make up the detector

  ! Note that we rely on recl being in bytes, if it is in words this will
  ! need to be changed to myrecl = 1
  integer, parameter :: myrecl = 4

  integer :: i, j, k, error, nx, ny, unit, unit2, tmp, tmp2, nframes, ios, ubnd, lbnd
  integer :: frame, filenum, radius, numpixels, neg, irec, minimum, start, finish

  real(kind=4) :: realtmp, header(2), buffer
  character(len=4) :: buf

  character(len=2000) :: fname, outname

  logical :: verbose = .FALSE., autoscale_min, autoscale_max, fixed_min, fixed_max
  logical :: fixed_nx, fixed_ny, have_dark, have_badpixel, usedark, add_noise=.FALSE.

  type (pnm_object) :: pnm

  real :: darkscale, threshold, noise

  integer, allocatable :: rawdata(:)
  logical, allocatable :: goodpixels(:,:)
  logical, allocatable :: frames(:)

  integer(i8_kind), dimension(:,:), allocatable :: data, newdata
  real(rd_kind), dimension(:,:), allocatable :: realdata, dark, medians
  logical, dimension(:,:), allocatable :: overexposed

  real(rd_kind) :: databuf(subwidth), databufnorm(subwidth), medianbuf, darkbuf(subwidth)
  
  type (binary_filehandle) :: fh

  type (mar_object) :: mar345

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  type (varying_string) :: myoptions(7)

  myoptions(1) = 'help'
  myoptions(2) = 'verbose'
  myoptions(3) = 'nx'
  myoptions(4) = 'ny'
  myoptions(5) = 'threshold'
  myoptions(6) = 'usedark'
  myoptions(7) = 'noise'

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

  if (num_args() < 2) then
     write(stderr,*) 'ERROR! Must supply a reference dark current file as first argument'
     write(stderr,*) '       and 1 or more dark current file(s) as subsequent arguments'
     call usage
     STOP
  end if

  verbose = option_exists('verbose')
  usedark = option_exists('usedark')

  ! See if we have specified a width
  if (option_exists('nx')) then
     ! Make sure we have a value
     if (.NOT. has_value('nx')) then
        write(stderr,*) 'Option nx must have a value!'
        call usage
        stop
     end if
     ! Get the number of pixels in the horizontal direction
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
     ! Get the number of pixels in the vertical direction
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

  ! See if we have specified a threshold for detecting aberrant peaks
  if (option_exists('threshold')) then
     ! Make sure we have a value
     if (.NOT. has_value('threshold')) then
        write(stderr,*) 'Option threshold must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     threshold = get_value('threshold')
  else
     threshold = 0.05
  end if

  ! See if we have specified we want to add some random noise to the background
  if (option_exists('noise')) then
     ! Make sure we have a value
     if (.NOT. has_value('noise')) then
        write(stderr,*) 'Option noise must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     noise = get_value('noise')
     add_noise = .TRUE.
     call rseed(123,456)
  end if

  ! Allocate the space required for the data
  allocate(realdata(nx,ny),data(nx,ny),rawdata(nx*ny),overexposed(nx,ny))

  ! See if we have specified a dark current image
  fname = next_arg()
  have_dark = .TRUE.

  unit = freeunit()
  open(unit=unit, file=trim(fname), form='UNFORMATTED', &
       action='READ', status='OLD', access='DIRECT', recl=myrecl) 
  
  ! Allocate the space required for the data
  allocate(dark(nx,ny))

  dark = 0
  irec = 3
     
  do i = 1, ny
     do j = 1, nx
        read(1, rec=irec) realtmp
        dark(j,i) = realtmp
        irec = irec + 1
     end do
  end do
  
  close(unit)
  
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

  allocate(medians(nx/subwidth,ny))

  ! unit = freeunit()
  ! open(unit=unit, file='medians.dat', form='FORMATTED', &
  !      action='WRITE', status='UNKNOWN')

  do i = 1, size(medians,2)
     do j = 1, size(medians,1)
        start = (j-1)*subwidth + 1
        finish = start + (subwidth-1)
        medians(j,i) = median(pack(dark(start:finish,i), mask=dark(start:finish,i)>0))
     end do
     ! write(unit,'(16F10.2)') medians(:,i)
  end do

  ! close(unit)

  ! Make a new pgm object -- assumes a top left origin in accordance
  ! with the mar format specification
  ! call new(pnm, nint(medians,i4_kind), origin='tl')
  ! pnm = nint(65535.*(medians/maxval(medians)),i4_kind)
  ! call write(pnm, 'medians.pgm')

  ! stop

  filenum = 0
     
  do while (have_args())

     filenum = filenum + 1

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

     read(unit, rec=1) header(1)
     read(unit, rec=2) header(2)

     irec=3

     do i = 1, ny
        do j = 1, nx
           read(unit, rec=irec) realtmp
           realdata(j,i) = realtmp
           irec=irec+1
        end do
     end do

     close(unit)

     if (usedark) then
        ! We will exclusively use the data from long time dark file, but
        ! scale it with the current dark file
        do i = 1, size(medians,2)
           INNER: do j = 1, size(medians,1)
              start = (j-1)*subwidth + 1
              finish = start + (subwidth-1)

              databuf = pack(realdata(start:finish,i), mask=.TRUE.)

              ! We won't do this is if more than half our pixels are missing
              if (count(databuf == 0) > subwidth/2) cycle INNER

              darkbuf = pack(dark(start:finish,i), mask=.TRUE.)

              ! We won't do this is if more than half of our long dark 
              ! pixels 
              if (count(darkbuf == 0) > subwidth/2) cycle INNER

              darkbuf = darkbuf/medians(j,i)

              ! Sort the data instead of calculating a median
              call sort(databuf)

              lbnd = 1
              do
                 if (databuf(lbnd) /= 0) exit
                 lbnd = lbnd + 1
              end do
              ubnd = subwidth
              do
                 medianbuf = databuf((ubnd-lbnd)/2 + lbnd)
                 if ((databuf(ubnd)-medianbuf)/medianbuf < threshold) exit
                 ubnd = ubnd - 1
                 ! Abandon this row if we have too few pixels left to calculate a meaningful median
                 if ((ubnd - lbnd) < subwidth/4) cycle INNER
              end do

              ! print *,i,j,lbnd,ubnd,medianbuf

              databuf = darkbuf*medianbuf
              realdata(start:finish,i) = databuf
           end do INNER
        end do
     else
        ! Here we'll use as much of the current dark file as possible
        do i = 1, size(medians,2)
           do j = 1, size(medians,1)
              start = (j-1)*subwidth + 1
              finish = start + (subwidth-1)
              ! print *,'j: ',j, start, finish
              databuf = pack(realdata(start:finish,i), mask=.TRUE.)
              if (count(databuf == 0) > 100) cycle
              ! print *,databuf
              darkbuf = pack(dark(start:finish,i), mask=.TRUE.)
              darkbuf = darkbuf/medians(j,i)
              if (count(darkbuf == 0) > 100) cycle
              if (all(darkbuf == 0)) cycle
              ! print *,darkbuf
              ! medians(j,i) = median(dark(j:j+(subwidth-1),i))
              medianbuf = median(pack(databuf, mask=(databuf>0)))
              databufnorm = databuf/medianbuf
              ! print *,'j: ',j, start, finish, medianbuf, count(abs(databufnorm - darkbuf) > threshold)
              ! print *,j,medianbuf
              ! print *,(abs(databufnorm - darkbuf) > threshold)
              ! print '(subwidth.2)',databufnorm - darkbuf
              where (abs(databufnorm - darkbuf) > threshold) databuf = darkbuf*medianbuf
              ! print *,databuf
              if (add_noise) databuf = databuf*(1 + rand_number()*noise)
              realdata(start:finish,i) = databuf
           end do
        end do
     end if
        
     ! Search for a ".tif" suffix .. if we find one than 
     ! replace it with "_clean.tif"
     i = index(fname,".tif")
     if (i == 0) then
        ! Couldn't find a matching suffix, so look for the end of the
        ! filename string in the character variable
        i = verify(fname," ",back=.TRUE.) + 1
     end if
     outname = fname
     write(outname(i:),'(A)') "_clean.tif"
     
     unit = freeunit()
     open(unit=unit, file=trim(outname), form='UNFORMATTED', &
         action='WRITE', status='UNKNOWN', access='DIRECT', recl=myrecl) 
     
     unit2 = freeunit()
     open(unit=unit2, file=trim(fname), form='UNFORMATTED', &
         action='READ', status='OLD', access='DIRECT', recl=myrecl) 
     
     read(unit2, rec=1) buffer
     write(unit, rec=1) buffer

     read(unit2, rec=2) buffer
     write(unit, rec=2) buffer

     irec=2

     do i = 1, ny
        do j = 1, nx
           irec = irec+1
           write(unit, rec=irec) real(realdata(j,i),rs_kind)
        end do
     end do

     do
        irec = irec+1
        ! read(unit2, rec=irec, end=1) buffer
        read(unit2, rec=irec, iostat=ios) buffer
        if (ios /= 0) exit
        write(unit, rec=irec) buffer
     end do

     close(unit)
     close(unit2)

     if (verbose) print *,'Wrote ',trim(outname)

  end do

  deallocate(data)

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Convert Perkin Elmer TIFF file to mar format'
    write(stderr,*)
    write(stderr,*) 'Usage: ni2pgm [--help] [--verbose] [--noscale] nifile(s)'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*) '  --verbose - verbose output'
    write(stderr,*) '  --noscale - do not rescale data to 65535'
    write(stderr,*)

  end subroutine usage

  subroutine read_header(fh)

    ! These Perkin Elmer files have an 122 byte header, which we read (and 
    ! discard) by reading a 61 block of 2 byte integers 

    type (binary_filehandle), intent(inout) :: fh

    integer :: header(61), magic, ifdoffset
    character(len=2) :: buf
    logical :: littleendian

    ! call read(fh, header,  2)
    call read(fh, buf)
    print *,buf

    if (buf == 'II') then
       littleendian=.true.
    else if (buf == 'MM') then
       littleendian=.false.
    else
       write(stderr,*) 'This is not a TIFF file! Bad endian indicator: ',buf
    end if
    call read(fh, magic, 2, littleendian=littleendian)
    print *,magic
    if (magic /= 42) then
       write(stderr,*) 'This is not a TIFF file! Bad endian magic number: ',magic
    end if
    call read(fh, ifdoffset, 4, littleendian=littleendian)
    print '(Z8)',ifdoffset
    ! call read(fh, magic, 2, littleendian=littleendian)
    ! print *,magic
           
  end subroutine read_header

  real function rand_number

    real :: x(1)

    call rannum(x, 1)

    rand_number = x(1)

  end function rand_number

end program darkcleaner
