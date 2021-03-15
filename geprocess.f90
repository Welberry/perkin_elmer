program geprocess

  use cmdline_arguments
  use string_functions, only: join, real
  use file_functions, only: stderr
  use pnm_class
  use precision, only: i8_kind, rd_kind 
  use binary_io, only: read, open, close, binary_filehandle
  use iso_varying_string
  use mar_class
  use statistics, only: mean

  implicit none

  integer, parameter :: dynamic_range = 16353

  integer :: i, j, k, error, nx, ny, unit, depth, tmp, tmp2, nframes
  integer :: frame, filenum, radius, numpixels

  character(len=2000) :: fname, outname

  logical :: verbose = .FALSE., autoscale, fixed_max, fixed_nx, fixed_ny
  logical :: fixed_depth, have_dark, have_gainmap, have_badpixel, average, smoothing, smoothing_dark

  type (pnm_object) :: pnm

  real :: maximum, darkscale, lambda, distance, phistart, phiinc

  integer(i8_kind) :: minimum

  integer, allocatable :: rawdata(:)
  logical, allocatable :: goodpixels(:,:)
  logical, allocatable :: frames(:)

  integer(i8_kind), dimension(:,:), allocatable :: data, dark, newdata
  real(rd_kind), dimension(:,:), allocatable :: gainmap
  logical, dimension(:,:), allocatable :: overexposed
  
  type (binary_filehandle) :: fh

  type (mar_object) :: mar345

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  type (varying_string) :: myoptions(19)

  myoptions(1) = 'help'
  myoptions(2) = 'average'
  myoptions(3) = 'verbose'
  myoptions(4) = 'depth'
  myoptions(5) = 'nx'
  myoptions(6) = 'ny'
  myoptions(7) = 'norm'
  myoptions(8) = 'dark'
  myoptions(9) = 'darkscale'
  myoptions(10) = 'nframes'
  myoptions(11) = 'wavelength'
  myoptions(12) = 'distance'
  myoptions(13) = 'phistart'
  myoptions(14) = 'phiinc'
  myoptions(15) = 'gainmap'
  myoptions(16) = 'badpixel'
  myoptions(17) = 'smooth'
  myoptions(18) = 'smoothdark'
  myoptions(19) = 'frame'

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
  average = option_exists('average')

  autoscale = .FALSE.
  fixed_max = .FALSE.

  ! See if we have specified a maximum we wish to scale the data to
  if (option_exists('norm')) then
     ! Make sure we have a value
     if (has_value('norm')) then
        ! Get the maximum value we will scale to
        maximum = get_value('norm')
        fixed_max = .TRUE.
        autoscale = .FALSE.
     else
        autoscale = .TRUE.
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

  ! Allocate the space required for the data
  allocate(data(nx,ny),rawdata(nx*ny),overexposed(nx,ny))

  ! See if we have specified a depth (bits per pixel)
  if (option_exists('depth')) then
     ! Make sure we have a value
     if (.NOT. has_value('depth')) then
        write(stderr,*) 'Option depth must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     depth = get_value('depth')
     fixed_depth = .TRUE.
  else
     depth = 2
     fixed_depth = .FALSE.
  end if

  ! See if we have the number of frames
  if (option_exists('nframes')) then
     ! Make sure we have a value
     if (.NOT. has_value('nframes')) then
        write(stderr,*) 'Option nframes must have a value!'
        call usage
        stop
     end if
     nframes = get_value('nframes')
  else
     nframes = 1
  end if

  allocate(frames(nframes))

  if (option_exists('frame')) then
     ! Make sure we have a value
     if (.NOT. has_value('frame')) then
        write(stderr,*) 'Option frame must have a value!'
        call usage
        stop
     end if
     frame = get_value('frame')
     frames = .FALSE.
     frames(frame) = .TRUE.
  else
     frames = .TRUE.
  end if

  if (option_exists('smooth')) then
     ! Make sure we have a value
     if (.NOT. has_value('smooth')) then
        write(stderr,*) 'Option smooth must have a value!'
        call usage
        stop
     end if
     radius = get_value('smooth')
     smoothing = .true.
  else
     smoothing = .false.
  end if

  if (option_exists('smoothdark')) then
     ! Make sure we have a value
     if (.NOT. has_value('smoothdark')) then
        write(stderr,*) 'Option smoothdark must have a value!'
        call usage
        stop
     end if
     radius = get_value('smoothdark')
     smoothing_dark = .true.
  else
     smoothing_dark = .false.
  end if

  if (fixed_nx) then
     if (fixed_ny) then
     else
        ny = nx
     end if
  else
     if (.not. fixed_ny) then
        write(stderr,*) 'Must specify at least nx or ny as command line options'
        call usage
        stop
     else
        nx = ny
     end if
  end if

  ! See if we have specified a gainmap image
  if (option_exists('gainmap')) then
     ! Make sure we have a value
     if (.NOT. has_value('gainmap')) then
        write(stderr,*) 'Option gainmap must specify a file name!'
        call usage
        stop
     end if
     fname = get_value('gainmap')
     have_gainmap = .TRUE.

     call open(fh, trim(fname))

     ! Read in the header and discard
     call read_header(fh)

     rawdata = 0
     
     call read(fh, rawdata,  4, littleendian=.true.)

     call close(fh)

     allocate(gainmap(nx,ny))

     gainmap = reshape(rawdata / mean(rawdata), (/nx, ny/))

  else
     have_gainmap = .FALSE.
  end if

  allocate(goodpixels(nx,ny))

  ! See if we have specified a badpixel image
  if (option_exists('badpixel')) then
     ! Make sure we have a value
     if (.NOT. has_value('badpixel')) then
        write(stderr,*) 'Option badpixel must specify a file name!'
        call usage
        stop
     end if
     fname = get_value('badpixel')
     have_badpixel = .TRUE.

     call open(fh, trim(fname))

     ! Read in the header and discard
     call read_header(fh)

     rawdata = 0
     
     call read(fh, rawdata,  2, littleendian=.true.)

     call close(fh)

     goodpixels = reshape(rawdata == 0, (/nx, ny/))
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
     call open(fh, trim(fname))

     ! Read in the header and discard
     call read_header(fh)

     ! Allocate the space required for the data
     allocate(dark(nx,ny))

     dark = 0
     
     do k = 1, nframes
        print '(A,I0)','Reading dark frame ',k
        do i = 1, ny
           do j = 1, nx
              call read(fh, tmp,  2, littleendian=.true.)
              dark(j,i) = dark (j,i) + tmp
              ! call read(fh, tmp2, 2, littleendian=.true.)
              ! dark(j,i) = ishft(tmp2,16) + tmp
           end do
        end do
     end do

     ! dark = nint(real(dark,rd_kind)/real(nframes,rd_kind))
     call close(fh)
     if (option_exists('darkscale')) then
        if (.NOT. has_value('darkscale')) then
           write(stderr,*) 'Option darkscale must specify a value!'
           call usage
           stop
        end if
        ! Get the maximum value we will scale to
        darkscale = get_value('darkscale')
        dark = nint(real(dark,rd_kind)*darkscale)
     end if
     if (smoothing_dark) then
        ! DON'T DO IT! THIS IS A BAD IDEA! The signal is so weak compared to the
        ! dark current that any smoothing makes a real mess of the output
        allocate(newdata(nx,ny))
        newdata = 0
        ! Look, we're cheating here alright?! We don't do the outer edge of the
        ! image and we're not that bothered as it is usually outside the detector
        ! area anyway. So there.
        do i = radius, ny-radius
           do j = radius, nx-radius
              if (.not. goodpixels(j,i)) cycle
              numpixels = count(goodpixels(j-radius:j+radius,i-radius:i+radius))
              if (numpixels == 0) cycle
              newdata(j,i) = sum(dark(j-radius:j+radius,i-radius:i+radius),mask=goodpixels(j-radius:j+radius,i-radius:i+radius))/numpixels
              ! Could conceivably make a new goodpixels mask by replacing some 
              ! missing pixels with the average of surrounding pixels. At this
              ! stage we won't "add" data
           end do
        end do
        dark = newdata
        deallocate(newdata)
     end if
     print *,maxval(dark),count(frames),nframes
     if (count(frames) /= nframes) dark = real(dark)*real(count(frames))/real(nframes)
     print *,maxval(dark),count(frames),nframes

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

  filenum = 0
     
  do while (have_args())

     filenum = filenum + 1

     ! The files are on the command line
     fname = next_arg()

     call open(fh, trim(fname))

     ! Read in the header and discard
     call read_header(fh)

     data = 0
     overexposed = .false.

     do k = 1, nframes
        ! frame = 0
        if (frames(k)) then
           print '(A,I0)','Reading frame ',k
           do i = 1, ny
              do j = 1, nx
                 call read(fh, tmp, depth, littleendian=.true.)
                 ! tmp = max(0,tmp - dark(j,i))
                 data(j,i) = data(j,i) + tmp
                 ! frame(j,i) = frame(j,i) + tmp
                 ! call read(fh, tmp,  2, littleendian=.true.)
                 ! call read(fh, tmp2, 2, littleendian=.true.)
                 ! data(j,i) = ishft(tmp2,16) + tmp
                 if (tmp >= dynamic_range) overexposed(j,i) = .TRUE.
              end do
           end do
        else
           print '(A,I0)','Skipping frame ',k
           call read(fh, rawdata,  depth, littleendian=.true.)
        end if
        ! if (have_dark) frame = frame - dark
        ! data = data + frame
     end do

     call close(fh)

     print '(I0,A)',count(overexposed),' overexposed pixels'

     if (have_dark) data = data - dark

     if (smoothing) then
        allocate(newdata(nx,ny))
        newdata = 0
        ! Look, we're cheating here alright?! We don't do the outer edge of the
        ! image and we're not that bothered as it is usually outside the detector
        ! area anyway. So there.
        do i = radius, ny-radius
           ! miny=max(1,i-1)
           ! maxy=max(y,i+1)
           do j = radius, nx-radius
              if (.not. goodpixels(j,i)) cycle
              numpixels = count(goodpixels(j-radius:j+radius,i-radius:i+radius))
              if (numpixels == 0) cycle
              newdata(j,i) = sum(data(j-radius:j+radius,i-radius:i+radius),mask=goodpixels(j-radius:j+radius,i-radius:i+radius))/numpixels
              ! Could conceivably make a new goodpixels mask by replacing some 
              ! missing pixels with the average of surrounding pixels. At this
              ! stage we won't "add" data
           end do
        end do
        data = newdata
        deallocate(newdata)
     end if

     if (autoscale) then
        minimum = minval(data,mask=goodpixels)
        data = data - minimum
        maximum = maxval(data,mask=goodpixels)
        print *,'scaling from ',minimum,' to ',minimum+maximum
        data = nint(real(data, rd_kind) * (real(2**16 - 1) / maximum) )
     else if (fixed_max) then
        print *,'scaling to ',maximum
        data = nint(real(data, rd_kind) * (real(2**16 - 1) / maximum) )
     end if

     where (data < 0) data = 0

     if (have_gainmap) data = nint(data*gainmap)

     if (average) data = nint(real(data, rd_kind) / real(nframes, rd_kind))

     if (have_badpixel) where (.not. goodpixels) data = 0

     where (overexposed) data = 999999

     call new(mar345, int(data), pixlength=200., pixheight=200., wavelength=lambda, distance=distance, phistart=phistart+phiinc*(filenum-1), phiend=phistart+phiinc*(filenum))

     ! Use the last input filename as a template for the output filename.
     ! Search for a ".***" suffix .. if we find one than 
     ! replace it with ".mar****"
     i = index(fname,".")
     if (i == 0) then
        ! Couldn't find a matching suffix, so look for the end of the
        ! filename string in the character variable
        i = verify(fname," ",back=.TRUE.) + 1
     end if
     outname = fname
     write(outname(i:),'(A,I0)') ".mar",nx
     
     if (verbose) print *,'Wrote ',trim(outname)

     ! Write the pgm file
     ! call write(pnm, trim(outname))
     call write(mar345, trim(outname))

  end do

  deallocate(data)

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Convert kuplot NIPL formatted files to pgm format'
    write(stderr,*)
    write(stderr,*) 'Usage: ni2pgm [--help] [--verbose] [--noscale] nifile(s)'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*) '  --verbose - verbose output'
    write(stderr,*) '  --noscale - do not rescale data to 65535'
    write(stderr,*)

  end subroutine usage

  subroutine read_header(fh)

    ! These GE files have an 8192 byte header, which we read (and discard)
    ! by reading a 2048 block of 4 byte integers 

    type (binary_filehandle), intent(inout) :: fh

    integer :: header(2048)

    call read(fh, header,  4)
           
  end subroutine read_header

end program geprocess
