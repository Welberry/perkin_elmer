program pestats

  !! Reports some stats on a Perkin Elmer tif file

  use iso_varying_string
  use cmdline_arguments
  use string_functions, only: join, real
  use file_functions, only: stderr, stdout, freeunit, open
  use sort_functions
  use variable_array
  use precision, only: i8_kind, rd_kind 

  implicit none

  !! $Log: $

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id:$"

  integer, parameter :: no_data = -9999

  integer :: nx, ny, npixels, nunique, error, ivalue, i, maskvalue

  character(len=500) :: fname, buffer, fmt, backgroundfile, scalefile

  real(rd_kind), allocatable :: data(:,:), data_1d(:)
  real(rd_kind), pointer     :: unique_values(:)
  integer, allocatable :: background(:,:)
  logical, allocatable :: datamask(:,:)

  real :: scale

  logical :: domin = .FALSE., domax = .FALSE., domean = .FALSE., domedian = .FALSE.
  logical :: have_scale = .FALSE., have_background = .FALSE., data_is_real

  type (varying_string) :: myoptions(12)

  logical :: verbose, quiet

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'min'
  myoptions(3) = 'max'
  myoptions(4) = 'mean'
  myoptions(5) = 'median'
  myoptions(6) = 'maskvalue'
  myoptions(7) = 'verbose'
  myoptions(8) = 'quiet'
  myoptions(9) = 'int'
  myoptions(10) = 'scale'
  myoptions(11) = 'nx'
  myoptions(12) = 'ny'

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

  data_is_real = .not. option_exists('int')

  if (option_exists('nx')) then
     if (.NOT. has_value('nx')) then
        write(stderr,*) 'Option nx must specify a filename!'
        call usage
        stop
     end if
     nx = get_value('nx')
     if (option_exists('ny')) then
        if (.NOT. has_value('ny')) then
           write(stderr,*) 'Option ny must specify a filename!'
           call usage
           stop
        end if
        ny = get_value('ny')
     else
        ny = nx
     end if
  else if (option_exists('ny')) then
     if (.NOT. has_value('ny')) then
        write(stderr,*) 'Option ny must specify a filename!'
        call usage
        stop
     end if
     ny = get_value('ny')
     nx = ny
  else
     nx = 2048
     ny = 2048
  end if

  if (option_exists('background')) then
     ! We have specified a backround file we will subtract
     ! before calculating statistics
     if (.NOT. has_value('background')) then
        write(stderr,*) 'Option background must specify a filename!'
        call usage
        stop
     end if
     backgroundfile = ""
     backgroundfile = get_value('background')
     have_background = .TRUE.
     ! call new(trim(backgroundfile),mar)
     if (verbose) print *,'Read background file'

     ! nx = size(mar)
     ! ny = size(mar)

     ! Allocate the space required for the data
     allocate(background(nx,ny))

     if (verbose) print *,'Retrieving background'
     ! background = get_data(mar)

     if (verbose) print *,'Destroying data'
     ! call destroy(mar)
  end if

  if (option_exists('scale')) then
     ! The user has specified a file containing scales which will 
     ! be applied to the background before subtracting it
     if (have_background) then
        if (.NOT. has_value('scale')) then
           write(stderr,*) 'Option scale must specify a filename!'
           call usage
           stop
        end if
        scalefile = get_value('scale')
        have_scale = .TRUE.
     else
        write(stderr,*) 'Option scale is only valid if a background file is specified'
     end if
  else
     scale = 1.0
     have_scale = .FALSE.
  end if

  if (num_args() < 1) then
     write(stderr,*) 'ERROR! Must supply perkin elmer data file(s) as command-line arguments'
     call usage
     STOP
  end if

  verbose = option_exists('verbose')
  quiet   = option_exists('quiet')

  if (quiet) verbose = .FALSE.

  ! Check what metrics we want to report
  domin    = option_exists('min')
  domax    = option_exists('max')
  domean   = option_exists('mean')
  domedian = option_exists('median')

  ! Our default is to report just median if no command line 
  ! options were used to specify what is required
  if (.NOT. any((/domin, domax, domedian, domean/))) then
     domedian = .TRUE.
  end if

  if (verbose) print *,'("domin = ",L," domax = ",L," domean = ",L," domedian = ",L)', domin, domax, domean, domedian

  ! See if we have specified a value that we will mask out and not
  ! count when calculating the requested statistical quantities
  if (option_exists('maskvalue')) then
     ! Make sure we have a value
     if (.NOT. has_value('maskvalue')) then
        write(stderr,*) 'Option maskvalue must have a value!'
        call usage
        stop
     end if
     ! Get the value
     maskvalue = get_value('maskvalue')
  else
     ! Default is to use the no data value
     maskvalue = no_data
  end if

  ! Allocate the space required for the data
  allocate(data(nx,ny),datamask(nx,ny))

  do while (have_args())

     ! The data files are on the command line
     fname = next_arg()

     if (.not. quiet) write(stdout,'(A,A)',advance='no') trim(fname), " : "

     call get_pe_data(trim(fname), data, data_is_real)
     if (verbose) print *,'Read data file'

     ! call mask_data(mar,maskvalue=maskvalue)

     ! We make the assumption that pixels with zero or negative 
     ! counts are ! outside the detector range and can be regarded 
     ! as having no data (this is mostly the area round the outside 
     ! of the circular image plate)
     where (data <= 0) data = maskvalue

     ! Make a logical mask which shows where we have legitimate
     ! data values
     datamask = (data /= maskvalue)

     ! Check if we have a background to remove before we 
     ! calculate the stats
     if (have_background) then
        if (have_scale) then
           scale = real(get_scale(trim(scalefile)))
           if (verbose) print *,"Scaling background to ",scale
        end if
        where (datamask) data = nint(real(data)-scale*real(background))
     end if

     if (verbose) print *,'Counting pixels'
     npixels = count(datamask)

     if (verbose) print *,'Doing stats'

     if (domin) then
        fmt = '(I0,F0.2)'
        if (.not. quiet) fmt = '("min = ",F0.2,X)'
        write(stdout,fmt,advance='no') minval(data,mask=datamask)
     end if
     if (domax) then
        fmt = '(I0,F0.2)'
        if (.not. quiet) fmt = '("max = ",F0.2,X)'
        write(stdout,fmt,advance='no') maxval(data,mask=datamask)
     end if
     if (domean) then
        fmt = '(F0.2,X)'
        if (.not. quiet) fmt = '("mean = ",F0.2,X)'
        write(stdout,fmt,advance='no') mean(data,datamask) ! sum(data)/npixels
     end if

     if (domedian) then

        allocate(data_1d(npixels))

        ! Make a one dimensional version of the data
        data_1d = pack(data,datamask)

        if (verbose) print *,'Sorting ...'

        ! Sort the pixel values
        call sort(data_1d)

        if (domedian) then
           fmt = '(F0.2,X)'
           if (.not. quiet) fmt = '("median = ",F0.2,X)'
           write(stdout,fmt,advance='no') data_1d(npixels/2)
        end if

        deallocate(data_1d)

     end if

     write(stdout,*)

  end do

  deallocate(data,datamask)

  if (have_background) deallocate(background)

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Usage: pestats --help --min --max --mean --median pefiles(s)'
    write(stderr,*)
    write(stderr,*) '  --help        - print this message.'
    write(stderr,*) '  --min         - calculate minimum'
    write(stderr,*) '  --max         - calculate maximum'
    write(stderr,*) '  --mean        - calculate mean'
    write(stderr,*) '  --median      - calculate median'
    write(stderr,*) '  --quiet       - output only the numbers'
    write(stderr,*) '  --int         - raw pe file stored as 32-bit integers'
!    write(stderr,*) '  --background  - subtract this background file before calculating stats'
    write(stderr,*) '  --scale       - file of scales (one per line) applied to the background'
    write(stderr,*) '                  before subtracting it'
    write(stderr,*)
    write(stderr,*) ' By default, median is output. Specifying an other statistic'
    write(stderr,*) ' will output *only* those statistics'
    write(stderr,*)

  end subroutine usage

  real function mean (array, mask)

    real(rd_kind), intent(in) :: array(:,:)
    logical, intent(in) :: mask(:,:)

    real(rd_kind) :: total
    integer(kind=8) :: npixels

    total = 0
    npixels = 0 
    do i = 1, size(array,2)
       total = total + sum(array(:,i),mask=mask(:,i))
       npixels = npixels + count(mask(:,i))
       ! print *,i,total
    end do
    mean = total / real(npixels,8)

  end function mean

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
          read(unit,*,end=1) scale
          exit
          ! We only get here if we are at the end of the file
1         rewind(unit)
       end do
    end if

  end function get_scale

  subroutine get_pe_data(filename, data, realdata)

    use binary_io, only: read, open, close, binary_filehandle

    character (len=*), intent(in) :: filename
    real(rd_kind), intent(inout)  :: data(:,:)
    logical :: realdata

    ! Local variables
    integer :: i, j, unit, irec, tmp, tmp2
    real(kind=4) :: realtmp

    ! Note that we rely on recl being in bytes, if it is in words this will
    ! need to be changed to myrecl = 1
    integer, parameter :: myrecl = 4
    
    type (binary_filehandle) :: fh

    if (realdata) then

       unit = freeunit()
       
       open(unit=unit, file=trim(fname), form='UNFORMATTED', &
            action='READ', status='OLD', access='DIRECT', recl=myrecl) 
       
       data = 0
       
       irec=3
       
       do i = 1, size(data,2)
          do j = 1, size(data,1)
             read(1, rec=irec) realtmp
             data(j,i) = realtmp
             irec=irec+1
          end do
       end do
       
       close(unit)

    else

       ! Data is stored as 32-bit unsigned integers

        call open(fh, trim(fname))
        
        do i = 1, ny
           do j = 1, nx
              ! Little endian, so we read in 2 lots of two bytes and 
              ! shift the second lot by 2 bytes and add together
              call read(fh, tmp,  2, littleendian=.true.)
              call read(fh, tmp2,  2, littleendian=.true.)
              data(j,i) = real(ishft(tmp2,16) + tmp)
           end do
        end do

        call close(fh)

    end if

  end subroutine get_pe_data

end program pestats
