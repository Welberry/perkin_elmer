program darkfilter

  ! This program is designed to clean up dark current files from the
  ! perkin elmer detector

  use cmdline_arguments
  use string_functions, only: join, real
  use file_functions, only: stderr, freeunit
  use pnm_class
  use precision, only: i8_kind, rs_kind, rd_kind, i4_kind
  use binary_io, only: read, open, close, binary_filehandle
  use iso_varying_string
  use mask_functions
  use fundamental_constants, only: pi
  use cluster_functions
  use variable_array, only: splice, push

  implicit none

  integer, parameter :: dynamic_range = 999999 !16353
  integer, parameter :: subwidth = 128 ! This is the width of the individual 
  ! elements which make up the detector

  ! Note that we rely on recl being in bytes, if it is in words this will
  ! need to be changed to myrecl = 1
  integer, parameter :: myrecl = 4

  integer :: i, j, k, error, nx, ny, unit, unit2, tmp, tmp2, nframes, ios, ubnd, lbnd
  integer :: frame, filenum, radius, numpixels, neg, irec, minimum, start, finish, nxny
  integer :: npixels, totalmasked

  real(kind=4) :: realtmp, header(2), buffer, magnification
  real :: x0, y0, x1, y1, xcentre, ycentre, minblobradius
  character(len=4) :: buf

  character(len=2000) :: fname, outname

  logical :: verbose = .FALSE., autoscale_min, autoscale_max, fixed_min, fixed_max, magnifying
  logical :: fixed_nx, fixed_ny, have_dark, have_badpixel, data_is_real, abs_difference
  logical :: blobbinating

  type (pnm_object) :: pnm

  real :: darkscale, threshold

  real(rd_kind), dimension(:,:), allocatable :: realdata, dark, diff
  integer, dimension(:,:), allocatable       :: clusters
  integer, dimension(:), allocatable         :: indices, onedclusters
  logical, dimension(:,:), allocatable       :: overexposed
  integer, dimension(:), pointer             :: clustertmp

  type (binary_filehandle) :: fh

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)
  type (varying_string) :: myoptions(9)

  myoptions(1) = 'help'
  myoptions(2) = 'quiet'
  myoptions(3) = 'nx'
  myoptions(4) = 'ny'
  myoptions(5) = 'threshold'
  myoptions(6) = 'int'
  myoptions(7) = 'abs'
  myoptions(8) = 'minblobradius'
  myoptions(9) = 'mag'

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

  data_is_real = .not. option_exists('int')
  verbose = .not. option_exists('quiet')
  abs_difference = option_exists('abs')

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

  nxny = nx*ny

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
     threshold = 0.01
  end if

  ! See if we have specified a maginification that we will use to 'inflate'
  ! the masks of our overloaded blobs
  if (option_exists('mag')) then
     ! Make sure we have a value
     if (.NOT. has_value('mag')) then
        write(stderr,*) 'Option mag must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     magnification = get_value('mag')
     magnifying = .true.
  else
     magnifying = .false.
  end if

  ! See if we have specified a minimum overloaded spot size
  if (option_exists('minblobradius')) then
     ! Make sure we have a value
     if (.NOT. has_value('minblobradius')) then
        write(stderr,*) 'Option minblobradius must have a value!'
        call usage
        stop
     end if
     ! Get the minimum radius of our overloaded spot before we
     ! will stick a tail on it and call it a weasel
     minblobradius = get_value('minblobradius')
     blobbinating = .true.
  else
     minblobradius = -1.
     blobbinating = .false.
  end if

  ! Allocate the space required for the data
  allocate(realdata(nx,ny),dark(nx,ny),diff(nx,ny))

  if (magnifying) allocate(clusters(nx,ny),onedclusters(nxny),indices(nxny))

  ! See if we have specified a dark current image
  fname = next_arg()
  have_dark = .TRUE.

  call get_pe_data(trim(fname), dark, data_is_real)

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

  filenum = 0

  do while (have_args())

     filenum = filenum + 1

     ! The files are on the command line
     fname = next_arg()

     ! The data is a tiff file, but we rely on the fact that the first 8 bytes
     ! are the standard tiff header 'II' followed by the number 42. This is 
     ! relying on little endian byte ordering. All bets are off it is big
     ! endian or this is run on a big endian machine.

     call get_pe_data(trim(fname), realdata, data_is_real)

     diff = realdata - dark
     if (abs_difference) then
        diff = abs(diff)
     else
        where (diff < 0.) diff = 0.
     end if

     where (diff /= 0.) diff = diff/dark

     ! Search for a ".tif" suffix .. if we find one than 
     ! replace it with "mask"
     i = index(fname,".dark.tif")
     if (i == 0) i = index(fname,".tif")
     if (i == 0) then
        ! Couldn't find a matching suffix, so look for the end of the
        ! filename string in the character variable
        i = verify(fname," ",back=.TRUE.) + 1
     end if
     outname = fname
     write(outname(i:),'(A)') ".mask"

     if (verbose) print *,'Masking ',count(diff>=threshold),' pixels at threshold ',threshold

     if (magnifying) then
        
        clusters = 0

        ! This function from the cluster_functions module will return an
        ! integer array where all the pixels which are 'true' in the input
        ! are labelled by cluster number. In this case all the true pixels
        ! will be those which exceed our mask value
        clusters = cluster_image( diff >= threshold )

        if (verbose) print *,'Setting up onedclusters ...'

        onedclusters = pack(clusters,mask=.TRUE.)

        if (verbose) print *,'Setting up indices ...'

        ! The much tidier inplace array notation was horrible slow for
        ! values of nxny >= 1000000. Use loop instead.
        ! indices = (/ (i, i=1,nxny) /)
        do i = 1, nxny
           indices(i) = i
        end do

        if (verbose) print *,'Masking data ...'

        print *,maxval(clusters)

        ! The clusters are numbered from 1 .. n, so the maximum value will
        ! be the number of different clusters
        do i = 1, maxval(clusters)

           ! Remove all the elements from the clustertmp array
           npixels = splice(clustertmp, 0)

           ! We make a temporary array which is just the indices of the
           ! cluster we are interested in (the push function returns the
           ! size of the array after the push, which is the number of 
           ! pixels in our cluster)
           npixels = push(clustertmp, pack(indices, mask=(onedclusters==i)))

           ! Find the centre of the cluster
           xcentre = sum(mod(clustertmp,nx))/npixels
           ycentre = sum((clustertmp/ny)+1)/npixels

           radius = sqrt(real(npixels) / pi)

           if (blobbinating .and. (radius > minblobradius)) then

              if (verbose) print '("Cluster ",I4," is centred at x=",I5," y=",I5," and contains ",I5," pixels")'&
                   ,i,xcentre,ny-ycentre,npixels

              ! Determine the 'new' radius of our grown blob -- the old radius
              ! multiplied by the magnification
              radius = magnification*radius

              ! Determine the 'bounding box' of our new blob -- make sure it
              ! stays within the edges of the data. This is not a great algorithm
              ! as we decrease the size of the bounding box, rather than 'clip'
              ! it later on -- which means the radius reduces also. Hmmmmm ...
              x0 = max(xcentre-radius,1.)
              y0 = max(ycentre-radius,1.)

              x1 = min(xcentre+radius,real(nx))
              y1 = min(ycentre+radius,real(ny))

              ! print *,'Masking section: ',x0,y0,x1,y1

              ! Replace a circular area bounded by (x0,y0) and (x1,y1) with the 
              ! current cluster id .i.e. crudely grow the cluster.
              clusters(x0:x1,y0:y1) = circular_mask(clusters(x0:x1,y0:y1),i)
           else
              if (verbose) print '("Cluster ",I4," is centred at x=",I5," y=",I5," and contains ",I5," pixels -- IGNORED")'&
                   ,i,xcentre,ycentre,npixels
           end if

        end do

        call write_mask(outname, clusters > 0)

     else

        call write_mask(outname, (diff>=threshold))

     end if

     if (verbose) print *,'Wrote ',trim(outname)

  end do

  deallocate(realdata, dark, diff)

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Create mask file for Perkin Elmer TIFF file by subtracting a known'
    write(stderr,*) '   good dark file from one contaminated with unerased bragg peaks.'
    write(stderr,*) '   Filter all pixels for which the ratio of their difference to the'
    write(stderr,*) '   good dark file are greater than threshold'
    write(stderr,*)
    write(stderr,*) 'Usage: darkfilter [--help] [--verbose] [--] nifile(s)'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*) '  --nx      - specify x dimension (default 2048)'
    write(stderr,*) '  --ny      - specify y dimension (default 2048)'
    write(stderr,*) '  --quiet   - less verbose output'
    write(stderr,*) '  --threshold - any pixels which are above this value will be masked'
    write(stderr,*) '  --abs     - filter on the absolute difference +ve and -ve'
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

  function circular_mask (square, mask_value) result(circle)
      
    ! A small square section of an image is the input
    integer, dimension(:,:), intent(in) :: square
    integer, intent(in)                 :: mask_value
    
    ! And a circular masked version of this is the output
    integer, dimension(size(square,1),size(square,2)) :: circle
    
    integer :: i, nx, ny, radius
       
    nx = size(square,1)
    ny = size(square,2)
    
    radius = (nx-1)/2
    
    circle = (nint(nx/2.) - spread((/ (i, i=1,nx) /),2,ny))**2 + (nint(ny/2.) - spread((/ (i, i=1,ny) /),1,ny) )**2

    totalmasked = totalmasked + count(circle > radius**2)
       
    where (circle <= radius**2) 
       circle = mask_value
    elsewhere
       circle = square
    end where
    
    ! do i = ny, 1, -1
    !    print '(100I5:)',distance(:,i)
    ! end do
    
  end function circular_mask

end program darkfilter
