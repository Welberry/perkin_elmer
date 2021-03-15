program find_bad_pixels

  use cmdline_arguments
  use pnm_class
  use statistics

  implicit none

  character(len=2000) :: fname

  type (pnm_object) :: pnm

  integer :: i, j, nx, ny, radius

  integer, allocatable :: data(:,:), output(:,:)

  fname = next_arg()

  call read(pnm,trim(fname))

  nx = size(pnm,1)
  ny = size(pnm,2)

  allocate(data(nx,ny), output(nx,ny))

  data = pnm

  print *,maxval(data),minval(data)

  radius = 4

  output = 0

  do j = 1+radius, ny-radius
     do i = 1+radius, nx-radius
        if (data(i,j) == 0) cycle
        if (abs(data(i,j) - median(pack(data(i-radius:i+radius, j-radius:j+radius),data(i-radius:i+radius,j-radius:j+radius)/=0))) > 20) then
           ! print *,i,j,data(i,j),median(pack(data(i-radius:i+radius, j-radius:j+radius),.TRUE.)),abs(data(i,j) - median(pack(data(i-radius:i+radius, j-radius:j+radius),.TRUE.)))
           output(i,j) = 1
        end if
     end do
  end do

  pnm = output
  
  call write(pnm,'bad_mask.pgm')

end program
