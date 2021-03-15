module cluster_functions

  implicit none

  ! Implements a weighted union-find search with path compression.

  ! Returns an array "pointer" (not really a pointer) which contains
  ! all the information about the clustering.

  public

  ! The 'empty' value in our pointer arrays
  integer, parameter :: empty = 0

contains

  integer recursive function find_root (ptr, site) result(root)

    integer, intent(inout) :: ptr(:)
    integer, intent(in)    :: site

    if (ptr(site) < 0) then 
       ! Found a root 
       root = site
    else
       ! Go looking for a root, and when you find it, set
       ! this site to it's value (path compression)
       ptr(site) = find_root(ptr,ptr(site))
       root = ptr(site)
       ! root = find_root(ptr,ptr(site))
    end if

    return

  end function find_root

  function cluster_image (sites) result(return_ptr)

    ! Input variables
    logical, dimension(:,:), intent(in)  :: sites

    ! Return type
    integer, dimension(size(sites,1),size(sites,2)) :: return_ptr

    ! This is the local 1D ptr array
    integer, dimension(size(sites)) :: ptr

    ! This is 1D version of the sites array
    logical, dimension(size(sites)) :: sites1d

    ! This is an array which will contain the locations of all the true
    ! values in sites
    integer, dimension(count(sites)) :: trues
    
    ! This is a temporary version of ptr which we will use when labelling
    ! the clusters at the end
    integer, dimension(size(sites)) :: ptrtmp

    integer :: big, nclusters, nx, ny, nxny, ntrue, i, j, r1, s1, r2, s2, nn(8)
    integer :: icounter, cluster_count

    nx = size(sites,1)
    ny = size(sites,2)
    nxny = nx * ny
    ntrue = size(trues)

    big = 0
    if (ntrue > 0) big = 1

    ! print *,'Packing sites1d'
    ! Make a one-dimensional version of sites
    sites1d = pack(sites, mask=.TRUE.)

    ! print *,'Packing trues'
    ! We make a sequential vector from 1 .. nxny and then pick out from
    ! it the points where our sites are .TRUE., thereby giving the index
    ! numbers of the true sites
    ! trues = pack( (/ (i, i=1,nxny) /), mask=sites1d) 
    icounter = 0
    do i = 1, nxny
       if (.NOT. sites1d(i)) cycle
       icounter = icounter + 1
       trues(icounter) = i 
    end do

    ! Initialise ptr array
    ptr = empty

    nclusters = ntrue

    ! print *,'Clustering true values'
    do i = 1, ntrue

       ! Initialise the first root, site and the count for this root
       r1 = trues(i)
       s1 = trues(i)
       ptr(s1) = -1

       nn = nearest_neighbours(s1,nx,nxny)

       ! Loop over the nearest-neighbour contacts
       do j = 1, 8

          s2 = nn(j)
          
          ! Not interested in sites outside our image
          if (s2 < 1) cycle

          ! Not interested in sites which aren't TRUE
          if (.NOT. sites1d(s2)) cycle

          ! This site has been assigned to a cluster already
          if (ptr(s2) /= empty) then

             ! Find to what cluster the second site belongs
             r2 = find_root(ptr,s2)

             ! If not the same, perform a weighted union
             if (r2 /= r1) then

                ! We have one less cluster
                nclusters = nclusters - 1

                if (ptr(r1) > ptr(r2)) then
                   ! The cluster rooted at r1 is smaller than that at r2
                   ! so increment the cluster counter at r2 
                   ptr(r2) = ptr(r2) + ptr(r1)
                   ! Point the site to the new root
                   ptr(r1) = r2
                   ! Make r1 point to r2 so we can have a one-liner at the
                   ! end of this if/then/else statement
                   r1 = r2
                else
                   ! The cluster rooted at r1 is larger than that at r2
                   ! so increment the cluster counter at r1 
                   ptr(r1) = ptr(r1) + ptr(r2)
                   ! Point the site to the new root
                   ptr(r2) = r1
                end if

                ! If necessary update the biggest cluster variable
                if (ptr(r1) < big) big = ptr(r1)

             end if
          end if
       end do
    end do

    big = -1*big

    ! print *,"Found ",nclusters," clusters, the biggest was ",big

    ptrtmp = ptr
    cluster_count = 0

    ! Swing through the pointer array assigning "cluster ids" to
    ! the sites, zero for no cluster id
    do i = 1, ntrue

       s1 = trues(i)

       ! Find the root for this site
       ! root = find_root(ptr,trues(i))
       r1 = find_root(ptr,s1)

       ! print *,i,root,ix,iy,mod(i,nx),i/nx + 1

       if (ptrtmp(r1) < 0) then
          cluster_count = cluster_count + 1
          ptrtmp(r1) = cluster_count
       end if

       ptrtmp(s1) = ptrtmp(r1)
       
    end do

    ! ptr = ptrtmp

    ! print *,'Labelled clusters 1 to ',cluster_count

    return_ptr = reshape(ptrtmp, shape(return_ptr))

  end function cluster_image

  function nearest_neighbours(site, nx, nxny) result(nn)

    integer, intent(in) :: site, nx, nxny

    integer, dimension(8) :: nn

    nn = 0

    nn(1) = site+nx
    nn(2) = site-nx

    ! Check if near right or left edge
    select case (mod(site,nx)) 
       case (0) ! right_edge
          nn(3:5) = 0
          nn(6) = site-1
          nn(7) = nn(1)-1
          nn(8) = nn(2)-1
       case (1) ! left_edge
          nn(3) = site+1
          nn(4) = nn(1)+1
          nn(5) = nn(2)+1
          nn(6:8) = 0
       case default
          nn(3) = site+1
          nn(4) = nn(1)+1
          nn(5) = nn(2)+1
          nn(6) = site-1
          nn(7) = nn(1)-1
          nn(8) = nn(2)-1
    end select

    where (nn > nxny) nn = 0
    where (nn < 0)    nn = 0

  end function nearest_neighbours
 
end module cluster_functions
