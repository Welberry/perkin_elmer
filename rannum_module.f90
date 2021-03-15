module rannum_module
  
  ! implicit none

  private

  !   This random number generator originally appeared in "Toward a Universal
  !   Random Number Generator" by George Marsaglia and Arif Zaman.
  !   Florida State University Report: FSU-SCRI-87-50 (1987)
  !   It was later modified by F. James and published in "A Review of Pseudo-
  !   random Number Generators"
  !   THIS IS THE BEST KNOWN RANDOM NUMBER GENERATOR AVAILABLE.
  !   (However, a newly discovered technique can yield
  !   a period of 10^600. But that is still in the development stage.)
  !   It passes ALL of the tests for random number generators and has a period
  !   of 2^144, is completely portable (gives bit identical results on all
  !   machines with at least 24-bit mantissas in the floating point
  !   representation).
  !   The algorithm is a combination of a Fibonacci sequence (with lags of 97
  !   and 33, and operation "subtraction plus one, modulo one") and an
  !   "arithmetic sequence" (using subtraction).
  !
  !   Use IJ = 1802 & KL = 9373 to test the random number generator. The
  !   subroutine RANMAR should be used to generate 20000 random numbers.
  !   Then display the next six random numbers generated multiplied by 4096*4096
  !   If the random number generator is working properly, the random numbers
  !   should be:
  !           6533892.0  14220222.0  7275067.0
  !           6172232.0  8354498.0   10633180.0

  public :: rannum, rseed
  
contains 
  
  subroutine rseed(ij,kl)
    
    ! this is the initialization routine for the random number 
    ! generator rannum() and must be called once prior to rannum().
    !
    ! note: the seed variables must have values between: 0 <= ij <= 31328
    !                                                    0 <= kl <= 30081
    real u(97), c, cd, cm
    integer i97, j97
      
    common /raset1/ u, c, cd, cm, i97, j97
    save
      
    if( ij .lt. 0  .or.  ij .gt. 31328  .or. kl .lt. 0  .or.  kl .gt. 30081 ) then
       print '(a)', ' the first random number seed must have a value between 0 and 31328'
       print '(a)',' the second seed must have a value between 0 and 30081'
       stop
    endif

    i = mod(ij/177, 177) + 2
    j = mod(ij    , 177) + 2
    k = mod(kl/169, 178) + 1
    l = mod(kl,     169) 
    
    do ii = 1, 97
       s = 0.0
       t = 0.5
       do jj = 1, 24
          m = mod(mod(i*j, 179)*k, 179)
          i = j
          j = k
          k = m
          l = mod(53*l+1, 169)
          if (mod(l*m, 64) .ge. 32) then
             s = s + t
          endif
          t = 0.5 * t
       end do
       u(ii) = s
    end do
    
    c = 362436.0 / 16777216.0
    cd = 7654321.0 / 16777216.0
    cm = 16777213.0 /16777216.0

    i97 = 97
    j97 = 33

    return
  end subroutine rseed

  subroutine rannum(rvec, len)
      
    ! this is a random number generator proposed by george marsaglia in 
    ! florida state university report: fsu-scri-87-50
    ! it was slightly modified by f. james to produce an array of 
    ! pseudorandom numbers.  
    
    real rvec(*)
    real u(97), c, cd, cm
    integer i97, j97
    integer ivec
      
    common /raset1/ u, c, cd, cm, i97, j97  
    save
 
    do ivec = 1, len
       uni = u(i97) - u(j97)
       if( uni .lt. 0.0 ) uni = uni + 1.0
       u(i97) = uni
       i97 = i97 - 1
       if(i97 .eq. 0) i97 = 97
       j97 = j97 - 1
       if(j97 .eq. 0) j97 = 97
       c = c - cd
       if( c .lt. 0.0 ) c = c + cm
       uni = uni - c
       if( uni .lt. 0.0 ) uni = uni + 1.0
       rvec(ivec) = uni
    end do
    return
  end subroutine rannum

  !----------------------------------
  !-----------------------gauss routine-------------------------
  !---------------------------------------------------------
  subroutine gaussm(ix,iy,s,am,v)
    real y(1)
    a=0.
    do i=1,12
       call rannum(y,1)
       a=a+y(1)
    end do
    v=(a-6.)*s+am
    return
  end subroutine gaussm

end module rannum_module
