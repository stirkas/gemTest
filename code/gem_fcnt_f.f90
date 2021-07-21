!-----------------------------------------------------------------------
!-----------------------------------------------------------------------


!   Functions, etc.
!   (things that don't change)

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!v
real function revers(num,n)
integer :: num,n,inum,iquot,irem
real :: rev,power
!    function to reverse the digits of num in base n.

rev = 0.d0
inum = num
power = 1.d0

!turn into do while
11 continue
iquot = int(inum/n)
irem = inum - n*iquot
power = power/n
rev = rev + irem*power
inum = iquot
if(inum.gt.0) goto 11

revers = rev
return
end function revers

!-----------------------------------------------------------------------


!this is the fortran portion that will call the c++ scrbes function
!program Csrcbes
  
!use iso_c_binding
!  
!    !The interface is what links in the c++ function. Biz, gam0, and gam1
!    !are all passed by pointer.
!  
!interface
!    subroutine srcbes(cbiz,cgam0,cgam1) bind(c)
!        use iso_c_binding
!        real (c_double) :: cbiz,cgam0,cgam1
!    end subroutine srcbes
!end interface
!  
    !Here we initialize the variables to be passed into the call
   ! REAL (8) ::biz,gam0,gam1 

    !call srcbes(biz,gam0,gam1)
  
    !end program Csrcbes
