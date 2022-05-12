program testcase
    use iso_c_binding

    implicit none

    real, dimension(100) :: test1
    real, dimension(100,100) :: test2
    real, dimension(100,100,100) :: test3

    interface 
        subroutine initptr(ptrtimer) bind(c, name = "initptr")
            use iso_c_binding
            type(c_ptr) :: ptrtimer
        end subroutine initptr

        subroutine savetimesf(ptrtimer) bind(c, name = "savetimesf")
            use iso_c_binding
            type(c_ptr) :: ptrtimer
        end subroutine savetimesf

        subroutine starttimerf(ptrtimer) bind(c, name = "starttimerf")
            use iso_c_binding
            type(c_ptr) :: ptrtimer
        end subroutine starttimerf

        subroutine stoptimerf(ptrtimer) bind(c, name = "stoptimerf")
            use iso_c_binding
            type(c_ptr) :: ptrtimer
        end subroutine stoptimerf

    end interface

    integer  i,j,k

    type(c_ptr) :: ptrtimer
    call initptr(ptrtimer)
    call starttimerf(ptrtimer)
    do i = 1, 100
        test1(i) = test1(i) + 1
    end do
    call stoptimerf(ptrtimer)
    call savetimesf(ptrtimer)
    
    do i = 1, 100
        do j = 1,100
            test2(i,j) = test2(i,j) + 1
        end do
    end do

    do i = 1, 100
        do j = 1,100
            do k = 1,100
                test3(i,j,k) = test3(i,j,k) + 1
            end do
        end do
    end do

end program testcase

    

