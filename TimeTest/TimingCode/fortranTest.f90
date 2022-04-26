program testcase
    use iso_c_binding

    real, dimension(100) :: test1
    real, dimension(100,100) :: test2
    real, dimension(100,100,100) :: test3

    interface 
        subroutine initclassptr(ptrtimer) bind(c, name = "initclassptr_")
            use iso_c_binding
            type(c_ptr) :: ptrtimer
        end subroutine initclassptr

        subroutine savetimesf(ptrtimer) bind(c, name = "savetimesf_")
            use iso_c_binding
            type(c_ptr) :: ptrtimer
        end subroutine savetimesf

        subroutine starttimerf(ptrtimer) bind(c, name = "starttimer_")
            use iso_c_binding
            type(c_ptr) :: ptrtimer
        end subroutine starttimerf

        subroutine stoptimerf(ptrtimer) bind(c, name = "stoptimer_")
            use iso_c_binding
            type(c_ptr) :: ptrtimer
        end subroutine stoptimerf

    end interface

    type(c_ptr) :: ptrtimer
    call initclassprt(ptrtimer)
    call StartTimer(ptrtimer)
    do i = 1, 100
        test1(i) = test1(i) + 1
    end do
    call stopTimerf(ptrtimer)
    call SaveTimesF(ptrtimer)
    
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

    

