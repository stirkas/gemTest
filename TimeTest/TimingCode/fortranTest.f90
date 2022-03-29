program testcase
    use iso_c_binding

    real, dimension(100) :: test1
    real, dimension(100,100) :: test2
    real, dimension(100,100,100) :: test3

    interface 
        subroutine initclassptr() bind(c, name = "initclassptr_")
            use iso_c_binding
        end subroutine initclassptr

        subroutine savetimesf() bind(c, name = "savetimesf_")
            use iso_c_binding
        end subroutine savetimesf

        subroutine starttimerf() bind(c, name = "starttimer_")
            use iso_c_binding
        end subroutine starttimerf

        subroutine endtimerf() bind(c, name = "endtimer_")
            use iso_c_binding
        end subroutine endtimerf

    end interface

    call StartTimer()
    do i = 1, 100
        test1(i) = test1(i) + 1
    end do
    call EndTimer()
    call SaveTimesF()
    
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

    

