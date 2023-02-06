module mod_tool
    ! 一些公共小程序

    implicit none

    contains

    integer function get_seed(randNum) result(seed)
        ! get a random number
        implicit none

        integer, optional, intent(in) :: randNum

        if (present(randNum)) then
           seed = randNum + 2*int(secnds(0.0))
        else
           seed = 7654321 + 2*int(secnds(0.0))
        end if

    end function get_seed

    !***********************************************************************
    ! This function generates Gaussian Random Deviates from uniform deviates.
    real function get_Gaussian_Random(idum) result(random)
        ! The function is from Press et al. (1992 p. 289).
        ! this is Box-Muller method !
        implicit none

        integer, intent(inout) :: idum

        real :: v1, v2
        real :: rsq 
        real :: factor
        real, save :: nextRandom = -999.
        logical, save :: isNextTime = .false.

        if (isNextTime) then
            random = nextRandom
            isNextTime = .false.
        else ! next time
            call get_random_in_a_circle(idum, v1, v2, rsq)
            factor = sqrt(-2.*log(rsq)/rsq)
            random = v2*factor
            ! for next time
            nextRandom = v1*factor
            isNextTime = .true.
        end if
    end function get_Gaussian_Random

    recursive subroutine get_random_in_a_circle(idum, v1, v2, rsq)

        integer, intent(inout) :: idum
        real, intent(out) ::  v1, v2
        real, intent(out) :: rsq

        real :: ran

        ! v1 = get_uniform_random(idum)
        ! v2 = get_uniform_random(idum)
        v1 = 2*ran(idum) - 1 ! 编译器自带的随机生成器
        v2 = 2*ran(idum) - 1
        rsq = v1**2 + v2**2

        ! see if they are in the unit circle
        if (rsq >= 1. .or. rsq <= 0) call get_random_in_a_circle(idum, v1, v2, rsq)    
    end subroutine get_random_in_a_circle

    real function get_uniform_random(idum) result(random)
        ! uniform random generator between -1 and 1
        implicit none
        
        integer, intent(inout) :: idum
        ! First proposed by Lewis, Goodman, and Miller in 1969
        integer, parameter :: IA=16807
        integer, parameter :: IM=2147483647
        integer, parameter :: IQ=127773
        integer, parameter :: IR=2836

        real, save :: am
        integer, save :: ix=-1, iy=-1, k

        if (idum <= 0 .or. iy < 0) then
            am=nearest(1.0, -1.0)/IM
            iy=ior(ieor(888889999,abs(idum)),1)
            ix=ieor(777755555, abs(idum))
            idum=abs(idum)+1
        end if
        ix=ieor(ix,ishft(ix,13))
        ix=ieor(ix,ishft(ix,-17))
        ix=ieor(ix,ishft(ix,5))

        k = iy/IQ
        iy = IA*(iy-k*IQ)-IR*k
        if (iy < 0) iy=iy+IM
        random = am*ior(iand(IM, ieor(ix,iy)),1)
        random = 2.*random-1.
    end function get_uniform_random

end module mod_tool
