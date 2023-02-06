module mod_constant
    ! 公用参数
    implicit none

    real, parameter ::  L = 0.0  ! Monin-Obukhov Length [m]
                                   ! L<0.0 unstable
                                   ! L>0.0 stable
                                   ! L=0.0 neutral
    real, parameter ::  k = 0.4  ! Von-Karman constant

    ! typical su, sv, sw in surface layer, obtained from kader and Yaglom (1990)
    ! The empirical parameters
    real, parameter :: SU = 2.70 ! u 
    real, parameter :: SV = 2.50 ! v
    real, parameter :: SW = 1.25 ! w

    ! the surface layer HEIGHT
    real, parameter :: SLH = 100. ! top of the surface layer for unstable [m]
    
    real, parameter :: PI = 4.0*atan(1.0)

    integer, parameter :: FLEN = 256 ! 文件名宽度: length of filename

end module mod_constant
