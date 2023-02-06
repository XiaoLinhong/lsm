module mod_structure
    ! 数据结构
    implicit none

    private

    public param ! parameters that is defined by users

    type param ! parameters that is defined by users

        logical :: debug = .true. ! information in running time

        ! parameters of underlying surface: 下垫面参数
        real :: z0    = 0.01 ! roughness length [m]
        real :: ustar = 0.5  ! friction velocity [m/s]
        ! parameters of source
        real :: zs = 0.1 ! source height

        ! box of simulation
        real    :: zMAX = 1000. ! top boundary [m]
        real    :: zMin = 0.1   ! ground height [m]
        real    :: dz   = 0.1   ! interval of z-axis [m]
        integer :: nz   = 100   ! number of z-axis grid

        real    :: yMAX = 40. ! maximum distance at y-axis [m]
        real    :: yMIN = 0.  ! start distance at y-axis [m]
        real    :: dy   = 0.2 ! interval of y-axis [m]
        integer :: ny   = 100 ! number of y-axis grid

        real    :: xMAX = 20. ! maximum distance at x-axis [m]
        real    :: xMIN = 0.  ! start distance at x-axis [m]
        real    :: dx   = 0.2 ! interval of x-axis [m]
        integer :: nx   = 100 ! number of x-axis grid

        real    :: volume = 0.01 ! volume of a grid [m_3]
       

        ! parameters of simulation
        integer :: nP = 1000000 ! number of particles

    end type param

end module mod_structure
