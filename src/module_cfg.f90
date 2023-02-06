module mod_cfg
    ! 读取配置文件
    use mod_constant, only : FLEN
    use mod_structure, only : param

    use flogger, only : log_notice, log_error

    implicit none

    private

    public get_cfg

    contains

    function get_input_name() result(input)
        ! 获取配置文件名称
        implicit none
        character(len=FLEN)   :: input     ! 配置文件名称
        integer               :: fnlen     ! The len of the name of input file. 
        integer               :: ierror    ! Flag for open file error 
        !  read the input file name
        if (command_argument_count() == 0) then
            input = 'namelist.input'
        else
            call get_command_argument(1, input, fnlen, ierror)
        end if
    end function get_input_name

    type(param) function get_cfg result(p)
        ! 获取配置参数
        implicit none

        ! open file
        integer  :: iHandle ! 句柄 
        integer  :: ierror ! Flag for open file error 
        character(len=FLEN) :: fileName ! 配置文件名称

        ! the parameters of underlying surface: 下垫面参数
        real :: z0    = 0.01 ! roughness length [m]
        real :: ustar = 0.5  ! friction velocity [m/s]
        ! the parameters of source
        real :: zs = 0.1 ! source height
        ! box of simulation
        real :: zMAX = 1000.  ! top boundary [m]
        real :: dz   = 0.1    ! interval of z-axis [m]
    
        real :: xMAX = 20. ! maximum distance at x-axis [m]
        real :: dx   = 0.2 ! interval of x-axis [m]

        real :: yMAX = 40. ! maximum distance at y-axis [m]
        real :: dy   = 0.2 ! interval of y-axis [m]

        NAMELIST /share/ z0, ustar, zs, zMAX, dz, xMAX, dx, yMAX, dy

        iHandle = 55
        fileName = get_input_name()
        open(iHandle, file=fileName, status='old', iostat=ierror, form='formatted')
        if(ierror /= 0) call log_error(trim(fileName) // ' not exist ...')

        ! share
        read(iHandle, nml=share)
        p%z0 = z0
        p%ustar = ustar
        p%zs = zs
        p%zMAX = zMAX
        p%dz = dz
        p%xMAX = xMAX
        p%dx = dx
        p%yMAX = yMAX
        p%dy = dy

        ! other
        p%nx = int((p%xMAX-p%zMIN)/dx)
        p%ny = int((p%yMAX-p%yMIN)/dy)
        p%nz = int((p%zMAX-p%zMIN)/dz)
        p%volume = dx*dy*dz

        ! 诊断
        if (p%debug) call log_notice(fileName)
        if (p%debug)    write(*, *) 'ustar', p%ustar
        if (p%debug)    write(*, *) 'height of source', p%zs
        if (p%debug) call log_notice('successfully read configuration file')
        close(iHandle)

    end function get_cfg

end module mod_cfg
