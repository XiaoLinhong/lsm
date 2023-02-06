program main
    !  Lagrangian Stochastic Model (LSM)
    ! constant
    use mod_constant, only : L, K, SU, SV, SW, SLH, PI

    ! namelist
    use mod_structure, only : param
    use mod_cfg, only : get_cfg

    ! random
    use mod_tool, only : get_seed
    use mod_tool, only : get_Gaussian_Random

    ! output
    use mod_ncio, only : write_3d_data

    ! 处理字符串
    use string, only : to_str

    ! 提示信息
    use flogger, only : log_warning, log_notice

    implicit none

    type(param) :: cfg ! 配置信息

    integer :: seed  ! random generator variables

    integer :: i ! value for loop 
    integer :: ii, jj, kk ! index of grid

    real :: x, y, z ! localtion of one particle

    real :: u, v, w ! wind of particle

    real :: uMean ! mean of u
    real :: sigma_u, sigma_v, sigma_w, uw ! wind statistics variables

    real :: TL !  Lagrangian time scale
    real :: Kh !  the diffusion coefficient
    real :: zForStability ! 区分surface layer和PBL
    real :: stab, psi, phi, phi_h ! stability

    ! generalized Langevin equation
    real :: dsigmau2_dz, dsigmav2_dz, dsigmaw2_dz    ! wind statistics variables
    real :: A, b, au, av, aw  ! LS model parameters

    real, allocatable, dimension(:, :, :) :: conc ! output concentration: nx, ny, nz

    real, allocatable, dimension(:) :: zAxis  !
    real, allocatable, dimension(:) :: yAxis  !
    real, allocatable, dimension(:) :: xAxis  !

    real :: tMAX = 0.1
    real :: dt

    ! get the parameters
    cfg = get_cfg()
    ! initialization for the conc array
    allocate( conc(cfg%nx, cfg%ny, cfg%nz) )
    conc = 0.0
    ! 输出坐标
    allocate(zAxis(cfg%nz))
    allocate(yAxis(cfg%ny))
    allocate(xAxis(cfg%nx))
    do i = 1, cfg%nz
       zAxis(i) =  cfg%zMin     + (i-1)*cfg%dz + cfg%dz/2
    end do 

    do i = 1, cfg%ny
       yAxis(i) =  -cfg%yMax/2. + (i-1)*cfg%dy + cfg%dy/2
    end do

    do i = 1, cfg%nx
       xAxis(i) =  cfg%xMin     + (i-1)*cfg%dx + cfg%dx/2
    end do 

    ! 用于生成随机数
    seed = get_seed()

    ! compute constant wind statistics
    uw = -cfg%ustar**2

    do i = 1, cfg%nP ! 每个粒子，单独模拟它的运动轨迹 
        if (mod(i, 100000) == 0) call log_notice( to_str(i, 10) )
        ! 粒子的初始位置
        x = 0.0
        y = 0.0
        z = cfg%zs ! 起始高位为排放源释放点
        
        ONE: do
            ! 与大气稳定度的相关参数计算
            dsigmau2_dz = 0.0
            dsigmav2_dz = 0.0
            dsigmaw2_dz = 0.0
            ! 在这之前，先求up,vp... 为什么要这么做？ 维纳过程，后面的状态和前面相关！
            stab = (1.0-3.0*z/9999.)**(1.0/3.0)
            sigma_u = su*cfg%ustar*stab 
            sigma_v = sv*cfg%ustar*stab
            sigma_w = sw*cfg%ustar*stab
            ! 随机扰动速度: 初始扰动数据
            u = get_Gaussian_Random(seed)*sigma_u
            v = get_Gaussian_Random(seed)*sigma_v
            w = get_Gaussian_Random(seed)*sigma_w

            zForStability = z
            if (L == 0.) then ! neutral
                stab = 1.0
                psi = 1.0
                phi_h = 1.0
            elseif (L > 0.0) then ! stable
                stab = 1.0
                psi = -5.0*z/L ! 和文献不一致 -4.7*z/L
                phi_h = 1.0+5.0*z/L ! 和公式上不一致 1.0 - 3.0*z/L
            else ! unstable
                if (z > SLH) then
                    zForStability = SLH
                    stab = (1.0-3.0*SLH/L)**(1.0/3.0)
                else
                    stab = (1.0-3.0*z/L)**(1.0/3.0)
                    dsigmau2_dz = -(2*SU*SU)*cfg%ustar*cfg%ustar/(stab*L) ! 这是什么变量?
                    dsigmav2_dz = -(2*SV*SV)*cfg%ustar*cfg%ustar/(stab*L)
                    dsigmaw2_dz = -(2*SW*SW)*cfg%ustar*cfg%ustar/(stab*L)
                end if
                !  dimensionless stability correction function(Stull, 1988)
                phi = (1.0-16.0*zForStability/L)**0.25
                psi = 2.0*log((1.0+phi)/2.0)+log((1.0+phi**2)/2)-2.0*atan(phi)+PI/2.0
                ! the stability correction(Hsieh et al., 2000)
                phi_h = 0.32*(0.037-zForStability/L)**(-1.0/3.0)                        
            end if

            uMEAN = cfg%ustar/K * (log(zForStability/cfg%z0)-psi) ! Monin Obukhov Similarity Theory

            !  (Hsieh and Katul, 1997)
            sigma_u = su*cfg%ustar*stab 
            sigma_v = sv*cfg%ustar*stab
            sigma_w = sw*cfg%ustar*stab

            ! 处理时间积分步长
            Kh = K*zForStability*cfg%ustar/phi ! diffusion coefficient
            TL = Kh/(sigma_w*sigma_w) !  K theory
            dt = min(0.02*TL, tMAX)
    
            ! 计算维纳过程的参数:  a flat homogeneous surface layer
            A = 2.0*((sigma_u*sigma_u)*(sigma_w*sigma_w)- uw*uw)
            b = sigma_w*sqrt(2.0/TL)
            au = (b*b)*(uw*w - u*sigma_w*sigma_w)/A + (sigma_w*sigma_w*dsigmau2_dz*u*w - uw*dsigmau2_dz*w*w )/A
            av = (-(b*b)*v)/2.0/sigma_v/sigma_v
            aw = (b*b)*(uw*u - w*sigma_u*sigma_u)/A + 0.5*dsigmaw2_dz + (-uw*dsigmaw2_dz*u*w + sigma_u*sigma_u*dsigmaw2_dz*w*w)/A
            ! 计算当前时刻轨迹点的速度: generalized Langevin equation
            u = u + au*dt + b*sqrt(dt)*get_Gaussian_Random(seed)
            v = v + av*dt + b*sqrt(dt)*get_Gaussian_Random(seed)
            w = w + aw*dt + b*sqrt(dt)*get_Gaussian_Random(seed)

            ! 计算当前时刻轨迹点的位置
            x = x + (uMEAN+u)*dt
            y = y + v*dt
            z = z + w*dt

            ! Concentration calculation, based on Kaplan and Dinar (1996) AE
            ii = int(x/cfg%dx)+1 ! 只往前走
            kk = int(z/cfg%dz)+1 ! 只往上走
            jj = int((y+cfg%yMax/2.0-cfg%dy/2.0)/cfg%dy)+1

            if (ii>0 .and. ii<=cfg%nx .and. jj>0 .and. jj<=cfg%ny .and. kk>0 .and. kk<=cfg%nz) then
                conc(ii, jj, kk) = conc(ii, jj, kk) + 1.0*dt/cfg%nP/cfg%volume
            end if
            
            ! boundary conditions
            ! The signs of u and w are reversed to preserve uw (Wilson and Sawford, (1996) BLM)
            if (z > cfg%zMAX) then ! 顶边界
                w = -w
                u = -u
                z = cfg%zMAX*2.0 - z
            end if

            if (z < cfg%zMIN) then ! 低边界
                w = -w
                u = -u
                z = cfg%zMIN*2.0 - z
            end if

            if (x > cfg%xMAX) exit ONE ! 到达x-axis 时，退出循环

            tMAX = 0.1/abs(w) ! 根据速度, 改变时间积分步长

        end do ONE
    end do

    ! 写出数据
    call write_3d_data('conc.nc',  conc, xAxis, yAxis, zAxis)

end program main
