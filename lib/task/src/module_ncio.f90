module mod_ncio
    ! 一些简单的NC IO程序
    use netcdf

    implicit none

    contains

    subroutine check(status, fileName, varName)
        ! 检查NC文件是否正确打开
        implicit none 

        integer, intent ( in) :: status
        character(*), intent(in), optional :: fileName
        character(*), intent(in), optional :: varName
    
        if(status /= nf90_noerr) then 
            write(*,'(A)') "ERROR : "//trim(adjustl(nf90_strerror(status)))
            if (present(fileName)) write(*, *) trim(fileName)
            if (present(varName)) write(*, *) trim(varName)
          stop 2
        end if
    end subroutine check 

    subroutine write_3d_data(filename, data3d, xAxis, yAxis, zAxis)
        ! 写出三维数据（中间缓冲数据）
        implicit none
        character(len=*), intent(in) :: filename
        real, dimension(:, :, :), intent(in) :: data3d ! nx, ny, nz
        real, dimension(:), intent(in) :: zAxis  !
        real, dimension(:), intent(in) :: yAxis  !
        real, dimension(:), intent(in) :: xAxis  !

        ! local variable
        integer :: ncid
        integer :: xID, yID, zID
        integer, dimension(4) :: varid ! nvar
    
        call check( nf90_create ( filename, nf90_clobber, ncid ) )
    
        ! define the dimensions
        call check( nf90_def_dim( ncid, "x" , size(data3d, 1), xID ) )
        call check( nf90_def_dim( ncid, "y" , size(data3d, 2), yID ) )
        call check( nf90_def_dim( ncid, "z" , size(data3d, 3), zID ) )
    
        ! define emission variables and variables attributes
        call check( nf90_def_var( ncid, 'x', NF90_FLOAT , (/ xID /), varid(1)  ) )
        call check( nf90_put_att( ncid, varid(1), "description" , 'x-Axis' ))
        call check( nf90_put_att( ncid, varid(1), "units"       ,  'm'))

        call check( nf90_def_var( ncid, 'y', NF90_FLOAT , (/ yID /), varid(2)  ) )
        call check( nf90_put_att( ncid, varid(2), "description" , 'y-Axis' ))
        call check( nf90_put_att( ncid, varid(2), "units"       ,  'm'))

        call check( nf90_def_var( ncid, 'z', NF90_FLOAT , (/ zID /), varid(3)  ) )
        call check( nf90_put_att( ncid, varid(3), "description" , 'z-Axis' ))
        call check( nf90_put_att( ncid, varid(3), "units"       ,  'm'))

        call check( nf90_def_var( ncid, 'conc', NF90_FLOAT , (/ xID, yID, zID /), varid(4)  ) )
        call check( nf90_put_att( ncid, varid(4), "description" , 'Concentration' ))
        call check( nf90_put_att( ncid, varid(4), "units"       ,  'ratio'))
    
        call check( nf90_enddef( ncid ) )

        ! 写出数据
        call check(nf90_put_var(ncid, varid(1), xAxis))
        call check(nf90_put_var(ncid, varid(2), yAxis))
        call check(nf90_put_var(ncid, varid(3), zAxis))
        call check(nf90_put_var(ncid, varid(4), data3d))

        call check(nf90_close(ncid))

    end subroutine write_3d_data

end module mod_ncio
