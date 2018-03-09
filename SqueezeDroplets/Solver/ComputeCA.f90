    subroutine ComputeCA(icord)
   
    implicit none
    include 'dpdflow.h'
    integer icord,num_cell,i,j,n,m,num_slice,n_cell,check_iter
    real * 8 wr_temp,CA_Drop_r_Z(nDpEnd-nWallAtom),Centre_X,Centre_Y,gap_slice,cost_value,cost_value_temp
    real * 8 a_drop,b_drop,a_drop_temp,b_drop_temp,alpha_GD,cost_bias,cost_temp,slope_temp,contactangles_GD
    real * 8 cost_value_top,cost_value_temp_top
    real * 8 a_drop_top,b_drop_top,a_drop_temp_top,b_drop_temp_top,alpha_GD_top,cost_bias_top,cost_temp_top,slope_temp_top,contactangles_GD_top
    real * 8 ,allocatable :: density_value(:,:),density_value_top(:,:)
    real * 8 CL_bottom,CL_top,CL_velocity_bottom,CL_velocity_top

    deltaQ = 0.7
    gap_slice = 0.7
    num_slice = 20
    num_cell = int(min(regionH(2)-1,regionH(1)-1)/deltaQ)+2
    ! bottom parameters
    cost_value = 0.0
    cost_value_temp = 0.0
    a_drop = 1.0
    b_drop = 1.0
    alpha_GD = 0.007
    !top parameters
    a_drop_top = 1.0
    b_drop_top = 1.0
    cost_value_top = 0.0
    cost_value_temp_top = 0.0
    alpha_GD_top = 0.007
    
    if(icord .eq. 0) then
        open(2017,file = './data/ComputeCA.plt')
        write(2017,'(''Variables= "Steps","ContactAngles","CL_bottom","CL_Velocity_bottom","ContactAngles_top","CL_top","CL_Velocity_top"'')')
        write(2017,'(''ZONE'')')
        write(2017,'(''F = POINT'')')
        
        open(2015,file='./data/cost_plt.plt')
        write(2015,'(''Variables= "Steps","cost"'')')
        write(2015,'(''ZONE'')')
        write(2015,'(''F = POINT'')')
        
        num_steps = 0
        Drop_bottom = 0.0
        Drop_top = 0.0
        allocate (each_slice_point(num_slice,2))
        allocate (num_particle(num_slice,num_cell))
        allocate(each_slice_point_top(num_slice,2))
        allocate(num_particle_top(num_slice,num_cell))
        each_slice_point(:,:) = 0.0
        num_particle(:,:) = 0
        each_slice_point_top(:,:) = 0.0
        num_particle_top(:,:) = 0
        CL_temp_bottom = 0.0
        CL_temp_top = 0.0
        
    elseif(icord .eq. 1) then
        num_steps = num_steps + 1
        Centre_X = 0.0
        Centre_Y = 0.0
        do n = nWallAtom + 1 , nDpEnd
            CA_Drop_r_Z(n - nWallAtom) = r(n,3)
            Centre_X = Centre_X + r(n,1)/(nDpEnd - nWallAtom)  ! find the centre of droplet 2017.8.30
            Centre_Y = Centre_Y + r(n,2)/(nDpEnd - nWallAtom)
        enddo
        do i = nWallAtom + 2 , nDpEnd    !使用选择排序对液滴粒子从小到大进行排序    2017.8.30
            do j = i , nWallAtom + 2 , -1
                if (CA_Drop_r_Z(j - nWallAtom) .lt. CA_Drop_r_Z(j - nWallAtom - 1)) then
                    wr_temp = CA_Drop_r_Z(j - nWallAtom)
                    CA_Drop_r_Z(j - nWallAtom) = CA_Drop_r_Z(j - nWallAtom - 1)
                    CA_Drop_r_Z(j - nWallAtom - 1) = wr_temp
                endif
            enddo
        enddo
        
        do n = 1 , 30
            Drop_bottom = Drop_bottom + CA_Drop_r_Z(n)/30.0 ! 取最低的50个粒子的z坐标平均值作为液滴的底
        enddo
        do n = nDpEnd - 30 + 1 - nWallAtom ,  nDpEnd - nWallAtom
            Drop_top = Drop_top + CA_Drop_r_Z(n)/30.0  ! 取最高的50个粒子的z坐标平均值作为液滴的顶
        enddo
        ! 处理从底部开始的粒子,自下往上
        do i = 1 , num_slice
            do  n = nWallAtom + 1, nDpEnd
                if(r(n,3) .ge. Drop_bottom/num_steps+(i-1)*gap_slice .and. r(n,3) .lt. Drop_bottom/num_steps+i*gap_slice) then
                     if (sqrt((r(n,1) - Centre_X)**2 + (r(n,2) - Centre_Y)**2) .le. 1.0) then
                         n_cell = 1
                     else
                         n_cell = int((sqrt((r(n,1) - Centre_X)**2 + (r(n,2) - Centre_Y)**2)-1.0)/deltaQ) + 2
                     endif
                     
                     if (n_cell .le. num_cell) then
                         num_particle(i,n_cell) = num_particle(i,n_cell) + 1
                         !write(*,*) "num_particle(",i,",",n_cell,") = ", num_particle(i,n_cell)
                     endif
                endif
            enddo
        enddo
        !处理从顶部开始的粒子,自上往下
        do i = 1 , num_slice
            do  n = nWallAtom + 1, nDpEnd
                if(r(n,3) .le. Drop_top/num_steps-(i-1)*gap_slice .and. r(n,3) .gt. Drop_top/num_steps-i*gap_slice) then
                     if (sqrt((r(n,1) - Centre_X)**2 + (r(n,2) - Centre_Y)**2) .le. 1.0) then
                         n_cell = 1
                     else
                         n_cell = int((sqrt((r(n,1) - Centre_X)**2 + (r(n,2) - Centre_Y)**2)-1.0)/deltaQ) + 2
                     endif
                     
                     if (n_cell .le. num_cell) then
                         num_particle_top(i,n_cell) = num_particle_top(i,n_cell) + 1
                         !write(*,*) "num_particle_top(",i,",",n_cell,") = ", num_particle_top(i,n_cell)
                     endif
                endif
            enddo
        enddo
        
    elseif(icord .eq. 2) then 
        allocate(density_value(num_slice,num_cell))
        allocate(density_value_top(num_slice,num_cell))
        density_value(:,:) = 0.0
        density_value_top(:,:) = 0.0
        contactangles_GD = 0.0
        ! 处理自下而上的粒子(计算每个格子的密度值)
        do i = 1 , num_slice
            if(num_particle(i,1) .gt. 0) then
                do j = 1 , num_cell
                    if(j .eq. 1) then
                        density_value(i,j) = (1.0*num_particle(i,j)/num_steps)/(gap_slice*pi*1.0**2)
                    else
                        density_value(i,j) = (1.0*num_particle(i,j)/num_steps)/(gap_slice*pi*(1.0+deltaQ*(j-1))**2 - gap_slice*pi*(1.0+deltaQ*(j-2))**2)
                    endif  
                    if(density_value(i,j) .le. density/20 .and. density_value(i,j+1) .le. density/20 .and. density_value(i,j+2) .le. density/20 ) then
                        each_slice_point(i,1) = 1.0 + deltaQ/20 + deltaQ*(j-2) !x方向坐标
                        each_slice_point(i,2) = Drop_bottom/num_steps + (i * gap_slice + (i - 1) * gap_slice)/2.0  ! z方向坐标
                        goto 91
                    endif
                enddo
            endif
91          continue
        enddo
        ! 处理自上而下的粒子
        do i = 1 , num_slice
            if(num_particle_top(i,1) .gt. 0) then
                do j = 1 , num_cell
                    if(j .eq. 1) then
                        density_value_top(i,j) = (1.0*num_particle_top(i,j)/num_steps)/(gap_slice*pi*1.0**2)
                    else
                        density_value_top(i,j) = (1.0*num_particle_top(i,j)/num_steps)/(gap_slice*pi*(1.0+deltaQ*(j-1))**2 - gap_slice*pi*(1.0+deltaQ*(j-2))**2)
                    endif  
                    if(density_value_top(i,j) .le. density/20 .and. density_value_top(i,j+1) .le. density/20 .and. density_value_top(i,j+2) .le. density/20 ) then
                        each_slice_point_top(i,1) = 1.0 + deltaQ/20 + deltaQ*(j-2) !x方向坐标
                        each_slice_point_top(i,2) = Drop_top/num_steps - (i * gap_slice + (i - 1) * gap_slice)/2.0  ! z方向坐标
                        goto 147
                    endif
                enddo
            endif
147         continue
        enddo
     
        open(2016,file = './data/Sp.plt')
        write(2016,'(''Variables= "Y","X"'')')
        write(2016,'(''ZONE'')')
        write(2016,'(''F = POINT'')')
        do i = 1, num_slice
            write(2016,'(3f16.6)') each_slice_point(i,2),each_slice_point(i,1)
        enddo
        write(2016,*) "--------Half of the high--------"
        do i = 1 , int((Drop_top/num_steps - Drop_bottom/num_steps)/2/gap_slice)
            write(2016,'(3f16.6)') each_slice_point(i,2),each_slice_point(i,1)
        enddo
        
        write(2016,*) "**********from top to bottom**************"
        do i = 1, num_slice
            write(2016,'(3f16.6)') each_slice_point_top(i,2),each_slice_point_top(i,1)
        enddo
        write(2016,*) "--------Half of the high--------"
        do i = 1 , int((Drop_top/num_steps - Drop_bottom/num_steps)/2/gap_slice)
            write(2016,'(3f16.6)') each_slice_point_top(i,2),each_slice_point_top(i,1)
        enddo
        
        CL_bottom = each_slice_point(1,1)
        CL_top = each_slice_point_top(1,1)
        
! 利用回归算法对液滴形状进行拟合----------x^2+(y-a)^2=b^2------------------------
        ! Contact Angles from bottom to top 
        cost_value = 0.0
        cost_temp = 0.0
        check_iter = 1
        m = 0
        do while (check_iter .eq. 1) 
            m = m + 1
            do i = 1 , int((Drop_top/num_steps - Drop_bottom/num_steps)/2/gap_slice)
                cost_value = cost_value + (((b_drop**2-(each_slice_point(i,2)-a_drop)**2-each_slice_point(i,1)**2))**2)/num_slice/4.0
            enddo
            do i = 1 , int((Drop_top/num_steps - Drop_bottom/num_steps)/2/gap_slice)  ! 利用梯度下降算法求得a,b的全局最优解
                a_drop_temp = a_drop
                b_drop_temp = b_drop
                a_drop = a_drop_temp - alpha_GD*(b_drop_temp**2-(each_slice_point(i,2)-a_drop_temp)**2-each_slice_point(i,1)**2)*(each_slice_point(i,2)-a_drop_temp)/num_slice
                b_drop = b_drop_temp - alpha_GD*((b_drop_temp**2-(each_slice_point(i,2)-a_drop_temp)**2-each_slice_point(i,1)**2))*b_drop_temp/num_slice
            enddo
            write(2015,'(i6,f16.6)') m,cost_value
            if (abs(cost_value - cost_temp) .le. 0.0000001) then
                check_iter = 0
            endif
            cost_bias = cost_value - cost_temp
            cost_temp = cost_value
            cost_value = 0.0
        enddo
        slope_temp = (a_drop - Drop_bottom/num_steps)/sqrt(b_drop**2 - (Drop_bottom/num_steps - a_drop)**2)
        contactangles_GD = 90 + atan(slope_temp)*180/pi
        
        ! Contact Angles from top to bottom
        cost_value_top = 0.0
        cost_temp_top = 0.0
        check_iter = 1
        m = 0
        do while (check_iter .eq. 1) 
            m = m + 1
            do i = 1 , int((Drop_top/num_steps - Drop_bottom/num_steps)/2/gap_slice)
                cost_value_top = cost_value_top + (((b_drop_top**2-(each_slice_point_top(i,2)-a_drop_top)**2-each_slice_point_top(i,1)**2))**2)/num_slice/4.0
            enddo
            do i = 1 , int((Drop_top/num_steps - Drop_bottom/num_steps)/2/gap_slice)  ! 利用梯度下降算法求得a,b的全局最优解
                a_drop_temp_top = a_drop_top
                b_drop_temp_top = b_drop_top
                a_drop_top = a_drop_temp_top - alpha_GD_top*(b_drop_temp_top**2-(each_slice_point_top(i,2)-a_drop_temp_top)**2-each_slice_point_top(i,1)**2)*(each_slice_point_top(i,2)-a_drop_temp_top)/num_slice
                b_drop_top = b_drop_temp_top - alpha_GD_top*((b_drop_temp_top**2-(each_slice_point_top(i,2)-a_drop_temp_top)**2-each_slice_point_top(i,1)**2))*b_drop_temp_top/num_slice
            enddo
            write(2015,'(i6,f16.6)') m,cost_value_top
            if (abs(cost_value_top - cost_temp_top) .le. 0.0000001) then
                check_iter = 0
            endif
            cost_bias_top = cost_value_top - cost_temp_top
            cost_temp_top = cost_value_top
            cost_value_top = 0.0
        enddo
        slope_temp_top = (a_drop_top - Drop_top/num_steps)/sqrt(b_drop_top**2 - (Drop_top/num_steps - a_drop_top)**2)
        if(slope_temp_top .le. 0) then
            contactangles_GD_top = 90 + abs(atan(slope_temp_top)*180/pi)
        else
            contactangles_GD_top = 90 - atan(slope_temp_top)*180/pi
        endif
        
        CL_velocity_bottom = (CL_bottom - CL_temp_bottom)/num_steps/deltaT
        CL_velocity_top = (CL_top - CL_temp_top)/num_steps/deltaT
        CL_temp_bottom = CL_bottom
        CL_temp_top = CL_top
        
        write(*,*) 'a = ',a_drop, '  ','b = ',b_drop,'a_top = ',a_drop_top, '  ','b_top = ',b_drop_top
        write(2016,'(4f16.6)') a_drop,b_drop,a_drop_top,b_drop_top
        close(2016)
        write(2017,'(i6,6f16.6)') stepCount, contactangles_GD, CL_bottom, CL_velocity_bottom, contactangles_GD_top, CL_top, CL_velocity_top
        
        Drop_bottom = 0.0
        Drop_top = 0.0
        num_steps = 0
        num_particle(:,:) = 0
        each_slice_point(:,:) = 0.0
        num_particle_top(:,:) = 0
        each_slice_point_top(:,:) = 0.0
        deallocate(density_value)
        deallocate(density_value_top)
    elseif(icord .eq. 3) then 
        deallocate(each_slice_point)
        deallocate(num_particle)
        deallocate(each_slice_point_top)
        deallocate(num_particle_top)
        close(2017)
        close(2015)
    endif

    return 
    end