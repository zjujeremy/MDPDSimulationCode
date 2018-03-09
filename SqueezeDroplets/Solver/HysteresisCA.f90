    subroutine HysteresisCA(icord)
   
    implicit none
    include 'dpdflow.h'
    
    integer icord, ifirst, i_ifirst, i, j, k, n, g_get_point_flog, g_get_point_flog_YDF, i_GD, i_RDF, j_RDF, g_m
    integer end_slice, start_slice
    integer g_slice_No, g_degree_No, g_cell_No, g_check_iter, g_num_cell_output
    real*8 g_Center_X, g_Center_Y, g_alpha_GD, g_contactangles_GD_line, g_each_CA(int(360.0/g_delta_degree)), g_each_CA_YDF(int(360.0/g_delta_degree)), g_each_CA_YDF_top(int(360.0/g_delta_degree)), g_each_CA_line(int(360.0/g_delta_degree))
    real*8 g_cost_value, g_cost_temp, g_a_drop, g_b_drop, g_a_drop_temp, g_b_drop_temp, g_slope_temp, g_contactangles_GD, g_a_drop_line, g_b_drop_line
    real*8 g_a_drop_YDF, g_b_drop_YDF, g_contactangles_GD_YDF, g_CA_temp
    real*8 g_a_drop_YDF_top, g_b_drop_YDF_top, g_contactangles_GD_YDF_top, g_CA_bottom, g_CA_top, g_CL_A_Velo_sum, g_CL_R_Velo_sum
    character*8 stc, i_stc
 
    if(icord .eq. 0) then 
        g_output_id = 339
        open(g_output_id, file = './data/CA_hysteresis.plt')    
        write(g_output_id, '(''TITLE = "contact angles in each degree(/30) VS steps"'')')
        write(g_output_id, '(''Variables= "Steps","D1","D2","D3","D4","D5","D6","D7","D8","D9","D10","D11","D12"'')') 
	    write(g_output_id, '(''ZONE'')')
	    write(g_output_id, '(''F = POINT'')')
        
        g_output_CL_id = 340
        open(g_output_CL_id, file = './data/CL_hysteresis.plt')    
        write(g_output_CL_id, '(''TITLE = "contact line in A&R VS steps"'')')
        write(g_output_CL_id, '(''Variables= "Steps","CL_Advance","Cl_Ad_velo","CL_Re","CL_Re_velo"'')') 
	    write(g_output_CL_id, '(''ZONE'')')
	    write(g_output_CL_id, '(''F = POINT'')')
        
        open(2018,file = './data/ComputeCA.plt')
        write(2018,'(''Variables= "Steps","ContactAngles_bottom","CL_bottom","CL_Velocity_bottom","CL_Velocity_bottom_sum","ContactAngles_top","CL_top","CL_Velocity_top","CL_Velocity_top_sum"'')')
        write(2018,'(''ZONE'')')
        write(2018,'(''F = POINT'')')
        
        g_num_steps = 0
        g_positionOver = 0
        g_Center_aver_X = 0.0
        g_Center_aver_Y = 0.0
        
        g_num_slice = int(region(3)/g_gap_slice)+1
        g_num_degree = int(360.0/g_delta_degree)
        g_num_cell = int((max(region(2), region(1))/g_first_deltaQ)**2)+1
        allocate(g_num_cell_particle(g_num_slice, g_num_degree, g_num_cell)) ! for RDF
        allocate(g_density_cell(g_num_slice, g_num_degree, g_num_cell))
        allocate(g_each_degree_point(g_num_degree, g_num_slice, 2))
        allocate(g_each_degree_point_YDF(g_num_degree, g_num_slice, 2))
        allocate(g_YDF_sum_coor_aver(g_num_slice, g_num_degree, g_num_cell))
        allocate(g_YDF_sum_num_parti(g_num_slice, g_num_degree, g_num_cell))
        allocate(g_YDF_sum_coor(g_num_slice, g_num_degree, g_num_cell))
        g_num_cell_particle(:, :, :) = 0
        g_density_cell(:, :, :) = 0.0
        g_each_degree_point(:, :, :) = 0.0
        g_each_degree_point_YDF(:, :, :) = 0.0
        g_YDF_sum_coor(:, :, :) = 0.0
        g_YDF_sum_coor_aver(:, :, :) = 0.0
        g_YDF_sum_num_parti(:, :, :) = 0
        
        call find_bottomAndtop
        
    elseif(icord .eq. 1) then 
        g_num_steps = g_num_steps + 1
        g_YDF_sum_num_parti(:, :, :) = 0
        g_YDF_sum_coor(:, :, :) = 0.0
        g_Center_X = 0.0
        g_Center_Y = 0.0
        do n = nWallAtom + 1 , nDpEnd
            g_Center_X = g_Center_X + r(n,1)/(nDpEnd - nWallAtom)  ! find the centre of droplet 2017.11.29
            g_Center_Y = g_Center_Y + r(n,2)/(nDpEnd - nWallAtom)
        enddo
        if(stepCount .le. g_outputFlu) then
            g_Center_aver_X = g_Center_aver_X + g_Center_X/(g_outputFlu-hysteresis_start)
            g_Center_aver_Y = g_Center_aver_Y + g_Center_Y/(g_outputFlu-hysteresis_start)
        else
            g_Center_aver_X = g_Center_aver_X + g_Center_X/g_outputFlu
            g_Center_aver_Y = g_Center_aver_Y + g_Center_Y/g_outputFlu
        endif
    
        do n = nWallAtom + 1, nDpEnd
            g_slice_No = int(abs(r(n, 3)-g_Drop_bottom)/g_gap_slice)+1
            if(r(n, 2) .ge. g_Center_Y) then
                g_CA_temp = acos((r(n, 1)-g_Center_X)/sqrt((r(n, 1)-g_Center_X)**2+(r(n, 2)-g_Center_Y)**2))*180.0/pi
                if(g_CA_temp .le. (g_delta_degree/2.0)) then
                    g_degree_No = 1
                elseif(g_CA_temp .gt. (180.0-g_delta_degree/2.0)) then
                    g_degree_No = int(180.0/g_delta_degree+1)
                else
                    g_degree_No = int((g_CA_temp-g_delta_degree/2.0)/g_delta_degree)+2
                endif
            elseif(r(n, 2) .lt. g_Center_Y) then
                g_CA_temp = acos(-(r(n, 1)-g_Center_X)/sqrt((r(n, 1)-g_Center_X)**2+(r(n, 2)-g_Center_Y)**2))*180.0/pi+180.0
                if(g_CA_temp .le. (180.0+g_delta_degree/2.0)) then
                    g_degree_No = int(180.0/g_delta_degree+1)
                elseif(g_CA_temp .gt. (360.0-g_delta_degree/2.0)) then
                    g_degree_No = 1
                else
                    g_degree_No = int((g_CA_temp-g_delta_degree/2.0)/g_delta_degree)+2
                endif
            endif
            g_cell_No = int((sqrt((r(n, 1)-g_Center_X)**2+(r(n, 2)-g_Center_Y)**2)/g_first_deltaQ)**2)+1
            if(g_cell_No .gt. g_num_cell) then
                g_positionOver = g_positionOver + 1
                write(*, *) "g_positionOver = ", g_positionOver
                goto 82
            endif
            g_num_cell_particle(g_slice_No, g_degree_No, g_cell_No) = g_num_cell_particle(g_slice_No, g_degree_No, g_cell_No) + 1 ! for RDF
            g_YDF_sum_coor(g_slice_No, g_degree_No, g_cell_No) = g_YDF_sum_coor(g_slice_No, g_degree_No, g_cell_No) + r(n, 3)     ! for YDF
            g_YDF_sum_num_parti(g_slice_No, g_degree_No, g_cell_No) = g_YDF_sum_num_parti(g_slice_No, g_degree_No, g_cell_No) + 1
82          continue
        enddo
        do i = 1, g_num_slice
            do j = 1, g_num_degree
                do k = 1, g_num_cell
                    if(g_YDF_sum_num_parti(i, j, k) .eq. 0) then
                        g_YDF_sum_coor_aver(i, j, k) = g_YDF_sum_coor_aver(i, j, k) + 0.0
                    else
                        g_YDF_sum_coor_aver(i, j, k) = g_YDF_sum_coor_aver(i, j, k) + g_YDF_sum_coor(i, j, k)/real(g_YDF_sum_num_parti(i, j, k))
                    endif
                enddo
            enddo
        enddo
        
    elseif(icord .eq. 2) then 
        call find_bottomAndtop
        
        write(stc,'(i8)')stepCount
        if(stepCount .lt. 10) then
            ifirst=index(stc,' ')+7
        elseif(stepCount .lt. 100) then
            ifirst=index(stc,' ')+6
        elseif(stepCount .lt. 1000) then
            ifirst=index(stc,' ')+5
        elseif(stepCount.lt.10000) then
            ifirst=index(stc,' ')+4
        elseif(stepCount.lt.100000) then
            ifirst=index(stc,' ')+3
        elseif(stepCount.lt.1000000) then
            ifirst=index(stc,' ')+2   
        elseif(stepCount.lt.10000000) then
            ifirst=index(stc,' ')+1 
        endif
        
        open(91, file = './data/CA_around_the_circumfrence/CA_degree'//stc(ifirst:)//'.plt')    
        write(91,'(''TITLE = "contact angles around the circumfrence"'')')
        write(91,'(''Variables= "CL/Circum", "Conatact Angles YDF", "Conatact Angles line", "Conatact Angles RDF"'')') 
	    write(91,'(''ZONE'')')
	    write(91,'(''F = POINT'')')
        
        do i = 1, g_num_degree
            do j = 1, g_num_slice
                g_get_point_flog = 1
                g_get_point_flog_YDF = 1
                do k = 1, g_num_cell-3
                    g_density_cell(j, i, k) = (1.0*g_num_cell_particle(j, i, k)/g_num_steps)/((g_gap_slice*pi*g_first_deltaQ**2)/g_num_degree)
                    g_density_cell(j, i, (k+1)) = (1.0*g_num_cell_particle(j, i, (k+1))/g_num_steps)/((g_gap_slice*pi*g_first_deltaQ**2)/g_num_degree)
                    g_density_cell(j, i, (k+2)) = (1.0*g_num_cell_particle(j, i, (k+2))/g_num_steps)/((g_gap_slice*pi*g_first_deltaQ**2)/g_num_degree)
                    g_density_cell(j, i, (k+3)) = (1.0*g_num_cell_particle(j, i, (k+3))/g_num_steps)/((g_gap_slice*pi*g_first_deltaQ**2)/g_num_degree)
                    g_YDF_sum_coor_aver(j, i, k) = g_YDF_sum_coor_aver(j, i, k)/real(g_num_steps)
                    if(k .eq. (g_num_cell-3)) then
                        g_YDF_sum_coor_aver(j, i, (k+1)) = g_YDF_sum_coor_aver(j, i, (k+1))/real(g_num_steps)
                        g_YDF_sum_coor_aver(j, i, (k+2)) = g_YDF_sum_coor_aver(j, i, (k+2))/real(g_num_steps)
                        g_YDF_sum_coor_aver(j, i, (k+3)) = g_YDF_sum_coor_aver(j, i, (k+3))/real(g_num_steps)
                    endif
                    !if(g_density_cell(j, i, k) .le. density/8.0 .and. g_density_cell(j, i, (k+1)) .le. density/8.0 .and. &
                    !         g_density_cell(j, i, (k+2)) .le. density/8.0 .and. g_density_cell(j, i, (k+3)) .le. density/8.0 .and. g_get_point_flog .eq. 1) then
                    !    g_each_degree_point(i, j, 1) = g_first_deltaQ*(sqrt(real(k))+sqrt(real(k-1)))/2.0
                    !    g_each_degree_point(i, j, 2) = g_Drop_bottom + (j-0.5)*g_gap_slice
                    !    g_get_point_flog = 0
                    !endif
                    if((g_get_point_flog_YDF .eq. 1 .and. k .lt. (g_num_cell-3) .and. g_YDF_sum_coor_aver(j, i, k) .gt. (j*g_gap_slice+g_Drop_bottom) .and. g_YDF_sum_coor_aver(j, i, (k+1))/real(g_num_steps) .gt. (j*g_gap_slice+g_Drop_bottom) .and. g_YDF_sum_coor_aver(j, i, (k+2))/real(g_num_steps) .gt. (j*g_gap_slice+g_Drop_bottom)) .or. &
                           (g_get_point_flog_YDF .eq. 1 .and. k .eq. (g_num_cell-3) .and. g_YDF_sum_coor_aver(j, i, k) .gt. (j*g_gap_slice+g_Drop_bottom) .and. g_YDF_sum_coor_aver(j, i, (k+1)) .gt. (j*g_gap_slice+g_Drop_bottom) .and. g_YDF_sum_coor_aver(j, i, (k+2)) .gt. (j*g_gap_slice+g_Drop_bottom))) then
                        g_each_degree_point_YDF(i, j, 1) = g_first_deltaQ*sqrt(real(k-1))
                        if(k .eq. 1) then 
                            g_each_degree_point_YDF(i, j, 2) = g_YDF_sum_coor_aver(j, i, k)
                        else
                            g_each_degree_point_YDF(i, j, 2) = g_YDF_sum_coor_aver(j, i, (k-1))
                        endif
                        g_get_point_flog_YDF = 0
                    endif
                enddo  
            enddo
            !if(i .eq. 1) then
            !    write(*, *) (g_each_degree_point_YDF(1, j_RDF, 2), j_RDF=1, 10)
            !    write(*, *) (g_each_degree_point_YDF(1, j_RDF, 1), j_RDF=1, 10)
            !endif
            
            ! output the RDF for the each degree
            write(i_stc,'(i8)') i
            if(i .lt. 10) then
                i_ifirst=index(i_stc,' ')+7
            elseif(i .lt. 100) then
                i_ifirst=index(i_stc,' ')+6
            elseif(i .lt. 1000) then
                i_ifirst=index(i_stc,' ')+5
            elseif(i .lt. 10000) then
                i_ifirst=index(i_stc,' ')+4
            elseif(i .lt. 100000) then
                i_ifirst=index(i_stc,' ')+3
            elseif(i .lt. 1000000) then
                i_ifirst=index(i_stc,' ')+2   
            elseif(i .lt. 10000000) then
                i_ifirst=index(i_stc,' ')+1 
            endif
            !if(mod((i-1), 3) .eq. 0) then
            if(i .eq. 1) then
                open(89, file = './data/Hysteresis_ContactAngles/RDF_'//stc(ifirst:)//'-'//i_stc(i_ifirst:)//'.plt')    
                write(89,'(''TITLE = "RDF In Each Degree"'')')
                write(89,'(''Variables= "Radius","RDF_1","RDF_2","RDF_3","RDF_4","RDF_5","RDF_6","RDF_7","RDF_8","RDF_9","RDF_10"'')') 
	            write(89,'(''ZONE'')')
	            write(89,'(''F = POINT'')')
                
                open(891, file = './data/YDF/YDF_'//stc(ifirst:)//'-'//i_stc(i_ifirst:)//'.plt')    
                write(891,'(''TITLE = "YDF In Each Degree"'')')
                write(891,'(''Variables= "Radius","YDF_1","YDF_2","YDF_3","YDF_4","YDF_5","YDF_6","YDF_7","YDF_8","YDF_9","YDF_10"'')') 
	            write(891,'(''ZONE'')')
	            write(891,'(''F = POINT'')')
                g_num_cell_output = int((2*RdsDp/g_first_deltaQ)**2)
                do i_RDF = 1, g_num_cell_output
                    write(89,'(11f16.6)') (sqrt(real(i_RDF))+sqrt(real(i_RDF-1)))*g_first_deltaQ/2, (g_density_cell(j_RDF, i, i_RDF), j_RDF=1, 10) 
                    write(891,'(11f16.6)') (sqrt(real(i_RDF))+sqrt(real(i_RDF-1)))*g_first_deltaQ/2, (g_YDF_sum_coor_aver(j_RDF, i, i_RDF), j_RDF=1, 10) 
                enddo
                close(89)
                close(891)
            endif

! 利用梯度下降算法求得a,b的全局最优解----------x^2+(y-a)^2=b^2------------------------
            !g_cost_value = 0.0
            !g_cost_temp = 0.0
            !g_check_iter = 1
            !g_m = 0
            !g_a_drop = 1.0
            !g_b_drop = 1.0
            !g_alpha_GD = 0.001
            !do while (g_check_iter .eq. 1) 
            !    g_m = g_m + 1
            !    do i_GD = 1, g_num_slice_for_GDalgo
            !        if(g_each_degree_point(i, i_GD, 1) .gt. g_first_deltaQ) then
            !            g_cost_value = g_cost_value + (((g_b_drop**2-(g_each_degree_point(i, i_GD, 2)-g_a_drop)**2-g_each_degree_point(i, i_GD, 1)**2))**2)/g_num_slice_for_GDalgo/4.0
            !        endif
            !    enddo
            !    do i_GD = 1, g_num_slice_for_GDalgo
            !        if(g_each_degree_point(i, i_GD, 1) .gt. g_first_deltaQ) then
            !            g_a_drop_temp = g_a_drop
            !            g_b_drop_temp = g_b_drop
            !            g_a_drop = g_a_drop_temp - g_alpha_GD*(g_b_drop_temp**2-(g_each_degree_point(i, i_GD, 2)-g_a_drop_temp)**2-g_each_degree_point(i, i_GD, 1)**2)*(g_each_degree_point(i, i_GD, 2)-g_a_drop_temp)/g_num_slice_for_GDalgo
            !            g_b_drop = g_b_drop_temp - g_alpha_GD*((g_b_drop_temp**2-(g_each_degree_point(i, i_GD, 2)-g_a_drop_temp)**2-g_each_degree_point(i, i_GD, 1)**2))*g_b_drop_temp/g_num_slice_for_GDalgo
            !        endif
            !    enddo
            !    if (abs(g_cost_value - g_cost_temp) .le. 0.0000001) then
            !        g_check_iter = 0
            !    endif
            !    g_cost_temp = g_cost_value
            !    g_cost_value = 0.0
            !enddo
            !g_slope_temp = (g_a_drop - g_Drop_bottom)/sqrt(g_b_drop**2 - (g_Drop_bottom - g_a_drop)**2)
            !g_contactangles_GD = 90 + atan(g_slope_temp)*180/pi 

! 利用梯度下降算法求得a,b的全局最优解----------x^2+(y-a)^2=b^2---------- 利用YDF取得的数据点--------------
            g_cost_value = 0.0
            g_cost_temp = 0.0
            g_check_iter = 1
            g_m = 0
            g_a_drop_YDF = -6.347737
            g_b_drop_YDF = 1.0
            g_alpha_GD = 0.001
            do while (g_check_iter .eq. 1) 
                g_m = g_m + 1
                do i_GD = 1, g_num_slice_for_GDalgo_YDF
                    if(g_each_degree_point_YDF(i, i_GD, 1) .gt. g_first_deltaQ) then
                        g_cost_value = g_cost_value + ((g_b_drop_YDF**2-(g_each_degree_point_YDF(i, i_GD, 2)-g_a_drop_YDF)**2-g_each_degree_point_YDF(i, i_GD, 1)**2)**2)/g_num_slice_for_GDalgo_YDF/4.0
                    endif
                enddo
                do i_GD = 1, g_num_slice_for_GDalgo_YDF
                    if(g_each_degree_point_YDF(i, i_GD, 1) .gt. g_first_deltaQ) then
                        g_a_drop_temp = g_a_drop_YDF
                        g_b_drop_temp = g_b_drop_YDF
                        g_a_drop_YDF = g_a_drop_temp - g_alpha_GD*(g_b_drop_temp**2-(g_each_degree_point_YDF(i, i_GD, 2)-g_a_drop_temp)**2-g_each_degree_point_YDF(i, i_GD, 1)**2)*(g_each_degree_point_YDF(i, i_GD, 2)-g_a_drop_temp)/g_num_slice_for_GDalgo_YDF
                        g_b_drop_YDF = g_b_drop_temp - g_alpha_GD*((g_b_drop_temp**2-(g_each_degree_point_YDF(i, i_GD, 2)-g_a_drop_temp)**2-g_each_degree_point_YDF(i, i_GD, 1)**2))*g_b_drop_temp/g_num_slice_for_GDalgo_YDF
                    endif
                enddo
                if (abs(g_cost_value - g_cost_temp) .le. 0.0000001) then
                    g_check_iter = 0
                endif
                g_cost_temp = g_cost_value
                g_cost_value = 0.0
            enddo
            !g_slope_temp = (g_a_drop_YDF - g_Drop_bottom)/sqrt(g_b_drop_YDF**2 - (g_Drop_bottom - g_a_drop_YDF)**2)
            g_slope_temp = (g_a_drop_YDF - initUcell(3)*gap(3)/2)/sqrt(g_b_drop_YDF**2 - (initUcell(3)*gap(3)/2 - g_a_drop_YDF)**2)
            g_contactangles_GD_YDF = 90 + atan(g_slope_temp)*180/pi
            if(isnan(g_contactangles_GD_YDF) .or. g_contactangles_GD_YDF .ge. 173) then 
                g_contactangles_GD_YDF = 180
            endif

            g_cost_value = 0.0
            g_cost_temp = 0.0
            g_check_iter = 1
            g_m = 0
            g_a_drop_YDF_top = -6.347737
            g_b_drop_YDF_top = 1.0
            g_alpha_GD = 0.001
            start_slice = int((g_Drop_top-g_Drop_bottom)/g_gap_slice)+1
            end_slice = start_slice - g_num_slice_for_GDalgo_YDF_top
            if(end_slice .lt. 0) then
                end_slice = 0
            endif
            do while (g_check_iter .eq. 1) 
                g_m = g_m + 1
                do i_GD = start_slice, end_slice, -1
                    !write(*, *) g_each_degree_point_YDF(i, i_GD, 1)
                    if(g_each_degree_point_YDF(i, i_GD, 1) .ge. g_first_deltaQ) then
                        g_cost_value = g_cost_value + ((g_b_drop_YDF_top**2-(g_each_degree_point_YDF(i, i_GD, 2)-g_a_drop_YDF_top)**2-g_each_degree_point_YDF(i, i_GD, 1)**2)**2)/g_num_slice_for_GDalgo_YDF_top/4.0
                    endif
                enddo
                !pause
                do i_GD = start_slice, end_slice, -1
                    if(g_each_degree_point_YDF(i, i_GD, 1) .ge. g_first_deltaQ) then
                        g_a_drop_temp = g_a_drop_YDF_top
                        g_b_drop_temp = g_b_drop_YDF_top
                        g_a_drop_YDF_top = g_a_drop_temp - g_alpha_GD*(g_b_drop_temp**2-(g_each_degree_point_YDF(i, i_GD, 2)-g_a_drop_temp)**2-g_each_degree_point_YDF(i, i_GD, 1)**2)*(g_each_degree_point_YDF(i, i_GD, 2)-g_a_drop_temp)/g_num_slice_for_GDalgo_YDF_top
                        g_b_drop_YDF_top = g_b_drop_temp - g_alpha_GD*((g_b_drop_temp**2-(g_each_degree_point_YDF(i, i_GD, 2)-g_a_drop_temp)**2-g_each_degree_point_YDF(i, i_GD, 1)**2))*g_b_drop_temp/g_num_slice_for_GDalgo_YDF_top
                    endif
                enddo
                if (abs(g_cost_value - g_cost_temp) .le. 0.0000001) then
                    g_check_iter = 0
                endif
                g_cost_temp = g_cost_value
                g_cost_value = 0.0
            enddo
            
            g_slope_temp = (g_a_drop_YDF_top - (g_Drop_top-g_gap_slice))/sqrt(g_b_drop_YDF_top**2 - ((g_Drop_top-g_gap_slice) - g_a_drop_YDF_top)**2)
            !write(*, *) "g_a_drop_YDF = ", g_a_drop_YDF, "g_b_drop_YDF = ", g_b_drop_YDF
            !write(*, *) "g_a_drop_YDF_top  = ", g_a_drop_YDF_top, "g_b_drop_YDF_top = ", g_b_drop_YDF_top
            !write(*, *) g_Drop_top, "  ", g_slope_temp
            !pause
            if(g_slope_temp .le. 0) then
                g_contactangles_GD_YDF_top = 90 + abs(atan(g_slope_temp)*180/pi)
            else
                g_contactangles_GD_YDF_top = 90 - atan(g_slope_temp)*180/pi
            endif
            if(isnan(g_contactangles_GD_YDF_top)) then
                g_contactangles_GD_YDF_top = 180
            endif
            
            
  ! 利用线性回归算法对求出液滴的切线----------y=ax+b-----a=a_drop_line  b=b_drop_line-------------------     
            !g_cost_value = 0.0
            !g_cost_temp = 0.0
            !g_check_iter = 1
            !g_m = 0
            !g_a_drop_line = 1.0
            !g_b_drop_line = 1.0
            !g_alpha_GD = 0.001
            !do while (g_check_iter .eq. 1) 
            !    g_m = g_m + 1
            !    do i_GD = 1, g_num_slice_for_GDalgo_line
            !        if(g_each_degree_point_YDF(i, i_GD, 1) .gt. g_first_deltaQ) then
            !            g_cost_value = g_cost_value + ((g_each_degree_point_YDF(i, i_GD, 1)*g_a_drop_line+g_b_drop_line-g_each_degree_point_YDF(i, i_GD, 2))**2)/g_num_slice_for_GDalgo_line/2.0
            !        endif
            !    enddo
            !    
            !    do i_GD = 1, g_num_slice_for_GDalgo_line
            !        if(g_each_degree_point_YDF(i, i_GD, 1) .gt. g_first_deltaQ) then
            !            g_a_drop_temp = g_a_drop_line
            !            g_b_drop_temp = g_b_drop_line
            !            g_a_drop_line = g_a_drop_temp - g_alpha_GD*(g_each_degree_point_YDF(i, i_GD, 1)*g_a_drop_line+g_b_drop_line-g_each_degree_point_YDF(i, i_GD, 2))*g_each_degree_point_YDF(i, i_GD, 1)/g_num_slice_for_GDalgo_line
            !            g_b_drop_line = g_b_drop_temp - g_alpha_GD*(g_each_degree_point_YDF(i, i_GD, 1)*g_a_drop_line+g_b_drop_line-g_each_degree_point_YDF(i, i_GD, 2))/g_num_slice_for_GDalgo_line
            !        endif 
            !    enddo
            !    if (abs(g_cost_value - g_cost_temp) .le. 0.0000001 .or. g_m .eq. 50000) then
            !        g_check_iter = 0
            !    endif
            !    g_cost_temp = g_cost_value
            !    g_cost_value = 0.0
            !enddo
            !g_contactangles_GD_line = 180 - atan(g_a_drop_line)*180/pi
            
            
            !g_each_CA(i) = g_contactangles_GD
            g_each_CA_YDF(i) = g_contactangles_GD_YDF
            g_each_CA_YDF_top(i) = g_contactangles_GD_YDF_top
            !g_each_CA_line(i) = g_contactangles_GD_line
            
            write(91,'(4f16.6)') (i-0.5)*g_delta_degree/360.0, g_contactangles_GD_YDF!, g_contactangles_GD_line, g_contactangles_GD
        enddo
        write(g_output_id, '(i6, 16f15.6)') stepCount, (g_each_CA_YDF(j_RDF), j_RDF=1, int(360.0/g_delta_degree), g_circle_output_delta)!, (g_each_CA(j_RDF), j_RDF=1, int(360.0/g_delta_degree), g_circle_output_delta)
        
        g_CL_A = 0.0  !contact line for bottom
        g_CL_R = 0.0  !contact line for top
        g_CA_bottom = 0.0
        g_CA_top = 0.0
        do i = 1, g_num_degree
            g_CL_A = g_CL_A + (g_each_degree_point_YDF(i, 1, 1) + g_Center_aver_X)/g_num_degree
            g_CL_R = g_CL_R + (g_each_degree_point_YDF(i, (g_Drop_top/g_gap_slice), 1) + g_Center_aver_X)/g_num_degree
            g_CA_bottom = g_CA_bottom + real((g_each_CA_YDF(i))/g_num_degree)
            g_CA_top = g_CA_top + real((g_each_CA_YDF_top(i))/g_num_degree)
        enddo
        
        if(stepCount .eq. g_outputFlu) then
            g_CL_A_prev = g_CL_A
            g_CL_R_prev = g_CL_R
        endif
        !g_CL_A = g_each_degree_point_YDF(1, 1, 1) + g_Center_aver_X
        !g_CL_R = -g_each_degree_point_YDF(int(180.0/g_delta_degree+1), 1, 1) + g_Center_aver_X
        g_CL_A_Velo = (g_CL_A - g_CL_A_prev)/(deltaT*g_outputFlu)
        g_CL_A_Velo_sum = g_CL_A_Velo
        g_CL_R_Velo = (g_Cl_R - g_Cl_R_prev)/(deltaT*g_outputFlu)
        g_CL_R_Velo_sum = sqrt(g_CL_R_Velo**2+wall_velo**2)
        write(g_output_CL_id, '(i6, 4f15.6)') stepCount, g_CL_A, g_CL_A_Velo, g_CL_R, g_CL_R_Velo
        g_CL_A_prev = g_CL_A
        g_CL_R_prev = g_CL_R
        write(2018, '(i6,8f16.6)') stepCount, g_CA_bottom, g_CL_A, g_CL_A_Velo, g_CL_A_Velo_sum, g_CA_top, g_CL_R, g_CL_R_Velo, g_CL_R_Velo_sum
        
        g_num_steps = 0
        g_Center_aver_X = 0.0
        g_Center_aver_Y = 0.0
        g_num_cell_particle(:, :, :) = 0
        g_YDF_sum_coor_aver(:, :, :) = 0.0
        g_density_cell(:, :, :) = 0.0
        g_each_degree_point(:, :, :) = 0.0
        g_each_degree_point_YDF(:, :, :) = 0.0
        
        !call find_bottomAndtop
        
        close(91)
    elseif(icord .eq. 3) then 
        deallocate(g_num_cell_particle)
        deallocate(g_density_cell)
        deallocate(g_each_degree_point)
        deallocate(g_each_degree_point_YDF)
        deallocate(g_YDF_sum_coor_aver)
        deallocate(g_YDF_sum_num_parti)
        deallocate(g_YDF_sum_coor)
        close(g_output_id)
        close(g_output_CL_id)
        close(2018)
    endif
    
    return 
    end
    
    
    subroutine find_bottomAndtop
   
    implicit none
    include 'dpdflow.h'
    
    integer n, i, j
    real*8 g_CA_drop_r_Z(nDpEnd-nWallAtom), g_temp
    
    g_Drop_bottom = 0.0
    g_Drop_top = 0.0
    do n = nWallAtom + 1 , nDpEnd
        g_CA_Drop_r_Z(n - nWallAtom) = r(n,3)
    enddo
    do i = nWallAtom + 2 , nDpEnd   
        do j = i , nWallAtom + 2 , -1
            if (g_CA_Drop_r_Z(j - nWallAtom) .lt. g_CA_Drop_r_Z(j - nWallAtom - 1)) then
                g_temp = g_CA_Drop_r_Z(j - nWallAtom)
                g_CA_Drop_r_Z(j - nWallAtom) = g_CA_Drop_r_Z(j - nWallAtom - 1)
                g_CA_Drop_r_Z(j - nWallAtom - 1) = g_temp
            else
                goto 34
            endif    
        enddo
34      continue  
    enddo
    do n = 1 , 10
        g_Drop_bottom = g_Drop_bottom + g_CA_Drop_r_Z(n)/10.0 
    enddo
    do n = nDpEnd - nWallAtom - 50 + 1  ,  nDpEnd - nWallAtom
        g_Drop_top = g_Drop_top + g_CA_Drop_r_Z(n)/50.0  
    enddo

    return
    end