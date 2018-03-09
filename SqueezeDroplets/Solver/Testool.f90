
    subroutine ouptParticleSituation
    ! To output the coordination, velocity and acceleration of every particle.
    
    implicit none
    include 'dpdflow.h'
    integer ifirst
    character*8 stc
    integer n,i
    
    write(stc,'(i8)')stepCount
    if(stepCount<10) then
        ifirst=index(stc,' ')+7
    elseif(stepCount<100) then
        ifirst=index(stc,' ')+6
    elseif(stepCount<1000) then
        ifirst=index(stc,' ')+5
    elseif(stepCount<10000) then
        ifirst=index(stc,' ')+4
    elseif(stepCount<100000) then
        ifirst=index(stc,' ')+3
    elseif(stepCount<1000000) then
        ifirst=index(stc,' ')+2   
    elseif(stepCount<10000000) then
        ifirst=index(stc,' ')+1 
    endif
    
    
	open(201, file = './instantState/Particles_'//stc(ifirst:)//'.plt')    
    write(201,'(''TITLE = "All Particles Coordinates"'')')
    write(201,'(''Variables= "X","Y","Z","Vx","Vy","Vz","ax","ay","az"'')')
	
    if (nWallAtom.ge. 1 .and. runId.ne.0.and. runId.ne.6) then
    !if (nWallAtom.ge. 1 .and. runId.eq.0.and. runId.ne.6) then
        write(201,'(''ZONE'')')
        write(201,'(''T="WallPariticles"'')')
        write(201,'(''I= '',i8,'',F = POINT'')') nWallAtom
        do n = 1, nWallAtom
	        write(201,'(9f15.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
        enddo
    endif

    if(nDpEnd>nWallAtom) then
        write(201,'(''ZONE'')')
        write(201,'(''T="BubblePariticles"'')') 
        write(201,'(''I= '',i8,'',F = POINT'')') nDpEnd-nWallAtom
       do n = nWallAtom + 1, nDpEnd
	      write(201,'(9f15.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
       enddo
    endif

    if(nChain .gt. 0) then
        write(201,'(''ZONE'')')
        write(201,'(''T="ChainPariticles"'')')
        write(201,'(''I= '',i8,'',F = POINT'')') nChainend-nDpEnd
       do n = nDpEnd + 1, nChainend
	      write(201,'(9f15.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
       enddo
       close(203)
    endif

 !   write(201,'(''ZONE'')')
 !   write(201,'(''T="FluidPariticles"'')')
 !   write(201,'(''I= '',i8,'',F = POINT'')') nAtom-nChainend
	!do n = nChainend + 1, nAtom
 !       write(201,'(9f15.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
 !   enddo


    write(201,'(''TEXT X=1 Y=90 T="   nAtom    nDp    Steps\\n'',3i8,'' "'')') &
                                     nAtom,nDpEnd-nWallAtom,stepCount
    close(201)
    write(*,*) 'We have output the all particles location | stepCount:',stepCount
 !pause
    
    end

!! Added by linyuqing

! --------------------------------------------------------------------------------------------------------------------
!subroutine ouptParticleGap(icord)
!    ! To output the coordination, velocity and acceleration of every particle.
!    
!    implicit none
!    include 'dpdflow.h'
!    integer n,i,icord
!    
!    if (icord .eq. 0) then
!	    open(202, file = './data/GapParticles.plt')    
!        write(202,'(''TITLE = "All Particles Coordinates"'')')
!        write(202,'(''Variables= "X","Y","Z","Vx","Vy","Vz","ax","ay","az"'')')
!    endif
!    
!    write(202,'(''ZONE'')')
!    write(202,'(''T="'',i8,''"'')') stepCount
!    write(202,'(''I= '',i8,'',F = POINT'')') nAtom
!	do n = 1, nAtom
!        write(202,'(9f10.4)') (r(n,i), i = 1,3),(rv(n,i), i = 1,3), (ra(n,i), i = 1,3)
!    enddo
!
!    write(202,'(''TEXT X=1 Y=90 T="   nAtom\\n'',i8,'' "'')') nAtom
!    
!    if (icord .eq. -1) close(202)
!
! !pause
!    
!    end
!
!!! Added by linyuqing
!
!! --------------------------------------------------------------------------------------------------------------------
!subroutine ouptSameParticlesAvg(icord)
!    ! To output the coordination, velocity and acceleration of every particle.
!    
!    implicit none
!    include 'dpdflow.h'
!    integer n,i,icord
!    real*8 Cr(3),Crv(3),Cra(3)
!    
!    if (icord .eq. 0) then 
!        open(205, file = './data/Fluid.plt')  
!        write(205,'(''Variables= "stepCount","CX","CY","CZ","CVx","CVy","CVz","Cax","Cay","Caz"'')')
!        write(205,'(''ZONE'')')
!        write(205,'(''F = POINT'')') 
!        if(nDp.eq. 1.and.nDpEnd>nWallAtom) then
!	        open(203, file = './data/Bubble.plt')  
!            write(203,'(''Variables= "stepCount","CX","CY","CZ","CVx","CVy","CVz","Cax","Cay","Caz"'')')
!            write(203,'(''ZONE'')')
!            write(203,'(''F = POINT'')') 
!        endif
!        if(nChain .gt. 0) then
!            open(204, file = './data/Chain.plt')  
!            write(204,'(''Variables= "stepCount","CX","CY","CZ","CVx","CVy","CVz","Cax","Cay","Caz"'')')
!            write(204,'(''ZONE'')')
!            write(204,'(''F = POINT'')') 
!        endif
!    endif
!    
!    do i =1,3
!        Cr(i)=0.
!        Crv(i)=0.
!        Cra(i)=0.
!    enddo
!	do n = nChainend+1, nAtom
!        do i =1,3
!            Cr(i)=Cr(i)+r(n,i)
!            Crv(i)=Crv(i)+rv(n,i)
!            Cra(i)=Cra(i)+ra(n,i)
!        enddo
!    enddo
!    do i =1,3
!        Cr(i)=Cr(i)/(nAtom-nChainend)
!        Crv(i)=Crv(i)/(nAtom-nChainend)
!        Cra(i)=Cra(i)/(nAtom-nChainend)
!    enddo
!    write(205,'(i8,9f10.4)') stepCount,(Cr(i), i = 1,3),(Crv(i), i = 1,3), (Cra(i), i = 1,3)
!    if (icord .eq. -1)  close(205)
!    
!    if(nDp.eq. 1.and.nDpEnd>nWallAtom) then
!        do i =1,3
!            Cr(i)=0.
!            Crv(i)=0.
!            Cra(i)=0.
!        enddo
!	    do n = nWallAtom+1, nDpEnd
!            do i =1,3
!                Cr(i)=Cr(i)+r(n,i)
!                Crv(i)=Crv(i)+rv(n,i)
!                Cra(i)=Cra(i)+ra(n,i)
!            enddo
!        enddo
!        do i =1,3
!            Cr(i)=Cr(i)/(nDpEnd-nWallAtom)
!            Crv(i)=Crv(i)/(nDpEnd-nWallAtom)
!            Cra(i)=Cra(i)/(nDpEnd-nWallAtom)
!        enddo
!        write(203,'(i8,9f10.4)') stepCount,(Cr(i), i = 1,3),(Crv(i), i = 1,3), (Cra(i), i = 1,3)
!        if (icord .eq. -1)  close(203)
!    endif
!
!    if(nChain .gt. 0) then
!        do i =1,3
!            Cr(i)=0.
!            Crv(i)=0.
!            Cra(i)=0.
!        enddo
!	    do n = nDpEnd+1, nChainend
!            do i =1,3
!                Cr(i)=Cr(i)+r(n,i)
!                Crv(i)=Crv(i)+rv(n,i)
!                Cra(i)=Cra(i)+ra(n,i)
!            enddo
!        enddo
!        do i =1,3
!            Cr(i)=Cr(i)/(nChainend-nDpEnd)
!            Crv(i)=Crv(i)/(nChainend-nDpEnd)
!            Cra(i)=Cra(i)/(nChainend-nDpEnd)
!        enddo
!        write(204,'(i8,9f10.4)') stepCount,(Cr(i), i = 1,3),(Crv(i), i = 1,3), (Cra(i), i = 1,3)
!        if (icord .eq. -1)  close(204)
!    endif  
!    
! !pause
!    
!    end
!
!    
    
    subroutine OutputDropState(icord)
    
    implicit none
    include 'dpdflow.h'
    
    integer State,n,i,nC,k,icord
    real*8 CDCoor(3),dC,dd
    
    integer nslice,ninner,Nth,Theta,checktemp
    real*8 bottomZ,slicegap,CircleCentX,CircleCentY,highZ,simuRad,ZZ,CA,Bottomtemp,tempCA,errorCA
    real*8 , allocatable :: vectorDropp(:,:),maxdistances(:,:),averagedistances(:)
    integer ifirst
    character*8 stc
    
    State = 1
    if(icord .eq. 0) then
        open(912, file = './data/DropState.plt', status = 'unknown') 
        write(912,'(''Variables= "stepCount","CDCoorZ","CDCoorY","CDCoorX","nC","State","Passk"'')')
        write(912,'(''ZONE'')')
        write(912,'(''F = POINT'')')
        
        open(913, file = './data/RedrawDrop_ContactAngles.plt', status = 'unknown') 
        write(913,'(''Variables= "stepCount","DropRadi","Radi(bottom)","HighZ","ContactAngles","CA_Bar"'')')
        write(913,'(''ZONE'')')
        write(913,'(''F = POINT'')')
        
        simuRadsum = 0.0
        averdistancessum = 0.0
        highZsum = 0.0
        CAsum = 0.0
        stepsum = 0
        CAmin = 180.0
        CAmax = 0.0
        CAerrorbar = 0.0
        CAtemp = 0.0
    elseif(icord .eq. 1) then

        if(mod(stepCount,100) .eq.0 ) then 
            write(stc,'(i8)')stepCount
            if(stepCount<10) then
                ifirst=index(stc,' ')+7
            elseif(stepCount<100) then
                ifirst=index(stc,' ')+6
            elseif(stepCount<1000) then
                ifirst=index(stc,' ')+5
            elseif(stepCount<10000) then
                ifirst=index(stc,' ')+4
            elseif(stepCount<100000) then
                ifirst=index(stc,' ')+3
            elseif(stepCount<1000000) then
                ifirst=index(stc,' ')+2   
            elseif(stepCount<10000000) then
                ifirst=index(stc,' ')+1 
            endif
        
            open(914, file = './SimuDrop/SimuDrop_'//stc(ifirst:)//'.plt')   
            write(914,'(''Variables= "radi","highz","stepCount"'')')
            write(914,'(''ZONE'')')
            write(914,'(''F = POINT'')')  
        endif
! output RedrawDrop_ContactAngles.plt Begin---------------------------------------------------------------------------------- 
        bottomZ = r(nWallAtom+1,3)
        do n = nWallAtom+2 ,nDpEnd
            if(r(n,3) .lt. bottomZ .and. r(n,3) .gt. r(1,3)) then
                bottomZ = r(n,3)
            endif
        enddo
        slicegap = 0.2
        nslice = int(2*RdsDp*2/slicegap)
        Bottomtemp = 1000.0   
        tempCA = 1000.0 
        errorCA = 100.0
        checktemp = 0
        allocate(averagedistances(nslice))
1993    continue
        do k = 1 , nslice
            averagedistances(k) = 0.0
            ninner = 0
            CircleCentX = 0
            CircleCentY = 0  
            allocate(vectorDropp(nAtom,5))
            do n = nWallAtom+1 ,nDpEnd
                if((r(n,3) .lt. bottomZ+k*slicegap) .and. (r(n,3) .ge. bottomZ+(k-1)*slicegap)) then  !Find out particles in this slice
                    ninner = ninner + 1
                    vectorDropp(ninner,3) = n     !record N th paticle
                    CircleCentX = r(n,1) + CircleCentX
                    CircleCentY = r(n,2) + CircleCentY
                endif
            enddo
            if(ninner .gt. 0) then
                CircleCentX = CircleCentX/ninner !Find out the center coordinate of this slice
                CircleCentY = CircleCentY/ninner
            endif   
            do n = 1 , ninner
                Nth = vectorDropp(n,3)
                vectorDropp(n,1) = r(Nth,1) - CircleCentX
                vectorDropp(n,2) = r(Nth,2) - CircleCentY
                vectorDropp(n,4) = sqrt(vectorDropp(n,1)*vectorDropp(n,1)+vectorDropp(n,2)*vectorDropp(n,2))  !The distances of particle and center
                if(r(Nth,2) .ge. CircleCentY) then   
                    vectorDropp(n,5) = acos(vectorDropp(n,1)/vectorDropp(n,4))*180.0/3.1415926535898    !The angles of the Drop particle's Vector and X-axis (бу) (.le. 180бу)
                else
                    vectorDropp(n,5) =360.0 - acos(vectorDropp(n,1)/vectorDropp(n,4))*180.0/3.1415926535898    ! ( .gt. 180бу )
                endif
            enddo
            allocate(maxdistances(nslice,100))
            do Theta = 1 , 50
                maxdistances(k,Theta) = 0
                do n = 1 , ninner
                    if(k .eq. 1 .and. vectorDropp(n,5) .ge. (Theta-1)*360.0/50.0 .and. vectorDropp(n,5) .lt. Theta*360.0/50.0 .and. vectorDropp(n,4) .ge. maxdistances(k,Theta) .and. vectorDropp(n,4) .le. Bottomtemp)then
                        maxdistances(k,Theta) = vectorDropp(n,4)
                    endif
                    if(k .ne. 1 .and. vectorDropp(n,5) .ge. (Theta-1)*360.0/50.0 .and. vectorDropp(n,5) .lt. Theta*360.0/50.0 .and. vectorDropp(n,4) .ge. maxdistances(k,Theta))then
                        maxdistances(k,Theta) = vectorDropp(n,4)
                    endif
                enddo
                averagedistances(k) = averagedistances(k) + maxdistances(k,Theta)/50.0
            enddo
            deallocate(vectorDropp,maxdistances)
            if(mod(stepCount,100) .eq.0 ) then
                write(914,'(2f12.6,i10)') averagedistances(k),((k - 1)*slicegap + k*slicegap)/2,stepCount
            endif
        enddo
         do k = 1 , nslice
            if(averagedistances(k) .eq. 0) then
                if(k .ne. 1)then
                    highZ = ((k - 1)*slicegap + (k - 2)*slicegap)/2
                else
                    write(*,*) "The Droplet out of the plant, the contact angles large than 180бу"
                    goto 2000
                endif    
                goto 1999
            endif
        enddo
1999    continue
        simuRad = sqrt(averagedistances(1)*averagedistances(1) + highZ*highZ)/2.0/cos(atan(averagedistances(1)/highZ))
        ZZ = highZ - simuRad
        CA = 90 + asin(ZZ/simuRad)*180.0/3.1415926535898
        Bottomtemp = averagedistances(1)
        errorCA = CA - tempCA
        !if(abs(errorCA) .ge. 1.0) then
        if(checktemp .eq. 0 .and. CA .lt. 90)then
            !tempCA = CA
            checktemp = 1
            goto 1993
        endif
        !if(CA .lt. CAmin) then
        !    CAmin = CA
        !endif
        !if(CA .gt. CAmax) then
        !    CAmax = CA
        !endif
        stepsum = stepsum + 1
        CAtemp(stepsum) = CA
        simuRadsum = simuRadsum + simuRad
        averdistancessum = averdistancessum + averagedistances(1)
        highZsum = highZsum + highZ 
        CAsum = CAsum + CA
2000    continue      

        if(mod(stepCount,100) .eq.0 ) then
            close(914)
        endif
        deallocate(averagedistances)
! output RedrawDrop_ContactAngles.plt End----------------------------------------------------------------------------------       

    elseif(icord .eq. 2) then
! output DropState.plt begin--------------------------------------------------
        do i = 1 , 3
            CDCoor(i) = 0
        enddo
        
        do i =1,3
            do n = nWallAtom+1, nDpEnd
                CDCoor(i) = CDCoor(i) + r(n,i)
            enddo
            CDCoor(i) = CDCoor(i)/(nDpEnd-nWallAtom)
        enddo
        
        nC= nWallAtom+1
        dC= 0.
        do i=1,3
            dC=dC+(r(nC,i)-CDCoor(i))**2
            dC=sqrt(dC)
        enddo
        do n = nWallAtom+2, nDpEnd !find out the center of the droplet
            dd=0.
            do i=1,3
                dd=dd+(r(n,i)-CDCoor(i))**2
                dd=sqrt(dd)
            enddo
            if(dd .lt. dC) then 
                dC=dd
                nC=n    ! record the particle that closest the center
            endif
        enddo
        
        k=0
        do n = nWallAtom+1, nDpEnd
            !if(r(n,3) .ge. (-regionH(3)/3) .or. abs(r(n,1)) .ge. (RdsDp*1.3) .or. abs(r(n,2)) .ge. (RdsDp*1.3) ) then
            if(abs(r(n,1)) .ge. 1.2*simuRadsum/stepsum .or. abs(r(n,2)) .ge. 1.2*simuRadsum/stepsum ) then            
                k=k+1
            endif
        enddo
        if(k .ge. int((nDpEnd-nWallAtom)/10)) then 
            State = 0
        endif
        write(912,'(i6,3f16.6,3i10)') stepCount,CDCoor(3),CDCoor(2),CDCoor(1),nC,State,k
! output DropState.plt end----------------------------------------------------------------------------------------------------------
        
!average the droplet's paraments begain--------------------------------------------------------------------------------------------       
        do i = 1 , stepsum
             CAerrorbar = CAerrorbar + (CAtemp(stepsum)-CAsum/stepsum)*(CAtemp(stepsum)-CAsum/stepsum)
        enddo
        CAerrorbar = CAerrorbar/stepsum
        CAerrorbar = sqrt(CAerrorbar)        
        write(913,'(i10,5f12.6)') stepCount,simuRadsum/stepsum,averdistancessum/stepsum,highZsum/stepsum,CAsum/stepsum,CAerrorbar
        simuRadsum = 0.0
        averdistancessum = 0.0
        highZsum = 0.0
        CAsum = 0.0
        stepsum = 0        
        CAmin = 180.0
        CAmax = 0.0
        CAerrorbar = 0.0
        CAtemp = 0.0
!average the droplet's paraments end--------------------------------------------------------------------------------------------        
    elseif(icord .eq. 3) then 
        close(912) 
        close(913)
    endif
    
    end
    

!! !! Added by linyuqing -------------------------------------------------------------------------------------------
!    
!! subroutine DiffusionCompute(icord)
!!     ! To output the coordination, velocity and acceleration of every particle.
!    
!!     implicit none
!!     include 'dpdflow.h'
!
!!     integer n,i,icord
!!     real * 8 Dmsd,time, DiffStartStep !Dgk,
!    
!!     real * 8, pointer::ZeroSate(:,:)
!!     common / Diff / ZeroSate
!!     common / Diffstep / DiffStartStep
!    
!! !    DiffintStep,DiffStep
!    
!    
!!     if (icord .eq. 0) then
!        
!!         allocate (ZeroSate(nAtom - nWallAtom,NDIM*2))
!!         do n = nWallAtom+1,nAtom
!!             do i= 1, NDIM
!!                 ZeroSate(n - nWallAtom,i) = r(n,i)+ ncc(n,i) * initUcell(i)*gap(i)
!! !                ZeroSate(n - nWallAtom,i+NDIM) = rv(n,i)
!!             enddo
!!         enddo
!!         DiffStartStep = stepCount
!!         time = 0.
!!         Dmsd = 0.
!! !        Dgk = 0.
!        
!!         open(401, file = './data/CG/Diffusivity.plt')    
!!         write(401,'(''TITLE = "Diffusivity('',i6,'')"'')') stepCount
!!         write(401,'(''Variables= "t","D(MSD)"'')')  !,"D(GK)"
!!         write(401,'(f10.4,2E15.7)') time,Dmsd !,Dgk
!!     endif
!
!!     if (icord .eq. 1) then
!!         time = (stepCount-DiffStartStep)*deltaT
!!         Dmsd = 0.
!! !        Dgk = 0.
!        
!!         do n = nWallAtom+1,nAtom
!!             do i= 1, NDIM
!!                 Dmsd = Dmsd + (r(n,i)+ ncc(n,i) * initUcell(i)*gap(i) - ZeroSate(n - nWallAtom,i))**2
!! !                Dgk = Dgk + ZeroSate(n - nWallAtom,i+NDIM)*rv(n,i)
!!             enddo
!!         enddo
!!         Dmsd = Dmsd /(nAtom - nWallAtom)
!! !        Dgk = Dgk /(nAtom - nWallAtom)
!!         write(401,'(f10.4,2E15.7)') time,Dmsd !,Dgk
!!     endif
!    
!!     if (icord .eq. 2) then
!!         deallocate (ZeroSate)
!!         close(401)
!!     endif
!    
!    
!!     end
!
!!! Added by linyuqing -------------------------------------------------------------------------------------   
!
!subroutine Barostat
!    ! To apply barostat to get constant Pressure
!    
!    implicit none
!    include 'dpdflow.h'
!
!    real*8 mu,P
!    
!    real*8 alpha,cc,d
!    
!    integer n,k,i,j,l,rhoSize(3),rhoSizeH(3),ic(3)
!    real*8 T,rho,rtemp,rhoV
!    
!    real*8, allocatable :: rhoGrid(:,:,:,:)
!    
!    alpha=0.1
!    cc=4.16
!    d=18.
!        
!	if(stepCount .eq. stepEquil) then
!    open(501, file = './data/Barostat/TrhoP.plt')    
!    write(501,'(''TITLE = "EOS"'')')
!    write(501,'(''Variables= "time","T","rho","P","mu","tau"'')')
!	endif
!    !get the liquid temperature    
!    T =0.0
!    do n = nChainend + 1, nAtom
!        do k =1,NDIM
!            T= T + rv(n,k)**2
!        enddo
!    enddo
!    T = T/(3.0*(nAtom - nChainend))
!    
!    !get the liquid density
!    do i=1,3
!    rhoSize(i) = int(region(i)/(2.*rCut))
!	if(mod(rhoSize(i),2) .eq. 0) rhoSize(i) = rhoSize(i) + 1
!	rhoSizeH(i) = (rhoSize(i) - 1)/2
!    enddo
!       
!    allocate(rhoGrid(-rhoSizeH(1):rhoSizeH(1),-rhoSizeH(2):rhoSizeH(2),-rhoSizeH(3):rhoSizeH(3),1:3))
!    
!    do i= -rhoSizeH(1),rhoSizeH(1),1
!    	do j= -rhoSizeH(2),rhoSizeH(2),1
!    		do k= -rhoSizeH(3),rhoSizeH(3),1
!	        do l =1,3
!	            rhoGrid(i,j,k,l) =0.0
!	        enddo
!	        enddo
!        enddo       
!    enddo
!    !pause
!    do n = nWallAtom + 1, nAtom
!        do i=1,3
!       		rtemp = (abs(r(n,i))**2)/(r(n,i)*(2.*rCut))
!       		rtemp = rtemp + (rtemp/(abs(rtemp)))*0.5
!       		ic(i)=int(rtemp)                    
!        enddo
!
!        if (abs(ic(1)) .le. rhoSizeH(1) .and. abs(ic(2)) .le. rhoSizeH(2) .and. abs(ic(3)) .le. rhoSizeH(3)) then
!	        rhoGrid(ic(1),ic(2),ic(3),1) =  rhoGrid(ic(1),ic(2),ic(3),1) + 1.
!	    
!	        if (n .le. nDpEnd) then
!	            rhoGrid(ic(1),ic(2),ic(3),2) =  rhoGrid(ic(1),ic(2),ic(3),2) + 1.
!	        elseif (n .le. nChainend) then
!	            rhoGrid(ic(1),ic(2),ic(3),3) =  rhoGrid(ic(1),ic(2),ic(3),3) + 1.
!	        endif 
!        endif      
!    enddo
!    l=0
!    rho= 0.
!    do i= -rhoSizeH(1)+1,rhoSizeH(1)-1,1
!    	do j= -rhoSizeH(2)+1,rhoSizeH(2)-1,1
!    		do k= -rhoSizeH(3)+1,rhoSizeH(3)-1,1
!	        if (rhoGrid(i,j,k,2) .eq. 0.0 .and. rhoGrid(i,j,k,3) .eq. 0.0 ) then
!	            rho = rho + rhoGrid(i,j,k,1)/((2.*rCut)**3)
!	            l=l+1
!	        endif 
!	        enddo
!        enddo       
!    enddo
!    !pause
!    if(l.gt.0) rho = rho /l
!
!    !! The Berendsen Barostat
!    P = 2*alpha*alphaB*(rCut2**4)*(rho**3)+(alpha*alphaf-2*alpha*alphaB*(rCut2**4)*cc)*(rho**2)+T*rho+2*alpha*alphaB*(rCut2**4)*d
!    print*,mu,rho,P,(region(i),i=1,3),(regionH(i),i=1,3),'1'
!    if(stepCount .eq.  (JStep + stepEquil)) P0 = JP0
!    mu=1
!    if (mod(stepCount,BStep) .eq. 0 .and. (rho0 .ge.0 .or. stepCount .ge.JStep+stepEquil)) then
!        !print*,deltaT,tau,rho,P,P0
!        mu= (1-((deltaT/tau)*(P0-P)))**(0.3333333333)
!        !
!        !do while (mu .le. 0.9)
!        !    if ((P0-P) .gt. 0) then
!        !        tau = tau *2
!        !    else
!        !        tau = tau /2
!        !    endif
!        !    print*, "the rise time (tau) is too small, change it to be: ", tau
!        !    mu= 1-((deltaT/tau)*(P0-P))!)**(0.3333333333)
!        !enddo
!        !
!        !do while (mu .gt. 1.1)
!        !    if ((P0-P) .gt. 0) then
!        !        tau = tau /2
!        !    else
!        !        tau = tau *2
!        !    endif
!        !    print*, "the rise time (tau) is too great, change it to be: ", tau
!        !    mu= 1-((deltaT/tau)*(P0-P))!)**(0.3333333333)
!        !enddo
!        
!        print*,mu,rho,P,(region(i),i=1,3),(regionH(i),i=1,3),'2'
!        
!        !refresh the region
!        do n = nChainend + 1, nAtom
!        	do i=1,3
!            if(abs(r(n,i)) .gt. (regionH(i)*lpercent)) r(n,i) = r(n,i)*mu
!            enddo
!        enddo
!        
!        do i=1,3
!	        region(i) = region(i)*mu
!	        regionH(i) = region(i)/2.
!
!	    	cells(i) = int(region(i) / rCut)
!	    	if(cells(i) .lt. 1) cells(i) = 1
!    	enddo
!        
!        binvolm = region(1)*region(2)*region(3)/hsize
!        
!        maxList =  nAtom+cells(1)*cells(2)*cells(3)
!    
!    endif
!    print*,mu,rho,P,(region(i),i=1,3),(regionH(i),i=1,3),'3'
!
!    write(501,'(f6.2,2f10.4,f16.6,f10.4,f10.4)') timeNow,T,rho,P,mu,tau
!
!    if(stepCount .ge. stepLimit) close(501)
!
!    deallocate(rhoGrid)
!    end
!
!!! Added by linyuqing   
!
!
!subroutine BubbleSize
!    ! To get the bubble radius and centre point
!
!    implicit none
!    include 'dpdflow.h'
!    
!    integer i,j,k,l,m,n,Npb
!    real*8 BCr(3),theta,Ixx,Iyy,Izz,Ixz,BRadius
!    real*8 S1,S2,S3
!    
!    real*8, allocatable :: BPr(:,:)
!    
!    Npb = nDpEnd - nWallAtom
!    
!    allocate(BPr(Npb,3))
!    
!    do i =1,3
!        BCr(i) = 0.
!    enddo 
!    
!    Ixx = 0.
!    Iyy = 0.
!    Izz = 0.
!    Ixz = 0.
!
!    do n = nWallAtom+1, nDpEnd
!        do i =1,3
!        BCr(i) = BCr(i) + r(n,i) + ncc(n,i)*region(i) 
!        BPr(n - nWallAtom,i) = r(n,i) + ncc(n,i)*region(i) 
!        enddo
!    enddo
!    do i =1,3
!        BCr(i) = BCr(i) / (nDpEnd - nWallAtom)
!        do n =1,Npb
!            BPr(n,i) = BPr(n,i) - BCr(i)
!        enddo
!    enddo   
!    
!    do n =1,Npb
!        Ixx = Ixx + BPr(n,1)**2
!        Iyy = Iyy + BPr(n,2)**2
!        Izz = Izz + BPr(n,3)**2
!        Ixz = Ixz + BPr(n,3)*BPr(n,1)
!    enddo
!       
!    theta = 0.5*atan(2*Ixz/(Ixx-Izz))
!
!    S1 = Ixx*cos(theta)**2 + Izz*sin(theta)**2 + 2*Ixz*sin(theta)*cos(theta)
!
!    S2 = Ixx*sin(theta)**2 + Izz*cos(theta)**2 - 2*Ixz*sin(theta)*cos(theta)
!
!    S3 = Iyy
!    
!    
!    BRadius = (125*S1*S2*S3/(Npb**3))**(0.0666666666666)
!    
!    if(stepCount .eq. stepEquil) then
!    open(601, file = './data/RP.plt')    
!    write(601,'(''TITLE = "RP"'')')
!    write(601,'(''Variables= "time","Radius","Cx","Cy","Cz"'')')
!    endif
!    
!    !print*, BRadius,(BCr(i),i=1,3)
!    
!    write(601,'(f6.2,4f10.4)') timeNow,BRadius,(BCr(i),i=1,3)
!    
!    if(stepCount .ge. stepLimit)  close(601)
!
!
!
!
!end
!!! Added by linyuqing   