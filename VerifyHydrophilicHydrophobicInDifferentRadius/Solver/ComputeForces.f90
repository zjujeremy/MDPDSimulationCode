!----------------------------------------------------------
!   Compute forces acting on particles
!----------------------------------------------------------

    subroutine ComputeForces

    implicit none
    include 'dpdflow.h'

!   Parameters
!   Locals

    integer c, i, j1, j2, k, m1, m2, m1X, m1Y, m1Z, m2X, m2Y, m2Z, n,  &
            offset, tag1, tag2, cellZ
    integer	iofX(14), iofY(14), iofZ(14)
    real*8  f, fr, fcVal, fdVal, frVal, ftot, weight, rij, weightc,    &
            weightd, weightr, rr, rdvij, uVal, fstl, rangls, dn,       &
            gamma1, cigama1
    real*8  ri(NDIM), dr(NDIM), invWid(NDIM), shift(NDIM)

    integer, allocatable :: cellList(:)
    real,    allocatable :: sqrDist(:)

    real*8 fCV, fDP, fRD

!    integer nbSP(MAXATOM), nbDP(MAXATOM)
!    real*8 pct(MAXATOM)

!---MDPD used variables-----------------------------------------------
    integer j, m, tag(MAXATOM, 50), ntag(MAXATOM)
    real*8  mrho(MAXATOM), ddr(MAXATOM, 50, NDIM), drij(MAXATOM, 50),  &
            fcRval, weightp, weightcd(MAXATOM,50), fCR
    data    pi / 3.14159265358979 /
!---------------------------------------------------------------------

    data    iofX / 0,1,1,0,-1,0,1,1,0,-1,-1,-1, 0, 1 /
    data    iofY / 0,0,1,1, 1,0,0,1,1, 1, 0,-1,-1,-1 /
    data    iofZ / 0,0,0,0, 0,1,1,1,1, 1, 1, 1, 1, 1 /

    atomID = 0
    allocate ( cellList(maxList+nAtomCopy) )
    allocate ( sqrDist(nAtom) )
    sqrDist = 100.

    do k = 1, NDIM
        invWid(k) = cells(k) / region(k)
    enddo

    if(runId .eq. 5) invWid(NDIM) = (cells(NDIM) - 1) / region(NDIM)
    
    do n = nAtom+nAtomCopy+1, maxList+nAtomCopy
        cellList(n) = 0
    enddo

    do n = nStartAtom, nAtom+nAtomCopy
        do k = 1, NDIM
            ra(n,k)   = 0.
            raCV(n,k) = 0.
            raCR(n,k) = 0.
            raDP(n,k) = 0.
            raRD(n,k) = 0.
            raSP(n,k) = 0.
            rforce(n,k) = 0.
            rforce(n,k+NDIM) = 0.
        enddo
        ntag(n) = 0
        mrho(n) = 0.

!       nbSP(n) = 0
!       nbDP(n) = 0
    enddo

    do n = nStartAtom, nAtom + nAtomCopy
        c = (int((r(n,3) + regionH(3))*invWid(3))*cells(2) +           &
             int((r(n,2) + regionH(2))*invWid(2)))*cells(1) +          &
             int((r(n,1) + regionH(1))*invWid(1)) + nAtom + nAtomCopy + 1

        if(c .lt. 1) then
            print*, 'in cf, stepCount = ', stepCount
            print*, ' c = ', c, '  n = ', n
            print*, 'r(n,*) = ', (r(n,k), k = 1, 3)
            pause
        endif

        if(c .gt. maxList + nAtomCopy) then
            print*, ' c = ', c, '  n = ', n, 'stepCount = ', stepCount
            print*, 'r(n) =  ', (r(n,k),  k = 1, 3)
            pause
        endif

        cellList(n) = cellList(c)
        cellList(c) = n
    enddo

    uSum = 0.
    virSum = 0.

    if(runId .eq. 5) then
        cellZ = cells(3) - 1
    else
        cellZ = cells(3)
    endif

	do m1Z = 1, cellZ
        do m1Y = 1, cells(2)
            do m1X = 1, cells(1)
                m1 = ((m1Z - 1)*cells(2) + m1Y - 1)*cells(1) + m1X + nAtom + nAtomCopy

                do offset = 1, 14
                    m2X = m1X + iofX(offset)
                    shift(1) = 0

                    if(m2X .gt. cells(1)) then
                        m2X = 1
                        shift(1) =  region(1)
                    elseif(m2X .eq. 0) then
                        m2X = cells(1)
                        shift(1) = -region(1)
                    endif

                    m2Y = m1Y + iofY(offset)
                    shift(2) = 0

                    if(m2Y .gt. cells(2)) then
                        m2Y = 1
                        shift(2) =  region(2)
                    elseif( m2Y .eq. 0) then
                        m2Y = cells(2)
                        shift(2) = -region(2)
                    endif

                    m2Z = m1Z + iofZ(offset)
                    shift(3) = 0

                    if(runId .eq. 1 .or. runId .eq. 5) then
                        if(m2Z .eq. 0 .or. m2Z .gt. cells(3)) goto 15
                    else
                        if(m2Z .gt. cells(3)) then
                            m2Z = 1
                            shift(3) =  region(3)
                        elseif(m2Z .eq. 0) then
                            m2Z = cells(3)
                            shift(3) = -region(3)
                        endif
                    endif

                    m2 = ((m2Z - 1)*cells(2) + m2Y - 1)*cells(1) + m2X + nAtom + nAtomCopy

                    j1 = cellList(m1)
                    do while(j1 .gt. 0)
                        j2 = cellList(m2)
                        do while(j2 .gt. 0)
                            if(m1 .ne. m2 .or. j2 .lt. j1) then

                                do k = 1, NDIM
                                    dr(k) = r(j1,k) - r(j2,k) - shift(k)
                                enddo

                                rr = dr(1)*dr(1) + dr(2)*dr(2) + dr(3)*dr(3)

                                if(rr .lt. rrCut) then
                                    rij = sqrt(rr)

                                    rdvij = 0.
                                    do k = 1, NDIM
                                        dr(k) = dr(k)/rij
                                        rdvij = rdvij +                &
                                            dr(k)*(rv(j1,k) - rv(j2,k))
                                    enddo 

                                    weight  = 1. - rij/rCut
                                    weightc = weight

! restrict the conservative cut radius equal to 1. (DP--7/12/11)
!                                   weightc = 1. - rij
!                                   weight  = 1. - rij/rcut

! inclusion region for repulsion potential in MDPD
                                    if(rij .lt. rCut2) then
! weight function related to local density
                                        weightp  = 15./(2.*pi*r3Cut2)* &
                                                    (1-rij/rCut2)**2
                                        ntag(j1) = ntag(j1) + 1
                                        ntag(j2) = ntag(j2) + 1
                                        tag(j1,ntag(j1)) = j2
                                        tag(j2,ntag(j2)) = j1
                                        weightcd(j1,ntag(j1)) = 1. -   &
                                                            rij/rCut2
                                        weightcd(j2,ntag(j2)) =        &
                                                   weightcd(j1,ntag(j1))
                                        do k = 1, NDIM
                                            ddr(j1,ntag(j1),k) =  dr(k)
                                            ddr(j2,ntag(j2),k) = -dr(k)
                                        enddo
                                        drij(j1,ntag(j1)) = rij
                                        drij(j2,ntag(j2)) = rij
! local density parameter in MDPD
                                        mrho(j1) = mrho(j1) + weightp
                                        mrho(j2) = mrho(j2) + weightp
                                    endif

                                    if(weightc .lt. 0.d0) weightc=0.d0
                                    if(weight  .lt. 0.d0) weight =0.d0

! original version of DPD: W^C(r) = W^R(r) = 1 - r, W^D(r) = (1-r)**2
                                    weightd = weight*weight
                                    weightr = weight

! modified version of DPD: W^C(r) = 1-r, W^R = (1-r)**1/4, W^D(r) = (1-r)**1/2
!                                   weightd = sqrt(weight)
!                                   weightr = sqrt(weightd)

                                    cigama1 = cigamaF
                                    gamma1  = gammaF

                                    tag1 = tagLE(j1)
                                    tag2 = tagLE(j2)

! both j1 and j2 are wall particles
                                    fcVal  = 0.
                                    uVal   = 0.

! one of j1 or j2 is a wall particle and another is a chain(or drop) bead or solvent particle
 			             if((tag1 .le. nWallAtom .and. tag2 .gt. nWallAtom) .or. &
			                (tag1 .gt. nWallAtom .and. tag2 .le. nWallAtom)) then
                                if(wr_mark(tag1) .eq. 1 .or. wr_mark(tag2) .eq. 1 ) then
    !                                        call BellShape(rij, fcVal)
                                    fcVal   =  alphaw_top*weightc
                                    uVal    = -alphaw_top*(rij - 0.5*rr)
                                    cigama1 =  cigamaw
                                    gamma1  =  gammaw

                                    k  = min(tag1,tag2)             ! wall particle (DP--21/11/11)
			                        i  = max(tag1,tag2)             ! solvent or bead particle
                                    dn = dr(1)*wn(k,1) +           &
                                            dr(2)*wn(k,2) +           &
                                            dr(3)*wn(k,3)

                                    if(abs(dn) .le. wLayer .and.   &
                                            rr .lt. sqrDist(i)) then
                                        sqrDist(i) = rr
                                        atomID(i)  = k
                                    endif
                                    goto 10                                                         
                                else
    !                                        call BellShape(rij, fcVal)
                                    fcVal   =  alphawf*weightc
                                    uVal    = -alphawf*(rij - 0.5*rr)
                                    cigama1 =  cigamaw
                                    gamma1  =  gammaw

                                    k  = min(tag1,tag2)             ! wall particle (DP--21/11/11)
			                        i  = max(tag1,tag2)             ! solvent or bead particle
                                    dn = dr(1)*wn(k,1) +           &
                                            dr(2)*wn(k,2) +           &
                                            dr(3)*wn(k,3)

                                    if(abs(dn) .le. wLayer .and.   &
                                            rr .lt. sqrDist(i)) then
                                        sqrDist(i) = rr
                                        atomID(i)  = k
                                    endif
                                    goto 10
                                endif
                          endif

! both j1 and j2 are droplet particle
			             if(tag1 .gt. nWallAtom .and. tag2 .gt. nWallAtom .and. &
				            tag1 .le. nDpEnd    .and. tag2 .le. nDpEnd) then

                            if(int((tag1-nWallAtom-1)/nPDP) .eq.  &
                               int((tag2-nWallAtom-1)/nPDP)) then
                                fcVal= alphaf*weightc
                                uVal =-alphaf*(rij - 0.5*rr)
!                               nbSP(j1) = nbSP(j1) + 1
!                               nbSP(j2) = nbSP(j2) + 1
                            else
!                               if(DpSign(j1) .eq. 1 .or. DpSign(j2) .eq. 1) then
!                                   fcVal =  alphaf*weightc
!                                   uVal  = -alphaf*(rij - 0.5*rr)
!                                   nbDP(j1) = nbDP(j1) + 1
!                                   nbDP(j2) = nbDP(j2) + 1
!                                else
                                    fcVal= alphaDD*weightc
                                    uVal =-alphaDD*(rij-0.5*rr)
!                                   nbDP(j1) = nbDP(j1) + 1
!                                   nbDP(j2) = nbDP(j2) + 1
!                                endif
                            endif

                            cigama1 = cigamaD
                            gamma1  = gammaD
                            goto 10
                        endif

! one of j1 and j2 is a drop particle and another is a solvent particle or chain particle
                        if((tag1 .gt. nDpEnd .and. tag2 .le. nDpEnd) .or.   &
                            (tag2 .gt. nDpEnd .and. tag1 .le. nDpEnd)) then

                            fcVal  = alphaFD*weightc
                            uVal   =-alphaFD*(rij - 0.5*rr)
                            cigama1= cigamaFD
                            gamma1 = gammaFD
                            
                            goto 10
                        endif

! both j1 and j2 are the beads of chains(or drops)
			            if(tag1 .gt. nDpEnd    .and. tag2 .gt. nDpEnd .and. &
				           tag1 .le. nChainend .and. tag2 .le. nChainend) then
 		  		            
                            if(int((tag1-nWallAtom)/ChainLen) .eq.  &
                               int((tag2-nWallAtom)/ChainLen)) then

!DP--21/12/2011                 if(abs(j1-j2) .le. 4) goto 10
                                
                                fcVal= alphapp*weightc
                                uVal =-alphapp*(rij-0.5*rr)
                             else
                                fcVal= alphaf *weightc
                                uVal =-alphaf *(rij-0.5*rr)
                             endif
                             goto 10
                        endif

! one of j1 and j2 is a bead of chains(or drops) and another is a solvent particle
			             if((tag1 .gt. nChainend .and. tag2 .le. nChainend) .or.   &
                            (tag2 .gt. nChainend .and. tag1 .le. nChainend)) then

                                        fcVal =  alphafp*weightc
                                        uVal  = -alphafp*(rij - 0.5*rr)
                                        goto 10
                        endif

! both j1 and j2 are solvent particles
    	                 if(tag1 .gt. nChainend .and. tag2 .gt. nChainend) then
                                        fcVal =  alphaf*weightc
                                        uVal  = -alphaf*(rij - 0.5*rr)
                         endif

10                                  fdVal=-gamma1*weightd*rdvij
                                    frVal= cigama1*weightr*rangls()*sdtinv
                                    ftot = fcVal + fdVal + frVal
                                    fstl = ftot*rij
! exclude the random force when computing stress (DP--7/12/11)
!                                   fstl = (fcVal + fdVal) * rij
                                    do k = 1, NDIM
!                                       f = ftot*dr(k)
                                        fCV = fcVal*dr(k)
                                        fDP = fdVal*dr(k)
                                        fRD = frVal*dr(k)

                                        fr = fstl*dr(k)*dr(k)

!                                        raCV(j1,k) = raCV(j1,k) + fCV / mass
!                                        raCV(j2,k) = raCV(j2,k) - fCV / mass

!                                        raDP(j1,k) = raDP(j1,k) + fDP / mass
!                                        raDP(j2,k) = raDP(j2,k) - fDP / mass

!                                        raRD(j1,k) = raRD(j1,k) + fRD / mass
!                                        raRD(j2,k) = raRD(j2,k) - fRD / mass
                                        
!                                        rforce(j1,k) = rforce(j1,k) + fr
!                                        rforce(j2,k) = rforce(j2,k) + fr

                                        raCV(tag1,k) = raCV(tag1,k) + fCV / mass
                                        raCV(tag2,k) = raCV(tag2,k) - fCV / mass

                                        raDP(tag1,k) = raDP(tag1,k) + fDP / mass
                                        raDP(tag2,k) = raDP(tag2,k) - fDP / mass

                                        raRD(tag1,k) = raRD(tag1,k) + fRD / mass
                                        raRD(tag2,k) = raRD(tag2,k) - fRD / mass
                                        
                                        rforce(tag1,k) = rforce(tag1,k) + fr
                                        rforce(tag2,k) = rforce(tag2,k) + fr
                                    enddo
                                    fr = fstl*dr(1)*dr(2)
!                                    rforce(j1,4) = rforce(j1,4) + fr
!                                    rforce(j2,4) = rforce(j2,4) + fr
		 	                        rforce(tag1,4) = rforce(tag1,4) + fr 
		 	                        rforce(tag2,4) = rforce(tag2,4) + fr

                                    fr = fstl*dr(1)*dr(3)
!                                    rforce(j1,5) = rforce(j1,5) + fr
!                                    rforce(j2,5) = rforce(j2,5) + fr
		 	                        rforce(tag1,5) = rforce(tag1,5) + fr 
		 	                        rforce(tag2,5) = rforce(tag2,5) + fr

                                    fr = fstl*dr(2)*dr(3)
!                                    rforce(j1,6) = rforce(j1,6) + fr
!                                    rforce(j2,6) = rforce(j2,6) + fr
		 	                        rforce(tag1,6) = rforce(tag1,6) + fr 
		 	                        rforce(tag2,6) = rforce(tag2,6) + fr

                                    uSum = uSum + uVal
                                    virSum = virSum + fcVal*rij
                                endif
                            endif
                            j2 = cellList(j2)
                        enddo
                        j1 = cellList(j1)
                    enddo
15              enddo
            enddo
        enddo
    enddo

   do n = nStartAtom, nAtom + nAtomCopy
       m = ntag(n)
       do i = 1, m
           j = tag(n,i)
           fcVal = alphawB*(mrho(n) + mrho(j)) * weightcd(n,i)          ! repulsive potential in MDPD
           fstl = fcVal * drij(n,i)
           do k = 1, NDIM
               fCR = fcVal * ddr(n,i,k)
               raCR(n,k) = raCR(n,k) + fCR / mass

               fr = fstl*ddr(n,i,k)*ddr(n,i,k)
               rforce(n,k) = rforce(n,k) + fr
           enddo
           fr = fstl*ddr(n,i,1)*ddr(n,i,2)
           rforce(n,4) = rforce(n,4) + fr

           fr = fstl*ddr(n,i,1)*ddr(n,i,3)
           rforce(n,5) = rforce(n,5) + fr

           fr = fstl*ddr(n,i,2)*ddr(n,i,3)
           rforce(n,6) = rforce(n,6) + fr
       enddo
   enddo

!! --- for Bell Shape potential case -----
!    do n = nWallAtom + 1, nAtom + nAtomCopy
!        m = ntag(n)
!        do i = 1, m
!            j = tag(n,i)
!            if(j .le. nWallAtom) goto 20
!            fcVal = alphaB*(mrho(n) + mrho(j)) * weightcd(n,i)      ! repulsive potential in MDPD
!            fstl  = fcVal * drij(n,i)
!            do k = 1, NDIM
!                fCR = fcVal * ddr(n,i,k)
!                raCR(n,k) = raCR(n,k) + fCR / mass
!
!                fr = fstl*ddr(n,i,k)*ddr(n,i,k)
!                rforce(n,k) = rforce(n,k) + fr
!            enddo
!            fr = fstl*ddr(n,i,1)*ddr(n,i,2)
!            rforce(n,4) = rforce(n,4) + fr
!
!            fr = fstl*ddr(n,i,1)*ddr(n,i,3)
!            rforce(n,5) = rforce(n,5) + fr
!
!            fr = fstl*ddr(n,i,2)*ddr(n,i,3)
!            rforce(n,6) = rforce(n,6) + fr
!20          continue
!        enddo
!    enddo

!    if(runID .eq. 3 .and. stepCount .gt. stepEquil .and.    &
!       mod((stepCount-stepEquil), 20) .eq. 0) then
!       do n = nWallAtom + 1, nDpEnd
!          pct(n) = real(nbDP(n)) / real(nbDP(n) + nbSP(n))
!          if(pct(n) .gt. cpct) DpSign(n) = 1
!       enddo
!    endif

!   print*, 'uSum = ', 0.5*uSUm, '  virSum = ', 0.5*virSum

    deallocate(cellList)
    deallocate(sqrDist)

    return
	end

!----------------------------------------------------------

    subroutine ComputeExternalForce

    implicit none
    include 'dpdflow.h'

!   Parameters
!   Locals

    integer	n, i, k
    real*8 CDp(10,NDIM), dCDp

    if((runId .eq. 1 .or. runId .eq. 4) .and.                        &
       gravField .ne. 0.                .and.                        &
       StepCount .ge. StepEquil) then
        do n = nWallAtom + 1, nAtom
!       do n = nChainEnd + 1, nAtom
	        ra(n,3) = ra(n,3) + gravField
	    enddo
    endif

    if(runId .eq. 2 .and. gravField .ne. 0. .and.  &
       StepCount .ge. StepEquil) then
        do n = nWallAtom + 1, nAtom
!       do n = nChainEnd + 1, nAtom
	        if(r(n,3) .gt. 0.) then
                ra(n,1) = ra(n,1) + gravField
            else
                ra(n,1) = ra(n,1) - gravField
            endif
	    enddo
    endif

    do i = 1, nDp
       do k = 1, NDIM
          do n = nWallAtom + (i-1)*nPDP + 1, nWallAtom + i*nPDP
             if(k .ne. NDIM) CDp(i,k) = CDp(i,k) + (r(n,k) + ncc(n,k)*region(k))
             if(k .eq. NDIM) CDp(i,k) = CDp(i,k) + (r(n,k) + ncc(n,k)*initUcell(k)*gap(k))
          enddo
          CDp(i,k) = CDp(i,k) / nPDP
       enddo
    enddo

!    dCDp = (CDp(1,1) - CDp(2,1)) / RdsDp
!
!    if(dCDp .le. 3.0 .and. dCDp .gt. -5.0) PstAdj = 0
!
    !----drop 1-------
!    do n = nWallAtom + 1, nWallAtom + nPDP
!       ra(n,2) = ra(n,2) - HPstAdj*(CDp(1,2)-CDp0(1,2))
!       if(PstAdj .eq. 1) then
!          ra(n,1) = ra(n,1) - HPstAdj*(CDp(1,1)-CDp0(1,1))
!          ra(n,3) = ra(n,3) - HPstAdj*(CDp(1,3)-CDp0(1,3))
!       endif
!    enddo
       
    !----drop 2-------
!    do n = nWallAtom + nPDP + 1, nDpEnd
!      ra(n,2) = ra(n,2) - HPstAdj*(CDp(2,2)-CDp0(2,2))
!       if(PstAdj .eq. 1) then
!!         ra(n,1) = ra(n,1) - HPstAdj*(CDp(2,1)-CDp0(2,1))
!          ra(n,3) = ra(n,3) - HPstAdj*(CDp(2,3)-CDp0(2,3))
!       endif
!    enddo

!   if(PstAdj .eq. 1) then
!   !----drop 1-------
!      do n = nWallAtom + 1, nWallAtom + nPDP
!         ra(n,1) = ra(n,1) - HPstAdj*(CDp(1,1)-CDp0(1,1))
!         ra(n,2) = ra(n,2) - HPstAdj*(CDp(1,2)-CDp0(1,2))
!         ra(n,3) = ra(n,3) - HPstAdj*(CDp(1,3)-CDp0(1,3))
!      enddo
!   !----drop 2-------
!      do n = nWallAtom + nPDP + 1, nDpEnd
!         ra(n,1) = ra(n,1) - HPstAdj*(CDp(2,1)-CDp0(2,1))
!         ra(n,2) = ra(n,2) - HPstAdj*(CDp(2,2)-CDp0(2,2))
!         ra(n,3) = ra(n,3) - HPstAdj*(CDp(2,3)-CDp0(2,3))
!      enddo
!   endif

	return

	end

!---------------------------------------------------------------------

    subroutine BellShape(rwf, fwf)

    implicit none
    include 'dpdflow.h'

    real*8 rwf, fwf, wca, wcr

    if(rwf .lt. rCut/3.) then
        wca = -18*(rCut/3.)**2/r3Cut + 12*(rCut/3.)/rrCut
    elseif(rwf .lt. 0.5*rCut) then
        wca = -18*rwf**2/r3Cut + 12*rwf/rrCut
    else
        wca =  6*rwf**2/r3Cut - 12*rwf/rrCut + 6/rCut
    endif

    if(rwf .lt. rCut2/3.) then
        wcr = -18*(rCut2/3.)**2/r3Cut2 + 12*(rCut2/3.)/rrCut2
    elseif(rwf .lt. 0.5*rCut2) then
        wcr = -18*rwf**2/r3Cut2 + 12*rwf/rrCut2
    elseif(rwf .le. rCut2) then
        wcr =  6*rwf**2/r3Cut2 - 12*rwf/rrCut2 + 6/rCut2
    endif

    fwf = alphaw*wca + alphawB*wcr

    end