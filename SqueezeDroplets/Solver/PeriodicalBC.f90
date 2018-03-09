!----------------------------------------------------------
!   Apply the periodical boundary condition on paritcles
!----------------------------------------------------------

	subroutine ApplyBoundaryCond

	implicit none
	include	'dpdflow.h'

!	parameters
!	Locals

	integer	k, n, m, vSign
	real*8 	e(NDIM)

!	skip over wall atoms

 	do n = nWallAtom+1, nAtom
       do k = 1, PNDIM
	      if(r(n,k) .gt. region(k) .or. r(n,k) .lt. -region(k)) then
	         print*, 'n = ', n, ' k = ', k, ' stepCount =', stepCount
	         print*, 'r(n,k) = ', r(n,k), '  region(k) = ', region(k)
	         print*, 'rv(n,k) = ', rv(n,k)
	         print*, 'molcule moves too far, reduce time step'
	         pause
	      endif

	      if(r(n,k) .ge. regionH(k)) then
	            r(n,k) = r(n,k) - region(k)
                ! record the times of the chain particle crossing the boundary (DP--22/11/11)
                ncc(n,k) = ncc(n,k) + 1
                
                if(runId .eq. 5 .and. stepCount .gt. stepEquil .and. k .eq. NDIM) then
!                    r(n,1) = r(n,1) - LEShiftX
                    r(n,1) = r(n,1) + ncc(n,1)*region(1) - LEShiftX
                    ncc(n,1) = anint(r(n,1)/region(1))       !!!!
                    r(n,1) = r(n,1) - anint(r(n,1)/region(1))*region(1)
                    rv(n,1) = rv(n,1) - LEShiftV
                    rvm(n,1) = rvm(n,1) - LEShiftV
                endif
                
	        elseif(r(n,k) .lt. - regionH(k)) then         
	            r(n,k) = r(n,k) + region(k)
                ncc(n,k) = ncc(n,k) - 1

                if(runId .eq. 5 .and. stepCount .gt. stepEquil .and. k .eq. NDIM) then
!                    r(n,1) = r(n,1) + LEShiftX
                    r(n,1) = r(n,1) + ncc(n,1)*region(1) + LEShiftX
                    ncc(n,1) = anint(r(n,1)/region(1))       !!!!
                    r(n,1) = r(n,1) - anint(r(n,1)/region(1))*region(1)
                    rv(n,1) = rv(n,1) + LEShiftV
                    rvm(n,1) = rvm(n,1) + LEShiftV
                endif

	      endif
	   enddo
	enddo

    if(runId .eq. 5) then
        nAtomCopy = 0
        do n = nWallAtom + 1, nAtom
            if(r(n,NDIM) .lt. -regionH(NDIM) + region(NDIM)/(cells(NDIM)-1)) then
                nAtomCopy = nAtomCopy + 1
                m = nAtom + nAtomCopy
                do k = 1, NDIM
                    r(m, k) = r(n,k)
                    rv(m,k) = rv(n,k)
                enddo
                tagLE(m) = n
                r(m, NDIM) = r(m, NDIM) + region(NDIM)
                r(m, 1) = r(m, 1) + LEShiftX
                r(m, 1) = r(m, 1) - anint(r(m, 1)/region(1))*region(1)
                rv(m,1) = rv(m,1) + LEShiftV
!           elseif(r(n,NDIM) .ge. regionH(NDIM) - region(NDIM)/(cells(NDIM)-1)) then
!!              if(r(n,NDIM) .ge. regionH(NDIM) - region(NDIM)/(cells(NDIM)-1)) then
!                   nAtomCopy = nAtomCopy + 1
!                   m = nAtom + nAtomCopy
!                   do k = 1, NDIM
!                       r(m, k) = r(n,k)
!                       rv(m,k) = rv(n,k)
!                   enddo
!                   tagLE(m) = n
!                   r(m, NDIM) = r(m, NDIM) - region(NDIM)
!                   r(m, 1) = r(m, 1) - LEShiftX
!                   r(m, 1) = r(m, 1) - anint(r(m, 1)/region(1))*region(1)
!                   rv(m, 1) = rv(m, 1) - LEShiftV
!               endif
            endif
        enddo
    endif
	return

	end
