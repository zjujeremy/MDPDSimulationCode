!-----------------------------------------------------
!   Extract a pipe from an original cuboid
!-----------------------------------------------------

    subroutine PipeConf

    implicit none
    include 'dpdflow.h'

    integer i, n, k
    integer pcnf
    integer nwall, nfluid
    real*8 prad, prdw, prad2, prdw2, pryz, pryz2
    real*8 trw(nWallAtom, NDIM), twn(nWallAtom, NDIM), trf(nAtom, NDIM), &
           trv(nAtom, NDIM), tra(nAtom, NDIM)

    pcnf = 101

    open (unit = pcnf, file = 'pipeconf.dat', status = 'unknown')

    read (pcnf, *) prad

    prdw = prad + 1.0*gap(1)
    
    prad2 = prad**2
    prdw2 = prdw**2
     
    nwall  = 0
    nfluid = 0

    do n = nStartAtom, nAtom
       pryz2 = r(n,2)**2 + r(n,3)**2
       if(pryz2 .ge. prad2 .and. pryz2 .le. prdw2) then
          nwall = nwall + 1
          pryz = sqrt(pryz2)
          do k = 1, NDIM
             trw(nwall,k) = r(n,k)
          enddo
          twn(nWall, 1) = 0.
          twn(nWall, 2) = -trw(nwall, 2)/pryz
          twn(nwall, 3) = -trw(nwall, 3)/pryz
       elseif(pryz2 .lt. prad2) then
          nfluid = nfluid + 1
          do k = 1, NDIM
             trf(nfluid,k) = r(n,k)
             trv(nfluid,k) = rv(n,k)
             tra(nfluid,k) = ra(n,k)
          enddo
       endif
     enddo

     nWallAtom = 0
     
     do n = 1, nwall
        nWallAtom = nWallAtom + 1
        do k = 1, NDIM
           r(nWallAtom, k)  = trw(n,k)
           wn(nWallAtom, k) = twn(n,k)
           rv(nWallAtom, k) = 0.
           ra(nWallAtom, k) = 0.
        enddo
     enddo

     nAtom = nWallAtom

     do n = 1, nfluid
        nAtom = nAtom + 1
        do k = 1, NDIM
           r(nAtom, k)  = trf(n,k)
           rv(nAtom, k) = trv(n,k)
           ra(nAtom, k) = tra(n,k)
        enddo
     enddo

     nStartAtom = 1
     PNDIM = 1

     nDpEnd = nWallAtom
     nChainend = nWallAtom

     write(30,'(//'' Coordinates of Wall Particles'')')
 	 write(30,'(2x, ''n'', 6x,''x'', 8x, ''y'', 8x, ''z'', 8x, &
		  ''nx'', 7x, ''ny'', 7x, ''nz'')')
	 do n = 1, nWallAtom
	    write(30,'(i6, 6f9.4)') n, (r(n,i), i = 1,3), (wn(n,i),i = 1,3)
	 enddo

 	 write(30,'(//'' Coordinates of Simple Particles''/)')
	 do n = nChainend + 1, nAtom
	    write(30,'(i6, 3f9.4)') n, (r(n,i), i = 1,3)
	 enddo

     end


