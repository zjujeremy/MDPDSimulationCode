!-----------------------------------------------------------------------
!	dpdpolym.h
!	include file of dpdpolym.f90
!-----------------------------------------------------------------------
!	Parameters to difine the size of the program

	  integer NDIM, MAXATOM, NHIST, MAXHSIZE, MAXWTM, MAXCHAIN
	  parameter (NDIM = 3, NHIST = 6,                                  &
     		     MAXATOM = 600000, MAXCHAIN = 100000,                  &
		         MAXHSIZE = 3000000, MAXWTM = 400000)

	  integer inpt, lcnf, lst1, lst2, lst3, lst4, lst5, forc
  	  integer moreCycles, nAtom, runId, stepAvg, stepCount, nWallAtom, &
	          nFreeAtom, nStartAtom, PNDIM, stepEquil, stepLimit, nPDP,&
			  randSeed, maxList, hSize, startSample, stepSample, nDp,  &
			  nChain, ChainLen, limitChainProps, stepChainProps,       &
			  countChainProps, nChainend, nChainConf, nDpEnd, stepFzn, &
			  stepStop, nAtomCopy
!             nflow, latticeMod, forceMod, stepAdjustTemp,             &
!             stepInitlzTemp
	  integer countGrid, limitGrid, stepGrid, nsolventend, PstAdj
	  integer atomID(MAXATOM), initUcell(NDIM), cells(NDIM), DpFzn(10),&
	          sizeHistGrid(NDIM), ncc(MAXATOM,NDIM), tagLE(MAXATOM),   &
              DpSign(MAXATOM)
	
	  real*8 deltaT, density, kinEnergy, pi, potEnergy, pressure, rCut,&
	         rrCut, sKinEnergy, sPressure, sTotEnergy, ssKinEnergy,    &
     		 ssPressure, ssTotEnergy, temperature, timeNow, gravField, &
     		 totEnergy, uSum, vMag, vSum, virSum, sInitKinEnergy,      &
			 wmingap, binvolm, alphaf, alphafp, alphapp, WLadjust,     &
		     cigamaw, gammaw, lambda, sdtinv, alphaw, alphawf, rcutw,  &
			 wLayer, aaDistSq, eeDistSq, radGyrSq, gMomRatio1, mass,   &
			 gMomRatio2, Hfene, rrfene, rmaxfene, reqfene, timeSteady, &
			 shearRate, RdsDp, alphafd, alphaDD, cigamaF, cigamaD,     &
			 gammaF, gammaD, gammaFD, cigamaFD,                        &
             HPstAdj, LEShiftX, LEShiftV, r3Cut, cpct
!            epslon, nEval, gravity(NDIM) 
	  real*8 r(MAXATOM,NDIM), rv(MAXATOM,NDIM), ra(MAXATOM,NDIM),      &
     		 rw(MAXWTM,NDIM), chainCentre(MAXCHAIN,NDIM),              &
			 sChain(MAXCHAIN), raCV(MAXATOM,NDIM), raDP(MAXATOM,NDIM), &
			 raRD(MAXATOM,NDIM), raSP(MAXATOM,NDIM), rvm(MAXATOM,NDIM)
	  real*8 region(NDIM), regionH(NDIM), histGrid(MAXHSIZE,NHIST),    &
		     strsGrid(MAXHSIZE,7), GridChainLen(MAXHSIZE,2), gap(NDIM),&
			 rforce(MAXATOM,6), wn(MAXWTM,NDIM), CDp0(10,NDIM)
!            profileV(MAXHSIZE), profileT(MAXHSIZE), flowvel(MAXHSIZE,3)
! ---- MDPD	------------------------------------------------------------	
      real*8 rCut2, alphaB, raCR(MAXATOM,NDIM), rrCut2, r3Cut2, alphawB
! ----------------------------------------------------------------------
!------ - Contact Angles----------------------------------------
		real * 8 simuRadsum, averdistancessum, highZsum, CAsum, CAmin, CAmax, CAtemp(5000), CAerrorbar
		integer stepsum
		

!--------------RDF para--------------------------------

		common / RDF / nS, nSB, RDFstepSample, BnGrid
		common / RDFGrid / BmGrid, StrGrid, StateGrid
		common / SPGrid / StateSpGrid, StrSpGrid, deltaQ, cellength
		integer RDFstepSample, BnGrid
		real * 8 deltaQ, cellength
		integer, pointer::nS(:, : ), nSB(:, : )
		real * 8, pointer:: BmGrid(:, : , : , : ), StrGrid(:, : , : , : ), StateGrid(:, : , : , : )
		real * 8, pointer::StateSpGrid(:, : ), StrSpGrid(:, : )
!--------------RDF para--------------------------------

		
      common / int1 / inpt, lcnf, lst1, lst2, lst3, lst4, lst5, forc
	  common / int2 / moreCycles, nAtom, nWallAtom, nFreeAtom, runId,  &
		              stepAvg, stepCount, nStartAtom, PNDIM, stepFzn,  &
     		          stepEquil, stepLimit, randSeed, maxList, hSize,  &
		              startSample, stepSample, nChain, ChainLen, nDp,  &
					  nChainend, nChainConf, nDpEnd, nPDP, stepStop,   &
                      nAtomCopy, tagLE
!                     latticeMod, forceMod, stepAdjustTemp, nflow,     &
!                     stepInitlzTemp 
	  common / int3 / countGrid, limitGrid, stepGrid, limitChainProps, &
		              stepChainProps, countChainProps, nsolventend,    &
					  PstAdj
	  common / int4 / atomID, initUcell, cells, sizeHistGrid, ncc,     &
	                  DpFzn, DpSign

	  common / real1 / deltaT, density, kinEnergy, pi, potEnergy, vMag,&
     		           pressure, rCut, rrCut, sKinEnergy, sPressure,   &
					   sTotEnergy, ssKinEnergy, ssPressure, virSum,    &
					   ssTotEnergy, temperature, timeNow, totEnergy,   &
					   uSum, vSum, gravField, sInitKinEnergy, WLadjust,&
		               wmingap, binvolm, alphaf, alphafp, alphapp,     &
					   cigamaw, gammaw, lambda, sdtinv, alphaw,        &
					   alphawf, rcutw, Hfene, rrfene, wLayer, aaDistSq,&
					   eeDistSq, radGyrSq, gMomRatio1, gMomRatio2,     &
					   rmaxfene, reqfene, timeSteady, RdsDp, alphaFD,  &
					   alphaDD, cigamaF, cigamaD, cigamaFD,            &
                       gammaF, gammaD, gammaFD, mass, HPstAdj,         &
                       r3Cut, cpct
!                      epslon, nEval, gravity 
	  common / real2 / r, rv, ra, rw, region, regionH, histGrid, gap,  &
	                   rforce, strsGrid, chainCentre, sChain,          &
					   GridChainLen, wn, shearRate, raCV, raDP, raRD,  &
					   raSP, CDp0, rvm
!                      profileV, profileT, flowvel
! ----- MDPD -----------------------------------------------------------
      common / realm / rCut2, alphaB, raCR, rrCut2, r3Cut2, alphawB
!-----------------------------------------------------------------------
!------ - Contact Angles----------------------------------------
	common / CAcode / simuRadsum, averdistancessum, highZsum, CAsum, CAmin, CAmax, CAtemp, CAerrorbar
	common / CAint / stepsum
!-------------------------------------------------------------- -

!-------- - squzee droplet par------------------------------------------ - -
	integer  wr_mark(MAXATOM)
	integer startSqueeze, calcu_step, delta_step
	real * 8 squeezeCa, gap_plane, wall_velo
	common / wallmark / wr_mark
	common / wallmarkk / startSqueeze, calcu_step, delta_step
	common / squeezePara / squeezeCa, gap_plane, wall_velo
!--------------------------------------compute contact angles---------

	integer num_steps
	real * 8  Drop_bottom, Drop_top, CL_temp_bottom, CL_temp_top
	real * 8, pointer::each_slice_point(:, : ), each_slice_point_top(:, : )
	integer, pointer::num_particle(:, : ), num_particle_top(:, : )
	common / computecontactangles_real8 / each_slice_point, each_slice_point_top
	common / computecontactangles_int / num_particle, num_particle_top
	common / computecontactangles__int / num_steps
	common / computecontactangles__real8 / Drop_bottom, Drop_top, CL_temp_bottom, CL_temp_top

	real * 8  alphaw_top
	common / drop_topwall_alpha / alphaw_top

	!---------------------------- - HysteresisCA parameter------------
	!hysteresis contact angles Parameters
	integer hysteresis_start, hysteresis_end, g_outputFlu, g_circle_output_delta, g_num_slice_for_GDalgo_YDF, g_num_slice_for_GDalgo_YDF_top, g_num_slice_for_GDalgo, g_num_slice_for_GDalgo_line
	real * 8 g_gap_slice, g_delta_degree, g_first_deltaQ

	integer g_num_steps, g_num_slice, g_num_degree, g_num_cell, g_output_id, g_output_id_YDF, g_positionOver, g_output_CL_id
	real * 8 g_Drop_bottom, g_Drop_top, g_CL_A, g_CL_A_prev, g_CL_A_Velo, g_CL_R, g_CL_R_prev, g_CL_R_Velo, g_Center_aver_X, g_Center_aver_Y
	integer, pointer::g_num_cell_particle(:, : , : ), g_YDF_sum_num_parti(:, : , : )
	real * 8, pointer::g_density_cell(:, : , : ), g_each_degree_point(:, : , : ), g_each_degree_point_YDF(:, : , : ), g_YDF_sum_coor_aver(:, : , : ), g_YDF_sum_coor(:, : , : )
	common / hysteresis_int_read / hysteresis_start, hysteresis_end, g_num_slice_for_GDalgo, g_num_slice_for_GDalgo_line, g_outputFlu, g_circle_output_delta, g_num_slice_for_GDalgo_YDF, g_num_slice_for_GDalgo_YDF_top
	common / hysteresis_real8_read / g_gap_slice, g_delta_degree, g_first_deltaQ
	common / hysteresis_int / g_num_steps, g_num_slice, g_num_degree, g_num_cell, g_output_id, g_output_id_YDF, g_positionOver, g_output_CL_id
	common / hysteresis_real8 / g_Drop_bottom, g_Drop_top, g_CL_A, g_CL_A_Velo, g_CL_R, g_CL_R_Velo, g_CL_A_prev, g_CL_R_prev, g_Center_aver_X, g_Center_aver_Y
	common / hysteresis_int_p / g_num_cell_particle, g_YDF_sum_num_parti
	common / hysteresis_real_p / g_density_cell, g_each_degree_point, g_each_degree_point_YDF, g_YDF_sum_coor_aver, g_YDF_sum_coor



