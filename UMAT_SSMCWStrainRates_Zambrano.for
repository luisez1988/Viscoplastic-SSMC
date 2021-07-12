!*******************************************************************************
!*******************************************************************************
!**			UMAT prepared by Luis Zambrano-Cruzatty based on Yerro 2015		  **
!**						email: luisez@vt.edu, luedzamb@espol.edu.ec			  **
!*******************************************************************************
!*******************************************************************************
! Copyright 2021 Luis Zambrano-Cruzatty and Alba Yerro

!Permission is hereby granted, free of charge, to any person obtaining a copy of 
!this software and associated documentation files (the "Software"), to deal in the 
!Software without restriction, including without limitation the rights to use, copy, 
!modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, 
!and to permit persons to whom the Software is furnished to do so, subject to the 
!following conditions:

!The above copyright notice and this permission notice shall be included in all copies
!or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
!INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
!PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
!FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
!OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
!DEALINGS IN THE SOFTWARE.	
	

Subroutine UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, &
				DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED,&
				CMNAME, NDI, NSHR, NTENS, &	NSTATEV, PROPS, NPROPS, COORDS, &
				PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, &
				LAYER, KSPT, KSTEP, KINC)

 !DEC$ ATTRIBUTES DLLEXPORT, ALIAS:"UMAT" :: UMAT
  INCLUDE 'ABA_PARAM.INC'

      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV), &
      DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS), &
      STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1), &
      PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

	  ! Arguments:
		!          I/O  Type
		!  PROPS    I   R()  : List with model parameters
		!  DSTRAN   I   R()  : Strain increment
		!  DTIME    I   R()  : Time increment
		!  DDSDDE   O   R(,) : Material stiffness matrix
		!  STRESS  I/O  R()  : stresses
		!  STATEV  I/O  R()  : state variables
		!

!---  Local variables
!
      Dimension :: DE(6,6), Sig(6), dEpsE(6), dEpsP(6), &
				   EpsE(6), EpsP(6), dEps(6), Eps(6), ERate(6), Erate0(6)
	  logical :: switch_smooth, switch_yield


		! Local variables:
		!
		!  DE        : Linear Elastic constitutive matrix
		!  dSig	     : Stress increment vector
		!  Sig	     : Stress vector
		!  dEpsE     : Elastic strain increment vector
		!  dEpsP     : Plastic strain increment vector
		!  dEps      : Total strain increment vector
		!  EpsE      : Elastic strain vector
		!  EpsP      : Plastic strain vector
		!  Eps       : Total strain vector
		!  EpsRate	 : Total strain rate tensor
		!

		  Rad  = 45d0 / datan(1d0)

	    !
		! SSMC w Strain rate model parameters
		!
		! Contents of PROPS(10) MCSS
		!  1 : G       shear modulus
		!  2 : ENU     Poisson's ratio
		!  3 : cp      peak cohesion
		!  4 : cr      residual cohesion
		!  5 : phip    peak friction angle
		!  6 : phir    residual friction angle
		!  7 : psip    peak dilation angle
		!  8 : psir    residual dilation angle
		!  9 : eta		   shape factor
		! 10 : alpha_G	   Shear modulus viscosity factor
		! 11 : alpha_K	   Bulk modulus viscosity factor
		! 12 : alpha_chi   friction and cohesion viscosity
		! 13 : alpha_beta  dilation angle viscosity
		! 14 : RefRate	   Reference strain rate
		! 15 : Switch_smooth Boolean switch for activating strain rate smoothing
		! 16 : N_S		   Degree of smoothening


		G_0      = PROPS(1)		  ! shear modulus
        ENU    = PROPS(2)         ! Poisson's ratio
        cp     = PROPS(3)         ! peak cohesion
        cr     = PROPS(4)         ! residual cohesion
        phip   = PROPS(5)/Rad     ! peak friction angle (rad)
        phir   = PROPS(6)/Rad     ! residual friction angle (rad)
        psip   = PROPS(7)/Rad     ! peak dilation angle (rad)
        psir   = PROPS(8)/Rad     ! residual dilation angle (rad)
        eta = PROPS(9)		      ! shape factor
		alpha_G	   = PROPS(10)	  ! Shear modulus viscosity factor
		alpha_K	   = PROPS(11)	  ! Bulk modulus viscosity factor
		if ((alpha_K==0.0d0).and.(alpha_G>0.0d0)) then ! in case user uses 0 it will be computed from nu
			alpha_K= alpha_G*2*(1+ENU)/(3*(1-2*ENU))
		endif
		alpha_chi  = PROPS(12)	  ! friction and cohesion viscosity
		if ((alpha_chi==0.0d0).and.(alpha_G>0.0d0)) then ! in case user uses 0 it will be computed from nu
			alpha_chi= 0.1*alpha_G
		endif
		alpha_beta = PROPS(13)    ! dilation angle viscosity
		if ((alpha_beta==0.0d0).and.(alpha_G>0.0d0)) then ! in case user uses 0 it will be computed from nu
			alpha_beta= alpha_G
		endif
		RefERate= PROPS(14)		  ! reference strain rate
		if (RefERate==0.0d0) then
			RefERate=2.7e-5
		endif
		call dbltobool(PROPS(15), switch_smooth)  ! switch for activating strain rate smoothening
		N_S=PROPS(16)             ! Degree of smoothening

        G   = STATEV(1)		      ! Previous shear modulus
		bK	 = STATEV(2)		  ! Previous bulk modulus
		c    = STATEV(3)          ! cohesion
        phi  = STATEV(4)          ! friction angle
        psi  = STATEV(5)          ! dilatancy angle
		call dbltobool(STATEV(6),switch_yield)       ! Point is plastifying
		Do i = 1,NTENS
        Erate0(i) = STATEV(6+i)   ! Previous strain rate
		end do
		Do i = 1,NTENS
        EpsP(i) = STATEV(12+i)	  ! Previous plastic strain
		end do
		N_i=STATEV(19)			  ! Current number of strain rate sums
		SUM_rate=STATEV(20)		  ! Current sum of strain rates
		!_____Error control state parameters__________________________________________________
		Error_Euler_max=0.0d0															      !
		Error_Yield_last=0.0d0															      !
		Error_Yield_max=0.0d0																  !
		!							                                                          !
		!_____________________________________________________________________________________!		




		X=1/DTIME
		if (DTIME==0.0d0) then
			ERate= 0.0d0    ! Current strain rate
		else
			ERate= X*DSTRAN ! Current strain rate
		end if



		bK_0= 2*G_0*(1+ENU)/(3*(1-2*ENU))

		IPL=0 !Variable required to check if the material is already in residual cond.

		!***********************************************************************************
		!Call the modified Euler algorithm
		call SSMC_Strain_Rate(NOEL,F1,F2,G_0, bK_0,  G, bK ,cp,cr,phip,phir,psip,psir,eta,c,phi,psi, &
			 stress,EpsP,DSTRAN,dEpsP,Sig, Erate0, ERate, RefERate, DTIME, alpha_G, &
			 alpha_K, alpha_chi, alpha_beta,  switch_smooth, N_S, N_i, SUM_rate, switch_yield ,IPL& !)
				,Error_Euler_max,  Error_Yield_last, Error_Yield_max)
		!************************************************************************************

		!************************************************************************************
		!Stress and state variables updating
		Do i=1,NTENS
          STRESS(i) = Sig(i)
        End Do

        STATEV(1) = G
        STATEV(2) = bK
        STATEV(3) = c
		STATEV(4) = phi
		STATEV(5) = psi
		STATEV(6)= logic2dbl(switch_yield)
		Do i = 1,NTENS
        STATEV(6+i) =Erate(i)   ! Current strain rate
		end do
        Do i = 1,NTENS
        STATEV(12+i) = EpsP(i)
		end do
		STATEV(19)=N_i
		STATEV(20)=SUM_rate
		!_____Error control state parameters__________________________________________________
		!     Comment if not wanted                                                           !
		!_____________________________________________________________________________________!	
		STATEV(21)=Error_Euler_max
		STATEV(22)=Error_Yield_last
		STATEV(23)=Error_Yield_max
		
		!************************************************************************************

		!************************************************************************************
		!Tangent stiffness matrix to be returned done by elastic stiffness
        F1  = bK+(4*G/3)
        F2  = bK-(2*G/3)
        DDSDDE = 0.0
        DDSDDE(1:3,1:3) = F2
        DDSDDE(1,1) = F1
        DDSDDE(2,2) = F1
        DDSDDE(3,3) = F1
        DDSDDE(4,4) = G
        DDSDDE(5,5) = G
        DDSDDE(6,6) = G
		!*************************************************************************************

		!End of UMAT
		Return
end subroutine UMAT





      Subroutine SSMC_Strain_Rate(IntGlo,D1,D2,G_0, K_0, G, K,cp,cr,phip,phir,psip,psir,eta,c,phi,psi, &
			                      Sig0,EpsP,DEps,dEpsP,SigC, Erate0, ERate, RefRate, Dtime,  &
			                      alpha_G, alpha_K, alpha_chi, alpha_beta,  switch_smooth, N_S, N_i, SUM_rate, &
								  switch_yield, IPL,Error_Euler_max,  Error_Yield_last, Error_Yield_max)
      !**********************************************************************
      !
      ! Elasto-plastic constitutive model with strain softening and strain rate effects,
	  ! based on the MOHR-COULOMB criterion (considering modifications of Abbo & Sloan (1995))
      ! and the CT. model proposed by Yerro (2015)
      ! Explicit MODIFIED EULER INTEGRATION SCHEME with automatic error control.
      ! Final correction of the yield surface drift (END OF STEP CORRECTION).
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i,j,n,it,m, MAXITER
      double precision :: F,F0,F2 !Evaluation of the Yield function
      double precision :: alpha, alpha1, alpha0 !Elastic Strain proportion
      double precision :: SSTOL !Tolerance Relative Error
      double precision :: YTOL !Tolerance Relative Error of the yield function evaluation
      double precision :: SPTOL !Tolerance Softening parameters
	  double precision :: LTOL ! Tolerance for unloading path
	  logical :: ApplyStrainRateUpdates !Check if the strain rate path crosses the reference line
	  !double precision :: PTOL ! Tolerance for parallelism between plastic and elastic strain increment
      double precision :: Rn !Relative error
      double precision :: T,DT,T1,beta,DTmin !Sub-stepping parameters
      double precision :: c1,phi1,psi1,c2,phi2,psi2
	  double precision :: cu, phiu, psiu, Gu, Ku
      double precision :: ctol,phitol,psitol !c,phi,psi tolerances
      double precision :: Dcr,Dphir,Dpsir !Difference between current and residual values
      double precision :: moduleEr,moduleSigDSig
      double precision :: EpsPEq,EpsPEq1,EpsPEq2 !Equivalent Plastic Deformation
      double precision :: DEpsPEq !Derivative Strain in function of Equivalent Plastic Deformation
	  double precision :: IErateI, IErate0I, Erate_dev, Erate_dev_0, Erate_vol, Erate_vol_0, dErate_eff ! 2 Norm of the strain rate tensor, time increment
	  double precision :: dGdErate, dKdErate, dcdErate, dphidErate, dpsidErate, kappa
	  double precision, dimension(2) :: DummyArg
      double precision, dimension(6) :: SigYield, SigYield2
      double precision, dimension(6) :: DSigPP,DSigP1,DSigP2
      double precision, dimension(6) :: DEpsPP,DEpsPP1,DEpsPP2
      double precision, dimension(6) :: DEpsS,DEpsSS,DSigE
      double precision, dimension(6) :: EpsP1,EpsP2
      double precision, dimension(6) :: dEpsEqdEpsP,dEpsEqdEpsP1
      double precision, dimension(6) :: sumSg,Er, dErate, dNratedErate, DErateS, ErateYield,DErateSS
	  !double precision, dimension(6) :: DIErateIDEpsRate  !Derivative of the norm of the strain rate respect to the strain rate
      double precision, dimension(3) :: dChisdEpsEq,dChisdEpsEq1 !Variation of softening parameters (c,phi,psi) in function of plastic strain
	  double precision, dimension(3) :: DChidErateef, DChidErateef1, DchiErate ! variation of state parameters with strain rate norm
	  double precision, dimension(3,6):: DChidErate
	  double precision, dimension(6,6):: DE
	  !double precision, dimension(6,6) :: Q ! rotational matrix to accommodate elastic and plastic strains
	  logical :: IsElasticUnloading = .false.
	  logical :: IsFailedStep = .false.
      !In variables
      integer, intent(in) :: IntGlo, N_S !Global ID of Gauss point or particle
      double precision, intent(in) :: G_0, K_0 !Elastic Parameters
	  double precision, intent(in) :: Dtime !Time step
      double precision, intent(in) :: cp,cr,phip,phir,psip,psir,eta !strain Softening parameter
	  double precision, intent(in) :: alpha_K, alpha_G, alpha_chi, alpha_beta, RefRate
      double precision, intent(in), dimension(6) :: Sig0 !Initial Stress
      double precision, intent(in), dimension(6) :: DEps !Incremental total strain
	  double precision, intent(in), dimension(6) :: Erate0   !Strain rate vector
      logical, intent(in):: switch_smooth
	  !Inout variables
      double precision, intent(inout):: D1,D2,c,phi,psi, G, K, SUM_rate !cohesion,friction angle and dilatancy angle
      double precision, intent(inout), dimension(6) :: EpsP !Accumulated Plastic Strain
      double precision, intent(inout), dimension(6) :: SigC !Final Stress
      double precision, intent(inout), dimension(6) :: ERate !Incremental Elastic Stress
	  double precision, intent(inout):: Error_Euler_max, Error_Yield_last, Error_Yield_max
      !Out variables
      integer, intent(out) :: IPL
      double precision, intent(out), dimension(6) :: DEpsP !Incremental plastic strain
	  logical, intent(inout):: 	switch_yield
	  integer, intent(inout):: N_i
	  !************************************************************************************************************
      !Initialization
      DEpsPEq = 0.0d0
      EpsPEq = 0.0d0
      SigYield = 0.0d0
      DEpsP = 0.0d0
      F = 0.0d0
      it = 0

	  if (G==0.0) then
		  G=G_0
	  endif

	  if (K==0.0) then
		  K=K_0
	  endif

	  if (c==0.0) then
		  c=cp
	  endif

	  if (phi==0.0) then
		  phi=phip
	  endif

	  if (psi==0.0) then
		  psi=psip
	  endif

	  !Tolerances
      SSTOL = 0.001d0 !Tolerance Relative Error (10-3 to 10-5)
      YTOL = 0.000000001d0 !Tolerance Error on the Yield surface (10-6 to 10-9)
      SPTOL = 0.0001d0 !Tolerance Softening Parameters (0.0001d0)
	  LTOL=0.01d0 !Tolerance for elastic unloading	  
	  MAXITER=20
      ctol = abs(cp-cr)*SPTOL
      phitol = abs(phip-phir)*SPTOL
      psitol = abs(psip-psir)*SPTOL
      DTmin = 0.000000001d0
	
	  call YieldFunctionValue(IntGlo,Sig0,c,phi,F0) ! initial yield function value
	  !Computes current effective strain rate tensor
	  call TwoNormTensor(Erate,6,IErateI)

	  ! Compute a smoothed strain rate if switch_smooth is true
	  if (switch_smooth) then
		  N_i=N_i+1
		  if (N_i<N_S) then !Not enough values
			  SUM_rate=SUM_rate+IErateI !Accumulate the strain rate
			  IErateI=SUM_rate/N_i !Takes the average
		  else !Enough values
			  call TwoNormTensor(Erate0,6,IErate0I) !previous step strain rate
			  SUM_rate=SUM_rate*(1.0-1.0/N_S) !Approximate sum without one term
			  SUM_rate=SUM_rate+IErateI
			  IErateI=SUM_rate/N_S !Averaged strain rate
		  endif
		  call TwoNormTensor(Erate,6,IErate0I)! Norm of the uncorrected strain rate
		  Erate=(IErateI/IErate0I)*Erate !Corrected strain rate tensor
	  endif

	  !Computes the increment of strain rate
	  dErate= Erate-Erate0
	!*****************************************************************************************
	!Store state variables
	Gu=G
	Ku=K
	phiu=phi
	psiu=psi
	cu=c
		 !Compute the norm of dErate
		 call CalculateEpsPEq(Erate,Erate_dev)
		 call CalculateEpsPEq(Erate0,Erate_dev_0)
		 !For the bulk modulus use the volumetric strain rate
		 Erate_vol=abs(Erate(1)+Erate(2)+Erate(3))
		 Erate_vol_0=abs(Erate0(1)+Erate0(2)+Erate0(3))
		 !call TwoNormTensor(Erate0,6,IErate0I) !previous step strain rate
		 dErate_eff=Erate_dev-Erate_dev_0
		 call check4crossing(Erate_dev_0, Erate_dev, dErate_eff, RefRate,ApplyStrainRateUpdates)
		 !Update G and K		 
		 if (ApplyStrainRateUpdates) then
			call UpdatePardue2StrainRate(alpha_G, Erate_dev_0, Erate_dev, dErate_eff, RefRate, G_0, Gu)			
		 end if
		 call CalculateEpsPEq(EpsP,EpsPEq) !Compute equivalent plastic strain
		 !Update state parameters
		 if (ApplyStrainRateUpdates) then
			 call UpdateSatePardue2StrainRate(EpsPEq, alpha_chi, alpha_beta, Erate_dev_0, Erate_dev,&
												dErate_eff, RefRate,cp, phip, psip, cr, phir, psir,&
											  eta, cu,phiu, psiu)
		 end if
		 dErate_eff=Erate_vol-Erate_vol_0
		 call check4crossing(Erate_vol_0, Erate_vol, dErate_eff, RefRate,ApplyStrainRateUpdates)
		 if (ApplyStrainRateUpdates) then
			call UpdatePardue2StrainRate(alpha_K, Erate_vol_0, Erate_vol, dErate_eff, RefRate, K_0, Ku)
		 endif


	 if (cu < cr.or.phiu < phir.or.psiu < psir) then !Check that the state parameters are never lower than the residual
        cu = max(c,cr)
        phiu = max(phi,phir)
        psiu = max(psi,psir)
	 end if
	!***********************************************************************************
	! Fill elastic material matrix
        D1  = Ku+(4*Gu/3)
        D2  = Ku-(2*Gu/3)
        DE  = 0.0
        DE(1:3,1:3) = D2
        DE(1,1) = D1
        DE(2,2) = D1
        DE(3,3) = D1
        DE(4,4) = Gu
        DE(5,5) = Gu
        DE(6,6) = Gu

	!***********************************************************************************
		! elastic stress increment
        Call MatVec( DE, 6, DEps, 6, DSigE)
        ! elastic stress
        Call AddVec( Sig0, DSigE, 1d0, 1d0, 6, SigC )
	!***********************************************************************************

	!***********************************************************************************
	  !Check for elasticity or plasticity loading
	  !Check the yield function value
      call YieldFunctionValue(IntGlo,SigC,cu,phiu,F)

      !If F<0 then the behavior is elastic --> Return
      if (F <= YTOL) then
		G=Gu
		K=Ku
		c=cu
		phi=phiu
		psi=psiu
        switch_yield=.false.
		IPL=0.0d0
        return
	  end if

      !If F>0, the behavior is elasto-plastic --> Continue
	  !***************************************************************************************

	  Dcr = abs(c - cr)
      Dphir = abs(phi - phir)
      Dpsir = abs(psi - psir)
      !Check if we are in residual conditions or in softening conditions
      if (Dcr <= ctol.and.Dphir <= phitol.and.Dpsir <= psitol) then
        IPL = 1 !IPL=1 Residual conditions --> no changes of the strength parameters
        c = cr
        phi = phir
        psi = psir
      else
        IPL = 2 !IPL=2 Softening Conditions --> changes of the strength parameters
	  end if

	  !****************************************************************************************
      !Determine the proportion (alpha) of the stress increment that lies within the yield function.
	  SigYield = Sig0 
	  DEpsS=DEps
	  ErateYield=Erate0
      if (F0 < -YTOL) then !In this Time increment there is part of elastic behavior
		  call NewtonRaphson(SigYield, phip, phir, cp, cr, psip, psir, G_0, K_0, &
							alpha_chi, alpha_G, alpha_K, alpha_beta, DEpsS, ErateYield, Erate &
							, Erate_dev_0, Erate_dev, Erate_vol, Erate_vol_0,refrate, phi, c, psi, K, G, alpha, &
							eta, EpsPEq, F0, MAXITER, YTOL)		  
	  else
		  !Determine Elastic Unloading path
		  call CheckElasticUnloading(IntGlo, Sig0, DSigE,c,phi,psi, IsElasticUnloading, LTOL)
		  if (IsElasticUnloading) then
			  call NewtonRaphson(SigYield, phip, phir, cp, cr, psip, psir, G_0, K_0, &
								alpha_chi, alpha_G, alpha_K, alpha_beta, DEpsS, ErateYield, Erate &
								, Erate_dev_0, Erate_dev, Erate_vol, Erate_vol_0, refrate, phi, c, psi, K, G, alpha, &
								eta, EpsPEq, F0, MAXITER, YTOL)
		  else
			  alpha = 0.0d0 !In this time increment all is plastic behavior
		  end if
	  end if

	  !****************************************************************************************


      !Determine the elastic portion of the stress increment
      !DSigE = alpha * DSigE !Correct Incremental Elastic Stress

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Determine the plastic portion of the stress increment.
      !The method used is the MODIFIED EULER INTEGRATION SCHEME with error control
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	  dErate=Erate-ErateYield
	  !Sub-stepping parameters
      T = 0.0d0
      DT = 1.0d0
	  m = 0 !Counter
      !Start the yielding
      Do while (T <= 1.0d0)

        Rn = 100
		!DrealTime=DT*Dtime ! Real time increment of sub-step
        call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)

		!****************************************************************************************
		!Calculate first and second estimate of stress increment and compute relative error
            !1)Calculation of the portion of the plastic strain increment (DEpsPP)
            DEpsSS = DT * DEpsS !Portion of the plastic strain increment
			DErateSS= DT * dErate
			!Calculate the new increments of strain rate
			call CalculateEpsPEq(ErateYield, Erate_dev_0)
			Erate_vol_0=abs(ErateYield(1)+ErateYield(2)+ErateYield(3))
			!update the strain rate
			ErateYield=ErateYield+DErateSS
			call CalculateEpsPEq(ErateYield, Erate_dev)
			Erate_vol=abs(ErateYield(1)+ErateYield(2)+ErateYield(3))
			dErate_eff=Erate_dev-Erate_dev_0
			call check4crossing(Erate_dev_0, Erate_dev, dErate_eff, RefRate,ApplyStrainRateUpdates)
            !Calculate a first estimate of the associated stress
            !hardening/softening parameter changes
            call PartialChiToPartialEpsilonEq(Erate_dev_0,alpha_chi, alpha_beta,eta,cp,cr,phip,phir,psip,psir, &
				EpsPEq,dChisdEpsEq,RefRate)
            call PartialEpsEqToPartialEpsP(EpsP,EpsPEq,dEpsEqdEpsP)

			!Actual increment of stress calculation
            call DetermineDSigAndDEpsP(IntGlo,G, K,G_0, K_0, D1, D2, alpha_k, alpha_G, alpha_chi, alpha_beta,&
							Erate_vol_0,Erate_vol, Erate_dev_0, Erate_dev, RefRate, c,phi,psi, cp, phip, psip, EpsPEq, eta,SigYield,&
							dEpsEqdEpsP,dChisdEpsEq, dErate_eff, DEpsSS,DSigP1,DEpsPP1, DchiErate)!First increment of stress
			!First plastic strain
            EpsP1 = EpsP + DEpsPP1
			call CalculateEpsPEq(EpsP1,EpsPEq1) !Determine Equivalent plastic Strain (EpsPEq)
            call PartialChiToPartialEpsilonEq(Erate_dev,alpha_chi, alpha_beta,eta,cp,cr,phip,phir,psip,psir, &
               EpsPEq1,dChisdEpsEq1, RefRate)
            call PartialEpsEqToPartialEpsP(EpsP1,EpsPEq1,dEpsEqdEpsP1)
			!Store initial state variables values
			c1 = c
			phi1 = phi
			psi1 = psi
			!Update the state variables according to the strain rate
			if (ApplyStrainRateUpdates) then
				c1=c1+DchiErate(1)
				phi1=phi1+DchiErate(2)
				psi1=psi1+DchiErate(3)
			end if
			call UpdateStateParameters(dChisdEpsEq1,dEpsEqdEpsP1,DEpsPP1,c1,phi1,psi1)
		!*****************************************************************************************
            !2)Calculate a second estimate of the associated stress
            !hardening/softening parameter changes
            SigYield2 = SigYield + DSigP1

            call DetermineDSigAndDEpsP(IntGlo,G, K,G_0, K_0, D1, D2, alpha_k, alpha_G, alpha_chi, alpha_beta,&
							Erate_vol_0,Erate_vol,Erate_dev_0, Erate_dev, RefRate, c1,phi1,psi1,cp, phip, psip, EpsPEq1, eta,SigYield2,&
							dEpsEqdEpsP1,dChisdEpsEq1, dErate_eff, DEpsSS,DSigP2,DEpsPP2, DchiErate)!First increment of stress
			!Second plastic strain
			EpsP2 = EpsP + DEpsPP2

            call CalculateEpsPEq(EpsP2,EpsPEq2) !Determine Equivalent plastic Strain (EpsPEq)
            call PartialChiToPartialEpsilonEq(Erate_dev,alpha_chi, alpha_beta,eta,cp,cr,phip,phir,psip,psir, &
               EpsPEq2,dChisdEpsEq1, RefRate)
            call PartialEpsEqToPartialEpsP(EpsP2,EpsPEq2,dEpsEqdEpsP1)
			!store initial state variables for second estimation
			c2 = c
			phi2 = phi
			psi2 = psi
			if (ApplyStrainRateUpdates) then
				c2=c2+DchiErate(1)
				phi2=phi2+DchiErate(2)
				psi2=psi2+DchiErate(3)
			end if

		   call UpdateStateParameters(dChisdEpsEq1,dEpsEqdEpsP1,DEpsPP2,c2,phi2,psi2)
			
		!****************************************************************************************
		   
		   !Increment os stress is computed as the averages of the two calculations
            DSigPP = 0.5d0 * (DSigP1 + DSigP2)



			!*****************************************************************************************
            !Calculation of the relative error
            Er = 0.5d0 * (DSigP1 - DSigP2)
            moduleEr = sqrt(Er(1)*Er(1)+Er(2)*Er(2)+ Er(3)*Er(3)+ Er(4)*Er(4)+Er(5)*Er(5)+Er(6)*Er(6))

            sumSg = SigYield + DSigPP
            moduleSigDSig = sqrt(sumSg(1)*sumSg(1) + sumSg(2)*sumSg(2) + sumSg(3)*sumSg(3)+ &
                                 sumSg(4)*sumSg(4) + sumSg(5)*sumSg(5) + sumSg(6)*sumSg(6))

            !Check the relative error (Rn) of the new stresses, with the defined tolerance (SSTOL)
            Rn = (moduleEr/moduleSigDSig)
			!store error
			if (Rn>Error_Euler_max) Error_Euler_max=Rn
			!************************************************************************************

			!4)If Rn>SSTOL, the loop is not finished and the sub-step is recalculated smaller
            if (Rn > SSTOL.and.m < 50) then
                beta = max (0.9d0*(sqrt(SSTOL/Rn)), 0.1d0)
                DT = max (DT*beta, DTmin)
                m = m + 1 !Update counter
				IsFailedStep= .true.
				ErateYield=ErateYield-DErateSS
			else
				 !Update the accumulated stresses, plastic strain and softening parameters
				SigYield = SigYield + DSigPP
				DEpsPP = 0.5d0 * (DEpsPP1 + DEpsPP2)
				DEpsP = DEpsP + DEpsPP
				EpsP = EpsP + DEpsPP
				call CalculateEpsPEq(EpsP,EpsPEq) !Determine Equivalent plastic Strain (EpsPEq)
				call PartialEpsEqToPartialEpsP(EpsP,EpsPEq,dEpsEqdEpsP)
				if (ApplyStrainRateUpdates) then
					call UpdatePardue2StrainRate(alpha_G, Erate_dev_0, Erate_dev, dErate_eff, RefRate, G_0, G)
					call UpdatePardue2StrainRate(alpha_K,Erate_vol_0, Erate_vol, dErate_eff, RefRate, K_0, K)
					c=c+DchiErate(1)
					phi=phi+DchiErate(2)
					psi=psi+DchiErate(3)
				end if
				call UpdateStateParameters(dChisdEpsEq,dEpsEqdEpsP,DEpsPP,c,phi,psi)
				D1  = K+(4*G/3)
				D2  = K-(2*G/3)
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				!!!!!!!!!!!!!!!!!!!!!!! END OF STEP CORRECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

				!**********************************************************************************
				!Check if we are on/in the yield surface, otherwise we are still outside (F>0)
				!and a correction is needed.
				call YieldFunctionValue(IntGlo,SigYield,c,phi,F)
				!store yield function errors
				Error_Yield_last=F
				if (abs(F)>abs(Error_Yield_max)) Error_Yield_max=F
				
				n=0 !Counter
				do while (abs(F) > YTOL.and.n < 10) !The correction is needed
					n = n + 1
					call CalculateEpsPEq(EpsP,EpsPEq)!Determine Equivalent plastic Strain (EpsPEq)
					call PartialChiToPartialEpsilonEq(Erate_dev,alpha_chi, alpha_beta,eta,cp,cr,phip,phir,psip,psir,&
						EpsPEq,dChisdEpsEq,RefRate)
					call PartialEpsEqToPartialEpsP(EpsP,EpsPEq,dEpsEqdEpsP)

					call EndOfStepCorrection(IntGlo,D1,D2,G,IPL,F,SigYield,dChisdEpsEq,dEpsEqdEpsP,EpsP, &
						c,phi,psi, DEpsPP,DEpsSS)
				end do
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!*********Psudo time updating*********************************************************
				beta = min (0.9d0*(sqrt(SSTOL/Rn)), 1.1d0)
				if (IsFailedStep) then !previous step failed
					IsFailedStep = .false.
					beta = min (beta, 1.0d0)
				end if
				!The substep is updated
				T1 = T + DT
				DT = beta * DT
				DT = max (DT, DTmin)
				DT = min (DT, 1.0d0-T1)

				!If T1>1 the calculation is finished
				If (T1 >= 1d0) then
					SigC = SigYield   !Determine Final stresses
					switch_yield=.true. ! Because it is yielding, no check on state variables is needed
					return
				end if
				T = T1
		!**********************************************************************************
		end if
      end do  !If T=1 the loop is finished
	  SigC = SigYield
	  switch_yield=.true. !Because it is yielding, no check on state variables is needed
        !end do
	end subroutine SSMC_Strain_Rate





	Subroutine CheckElasticUnloading(IntGlo, Sig, DSig,c,phi,psi, IsElasticUnloading, LTOL)
	!**********************************************************************
	!
	! It returs true if the stress path correspond to elastic unloading
	!
	!**********************************************************************

	implicit none

	!Local Variables
	double precision ::  p,J,Lode,S3TA, NormDSig, NormDFDSig, DSigdpDFDSig, CTeta
	double precision, dimension(6) :: DFDsig, dPdSig
	!In Variables
	double precision, intent(in), dimension(6) :: Sig, DSig
	double precision, intent(in) :: c, phi, psi, LTOL
	integer, intent(in) :: IntGlo       !Global ID of Gauss point or particle
	!Out Variables
	Logical, intent(out) :: IsElasticUnloading

	  !Calculation of the invariants (p',J,Lode)
      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
	  !Calulate yield function derivative with respect to stress
	  call PartialForGtoPartialSigma(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,dPdSig)
	  !calculate the 2-norms of the tensors
	  call TwoNormTensor(Dsig, 6, NormDSig)
	  call TwoNormTensor(DFDSig, 6, NormDFDSig)
	  !Calculate Dsig:DFDSig
	  call MatAdpMatB(DSig, DFDSig, 6, DSigdpDFDSig)

	  !Cos of angle between stress path and normal vector to yield surface
	  CTeta=DSigdpDFDSig/(NormDSig*NormDFDSig)

	  If (CTeta < -LTOL) then !Elastic Unloading
		  IsElasticUnloading=.true.
	  else
		  IsElasticUnloading=.false.
	  end if
	end subroutine CheckElasticUnloading
	
	
	Subroutine NewtonRaphson(Sig0, phip, phir, cp, cr, psip, psir, G_0, K_0, &
							alpha_chi, alpha_G, alpha_K, alpha_psi, dEps, Erate_0, Erate &
							, IErate_0I, IErateI, Erate_vol, Erate_vol_0, refrate, phi, c, psi, K, G, alpha, &
							eta, Epspq, F0, MAXITER, FTOL)
	!_____________________________________________________________________________________
	!  This subroutine determine the elastic proportion with consistent state and elastic
	!  parameters
	!_____________________________________________________________________________________
	implicit none
	
	!Input variables
	double precision, dimension(6), intent(in):: Erate
	double precision, intent(in):: phip, phir, cp, cr, psip, psir, G_0, K_0, &
								  alpha_chi, alpha_G, alpha_K, alpha_psi, refrate, &
								  IErateI, IErate_0I,Erate_vol, Erate_vol_0, FTOL, Epspq, eta
	integer, intent(in):: MAXITER
	
	!output variables
	double precision, dimension(6), intent(inout):: Sig0, dEps, Erate_0
	double precision, intent(inout):: phi, c, psi, K, G, alpha, F0
	
	!local variables
	double precision:: FT, dG, dK, dc, dphi, dpsi, K_u, G_u, phi_u, c_u, psi_u, &
					   derate_alpha(6), dEps_alpha(6), Erate_alpha(6), IErateI_alpha, dErate_eff, &
						dSig_vel(6), DE(6,6), D1, D2, SIg_alpha(6), p, J, Lode, S3TA, DFDSig(6), &
						dPPdSig(6), dFdChis(2), F_prime, aux(6,6), dDdK(6,6), dDdG(6,6), dSigdAlpha(6), &
						Erate_vol_alpha
	
	integer:: n, I
	logical:: ApplyStrainRateUpdates
	
	FT=1000
	alpha=1.0d0
	n=0
	
	!Compute change due to strain rate
	dErate_eff=IErateI-IErate_0I
	call check4crossing(IErate_0I, IErateI, dErate_eff, RefRate,ApplyStrainRateUpdates)
	if (ApplyStrainRateUpdates) then !Check if the strain rates are enough to cause updating
		dG=G_0*alpha_G*log10(IErateI/IErate_0I)
		dphi=phip*alpha_chi*log10(IErateI/IErate_0I)*exp(-eta*Epspq)
		dc=cp*alpha_chi*log10(IErateI/IErate_0I)*exp(-eta*Epspq)
		dpsi=psip*alpha_psi*log10(IErateI/IErate_0I)*exp(-eta*Epspq)
	else 
		dG=0.0d0
		dphi=0.0d0
		dc=0.0d0
		dpsi=0.0d0	
	endif
	
	dErate_eff=Erate_vol-Erate_vol_0
	call check4crossing(Erate_vol_0, Erate_vol, dErate_eff, RefRate,ApplyStrainRateUpdates)
	if (ApplyStrainRateUpdates) then
		dK=K_0*alpha_K*log10(Erate_vol/Erate_vol_0)
	else
		dK=0.0d0
	endif
	
	do while ((abs(FT)>=FTOL).and.(n<=MAXITER))
		!____ store variables_______________________________________________________________________
		K_u=K
		G_u=G
		phi_u=phi
		c_u=c
		psi_u=psi
		F0=FT
		n=n+1
		
		!___Compute trial strain and strain rate____________________________________________________
		derate_alpha=alpha*(Erate-Erate_0)
		dEps_alpha=alpha*dEps
		Erate_alpha=Erate_0+derate_alpha
		call CalculateEpsPEq(Erate_alpha,IErateI_alpha)
		dErate_eff=IErateI_alpha-IErate_0I
		
		!__ Update state and elastic parameters_____________________________________________________
		call check4crossing(IErate_0I, IErateI_alpha, dErate_eff, RefRate,ApplyStrainRateUpdates)
		
		if (ApplyStrainRateUpdates) then !Update state parameters due to strain rate
			call UpdatePardue2StrainRate(alpha_G, IErate_0I, IErateI_alpha, dErate_eff, RefRate, G_0, G_u)			
			call UpdateSatePardue2StrainRate(Epspq, alpha_chi, alpha_psi, IErate_0I, IErateI_alpha, dErate_eff, RefRate,&
												cp, phip, psip,cr, phir, psir, eta, c_u,phi_u, psi_u)
			
		endif
		!compute volumetric strain rates
		Erate_vol_alpha=abs(Erate_alpha(1)+Erate_alpha(2)+Erate_alpha(3))
		dErate_eff=Erate_vol_alpha-Erate_vol_0
		call check4crossing(Erate_vol_0, Erate_vol_alpha, dErate_eff, RefRate,ApplyStrainRateUpdates)
		if (ApplyStrainRateUpdates) then 
			call UpdatePardue2StrainRate(alpha_K,Erate_vol_0, Erate_vol_alpha, dErate_eff, RefRate, K_0, K_u)
		endif
		
		!__ Update trial stress________________________________________________________________________
		! Fill elastic material matrix
        D1  = K_u+(4*G_u/3)
        D2  = K_u-(2*G_u/3)
        DE  = 0.0
        DE(1:3,1:3) = D2
        DE(1,1) = D1
        DE(2,2) = D1
        DE(3,3) = D1
        DE(4,4) = G_u 
        DE(5,5) = G_u
        DE(6,6) = G_u
		! elastic stress increment
        Call MatVec( DE, 6, dEps_alpha, 6, dSig_vel)
        ! elastic stress
        Call AddVec( Sig0, dSig_vel, 1d0, 1d0, 6, Sig_alpha )
		
		!__ Evaluate Yield function and derivatives____________________________________________________
		call YieldFunctionValue(1,Sig_alpha,c_u,phi_u,FT)
		call CalculateInvariants(1,Sig_alpha,p,J,Lode,S3TA)
		call PartialForGtoPartialSigma(Sig_alpha,p,J,Lode,S3TA,c_u,phi_u,psi_u,DFDSig,dPPdSig)
		call PartialFtoPartialChi(p,J,Lode,S3TA,c_u,phi_u,dFdChis)
		
		!__ Update alpha________________________________________________________________________________
		F_prime=dFdChis(1)*dc+dFdChis(2)*dphi
		dDdK=0.0
		dDdK(1:3,1:3) = 1.0
		dDdG=0.0
		dDdG(1:3,1:3)=-2./3.
		dDdG(1,1)=4./3.
		dDdG(2,2)=4./3.
		dDdG(3,3)=4./3.
		dDdG(4,4)=1.0
		dDdG(5,5)=1.0
		dDdG(5,5)=1.0
		aux=alpha* (dK*dDdK+dG*dDdG)
		aux=DE+aux
		call MatVec(aux, 6, dEps, 6, dSigdAlpha)
		do I=1,6 !dot product
			F_prime=F_prime+dFdSig(I)*dSigdAlpha(I)
		enddo
		alpha=alpha-FT/F_prime
	enddo
	K=K_u
	G=G_u
	phi=phi_u
	c=c_u
	psi=psi_u
	Sig0=Sig_alpha
	dEps=(1-alpha)*dEps
	Erate_0=Erate_0+(1-alpha)*(Erate-Erate_0)	
	end subroutine NewtonRaphson
	


 Subroutine YieldFunctionValue(IntGlo,Sig,c,phi,F)
      !**********************************************************************
      !
      ! In this subroutine the yield function evaluated is a smooth hyperbolic approximation to the
      ! Mohr-Coulomb yield criterion (Abbo and Sloan, 1995).
      !
      ! The edges of the hexagonal pyramid and the tip have been smoothed.
      ! There are two parameters aSmooth (smoothes the tip) and ATTRAN(smoothes the edges)
      ! In this case aSmooth=0.0005*c*cot(phi) and LodeT=25�.
      ! If aSmooth=0 and LodeT=30� the original Mohr-Coulomb is obtained.
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision ::  p,J,Lode,S3TA !Invariants
      double precision ::  COH, SPHI, CPHI, COTPHI, STA, CTA, K, aSmooth, ASPHI2, SGN, A, B
      double precision, parameter :: C00001 = 1.0d0 !Parameters
      double precision, parameter :: C00003 = 3.0d0
      double precision, parameter :: C00P50 = 0.0005d0
      double precision, parameter :: C00000 = 0.0d0
      double precision, parameter :: C00IR3 = 0.577350269189626d0
      double precision, parameter :: C000P1 = 0.00000000001D0
      !Constants for rounded K function (for LodeT=25)
      !double precision, parameter :: A1 = 1.432052062044227d0
      !double precision, parameter :: A2 = 0.406941858374615d0
      !double precision, parameter :: B1 = 0.544290524902313d0
      !double precision, parameter :: B2 = 0.673903324498392d0
      !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !double precision, parameter :: A1 = -138300705.446275
      !double precision, parameter :: A2 = -138300706.472675
      !double precision, parameter :: B1 = -138300706.3123
      !double precision, parameter :: B2 = 0.192450089729875
      !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
      !In variables
      double precision, intent(in), dimension(6) :: Sig
      double precision, intent(in) :: c,phi
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle

      !Out variables
      double precision, intent(out) :: F

      F = C00000

      !Calculation of the invariants (p',J,Lode)
      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!Evaluation of the yield function with Smoothing!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !Material parameters
      COH = c     !Cohesion
      SPHI = sin(phi)
      CPHI = cos(phi)
      COTPHI = CPHI/SPHI
      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
      if (abs(phi) == C00000) then
            ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
      end if

      !Calculate K function
      if (abs(Lode) < ATTRAN) then
        STA = sin(Lode)
        CTA = cos(Lode)
        K = CTA - STA*SPHI*C00IR3
      else
        SGN = SIGN(C00001,Lode)
        A = A1 + A2*SGN*SPHI
        B = B1*SGN + B2*SPHI
        K = A - B*S3TA
      end if

      !Calculate value of Hyperbolic Yield function
      F = p*SPHI + sqrt(J*J*K*K+ASPHI2) - COH*CPHI

	end subroutine YieldFunctionValue


	Subroutine PartialForGtoPartialSigma(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,dPPdSig)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the yield function (F) and the plastic potencial punction (P).
      ! Based on Abbo & Sloan (1995)
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i
      double precision :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B, &
                           D, aSmooth, ASPHI2, SGN, T3TA, C3TA, J2, psi2
      double precision ::   K, dKdLode
      double precision :: SPSI, CPSI, TPSI, COTPSI, ASPSI2
      double precision :: i1, i2, Sx, Sy, Sz
      double precision :: DFDp,DFDJ,DFDLode !Derivatives F respect Invariants
      double precision :: DPDp,DPDJ,DPDLode !Derivatives P respect Invariants
      double precision :: C1, C2, C3
      double precision, dimension(6):: DpDSig,DJDSig,DJ3DSig !Derivatives Invariants

      double precision, parameter :: C00001 = 1.0D0 !Parameters
      double precision, parameter :: C000P5 = 0.5D0
      double precision, parameter :: C00P50 = 0.0005D0
      double precision, parameter :: C00000 = 0.0D0
      double precision, parameter :: C00003 = 3.0D0
      double precision, parameter :: C00004 = 4.0D0
      double precision, parameter :: C00002 = 2.0D0
      double precision, parameter :: CP3333 = 0.333333333333333D0
      double precision, parameter :: C00IR3 = 0.577350269189626D0
      double precision, parameter :: C0R3I2 = 0.866025403784439D0
      double precision, parameter :: C000P1 = 0.000000000000001D0
      double precision, parameter :: J0 = 0.001D0
      !Constants for rounded K function (for LodeT=25)
      !double precision, parameter :: A1 = 1.432052062044227d0
      !double precision, parameter :: A2 = 0.406941858374615d0
      !double precision, parameter :: B1 = 0.544290524902313d0
      !double precision, parameter :: B2 = 0.673903324498392d0
      !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !double precision, parameter :: A1 = -138300705.446275
      !double precision, parameter :: A2 = -138300706.472675
      !double precision, parameter :: B1 = -138300706.3123
      !double precision, parameter :: B2 = 0.192450089729875
      !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians
      !In variables
      double precision, intent(in) ::  c,phi,psi !Soft Parameters
      double precision, intent(in), dimension(6) :: Sig
      !Out variables
      double precision, intent(out), dimension(6) :: DFDSig, dPPdSig !Derivatives respect Sigma
      !Inout variables
      double precision, intent(inout) :: p,J,Lode,S3TA !Invariants

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! DFDSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !Material parameters
      COH = c !Cohesion
      SPHI = sin(phi)
      CPHI = cos(phi)
      COTPHI = CPHI/SPHI
      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
      if (abs(phi) == C00000) then
        ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
      end if

      if (J == C00000) then
        J2 = C000P1
        J = sqrt(J2)
      else
        J2 = J*J
      end if

      CTA = cos(Lode)
      C3TA = CTA*(C00004*CTA*CTA-C00003)
      T3TA = S3TA/C3TA

      !Calculate K function and its derivative
      if (abs(Lode) < ATTRAN) then
        STA = S3TA/(C00004*CTA*CTA-C00001)
        K = CTA - STA*SPHI*C00IR3
        dKdLode =  - STA - C00IR3*SPHI*CTA
      else
        SGN = SIGN(C00001,Lode) ! It puts the Lode's sign to the number 1
        A = A1 + A2*SGN*SPHI
        B = B1*SGN + B2*SPHI
        K = A - B*S3TA
        dKdLode = - C00003*B*C3TA
      end if

      !Derivative Dp/DSig
      DpDSig(1) = CP3333
      DpDSig(2) = CP3333
      DpDSig(3) = CP3333
      DpDSig(4) = C00000
      DpDSig(5) = C00000
      DpDSig(6) = C00000

      !Derivative DJ/DSig
      i1 = C000P5/J
      if (J < 0.0001) then
        i1 = 0.0d0
      end if
      Sx = Sig(1)-p
      Sy = Sig(2)-p
      Sz = Sig(3)-p

      DJDSig(1) = i1 * Sx
      DJDSig(2) = i1 * Sy
      DJDSig(3) = i1 * Sz
      DJDSig(4) = i1 * C00002 * Sig(4)
      DJDSig(5) = i1 * C00002 * Sig(5)
      DJDSig(6) = i1 * C00002 * Sig(6)

      !Derivative DJ3/DSig
      i2 = CP3333*J*J
      DJ3DSig(1) = (Sy*Sz - Sig(5)*Sig(5) + i2)
      DJ3DSig(2) = (Sx*Sz - Sig(6)*Sig(6) + i2)
      DJ3DSig(3) = (Sx*Sy - Sig(4)*Sig(4) + i2)
      DJ3DSig(4) = C00002*(Sig(5)*Sig(6) - Sz*Sig(4))
      DJ3DSig(5) = C00002*(Sig(6)*Sig(4) - Sx*Sig(5))
      DJ3DSig(6) = C00002*(Sig(4)*Sig(5) - Sy*Sig(6))

      D = J*K/(sqrt(J2*K*K + ASPHI2))

      !C1F
      C1 = SPHI
      !C2F
      C2 = D*K - T3TA*D*dKdLode
      !C3F
      C3 = -C0R3I2*dKdLode*D/(J2*C3TA)

      !DFDSig!
      do i=1,6
        DFDSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!! dPdSig = DFDSig (if associated Flow Rule)  !!!!!!!!!!!!!!!!!!!!!!
      !!!!! or
      !!!!! dPdSig = C1*DPDSig + C2*DJDSig + C3*DJ3DSig  !!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (abs(J) < J0) then
        psi2 = phi - abs(J)*(phi - psi)/J0
      else
        psi2 = psi
      end if

      if (phi == psi2) then !If Associated Flow Rule, then dPdSig = DFDSig
        dPPdSig = DFDSig

      else !If Non-Associated Flow Rule, then calculate...
        !Material parameters
        SPSI = sin(psi2)
        CPSI = cos(psi2)
        COTPSI = CPSI/SPSI
        aSmooth = C00P50*COH*COTPSI !Smoothing parameter
        ASPSI2 = aSmooth*aSmooth*SPSI*SPSI
        if (abs(psi2) == C00000)then
            ASPSI2 = C00000
        end if

        !Calculate K function and its derivative
        if (abs(Lode) <= ATTRAN) then
            K = CTA - STA*SPSI*C00IR3
            dKdLode = - STA - C00IR3*SPSI*CTA
        else
            A = A1 + A2*SGN*SPSI
            B = B1*SGN + B2*SPSI
            K = A - B*S3TA
            dKdLode = - C00003*B*C3TA
        end if

        D = J*K/(sqrt(J*J*K*K + ASPSI2))

        !C1F
        C1 = SPSI
        !C2F
        C2 = D*K - T3TA*D*dKdLode
        !C3F
        C3 = -C0R3I2*dKdLode*D/(J2*C3TA)

        !dPdSig
        do i=1,6
            dPPdSig(i) = C1*DpDSig(i) + C2*DJDSig(i) + C3*DJ3DSig(i)
        end do

      end if

	end subroutine PartialForGtoPartialSigma

	Subroutine UpdateStateParameters(dChidEpseq,dEpsEqdEpsp,dEpsP,c, phi, psi )
      !**********************************************************************
      !
      ! Calculation of strength parameters (c, phi, psi)
      !
      !**********************************************************************

      implicit none
	  !input variables
	  double precision, dimension(3), intent(in):: dChidEpseq
	  double precision, dimension(6), intent(in):: dEpsEqdEpsp, dEpsP
	  !output variables
	  double precision, intent(inout):: c, phi, psi
	  !Local Variables
	  double precision, dimension(3,6):: dChidEpsp
	  double precision, dimension(3):: dChi
	  integer:: i, j



	  !Compute dChidEpsp=dChidEpseq*dEpsEqdEpsp
	  do i=1, 3
		  dChi(i)=0.0d0
		  do j=1, 6
			  dChidEpsp(i,j)=dChidEpseq(i)*dEpsEqdEpsp(j)
		  end do
	  end do

	  !Compute increment due to plastic strain
	  do i=1, 3
		  do j=1,6
			  dChi(i)=dChi(i)+dChidEpsp(i,j)*dEpsP(j)
		  end do
	  end do

	  !update state parameters
	  c=c+dChi(1)
	  phi=phi+dChi(2)
	  psi=psi+dChi(3)
      end subroutine UpdateStateParameters

	subroutine UpdateSatePardue2StrainRate(EpsPEq, alpha_chi, alpha_beta, IErate0I, IErateI, dErate_eff, RefRate, &
	                                         cp, phip, psip, cr, phir, psir, eta, c, phi, psi)

	!***************************************************************************
	! Updates the state parameters
	!
	!***************************************************************************
	!Input parameters
	double precision, intent(in):: EpsPEq, alpha_chi, alpha_beta, IErateI, dErate_eff, IErate0I, RefRate
	double precision, intent(in):: cp, phip, psip, eta, cr, phir, psir
	!Output variables
	double precision, intent(inout):: c,phi, psi
	!local variables
	integer :: i, n
	double precision:: cp_rate, phi_rate, psi_rate

	!Update c
	cp_rate=cp*(1+alpha_chi*log10(IErateI/RefRate))
	c=cr+(cp_rate-cr)*exp(-eta*EpsPEq)

	!Update phi
	phi_rate=phip*(1+alpha_chi*log10(IErateI/RefRate))
	phi=phir+(phi_rate-phir)*exp(-eta*EpsPEq)
	!Update psi
	psi_rate=psip*(1+alpha_beta*log10(IErateI/RefRate))
	psi=psir+(psi_rate-psir)*exp(-eta*EpsPEq)
	end subroutine UpdateSatePardue2StrainRate
	
	
		subroutine UpdatePardue2StrainRate(alpha,IErate0I, IErateI, dNErate, Refrate, Par_0, Par)
	!***************************************************************************
	! Updates parameters using the logarithmic law
	!
	!***************************************************************************
	implicit none
	!input variables
	double precision, intent(in):: alpha, Refrate, dNErate, Par_0
	!inout variables
	double precision, intent(inout):: Par,  IErateI, IErate0I
	!local variables
	integer :: i, n, usecase
	double precision:: eratio, dPar, K
	
	Par=Par_0*(1.0d0+alpha*log10(IErateI/refrate))	

	end subroutine UpdatePardue2StrainRate


	subroutine PartialChiToPartialEpsilonEq(NEpsRate, K, Kpsi, factor,cp,cr,phip,phir,psip,psir, &
		EpsPEq,dChisdEpsEq, RefRate)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the state parameters with respect
      ! the equivalent plastic shear strain
      !
      !**********************************************************************

      implicit none

      !Local Variables
	  double precision :: Tpsip, Tphip, phipc, cpc, psipc, aux
	  !In Variables
      double precision, intent(in) :: NEpsRate, EpsPEq, K, Kpsi, RefRate
      double precision, intent(in) :: factor,cp,cr,phip,phir,psip,psir
      !Out Variables
      double precision, intent(out), dimension(3):: dChisdEpsEq

	  if (NEpsRate>RefRate) then
	  !Increase in peak friction angle due to strain rates:
	  Tphip=phip*(1+K*log10(NEpsRate/RefRate))
	  Tpsip=psip*(1+Kpsi*log10(NEpsRate/RefRate))
	  !Increase in Cohesion due to strain rates
	  cpc=cp*(1+K*log10(NEpsRate/RefRate))
	  else
		Tphip=phip
		Tpsip=psip
		cpc=cp
	  end if
      !Derivative Cohesion respect Equivalent Plastic Strain (Dc/DPEq)
      dChisdEpsEq(1) = -factor * (cpc - cr) * (exp(-factor*EpsPEq))
      !Derivative Friction angle respect Equivalent Plastic Strain (Dphi/DPEq)
      dChisdEpsEq(2) = -factor * (Tphip - phir) * (exp(-factor*EpsPEq))
      !Derivative Dilatancy angle respect Equivalent Plastic Strain (Dpsi/DPEq)
      dChisdEpsEq(3) = -factor * (Tpsip - psir) * (exp(-factor*EpsPEq))

	end subroutine PartialChiToPartialEpsilonEq


	subroutine PartialEpsEqToPartialEpsP(EpsP,EpsPEq,dEpsEqdEpsP)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the equivalent plastic shear strain
      ! with respect the plastic strain
      !
      !**********************************************************************

      implicit none

      !Local Variables
      double precision :: k1, k2, k3
      double precision :: EpsPM
      double precision, dimension(3) :: EpsDev
      !In Variables
      double precision, intent(in), dimension(6) :: EpsP
      double precision, intent(in) :: EpsPEq
      !Out Variables
      double precision, intent(out), dimension(6):: dEpsEqdEpsP

      k1 = 2.0d0/(3.0d0*EpsPEq)
      if (EpsPEq < 0.00000000001d0) then
        k1 = 0.0d0
      end if
      k2 = k1 * 1.0d0/3.0d0
      k3 = k1 * 2.0d0

      EpsPM = k2 * (EpsP(1) + EpsP(2) + EpsP(3))
      EpsDev(1) = EpsP(1)-EpsPM
      EpsDev(2) = EpsP(2)-EpsPM
      EpsDev(3) = EpsP(3)-EpsPM

      dEpsEqdEpsP(1) = k2 * ( 2.0d0*EpsDev(1) - EpsDev(2) - EpsDev(3))
      dEpsEqdEpsP(2) = k2 * (-EpsDev(1) + 2.0d0*EpsDev(2) - EpsDev(3))
      dEpsEqdEpsP(3) = k2 * (-EpsDev(1) - EpsDev(2) + 2.0d0*EpsDev(3))
      dEpsEqdEpsP(4) = k2 * EpsP(4)
      dEpsEqdEpsP(5) = k2 * EpsP(5)
      dEpsEqdEpsP(6) = k2 * EpsP(6)

	end subroutine PartialEpsEqToPartialEpsP


  Subroutine PartialFtoPartialChi(p,J,Lode,S3TA,c,phi,dFdChis)
      !**********************************************************************
      !
      ! Calculation of the derivatives of the yield function (F) with respect the state parameters
      ! The state parameters are: cohesion (COH) and friction angle (PHI)
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision :: COH, SPHI, CPHI, TPHI, COTPHI, STA, CTA, A, B, &
                          Denom, Num, aSmooth, ASPHI2, SGN
      double precision :: K, dKdPhi, dadc, dadPhi
      double precision, parameter :: C00001 = 1.0D0 !Parameters
      double precision, parameter :: C00P50 = 0.0005D0
      double precision, parameter :: C00000 = 0.0D0
      double precision, parameter :: C00003 = 3.0D0
      double precision, parameter :: C00002 = 2.0D0
      double precision, parameter :: C00IR3 = 0.577350269189626D0
      double precision, parameter :: C000P1 = 0.00000000001D0
      !Constants for rounded K function (for LodeT=25)
      !double precision, parameter :: A1 = 1.432052062044227d0
      !double precision, parameter :: A2 = 0.406941858374615d0
      !double precision, parameter :: B1 = 0.544290524902313d0
      !double precision, parameter :: B2 = 0.673903324498392d0
      !double precision, parameter :: ATTRAN = 0.436332312998582d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=29.5)
      double precision, parameter :: A1 = 7.138654723242414d0
      double precision, parameter :: A2 = 6.112267270920612d0
      double precision, parameter :: B1 = 6.270447753139589d0
      double precision, parameter :: B2 = 6.398760841429403d0
      double precision, parameter :: ATTRAN = 0.514872129338327d0 !Smoothing parameter: LodeT in radians
      !Constants for rounded K function (for LodeT=30)
      !double precision, parameter :: A1 = -138300705.446275
      !double precision, parameter :: A2 = -138300706.472675
      !double precision, parameter :: B1 = -138300706.3123
      !double precision, parameter :: B2 = 0.192450089729875
      !double precision, parameter :: ATTRAN = 0.523598776 !Smoothing parameter: LodeT in radians

      !In variables
      double precision, intent(in) :: p,J,Lode,S3TA !Invariants
      double precision, intent(in) :: c,phi !Soft Parameters
      !Out variables
      double precision, intent(out), dimension(2) :: dFdChis !Derivatives respect Soft Parameters


      !Material parameters
      COH = c !Cohesion
      SPHI = sin(phi)
      CPHI = cos(phi)
      COTPHI = CPHI/SPHI

      !Calculate aSmooth and its derivatives
      if (abs(phi) == C00000) then
        COTPHI = C00000
        dadc = C00000
        dadPhi = C00000
      else
        dadc = C00P50*CPHI/SPHI
        dadPhi = - C00P50*COH/(SPHI*SPHI)
      end if
      aSmooth = C00P50*COH*COTPHI !Smoothing parameter
      ASPHI2 = aSmooth*aSmooth*SPHI*SPHI
      if (abs(phi) == C00000) then
       ASPHI2 = C00P50*C00P50*COH*COH*CPHI*CPHI
      end if

      !Calculate K function and its derivatives
      if (abs(Lode) <= ATTRAN) then
        STA = sin(Lode)
        CTA = cos(Lode)
        K = CTA - STA*SPHI*C00IR3
        dKdPhi = - C00IR3*CPHI*STA
      else
        SGN = SIGN(C00001,Lode) !It puts the Lode's sign to the number 1
        A = A1 + A2*SGN*SPHI
        B = B1*SGN + B2*SPHI
        K = A - B*S3TA
        dKdPhi = A2*SGN*CPHI - B2*CPHI*S3TA
      end if

      !Operating..
      Denom = (sqrt(J*J*K*K + ASPHI2))
      Num =  J*J*K*dKdPhi + aSmooth*SPHI*SPHI*dadPhi + aSmooth*aSmooth*SPHI*CPHI

      !Derivative DF/Dc
      dFdChis(1) = aSmooth*SPHI*SPHI*dadc/Denom - CPHI

      !Derivative DF/Dphi
      dFdChis(2) = p*CPHI + Num/Denom + COH*SPHI

      if (J <= C00000) then
        dFdChis(1) = - CPHI
        dFdChis(2) = p*CPHI + COH*SPHI
      end if

      end subroutine PartialFtoPartialChi


  Subroutine  DetermineDSigAndDEpsP(IntGlo,G,K,G_0, K_0,D1,D2, alpha_K, alpha_G, alpha_chi, alpha_beta &
			       ,Erate_vol0, Erate_vol, IErate0I, IErateI, RefRate,c,phi,psi, cp, phip, psip, EpsPEq, eta, Sig, &
					dEpsEqdEpsP ,dChisdEpsEq,dErate_eff, DEps,DSig,DEpsP, DchiErate)
      !**********************************************************************
      !
      ! Calculation of the stress increment and plastic strain increment
      !
      !         dSig = Dep * dEps
      !         dEpsP = Lambda * DPDSig
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i
      double precision :: A,Ai,Denom,LambdaNum,Lambda, aux
	  double precision :: p,J,Lode,S3TA, Gi, Ki !Invariants
      double precision, dimension(6,6) :: DE
      double precision, dimension(6) :: DummyVec, DummyVec2, Sigel
      double precision, dimension(6) :: dPdSig !Derivatives Plastic potential respect net stress
      double precision, dimension(6) :: DFDSig !Derivatives Yield function respect net stress
      double precision, dimension(2) :: dFdChis !Derivatives Yield function respect State Parameters
	  logical::ApplyStrainRateUpdates
      !In Variables
      double precision, intent(in) :: c,phi,psi, IErate0I, IErateI,Erate_vol0, Erate_vol, RefRate, G_0, K_0 !Softening parameters
	  double precision, intent(in) :: cp, phip, psip, EpsPEq, eta
	  double precision, intent(in) ::  alpha_K, alpha_G, alpha_chi, alpha_beta, dErate_eff
	  double precision, intent(inout) :: D1,D2,G,K !Elastic parameters
      double precision, intent(in), dimension(6):: dEpsEqdEpsP
      double precision, intent(in), dimension(6) :: Sig
      double precision, intent(in), dimension(3) :: dChisdEpsEq  !Derivatives respect Equivalent Plastic Strain
	  double precision, intent(in), dimension(6) :: DEps

      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
      !Out Variables
      double precision, intent(out), dimension(6) :: DSig
      double precision, intent(inout), dimension(6) :: DEpsP
	  double precision, intent(inout), dimension(3) :: DchiErate
		
      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
      call PartialForGtoPartialSigma(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,dPdSig)
      call PartialFtoPartialChi(p,J,Lode,S3TA,c,phi,dFdChis)

      !Parameter A (H = -A --> A>0 softening / A<0 hardening)
      A = 0.0d0
      Ai = (dFdChis(1)*dChisdEpsEq(1) + dFdChis(2)*dChisdEpsEq(2))
      do i=1,6
      A = A + Ai * dEpsEqdEpsP(i) * dPdSig(i)
	  end do
	  
	  Gi=G
	  Ki=K

	  !Parameter Hr and updating the elastic parameters
	  call check4crossing(IErate0I, IErateI, dErate_eff, RefRate,ApplyStrainRateUpdates)
	  if (ApplyStrainRateUpdates) then
	  !Elastic parameters updating with the strain rates
		  DchiErate=0.0d0
		call UpdatePardue2StrainRate(alpha_G, IErate0I, IErateI, dErate_eff, RefRate, G_0, Gi)
		
		DchiErate(1)=cp*alpha_chi*log10(IErateI/IErate0I)
		DchiErate(2)=phip*alpha_chi*log10(IErateI/IErate0I)
		DchiErate(3)=psip*alpha_beta*log10(IErateI/IErate0I)
	  end if
	  aux=Erate_vol-ERate_vol0
	  call check4crossing(ERate_vol0, Erate_vol, aux, RefRate,ApplyStrainRateUpdates)
	  if (ApplyStrainRateUpdates) then
		   call UpdatePardue2StrainRate(alpha_K,ERate_vol0, Erate_vol, dErate_eff, RefRate, K_0, Ki)
	  endif  
	 
		
		! Fill elastic material matrix
        D1  = Ki+(4*Gi/3)
        D2  = Ki-(2*Gi/3)
        DE  = 0.0
        DE(1:3,1:3) = D2
        DE(1,1) = D1
        DE(2,2) = D1
        DE(3,3) = D1
        DE(4,4) = Gi !Ask about value. Why is this not 2?
        DE(5,5) = Gi
        DE(6,6) = Gi

	  !Compute Sigel
	  Call MatVec( DE, 6, DEps, 6, Sigel)
	  !Computes D*dPdSigma
	  Call MatVec( DE, 6, dPdSig, 6, DummyVec2)

	  !Computes the numerators of the Lambda equation
	  call DotProduct(DFDSig, Sigel, 6, LambdaNum)
		if (ApplyStrainRateUpdates) then
			do i=1, 2
				LambdaNum=LambdaNum+dFdChis(i)*DChiErate(i)		
			enddo		
		end if
	  !Computes the denominator of the Lambda equation
	  call DotProduct(DFDSig, DummyVec2, 6, denom)
	  denom=denom+A

	  !Computes lambda
	  Lambda=LambdaNum/denom

	  !Computes Plastic strain increment Depsp=lambd*dP/dSigma
	  dEpsP=Lambda* dPdSig
	  
	  call MatVec( DE, 6, DEpsP, 6, DummyVec)
	  DSig=Sigel-DummyVec
	end subroutine DetermineDSigAndDEpsP

	subroutine EndOfStepCorrection(IntGlo,D1,D2,GG,IPL,F,Sig,dChisdEpsEq,dEpsEqdEpsP,EpsP,c,phi,psi,&
				DEpsP,DEps)
      !**********************************************************************
      !
      ! Final correction of the yield surface drift (END OF STEP CORRECTION).
      ! The stresses, the plastic strain and the strength parameters are corrected.
      !
      !**********************************************************************

      implicit none

      !Local variables
      integer :: i,k
      double precision :: p,J,Lode,S3TA !Invariants
      double precision :: Lambda,param, param2,c2,phi2,psi2,F2
      double precision :: Denom,A,Ai,Fact
      double precision, dimension(2) :: dFdChis
      double precision, dimension(6) :: dPdSig,DFDSig,Sig2,EpsP2
      double precision, dimension(6) :: Denom1

      double precision, dimension(3) :: Dh
      !In Variables
      integer, intent(in) :: IntGlo,IPL !Global ID of Gauss point or particle
      double precision, intent(in):: D1,D2,GG
      double precision, intent(in), dimension(3) :: dChisdEpsEq !Derivatives respect Equivalent Plastic Strain
      double precision, intent(in), dimension(6) :: dEpsEqdEpsP !Derivatives respect Equivalent Plastic Strain
	  double precision, intent(in), dimension(6) :: DEps
      !InOut Variables
      double precision, intent(inout):: c,phi,psi
      double precision, intent(inout), dimension(6) :: Sig
      double precision, intent(inout), dimension(6) :: EpsP, DEpsP
      double precision, intent(inout):: F

      call CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
      call PartialForGtoPartialSigma(Sig,p,J,Lode,S3TA,c,phi,psi,DFDSig,dPdSig)
      call PartialFtoPartialChi(p,J,Lode,S3TA,c,phi,dFdChis)

      !Parameter A (hardening/softening parameter)
      A = 0.0d0
      Ai = (dFdChis(1)*dChisdEpsEq(1) + dFdChis(2)*dChisdEpsEq(2))
      do i=1,6
        A = A + Ai * dEpsEqdEpsP(i) * dPdSig(i)
	  end do

	  Denom1(1) = dPdSig(1)*D1 + dPdSig(2)*D2 + dPdSig(3)*D2
      Denom1(2) = dPdSig(1)*D2 + dPdSig(2)*D1 + dPdSig(3)*D2
      Denom1(3) = dPdSig(1)*D2 + dPdSig(2)*D2 + dPdSig(3)*D1
      Denom1(4) = dPdSig(4)*GG
      Denom1(5) = dPdSig(5)*GG
      Denom1(6) = dPdSig(6)*GG

      Denom = Denom1(1)*DFDSig(1) + Denom1(2)*DFDSig(2) + &
             Denom1(3)*DFDSig(3) + Denom1(4)*DFDSig(4) + &
              Denom1(5)*DFDSig(5) + Denom1(6)*DFDSig(6) - A
	  Fact = 1d0/Denom

      Lambda = F*Fact !factor correction

      Sig2 = Sig - Lambda * Denom1 ! Sig2 = Sig + fact * Denom1 Stress corrected
      DEpsP = Lambda * dPdSig
      EpsP2 = EpsP + DEpsP

      if (IPL == 1)then
        Dh = 0.0d0
      else
        param = dEpsEqdEpsP(1) * DEpsP(1) + dEpsEqdEpsP(2) * DEpsP(2) + dEpsEqdEpsP(3) * DEpsP(3) + &
               dEpsEqdEpsP(4) * DEpsP(4) + dEpsEqdEpsP(5) * DEpsP(5) + dEpsEqdEpsP(6) * DEpsP(6)
		Dh(1) = dChisdEpsEq(1)*param
        Dh(2) = dChisdEpsEq(2)*param
        Dh(3) = dChisdEpsEq(3)*param

		Dh(1) = min (Dh(1), 0.0d0)
        Dh(2) = min (Dh(2), 0.0d0)
        Dh(3) = min (Dh(3), 0.0d0)
      end if

      c2 = c + Dh(1)
      phi2 = phi + Dh(2)
      psi2 = psi + Dh(3)

      call YieldFunctionValue(IntGlo,Sig2,c2,phi2,F2)

      if ((abs(F2) > abs(F)).or.(Denom == 0.0d0)) then !NormalCorrectionScheme
        Denom = 0.0d0
        Denom = DFDSig(1)*DFDSig(1) + DFDSig(2)*DFDSig(2) + &
                 DFDSig(3)*DFDSig(3) + DFDSig(4)*DFDSig(4) + &
                 DFDSig(5)*DFDSig(5) + DFDSig(6)*DFDSig(6)

        Lambda = F/Denom
        Sig = Sig - Lambda * DFDSig
        DEpsP = Lambda * dPdSig
        EpsP = EpsP + DEpsP
        call YieldFunctionValue(IntGlo,Sig,c,phi,F)
      else
        Sig = Sig2
        EpsP = EpsP2
        c = c2
        phi = phi2
        psi = psi2
        F = F2
      end if

	end subroutine EndOfStepCorrection

	Subroutine CalculateInvariants(IntGlo,Sig,p,J,Lode,S3TA)
      !**********************************************************************
      !
      ! Calcuation of the invariants (defined as Abbo & Sloan (1995))
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision :: Sx,Sy,Sz,SqTxy,SqTyz,SqTxz,suma,h1,h2,J2,J3
      double precision, parameter :: C00000 = 0.0D0
      double precision, parameter :: C00001 = 1.0D0
      double precision, parameter :: C00P16 = 0.166666666666666D0
      double precision, parameter :: C00002 = 2.0D0
      double precision, parameter :: C00003 = 3.0D0
      double precision, parameter :: CP3333 = 0.333333333333333D0
      double precision, parameter :: C00IR3 = 0.577350269189626D0
      double precision, parameter :: TINY = 0.000000000000001D0
      !In variables
      double precision, intent(in), dimension(6) :: Sig
      integer, intent(in) :: IntGlo !Global ID of Gauss point or particle
      !Out variables
      double precision, intent(out) :: p,J,Lode,S3TA !Invariants

      p = C00000
      J = C00000
      Lode = C00000

      !Mean stress (p)
      p = CP3333 * (Sig(1) + Sig(2) + Sig(3))

      !Deviatoric stress (J)
      Sx = Sig(1) - p
      Sy = Sig(2) - p
      Sz = Sig(3) - p
      suma = (Sig(1)-Sig(2))*(Sig(1)-Sig(2))+(Sig(1)-Sig(3))*(Sig(1)-Sig(3))+(Sig(2)-Sig(3))*(Sig(2)-Sig(3))
      SqTxy =  Sig(4) * Sig(4)
      SqTyz =  Sig(5) * Sig(5)
      SqTxz =  Sig(6) * Sig(6)

      J2 = C00P16 * suma + SqTxy + SqTyz + SqTxz
      J3 = Sx*Sy*Sz + C00002 * Sig(4)*Sig(5)*Sig(6) - Sx*SqTyz - Sy*SqTxz - Sz*SqTxy
      J = SQRT(J2)

      !Lode's angle (Lode)
      if (J2 > C00000) then

        h1 = -C00003/(C00002*C00IR3)
        h2 = J3/(J*J*J)
        S3TA = h1*h2
        if (S3TA < -C00001) then
            S3TA = -C00001
        else if (S3TA > C00001) then
            S3TA = C00001
      end if
        Lode = CP3333*asin(S3TA)
      else  !Special case of zero deviatoric stress
        Lode = C00000
        S3TA = C00000
      end if

	end subroutine CalculateInvariants


	subroutine CalculateEpsPEq(EpsP,EpsPEq)
      !**********************************************************************
      !
      ! Calculation of the equivalent plastic shear strain
      !
      !**********************************************************************

      implicit none

      !Local variables
      double precision:: EpsPM, C1, C2
      double precision, dimension(3) :: EpsDev
      !In variables
      double precision, intent(in), dimension(6) :: EpsP
      !Out variables
      double precision, intent(out) :: EpsPEq

      !EpsPEq = ((2/3)ep:ep)^(1/2), ep is the deviatoric plastic strain

      EpsPM = (1.0d0/3.0d0) * (EpsP(1) + EpsP(2) + EpsP(3))
      EpsDev(1) = EpsP(1)-EpsPM
      EpsDev(2) = EpsP(2)-EpsPM
      EpsDev(3) = EpsP(3)-EpsPM
      C1 = 2.0d0/3.0d0
      C2 = C1 * 0.5d0

      EpsPEq = sqrt(C1*EpsDev(1)*EpsDev(1) + C1*EpsDev(2)*EpsDev(2) +  C1*EpsDev(3)*EpsDev(3) + &
                     C2*EpsP(4)*EpsP(4) + C2*EpsP(5)*EpsP(5) + C2*EpsP(6)*EpsP(6))

	end subroutine CalculateEpsPEq

	subroutine ComputeRotationAtoB(DEpsel,DEpsP,kappa,NDEpsel,NDEpsP,PTOL,Q,cost)
	!**********************************************************************
	!
	! Calculation of the rotation matrix between DEpsP and DEpsel
	!
	!**********************************************************************

	implicit none
	!Local variables
	integer :: i, j
	double precision, dimension(6):: x, u, v
	double precision, dimension(6,6):: Id, uuT, vvT, Rot
	double precision, dimension(6,2):: uv
	double precision, dimension(2,6):: uvT
	double precision, dimension(2,2):: Rt
	double precision::normu, normv, Sint, Aux
	!In variables
	double precision, intent(in), dimension(6) ::DEpsel, DEpsP
	double precision, intent(in):: kappa, NDEpsel, NDEpsP, PTOL
	!Out variables
	double precision, intent(out)::Cost
	double precision, intent(out), dimension (6,6):: Q

	Call identitymatrix(6, Id)
	!cheack angle
	call DotProduct(DEpsel,DEpsP,6,Cost)
	Cost=Cost/(NDEpsel*NDEpsP)

	if (abs(Cost)<(1-PTOL)) then !Elastic and plastic strain increments are not parallel
		u= kappa * DEpsP
		call NormVec(u,6,normu) !determine the nirm of u
		u=(1/normu)*u !unit vector in the DEpsP direction

		call DotProduct(u,DEpsel,6,Aux) !u.DEpsel
		v=Aux*u !projection of DEpsel into u
		v=DEpsel-v !ortoghonal to u
		call NormVec(v,6,normv) ! norm of v
		v= (1/normv)*v ! unit vector orthogonal to DEpsP

		!Ensamble 2D rot matrix
		Sint=sqrt(1-Cost*Cost)
		Rt(1,1)=Cost
		Rt(2,2)=Cost
		Rt(1,2)=-Sint
		RT(2,1)=Sint

		!Proceed to apply Gram-Shmidt and rotartion
		call VecAOVecB(u,u,6,uuT) !uuT
		call VecAOVecB(v,v,6,vvT) !vvT

		do i=1,6
			uv(i,1)=u(i)
			uv(i,2)=v(i)
			uvT(1,i)=u(i)
			uvT(2,i)=v(i)
		end do

		call MatAMatB(uv,Rt,6,2,2,uv) !uv*Rt (6x2)(2x2)
		call MatAMatB(uv,uvT,6,2,6,Rot) !uv*Rt*uvT (6x2)(2x6)

		Q=Id-uuT-vvT+Rot !orthogonal rotational matrix

	else ! Paralelism can be assumed
		Q=Id
	end if


	end subroutine ComputeRotationAtoB

	Subroutine MatVec(xMat,IM,Vec,N,VecR)
!***********************************************************************
!
!     Calculate VecR = xMat*Vec
!
! I   xMat  : (Square) Matrix (IM,*)
! I   Vec   : Vector
! I   N     : Number of rows/colums
! O   VecR  : Resulting vector
!
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat(IM,*),Vec(*),VecR(*)
!***********************************************************************
      Do I=1,N
        X=0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
      Return
	End   Subroutine MatVec



	 Subroutine AddVec(Vec1,Vec2,R1,R2,N,VecR)
!***********************************************************************
!
!     Calculate VecR() = R1*Vec1()+R2*Vec2()
!
! I   Vec1,
! I   Vec2  : Vectors
! I   R1,R2 : Multipliers
! I   N     : Number of rows
! O   VecR  : Resulting vector
!
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
!***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
	End Subroutine AddVec

	Subroutine TwoNormTensor(Tensor, N, TwoNorm)
!***********************************************************************
!
!     Calculate 2NormTensor = sqrt(A:A)
!
! I   Tensor  : (Square or vecetor of dimension N)
! I   N     :   Number of elments
! O   2Norm : Resulting norm
!
!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension Tensor(N)
!***********************************************************************
	  X=N/2
      TwoNorm=0.0d0
	  Do I=1,X
		  TwoNorm=TwoNorm+Tensor(I)*Tensor(I)
	  end Do
	  Do I=X+1,N
		  TwoNorm=TwoNorm+0.5*(Tensor(I)*Tensor(I))!The convention in UMAT is to use engineering shear strains
	  end do
	  TwoNorm=sqrt(TwoNorm)

	end subroutine TwoNormTensor

		Subroutine MatAdpMatB(xMatA, xMatB,N, DProduct)
!***********************************************************************
!
!     Calculate the scalar produc of (A:B)
!
! I   xMatA xMATB  : (Tensors N)
! I   N     :   Number of elments
! O   DProduct : Resulting scalar produc
!
!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension xMatA(N), xMatB(N)
!***********************************************************************
	  X=N/2
      DProduct=0.0d0
	  Do I=1,X
		  DProduct=DProduct+xMatA(I)*xMatB(I)
	  end Do
	  Do I=X+1,N
		  DProduct=DProduct+2*xMatA(I)*xMatB(I)
	  end do

	end subroutine MatAdpMatB

	Subroutine identitymatrix(N,dIMat)
	!***********************************************************************
	!
	!     Create the Identity matrix
	!
	! I   N     :   Dimension
	! O   IMat : Resulting Identity matix
	!
	!***********************************************************************

	  Implicit Double Precision (A-H,O-Z)
      Dimension dIMat(N,N)

		 dIMat = 0.0d0                           ! Initialize the array.
		forall(j = 1:N) dIMat(j,j) = 1.0d0    ! Set the diagonal.

	end Subroutine identitymatrix

	Subroutine MatAMatB(xMatA, xMatB,ir1,ic1,ic2, AB)
!***********************************************************************
!
!     Calculate A(nxm)B(mxp) (AB)
!
! I   MatA MATB  : matrices
! I   r1     :   Number of rows in A
! I   c1     :   Number of columns in A
! I   c2     :   Number of columns in B
! O   AB : Resulting matrix
!
!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension xMatA(ir1,ic1), xMatB(ic1,ic2), AB(ir1,ic2)
!***********************************************************************


	  Do I=1,ir1
		  Do J=1,ic2
			  X=0
			  Do K=1,ic1
				  X=X+xMatA(I,K)*xMatB(K,J)
			  end do
			  AB(I,J)=X
		  end do
	  end do

	end subroutine MatAMatB

  Subroutine VecAOVecB(VecA, VecB,N, AOB)
!***********************************************************************
!
!     Calculate the inner product A(Nx1)B(1xN) (AOXB)
!
! I   VecA VecB  : Vectors
! I   N     :   Dimension
! O   AOB : Resulting product
!
!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension VecA(N), VecB(N), AOB(N,N)
!***********************************************************************

	  Do I=1,N
		  Do J=1,N
			  AOB(I,J)=VecA(I)*VecB(J)
		  end do
	  end do

	end subroutine VecAOVecB


  Subroutine DotProduct(VecA, VecB,N, Dp)
!***********************************************************************
!
!     Calculate the dot product of A(Nx1) and B(1xN)
!
! I   VecA VecB  : Vectors
! I   N     :   Dimension
! O   Dp : Dot product
!
!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension VecA(N), VecB(N)
!***********************************************************************
	  Dp=0.0d0
	  Do I=1,N
			  Dp=Dp+VecA(I)*VecB(I)
	  end do

	end subroutine DotProduct

Subroutine NormVec(VecA,N, xNorm)
!***********************************************************************
!
!     Calculate the norm of A
!
! I   VecA  : Vector
! I   N     :   Dimension
! O   Norm : norm of A
!
!***********************************************************************
	  Implicit Double Precision (A-H,O-Z)
      Dimension VecA(N)
	  double precision, intent(out):: xNorm
!***********************************************************************
	  Norm=0.0d0
	  call DotProduct(VecA, VecA, N, xNorm)
	  xNorm=sqrt(xNorm)

	end subroutine NormVec


	subroutine dbltobool(A,B)
	!******************************************************************
	! Takes a double which values are either 1.0 or 0.0 and returns a *
	! Boolean
	!******************************************************************
	implicit none
	double precision, intent(in):: A
	logical, intent(out):: B
	if (A<1.0) then
		B=.false.
	else
		B=.true.
	endif
	end subroutine dbltobool

	function logic2dbl(a)
	  logical, intent(in) :: a

	  if (a) then
		logic2dbl = 1.d0
	  else
		logic2dbl = 0.d0
	  end if
	end function logic2dbl
	
	Subroutine check4crossing(IErate0I, IErateI, dErate_eff,RateRef, Apply)
		implicit none
		double precision, intent(inout):: IErate0I, IErateI, dErate_eff, RateRef
		logical:: cond1, cond2, cond3
		logical, intent(out)::Apply
		  Apply=.false.
		  if(IErate0I<=RateRef) IErate0I=RateRef
		  if (IErateI<=RateRef) IErateI=RateRef
		  cond1=(IErate0I<=RateRef).and.(IErateI>RateRef)
		  cond2=(IErate0I>RateRef).and.(IErateI<=RateRef)
		  dErate_eff=IErateI-IErate0I
		  cond3=(IErate0I>=RateRef).and.(IErateI>=RateRef)
		  if (cond1.or.cond2.or.cond3) Apply=.true.
	end Subroutine check4crossing
