c
************************************************************************
c
c           SUBROUTINE  TSCHEM
c
c   This subroutine manages the whole integration scheme.
c   The following equations are solved:          
c   
c    ~~     n
c   Q  -  Q                n         n       n-1   alp       2  ~~   n    alp*Ek   ~~   n
c  --------- = -alp*grad (P ) + gam*H + rho*H   + ----- nabla ( Q + Q ) - ------ ( Q + Q )
c    d t                                          2 Re                      2
c
c          i                           i               i
c   where H  are the nonlinear terms, P  the pressure Q  the velocities
c       ~~
c   and Q  the provisional non solenoidal velocity field.
c   The superscripts (~~, n, n-1) indicate the time step level. 
c                        n
c   The nonlinear terms H  are computed in the routines HDNL*, while
c   in the routines INVTR* are computed the remaining terms, updated
c   the non linear terms and inverted the equation to find the provisional
c   field at the new time step.
c       ~~
c   The Q  velocity field is projected onto a solenoidal field by a 
c   scalar Phi computed through the equation
c
c                         2            1          ~~
c                    nabla (Phi ) =  ------ div ( Q  )
c                                    alp dt
c
c   The right hand side of this equation is computed in the routine
c   DIVG, while the equation is solved in PHCALC.
c
c   In the routine UPDVP the solenoidal velocity field at the new time
c   step is then computed through
c
c                n+1  ~~
c               Q   = Q  - alt*dt grad (Phi)
c
c   Finally in the routine PRCALC is updated the pressure field
c
c                n+1   n        alp dt      2
c               P   = P + Phi - ------ nabla (Phi)
c                                2 Re
c
c   When the scalar field is computed (density, concentration,
c   temperature) the routines HDNLRO and INVTRRO are used. The same
c   strategy at the velocity field is used, except that the scalar
c   field does not need any correction.
c
c   All variables are located on a staggered grid with the velocities
c   on the faces of the computational cell and all the scalars at the
c   centre. This is important when terms belonging to different equations
c   are avaluated.
c
c   Further details of the scheme can be found in the paper
c   "A finite-difference scheme for three-dimensional incompressible
c    flows in cylindrical coordinates" by R. Verzicco and P. Orlandi
c    J. of Comp. Phys. 1996.
c
c
      subroutine tschem(time)
      include'param.f'
c 
c     REAL etime_,t(2),timef
c     REAL start,finish
      integer start,finish,rate
c     REAL q1s(m1,m3),q3s(m1,m3)
c
c   TIME INTEGRATION : implicit viscous, 3rd order RK (Adams Bashfort)  
c                                                                       
c       ttt_00 = etime_(t)
c       et_00  = timef()
      call system_clock(start, rate)
      do 2000 ns=1,nsst                                                 
c       ttt_in = etime_(t)
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)
c       write(6,*) ' ---------------------------------------- '
c       write(6,*) ' T = ',time,' NTIME = ',ntime,' NS = ',ns
c       timhd0 = etime_(t)
        call hdnl1                 !! NON LINEAR TERMS
        call hdnl2                 !!  "    "      "
        call hdnl3                 !!  "    "      "
c       call hdnlro                !!  "    "      "
        call invtr1                !! DQ1HAT=Q1HAT-Q1(N)
        call invtr2                !! DQ2HAT=Q2HAT-Q2(N)
        call invtr3                !! DQ3HAT=Q3HAT-Q3(N)
c       call invtrro               !! 
c
c       tiphc0  = etime_(t)
        call divg           
c       tiphc  = etime_(t)
c       tiphc0  = etime_(t)
        call phcalc                !! PRESSURE (FFT & TRID)
c       tiphc  = etime_(t)
c       tiphc0  = etime_(t)
        call updvp                 !! SOLENOIDAL VEL FIELD
c       tiphc  = etime_(t)
c     qmax = 0.
c     qtot = 0.
c     call divgck(qmax,qtot)
c       tiphc0  = etime_(t)
        call prcalc                !! PRESSURE FIELD
c       write(6,*) ' Total_Time = ',ttt_fi-ttt_in
 2000 continue                                                          
c       ttt_11 = etime_(t)
c       et_11 = timef()
        call system_clock(finish, rate)
        write(6,*) ' ---------------------------------------- '
        write(6,*) ' T = ',time,' NTIME = ',ntime,' DT = ',dt
c       write(6,*) '      Total_Iter_Time = ',finish-start
        write(6,*) '      Total_Iter_Time = ',(finish-start)/rate
c       write(6,*) ' Total_Elap_Time = ',(et_11-et_00)*0.001
      return                                                            
      end                                                               
c
