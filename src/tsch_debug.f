c
************************************************************************
c
c           SUBROUTINE  TSCHEM
c
c   This subroutine manages the whole integration scheme.
c   The following equations are solved:          
c   
c    ~~     n
c   Q  -  Q                n         n       n-1   alp       2  ~~   n 
c  --------- = -alp*grad (P ) + gam*H + rho*H   + ----- nabla ( Q + Q )
c    d t                                          2 Re
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
      REAL etime_,t(2),timef
      REAL start,finish
      REAL q1s(m1,m3),q3s(m1,m3) 
c
c   TIME INTEGRATION : implicit viscous, 3rd order RK (Adams Bashfort)  
c                                                                       
c       ttt_00 = etime_(t)
c       et_00  = timef()
      call cpu_time(start)
      do 2000 ns=1,nsst                                                 
c       ttt_in = etime_(t)
        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)
        write(6,*) ' ---------------------------------------- '
        write(6,*) ' T = ',time,' NTIME = ',ntime,' NS = ',ns
c
c       Axis Homogeneization
c
        k=1    ! NORTH POLE
        do j=1,n2m
         uaxN(1,j) = 0.
         uaxN(2,j) = 0.
         uaxN(3,j) = 0.
         do i=1,n1m
          is=isym(i)
          jm=jmv(j)
          utax = 0.25*(q1(i,j,1)-q1(is,j,1) 
     %               +q1(i,jm,1)-q1(is,jm,1))
          urax = 0.5*(q2(i,j,1)+q2(is,j,1))/(rc(j)**2)
          upax = 0.25*(q3(i,j,2)-q3(is,j,2)
     %                +q3(i,jm,2)-q3(is,jm,2))/sin(zz(2))
          uux = -utax*sin(thetac(i))+upax*cos(thetac(i))
          uuy = +utax*cos(thetac(i))+upax*sin(thetac(i))
          uuz = +urax
          uaxN(1,j) = uaxN(1,j) + uux
          uaxN(2,j) = uaxN(2,j) + uuy
          uaxN(3,j) = uaxN(3,j) + uuz
         end do
         uaxN(1,j) = uaxN(1,j)/float(n1m)
         uaxN(2,j) = uaxN(2,j)/float(n1m)
         uaxN(3,j) = uaxN(3,j)/float(n1m)
        end do
c
        k=N3!(m)    ! SOUTH POLE
        do j=1,n2m
         uaxS(1,j) = 0.
         uaxS(2,j) = 0.
         uaxS(3,j) = 0.
         do i=1,n1m
          is=isym(i)
          jm=jmv(j)
          utax = 0.25*(q1(i,j,n3m)-q1(is,j,n3m) 
     %               +q1(i,jm,n3m)-q1(is,jm,n3m))
          urax = 0.5*(q2(i,j,n3m)+q2(is,j,n3m))/(rc(j)**2)
          upax = 0.25*(q3(i,j,n3m)-q3(is,j,n3m)
     %                +q3(i,jm,n3m)-q3(is,jm,n3m))/sin(zz(n3m))
c         utax = 0.5*(q1(i,j,n3m)-q1(is,j,n3m))
c         urax = 0.5*(q2(i,j,n3m)+q2(is,j,n3m))/(rc(j)**2)
c         upax = 0.5*(q3(i,j,n3m)-q3(is,j,n3m))/sin(zz(n3m))
          uux = -utax*sin(thetac(i))-upax*cos(thetac(i))
          uuy = +utax*cos(thetac(i))-upax*sin(thetac(i))
          uuz = -urax
          uaxS(1,j) = uaxS(1,j) + uux
          uaxS(2,j) = uaxS(2,j) + uuy
          uaxS(3,j) = uaxS(3,j) + uuz
         end do
         uaxS(1,j) = uaxS(1,j)/float(n1m)
         uaxS(2,j) = uaxS(2,j)/float(n1m)
         uaxS(3,j) = uaxS(3,j)/float(n1m)
        end do
c
c
c    components in spherical coordinates
c
c    North pole axis
c
        do j=1,n2m
         jp=jpv(j)
         do i=1,n1m
          qaxN(1,i,j)=-0.5*(uaxN(1,j)+uaxN(1,jp))*sin(thetac(i))
     %                +0.5*(uaxN(2,j)+uaxN(2,jp))*cos(thetac(i))
          qaxN(2,i,j)=uaxN(3,j)
          qaxN(3,i,j)= 0.5*(uaxN(1,j)+uaxN(1,jp))*cos(thetam(i))
     %                +0.5*(uaxN(2,j)+uaxN(2,jp))*sin(thetam(i))
         end do
        end do
c
c    South pole axis
c
        do j=1,n2m
         jp=jpv(j)
         do i=1,n1m
          qaxS(1,i,j)=-0.5*(uaxS(1,j)+uaxS(1,jp))*sin(thetac(i))
     %                +0.5*(uaxS(2,j)+uaxS(2,jp))*cos(thetac(i))
          qaxS(2,i,j)=-uaxS(3,j)
          qaxS(3,i,j)=-0.5*(uaxS(1,j)+uaxS(1,jp))*cos(thetam(i))
     %                -0.5*(uaxS(2,j)+uaxS(2,jp))*sin(thetam(i))
         end do
        end do
c       et0=timef()
c       timhd0 = etime_(t)
        call hdnl1                    !! NON LINEAR TERMS
        call hdnl2                 !!  "    "      "
        call hdnl3                          !!  "    "      "
        call hdnl3                          !!  "    "      "
        call hdnlro                          !!  "    "      "
      do l=1,ndv
        vmax(l)=0.
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
              if(l.eq.1) vca=abs(dq(i,j,k))
              if(l.eq.2) vca=abs(dph(i,j,k))
              if(l.eq.3) vca=abs(qcap(i,j,k))
              if(vca.ge.vmax(l)) then
                vmax(l)=vca
                imaxv(l)=i
                jmaxv(l)=j
                kmaxv(l)=k
              end if
            end do
          end do
        end do
      end do
      write(6,*) ' hdnl1 ',vmax(1),imaxv(1),jmaxv(1),kmaxv(1)
      write(6,*) ' hdnl2 ',vmax(2),imaxv(2),jmaxv(2),kmaxv(2)
      write(6,*) ' hdnl3 ',vmax(3),imaxv(3),jmaxv(3),kmaxv(3)
        call invtr1                        !! DQ1HAT=Q1HAT-Q1(N)
        call invtr2                      !! DQ2HAT=Q2HAT-Q2(N)
        call invtr3                      !! DQ3HAT=Q3HAT-Q3(N)
        call invtrro                      !! DQ3HAT=Q3HAT-Q3(N)
c
c       timin3 = etime_(t)
c       write(6,*) ' t_in  = ',timin3-timin0
      do l=1,ndv
        vmax(l)=0.
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
              if(l.eq.1) vca=abs(q1(i,j,k))
              if(l.eq.2) vca=abs(q2(i,j,k))
              if(l.eq.3) vca=abs(q3(i,j,k))
              if(vca.ge.vmax(l)) then
                vmax(l)=vca
                imaxv(l)=i
                jmaxv(l)=j
                kmaxv(l)=k
              end if
            end do
          end do
        end do
      end do
      write(6,*) ' invtr1 ',vmax(1),imaxv(1),jmaxv(1),kmaxv(1)
      write(6,*) ' invtr2 ',vmax(2),imaxv(2),jmaxv(2),kmaxv(2)
      write(6,*) ' invtr3 ',vmax(3),imaxv(3),jmaxv(3),kmaxv(3)
      demax = 0.
      demin = 1.
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
              vca=dens(i,j,k)
              if(vca.ge.demax) then
                demax=vca
                imaxv(l)=i
                jmaxv(l)=j
                kmaxv(l)=k
              end if
              if(vca.le.demin) demin=vca
            end do
          end do
        end do
      write(6,*) ' invtrro ',demax,imaxv(1),jmaxv(1),kmaxv(1)
      write(6,*) ' invtrro ',demin
c       tiphc0  = etime_(t)
        call divg           
c       tiphc  = etime_(t)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c       write(6,*) ' t_dvg = ',tiphc-tiphc0
        vmax(1)=0.
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
               vca=abs(dph(i,j,k))
              if(vca.ge.vmax(1)) then
                vmax(1)=vca
                imaxv(1)=i
                jmaxv(1)=j
                kmaxv(1)=k
              end if
            end do
          end do
        end do
      write(6,*) ' divg ',vmax(1),imaxv(1),jmaxv(1),kmaxv(1)
c       et0=timef()
c       tiphc0  = etime_(t)
        call phcalc          !! PRESSURE (FFT & TRID)
c       tiphc  = etime_(t)
c       et1=timef()
c       write(6,*) ' Et_phc  = ',et1-et0
        vmax(1)=0.
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
               vca=abs(dph(i,j,k))
              if(vca.ge.vmax(1)) then
                vmax(1)=vca
                imaxv(1)=i
                jmaxv(1)=j
                kmaxv(1)=k
              end if
            end do
          end do
        end do
      write(6,*) ' phcal ',vmax(1),imaxv(1),jmaxv(1),kmaxv(1)
c       write(6,*) ' t_phc = ',(tiphc-tiphc0)*0.001
c       tiphc0  = etime_(t)
        call updvp                 !! SOLENOIDAL VEL FIELD
c       tiphc  = etime_(t)
c       write(6,*) ' t_upd = ',tiphc-tiphc0
      do l=1,ndv
        vmax(l)=0.
        do k=1,n3m
          do j=1,n2m
            do i=1,n1m
              if(l.eq.1) vca=abs(q1(i,j,k))
              if(l.eq.2) vca=abs(q2(i,j,k))
              if(l.eq.3) vca=abs(q3(i,j,k))
              if(vca.ge.vmax(l)) then
                vmax(l)=vca
                imaxv(l)=i
                jmaxv(l)=j
                kmaxv(l)=k
              end if
            end do
          end do
        end do
      end do
      write(6,*) ' updv1 ',vmax(1),imaxv(1),jmaxv(1),kmaxv(1)
      write(6,*) ' updv2 ',vmax(2),imaxv(2),jmaxv(2),kmaxv(2)
      write(6,*) ' updv3 ',vmax(3),imaxv(3),jmaxv(3),kmaxv(3)
      qmax = 0.
      qtot = 0.
      call divgck(qmax,qtot)
c       vmax(1)=0.
c       do k=1,n3m
c         do j=1,n2m
c           do i=1,n1m
c              vca=abs(dq(i,j,k))
c             if(vca.ge.vmax(1)) then
c               vmax(1)=vca
c               imaxv(1)=i
c               jmaxv(1)=j
c               kmaxv(1)=k
c             end if
c           end do
c         end do
c       end do
c     write(6,*) ' updvdq ',vmax(1),imaxv(1),jmaxv(1),kmaxv(1)
c       tiphc0  = etime_(t)
        call prcalc                         !! PRESSURE FIELD
c       tiphc  = etime_(t)
c       write(6,*) ' t_prc = ',tiphc-tiphc0
c       ttt_fi = etime_(t)
c       write(6,*) ' Total_Time = ',ttt_fi-ttt_in
 2000 continue                                                          
      call cpu_time(finish)
c       ttt_11 = etime_(t)
c       et_11 = timef()
        write(6,*) ' ---------------------------------------- '
        write(6,*) ' T = ',time,' NTIME = ',ntime,' DT = ',dt
        write(6,*) '      Total_Iter_Time = ',finish-start
c       write(6,*) ' Total_Elap_Time = ',(et_11-et_00)*0.001
      return                                                            
      end                                                               
c
