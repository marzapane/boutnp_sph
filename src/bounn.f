c     this code is made for simulating three-dimensional flows in polar
c     ( cilindrical ) coordinates.
c     boundary condition are slip-walls ( =0 ) or non-slip wall ( =1)
c     in the radial and the vertical direction.
c
c     this code allows , by changing suitable indices , to introduce
c     the buoiancy term ( ibuo=1 ) the coriolis terms ( icorio=1 )
c     and to consider the transport of a scalar ( isca=1 )
c     if the scalar is passiv, the term for the centrifug force can
c     be included in the pressure so the index ifugo can be set =0,
c     otherwise it has to be set =1 (if the scalar is not a passiv one)
c
c      navier-stokes equations are solved by a fractional step method
c     ( Kim and Moin ) with the pressure in the first step.
c
c     The time advancement of the solution is obtained by a
c     Runge-Kutta 3rd order low storage scheme (Wray) or a 2nd
c     order Adams-Bashfort scheme.
c
c     The Poisson  equation for the pressure is solved directly
c     introducing FFT in the azimutal and vertical direction.
c     Because of the presence of the walls in the vertical direction
c     the cosFFT in that direction is required
c
c                 Roberto Verzicco and Paolo Orlandi
c                 dipartimento di meccanica ed aeronautica
c                 universita' la sapienza di roma
c
c
c
c     All variables are calculated in a staggered grid:
c
c        Instead of velocities, "flux" variables are introduced to avoid
c        the singularity for r=0
c        q1=vtheta, q2= r*vr, q3= vz
c
c        dq1,dq2,dq3 : velocity correction
c
c        qcap :divergence of the  non free divergent velocity field
c
c        non linear terms:
c
c        ru1, ru2, ru3, ruro : old step
c        h1,h2,h3, hro : new step
c        dens : density
c        pr : pressure
c        dph : pressure correction
c       pressure solver is in the package press.f
c       non linear terms are calculated in hdnl routines
c       the invertions of momentum equations is performed in invtr
c       routines
c
      program papero
      include'param.f'
      character*4 dummy
c
c     timings
      real timef
      real tot_time(2)    ! total time
c
c
c     tot_time(1)=timef()
c
      open(unit=15,file='bou.in',status='old')
        read(15,301) dummy
        read(15,*) n1,n2,n3,nsst,nwrit,nread
        read(15,301) dummy
        read(15,*) n1p,n2p,n3p
        read(15,301) dummy
        read(15,*) ntst,tprint,tpin,tmax,ireset
        read(15,301) dummy
        read(15,*) alx3,istr3,str3,rmed31,etdp3,strb3
        read(15,301) dummy
        read(15,*)strr,rext,rint,istr,rmed,etdp,strb
        read(15,301) dummy
        read(15,*) ray,ri,ros,pra,dt,resid,cflmax
        read(15,301) dummy
        read(15,*) inslwri,inslwr,ifugo,icorio,ibuo,isca,ivisc
        read(15,301) dummy
        read(15,*) epsil,lamb,rper,nsect
        read(15,301) dummy
        read(15,*) istat,timeav,timerms,iresta
        read(15,301) dummy
        read(15,*) nlev,dism,nn
        read(15,301) dummy
        read(15,*) idtv,dtmax,cfllim
        read(15,301) dummy
        read(15,*) nini,nfin,nstri
        read(15,301) dummy
        read(15,*) nson
        read(15,301) dummy
        read(15,*) indexForInitial
        read(15,301) dummy
        read(15,*) iforc
301     format(a4)
      close(15)
      r0 = rext
      if(nson.ne.0) then
      open(unit=15,file='sonda.in',status='old')
        do i = 1,nson
          read(15,*) (coson(i,j),j=1,3)
        end do
      close(15)
      end if
c
c
c     DEFINITIONS FOR THE NATURAL CONVECTION
c
      ren = sqrt(ray/pra)
      pec = sqrt(pra*ray)
c                          
c
c
      call openfi
c
      pi=2.*asin(1.)
      tfini=dt*ntst
      tfin=tfini
      n1m=n1-1
      if(n1.eq.1) n1m = 1
      n2m=n2-1
      n3m=n3-1
      n3mh=n3m/2+1
      if(n1m.eq.1) then
        iaxsy=0
      else
        iaxsy=1
      endif
c
      write(6,112)2.*r0/alx3
      write(32,112)2.*r0/alx3
  112 format(//,20x,'R A Y B E N ',//,10x,
     % 'Cilindro a rapporto d''aspetto:  D/H = ',f4.2)
      write(6,202) ray,pra
      write(32,202) ray,pra
  202 format(/,5x,'Parametri: ',' Ray=',e9.3,' Pra= ',e9.3)
      if(idtv.eq.1) then
         write(6,204) cflmax
         write(32,204) cflmax
  204 format(/,5x,'calcolo effettuato con dt variabile e cfl= ',
     % e10.4,/ )
      else
         write(6,205) dtmax,cfllim
         write(32,205) dtmax,cfllim
  205 format(/,5x,'calcolo effettuato con dt= ',e10.4,' e cflmax=
     %  ',e10.4,/ )
      endif
c
c     assign coefficients for time marching schemes
c
      if(nsst.gt.1) then
        gam(1)=8./15.
        gam(2)=5./12.
        gam(3)=3./4.
        rom(1)=0.
        rom(2)=-17./60.
        rom(3)=-5./12.
        write(6,100)
  100   format(/,5x,'The time scheme is a III order Runge-Kutta')
        write(6,101) (gam(n),n=1,nsst),(rom(n),n=1,nsst)
  101   format(/,5x,'gam= ',3f8.3,4x,'ro= ',3f8.3)
      else
        gam(1)=1.5
        gam(2)=0.
        gam(3)=0.
        rom(1)=-0.5
        rom(2)=0.
        rom(3)=0.
        write(6,110)
  110   format(/,5x,'The time scheme is the Adams-Bashfort')
        write(6,111) gam(1),rom(1)
  111   format(/,5x,'gam= ',f8.3,4x,'ro= ',f8.3)
      endif
      do 10 ns=1,nsst
        alm(ns)=(gam(ns)+rom(ns))
   10 continue
c
c     the solution of the problem starts
c
      call gcurv
c
c
c     tot_time(2)=timef()
      write(6,*) ' Total Time (all simulation) ='
     &             ,(tot_time(2)-tot_time(1))*0.001
c
      stop
      end
c
c***********************************************************************
      subroutine openfi
      include'param.f'
      character*14 filth,filte,filba,filen,filso,filve,filnu,filbaper
      common/names/filth,filte,filba,filen,filso,filve,filnu,filbaper
c     open(46,file='nfbou',status='old')
c     read(46,'(a)')filth
c     read(46,'(a)')filte
c     read(46,'(a)')filba
c     read(46,'(a)')filen
c     read(46,'(a)')filso
      filth = 'bou.out'
      filte = 'boute.out'
c     filba = 'bouba.out'
      filen = 'bouen.out'
c     filso = 'paso'
c     filbaper = 'filbaper.out'
      open(32,file=filth,status='unknown')
c     open(34,file=filba,status='unknown')
c     open(35,file=filbaper,status='unknown')
c     open(39,file=filen,status='unknown')
      open(97,file=filte,status='unknown')
c     filnu='nusse.out'
      filve='maxvel.out'
c     open(95,file=filnu,status='unknown')
      open(96,file=filve,status='unknown')
      rewind(32)
c     rewind(34)
c     rewind(35)
c     rewind(39)
      rewind(95)
      rewind(97)
c     close(46)
      return
      end
c
c
c***********************************************************************
      subroutine divgck(qmax,qtot)
      include'param.f'
C     real   usrnu3,udx2,fac2,dqcap
c
      dvol= 1./(dx1*dx2*dx3)
      qtot=0.
      qmax=0.
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,ipv,usrm,usinpm)
!$OMP$ SHARED(q1,q2,q3,dx1,udx2m,udx3m,rm)
!$OMP$ PRIVATE(kp,jp,ip,dqcap,ic,jc,kc,fact)
!$OMP$ REDUCTION(max:qmax)
!$OMP$ REDUCTION(+:qtot)
      do 11 kc=1,n3m
        kp=kc+1
        do 11 jc=1,n2m
          jp=jc+1
          fact = rm(jc)*usinpm(kc)
            do 11 ic=1,n1m
              ip=ipv(ic)
              dqcap =(
     %               (q2(ic,jp,kc)-q2(ic,jc,kc))*udx2m(jc)
     %              +((q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
     %              + (q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc))*fact
     %               )
              qmax = max(abs(dqcap),qmax)
c             if(abs(dqcap).gt.qmax) then
c               imxq=ic
c               jmxq=jc
c               kmxq=kc
c               qmax = abs(dqcap)
c             endif
              qtot=qtot+dqcap*rm(jc)*g2rm(jc)*g3rm(kc)
   11 continue
!$OMP  END PARALLEL DO
      qqmax=qmax
      qtot=qtot*dvol
      write(6,*) ' QMAX = ',qmax
c     write(6,*) ' QMAX = ',qmax,' @ i,j,k = ',imxq,jmxq,kmxq
      write(6,*) ' QTOT = ',qtot
      return
      end
c
c***********************************************************************
      subroutine meshes
      include'param.f'
      real   pi
      pi=2.d0*asin(1.d0)
      if(nsect.eq.0) then
      dx1=2.*pi/float(n1m)
        else
      dx1=2.*pi/(float(n1m)*float(lamb))
        endif
      dx2=1.d0/float(n2m)
!     dx3=pi/(float(n3m))
      dx3=1.d0/(float(n3m))
      write(6,*)'dx1, dx2, dx3  ',dx1,dx2,dx3
      dx1=1.d0/dx1
      dx2=1.d0/dx2
      dx3=1.d0/dx3
      dx1q=dx1*dx1
      dx2q=dx2*dx2
      dx3q=dx3*dx3
      return
      end
c
c***********************************************************************
      subroutine vmaxv
      include'param.f'
c
        vmax1=0.
        vmax2=0.
        vmax3=0.
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,q1,q2,q3)
!$OMP$ PRIVATE(i,j,k)
!$OMP$ REDUCTION(max:vmax1)
!$OMP$ REDUCTION(max:vmax2)
!$OMP$ REDUCTION(max:vmax3)
        do 312 k=1,n3m
          do 312 j=2,n2m
            do 312 i=1,n1m
              vmax1 = max(vmax1,abs(q1(i,j,k)))
              vmax2 = max(vmax2,abs(q2(i,j,k)))
              vmax3 = max(vmax3,abs(q3(i,j,k)))
  312   continue
!$OMP  END PARALLEL DO
        vmax(1) = vmax1
        vmax(2) = vmax2
        vmax(3) = vmax3
        return
      end
c
c
c**************************************************************
c
      subroutine cordin
      include'param.f'
      real   eta(m2),etaz(m3),etarm(2*m2)
      real   xt31(m3),xt32(m3)
      real   tstr,tstr2,tstr3,tstr4
      real   vk(m3),rk(m2)
c     real   rint,rext,x3,alx3
c     real   xn11,r0,x2dp,rcdp,xx1,xx2
      pi=2.d0*asin(1.d0)
c
c     AZIMUTHAL COORDINATE DEFINITION
c
      do i=1,n1
        thetac(i)= float(i-1)/dx1
      end do
      do i=1,n1m
        thetam(i)= (float(i-1)+0.5d0)/dx1
      end do
c
c     RADIAL COORDINATE DEFINITION
c
      tstr= tanh(strr)
      tstr2=tanh(strr*0.5d0)
      tstr3=tanh(strr*(etdp-0.5d0))
      tstr4=tanh(strb*(etdp-1.d0))
      xn11= rmed * tanh(strr*0.5d0)/tstr3
      rext = r0
c
c     UNIFORM GRID
c
      if (istr.eq.-1) then
        do  j=1,n2
         x2=float(j-1)/float(n2m)
         rc(j)= rint + x2*(rext-rint)
        end do
      end if
c
c     CLUSTERING AT THE AXIS
c
      if (istr.eq.0) then
        tstr2=tanh(strr)
        do  j=1,n2
          x2=float(j-1)/float(n2m)
          eta(j)=1.d0 + tanh(strr*(x2-1.d0))/tstr2
          rc(j)=rint + eta(j)*(rext-rint)
        enddo
      end if
c
c     CLUSTERING AT THE EXTERNAL RADIAL WALL
c
      if (istr.eq.1) then
        do j=1,n2
          x2dp=float(j-1)/float(n2m)
          rcdp=tanh(strr*x2dp)/tstr
          rc(j)=rint+rcdp*(rext-rint)
        end do
      end if
c
c     VARIABLE CLUSTERING
c
      if (istr.eq.2) then
        do j=1,n2
          x2=(j-1)/dx2
          x22 = 0.5d0 + x2*0.5d0
          xx1= rmed * tanh(strr*(x22- 0.5d0))/tstr3
          xx2= 1./xn11 + (1.d0 - 1.d0/xn11)* tanh(strb*(x22-1.d0))/tstr4
          eta(j)= xx1*xx2
          rc(j)=rint+eta(j)*(rext-rint)
        end do
      end if
      rc(1)=rint
C
C    MESH FROM AN EXTERNAL INPUT FILE
C
      if (istr.eq.3) then
        open(66,file='gridr.data',form='formatted')
        do j=1,n2
          read(66,*) nx2,rc(j)
        end do
        close(66)
      end if
C
C     MODIFIED CHANNEL TYPE CLUSTERING
C
      if (istr.eq.4) then
      nclip = 1  ! 7 used at Ra=2x10^11
c     nclip = 14
      n2mo = n2+nclip+nclip
      do j=1,n2mo
        etarm(j)=+cos(pi*(float(j)-0.5)/float(n2mo))
        write(200,*) j,etarm(j)
      end do
      close(200)
      do j=1,n2
        eta(j)=etarm(j+nclip)
      end do
      delet = eta(1)-eta(n2)
      etain = eta(1)
      do j=1,n2
        eta(j)=eta(j)/(0.5*delet)
      end do
      rc(1) = rint
      do j=2,n2m
        rc(j) = rint+(rext-rint)*(1.-eta(j))*0.5
      end do
      rc(n2) = rext
      end if
C
C     CHANNEL TYPE CLUSTERING
C
      if (istr.eq.5) then
      tstr=tanh(strr*0.5)
      do j=1,n2
        x2=float(j-1)/float(n2m)
        eta(j)=0.5d0*(1.d0+tanh(strr*(x2-0.5d0))/tstr)
      end do
      do j=1,n2
c       zz(k)=(-0.5+etaz(k))*2.
        rc(j)=rint+(rext-rint)*eta(j)
      end do
      end if
c
c
c     STAGGERED COORDINATES AND
c     METRIC QUANTITIES
c
      do j=1,n2m
        rm(j)=(rc(j)+rc(j+1))*0.5d0
        g2rm(j)=(rc(j+1)-rc(j))*dx2
      end do
      do j=2,n2m
        g2rc(j)=(rc(j+1)-rc(j-1))*dx2*0.5d0
      end do
      g2rc(1)=(rc(2)-rc(1))*dx2
      g2rc(n2)=(rc(n2)-rc(n2m))*dx2
c
c     AXIAL COORDINATE DEFINITION
c
c
      tstr3=tanh(str3)
c
c     UNIFORM GRID
c
c       write(6,*) istr3,alx3
      if (istr3.lt.0) then
        do k=1,n3
          zz(k)= alx3*float(k-1)/dx3
        end do
        do k=1,n3m
          zm(k)= alx3*(float(k-1)+0.5d0)/dx3
        end do

      endif
c
c     CLUSTERING AT THE SOUTH WALL
c
      if (istr3.eq.0) then
        do k=1,n3
          x3=float(k-1)/float(n3m)
          etaz(k)=alx3*(1.d0+tanh(str3*(x3-1.d0))/tstr3)
          zz(k)=etaz(k)
        enddo
      endif
c
c     CLUSTERING AT THE NORTH WALL
c
      if (istr3.eq.1) then
        do k=1,n3
          x3=float(k-1)/float(n3m)
          etaz(k)=alx3*tanh(str3*x3)/tstr3
          zz(k)=etaz(k)
        enddo
      endif
c
c     VARIABLE CLUSTERING
c
      if (istr3.eq.2) then
        dx3t=n3m
        tstr3=tanh(str3*etdp3)
        do k=1,n3
          x3=float(k-1)/dx3t
          xt31(k)=rmed31/alx3*tanh(str3*x3)/tstr3
        end do
        do k=1,n3
          x3=(k-1)/dx3t
          xt32(k)=1./xt31(n3)+(1.d0-1.d0/xt31(n3))
     %    *tanh(strb3*(x3-1.d0))/tanh(strb3*(etdp3-1.d0))
          etaz(k)= xt31(k)*xt32(k)
          zz(k)=etaz(k)*alx3
        end do
      end if
C
C     CHANNEL TYPE CLUSTERING
C
      if (istr3.eq.5) then
      tstr3=tanh(str3*0.5)
      do k=1,n3
        x3=float(k-1)/float(n3m)
        etaz(k)=0.5d0*(1.d0+tanh(str3*(x3-0.5d0))/tstr3)
      end do
      do k=1,n3
c       zz(k)=(-0.5+etaz(k))*2.
        zz(k)=+etaz(k)
      end do
      end if
C
C     CHEBYCHEV TYPE CLUSTERING
C
      if (istr3.eq.7) then
      do k=1,n3
        etaz(k)=+cos(pi*(float(k)-0.5)/float(n3))
        write(200,*) k,etaz(k)
      end do
      close(200)
      zz(1) = 0.
      do k=2,n3m
        zz(k) = alx3*(1.-etaz(k))*0.5
      end do
      zz(n3) = alx3
      end if
c
c    CLUSTERED AT THE CENTER
c
      if (istr3 .EQ. 8) then
      do k=1,n3
          x3=float(k-1)/float(n3m)
           if (x3 .LE. 0.5d0) then
           etaz(k)=0.5d0*tanh(str3*x3)/tanh(str3*0.5d0)
           zz(k)=alx3*etaz(k)
           else
              etaz(k)=1.d0+0.5d0*tanh(str3*(x3-1))/tanh(str3*0.5d0)
              zz(k)=alx3*etaz(k)
           end if
c        write(6,*) k,zz(k)
       end do
       end if
c      pause
C
C    MESH FROM AN EXTERNAL INPUT FILE
C
      if (istr3.eq.3) then
        open(66,file='gridp.data',form='formatted')
        do k=1,n3
          read(66,*) nx3,zz(k)
        end do
        close(66)
      end if
c
c     STAGGERED COORDINATES AND
c     METRIC QUANTITIES
c
      do k=1,n3m
        zm(k)=(zz(k)+zz(k+1))*0.5d0
        g3rm(k)=(zz(k+1)-zz(k))*dx3
      enddo
      do k=2,n3m
        g3rc(k)=(zz(k+1)-zz(k-1))*dx3*0.5d0
      enddo
      g3rc(1)=(zz(2)-zz(1))*dx3
      g3rc(n3)= (zz(n3)-zz(n3m))*dx3
c
c     WRITE GRID INFORMATION
c
c     open(unit=98,file='spacor.out',status='unknown')
c     do i=1,n1
c       write(98,345) i,thetac(i),thetam(i)
c     end do
c     close(98)
c     open(unit=98,file='radcor.out',status='unknown')
c     do j=1,n2
c       write(98,345) j,rc(j),rm(j),g2rc(j),g2rm(j)
c     end do
c     close(98)
c     open(unit=78,file='axicor.out',status='unknown')
c     do k=1,n3
c       write(78,345) k,zz(k),zm(k),g3rc(k),g3rm(k)
c     end do
c     close(78)
c345  format(i4,4(2x,e13.5))
c
c     QUANTITIES FOR DERIVATIVES
c
c     open(unit=78,file='fact3.out',status='unknown')
       usinpm(1) = 1./sin(zm(1))
       utanpm(1) = 1./tan(zm(1))
       sinpm(1) = sin(zm(1))
       sinpc(1) = sin(zz(1))
       cospm(1) = cos(zm(1))
       cospc(1) = cos(zz(1))
c      write(999,*) zm(1),sin(zm(1)),cos(zm(1)),tan(zm(1))
      do k=2,n3m
       utanpc(k) = 1./tan(zz(k))
       utanpm(k) = 1./tan(zm(k))
       usinpm(k) = 1./sin(zm(k))
       usinpc(k) = 1./sin(zz(k))
       sinpm(k) = sin(zm(k))
       sinpc(k) = sin(zz(k))
       cospm(k) = cos(zm(k))
       cospc(k) = cos(zz(k))
c      write(999,*) zm(k),sin(zm(k)),cos(zm(k)),tan(zm(k))
      end do
c     close(999)
       sinpc(n3) = sin(zz(n3))
       cospc(n3) = cos(zz(n3))
      do k=1,n3m
        udx3m(k) = dx3/g3rm(k)
        udx3c(k) = dx3/g3rc(k)
        udx3vis2(k) = 2.*dx3/g3rm(k)/ren
c       write(78,*) k,udx3m(k),udx3c(k)
      end do
      udx3c(n3) = dx3/g3rc(n3)
c       write(78,*) n3,udx3m(n3m),udx3c(n3)
c       close(78)
c     open(unit=78,file='fact2m.out',status='unknown')
      do j=1,n2m
        rmsg2rm(j) = rm(j)**2/g2rm(j)
        usg2rc(j) = 1.d0/g2rc(j)
        udx2rm(j) = dx2/(g2rm(j)*rm(j))
        udx2rmq(j) = dx2/(g2rm(j)*rm(j)**2)
        udx2m(j) = dx2/g2rm(j)
        udx2c(j) = dx2/g2rc(j)
        usrm(j) = 1.d0/rm(j)
        usrmq(j) = 1.d0/rm(j)**2
        udx1qm(j) = dx1q/rm(j)**2
c       write(78,*) j,udx2rm(j),udx2m(j),udx2c(j),
c    %               usrm(j),udx1qm(j)
      end do
      usg2rc(n2) = 1.d0/g2rc(n2)
      udx2c(n2) = dx2/g2rc(n2)
c       write(78,*) n2,udx2c(j)
c       close(78)
c     open(unit=78,file='fact2c.out',status='unknown')
        usrc(1) = 1.d0
        udx1qc(1) = 1.d0
      do j=2,n2
        usrc(j) = 1.d0/rc(j)
        udx1qc(j) = dx1q/rc(j)**2
c       write(78,*) j,udx1qc(j),usrc(j)
      end do
c       close(78)
      do k=1,n3m
      do j=1,n2m
       udx1vis1(j,k)=2./(ren*rm(j)**2*sin(zm(k)))
       udx1vis3(j,k)=2.*cos(zm(k))/(ren*(rm(j)*sin(zm(k)))**2)
      end do
      end do
c
c   viscous indices
c
c------------------------------------------
      do k=1,n3
         vk(k) = 1.
      end do
      small = 1. !0.0002
      angin=0.25 !zz(4)
      angfi=zz(n3)-0.25 !zz(n3-3)
      do k=1,n3
        if((zz(k).ge.angin).and.(zz(k).le.angfi)) then
         vk(k) = 1.
        end if
        if(zz(k).le.angin) then
         vk(k) = small + (1.-small)*(
     %         3.*(zz(k)/angin)**2-2.*(zz(k)/angin)**3)
        end if
        if(zz(k).gt.angfi) then
         vk(k) = small + (1.-small)*(
     %          3.*((zz(n3)-zz(k))/(zz(n3)-angfi))**2
     %         -2.*((zz(n3)-zz(k))/(zz(n3)-angfi))**3)
        end if
c       write(695,*) zz(k),vk(k)
      end do
c     close(695)
      do j=1,n2
        rk(j) = 1.
      end do
      small = 1. !0.0002
      ragin=0.25
      do j=1,n2m
        if(rc(j).ge.ragin) then
          rk(j) = 1.
        end if
        if(rc(j).lt.ragin) then
          rk(j) = small + (1.-small)*(
     %          3.*(rc(j)/ragin)**2-2.*(rc(j)/ragin)**3)
        end if
c       write(694,*) rc(j),rk(j)
      end do
c     close(694)
c------------------------------------------
      do k=1,n3
      do j=1,n2
       fvisc(j,k) = 1.
       fvisc1(j,k) = 1.
       fvisc2(j,k) = 1.
       fvisc3(j,k) = 1.
       fvisc1e(j,k) = 1.
       fvisc2e(j,k) = 1.
       fvisc3e(j,k) = 1.
      end do
      end do
      do k=1,n3
c      fvisc1(1,k) = small**2 !vk(k)*rk(1)
c      fvisc2(1,k) = small**2 !vk(k)*rk(1)
c      fvisc3(1,k) = small**2 !vk(k)*rk(1)
c      fvisc(1,k) = small**2 !vk(k)*rk(1)
      do j=1,n2
       if((zz(k).le.angin).or.(zz(k).ge.angfi)) then
       if(rc(j).le.ragin) then
       fvisc1(j,k) = vk(k)*rk(j)
       fvisc2(j,k) = vk(k)*rk(j)
       fvisc3(j,k) = vk(k)*rk(j)
       if(j.eq.1) fvisc3(j,k) = 0.
       fvisc(j,k) =  vk(k)*rk(j)
       end if
       end if
c      fvisc1(j,1) = 0.01
c      fvisc2(j,1) = 0.01
c      fvisc3(j,2) = 0.01
c      fvisc1(j,2) = 0.1
c      fvisc2(j,2) = 0.1
c      fvisc3(j,3) = 0.1
c      fvisc1(j,n3m) = 0.01
c      fvisc2(j,n3m) = 0.01
c      fvisc3(j,n3m) = 0.01
c      fvisc1(j,n3m-1) = 0.1
c      fvisc2(j,n3m-1) = 0.1
c      fvisc3(j,n3m-1) = 0.1
      end do
      end do
      do k=1,n3
      do j=1,n2
       fvisc1e(j,k) = fvisc1(j,k)   
       fvisc2e(j,k) = fvisc2(j,k)   
       fvisc3e(j,k) = fvisc3(j,k)   
      end do
      end do
      return
      end
c
c
c***********************************************************************
      subroutine continua(time)
      include'param.f'
      character*70 filcnw
C     filcnw = './dati/continua.dat'
c     filcnw = '/scratch/verzicco/barun/dati/continua.dat'
C     open(13,file=filcnw,form='unformatted',status='unknown')
C       rewind(13)
C       write(13)n1,n2,n3
C       write(13) epsil,lamb,re,time
C       write(13) (((dens(i,j,k),i=1,n1),j=1,n2),k=1,n3),
C    %          (((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3),
C    %          (((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3),
C    %          (((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c    %         ,(((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
C     close(13)
      filcnw = './dati/continua.dat'
c     filcnw = '/scratch/verzicco/barun/dati/continua.dat'
      open(13,file=filcnw,form='unformatted',status='unknown')
        rewind(13)
        write(13)n1,n2,n3
        write(13) epsil,lamb,re,time
        write(13) (((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3)
        write(13) (((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3)
        write(13) (((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
        write(13) (((dens(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      close(13)
c       open(13,file='densin.dat',form='unformatted',status='unknown')
c       rewind(13)
c       write(13)n1,n2,n3
c       write(6,*) ' n1,n2,n3 ',n1,n2,n3
c       write(13) appa,na,appa,time
c       write(6,*) appa,na,appa,time
c       write(13)
c     close(13)
c     open(13,file='q1in.dat',form='unformatted',status='unknown')
c       rewind(13)
c       write(13)n1,n2,n3
c       write(13) appa,na,appa,time
c       write(13)
c    %          (((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     close(13)
c     open(13,file='q2in.dat',form='unformatted',status='unknown')
c       rewind(13)
c       write(13)n1,n2,n3
c       write(13) appa,na,appa,time
c       write(13)
c    %          (((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     close(13)
c     open(13,file='q3in.dat',form='unformatted',status='unknown')
c       rewind(13)
c       write(13)n1,n2,n3
c       write(13) appa,na,appa,time
c       write(13)
c    %          (((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     close(13)
      return
      end
c
c***********************************************************************
      subroutine inirea(ntii,time)
      include'param.f'
      character*70 filcnw
c     real appa
C     filcnw = './dati/continua.dat'
c     filcnw = '/scratch/verzicco/barun/dati/continua.dat'
C     open(13,file=filcnw,form='unformatted',status='unknown')
C       rewind(13)
C       read(13)n1,n2,n3
C       write(6,*)' ===>  ',n1,n2,n3
C       read(13) appa,na,appa,time
C       write(6,*) ' ===>  ',appa,na,appa,time
C       read(13)
C    %          (((dens(i,j,k),i=1,n1),j=1,n2),k=1,n3),
C    %          (((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3),
C    %          (((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3),
C    %          (((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c    %         ,(((pr(i,j,k),i=1,n1),j=1,n2),k=1,n3)
C     close(13)
      filcnw = './dati/continua.dat'
      open(13,file=filcnw,form='unformatted',status='unknown')
        rewind(13)
        read(13)n1,n2,n3
        read(13) epsil,lamb,re,time
        read(13) (((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3)
        read(13) (((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3)
        read(13) (((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
        read(13) (((dens(i,j,k),i=1,n1),j=1,n2),k=1,n3)
      close(13)
c     open(13,file='densin.dat',form='unformatted',status='unknown')
c       rewind(13)
c       read(13)n1,n2,n3
c       write(6,*) ' n1,n2,n3 ',n1,n2,n3
c       read(13) appa,na,appa,time
c       write(6,*) appa,na,appa,time
c       read(13)
c     close(13)
c     open(13,file='q1in.dat',form='unformatted',status='unknown')
c       rewind(13)
c       read(13)n1,n2,n3
c       read(13) appa,na,appa,time
c       read(13)
c    %          (((q1(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     close(13)
c     open(13,file='q2in.dat',form='unformatted',status='unknown')
c       rewind(13)
c       read(13)n1,n2,n3
c       read(13) appa,na,appa,time
c       read(13)
c    %          (((q2(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     close(13)
c     open(13,file='q3in.dat',form='unformatted',status='unknown')
c       rewind(13)
c       read(13)n1,n2,n3
c       read(13) appa,na,appa,time
c       read(13)
c    %          (((q3(i,j,k),i=1,n1),j=1,n2),k=1,n3)
c     close(13)
c
c    Reconstruction of the q2 velocity component
c    by the free-divergence condition
c
C     do k=1,n3m
C       usrnu3=dx3/g3rm(k)
C       kp=k+1
C       do i=1,n1m
C         ip=ipv(i)
C         q2(i,1,k) = 0.
C         do j=2,n2m
C           jm=j-1
C           udx2=dx2/g2rm(jm)
C           usurm = 1./rm(jm)
C           q2(i,j,k)=q2(i,jm,k)-(
C    %             (q1(ip,jm,k)-q1(i,jm,k))*dx1*usurm
C    %            +(q3(i,jm,kp)-q3(i,jm,k))*usrnu3
C    %                                       )/(udx2*usurm)
C         end do
C       end do
C     end do
      if (ireset.eq.1) then
      ihist=0
      ntii=0
      time=0.
      endif
      return
      end
c
c***********************************************************************
      subroutine initia
      include'param.f'
!$OMP  PARALLEL DO
!$OMP$ SHARED(pr,dens,q1,q2,q3,dq,ru1,ru2,ru3,ruro,hro,rhs,dph,qcap)
!$OMP$ PRIVATE(i,j,k,iadd)
      do 4 k=1,m3
        do 4 j=1,m2
          do 4 i=1,m1
            pr(i,j,k)=0.
            q1(i,j,k)=0.
            q2(i,j,k)=0.
            q3(i,j,k)=0.
            dq(i,j,k)=0.
            ru1(i,j,k)=0.
            ru2(i,j,k)=0.
            ru3(i,j,k)=0.
            iadd = i + (j-1)*m1+(k-1)*m1*m2
            rhs(iadd) = 0.
            dph(i,j,k) = 0.
            qcap(i,j,k) = 0.
            dens(i,j,k)=0.
            ruro(i,j,k)=0.
            hro(i,j,k)=0.
    4 continue
!$OMP  END PARALLEL DO
      do j=1,m2
        amscj(j) = 0.
        acscj(j) = 0.
        apscj(j) = 0.
        rc(j) = 0.
        rm(j) = 0.
        g2rc(j) = 0.
        g2rm(j) = 0.
      end do
      do k=1,m3
        am3ssk(k) = 0.
        ac3ssk(k) = 0.
        ap3ssk(k) = 0.
        zz(k) = 0.
        zm(k) = 0.
        g3rc(k) = 0.
        g3rm(k) = 0.
      end do
      return
      end
c
c***********************************************************************
      subroutine outpf2(time)
      include'param.f'
      character*70 namfi3
      character*5 ipfi
      real vx(m1,m2,m3),vy(m1,m2,m3),vz(m1,m2,m3)
cccccccccccccccccccccc
      go to 3476
c     Spherical components
      do kc=2,n3m
       km=kmv(kc)
       kp=kpv(kc)
        do jc=2,n2
        jm=jmv(jc)
        jp=jpv(jc)
         do ic=1,n1m
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         isy = isym(ic)
c Cartesian coordinates
          xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
          yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
          zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vx(ic,jc,kc) = 0.25*(q1(ic,jc,kc)+q1(ic,jm,kc)
     %              +q1(ic,jc,km)+q1(ic,jm,km))
          vy(ic,jc,kc) = 0.25*(q2(ic,jc,kc)+q2(im,jc,kc)
     %               +q2(ic,jc,km)+q2(im,jc,km))/rc(jc)**2
          vz(ic,jc,kc) = 0.25*(q3(ic,jc,kc)+q3(im,jc,kc)
     %              +q3(ic,jm,kc)+q3(im,jm,kc))/sin(zz(kc)) ! C L secondo me qui va messo KC, JM
            end do
          end do
        end do
c variables at the sphere centre
      do kc=2,n3m
       km=kmv(kc)
       kp=kpv(kc)
        jc=1
        jm=jmv(jc)
        jp=jpv(jc)
         do ic=1,n1m
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         ksyy = ksym(km)
         ksyx = ksyc(kc)
         ksymc = ksyc(km)
         isy = isym(ic)
         isyx=isym(im)
          xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
          yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
          zzz = rc(jc)*cos(zz(kc))
          vx(ic,jc,kc) = 0.25*(q1(ic,jc,kc)-q1(isy,jc,ksy)
     %              +q1(ic,jc,km)-q1(isy,jc,ksyy))
          vy(ic,jc,kc) = 0.125*(q2(ic,jp,kc)-q2(isy,jp,ksy)
     %               +q2(im,jp,kc)-q2(isyx,jp,ksy)
     %               +q2(ic,jp,km)-q2(isy,jp,ksyy)
     %               +q2(im,jp,km)-q2(isyx,jp,ksyy))/rc(jp)**2
          vz(ic,jc,kc) = 0.25*(q3(ic,jc,kc)+q3(isy,jc,ksyx)
     %              +q3(im,jc,kc)+q3(isyx,jc,ksyx))/sin(zz(kc))
            end do
        end do
c North pole
       kc=1
       kp=kpv(kc)
        do jc=2,n2
        jm=jmv(jc)
        jp=jpv(jc)
         do ic=1,n1m
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         ksc = ksyc(kc)
         isy = isym(ic)
         isyx = isym(im)
c Cartesian coordinates
          xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
          yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
          zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vx(ic,jc,kc) = 0.25*(q1(ic,jc,kc)-q1(isy,jc,kc)
     %              + q1(ic,jm,kc)-q1(isy,jm,kc))
          vy(ic,jc,kc) = 0.25*(q2(ic,jc,kc)+q2(isy,jc,kc)
     %          +q2(im,jc,kc)+q2(isyx,jc,kc) )/rc(jc)**2
          vz(ic,jc,kc) = 0.125*(q3(ic,jc,kp)-q3(isy,jc,kp)
     %     +q3(im,jc,kp)-q3(isyx,jc,kp)
     %     +q3(ic,jm,kp)-q3(isy,jm,kp)
     %     +q3(im,jm,kp)-q3(isyx,jm,kp)  )/sin(zz(kp))
            end do
        end do
c South pole
       kc=n3
       km=n3m
        do jc=2,n2
        jm=jmv(jc)
        jp=jpv(jc)
         do ic=1,n1m
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         ksc = ksyc(kc)
         isy = isym(ic)
         isyx = isym(im)
c Cartesian coordinates
          xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
          yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
          zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vx(ic,jc,kc) = 0.25*(q1(ic,jc,km)-q1(isy,jc,km)
     %              + q1(ic,jm,km)-q1(isy,jm,km))
          vy(ic,jc,kc) = 0.25*(q2(ic,jc,km)+q2(isy,jc,km)
     %          +q2(im,jc,km)+q2(isyx,jc,km) )/rc(jc)**2

          vz(ic,jc,kc) = 0.125*(q3(ic,jc,km)-q3(isy,jc,km)
     %     +q3(im,jc,km)-q3(isyx,jc,km)
     %     +q3(ic,jm,km)-q3(isy,jm,km)
     %     +q3(im,jm,km)-q3(isyx,jm,km)  )/sin(zz(km))
            end do
        end do
!        pause
         do ic=1,n1
      vx(ic,1,1) = vx(ic,1,2)
      vy(ic,1,1) = vy(ic,1,2)
      vz(ic,1,1) = vz(ic,1,2)
      vx(ic,1,n3) = vx(ic,1,n3m)
      vy(ic,1,n3) = vy(ic,1,n3m)
      vz(ic,1,n3) = vz(ic,1,n3m)
        end do
      do kc=1,n3
       do jc=1,n2
        pr(n1,jc,kc) = pr(1,jc,kc)
        vx(n1,jc,kc) = vx(1,jc,kc)
        vy(n1,jc,kc) = vy(1,jc,kc)
        vz(n1,jc,kc) = vz(1,jc,kc)
       end do
      end do
      itime=1000.*time !+0.5
      write(ipfi,99) itime
      namfi3='./dati/bosp'//ipfi//'.q'
c     namfi3='/scratch/verzicco/barun/dati/bou'//ipfi//'.dat'
      write(6,201)time,namfi3
      n1pp=(n1-1)/n1p!+1
      n2pp=(n2-1)/n2p!+1
      n3pp=(n3-1)/n3p!+1
      open(99,file=namfi3,form='unformatted',status='unknown')
        rewind(99)
        write(99) 1
        write(99) n1,n2,n3
        write(99) epsil,epsil,re,time
        write(99)
c    %   ((((1.),
     %   ((((sqrt(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vx(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vy(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vz(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   ((((pr(i,j,k)),
     %   i=1,n1),j=1,n2),k=1,n3)
      close(99)
cccccccccccccccccccccc
 3476  continue
c     Cartesian components
      do kc=2,n3m
       km=kmv(kc)
       kp=kpv(kc)
        do jc=2,n2
        jm=jmv(jc)
        jp=jpv(jc)
         do ic=1,n1m
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         isy = isym(ic)
c Cartesian coordinates
          xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
          yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
          zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vt = 0.25*(q1(ic,jc,kc)+q1(ic,jm,kc)
     %              +q1(ic,jc,km)+q1(ic,jm,km))
          vr = 0.25*(q2(ic,jc,kc)+q2(im,jc,kc)
     %               +q2(ic,jc,km)+q2(im,jc,km))/rc(jc)**2
          vp = 0.25*(q3(ic,jc,kc)+q3(im,jc,kc)
     %              +q3(ic,jm,kc)+q3(im,jm,kc))/sin(zz(kc)) ! C L secondo me qui va messo KC, JM
          ct = cos(thetac(ic))
          st = sin(thetac(ic))
          cp = cos(zz(kc))
          sp = sin(zz(kc))
c Cartesian velocities
      vx(ic,jc,kc) = vr*ct*sp-vt*st*sp+vp*ct*cp
      vy(ic,jc,kc) = vr*st*sp+vt*ct*sp+vp*st*cp
      vz(ic,jc,kc) = vr*cp-vp*sp
            end do
          end do
        end do
c variables at the sphere centre
      do kc=2,n3m
       km=kmv(kc)
       kp=kpv(kc)
        jc=1
        jm=jmv(jc)
        jp=jpv(jc)
         do ic=1,n1m
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         ksyy = ksym(km)
         ksyx = ksyc(kc)
         ksymc = ksyc(km)
         isy = isym(ic)
         isyx=isym(im)
          xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
          yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
          zzz = rc(jc)*cos(zz(kc))
          vt = 0.25*(q1(ic,jc,kc)-q1(isy,jc,ksy)
     %              +q1(ic,jc,km)-q1(isy,jc,ksyy))
          vr = 0.125*(q2(ic,jp,kc)-q2(isy,jp,ksy)
     %               +q2(im,jp,kc)-q2(isyx,jp,ksy)
     %               +q2(ic,jp,km)-q2(isy,jp,ksyy)
     %               +q2(im,jp,km)-q2(isyx,jp,ksyy))/rc(jp)**2
          vp = 0.25*(q3(ic,jc,kc)+q3(isy,jc,ksyx)
     %              +q3(im,jc,kc)+q3(isyx,jc,ksyx))/sin(zz(kc))
          ct = cos(thetac(ic))
          st = sin(thetac(ic))
          cp = cos(zz(kc))
          sp = sin(zz(kc))
      vx(ic,jc,kc) = vr*ct*sp-vt*st*sp+vp*ct*cp
      vy(ic,jc,kc) = vr*st*sp+vt*ct*sp+vp*st*cp
      vz(ic,jc,kc) = vr*cp-vp*sp
            end do
        end do
c North pole
       kc=1
       kp=kpv(kc)
        do jc=2,n2
        jm=jmv(jc)
        jp=jpv(jc)
         do ic=1,n1m
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         ksc = ksyc(kc)
         isy = isym(ic)
         isyx = isym(im)
c Cartesian coordinates
          xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
          yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
          zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vt = 0.25*(q1(ic,jc,kc)-q1(isy,jc,kc)
     %              + q1(ic,jm,kc)-q1(isy,jm,kc))
          vr = 0.25*(q2(ic,jc,kc)+q2(isy,jc,kc)
     %          +q2(im,jc,kc)+q2(isyx,jc,kc) )/rc(jc)**2
          vp = 0.125*(q3(ic,jc,kp)-q3(isy,jc,kp)
     %     +q3(im,jc,kp)-q3(isyx,jc,kp)
     %     +q3(ic,jm,kp)-q3(isy,jm,kp)
     %     +q3(im,jm,kp)-q3(isyx,jm,kp)  )/sin(zz(kp))
          ct = cos(thetac(ic))
          st = sin(thetac(ic))
          cp = cos(zz(kc))
          sp = sin(zz(kc))
c Cartesian velocities
      vx(ic,jc,kc) = vr*ct*sp-vt*st*sp+vp*ct*cp
      vy(ic,jc,kc) = vr*st*sp+vt*ct*sp+vp*st*cp
      vz(ic,jc,kc) = vr*cp-vp*sp
            end do
        end do
c South pole
       kc=n3
       km=n3m
        do jc=2,n2
        jm=jmv(jc)
        jp=jpv(jc)
         do ic=1,n1m
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         ksc = ksyc(kc)
         isy = isym(ic)
         isyx = isym(im)
c Cartesian coordinates
          xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
          yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
          zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vt = 0.25*(q1(ic,jc,km)-q1(isy,jc,km)
     %              + q1(ic,jm,km)-q1(isy,jm,km))
          vr = 0.25*(q2(ic,jc,km)+q2(isy,jc,km)
     %          +q2(im,jc,km)+q2(isyx,jc,km) )/rc(jc)**2

          vp = 0.125*(q3(ic,jc,km)-q3(isy,jc,km)
     %     +q3(im,jc,km)-q3(isyx,jc,km)
     %     +q3(ic,jm,km)-q3(isy,jm,km)
     %     +q3(im,jm,km)-q3(isyx,jm,km)  )/sin(zz(km))

          ct = cos(thetac(ic))
          st = sin(thetac(ic))
          cp = cos(zz(kc))
          sp = sin(zz(kc))
c Cartesian velocities
      vx(ic,jc,kc) = vr*ct*sp-vt*st*sp+vp*ct*cp
      vy(ic,jc,kc) = vr*st*sp+vt*ct*sp+vp*st*cp
      vz(ic,jc,kc) = vr*cp-vp*sp

C          write(6, *) 'I=',ic,' J= ',jc,'K= ',kc
C          write(6,*) q1(ic,jc,km)
C          write(6,*) q2(ic,jc,km)
C          write(6,*) q3(ic,jc,km)
C          write(6,*) vt, vr, vp
C          write(6,*) vx(ic,jc,kc)
C          write(6,*) vy(ic,jc,kc)
C          write(6,*) vz(ic,jc,kc)
C          write(6,*) sin(zz(km))
C          pause
            end do
        end do
!        pause
         do ic=1,n1
      vx(ic,1,1) = vx(ic,1,2)
      vy(ic,1,1) = vy(ic,1,2)
      vz(ic,1,1) = vz(ic,1,2)
      vx(ic,1,n3) = vx(ic,1,n3m)
      vy(ic,1,n3) = vy(ic,1,n3m)
      vz(ic,1,n3) = vz(ic,1,n3m)
        end do
      do kc=1,n3
       do jc=1,n2
        pr(n1,jc,kc) = pr(1,jc,kc)
        vx(n1,jc,kc) = vx(1,jc,kc)
        vy(n1,jc,kc) = vy(1,jc,kc)
        vz(n1,jc,kc) = vz(1,jc,kc)
       end do
      end do
cccccccccccccccccccccc
c     do k=1,n3m
c      do j=1,n2m
c       do i=1,n1m
c        ip = ipv(i)
c        ct= cos(thetam(i))
c        st= sin(thetam(i))
c        cp= cos(zm(k))
c        sp= sin(zm(k))
c        vx(i,j,k) = 0.5*(q2(i,j,k)+q2(i,j+1,k))*ct*sp/rm(j)**2-
c    %               0.5*(q1(i,j,k)+q1(ip,j,k))*st*sp*rm(j)+
c    %               0.5*(q3(i,j,k)+q3(i,j,k+1))*ct*cp*rm(j)/sp
c        vy(i,j,k) = 0.5*(q2(i,j,k)+q2(i,j+1,k))*st*sp/rm(j)**2+
c    %               0.5*(q1(i,j,k)+q1(ip,j,k))*ct*sp*rm(j)+
c    %               0.5*(q3(i,j,k)+q3(i,j,k+1))*st*cp*rm(j)/sp
c        vz(i,j,k) = 0.5*(q2(i,j,k)+q2(i,j+1,k))*cp/rm(j)**2-
c    %               0.5*(q3(i,j,k)+q3(i,j,k+1))*sp*rm(j)/sp
c       end do
c      end do
c     end do
      itime=10.*time !+0.5
      write(ipfi,99) itime
   99 format(i5.5)
      namfi3='./dati/boum'//ipfi//'.q'
c     namfi3='/scratch/verzicco/barun/dati/bou'//ipfi//'.dat'
      write(6,201)time,namfi3
  201 format(10x,'At t=',e10.3,' write data on ',a20)
      n1pp=(n1-1)/n1p!+1
      n2pp=(n2-1)/n2p!+1
      n3pp=(n3-1)/n3p!+1
      open(99,file=namfi3,form='unformatted',status='unknown')
        rewind(99)
        write(99) 1
        write(99) n1,n2,n3
        write(99) epsil,epsil,re,time
        write(99)
c    %   ((((1.),
     %   ((((sqrt(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vx(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vy(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vz(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
c    %   ((((0.5*(q1(i,j,k)+q1(i+1,j,k))),
c    %   i=1,n1m),j=1,n2m),k=1,n3m),
c    %   ((((0.5*(q2(i,j,k)+q2(i,j+1,k))),
c    %   i=1,n1m),j=1,n2m),k=1,n3m),
c    %   ((((0.5*(q3(i,j,k)+q3(i,j,k+1))),
c    %   i=1,n1m),j=1,n2m),k=1,n3m),
     %   ((((pr(i,j,k)),
     %   i=1,n1),j=1,n2),k=1,n3)
      close(99)
C      namfi3='SXY_'//ipfi//'.dat'
C      open(69,file=namfi3,form='formatted',status='unknown')
C      rewind(69)
C      write(69,*)'Variables = X, Y, Z, Vt, Vr, Vp, Pr'
C      write(69,*)'ZONE F=POINT, I= ',n1,', J= ',n2m, ', K= ',1
       kc=n3m/2+1
       km=kmv(kc)
       kp=kpv(kc)
        do jc=1,n2m
        jm=jmv(jc)
        jp=jpv(jc)
         do ii=1,n1
         ic = ii
         if(ii.eq.n1) ic=1
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         isy = isym(ic)
          xxx = rm(jc)*cos(thetam(ic))*sin(zm(kc))
          yyy = rm(jc)*sin(thetam(ic))*sin(zm(kc))
          zzz = rm(jc)*cos(zm(kc))
          vt = 0.5*(q1(ic,jc,kc)+q1(ip,jc,kc))
          if(jc.eq.1) then
           q2sy = 0.5*(q2(ic,jp,kc)-q2(isy,jp,ksy))/rc(jp)**2
           vr = 0.5*(q2sy+q2(ic,jp,kc)/rc(jp)**2)
          else
           vr = 0.5*(q2(ic,jc,kc)/rc(jc)**2+q2(ic,jp,kc)/rc(jp)**2)
          end if
          if((kc.ne.1).and.(kc.ne.n3m)) then
          vp = 0.5*(q3(ic,jc,kc)/sin(zz(kc))+q3(ic,jc,kp)/sin(zz(kp)))
          else if(kc.eq.1) then
C L           ksy = ksyc(kp)
           ksy=kp
           q3sy = 0.5*(q3(ic,jc,kp)/sin(zz(kp))-
     %                 q3(isy,jc,ksy)/sin(zz(ksy)))
          vp = 0.5*(q3sy+q3(ic,jc,kp)/sin(zz(kp)))
          else if(kc.eq.n3m) then
C L           ksy = ksyc(kc)
           ksy=kc
           q3sy = 0.5*(q3(ic,jc,kc)/sin(zz(kc))-
     %                 q3(isy,jc,ksy)/sin(zz(ksy)))
          vp = 0.5*(q3sy+q3(ic,jc,kc)/sin(zz(kc)))
          end if
          ct = cos(thetam(ic))
          st = sin(thetam(ic))
          cp = cos(zm(kc))
          sp = sin(zm(kc))
          vvx = vr*ct*sp-vt*rm(jc)*st*sp+vp*rm(jc)*ct*cp
          vvy = vr*st*sp+vt*rm(jc)*ct*sp+vp*rm(jc)*st*cp
          vvz = vr*cp-vp*rm(jc)*sp
C                      write(69,*)xxx,yyy,zzz
C     %           ,vt
C     %           ,vr
C     %           ,vp
C     %           ,sqrt(vt**2+vr**2+vp**2)
c    %           ,vvx
c    %           ,vvy
c    %           ,vvz
c    %           ,sqrt(vvx**2+vvy**2+vvz**2)
            end do
          end do
c       do jc=1,n2m
c       jm=jmv(jc)
c       jp=jpv(jc)
c        do ii=1,n1
c        ic = ii
c        if(ii.eq.n1) ic=1
c        im=imv(ic)
c        ip=ipv(ic)
c         xxx = rm(jc)*cos(thetam(ic))*sin(zm(kc))
c         yyy = rm(jc)*sin(thetam(ic))*sin(zm(kc))
c         zzz = rm(jc)*cos(zm(kc))
c                     write(69,*)xxx,yyy,zzz
c    %           ,0.5*(q1(ic,jc,kc)+q1(ip,jc,kc))
c    %           ,0.5*(q2(ic,jc,kc)+q2(ic,jp,kc))/rm(jc)**2
c    %           ,0.5*(q3(ic,jc,kc)+q3(ic,jc,kp))/sin(zm(kc))
c    %           ,pr(ic,jc,kc)
c           end do
c         end do
c     do kc=1,n3m
c      km=kmv(kc)
c      kp=kpv(kc)
c      do jc=1,n2m
c       jm=jmv(jc)
c       jp=jpv(jc)
c        ic=n1m/2+1
c        im=imv(ic)
c        ip=ipv(ic)
c         xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
c         yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
c         zzz = rc(jc)*cos(zz(kc))
c                     write(69,*)xxx,yyy,zzz
c    %           ,0.25*(q1(ic,jc,kc)+q1(ic,jm,kc)+
c    %                  q1(ic,jm,km)+q1(ic,jc,km))
c    %           ,0.25*(q2(ic,jc,kc)+q2(im,jc,kc)+
c    %                  q2(im,jc,km)+q2(ic,jc,km))
c    %           ,0.25*(q3(ic,jc,kc)+q3(im,jc,kc)+
c    %                  q3(im,jm,kc)+q3(ic,jm,kc))
c    %           ,0.125*(pr(ic,jc,kc)+pr(im,jc,kc)+
c    %                   pr(im,jm,kc)+pr(ic,jm,kc)+
c    %                   pr(ic,jc,km)+pr(im,jc,km)+
c    %                   pr(im,jm,km)+pr(ic,jm,km))
c           end do
c         end do
C!      close(69)
      namfi3='CXY_'//ipfi//'.dat'
C!      open(69,file=namfi3,form='formatted',status='unknown')
C!      rewind(69)
C!      write(69,*)'Variables = X, Y, Z, Vx, Vy, Vz, Pr'
C!      write(69,*)'ZONE F=POINT, I= ',n1,', J= ',n2m, ', K= ',1
       kc=n3m/2+1
       km=kmv(kc)
       kp=kpv(kc)
        do jc=1,n2m
        jm=jmv(jc)
        jp=jpv(jc)
         do ii=1,n1
         ic = ii
         if(ii.eq.n1) ic=1
         im=imv(ic)
         ip=ipv(ic)
         ksy = ksym(kc)
         isy = isym(ic)
          xxx = rm(jc)*cos(thetam(ic))*sin(zm(kc))
          yyy = rm(jc)*sin(thetam(ic))*sin(zm(kc))
          zzz = rm(jc)*cos(zm(kc))
          vt = 0.5*(q1(ic,jc,kc)+q1(ip,jc,kc))
          if(jc.eq.1) then
           q2sy = 0.5*(q2(ic,jp,kc)-q2(isy,jp,ksy))/rc(jp)**2
           vr = 0.5*(q2sy+q2(ic,jp,kc)/rc(jp)**2)
          else
           vr = 0.5*(q2(ic,jc,kc)/rc(jc)**2+q2(ic,jp,kc)/rc(jp)**2)
          end if
          if((kc.ne.1).and.(kc.ne.n3m)) then
          vp = 0.5*(q3(ic,jc,kc)/sin(zz(kc))+q3(ic,jc,kp)/sin(zz(kp)))
          else if(kc.eq.1) then
C L           ksy = ksyc(kp)
           ksy = kp
           q3sy = 0.5*(q3(ic,jc,kp)/sin(zz(kp))-
     %                 q3(isy,jc,ksy)/sin(zz(ksy)))
          vp = 0.5*(q3sy+q3(ic,jc,kp)/sin(zz(kp)))
          else if(kc.eq.n3m) then
C L           ksy = ksyc(kc)
           ksy=kc
           q3sy = 0.5*(q3(ic,jc,kc)/sin(zz(kc))-
     %                 q3(isy,jc,ksy)/sin(zz(ksy)))
          vp = 0.5*(q3sy+q3(ic,jc,kc)/sin(zz(kc)))
          end if
          ct = cos(thetam(ic))
          st = sin(thetam(ic))
          cp = cos(zm(kc))
          sp = sin(zm(kc))
C!                      write(69,*)xxx,yyy,zzz
C!     %           ,vr*ct*sp-vt*rm(jc)*st*sp+vp*rm(jc)*ct*cp
C!     %           ,vr*st*sp+vt*rm(jc)*ct*sp+vp*rm(jc)*st*cp
C!     %           ,vr*cp-vp*rm(jc)*sp
C!     %           ,pr(ic,jc,kc)
            end do
          end do
C!      close(69)
      return
      end
c
c
c
c***********************************************************************
c                                                                      *
c   ********************* subr fftv                                    *
c   this subroutine initialize the temperton fft                       *
c                                                                      *
c***********************************************************************
      subroutine fftv
      include 'param.f'
      n1mh=n1m/2
      n1mp=n1mh+1
      nv1fft=n1m
      call fftfax(nv1fft,ifxv,trigxv)
      return
      end
c***********************************************************************
c                                                                      *
c   ************** WSON ************************                       *
c                                                                      *
c***********************************************************************
      subroutine wson(time)
      include 'param.f'
      character*3 ipfis
      character*70 namson
      character*14 filth,filte,filba,filen,filso,filve,filnu,filbaper
      common/names/filth,filte,filba,filen,filso,filve,filnu,filbaper
c
        nul=150
      do i = 1,nson
        ir = mason(i,1)
        jr = mason(i,2)
        kr = mason(i,3)
        irm = imv(ir)
        jrm = jmv(jr)
        krm = kmv(kr)
        krp = kpv(kr)
         if(jr.eq.1) then
           rrr= 1.
         else
           rrr=rc(jr)
         endif
        udx3=dx3/(g3rc(kr)*rrr)
        udx2=dx2/g2rc(jr)
c
c       azimuthal vorticity
c
              dq2x3=(q2(ir,jr,kr)-q2(ir,jr,krm))*udx3
              dq3x2=(q3(ir,jr,kr)-q3(ir,jrm,kr))*udx2
              voraz=dq2x3-dq3x2
c
c       radial vorticity
c
              dq3x1=(q3(ir,jr,kr)-q3(irm,jr,kr))*dx1/rm(jr)*iaxsy
              dq1x3=(q1(ir,jr,kr)-q1(ir,jr,krm))*dx3/g3rc(kr)
              vorr=dq3x1-dq1x3
c
c       axial vorticity
c
              dq1x2=(q1(ir,jr,kr)*rm(jr)-q1(ir,jrm,kr)*rm(jrm))
     %           *udx2
              dq2x1=(q2(ir,jr,kr)-q2(irm,jr,kr))*dx1/rrr*iaxsy
              vorz=(dq1x2-dq2x1)/rrr
        delq3=+q3(ir,jr,krp)-q3(ir,jr,krm)
        q2sur = q2(ir,jr,kr)/rrr
        write(ipfis,99) nul
   99 format(i3.3)
      namson='son'//ipfis//'.out'
c     namson='/scratch/verzicco/barun/son'//ipfis//'.out'
      if (time.gt.0.) then
            open(99,file=namson,status='old',position='append')
      else
            open(99,file=namson,status='replace')
      end if
c     if(ntime.eq.1) go to 30
c20   read(99,100,end=10) aa,aa,aa,aa,aa,aa,aa,aa,aa,aa,aa
c     go to 20
c10   continue
c     backspace(99)
 30   write(99,100)
     %   time,q1(ir,jr,kr),q2sur,q3(ir,jr,kr),
     %   pr(ir,jr,kr),delq3,voraz,vorr,vorz
      close(99)
        nul=nul+1
      end do
 100  format(2x,e14.7,2x,10(e12.5,2x))
      pi = 2.*asin(1.)
      return
      end
c***********************************************************************
c                                                                      *
c   ************** INDSON ************************                     *
c                                                                      *
c***********************************************************************
      subroutine indson
c
c     This routnine finds the indices in the computational grid
c     close to the physical location of the probe.
c     This is necessary since the location of the probe is given
c     in terms of physical coordinates (sonda.in)
c
      include 'param.f'
c
        pi = 2. * asin(1.)
      do ll = 1,nson
         ierr=0
         if((coson(ll,2).lt.0.).or.(coson(ll,2).gt.r0)) ierr=1
         if((coson(ll,3).lt.0.).or.(coson(ll,3).gt.alx3)) ierr=1
      if(ierr.eq.1) then     !{ this finishes at the end of the routine
         write(6,*) 'WARNING: la sonda ',ll,' e'' fuori dal dominio '
         write(32,*) 'WARNING: la sonda ',ll,' e'' fuori dal dominio '
      else
c
c     SEARCHING FOR THE INDICES IN THE GRID CLOSE TO THE PROBE
c
         do j=1,n2m
            jp = j+1
            if((rc(j).le.coson(ll,2)).and.(rc(jp).gt.coson(ll,2))) then
               mason(ll,2)=j
            end if
         end do
         do k=1,n3m
            kp = k+1
            if((zz(k).le.coson(ll,3)).and.(zz(kp).gt.coson(ll,3))) then
               mason(ll,3)=k
            end if
         end do
CC       cosnn = coson(ll,1)/ (2. * pi)
CC       ncosnn = cosnn
c
c        THIS COORDINATE IS IN BETWEEN 0 AND 2*PI (just in case!)
c
CC       cosnn = cosnn - ncosnn
         do i=1,n1m
         if(n1m.gt.1) then
            ip = ipv(i)
            if(coson(ll,1).gt.1.e-02) then
               if((thetac(i).le.coson(ll,1)).and.
     %            (thetac(ip).gt.coson(ll,1))) then
                  mason(ll,1)=i
               end if
            else
               mason(ll,1)=1
            end if
         else
            mason(ll,1)=1
         end if
         end do
      end if    !  here }
      end do
c
      return
      end
c
c***********************************************************************
c                                                                      *
c   ******* ENERC : compute energy spectrum (in a cylinder) *******    *
c                                                                      *
c***********************************************************************
      subroutine enerc(time)
      include 'param.f'
      dimension ener(m1),rtr(m2),ffat(m1)
      real   xxr(m1m+2,m3m),wwork(m1m+1,m3m)
      character*14 filth,filte,filba,filen,filso,filve,filnu,filbaper
      common/names/filth,filte,filba,filen,filso,filve,filnu,filbaper
c
c***********************************************************************
c
      call fftv
      dvol = 1./(dx2 * dx3)
      pi = 2. * asin (1.)
      n1mh=n1m/2+1
      n1md=n1m+2
      ffat(1) = pi
      do i = 2,n1mh
        ffat(i) = 2.*pi
      end do
      rtr(1)=1.
      do 111 j=2,n2
        rtr(j)=1./rc(j)
 111  continue
      do 1 j=1,n2
        do 11 k=1,n3m
          xxr(1,k)=q2(n1m,j,k)*rtr(j)
          do 21 i=1,n1m
            is=i+1
            xxr(is,k)=q2(i,j,k)*rtr(j)
   21     continue
          xxr(n1md,k)=q2(1,j,k)*rtr(j)
   11   continue
c
c     fft applied to the (q_2) by fft99
c   from physical to wave number space
c
        call fft99(xxr,wwork,trigxv,ifxv,1,m1+1,n1m,n3m,-1)
        do 31 i=1,n1mh
          ip=2*i
          id=2*i-1
          do 41 k=1,n3m
            iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
            rhs(iadd )=real( cmplx( xxr(id,k),xxr(ip,k))*
     1          conjg(cmplx( xxr(id,k),xxr(ip,k)) )  )
  41      continue
  31    continue
   1  continue
c
c
      do 2 j=1,n2m
        do 12 k=1,n3m
          xxr(1,k)=q1(n1m,j,k)
          do 22 i=1,n1m
            is=i+1
            xxr(is,k)=q1(i,j,k)
   22     continue
          xxr(n1md,k)=q1(1,j,k)
   12   continue
c
c     fft applied to the (q_1) by fft99
c   from physical to wave number space
c
        call fft99(xxr,wwork,trigxv,ifxv,1,m1+1,n1m,n3m,-1)
        do 32 i=1,n1mh
          ip=2*i
          id=2*i-1
          do 42 k=1,n3m
            iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
            rhs(iadd )=rhs(iadd )+real( cmplx( xxr(id,k),xxr(ip,k))*
     1                conjg(cmplx( xxr(id,k),xxr(ip,k)) )  )
  42        continue
  32    continue
   2  continue
c
c
      do 3 j=1,n2m
        do 13 k=1,n3m
          xxr(1,k)=q3(n1m,j,k)
          do 23 i=1,n1m
            is=i+1
            xxr(is,k)=q3(i,j,k)
   23     continue
          xxr(n1md,k)=q3(1,j,k)
   13   continue
c
c     fft applied to the (q_3) by fft99
c   from physical to wave number space
c
        call fft99(xxr,wwork,trigxv,ifxv,1,m1+1,n1m,n3m,-1)
        do 33 i=1,n1mh
          ip=2*i
          id=2*i-1
          do 43 k=1,n3m
            iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
            rhs(iadd )=rhs(iadd )+real( cmplx( xxr(id,k),xxr(ip,k))*
     1                conjg(cmplx( xxr(id,k),xxr(ip,k)) )  )
  43      continue
  33    continue
   3  continue
c
c
c***********************************************************************
c***********************************************************************
      do 116 i=1,n1mh
        ener(i)=0.
        do 116 k=1,n3m
          do 116 j=1,n2m
            iadd=i+(j-1)*n1m+(k-1)*n1m*n2m
            ener(i)=ener(i)
     %             +ffat(i)*rhs(iadd)*rm(j)*g2rm(j)*g3rm(k)*dvol
 116  continue
      open(39,file=filen,status='unknown')
      write(39,546) time,(ener(i),i=nini,nfin,nstri)
 546  format(16(1x,e14.6))
      close(39)
      return
      end
c***********************************************************************
      subroutine cfl(cflm)
      include'param.f'
c
      cflm=0.
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m)
!$OMP$ SHARED(q2,q3,ipv,udx2rm,udx3m)
!$OMP$ PRIVATE(kp,jp,ip,qcf,i,j,k)
!$OMP$ REDUCTION(max:cflm)
      do 7 k=1,n3m
        kp=k+1
        do 7 j=1,n2m
          jp=j+1
          do 7 i=1,n1m
            ip=ipv(i)
            qcf=((
     %           +abs((q2(i,j,k)+q2(i,jp,k))*0.5*udx2rm(j)))
     %           +abs((q3(i,j,k)+q3(i,j,kp))*0.5*udx3m(k)))
            cflm = max(cflm,qcf)
  7   continue
!$OMP  END PARALLEL DO
      return
      end
c
c***********************************************************************
      subroutine outpf(time)
      include'param.f'
      character*70 namfi3
      character*5 ipfi
      dimension vx(m1,m2,m3),vy(m1,m2,m3),vz(m1,m2,m3)
      dimension prc(m1,m2,m3)
      dimension temp(m1,m2,m3)
cccccccccccccccccccccc
c     go to 7766
c     Spherical components

      do kc=2,n3m
c       km=kmv(kc)
        km=kc-1
c       kp=kpv(kc)
        do jc=2,n2
c         jm=jmv(jc)
          jm=jc-1
c         jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
c           ip=ipv(ic)
c           ksy = ksym(kc)
c           isy = isym(ic)
c Cartesian coordinates
c           xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
c           yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
c           zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
            vt = 0.25*(q1(ic,jc,kc)+q1(ic,jm,kc)
     %                +q1(ic,jc,km)+q1(ic,jm,km))
            vr = 0.25*(q2(ic,jc,kc)+q2(im,jc,kc)
     %                 +q2(ic,jc,km)+q2(im,jc,km))/rc(jc)**2
            vp = 0.25*(q3(ic,jc,kc)+q3(im,jc,kc)
     %                +q3(ic,jm,kc)+q3(im,jm,kc))/sin(zz(kc))
            vx(ic,jc,kc) = vt
            vy(ic,jc,kc) = vr
            vz(ic,jc,kc) = vp
            prc(ic,jc,kc)= 0.125*(pr(ic,jc,kc)+pr(ic,jm,kc)
     %                           +pr(ic,jc,km)+pr(ic,jm,km) 
     %                           +pr(im,jc,kc)+pr(im,jm,kc)
     %                           +pr(im,jc,km)+pr(im,jm,km) )
            temp(ic,jc,kc)= 0.125*(dens(ic,jc,kc)+dens(ic,jm,kc)
     %                           +dens(ic,jc,km)+dens(ic,jm,km) 
     %                           +dens(im,jc,kc)+dens(im,jm,kc)
     %                           +dens(im,jc,km)+dens(im,jm,km) )
          end do
        end do
      end do
c North pole
      kc=1
      kp=kpv(kc)
      do jc=2,n2
c       jm=jmv(jc)
        jm=jc-1
c       jp=jpv(jc)
        do ic=1,n1m
c         im=imv(ic)
          ip=ipv(ic)
c         ksy = ksym(kc)
c         ksc = ksyc(kc)
          isy = isym(ic)
c         isyx = isym(im)
c Cartesian coordinates
c         xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
c         yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
c         zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vt = 0.25*(q1(ic,jc,kc)-q1(isy,jc,kc)
     %              + q1(ic,jm,kc)-q1(isy,jm,kc))
          vr = 0.5*(q2(ic,jc,kc)+q2(isy,jc,kc))/rc(jc)**2
c         vr = 0.25*(q2(ic,jc,kc)+q2(isy,jc,kc)
c    %          +q2(im,jc,kc)+q2(isyx,jc,kc) )/rc(jc)**2
          vp = 0.25*(q3(ic,jc,kp)-q3(isy,jc,kp)
     %              +q3(ic,jm,kp)-q3(isy,jm,kp))/sin(zz(kp))
c         vp = 0.125*(q3(ic,jc,kp)-q3(isy,jc,kp)
c    %     +q3(im,jc,kp)-q3(isyx,jc,kp)
c    %     +q3(ic,jm,kp)-q3(isy,jm,kp)
c    %     +q3(im,jm,kp)-q3(isyx,jm,kp)  )/sin(zz(kp))
          vx(ic,jc,kc) = vt
          vy(ic,jc,kc) = vr
          vz(ic,jc,kc) = vp
          prc(ic,jc,kc) = 0.25*(pr(ic,jc,kc)+pr(isy,jc,kc)
     %                         +pr(ic,jm,kc)+pr(isy,jm,kc))
          temp(ic,jc,kc) = 0.25*(dens(ic,jc,kc)+dens(isy,jc,kc)
     %                         +dens(ic,jm,kc)+dens(isy,jm,kc))
        end do
      end do
c South pole
      kc=n3
      km=n3m
      do jc=2,n2
c       jm=jmv(jc)
        jm=jc-1
c       jp=jpv(jc)
        do ic=1,n1m
c         im=imv(ic)
c         ip=ipv(ic)
c         ksy = ksym(kc)
c         ksc = ksyc(kc)
          isy = isym(ic)
c         isyx = isym(im)
c Cartesian coordinates
c         xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
c         yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
c         zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vt = 0.25*(q1(ic,jc,kc)-q1(isy,jc,kc)
     %              + q1(ic,jm,kc)-q1(isy,jm,kc))
          vr = 0.5*(q2(ic,jc,km)+q2(isy,jc,km))/rc(jc)**2
c         vr = 0.25*(q2(ic,jc,km)+q2(isy,jc,km)
c    %          +q2(im,jc,km)+q2(isyx,jc,km) )/rc(jc)**2
          vp = 0.25*(q3(ic,jc,km)-q3(isy,jc,km)
     %              +q3(ic,jm,km)-q3(isy,jm,km))/sin(zz(km))
c         vp = 0.125*(q3(ic,jc,km)-q3(isy,jc,km)
c    %     +q3(im,jc,km)-q3(isyx,jc,km)
c    %     +q3(ic,jm,km)-q3(isy,jm,km)
c    %     +q3(im,jm,km)-q3(isyx,jm,km)  )/sin(zz(km))
          vx(ic,jc,kc) = vt
          vy(ic,jc,kc) = vr
          vz(ic,jc,kc) = vp
          prc(ic,jc,kc) = 0.25*(pr(ic,jc,km)+pr(isy,jc,km)
     %                         +pr(ic,jm,km)+pr(isy,jm,km))
          temp(ic,jc,kc) = 0.25*(dens(ic,jc,km)+dens(isy,jc,km)
     %                         +dens(ic,jm,km)+dens(isy,jm,km))
        end do
      end do
c variables at the sphere centre
      do kc=2,n3m
c       km=kmv(kc)
        km=kc-1
c       kp=kpv(kc)
        jc=1
        jm=jmv(jc)
        jp=jpv(jc)
        do ic=1,n1m
          im=imv(ic)
c         ip=ipv(ic)
          ksy = ksym(kc)
          ksyy = ksym(km)
          ksyx = ksyc(kc)
c         ksymc = ksyc(km)
          isy = isym(ic)
          isyx=isym(im)
c         xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
c         yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
c         zzz = rc(jc)*cos(zz(kc))
          vt = 0.25*(q1(ic,jc,kc)-q1(isy,jc,ksy)
     %              +q1(ic,jc,km)-q1(isy,jc,ksyy))
          vr = 0.125*(q2(ic,jp,kc)-q2(isy,jp,ksy)
     %               +q2(im,jp,kc)-q2(isyx,jp,ksy)
     %               +q2(ic,jp,km)-q2(isy,jp,ksyy)
     %               +q2(im,jp,km)-q2(isyx,jp,ksyy))/rc(jp)**2
          vp = 0.25*(q3(ic,jc,kc)+q3(isy,jc,ksyx)
     %              +q3(im,jc,kc)+q3(isyx,jc,ksyx))/sin(zz(kc))
c Cartesian velocities
          vx(ic,jc,kc) = vt
          vy(ic,jc,kc) = vr
          vz(ic,jc,kc) = vp
          if(inslwri.eq.0) then
            vx(ic,1,1) = 0.
            vy(ic,1,1) = 0.
            vz(ic,1,1) = 0.
          end if
          prc(ic,jc,kc) = 0.25*(pr(ic,jc,kc)+pr(isy,jc,ksy)
     %                         +pr(ic,jm,kc)+pr(isy,jm,ksy))
          temp(ic,jc,kc) = 1. !0.25*(dens(ic,jc,kc)+dens(isy,jc,ksy)
c    %                         +dens(ic,jm,kc)+dens(isy,jm,ksy))
        end do
      end do
      do ic=1,n1m
        if(inslwri.eq.1) then
          vx(ic,1,1) = vx(ic,1,2)
          vy(ic,1,1) = vy(ic,1,2)
          vz(ic,1,1) = vz(ic,1,2)
        else
          vx(ic,1,1) = 0.
          vy(ic,1,1) = 0.
          vz(ic,1,1) = 0.
        end if
        if(inslwr.eq.1) then
          vx(ic,n2,1) = vx(ic,n2m,2)
          vy(ic,n2,1) = vy(ic,n2m,2)
          vz(ic,n2,1) = vz(ic,n2m,2)
        else
          vx(ic,1,1) = 0.
          vy(ic,1,1) = 0.
          vz(ic,1,1) = 0.
        end if
        prc(ic,1,1) = prc(ic,1,2)
        temp(ic,1,1) = temp(ic,1,2)
        vx(ic,1,n3) = vx(ic,1,n3m)
        vy(ic,1,n3) = vy(ic,1,n3m)
        vz(ic,1,n3) = vz(ic,1,n3m)
        prc(ic,1,n3) = prc(ic,1,n3m)
        temp(ic,1,n3) = temp(ic,1,n3m)
      end do
      do kc=1,n3
        do ic=1,n1
          pr(n1,jc,kc) = pr(1,jc,kc)
          if(inslwr.eq.0) then
            vx(n1,jc,kc) = 0.
            vy(n1,jc,kc) = 0.
            vz(n1,jc,kc) = 0.
          end if
          temp(ic,n2,kc) = 0.
        end do
      end do
      do kc=1,n3
        do jc=1,n2
          pr(n1,jc,kc) = pr(1,jc,kc)
          vx(n1,jc,kc) = vx(1,jc,kc)
          vy(n1,jc,kc) = vy(1,jc,kc)
          vz(n1,jc,kc) = vz(1,jc,kc)
          prc(n1,jc,kc) = prc(1,jc,kc)
          temp(n1,jc,kc) = temp(1,jc,kc)
        end do
      end do
      go to 7755
c
c   Omogeneizzazione all'asse
c

c     
      do j=1,n2
        vxas = 0.
        vyas = 0.
        vzas = 0.
        do i=1,n1
          vxas = vxas + vx(i,j,n3m)
          vyas = vyas + vy(i,j,n3m)
          vzas = vzas + vz(i,j,n3m)
        end do
        vxas = vxas/float(n1)
        vyas = vyas/float(n1)
        vzas = vzas/float(n1)
        do i=1,n1
          vx(i,j,n3) = vxas
          vy(i,j,n3) = vyas
          vz(i,j,n3) = vzas
        end do
      end do
c     
      do j=1,n2
        vxas = 0.
        vyas = 0.
        vzas = 0.
        do i=1,n1
          vxas = vxas + vx(i,j,2)
          vyas = vyas + vy(i,j,2)
          vzas = vzas + vz(i,j,2)
        end do
        vxas = vxas/float(n1)
        vyas = vyas/float(n1)
        vzas = vzas/float(n1)
        do i=1,n1
          vx(i,j,1) = vxas
          vy(i,j,1) = vyas
          vz(i,j,1) = vzas
        end do
      end do
c     
      vxce = 0.
      vyce = 0.
      vzce = 0.
      do k=1,n3
        do i=1,n1
          vxce = vxce + vx(i,1,k)
          vyce = vyce + vy(i,1,k)
          vzce = vzce + vz(i,1,k)
        end do
      end do
      vxce = vxce/float(n1*n3)
      vyce = vyce/float(n1*n3)
      vzce = vzce/float(n1*n3)
      do k=1,n3
        do i=1,n1
          vx(i,1,k) = vxce
          vy(i,1,k) = vyce
          vz(i,1,k) = vzce
        end do
      end do
      if(rc(1).le.1.e-3) then
        jc=1
        pmcen = 0.
        do kc=1,n3m
          do ic=1,n1m
            pmcen = pmcen + pr(ic,jc,kc)
          end do
        end do
        pmcen = pmcen/float(n1m*n3m)
        do kc=1,n3 
          do ic=1,n1 
            prc(ic,jc,kc) = pmcen
          end do
        end do
      end if
 7755 continue
      itime=10.*time !+0.5
      write(ipfi,99) itime
      namfi3='./dati/spbou'//ipfi//'.q'
c     namfi3='/scratch/verzicco/barun/dati/bou'//ipfi//'.dat'
      write(6,201)time,namfi3
      n1pp=(n1-1)/n1p!+1
      n2pp=(n2-1)/n2p!+1
      n3pp=(n3-1)/n3p!+1
      open(99,file=namfi3,form='unformatted',status='unknown')
      rewind(99)
      write(99) 1
      write(99) n1,n2,n3
      write(99) epsil,epsil,re,time
      write(99)
c    %   ((((1.),
     %   ((((sqrt(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vx(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vy(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vz(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   ((((temp(i,j,k)),
     %   i=1,n1),j=1,n2),k=1,n3)
      close(99)

************************************************************************

c7766 continue

c     Cartesian components
      do kc=2,n3m
c       km=kmv(kc)
        km=kc-1
c       kp=kpv(kc)
        do jc=2,n2
c         jm=jmv(jc)
          jm=jc-1
c         jp=jpv(jc)
          do ic=1,n1m
            im=imv(ic)
c           ip=ipv(ic)
c           ksy = ksym(kc)
c           isy = isym(ic)
c Cartesian coordinates
c           xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
c           yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
c           zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
            vt = 0.25*(q1(ic,jc,kc)+q1(ic,jm,kc)
     %                +q1(ic,jc,km)+q1(ic,jm,km))
            vr = 0.25*(q2(ic,jc,kc)+q2(im,jc,kc)
     %                 +q2(ic,jc,km)+q2(im,jc,km))/rc(jc)**2
            vp = 0.25*(q3(ic,jc,kc)+q3(im,jc,kc)
     %                +q3(ic,jm,kc)+q3(im,jm,kc))/sin(zz(kc))
            ct = cos(thetac(ic))
            st = sin(thetac(ic))
            cp = cos(zz(kc))
            sp = sin(zz(kc))
c           vx(ic,jc,kc) = vr*ct*sp-vt*st*sp+vp*ct*cp
c           vy(ic,jc,kc) = vr*st*sp+vt*ct*sp+vp*st*cp
c           vz(ic,jc,kc) = vr*cp-vp*sp
            vx(ic,jc,kc) = vr*ct*sp-vt*st+vp*ct*cp
            vy(ic,jc,kc) = vr*st*sp+vt*ct+vp*st*cp
            vz(ic,jc,kc) = vr*cp-vp*sp
            prc(ic,jc,kc)= 0.125*(pr(ic,jc,kc)+pr(ic,jm,kc)
     %                           +pr(ic,jc,km)+pr(ic,jm,km) 
     %                           +pr(im,jc,kc)+pr(im,jm,kc)
     %                           +pr(im,jc,km)+pr(im,jm,km) )
            temp(ic,jc,kc)= 0.125*(dens(ic,jc,kc)+dens(ic,jm,kc)
     %                           +dens(ic,jc,km)+dens(ic,jm,km) 
     %                           +dens(im,jc,kc)+dens(im,jm,kc)
     %                           +dens(im,jc,km)+dens(im,jm,km) )
          end do
        end do
      end do
c variables at the sphere centre
      do kc=2,n3m
c       km=kmv(kc)
        km=kc-1
c       kp=kpv(kc)
        jc=1
        jm=jmv(jc)
        jp=jpv(jc)
        do ic=1,n1m
          im=imv(ic)
c         ip=ipv(ic)
          ksy = ksym(kc)
          ksyy = ksym(km)
          ksyx = ksyc(kc)
c         ksymc = ksyc(km)
          isy = isym(ic)
          isyx=isym(im)
c         xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
c         yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
c         zzz = rc(jc)*cos(zz(kc))
          vt = 0.25*(q1(ic,jc,kc)-q1(isy,jc,ksy)
     %              +q1(ic,jc,km)-q1(isy,jc,ksyy))
          vr = 0.125*(q2(ic,jp,kc)-q2(isy,jp,ksy)
     %               +q2(im,jp,kc)-q2(isyx,jp,ksy)
     %               +q2(ic,jp,km)-q2(isy,jp,ksyy)
     %               +q2(im,jp,km)-q2(isyx,jp,ksyy))/rc(jp)**2
          vp = 0.25*(q3(ic,jc,kc)+q3(isy,jc,ksyx)
     %              +q3(im,jc,kc)+q3(isyx,jc,ksyx))/sin(zz(kc))
          ct = cos(thetac(ic))
          st = sin(thetac(ic))
          cp = cos(zz(kc))
          sp = sin(zz(kc))
c Cartesian velocities
c         vx(ic,jc,kc) = vr*ct*sp-vt*st*sp+vp*ct*cp
c         vy(ic,jc,kc) = vr*st*sp+vt*ct*sp+vp*st*cp
c         vz(ic,jc,kc) = vr*cp-vp*sp
          vx(ic,jc,kc) = vr*ct*sp-vt*st+vp*ct*cp
          vy(ic,jc,kc) = vr*st*sp+vt*ct+vp*st*cp
          vz(ic,jc,kc) = vr*cp-vp*sp
          prc(ic,jc,kc) = 0.25*(pr(ic,jc,kc)+pr(isy,jc,ksy)
     %                         +pr(ic,jm,kc)+pr(isy,jm,ksy))
          temp(ic,jc,kc) = 0.25*(dens(ic,jc,kc)+dens(isy,jc,ksy)
     %                         +dens(ic,jm,kc)+dens(isy,jm,ksy))
        end do
      end do
c North pole
      kc=1
      kp=kpv(kc)
      do jc=2,n2
c       jm=jmv(jc)
        jm=jc-1
c       jp=jpv(jc)
        do ic=1,n1m
          im=imv(ic)
c         ip=ipv(ic)
c         ksy = ksym(kc)
c         ksc = ksyc(kc)
          isy = isym(ic)
          isyx = isym(im)
c Cartesian coordinates
c         xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
c         yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
c         zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vt = 0.25*(q1(ic,jc,kc)-q1(isy,jc,kc)
     %              + q1(ic,jm,kc)-q1(isy,jm,kc))
          vr = 0.25*(q2(ic,jc,kc)+q2(isy,jc,kc)
     %          +q2(im,jc,kc)+q2(isyx,jc,kc) )/rc(jc)**2
          vp = 0.125*(q3(ic,jc,kp)-q3(isy,jc,kp)
     %       +q3(im,jc,kp)-q3(isyx,jc,kp)
     %       +q3(ic,jm,kp)-q3(isy,jm,kp)
     %       +q3(im,jm,kp)-q3(isyx,jm,kp) )/sin(zz(kp))
          ct = cos(thetac(ic))
          st = sin(thetac(ic))
          cp = cos(zz(kc))
          sp = sin(zz(kc))
c         vx(ic,jc,kc) = vr*ct*sp-vt*st*sp+vp*ct*cp
c         vy(ic,jc,kc) = vr*st*sp+vt*ct*sp+vp*st*cp
c         vz(ic,jc,kc) = vr*cp-vp*sp
          vx(ic,jc,kc) = vr*ct*sp-vt*st+vp*ct*cp
          vy(ic,jc,kc) = vr*st*sp+vt*ct+vp*st*cp
          vz(ic,jc,kc) = vr*cp-vp*sp
          prc(ic,jc,kc) = 0.25*(pr(ic,jc,kc)+pr(isy,jc,kc)
     %                         +pr(ic,jm,kc)+pr(isy,jm,kc))
          temp(ic,jc,kc) = 0.25*(dens(ic,jc,kc)+dens(isy,jc,kc)
     %                         +dens(ic,jm,kc)+dens(isy,jm,kc))
        end do
      end do
c South pole
      kc=n3
      km=n3m
      do jc=2,n2
c       jm=jmv(jc)
c       if(jc.eq.2) jm=2  !because for j=1 we haven't defined the variables jet
        jm = jc-1
c       jp=jpv(jc)
        do ic=1,n1m
          im=imv(ic)
c         ip=ipv(ic)
c         ksy = ksym(kc)
c         ksc = ksyc(kc)
          isy = isym(ic)
          isyx = isym(im)
c Cartesian coordinates
c         xxx = rc(jc)*cos(thetac(ic))*sin(zz(kc))
c         yyy = rc(jc)*sin(thetac(ic))*sin(zz(kc))
c         zzz = rc(jc)*cos(zz(kc))
c velocities at the cell nodes
          vt = 0.25*(q1(ic,jc,kc)-q1(isy,jc,kc)
     %              + q1(ic,jm,kc)-q1(isy,jm,kc))
          vr = 0.25*(q2(ic,jc,km)+q2(isy,jc,km)
     %          +q2(im,jc,km)+q2(isyx,jc,km) )/rc(jc)**2

          vp = 0.125*(q3(ic,jc,km)-q3(isy,jc,km)
     %     +q3(im,jc,km)-q3(isyx,jc,km)
     %     +q3(ic,jm,km)-q3(isy,jm,km)
     %     +q3(im,jm,km)-q3(isyx,jm,km)  )/sin(zz(km))

          ct = cos(thetac(ic))
          st = sin(thetac(ic))
          cp = cos(zz(kc))
          sp = sin(zz(kc))
c         vx(ic,jc,kc) = vr*ct*sp-vt*st*sp+vp*ct*cp
c         vy(ic,jc,kc) = vr*st*sp+vt*ct*sp+vp*st*cp
c         vz(ic,jc,kc) = vr*cp-vp*sp
          vx(ic,jc,kc) = vr*ct*sp-vt*st+vp*ct*cp
          vy(ic,jc,kc) = vr*st*sp+vt*ct+vp*st*cp
          vz(ic,jc,kc) = vr*cp-vp*sp
          prc(ic,jc,kc) = 0.25*(pr(ic,jc,km)+pr(isy,jc,km)
     %                         +pr(ic,jm,km)+pr(isy,jm,km))
          temp(ic,jc,kc) = 0.25*(dens(ic,jc,km)+dens(isy,jc,km)
     %                         +dens(ic,jm,km)+dens(isy,jm,km))
        end do
      end do
      do ic=1,n1
        vx(ic,1,1) = vx(ic,1,2)
        vy(ic,1,1) = vy(ic,1,2)
        vz(ic,1,1) = vz(ic,1,2)
        vx(ic,1,n3) = vx(ic,1,n3m)
        vy(ic,1,n3) = vy(ic,1,n3m)
        vz(ic,1,n3) = vz(ic,1,n3m)
        prc(ic,1,n3) = prc(ic,1,n3m)
        temp(ic,1,n3) = temp(ic,1,n3m)
      end do
      do kc=1,n3
        do jc=1,n2
          prc(n1,jc,kc) = prc(1,jc,kc)
          vx(n1,jc,kc) = vx(1,jc,kc)
          vy(n1,jc,kc) = vy(1,jc,kc)
          vz(n1,jc,kc) = vz(1,jc,kc)
          pr(n1,jc,kc) = pr(1,jc,kc)
          temp(n1,jc,kc) = temp(1,jc,kc)
        end do
      end do
c
c   Omogeneizzazione all'asse
c

c     
      do j=1,n2
        vxas = 0.
        vyas = 0.
        vzas = 0.
        do i=1,n1
          vxas = vxas + vx(i,j,n3m)
          vyas = vyas + vy(i,j,n3m)
          vzas = vzas + vz(i,j,n3m)
        end do
        vxas = vxas/float(n1)
        vyas = vyas/float(n1)
        vzas = vzas/float(n1)
        do i=1,n1
          vx(i,j,n3) = vxas
          vy(i,j,n3) = vyas
          vz(i,j,n3) = vzas
        end do
      end do
      
      do j=1,n2
        vxas = 0.
        vyas = 0.
        vzas = 0.
        do i=1,n1
          vxas = vxas + vx(i,j,2)
          vyas = vyas + vy(i,j,2)
          vzas = vzas + vz(i,j,2)
        end do
        vxas = vxas/float(n1)
        vyas = vyas/float(n1)
        vzas = vzas/float(n1)
        do i=1,n1
          vx(i,j,1) = vxas
          vy(i,j,1) = vyas
          vz(i,j,1) = vzas
        end do
      end do
c     
      vxce = 0.
      vyce = 0.
      vzce = 0.
      do k=1,n3
        do i=1,n1
          vxce = vxce + vx(i,1,k)
          vyce = vyce + vy(i,1,k)
          vzce = vzce + vz(i,1,k)
        end do
      end do
      vxce = vxce/float(n1*n3)
      vyce = vyce/float(n1*n3)
      vzce = vzce/float(n1*n3)
      do k=1,n3
        do i=1,n1
          vx(i,1,k) = vxce
          vy(i,1,k) = vyce
          vz(i,1,k) = vzce
        end do
      end do
      if(rc(1).le.1.e-3) then
        jc=1
        pmcen = 0.
        do kc=1,n3m
          do ic=1,n1m
            pmcen = pmcen + pr(ic,jc,kc)
          end do
        end do
        pmcen = pmcen/float(n1m*n3m)
        do kc=1,n3 
          do ic=1,n1 
            prc(ic,jc,kc) = pmcen
          end do
        end do
      end if
cccccccccccccccccccccc
c     do k=1,n3m
c       do j=1,n2m
c         do i=1,n1m
c           ip = ipv(i)
c           ct= cos(thetam(i))
c           st= sin(thetam(i))
c           cp= cos(zm(k))
c           sp= sin(zm(k))
c           vx(i,j,k) = 0.5*(q2(i,j,k)+q2(i,j+1,k))*ct*sp/rm(j)**2-
c    %                  0.5*(q1(i,j,k)+q1(ip,j,k))*st*sp*rm(j)+
c    %                  0.5*(q3(i,j,k)+q3(i,j,k+1))*ct*cp*rm(j)/sp
c           vy(i,j,k) = 0.5*(q2(i,j,k)+q2(i,j+1,k))*st*sp/rm(j)**2+
c    %                  0.5*(q1(i,j,k)+q1(ip,j,k))*ct*sp*rm(j)+
c    %                  0.5*(q3(i,j,k)+q3(i,j,k+1))*st*cp*rm(j)/sp
c           vz(i,j,k) = 0.5*(q2(i,j,k)+q2(i,j+1,k))*cp/rm(j)**2-
c    %                  0.5*(q3(i,j,k)+q3(i,j,k+1))*sp*rm(j)/sp
c         end do
c       end do
c     end do
      itime=10.*time !+0.5
      write(ipfi,99) itime
   99 format(i5.5)
      namfi3='./dati/boum'//ipfi//'.q'
c     namfi3='/scratch/verzicco/barun/dati/bou'//ipfi//'.dat'
      write(6,201)time,namfi3
  201 format(10x,'At t=',e10.3,' write data on ',a20)
      n1pp=(n1-1)/n1p!+1
      n2pp=(n2-1)/n2p!+1
      n3pp=(n3-1)/n3p!+1
      open(99,file=namfi3,form='unformatted',status='unknown')
      rewind(99)
      write(99) 1
      write(99) n1,n2,n3
      write(99) epsil,epsil,re,time
      write(99)
c    %   ((((sqrt(vx(i,j,k)**2+vy(i,j,k)**2+vz(i,j,k)**2)),
     %   ((((1.),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vx(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vy(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   (((((vz(i,j,k))),
     %   i=1,n1),j=1,n2),k=1,n3),
c    %   ((((0.5*(q1(i,j,k)+q1(i+1,j,k))),
c    %   i=1,n1m),j=1,n2m),k=1,n3m),
c    %   ((((0.5*(q2(i,j,k)+q2(i,j+1,k))),
c    %   i=1,n1m),j=1,n2m),k=1,n3m),
c    %   ((((0.5*(q3(i,j,k)+q3(i,j,k+1))),
c    %   i=1,n1m),j=1,n2m),k=1,n3m),
     %   ((((temp(i,j,k)),
!     %   ((((pr(i,j,k)),
     %   i=1,n1),j=1,n2),k=1,n3)
      close(99)

      return
      end
c
c***********************************************************************
      subroutine diss(time)
      include'param.f'
c
      pi=2.*asin(1.)
      vol = 4.*pi/3.*(rc(n2)**3-rc(1)**3)
c     vol = 0.
c     do kc=1,n3m
c      do jc=1,n2m
c       do ic=1,n1m
c        vol = vol + 
c    %    g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
c       end do
c      end do
c     end do
c
      disski = 0.
c     go to 4343
c
c   ! vorticita' longitudinale !   
      do kc=2,n3m-1
        km=kmv(kc)           
        kp=kpv(kc)           
        kpp=kc+1
        do jc=2,n2m-1     
          jm=jmv(jc)       
          jp=jpv(jc)       
          jpp=jc+1
          udx3=dx3/(g3rc(kc))
          udx2=dx2/g2rc(jc)
          udx3p=dx3/(g3rc(kp))
          udx3pp=dx3/(g3rc(kpp))
          udx2p=dx2/g2rc(jpp)
            do ic=1,n1m 
              dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*udx3/rc(jc)**2
              dq3x2=(rm(jc)*q3(ic,jc,kc)-rm(jm)*q3(ic,jm,kc))*udx2    
     %                /sinpc(kc)
              voraz=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kc)-q2(ic,jpp,km))*udx3/rc(jpp)**2
              dq3x2=(rm(jp)*q3(ic,jp,kc)-rm(jc)*q3(ic,jc,kc))*udx2p   
     %                /sinpc(kc)
              vorazj=(dq3x2-dq2x3)/rc(jpp)
              dq2x3=(q2(ic,jc,kp)-q2(ic,jc,kc))*udx3p/rc(jc)**2
              dq3x2=(rm(jc)*q3(ic,jc,kpp)-rm(jm)*q3(ic,jm,kpp))*udx2    
     %                /sinpc(kpp)
              vorazk=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kp)-q2(ic,jpp,kc))*udx3pp/rc(jpp)**2
              dq3x2=(rm(jp)*q3(ic,jp,kpp)-rm(jc)*q3(ic,jc,kpp))*udx2p   
     %                /sinpc(kpp)
              vorazjk=(dq3x2-dq2x3)/rc(jpp)
              om1=0.25*(voraz+vorazj+vorazk+vorazjk)
c             write(6,*) ic,jc,kc,voraz,vorazj,vorazk,vorazjk
              disski = disski +om1**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
            end do
          end do
c Inner shell
          jc=1
          jm=jmv(jc)       
          jp=jpv(jc)       
          jpp=jc+1
          udx3=dx3/(g3rc(kc))
          udx2=dx2/g2rc(jc)
          udx3p=dx3/(g3rc(kp))
          udx3pp=dx3/(g3rc(kpp))
          udx2p=dx2/g2rc(jpp)
            do ic=1,n1m 
              dq2x3=0.
              dq3x2=(rm(jc)*q3(ic,jc,kc)-0.)*2.*udx2/sinpc(kc)
              voraz=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kc)-q2(ic,jpp,km))*udx3/rc(jpp)**2
              dq3x2=(rm(jp)*q3(ic,jp,kc)-rm(jc)*q3(ic,jc,kc))*udx2p   
     %                /sinpc(kc)
              vorazj=(dq3x2-dq2x3)/rc(jpp)
              dq2x3=0.
              dq3x2=(rm(jc)*q3(ic,jc,kpp)-0.)*2.*udx2/sinpc(kpp)
              vorazk=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kp)-q2(ic,jpp,kc))*udx3pp/rc(jpp)**2
              dq3x2=(rm(jp)*q3(ic,jp,kpp)-rm(jc)*q3(ic,jc,kpp))*udx2p   
     %                /sinpc(kpp)
              vorazjk=(dq3x2-dq2x3)/rc(jpp)
              om1=0.25*(voraz+vorazj+vorazk+vorazjk)
c             write(6,*) ic,jc,kc,voraz,vorazj,vorazk,vorazjk
              disski = disski +om1**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
            end do
c Outer shell
          jc=n2m
          jm=jmv(jc)       
          jp=jpv(jc)       
          jpp=jc+1
          udx3=dx3/(g3rc(kc))
          udx2=dx2/g2rc(jc)
          udx3p=dx3/(g3rc(kp))
          udx3pp=dx3/(g3rc(kpp))
          udx2p=dx2/g2rc(jpp)
            do ic=1,n1m 
              dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*udx3/rc(jc)**2
              dq3x2=(rm(jc)*q3(ic,jc,kc)-rm(jm)*q3(ic,jm,kc))*udx2    
     %                /sinpc(kc)
              voraz=(dq3x2-dq2x3)/rc(jc)
              dq2x3=0.
              dq3x2=(0.-rm(jc)*q3(ic,jc,kc))*2.*udx2p/sinpc(kc)
              vorazj=(dq3x2-dq2x3)/rc(jpp)
              dq2x3=(q2(ic,jc,kp)-q2(ic,jc,kc))*udx3/rc(jc)**2
              dq3x2=(rm(jc)*q3(ic,jc,kpp)-rm(jm)*q3(ic,jm,kpp))*udx2    
     %                /sinpc(kpp)
              vorazk=(dq3x2-dq2x3)/rc(jc)
              dq2x3=0.
              dq3x2=(0.-rm(jc)*q3(ic,jc,kpp))*2.*udx2p/sinpc(kpp)
              vorazjk=(dq3x2-dq2x3)/rc(jpp)
              om1=0.25*(voraz+vorazj+vorazk+vorazjk)
c             write(6,*) ic,jc,kc,voraz,vorazj,vorazk,vorazjk
              disski = disski +om1**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
            end do
c
        end do
c  
c  North pole
        kc=1
        km=kmv(kc)           
        kp=kpv(kc)           
        kpp=kc+1
        do jc=2,n2m-1     
          jm=jmv(jc)       
          jp=jpv(jc)       
          jpp=jc+1
          udx3=dx3/(g3rc(kc))
          udx2=dx2/g2rc(jc)
          udx3p=dx3/(g3rc(kp))
          udx3pp=dx3/(g3rc(kpp))
          udx2p=dx2/g2rc(jpp)
            do ic=1,n1m 
              isy=isym(ic)
              dq2x3=(q2(ic,jc,kc)-q2(isy,jc,kc))*udx3/rc(jc)**2
              dq3x2=(rm(jc)*(q3(ic,jc,kp)-q3(isy,jc,kp))*0.5
     %              -rm(jm)*(q3(ic,jm,kp)-q3(isy,jm,kp))*0.5)*udx2    
     %                /sinpc(kp)
              voraz=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kc)-q2(isy,jpp,kc))*udx3p/rc(jpp)**2
              dq3x2=(rm(jp)*(q3(ic,jp,kp)-q3(isy,jp,kp))*0.5
     %              -rm(jc)*(q3(ic,jc,kp)-q3(isy,jc,kp))*0.5)*udx2p   
     %                /sinpc(kp)
              vorazj=(dq3x2-dq2x3)/rc(jpp)
              dq2x3=(q2(ic,jc,kp)-q2(ic,jc,kc))*udx3/rc(jc)**2
              dq3x2=(rm(jc)*q3(ic,jc,kpp)-rm(jm)*q3(ic,jm,kpp))*udx2    
     %                /sinpc(kp)
              vorazk=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kp)-q2(ic,jpp,kc))*udx3pp/rc(jpp)**2
              dq3x2=(rm(jp)*q3(ic,jp,kpp)-rm(jc)*q3(ic,jc,kpp))*udx2p   
     %                /sinpc(kp)
              vorazjk=(dq3x2-dq2x3)/rc(jpp)
              om1=0.25*(voraz+vorazj+vorazk+vorazjk)
c
              disski = disski +om1**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
            end do
          end do
c  Inner shell
          jc=1     
          jm=jmv(jc)       
          jp=jpv(jc)       
          jpp=jc+1
          udx3=dx3/(g3rc(kc))
          udx2=dx2/g2rc(jc)
          udx3p=dx3/(g3rc(kp))
          udx3pp=dx3/(g3rc(kpp))
          udx2p=dx2/g2rc(jpp)
            do ic=1,n1m 
              isy=isym(ic)
              dq2x3=0.
              dq3x2=(rm(jc)*(q3(ic,jc,kp)-q3(isy,jc,kp))*0.5
     %              -0.                                     )*2.*udx2
     %                /sinpc(kp)
              voraz=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kc)-q2(isy,jpp,kc))*udx3/rc(jpp)**2
              dq3x2=(rm(jp)*(q3(ic,jp,kp)-q3(isy,jp,kp))*0.5
     %              -rm(jc)*(q3(ic,jc,kp)-q3(isy,jc,kp))*0.5)*udx2p   
     %                /sinpc(kp)
              vorazj=(dq3x2-dq2x3)/rc(jpp)
              dq2x3=0.
              dq3x2=(rm(jc)*q3(ic,jc,kpp)-0.)*2.*udx2    
     %                /sinpc(kp)
              vorazk=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kp)-q2(ic,jpp,kc))*udx3pp/rc(jpp)**2
              dq3x2=(rm(jp)*q3(ic,jp,kpp)-rm(jc)*q3(ic,jc,kpp))*udx2p   
     %                /sinpc(kp)
              vorazjk=(dq3x2-dq2x3)/rc(jpp)
              om1=0.25*(voraz+vorazj+vorazk+vorazjk)
c
              disski = disski +om1**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
            end do
c Outer shell
          jc=n2m     
          jm=jmv(jc)       
          jp=jpv(jc)       
          jpp=jc+1
          udx3=dx3/(g3rc(kc))
          udx2=dx2/g2rc(jc)
          udx3p=dx3/(g3rc(kp))
          udx3pp=dx3/(g3rc(kpp))
          udx2p=dx2/g2rc(jpp)
            do ic=1,n1m 
              isy=isym(ic)
              dq2x3=(q2(ic,jc,kc)-q2(isy,jc,kc))*udx3/rc(jc)**2
              dq3x2=(rm(jc)*(q3(ic,jc,kp)-q3(isy,jc,kp))*0.5
     %              -rm(jm)*(q3(ic,jm,kp)-q3(isy,jm,kp))*0.5)*udx2    
     %                /sinpc(kp)
              voraz=(dq3x2-dq2x3)/rc(jc)
              dq2x3=0.
              dq3x2=(0.
     %              -rm(jc)*(q3(ic,jc,kp)-q3(isy,jc,kp))*0.5)*2.*udx2p
     %                /sinpc(kp)
              vorazj=(dq3x2-dq2x3)/rc(jpp)
              dq2x3=(q2(ic,jc,kp)-q2(ic,jc,kc))*udx3p/rc(jc)**2
              dq3x2=(rm(jc)*q3(ic,jc,kpp)-rm(jm)*q3(ic,jm,kpp))*udx2    
     %                /sinpc(kp)
              vorazk=(dq3x2-dq2x3)/rc(jc)
              dq2x3=0.
              dq3x2=(0.            -rm(jc)*q3(ic,jc,kpp))*2.*udx2p   
     %                /sinpc(kp)
              vorazjk=(dq3x2-dq2x3)/rc(jpp)
              om1=0.25*(voraz+vorazj+vorazk+vorazjk)
c
              disski = disski +om1**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
            end do
c South pole
        kc=n3m
        km=kmv(kc)           
        kp=kpv(kc)           
        kpp=kc+1
        do jc=2,n2m-1     
          jm=jmv(jc)       
          jp=jpv(jc)       
          jpp=jc+1
          udx3=dx3/(g3rc(kc))
          udx2=dx2/g2rc(jc)
          udx3p=dx3/(g3rc(kpp))
          udx3pp=dx3/(g3rc(kpp))
          udx2p=dx2/g2rc(jpp)
            do ic=1,n1m 
              isy=isym(ic)
              dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*udx3/rc(jc)**2
              dq3x2=(rm(jc)*q3(ic,jc,kc)-rm(jm)*q3(ic,jm,kc))*udx2    
     %                /sinpc(kc)
              voraz=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kc)-q2(ic,jpp,km))*udx3/rc(jpp)**2
              dq3x2=(rm(jp)*q3(ic,jp,kc)-rm(jc)*q3(ic,jc,kc))*udx2p   
     %                /sinpc(kc)
              vorazj=(dq3x2-dq2x3)/rc(jpp)
              dq2x3=(q2(isy,jc,kc)-q2(ic,jc,kc))*udx3p/rc(jc)**2
              dq3x2=(rm(jc)*(q3(ic,jc,kc)-q3(isy,jc,kc))*0.5
     %              -rm(jm)*(q3(ic,jm,kc)-q3(isy,jm,kc))*0.5)*udx2    
     %                /sinpc(kc)
              vorazk=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kp)-q2(ic,jpp,kc))*udx3pp/rc(jpp)**2
              dq3x2=(rm(jp)*(q3(ic,jp,kc)-q3(isy,jp,kc))*0.5
     %              -rm(jc)*(q3(ic,jc,kc)-q3(isy,jc,kc))*0.5)*udx2p   
     %                /sinpc(kc)
              vorazjk=(dq3x2-dq2x3)/rc(jpp)
              om1=0.25*(voraz+vorazj+vorazk+vorazjk)
c             write(6,*) ic,jc,kc,voraz,vorazj,vorazk,vorazjk
              disski = disski +om1**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
            end do
          end do
c Inner shell
          jc=1
          jm=jmv(jc)       
          jp=jpv(jc)       
          jpp=jc+1
          udx3=dx3/(g3rc(kc))
          udx2=dx2/g2rc(jc)
          udx3p=dx3/(g3rc(kp))
          udx3pp=dx3/(g3rc(kpp))
          udx2p=dx2/g2rc(jpp)
            do ic=1,n1m 
              isy=isym(ic)
              dq2x3=0.
              dq3x2=(rm(jc)*q3(ic,jc,kc)-0.)*2.*udx2/sinpc(kc)
              voraz=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(ic,jpp,kc)-q2(ic,jpp,km))*udx3p/rc(jpp)**2
              dq3x2=(rm(jp)*q3(ic,jp,kc)-rm(jc)*q3(ic,jc,kc))*udx2p   
     %                /sinpc(kc)
              vorazj=(dq3x2-dq2x3)/rc(jpp)
              dq2x3=0.
              dq3x2=(rm(jc)*(q3(ic,jc,kc)-q3(isy,jc,kc))*0.5
     %               -0.)*2.*udx2/sinpc(kp)
              vorazk=(dq3x2-dq2x3)/rc(jc)
              dq2x3=(q2(isy,jpp,kc)-q2(ic,jpp,kc))*udx3pp/rc(jpp)**2
              dq3x2=(rm(jp)*(q3(ic,jp,kc)-q3(isy,jp,kc))*0.5
     %              -rm(jc)*(q3(ic,jc,kc)-q3(isy,jc,kc))*0.5)*udx2p   
     %                /sinpc(kc)
              vorazjk=(dq3x2-dq2x3)/rc(jpp)
              om1=0.25*(voraz+vorazj+vorazk+vorazjk)
c             write(6,*) ic,jc,kc,voraz,vorazj,vorazk,vorazjk
              disski = disski +om1**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
            end do
c Outer shell
          jc=n2m
          jm=jmv(jc)       
          jp=jpv(jc)       
          jpp=jc+1
          udx3=dx3/(g3rc(kc))
          udx2=dx2/g2rc(jc)
          udx3p=dx3/(g3rc(kp))
          udx3pp=dx3/(g3rc(kpp))
          udx2p=dx2/g2rc(jpp)
            do ic=1,n1m 
              isy=isym(ic)
              dq2x3=(q2(ic,jc,kc)-q2(ic,jc,km))*udx3/rc(jc)**2
              dq3x2=(rm(jc)*q3(ic,jc,kc)-rm(jm)*q3(ic,jm,kc))*udx2    
     %                /sinpc(kc)
              voraz=(dq3x2-dq2x3)/rc(jc)
              dq2x3=0.
              dq3x2=(0.-rm(jc)*q3(ic,jc,kc))*2.*udx2p/sinpc(kc)
              vorazj=(dq3x2-dq2x3)/rc(jpp)
              dq2x3=(q2(isy,jc,kc)-q2(ic,jc,kc))*udx3/rc(jc)**2
              dq3x2=(rm(jc)*(q3(ic,jc,kc)-q3(isy,jc,kc))*0.5
     %              -rm(jm)*(q3(ic,jm,kc)-q3(isy,jc,kc))*0.5)*udx2    
     %                /sinpc(kc)
              vorazk=(dq3x2-dq2x3)/rc(jc)
              dq2x3=0.
              dq3x2=(0.-rm(jc)*(q3(ic,jc,kc)-q3(isy,jc,kc))*0.5)
     %              *2.*udx2p/sinpc(kc)
              vorazjk=(dq3x2-dq2x3)/rc(jpp)
              om1=0.25*(voraz+vorazj+vorazk+vorazjk)
c             write(6,*) ic,jc,kc,voraz,vorazj,vorazk,vorazjk
              disski = disski +om1**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
            end do
c
c   ! vorticita' radiale !  
      do 4 kc=2,n3m-1
        km=kmv(kc)                         
        kp=kpv(kc)                         
        kpp=kc+1
        udx3=dx3/(g3rc(kc))
        udx3p=dx3/(g3rc(kpp))
        do 4 jc=1,n2m                     
         fact = 1./(rm(jc)*sinpc(kc))
         factk= 1./(rm(jc)*sinpc(kp))
          do 4 ic=1,n1m                    
           im=imv(ic)                    
           ip=ipv(ic)                    
           dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/sinpc(kc)
           dq1x3=(q1(ic,jc,kc)*sinpm(kc)-q1(ic,jc,km)*sinpm(km))*udx3
           vorr=(dq1x3-dq3x1)*fact
           dq3x1=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1/sinpc(kc)
           dq1x3=(q1(ip,jc,kc)*sinpm(kc)-q1(ip,jc,km)*sinpm(km))*udx3
           vorri=(dq1x3-dq3x1)*fact
           dq3x1=(q3(ic,jc,kpp)-q3(im,jc,kpp))*dx1/sinpc(kp)
           dq1x3=(q1(ic,jc,kp)*sinpm(kp)-q1(ic,jc,kc)*sinpm(kc))*udx3p
           vorrk=(dq1x3-dq3x1)*factk
           dq3x1=(q3(ip,jc,kpp)-q3(ic,jc,kpp))*dx1/sinpc(kp)
           dq1x3=(q1(ip,jc,kp)*sinpm(kp)-q1(ip,jc,kc)*sinpm(kc))*udx3p
           vorrik=(dq1x3-dq3x1)*factk
           om2=0.25*(vorr+vorri+vorrk+vorrik)
              disski = disski +om2**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
    4 continue                   
c   North pole
      kc=1
        km=kmv(kc)                         
        kp=kpv(kc)                         
        kpp=kc+1
        udx3=dx3/(g3rc(kc))
        udx3p=dx3/(g3rc(kpp))
        do 44 jc=1,n2m                     
         fact = 1./(rm(jc)*sinpc(kp))
         factk= 1./(rm(jc)*sinpc(kp))
          do 44 ic=1,n1m                    
           im=imv(ic)                    
           ip=ipv(ic)                    
           isy=isym(ic)
           ism=isym(im)
           isp=isym(ip)
           dq3x1=((q3(ic,jc,kp)-q3(isy,jc,kp))*0.5
     %           -(q3(im,jc,kp)-q3(ism,jc,kp))*0.5)*dx1/sinpc(kp)
           dq1x3=(q1(ic,jc,kc)*sinpm(kc)+q1(isy,jc,kc)*sinpm(kc))*udx3
           vorr=(dq1x3-dq3x1)*fact
           dq3x1=((q3(ip,jc,kp)-q3(isp,jc,kp))*0.5
     %           -(q3(ic,jc,kp)-q3(isy,jc,kp))*0.5)*dx1/sinpc(kp)
           dq1x3=(q1(ip,jc,kc)*sinpm(kc)+q1(isp,jc,kc)*sinpm(kc))*udx3
           vorri=(dq1x3-dq3x1)*fact
           dq3x1=(q3(ic,jc,kpp)-q3(im,jc,kpp))*dx1/sinpc(kp)
           dq1x3=(q1(ic,jc,kp)*sinpm(kp)-q1(ic,jc,kc)*sinpm(kc))*udx3p
           vorrk=(dq1x3-dq3x1)*factk
           dq3x1=(q3(ip,jc,kpp)-q3(ic,jc,kpp))*dx1/sinpc(kp)
           dq1x3=(q1(ip,jc,kp)*sinpm(kp)-q1(ip,jc,kc)*sinpm(kc))*udx3p
           vorrik=(dq1x3-dq3x1)*factk
           om2=0.25*(vorr+vorri+vorrk+vorrik)
              disski = disski +om2**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
   44 continue                   
c  South pole
        kc=n3m
        km=kmv(kc)                         
        kp=kpv(kc)                         
        kpp=kc+1
        udx3=dx3/(g3rc(kc))
        udx3p=dx3/(g3rc(kpp))
        do 45 jc=1,n2m                     
         fact = 1./(rm(jc)*sinpc(kc))
         factk= 1./(rm(jc)*sinpc(kc))
          do 45 ic=1,n1m                    
           im=imv(ic)                    
           ip=ipv(ic)                    
           isy=isym(ic)
           ism=isym(im)
           isp=isym(ip)
           dq3x1=(q3(ic,jc,kc)-q3(im,jc,kc))*dx1/sinpc(kc)
           dq1x3=(q1(ic,jc,kc)*sinpm(kc)-q1(ic,jc,km)*sinpm(km))*udx3
           vorr=(dq1x3-dq3x1)*fact
           dq3x1=(q3(ip,jc,kc)-q3(ic,jc,kc))*dx1/sinpc(kc)
           dq1x3=(q1(ip,jc,kc)*sinpm(kc)-q1(ip,jc,km)*sinpm(km))*udx3
           vorri=(dq1x3-dq3x1)*fact
           dq3x1=((q3(ic,jc,kc)-q3(isy,jc,kc))*0.5
     %           -(q3(im,jc,kc)-q3(ism,jc,kc))*0.5)*dx1/sinpc(kc)
           dq1x3=(-q1(isy,jc,kc)-q1(ic,jc,kc))*sinpm(kc)*udx3p
           vorrk=(dq1x3-dq3x1)*factk
           dq3x1=((q3(ip,jc,kc)-q3(isp,jc,kc))*0.5
     %           -(q3(ic,jc,kc)-q3(isy,jc,kc))*0.5)*dx1/sinpc(kc)
           dq1x3=((q1(ip,jc,kc)-q1(isp,jc,kc))*0.5*sinpm(kc)
     %            -q1(ip,jc,kc)*sinpm(kc))*udx3p
           vorrik=(dq1x3-dq3x1)*factk
           om2=0.25*(vorr+vorri+vorrk+vorrik)
              disski = disski +om2**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
   45 continue                   
c  
c   ! vorticita' latitude !
      do kc=1,n3m                         
        do jc=2,n2m-1
          jm=jmv(jc)                           
          jp=jpv(jc)                           
          jpp=jp+1
          udx2 = dx2/g2rc(jc)                               
          udx2p= dx2/g2rc(jp)                               
            do ic=1,n1m                       
            im=imv(ic)                           
            ip=ipv(ic)                           
            dq1x2=(q1(ic,jc,kc)*rm(jc)-q1(ic,jm,kc)*rm(jm)) 
     %           *udx2
            dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1/(sinpm(kc)*rc(jc)**2)
            vorz=(dq2x1-dq1x2)/rc(jc) 
            dq1x2=(q1(ip,jc,kc)*rm(jc)-q1(ip,jm,kc)*rm(jm)) 
     %           *udx2
            dq2x1=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1/(sinpm(kc)*rc(jc)**2)
            vorzi=(dq2x1-dq1x2)/rc(jc) 
            dq1x2=(q1(ic,jp,kc)*rm(jp)-q1(ic,jc,kc)*rm(jc)) 
     %           *udx2p
            dq2x1=(q2(ic,jpp,kc)-q2(im,jpp,kc))*dx1
     %           /(sinpm(kc)*rc(jpp)**2)
            vorzj=(dq2x1-dq1x2)/rc(jpp) 
            dq1x2=(q1(ip,jp,kc)*rm(jp)-q1(ip,jc,kc)*rm(jc)) 
     %           *udx2p
            dq2x1=(q2(ip,jpp,kc)-q2(ic,jpp,kc))*dx1
     %           /(sinpm(kc)*rc(jpp)**2)
            vorzij=(dq2x1-dq1x2)/rc(jpp) 
            om3=0.25*(vorz+vorzi+vorzj+vorzij)
              disski = disski +om3**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
       end do
c Inner shell
         jc=1
          jm=jmv(jc)                           
          jp=jpv(jc)                           
          jpp=jp+1
          udx2 = dx2/g2rc(jc)                               
          udx2p= dx2/g2rc(jp)                               
            do ic=1,n1m                       
            im=imv(ic)                           
            ip=ipv(ic)                           
            dq1x2=(q1(ic,jc,kc)*rm(jc)-0.)*2.*udx2
            dq2x1=0.
            vorz=(dq2x1-dq1x2)/rc(jc) 
            dq1x2=(q1(ip,jc,kc)*rm(jc)-0.)*2.*udx2
            dq2x1=0.
            vorzi=(dq2x1-dq1x2)/rc(jc) 
            dq1x2=(q1(ic,jp,kc)*rm(jp)-q1(ic,jc,kc)*rm(jc)) 
     %           *udx2p
            dq2x1=(q2(ic,jpp,kc)-q2(im,jpp,kc))*dx1
     %           /(sinpm(kc)*rc(jpp)**2)
            vorzj=(dq2x1-dq1x2)/rc(jpp) 
            dq1x2=(q1(ip,jp,kc)*rm(jp)-q1(ip,jc,kc)*rm(jc)) 
     %           *udx2p
            dq2x1=(q2(ip,jpp,kc)-q2(ic,jpp,kc))*dx1
     %           /(sinpm(kc)*rc(jpp)**2)
            vorzij=(dq2x1-dq1x2)/rc(jpp) 
            om3=0.25*(vorz+vorzi+vorzj+vorzij)
              disski = disski +om3**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
c Outer shell
          jc=n2m
          jm=jmv(jc)                           
          jp=jpv(jc)                           
          jpp=jp+1
          udx2 = dx2/g2rc(jc)                               
          udx2p= dx2/g2rc(jp)                               
            do ic=1,n1m                       
            im=imv(ic)                           
            ip=ipv(ic)                           
            dq1x2=(q1(ic,jc,kc)*rm(jc)-q1(ic,jm,kc)*rm(jm)) 
     %           *udx2
            dq2x1=(q2(ic,jc,kc)-q2(im,jc,kc))*dx1/(sinpm(kc)*rc(jc)**2)
            vorz=(dq2x1-dq1x2)/rc(jc) 
            dq1x2=(q1(ip,jc,kc)*rm(jc)-q1(ip,jm,kc)*rm(jm)) 
     %           *udx2
            dq2x1=(q2(ip,jc,kc)-q2(ic,jc,kc))*dx1/(sinpm(kc)*rc(jc)**2)
            vorzi=(dq2x1-dq1x2)/rc(jc) 
            dq1x2=(0.-q1(ic,jc,kc)*rm(jc))*2.*udx2p
            dq2x1=0.
            vorzj=(dq2x1-dq1x2)/rc(jpp) 
            dq1x2=(0.-q1(ip,jc,kc)*rm(jc))*2.*udx2p
            dq2x1=0.
            vorzij=(dq2x1-dq1x2)/rc(jpp) 
            om3=0.25*(vorz+vorzi+vorzj+vorzij)
              disski = disski +om3**2 
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
c
      end do
c     write(6,*) ' E ',disski
c
      disski = disski/vol
c
c     Temperature gradient
c
      disste = 0.
      do kc=2,n3m-1
       kp=kpv(kc)
       km=kmv(kc)
       do jc=2,n2m-1
        jp=jpv(jc)
        jm=jmv(jc)
        udx2p=dx2/g2rc(jp)
        udx2=dx2/g2rc(jc)
        udx3p=dx3/(g3rc(kp)*rm(jc))
        udx3=dx3/(g3rc(kc)*rm(jc))
        fact=1./(rm(jc)*sinpm(kc))
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)
         dtdth = (dens(ip,jc,kc)-dens(im,jc,kc))*0.5*dx1*fact
         dtdr  = ((dens(ic,jp,kc)-dens(ic,jc,kc))*udx2p+
     %            (dens(ic,jc,kc)-dens(ic,jm,kc))*udx2)*0.5
         dtdp  = ((dens(ic,jc,kp)-dens(ic,jc,kc))*udx3p+
     %            (dens(ic,jc,kc)-dens(ic,jc,km))*udx3)*0.5
         disste = disste +
     %   (dtdth**2+dtdr**2+dtdp**2)
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
       end do
      end do
c North pole
      kc=1
       kp=kpv(kc)
       km=kmv(kc)
       do jc=2,n2m-1
        jp=jpv(jc)
        jm=jmv(jc)
        udx2p=dx2/g2rc(jp)
        udx2=dx2/g2rc(jc)
        udx3p=dx3/(g3rc(kp)*rm(jc))
        udx3=dx3/(g3rc(kc)*rm(jc))
        fact=1./(rm(jc)*sinpm(kc))
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)
         isy=isym(ic)
         dtdth = (dens(ip,jc,kc)-dens(im,jc,kc))*0.5*dx1*fact
         dtdr  = ((dens(ic,jp,kc)-dens(ic,jc,kc))*udx2p+
     %            (dens(ic,jc,kc)-dens(ic,jm,kc))*udx2)*0.5
         dtdp  = ((dens(ic,jc,kp)-dens(ic,jc,kc))*udx3p+
     %            (dens(ic,jc,kc)-dens(isy,jc,kc))*udx3)*0.5
         disste = disste +
     %   (dtdth**2+dtdr**2+dtdp**2)
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
       end do
c South pole
      kc=n3m
       kp=kpv(kc)
       kpp=kc+1
       km=kmv(kc)
       do jc=2,n2m-1
        jp=jpv(jc)
        jm=jmv(jc)
        udx2p=dx2/g2rc(jp)
        udx2=dx2/g2rc(jc)
        udx3p=dx3/(g3rc(kpp)*rm(jc))
        udx3=dx3/(g3rc(kc)*rm(jc))
        fact=1./(rm(jc)*sinpm(kc))
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)
         isy=isym(ic)
         dtdth = (dens(ip,jc,kc)-dens(im,jc,kc))*0.5*dx1*fact
         dtdr  = ((dens(ic,jp,kc)-dens(ic,jc,kc))*udx2p+
     %            (dens(ic,jc,kc)-dens(ic,jm,kc))*udx2)*0.5
         dtdp  = ((dens(isy,jc,kc)-dens(ic,jc,kc))*udx3p+
     %            (dens(ic,jc,kc)-dens(ic,jc,km))*udx3)*0.5
         disste = disste +
     %   (dtdth**2+dtdr**2+dtdp**2)
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
       end do
c inner shell
      do kc=2,n3m-1
       kp=kpv(kc)
       km=kmv(kc)
        jc=1
        jp=jpv(jc)
        jm=jmv(jc)
        udx3p=dx3/(g3rc(kp)*rm(jc))
        udx3=dx3/(g3rc(kc)*rm(jc))
        fact=1./(rm(jc)*sinpm(kc))
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)
         dtdth = (dens(ip,jc,kc)-dens(im,jc,kc))*0.5*dx1*fact
         dtdr  =((dens(ic,jp,kc)-dens(ic,jc,kc))*dx2/g2rc(2)
     %          +(dens(ic,jc,kc)-denbs(ic,kc))*dx2/g2rc(1)*2.)*0.5
         dtdp  = ((dens(ic,jc,kp)-dens(ic,jc,kc))*udx3p+
     %            (dens(ic,jc,kc)-dens(ic,jc,km))*udx3)*0.5
         disste = disste +
     %   (dtdth**2+dtdr**2+dtdp**2)
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
      end do
c inner shell North
       kc=1
       kp=kpv(kc)
       km=kmv(kc)
        jc=1
        jp=jpv(jc)
        jm=jmv(jc)
        udx3p=dx3/(g3rc(kp)*rm(jc))
        udx3=dx3/(g3rc(kc)*rm(jc))
        fact=1./(rm(jc)*sinpm(kc))
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)
         isy=isym(ic)
         dtdth = (dens(ip,jc,kc)-dens(im,jc,kc))*0.5*dx1*fact
         dtdr  =((dens(ic,jp,kc)-dens(ic,jc,kc))*dx2/g2rc(2)
     %          +(dens(ic,jc,kc)-denbs(ic,kc))*dx2/g2rc(1)*2.)*0.5
         dtdp  = ((dens(ic,jc,kp)-dens(ic,jc,kc))*udx3p+
     %            (dens(ic,jc,kc)-dens(isy,jc,kc))*udx3)*0.5
         disste = disste +
     %   (dtdth**2+dtdr**2+dtdp**2)
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
c inner shell South
       kc=n3m
       kp=kpv(kc)
       km=kmv(kc)
       kpp=kc+1
        jc=1
        jp=jpv(jc)
        jm=jmv(jc)
        udx3p=dx3/(g3rc(kpp)*rm(jc))
        udx3=dx3/(g3rc(kc)*rm(jc))
        fact=1./(rm(jc)*sinpm(kc))
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)
         isy=isym(ic)
         dtdth = (dens(ip,jc,kc)-dens(im,jc,kc))*0.5*dx1*fact
         dtdr  =((dens(ic,jp,kc)-dens(ic,jc,kc))*dx2/g2rc(2)
     %          +(dens(ic,jc,kc)-denbs(ic,kc))*dx2/g2rc(1)*2.)*0.5
         dtdp  = ((dens(isy,jc,kc)-dens(ic,jc,kc))*udx3p+
     %            (dens(ic,jc,kc)-dens(ic,jc,km))*udx3)*0.5
         disste = disste +
     %   (dtdth**2+dtdr**2+dtdp**2)
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
c outer shell
      do kc=2,n3m-1
       kp=kpv(kc)
       km=kmv(kc)
        jc=n2m
        jp=jpv(jc)
        jm=jmv(jc)
        udx2=0.5*dx2/g2rm(jc)
        udx3p=dx3/(g3rc(kp)*rm(jc))
        udx3=dx3/(g3rc(kc)*rm(jc))
        fact=1./(rm(jc)*sinpm(kc))
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)
         dtdth = (dens(ip,jc,kc)-dens(im,jc,kc))*0.5*dx1*fact
         dtdr  =((dens(ic,jc,kc)-dens(ic,jm,kc))*dx2/g2rc(jc)
     %          +(denbn(ic,kc)-dens(ic,jc,kc))*dx2/g2rc(n2)*2.)*0.5
         dtdp  = ((dens(ic,jc,kp)-dens(ic,jc,kc))*udx3p+
     %            (dens(ic,jc,kc)-dens(ic,jc,km))*udx3)*0.5
         disste = disste +
     %   (dtdth**2+dtdr**2+dtdp**2)
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
      end do
c outer shell North
       kc=1
       kp=kpv(kc)
       km=kmv(kc)
        jc=n2m
        jp=jpv(jc)
        jm=jmv(jc)
        udx2=0.5*dx2/g2rm(jc)
        udx3p=dx3/(g3rc(kp)*rm(jc))
        udx3=dx3/(g3rc(kc)*rm(jc))
        fact=1./(rm(jc)*sinpm(kc))
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)
         isy=isym(ic)
         dtdth = (dens(ip,jc,kc)-dens(im,jc,kc))*0.5*dx1*fact
         dtdr  =((dens(ic,jc,kc)-dens(ic,jm,kc))*dx2/g2rc(jc)
     %          +(denbn(ic,kc)-dens(ic,jc,kc))*dx2/g2rc(n2)*2.)*0.5
         dtdp  = ((dens(ic,jc,kp)-dens(ic,jc,kc))*udx3p+
     %            (dens(ic,jc,kc)-dens(isy,jc,kc))*udx3)*0.5
         disste = disste +
     %   (dtdth**2+dtdr**2+dtdp**2)
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
c outer shell South
       kc=n3m
       kp=kpv(kc)
       km=kmv(kc)
       kpp=kc+1
        jc=n2m
        jp=jpv(jc)
        jm=jmv(jc)
        udx2=0.5*dx2/g2rm(jc)
        udx3p=dx3/(g3rc(kpp)*rm(jc))
        udx3=dx3/(g3rc(kc)*rm(jc))
        fact=1./(rm(jc)*sinpm(kc))
        do ic=1,n1m
         ip=ipv(ic)
         im=imv(ic)
         isy=isym(ic)
         dtdth = (dens(ip,jc,kc)-dens(im,jc,kc))*0.5*dx1*fact
         dtdr  =((dens(ic,jc,kc)-dens(ic,jm,kc))*dx2/g2rc(jc)
     %          +(denbn(ic,kc)-dens(ic,jc,kc))*dx2/g2rc(n2)*2.)*0.5
         dtdp  = ((dens(isy,jc,kc)-dens(ic,jc,kc))*udx3p+
     %            (dens(ic,jc,kc)-dens(ic,jc,km))*udx3)*0.5
         disste = disste +
     %   (dtdth**2+dtdr**2+dtdp**2)
     %   *g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
        end do
c
      disste = disste/vol
c 
       open(38,file='dissipa.out',status='unknown',position='append')
       write(38,540) time, disski, disste
c      pause
 540   format(3(1x,e14.6))
       close(38)
c
      return
      end
c
c
c***********************************************************************
      subroutine outh(time,cflm)
      include'param.f'
      character*14 filth,filte,filba,filen,filso,filve,filnu,filbaper
      common/names/filth,filte,filba,filen,filso,filve,filnu,filbaper
c
      dpi=2.*2.*asin(1.)
      write(6,159)ntime,time,
     % vmax(1),imaxv(1),jmaxv(1),kmaxv(1),
     % vmax(2),imaxv(2),jmaxv(2),kmaxv(2),
     % vmax(3),imaxv(3),jmaxv(3),kmaxv(3),
     % cflm,qqmax,imxq,jmxq,kmxq,densm,denmax,denmin
c      open(32,file=filth,status='unknown')
      write(32,159)ntime,time,
     % vmax(1),imaxv(1),jmaxv(1),kmaxv(1),
     % vmax(2),imaxv(2),jmaxv(2),kmaxv(2),
     % vmax(3),imaxv(3),jmaxv(3),kmaxv(3),
     % cflm,qqmax,imxq,jmxq,kmxq,densm,denmax,denmin
c      close(32)
  159 format(1x,i5,2x,e10.4,4x,/,4x,
     % 3(e9.3,1x,i3,1x,i3,1x,i3,1x),/,4x,e9.3,
     % 2x,e9.3,1x,3(i3,1x),3(e12.6,1x))
c
c      open(97,file=filte,status='unknown',position='append')
       write(97,546) time, densm, anusslow, anussupp
 546   format(4(1x,e14.6))
c      close(97)
       eene = 0.
       do k=2,n3m-1
        kp=kpv(k)
        do j=1,n2m
         jp=jpv(j)
         do i=1,n1m
          ip = ipv(i)
          eene = eene + 0.5*(
     %         ((q1(i,j,k)+q1(ip,j,k))*0.5)**2+
     %         ((q2(i,j,k)+q2(i,jp,k))*0.5/rm(j)**2)**2+   ! senza **2 su rm(j)?
     %         ((q3(i,j,k)+q3(i,j,kp))*0.5/sinpm(k))**2)
     %          *g2rm(j)*g3rm(k)/(dx1*dx2*dx3)*rm(j)**2*sinpm(k)
          end do
        end do
      end do
C         attorno ai poli.

C           north pole

        k= 2
        kp=kpv(k)
        do j=1,n2m
         jp=jpv(j)
         do i=1,n1m
          ip = ipv(i)
            isy = isym(i)
      q3s = 0.5*(q3(i,j,kp) - q3(isy,j,kp))

          eene = eene + 0.5*(
     %         ((q1(i,j,k)+q1(ip,j,k))*0.5)**2+
     %         ((q2(i,j,k)+q2(i,jp,k))*0.5/rm(j)**2)**2+   ! senza **2 su rm(j)?
     %         ((q3s+q3(i,j,kp))*0.5/sinpm(k))**2)
     %          *g2rm(j)*g3rm(k)/(dx1*dx2*dx3)*rm(j)**2*sinpm(k)
          end do
        end do


C           south pole

        k= n3m
        kp=kpv(k)
        do j=1,n2m
         jp=jpv(j)
         do i=1,n1m
          ip = ipv(i)

            isy = isym(i)
        q3s = 0.5*(q3(i,j,k) - q3(isy,j,k))
          eene = eene + 0.5*(
     %         ((q1(i,j,k)+q1(ip,j,k))*0.5)**2+
     %         ((q2(i,j,k)+q2(i,jp,k))*0.5/rm(j)**2)**2+   ! senza **2 su rm(j)?
     %         ((q3(i,j,k)+q3s)*0.5/sinpm(k))**2)
     %          *g2rm(j)*g3rm(k)/(dx1*dx2*dx3)*rm(j)**2*sinpm(k)
          end do
        end do


       open(38,file='enerv.out',status='unknown',position='append')
       write(38,540) time, eene
       write(*,*) ' T = ',time,' E = ', eene
c      pause
 540   format(2(1x,e14.6))
       close(38)
c
c      open(96,file=filve,status='unknown')
       write(96,*) time, vmax(1), vmax(2), vmax(3)
c      close(96)
      return
      end
c
c
c***********************************************************************
      subroutine densmc
      include'param.f'
      real   dden,vol
      dden = 0.0d0
      vol = 0.0d0
c
      densm=0.
      denmax=-10.*10.
      denmin=10.*10.
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,rm,g2rm,g3rm)
!$OMP$ SHARED(dens,sinpm,dx1,dx2,dx3)
!$OMP$ PRIVATE(fac2,ic,jc,kc)
!$OMP$ REDUCTION(max:denmax)
!$OMP$ REDUCTION(min:denmin)
!$OMP$ REDUCTION(+:dden)
!$OMP$ REDUCTION(+:vol)
      do 11 kc=1,n3m
        do 11 jc=1,n2m
          fac2 = g2rm(jc)*g3rm(kc)/(dx1*dx2*dx3)*rm(jc)**2*sinpm(kc)
          do 11 ic=1,n1m
            denmax = max(denmax,dens(ic,jc,kc))
            denmin = min(denmin,dens(ic,jc,kc))
C           if(dens(ic,jc,kc).gt.denmax) denmax=dens(ic,jc,kc)
C           if(dens(ic,jc,kc).lt.denmin) denmin=dens(ic,jc,kc)
            dden=dden+(dens(ic,jc,kc)*fac2)
            vol=vol+(fac2)
   11 continue
!$OMP  END PARALLEL DO
      densm=(dden/vol)
c      write(6,*) 'in densmc densm,densmin,denmax',densm,denmin,denmax
c
c     COMPUTATION OF THE NUSSELT NUMBER AT THE 
c     INNER AND OUTER SPHERES
c
      anusslow = 0.
      anussupp = 0.
      del1 = rm(1)-rc(1)
      del  = rm(2)-rc(1)
      udel1q = 1./del1**2
      udelq = 1./del**2
      fcder = 1./(1./del1 - 1./del)
      del1n = -rm(n2m)+rc(n2)
      deln  = -rm(n2m-1)+rc(n2)
      udel1qn = 1./del1n**2
      udelqn = 1./deln**2
      fcdern = 1./(1./del1n - 1./deln)
      pi = 2. * asin(1.)
      surflow = 4.*pi * rint**2
      surfupp = 4.*pi * rext**2
      do k=1,n3m
       faclow = g3rm(k)/(dx1*dx3)*rc(1)**2*sinpm(k)
       facupp = g3rm(k)/(dx1*dx3)*rc(n2)**2*sinpm(k)
       do i=1,n1m
           grtlow = (dens(i,1,k)-denbs(i,k))/(rm(1)-rc(1))
c          grtlow = (dens(i,1,k)*udel1q-dens(i,2,k)*udelq
c    %              -denbs(i,k)*(udel1q-udelq) )
c    %              *fcder
           grtupp = (denbn(i,k)-dens(i,n2m,k))/(rc(n2)-rm(n2m))
c          grtupp = (dens(i,n2m,k)*udel1qn-dens(i,n2m-1,k)*udelqn
c    %              -denbn(i,k)*(udel1qn-udelqn) )
c    %              *fcdern
           anusslow = anusslow + grtlow*faclow
           anussupp = anussupp + grtupp*facupp
         end do
      end do
      anusslow = anusslow / surflow
      anussupp = anussupp / surfupp
      return
      end
