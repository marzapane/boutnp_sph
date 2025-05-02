c***********************************************************************
c     particular care is necessary for the metric quantities in the
c     non-linear terms
c***********************************************************************
      subroutine hdnlro
      include 'param.f'
c
c     h term for the q1 momentum equation at i,j+1/2,k+1/2
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,kmv,jpv,jmv,ipv,imv)
!$OMP$ SHARED(rm,zm,g2rm,g3rm)
!$OMP$ SHARED(usrm,usinpm,utanpm)
!$OMP$ SHARED(dx1,dx2,dx3)
!$OMP$ SHARED(hro,q1,q2,q3,dens)
!$OMP$ PRIVATE(kc,kmm,kpp,kp,jc,jmm,jpp,jp,jm,ic,ip,im)
!$OMP$ PRIVATE(udx1nl,fact2,fact3)
!$OMP$ PRIVATE(h11,h12,h13,hq1)
      do 10 kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 10 jc=1,n2m
      jmm=jmv(jc)
      jpp=jpv(jc)
      jp=jc+1
      jm=jc-1
      udx1nl = usrm(jc)*usinpm(kc)*dx1*0.25
      fact2  = usrm(jc)*usinpm(kc)
      fact3  = usrm(jc)*usinpm(kc)*utanpm(kc)
      do 10 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
c
c                 1      d  T q_t
c            ---------- -----------  azimuthal non-linear term
c             r sin phi  d   t
c
      h11=( q1(ip,jc,kc)*(dens(ip,jc,kc)+dens(ic,jc,kc))
     %     -q1(ic,jc,kc)*(dens(im,jc,kc)+dens(ic,jc,kc))
     %    )*dx1*0.5/(rm(jc)*sin(zm(kc)))
C    %    )*usrm(jc)*usinpm(kc)*dx1*0.5
c
c             1    d  T q_r
c            ---  ----------- 
c             r^2   d   r
c
c
      h12 =( q2(ic,jp,kc)*(dens(ic,jpp,kc)+dens(ic,jc,kc))
     %      -q2(ic,jc,kc)*(dens(ic,jc,kc)+dens(ic,jmm,kc))
     %     )*dx2/(g2rm(jc)*rm(jc)**2)*0.5 !udx2rmq(jc)*0.5
c
c   q1 q3 term
c
c
c          1      d  T q_p
c       -------  -----------
c      r sin phi  d   p
c
      h13=(q3(ic,jc,kp)*(dens(ic,jc,kpp)+dens(ic,jc,kc))
     %    -q3(ic,jc,kc)*(dens(ic,jc,kc)+dens(ic,jc,kmm))
     %    )*dx3/g3rm(kc)*0.5*usrm(jc)*usinpm(kc)  !udx3m(kc)*0.5*fact2
c
c
      hq1=h11+h12+h13
c
c
      hro(ic,jc,kc)=-hq1
   10 continue
!$OMP  END PARALLEL DO
c
      return
      end
c
c***********************************************************************
c     particular care is necessary for the metric quantities in the
c     non-linear terms
c***********************************************************************
      subroutine hdnl1
      include 'param.f'
c
c     h term for the q1 momentum equation at i,j+1/2,k+1/2
c
      integer om_th,om_z,ampl
      real ph_th,ph_z,forc_mod,rand,pi
      pi = 2. * asin (1.)
c     definition of the random coeff. for the forcing
      rand = 1.d0
      call random_number(rand)
      ampl = nint(n1/32.)
      om_th = n1/4 - ampl/2 + nint(rand * ampl)
      call random_number(rand)
      ampl = nint(n3/32.)
      om_z = n3/4 - ampl/2 + nint(rand * ampl/2)*2
      call random_number(rand)
      ph_th = rand * 2*pi
      call random_number(rand)
      ph_z = rand * 2*pi
      call random_number(rand)
      forc_mod = forc_scale * rand
      r_max = rc(n2)
      r_min = rc(1)

      ivisc1=1
      small = 1.e-10
      icon = 0
      usros = 1.d0/ros
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,kmv,jpv,jmv,ipv,imv)
!$OMP$ SHARED(rm,rc,thetac,zm,g2rm,g3rm,isym,ksym,ksyc)
!$OMP$ SHARED(usrm,usrc,usinpm,usinpc,sinpm,sinpc,cospm)
!$OMP$ SHARED(utanpm,qaxN,qaxS)
!$OMP$ SHARED(dq,q1,q2,q3,dens)
!$OMP$ SHARED(dx1,dx2,dx3,udx3m,udx1vis1,udx1vis3,pi)
!$OMP$ SHARED(fvisc,ivisc1,usros,icorio,ibuo,iforc)
!$OMP$ SHARED(forc_mod,om_th,om_z,ph_th,ph_z,r_min,r_max)
!$OMP$ PRIVATE(kc,kmm,kpp,kp,jc,jmm,jpp,jp,jm,ic,ip,im)
!$OMP$ PRIVATE(udx1nl,fact2,fact3,insy,insym,knsy)
!$OMP$ PRIVATE(ksm,ksy,ism,ismm,q3sy,q3sym,q3s,q3sm)
!$OMP$ PRIVATE(q2e,q2w,d11q2e,q3e,q3w,d11q3e)
!$OMP$ PRIVATE(h11,h12,h13,h12n,h13n,hq1,coriol,densit)
!$OMP$ PRIVATE(forc,om_tape)
      do 10 kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      do 10 jc=1,n2m
      jmm=jmv(jc)
      jpp=jpv(jc)
      jp=jc+1
      jm=jc-1
      udx1nl = usrm(jc)*usinpm(kc)*dx1*0.25
      fact2  = usrm(jc)*usinpm(kc)
      fact3  = usrm(jc)*usinpm(kc)*utanpm(kc)
      do 10 ic=1,n1m
      ip=ipv(ic)
      im=imv(ic)
      insy = isym(ic)
      insym = isym(im)
      knsy = ksym(kc)
c
c                 1      d  q_t q_t
c            ---------- -----------  azimuthal non-linear term
c             r sin phi  d   t
c
      h11=( (q1(ip,jc,kc)+q1(ic,jc,kc))*(q1(ip,jc,kc)+q1(ic,jc,kc))
     %     -(q1(im,jc,kc)+q1(ic,jc,kc))*(q1(im,jc,kc)+q1(ic,jc,kc))
     %    )*dx1*0.25/(rm(jc)*sin(zm(kc)))
C    %    )*usrm(jc)*usinpm(kc)*dx1*0.25
c
c             1    d  q_t q_r
c            ---  ----------- 
c             r^2   d   r
c
c
      h12 =( (q2(ic,jp,kc)+q2(im,jp,kc))*(q1(ic,jpp,kc)+q1(ic,jc,kc))
     %      -(q2(ic,jc,kc)+q2(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jmm,kc))
     %     )*dx2/(g2rm(jc)*rm(jc)**2)*0.25 !udx2rmq(jc)*0.25
c
c   q1 q3 term
c
c
c          1      d  q_t q_p
c       -------  -----------
c      r sin phi  d   p
c
      h13=((q3(ic,jc,kp)+q3(im,jc,kp))*(q1(ic,jc,kpp)+q1(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(im,jc,kc))*(q1(ic,jc,kc)+q1(ic,jc,kmm))
     %    )*dx3/g3rm(kc)*0.25*usrm(jc)*usinpm(kc)  !udx3m(kc)*0.25*fact2
c
c     non differentiated terms
c
c         q_t q_r
c        --------    
c          r^3
c
      h12n=  q1(ic,jc,kc)/rm(jc) * (
     %                  (q2(ic,jp,kc)+q2(im,jp,kc))/(rc(jp)**2)
     %                 +(q2(ic,jc,kc)+q2(im,jc,kc))/(rc(jc)**2) )*0.25
c
c             q_t q_p
c        ----------------
c       r sin(phi)tan(phi)
c
      if((kc.ne.1).and.(kc.ne.n3m)) then
       h13n=  q1(ic,jc,kc)*usrm(jc)*utanpm(kc) * (
     %          (q3(im,jc,kc)*usinpc(kc)+q3(im,jc,kp)*usinpc(kp))
     %         +(q3(ic,jc,kc)*usinpc(kc)+q3(ic,jc,kp)*usinpc(kp)) )*0.25
      else if(kc.eq.1) then 
      ksm = ksyc(kp)
      ism = isym(ic)
      ismm = isym(im)
c  qaxS(3,m1,m2)
      q3sy = qaxN(3,ic,jc)*sinpc(kp) ! (q3(ic,jc,kp)-q3(ism,jc,kp))*0.5
      q3sym = qaxN(3,im,jc)*sinpc(kp) ! (q3(im,jc,kp)-q3(ismm,jc,kp))*0.5
      h13n=  q1(ic,jc,kc)*usrm(jc)*utanpm(kc) * (
     %                  (q3sy+q3(im,jc,kp))*usinpc(kp)
     %                 +(q3sym+q3(ic,jc,kp))*usinpc(kp) )*0.25
      else if(kc.eq.n3m) then 
      ksm = ksyc(kc)
      ism = isym(ic)
      ismm = isym(im)
      q3sy = qaxS(3,ic,jc)*sinpc(kc) ! (q3(ic,jc,kc)-q3(ism,jc,kc))*0.5
      q3sym = qaxS(3,im,jc)*sinpc(kc) ! (q3(im,jc,kc)-q3(ismm,jc,kc))*0.5
      h13n=  q1(ic,jc,kc)*usrm(jc)*utanpm(kc) * (
     %                  (q3sy+q3(im,jc,kc))*usinpc(kc)
     %                 +(q3sym+q3(ic,jc,kc))*usinpc(kc) )*0.25
      end if
c     h13n = 0.
c
      hq1=h11+h12+h13+h12n+h13n
c
c                  2        d  q_r / r^2
c            -------------  -----------
c            Re r^2 sin phi  d   t
c
       q2e=(q2(ic,jp,kc)*usrc(jp)**2+q2(ic,jc,kc)*usrc(jc)**2)*0.5
       q2w=(q2(im,jp,kc)*usrc(jp)**2+q2(im,jc,kc)*usrc(jc)**2)*0.5
       d11q2e=(q2e-q2w)*udx1vis1(jc,kc)*float(ivisc1)*dx1*fvisc(jc,kc)
c
c            2 cos phi        d  q_p / sin phi
c            ------          ---------------
c            Re r sin^2 phi       d   t
c
      if((kc.ne.1).and.(kc.ne.n3m)) then
      q3e=(q3(ic,jc,kp)*usinpc(kp)+q3(ic,jc,kc)*usinpc(kc))*0.5
      q3w=(q3(im,jc,kp)*usinpc(kp)+q3(im,jc,kc)*usinpc(kc))*0.5
      d11q3e=(q3e-q3w)*udx1vis3(jc,kc)*float(ivisc1)*dx1*fvisc(jc,kc)
      else if(kc.eq.1) then 
c     isy = isym(ic)
      ksy = ksyc(kp)
      q3s = qaxN(3,ic,jc)*sinpc(kp) !0.5*(q3(ic,jc,kp)-q3(isy,jc,kp))
c     isy = isym(im)
      q3sm= qaxN(3,im,jc)*sinpc(kp) !0.5*(q3(im,jc,kp)-q3(isy,jc,kp))
      q3e=(q3(ic,jc,kp)+q3s)*usinpc(kp)*0.5
      q3w=(q3(im,jc,kp)+q3sm)*usinpc(kp)*0.5
      d11q3e=(q3e-q3w)*udx1vis3(jc,kc)*float(ivisc1)*dx1*fvisc(jc,kc)
      else if(kc.eq.n3m) then 
c     isy = isym(ic)
      ksy = ksyc(kc)
      q3s = qaxS(3,ic,jc)*sinpc(kc) ! 0.5*(q3(ic,jc,kc)-q3(isy,jc,kc))
c     isy = isym(im)
      q3sm= qaxS(3,im,jc)*sinpc(kc) ! 0.5*(q3(im,jc,kc)-q3(isy,jc,kc))
      q3e=(q3s+q3(ic,jc,kc))*usinpc(kc)*0.5
      q3w=(q3sm+q3(im,jc,kc))*usinpc(kc)*0.5
      d11q3e=(q3e-q3w)*udx1vis3(jc,kc)*float(ivisc1)*dx1*fvisc(jc,kc)
      end if
c
c   coriolis term
c
c
c              1    ( q_r          q_p        )
c           - ----- (-----sin(p)+ ------cos(p))  
c              Ro   ( r^2         sin(p)      )
c
      coriol= - (
     %          (q2(ic,jc,kc)+q2(im,jc,kc) +
     %           q2(ic,jp,kc)+q2(im,jp,kc))*usrm(jc)**2*sinpm(kc)
     %         +(q3(ic,jc,kc)+q3(im,jc,kc) +
     %           q3(ic,jc,kp)+q3(im,jc,kp))*usinpm(kc)*cospm(kc)
     %           )*0.25*usros
     %          *float(icorio)
c
c
c
c  add the buoiancy term  
c
      densit = 0.5*(dens(ic,jc,kc)+dens(im,jc,kc))*float(ibuo)

      dq(ic,jc,kc)=d11q2e+d11q3e-hq1+coriol
     %             -densit*(sin(thetac(ic))) !gravita'' diretta lungo x 
c
c     add the random forcing term
c
      if(iforc.eq.1) then
        om_tape = nint(om_th*sin(zm(kc)))
        forc = forc_mod
     %        * sin(om_tape * thetac(ic) + ph_th)
     %        * (rm(jc)-r_min)/(r_max-r_min)
     %        * cos(om_z * zm(kc) + ph_z)
        dq(ic,jc,kc) = dq(ic,jc,kc) + forc !simple oscillating random forcing
      end if

   10 continue
!$OMP  END PARALLEL DO
c
c     close(888)
c     pause
      return
      end
c
c***********************************************************************
      subroutine hdnl2
      include 'param.f'
c
c     h term for the q2 momentum equation at i+1/2,j,k+1/2
c
      integer om_th,om_z,ampl
      real ph_th,ph_z,forc_mod,rand,pi
      pi = 2. * asin (1.)
c     definition of the random coeff. for the forcing
      rand = 1.d0
      call random_number(rand)
      ampl = nint(n1/32.)
      om_th = n1/4 - ampl/2 + nint(rand * ampl)
      call random_number(rand)
      ampl = nint(n3/32.)
      om_z = n3/4 - ampl/2 + nint(rand * ampl/2)*2
      call random_number(rand)
      ph_th = rand * 2*pi
      call random_number(rand)
      ph_z = rand * 2*pi
      call random_number(rand)
      forc_mod = forc_scale * rand
      r_max = rc(n2)
      r_min = rc(1)

      ivisc2=1
      icon = 0
      usros = 1.d0/ros
      sudx1vis2 = 2.*dx1/ren
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,kmv,ipv,imv)
!$OMP$ SHARED(rm,rc,thetac,zm,g3rm,isym,ksym,ksyc)
!$OMP$ SHARED(usrc,usinpm,usinpc,sinpm,sinpc)
!$OMP$ SHARED(qaxN,qaxS,ren)
!$OMP$ SHARED(dph,q1,q2,q3,dens)
!$OMP$ SHARED(dx1,dx3,udx2c,udx3m,pi)
!$OMP$ SHARED(fvisc,ivisc2,usros,icorio,ibuo,iforc)
!$OMP$ SHARED(forc_mod,om_th,om_z,ph_th,ph_z,r_min,r_max)
!$OMP$ PRIVATE(kc,kmm,kpp,kp,jc,jp,jm,ic,imm,ipp)
!$OMP$ PRIVATE(udx1,fact,fact2,fact3,insy,knsy)
!$OMP$ PRIVATE(ksm,ism,q3sy,q3sym,q3s,q3n)
!$OMP$ PRIVATE(q1e,q1w,d11q1v,d33q3v)
!$OMP$ PRIVATE(h21,h22,h23,h11n,h33n,hq2,coriol,densit)
!$OMP$ PRIVATE(forc,om_tape)
      do 20 kc=1,n3m
      kmm=kmv(kc)
      kpp=kpv(kc)
      kp=kc+1
      fact2 = usinpm(kc)**2
      fact3 = usinpm(kc)
      do 20 jc=2,n2m
      jm=jc-1
      jp=jc+1
      udx1 = dx1*0.25
      fact = usrc(jc)*usinpm(kc)
      do 20 ic=1,n1m
      imm=imv(ic)
      ipp=ipv(ic)
      insy = isym(ic)
      knsy = ksym(kc)
c
c             1          d  q_t q_r
c           ---------  -----------
c           r sin phi    d   t
c
      h21=((q1(ipp,jc,kc)+q1(ipp,jm,kc))
     %     *(q2(ipp,jc,kc)+q2(ic,jc,kc))
     %     -(q1(ic,jc,kc)+q1(ic,jm,kc))
     %     *(q2(ic,jc,kc)+q2(imm,jc,kc)))*
     %     dx1*0.25*usrc(jc)*usinpm(kc)
c
c                 d  q_r q_r / r^2
c                ---------------
c                 d   r
c
      h22=( (q2(ic,jp,kc)*usrc(jp)**2+q2(ic,jc,kc)*usrc(jc)**2)
     %     *(q2(ic,jp,kc)+q2(ic,jc,kc))
     %     -(q2(ic,jc,kc)*usrc(jc)**2+q2(ic,jm,kc)*usrc(jm)**2)
     %     *(q2(ic,jc,kc)+q2(ic,jm,kc))
     %    )*udx2c(jc)*0.25
c
c             1          d  q_p q_r
c           ---------  -----------
c           r sin phi    d   p
c
      h23=((q3(ic,jc,kp)+q3(ic,jm,kp))*(q2(ic,jc,kpp)+q2(ic,jc,kc))
     %    -(q3(ic,jc,kc)+q3(ic,jm,kc))*(q2(ic,jc,kc)+q2(ic,jc,kmm))
     %    )*udx3m(kc)*0.25*usrc(jc)*usinpm(kc)
c
c   terms not differentiated
c
c             r*( q_t )^2
c
      h11n=((q1(ipp,jc,kc)+q1(ipp,jm,kc)
     %       +q1(ic,jc,kc)+q1(ic,jm,kc))*0.25)**2*rc(jc)
c
c
c
c             r*( q_p )^2/sin^2 phi
c
      if((kc.ne.1).and.(kc.ne.n3m)) then
      h33n=(((q3(ic,jc,kp)+q3(ic,jm,kp))*usinpc(kp)
     %       +(q3(ic,jc,kc)+q3(ic,jm,kc))*usinpc(kc))*0.25)**2*rc(jc)
      else if(kc.eq.1) then 
      ism = isym(ic)
      ksm = ksyc(kp)
      q3sy = qaxN(3,ic,jc)*sinpc(kp) !0.5*(q3(ic,jc,kp)-q3(isy,jc,kp))
c     q3sy = (q3(ic,jc,kp)-q3(ism,jc,kp))*0.5
      q3sym = qaxN(3,ic,jm)*sinpc(kp) !0.5*(q3(ic,jc,kp)-q3(isy,jc,kp))
c     q3sym = (q3(ic,jm,kp)-q3(ism,jm,kp))*0.5
      h33n=(((q3(ic,jc,kp)+q3(ic,jm,kp)+q3sy+q3sym)*usinpc(kp)
     %       )*0.25)**2*rc(jc)
      else if(kc.eq.n3m) then 
      ksm = ksyc(kc)
      ism = isym(ic)
      q3sy = qaxS(3,ic,jc)*sinpc(kc) !0.5*(q3(ic,jc,kp)-q3(isy,jc,kp))
c     q3sy = (q3(ic,jc,kc)-q3(ism,jc,kc))*0.5
      q3sym = qaxS(3,ic,jm)*sinpc(kc) !0.5*(q3(ic,jc,kp)-q3(isy,jc,kp))
c     q3sym = (q3(ic,jm,kc)-q3(ism,jm,kc))*0.5
      h33n=(((q3(ic,jc,kc)+q3(ic,jm,kc)+q3sy+q3sym)*usinpc(kc)
     %       )*0.25)**2*rc(jc)
      end if
c     h33n = 0.
c
      hq2=h21+h22+h23-h11n-h33n
c
c   off diagonal viscous terms
c
c
c                 2        d  q_t
c          - ----------  --------
c            Re sin phi    d   t
c
      q1e=(q1(ipp,jc,kc)+q1(ipp,jm,kc))*0.5
      q1w=(q1(ic,jc,kc)+q1(ic,jm,kc))*0.5
      d11q1v=-(q1e-q1w)*2.*dx1/ren*usinpm(kc)*float(ivisc2)*fvisc(jc,kc)
c
c                 2        d  q_p
c          - ----------  --------
c            Re sin phi    d   p
c
      q3s=(q3(ic,jc,kp)+q3(ic,jm,kp))*0.5
      q3n=(q3(ic,jc,kc)+q3(ic,jm,kc))*0.5
      d33q3v=-(q3s-q3n)*2.*dx3/g3rm(kc)/ren*usinpm(kc)
     %               *float(ivisc2)*fvisc(jc,kc)
c
c
c   coriolis term
c
c
c              1
c             -----   r^2 sin(p) q_t 
c              Ro
c
      coriol=(q1(ic,jc,kc)+q1(ipp,jc,kc)+q1(ic,jm,kc)+
     %        q1(ipp,jm,kc))*0.25*rc(jc)**2*sinpm(kc)
     %        *usros*float(icorio)
c
c
c  add the buoiancy term  
c
      densit = 0.5*(dens(ic,jc,kc)+dens(ic,jm,kc))*rc(jc)**2*float(ibuo)

      dph(ic,jc,kc)=d11q1v+d33q3v-hq2+coriol
     %               +densit*(sinpm(kc)*cos(thetam(ic))) !gravita'' allimeata in x
c    %               +densit*(-cospm(kc))         !gravita'' verticale
c    %               +densit*(rc(n2)/rc(jc))**2  !tutta la massa al centro
c    %               +densit                     !gravita' costante
c    %               +densit*(rc(jc)/rc(n2))     !self-gravitating sphere

c
c     add the random forcing term
c
      if(iforc.eq.1) then
        om_tape = nint(om_th*sin(zm(kc)))
        forc = forc_mod * rc(jc)**2
     %        * sin(om_tape * thetam(ic) + ph_th)
     %        * (rc(jc)-r_min)/(r_max-r_min)
     %        * cos(om_z * zm(kc) + ph_z)
        dph(ic,jc,kc) = dph(ic,jc,kc) + forc !simple oscillating random forcing
      end if

   20 continue
!$OMP  END PARALLEL DO
c
      return
      end
c
c***********************************************************************
      subroutine hdnl3
      include 'param.f'
c
c     h term for the q3 momentum equation at i+1/2,j+1/2,k
c
      integer om_th,om_z,ampl
      real ph_th,ph_z,forc_mod,rand
      pi = 2. * asin (1.)
c     definition of the random coeff. for the forcing
      rand = 1.d0
      call random_number(rand)
      ampl = nint(n1/32.)
      om_th = n1/4 - ampl/2 + nint(rand * ampl)
      call random_number(rand)
      ampl = nint(n3/32.)
      om_z = n3/4 - ampl/2 + nint(rand * ampl/2)*2
      call random_number(rand)
      ph_th = rand * 2*pi
      call random_number(rand)
      ph_z = rand * 2*pi
      call random_number(rand)
      forc_mod = forc_scale * rand
      r_max = rc(n2)
      r_min = rc(1)

      ivisc3 = 1
      small = 1.e-10
      icon = 0
      usros = 1.d0/ros
      sudx1vis2 = 2.*dx1/ren
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,jpv,jmv,ipv,imv)
!$OMP$ SHARED(rm,thetac,thetam,zz,zm,g2rm,g3rc,isym)
!$OMP$ SHARED(usrm,usrmq,usrc,usinpc,sinpc,cospc)
!$OMP$ SHARED(utanpc,qaxN,qaxS,ren)
!$OMP$ SHARED(qcap,q1,q2,q3,dens)
!$OMP$ SHARED(dx1,dx2,dx3,pi)
!$OMP$ SHARED(fvisc,ivisc3,usros,icorio,ibuo,iforc)
!$OMP$ SHARED(forc_mod,om_th,om_z,ph_th,ph_z,r_min,r_max)
!$OMP$ PRIVATE(kc,km,kp,jc,jmm,jpp,jp,ic,ipp,imm)
!$OMP$ PRIVATE(sudx1vis2,sudx1vis3,udx1,fact2,fact3,insy)
!$OMP$ PRIVATE(q3km,q3kp,q1e,q1w,d11q1v,q2s,q2n,d11q2v)
!$OMP$ PRIVATE(h31,h32,h33,h23n,h11n,hq3,coriol,densit)
!$OMP$ PRIVATE(forc,om_tape)
      do 30 kc=2,n3m
c     km=kmv(kc)
      km=kc-1
      kp=kc+1
      sudx1vis3 = 2.*dx3/g3rc(kc)/ren
      do 30 jc=1,n2m
      jp=jc+1
      jmm=jmv(jc)
      jpp=jpv(jc)
      udx1=dx1*usrm(jc)*usinpc(kc)*0.25
      fact3 = utanpc(kc)*usrmq(jc)
      fact2 = sinpc(kc)*usrmq(jc)
      do 30 ic=1,n1m
      imm=imv(ic)
      ipp=ipv(ic)
c
c                1        d  q_p q_t
c           ----------   -----------
c            r sin phi     d   t
c
      h31=((q1(ipp,jc,kc)+q1(ipp,jc,km))
     %        *(q3(ipp,jc,kc)+q3(ic,jc,kc))
     %        -(q1(ic,jc,kc)+q1(ic,jc,km))
     %        *(q3(ic,jc,kc)+q3(imm,jc,kc)))*
c    %        dx1*usrm(jc)*usinpc(kc)*0.25 
     %        dx1/(rm(jc)*sin(zz(kc)))*0.25 
c
c            1    d  q_p q_r
c           ---  -----------
c            r^2  d   r
c
      h32=(((q2(ic,jp,kc)+q2(ic,jp,km))
     %     *(q3(ic,jpp,kc)+q3(ic,jc,kc)))
     %    -((q2(ic,jc,kc)+q2(ic,jc,km))
     %     *(q3(ic,jc,kc)+q3(ic,jmm,kc))))*
     %     dx2/(g2rm(jc)*rm(jc)**2)*0.25
c
c            1    d  q_p q_p/sin phi
c          ----  -----------
c            r    d   p
c
      if((kc.ne.2).and.(kc.ne.n3m)) then
      h33=((q3(ic,jc,kp)+q3(ic,jc,kc))
     %    *(q3(ic,jc,kp)/sin(zz(kp))+q3(ic,jc,kc)/sin(zz(kc)))
     %    -(q3(ic,jc,kc)+q3(ic,jc,km))
     %    *(q3(ic,jc,kc)/sin(zz(kc))+q3(ic,jc,km)/sin(zz(km)))
     %    )*dx3/g3rc(kc)*usrm(jc)*0.25
c     h33=((q3(ic,jc,kp)+q3(ic,jc,kc))
c    %    *(q3(ic,jc,kp)*usinpc(kp)+q3(ic,jc,kc)*usinpc(kc))
c    %    -(q3(ic,jc,kc)+q3(ic,jc,km))
c    %    *(q3(ic,jc,kc)*usinpc(kc)+q3(ic,jc,km)*usinpc(km))
c    %    )*udx3c(kc)*usrm(jc)*0.25
     
      else if(kc.eq.2) then 
      insy = isym(ic)
c     q3km =(q3(ic,jc,kc)-q3(insy,jc,kc))*0.5
      q3km = qaxN(3,ic,jc)*sin(zz(kc))
      h33=((q3(ic,jc,kp)+q3(ic,jc,kc))
     %    *(q3(ic,jc,kp)/sin(zz(kp))+q3(ic,jc,kc)/sin(zz(kc)))
     %    -(q3(ic,jc,kc)+q3(ic,jc,km))
     %    *(q3(ic,jc,kc)+q3km)/sin(zz(kc))
     %    )*dx3/g3rc(kc)/rm(jc)*0.25
      else if(kc.eq.n3m) then 
      insy = isym(ic)
c     q3kp=(q3(ic,jc,kc)-q3(insy,jc,kc))*0.5
      q3kp = qaxS(3,ic,jc)*sin(zz(kc))
      h33=((q3(ic,jc,kp)+q3(ic,jc,kc))
     %    *(q3kp+q3(ic,jc,kc))/sin(zz(kc))
     %    -(q3(ic,jc,kc)+q3(ic,jc,km))
     %    *(q3(ic,jc,kc)/sin(zz(kc))+q3(ic,jc,km)/sin(zz(km)))
     %    )*dx3/g3rc(kc)/rm(jc)*0.25
      end if
c
c   terms not differentiated
c
c              q_p q_r
c              -------
c               r^3
c
      h23n=q3(ic,jc,kc)*((q2(ic,jp,kc)+q2(ic,jp,km))*usrc(jp)**2
     %       +(q2(ic,jc,kc)+q2(ic,jc,km))*usrc(jc)**2)*0.25*usrm(jc)
c
c              q_t^2 cos phi
c            - -------------
c                   r
c
      h11n= ((q1(ipp,jc,kc)+q1(ipp,jc,km)
     %       +q1(ic,jc,kc)+q1(ic,jc,km))*0.25)**2*usrm(jc)*cospc(kc)
c
      hq3=h31+h32+h33-h11n+h23n
c
c   off diagonal viscous terms
c
c
c            2 cos phi        d  q_t
c          - ----------     --------
c            Re r^2 sin phi   d   t
c
      q1e=(q1(ipp,jc,kc)+q1(ipp,jc,km))*0.5
      q1w=(q1(ic,jc,kc)+q1(ic,jc,km))*0.5
      d11q1v=-(q1e-q1w)*sudx1vis2*fact3*float(ivisc3)*fvisc(jc,kc)
c
c            2 sin phi   d  q_r/r^2
c            ----------  --------
c            Re r^2      d   p
c
      q2s=(q2(ic,jc,kc)*usrc(jc)**2+q2(ic,jp,kc)*usrc(jp)**2)*0.5
      q2n=(q2(ic,jc,km)*usrc(jc)**2+q2(ic,jp,km)*usrc(jp)**2)*0.5
      d11q2v= (q2s-q2n)*sudx1vis3*fact2*float(ivisc3)*fvisc(jc,kc)
c
c   coriolis term
c
c
c              1
c             -----   sin(p) cos(p) q_t 
c              Ro
c
      coriol=(q1(ic,jc,kc)+q1(ipp,jc,kc)+q1(ic,jc,km)+
     %        q1(ipp,jc,km))*0.25*cospc(kc)*sinpc(kc)
     %        *usros*float(icorio)
c
c
c
c  add the buoiancy term  
c
      densit = 0.5*(dens(ic,jc,kc)+dens(ic,jc,km))*sinpc(kc)*float(ibuo)

      qcap(ic,jc,kc)=d11q1v+d11q2v-hq3+coriol
     %               +densit*(+cospc(kc)*cos(thetam(ic))) !gravita'' allimeata in x
c    %               +densit*(+sinpc(kc))         !gravita'' verticale
c
c     add the randomforcing term
c
      if(iforc.eq.1) then
        om_tape = nint(om_th*sin(zz(kc)))
        forc = forc_mod * sin(zz(kc))
     %        * sin(om_tape * thetam(ic) + ph_th)
     %        * (rm(jc)-r_min)/(r_max-r_min)
     %        * cos(om_z * zz(kc) + ph_z)
        qcap(ic,jc,kc) = qcap(ic,jc,kc) + forc
      end if
     
c
   30 continue
!$OMP  END PARALLEL DO
c
      return
      end
c
