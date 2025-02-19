************************************************************************
      subroutine gcurv
      include'param.f'
c
c     Code for the computation of three-dimensional incompressible flows
c     in cylindrical polar coordinates.
c
c     This code solves flow fields bounded in the x3 (axial) and
c     x2 (radial) directions. All boundaries can be no-slip or free-slip
c     by setting the appropriate indices in the input file.
c     The geometry  (a cylindrical can) is given in the subroutine cordi.
c     The discretization is uniform in the axial (3) and azimuthal (1)
c     directions, while can be non-uniform in the radial direction.
c
c     The equations for the following variables
c
c      q1=v(theta)    q2=v(r)*r     q3=v(zeta)
c
c     are discretized by finite-difference schemes.
c     The introduction of the variable q2 is necessary to avoid the
c     problem of the singularity of the equations at the axis of symmetry
c     (r=0).
c
c     All spatial derivatives are discretized by central second-order
c     accurate finite-difference schemes including the non linear terms.
c
c     The nonlinear terms are treated explicitly, while the viscous terms
c     are computed implicitly. This would lead to the inversion of a
c     large banded matrix, that however is avoided by introducing
c     a factored scheme bringing to the solution of three tridiagonal
c     matrices for each velocity component (subroutine INVTR*).
c
c     In time a fractional-step procedure is used in the version of
c     Nagi Mansour introducing the pressure in the first step.
c
c     The non-linear terms and the cross derivatives of the viscous terms
c     are discretized by explicit  Adams-Bashfort or 3rd order Runge-Kutta
c     method (A. Wray, personal communication).
c
c     The scalar quantity Phi, which projects the provisional velocity field
c     onto  a divergence free field, is solved by a direct method.
c     For the axial and azimuthal directions modified wave numbers coupled
c     with trigonometric expansions (FFTs) are used. The equation is then
c     solved by simply inverting a tridiagonal matrix for the radial direction.
c     No explicit boundary conditions are necessary for this Poisson equation.
c
c     Other details of the scheme are given in the introduction of the
c     subroutine TSCHEM
c
c     timings
c     real timef
      real init_time(2)         ! initialization time
      real iter_time(2)         ! iteration time
      character*25 namfi
      dimension xxx(m1,m2,m3),yyy(m1,m2,m3),zzz(m1,m2,m3)
c
      pi=2.*asin(1.)
c
c
c     init_time(1)=timef()
c
      call initia
C**************************
      imovie = 0
      tframe = 0.05
C**************************
c     fvisc = 1.000 
c     fvisc1 = 1.0000 
c     fvisc2 = 1.00
c     fvisc3 = 0.0000 
c
c     grid information
c
      call meshes
      call indic
      call cordin
      write(6,754)n1,n2,n3
      write(32,754)n1,n2,n3
  754 format(/,5x,'numero di punti griglia: ',' n1= ',i3,' n2= ',i3,
     % ' n3= ',i3/)
      if(iaxsy.eq.0) write(6,856)
      if(iaxsy.eq.0) write(32,856)
  856 format(/,5x,'calcolo con condizione di assialsimmetria   '
     % ,i3/)
      write(6,755) dx1,dx2,dx3,dt,ntst
      write(32,755) dx1,dx2,dx3,dt,ntst
  755 format(/,2x,'dx1=',e10.3,' dx2=',e10.3,' dx3=',e10.3,'  dt=',e10.3,
     % '  ntst=',i5,/)
      time=0.
      ntii=0
      if(nson.ne.0) then
        call indson
        do i = 1,nson
          write(6,123) i,(mason(i,j),j=1,3)
          write(32,123) i,(mason(i,j),j=1,3)
        end do
  123   format(2x,'sonda n= ',i3,' : i=',i3,', j=',i3,', k=',i3)
        do i = 1,nson
          write(6,923) i,(coson(i,j),j=1,3)
          write(32,923) i,(coson(i,j),j=1,3)
        end do
      endif
  923 format(2x,'sonda n= ',i3,' : th=',f8.2,', r=',f8.2,', x=',f8.2)
      write(6,*)
      write(32,*)
      write(6,*) ' AAA '
C    PRINT 3D coordinate file
      write(*,*)'-------- Grid file -----------'
      do kc=1,n3
        do jc=1,n2
          do ii=1,n1
c           ic = ii
c           if(ii.eq.n1) ic=1
c           xxx(ic,jc,kc) = rm(jc)*cos(thetam(ic))*sin(zm(kc))
c           yyy(ic,jc,kc) = rm(jc)*sin(thetam(ic))*sin(zm(kc))
c           zzz(ic,jc,kc) = rm(jc)*cos(zm(kc))
c           ic = ii
c           if(ii.eq.n1) ic=1
            xxx(ii,jc,kc) = rc(jc)*cos(thetac(ii))*sin(zz(kc))
            yyy(ii,jc,kc) = rc(jc)*sin(thetac(ii))*sin(zz(kc))
            zzz(ii,jc,kc) = rc(jc)*cos(zz(kc))
          end do
        end do
      end do

      namfi='grid.xyz'

      open(99,file=namfi,form='unformatted')
      write(99) 1
      write(99) n1,n2,n3
      write(99)
     %   ((((xxx(i,j,k)),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   ((((yyy(i,j,k)),
     %   i=1,n1),j=1,n2),k=1,n3),
     %   ((((zzz(i,j,k)),
     %   i=1,n1),j=1,n2),k=1,n3)
      close(99)
c
c   read or create initial fields
c
      call phini
      write(6,*) ' AAA2 '
c
      do 22 l=1,ndv
        vmax(l)=0.
   22 continue
c
      call densbo
c
c
c      create the initial conditions
c
      if(nread.eq.0) then
        write(6,*)' nread=0 ---> cond. iniz. create nel code'
        write(32,*)' nread=0 ---> cond. iniz. create nel code'
        ntii=0
        ntime=0
        time=0.
        cflm=0.
**** SCRIVO CONDIZIONI INIZIALI***

c     indexForInitial = 4
c index = 1 -> ball
c index = 2 -> Hill vortex
c index = 3 -> null (with random pert.)
c index = 4 -> homogeneous velocity
        if ( indexForInitial .eq. 1 ) then
          call inqpr
        elseif ( indexForInitial .eq. 2 ) then
          call inqHill
        elseif ( indexForInitial .eq. 3 ) then
          call inqRB
        elseif ( indexForInitial .eq. 4 ) then
          call inqHR
        end if
          write(6,*) ' AAA3 '
        call divgck(dmax,dtot)
          write(6,*) ' AAA4 '
        write(6,*)' initial divg dmax,dtot  ',dmax,dtot
        write(32,*)' initial divg dmax,dtot  ',dmax,dtot
        call outpf(time)
      else
        write(6,*)' nread=1 ---> cond. iniz. lette da file'
        write(32,*)' nread=1 ---> cond. iniz. lette da file'
        call inirea(ntii,time)
cc      call divgck(dmax,dtot)
        write(6,*)' initial divg dmax,dtot  ',dmax,dtot
        write(32,*)' initial divg dmax,dtot  ',dmax,dtot
c       call cfl(cflm)
c       call vmaxv
c       call outh(time,dt)
c       if(iaxsy.ne.0) call enerc(time)
c       if(nson.ne.0) call wson(time)
        call init_spec
        write(6,*)' spectrum computation initiated'
        write(32,*)' spectrum computation initiated'
      endif
c
      ntstf=ntii+ntst
      ntii=ntii+1
      write(6,711)tprint,ntii,ntstf,tpin
      write(32,711)tprint,ntii,ntstf,tpin
  711 format(3x,'check in cond : tprint =',f8.4,'  ntii =',i8,
     1       '  ntstf =',i8,2x,'tpin =',f8.4//)
      if(idtv.eq.1) then
        write(6,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),
     %  jmaxv(1),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),
     %  imaxv(3),jmaxv(3),kmaxv(3),dt,qmax,iqm,jqm,kqm,densm,densmax,
     %  densmin'
        write(32,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),
     %  jmaxv(1),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),
     %  imaxv(3),jmaxv(3),kmaxv(3),dt,qmax,iqm,jqm,kqm,densm,densmax,
     %  densmin'
c       call outh(time,dt)
c       if(iaxsy.ne.0) call enerc(time)
c       if(nson.ne.0) call wson(time)
      else
        write(6,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),
     %  jmaxv(1),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),
     %  imaxv(3),jmaxv(3),kmaxv(3),cflm,qmax,iqm,jqm,kqm,densm,densmax,
     %  densmin'
        write(32,*)'ntime,time,voamax,vormax,vozmax,vmax(1),imax(i),
     %  jmaxv(1),kmaxv(1),vmax(2),imaxv(2),jmaxv(2),kmaxv(2),vmax(3),
     %  imaxv(3),jmaxv(3),kmaxv(3),cflm,qmax,iqm,jqm,kqm,densm,densmax,
     %  densmin'
        cflm=cflm*dt
c       call outh(time,cflm)
c       if(iaxsy.ne.0) call enerc(time)
c       if(nson.ne.0) call wson(time)
      endif
c
      call coetar
      write(6,*) ' AAA5 '
c
c     init_time(2)=timef()
      write(6,*) '       Init Time = '
     %             ,(init_time(2)-init_time(1))*0.001
c


c  ********* starts the time dependent calculation ***
c

c     initialize the sonda files
      if(nson.ne.0) call wson(time)

      do 350 ntime=ntii,ntstf
c
c     the calculation stops if the velocities are diverging for numerical
c     stability conditions (courant number restrictions)
c
c           iter_time(1)=timef()
c
        call cfl(cflm)
        if(idtv.eq.1) then
          if(ntime.ne.1) then
            dt=cflmax/cflm
            if(dt.gt.dtmax) dt=dtmax
          endif
          if(dt.lt.1E-16) go to 166
          else
            cflm=cflm*dt
          if(cflm.gt.cfllim) go to 165
        endif
        beta=dt/ren*0.5
        call tschem(time)
c       call outpf(time)
c       pause
        time=time+dt
        if(ntime.eq.1) go to 306
        if(nson.ne.0) then
          if(mod(time,tpin/4.).lt.dt) call wson(time)
        endif
        if(mod(time,tpin).lt.dt) go to 306
        go to 305
  306   continue
        call vmaxv
        if(vmax(1).gt.1000.and.vmax(2).gt.1000) go to 266
        call cfl(cflm)
        call divgck(dmax,dtot)
        if(dtot .ne. dtot) go to 366
        call densmc
        if(idtv.eq.1) then
          call outh(time,dt)
          call diss(time)
CC        if(iaxsy.ne.0) call enerc(time)
CRV       if(nson.ne.0) call wson(time)
        else
          cflm=cflm*dt
          call outh(time,cflm)
          call diss(time)
CC        if(iaxsy.ne.0) call enerc(time)
CRV       if(nson.ne.0) call wson(time)
        endif
        if(nread.eq.1) call spect
        if(dmax.gt.resid) go to 169
c
  305   continue
c
c     write the flow field
c
        if(mod(time,tprint/2.).lt.dt) then
          if(nwrit.eq.1) then
            call continua(time)
          endif
        endif
        if(mod(time,tprint).lt.dt) then
          call outpf(time)
          if(nread.eq.1) call print_spec(time)
        endif
        if(time.gt.tmax) go to 167
c
c       iter_time(2)=timef()
c       write(6,*) '  Iteration Time = '
c    %             ,(iter_time(2)-iter_time(1))*0.001
c
  350 continue
c
c     write the restart file anyway before
c     exiting
c
      call continua(time)
      go to 167
  165 continue
      write(6,164)
      write(32,164)
  164 format(10x,'too large CFL')
      go to 167
  166 continue
      write(6,168) dt
      write(32,168) dt
  168 format(10x,'too small dt = ',e14.7)
      go to 167
  266 continue
      write(6,268)
      write(32,268)
  268 format(10x,'velocities diverged')
      go to 167
  366 continue
      write(6,368)
      write(32,368)
  368 format(10x,'NaN appeared')
      go to 167
  169 continue
      write(6,178) dmax,imxq,jmxq,kmxq
      write(32,178) dmax,imxq,jmxq,kmxq
  178 format(10x,'too large local residue for mass conservation : '
     1       ,e12.5,' at ',3(1x,i3) )
  167 continue
c
c     iter_time(2)=timef()
      write(6,*) '  Iteration Time = '
     %             ,(iter_time(2)-iter_time(1))*0.001
c
      return
      end
c
************************************************************************
c
c     in this subroutine the indices ip,im,jp,jm,kp,km are calculated.
c     the indices jpc,jup,jmc,jum,kpc,kup,kmc,kum are indices
c     for the walls
c
      subroutine indic
      include'param.f'
c
c
c   azimuthal direction
c
      do 1 ic=1,n1m
        imv(ic)=ic-1
        if(ic.eq.1) imv(ic)=n1m
        ipv(ic)=ic+1
        if(ic.eq.n1m) ipv(ic)=1
    1 continue
c
c    symmetrical index for the azimuthal location
c    i + pi
c
c     write(6,*) ' I       ISYM   '
      do i=1,n1m
       isym(i) = i + n1m/2
       if(isym(i).gt.n1m) isym(i) = isym(i) - n1m
c      write(6,*) i,isym(i)
      enddo
c     pause
c
c   streamwise direction
c
      do 4 kc=1,n3m
        kmv(kc)=kc-1
        kpv(kc)=kc+1
        if(kc.eq.1) kmv(kc)=kc
        if(kc.eq.n3m) kpv(kc)=kc
    4 continue
c
c    symmetrical index for the latitude location
c     pi - k
c
c     write(6,*) ' K       KSYM   KSYC '
      do k=1,n3m
       ksym(k) = n3m - k + 1
       ksyc(k) = n3 - k + 1
c     write(6,*) k,ksym(k), ksyc(k)
      enddo
       k=n3
       ksyc(k) = n3 - k + 1
c     write(6,*) k,' -- ', ksyc(k)
c     pause

c
c     direction normal to non-slip walls
c
      do 150 kc=1,n3m
        kpc(kc)=kpv(kc)-kc
        kmc(kc)=kc-kmv(kc)
        kup(kc)=1-kpc(kc)
        kum(kc)=1-kmc(kc)
  150 continue
c
c   indices for radial direction
c
      do 3 jc=1,n2m
        jmv(jc)=jc-1
        jpv(jc)=jc+1
        if(jc.eq.1) jmv(jc)=jc
        if(jc.eq.n2m) jpv(jc)=jc
    3 continue
c
c   indices for radial direction
c
      do 15 jc=1,n2m
        jpc(jc)=jpv(jc)-jc
        jmc(jc)=jc-jmv(jc)
        jup(jc)=1-jpc(jc)
        jum(jc)=1-jmc(jc)
   15 continue
c
      return
      end
c
************************************************************************
c
c     This subroutine computes the div(q^).
c     q^ is the non-solenoidal field containing the modified variables
c     q1=v(theta)    q2=v(r)*r^2     q3=v(phi)*sin phi
c     Note also that instead of computing the raw divergence it is computed
c     the divergence multiplied by r^2 for numerical/computational reasons.
c
c     d q_r      r     d q_p      r    d q_t
c     ----- + -------- ----- + ------- -----
c     d r     sin phi  d p     sin phi  d t
c
      subroutine divg
      include'param.f'
c
      usdtal = 1.d0/(dt*al)
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,ipv,usrm,rm,usinpm)
!$OMP$ SHARED(q1,q2,q3,dph,dx1,udx2m,udx3m,usdtal)
!$OMP$ PRIVATE(kp,jp,ip,ic,jc,kc,fact)
      do 11 kc=1,n3m
        kp=kc+1
        do 11 jc=1,n2m
          jp=jc+1
          fact = rm(jc)*usinpm(kc)
          do 11 ic=1,n1m
            ip=ipv(ic)
            dph(ic,jc,kc)=(
     %             (q2(ic,jp,kc)-q2(ic,jc,kc))*udx2m(jc)
     %            +((q1(ip,jc,kc)-q1(ic,jc,kc))*dx1
     %            + (q3(ic,jc,kp)-q3(ic,jc,kc))*udx3m(kc))*fact
     %                     )*usdtal
   11 continue
!$OMP  END PARALLEL DO
      return
      end
c
c
c************************************************************************
c
c     This subroutine projects q^ onto a solenoidal velocity field.
c
c       q^(n+1)=q^-grad(dph)*dt*al
c
      subroutine updvp
      include'param.f'
c
c  ***********  compute the q1 velocity component
c               v1dgf=component 1 of grad(dph)
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,imv)
!$OMP$ SHARED(q1,dph,usrm,al,dt,dx1,usinpm)
!$OMP$ PRIVATE(im,usurm,ic,jc,kc)
      do 1 kc=1,n3m
        do 1 jc=1,n2m
          usurm = al*dt*dx1*usrm(jc)*usinpm(kc)
          do 1 ic=1,n1m
            im=imv(ic)
            q1(ic,jc,kc)=q1(ic,jc,kc)-
     %      (dph(ic,jc,kc)-dph(im,jc,kc))*usurm
    1 continue
!$OMP END PARALLEL DO
c
c
c  ***********  compute the q2 velocity component
c               v2dgf=component 2 of grad(dph)
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,udx2c)
!$OMP$ SHARED(q2,dph,rc,al,dt,n2)
!$OMP$ PRIVATE(jm,usurm,ic,jc,kc)
      do 2 kc=1,n3m
        do 2 jc=2,n2m
          jm=jc-1
          usurm = al*dt*rc(jc)**2*udx2c(jc)
          do 2 ic=1,n1m
            q2(ic,jc,kc)=q2(ic,jc,kc)-
     %      (dph(ic,jc,kc)-dph(ic,jm,kc))*usurm
    2 continue
!$OMP END PARALLEL DO
      do 21 kc=1,n3m
        do 21 ic=1,n1m
          q2(ic,1,kc)=0.
          q2(ic,n2,kc)=0.
   21 continue
c
c  ***********  compute the q3 velocity component
c               q3 is the cartesian component
c               v3dgf=component 3 of grad(dph)
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,imv,udx3c)
!$OMP$ SHARED(q3,dph,al,dt,usrm,sinpc)
!$OMP$ PRIVATE(km,usurm,ic,jc,kc)
      do 5 kc=2,n3m
        km=kc-1
        do 5 jc=1,n2m
          usurm = al*dt*udx3c(kc)*usrm(jc)*sinpc(kc)
          do 5 ic=1,n1m
            q3(ic,jc,kc)=q3(ic,jc,kc)-
     %      (dph(ic,jc,kc)-dph(ic,jc,km))*usurm
    5 continue
!$OMP END PARALLEL DO
      do 210 jc=1,n2m
        do 210 ic=1,n1m
          q3(ic,jc,1)=0.
          q3(ic,jc,n3)=0.
  210 continue
      return
      end
c
c
c
c
************************************************************************
c   this subroutine performs the calculation of the pressure.
c   this depends on the fractional step
c
      subroutine prcalc
      include'param.f'
c
c     the pressure is evaluated at the center of the box.
c
c     p^{n+1} = p^{n} + phi^{n+1} - b * Nabla^2 phi^{n+1}
c
      be=al*beta
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,kpv,kmv,jmv,jpv,imv,ipv)
!$OMP$ SHARED(pr,dph,be,dx1q,usinpm)
!$OMP$ SHARED(apphk,acphk,amphk,apphj,acphj,amphj)
!$OMP$ PRIVATE(kp,km,jm,jp,im,ip,ic,jc,kc)
      do 1 kc=1,n3m
        kp=kpv(kc)
        km=kmv(kc)
        do 1 jc=1,n2m
          jm=jmv(jc)
          jp=jpv(jc)
            do 1 ic=1,n1m
              ip=ipv(ic)
              im=imv(ic)
              pr(ic,jc,kc)=pr(ic,jc,kc)+dph(ic,jc,kc)
     %                  -be*(
     %                       (dph(ip,jc,kc)-2.*dph(ic,jc,kc)
     %                       +dph(im,jc,kc))*dx1q*usinpm(kc)**2
     %                      +(dph(ic,jc,kp)*apphk(kc)
     %                       +dph(ic,jc,kc)*acphk(1,kc)
     %                       +dph(ic,jc,km)*amphk(kc))
     %                      +(dph(ic,jp,kc)*apphj(jc)
     %                       +dph(ic,jc,kc)*acphj(jc)
     %                       +dph(ic,jm,kc)*amphj(jc))
     %                      )*usrm(jc)**2 !*float(ivisc)
    1 continue
!$OMP  END PARALLEL DO
c      if(rc(1).le.1.e-3) then
c      jc=1
c      pmcen = 0.
c      do kc=1,n3m
c       do ic=1,n1m
c        pmcen = pmcen + pr(ic,jc,kc)
c       end do
c      end do
c      pmcen = pmcen/float(n1m*n3m)
c      do kc=1,n3m
c       do ic=1,n1m
c        pr(ic,jc,kc) = pmcen
c       end do
c      end do
c      end if
c
      pmed = pr(1,n2m,n3m)
c     pmed = 0.
!$OMP  PARALLEL DO
!$OMP$ SHARED(pr,pmed,n1m,n2m,n3m)
!$OMP$ PRIVATE(i,j,k)
      do i=1,n1m
        do j=1,n2m
          do k=1,n3m
            pr(i,j,k) = pr(i,j,k) - pmed
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c         k=1
c         do j=1,n2m
c           do i=1,n1m/2
c            is=isym(i)
c            app = 0.5*(pr(i,j,k)+pr(is,j,k))
c            pr(i,j,k)=app
c            pr(is,j,k)=app
c           end do
c         end do
c         k=n3m
c         do j=1,n2m
c           do i=1,n1m/2
c            is=isym(i)
c            app = 0.5*(pr(i,j,k)+pr(is,j,k))
c            pr(i,j,k)=app
c            pr(is,j,k)=app
c           end do
c         end do
c      kc=1
c      do jc=1,n2m
c       pmeax = 0.
c       do ic=1,n1m
c        pmeax = pmeax + pr(ic,jc,kc)
c       end do
c       pmeax = pmeax/float(n1m)
c       do ic=1,n1m
c        pr(ic,jc,kc) = pmeax
c       end do
c      end do
c      kc=n3m
c      do jc=1,n2m
c       pmeax = 0.
c       do ic=1,n1m
c        pmeax = pmeax + pr(ic,jc,kc)
c       end do
c       pmeax = pmeax/float(n1m)
c       do ic=1,n1m
c        pr(ic,jc,kc) = pmeax
c       end do
c      end do
      return
      end
c
c
c***********************************************************************
c                                                                      *
c                       CONDIZIONI INIZIALI                            *
c                                                                      *
c***********************************************************************
      subroutine inqpr
      include'param.f'
      pi=2.*asin(1.)
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,q1,q2,q3)
!$OMP$ PRIVATE(i,j,k)
      do 80 k=1,n3
         do 80 j=1,n2
            do 80 i=1,n1
              q1(i,j,k)=0.
              q2(i,j,k)=0.
              q3(i,j,k)=0.
  80  continue
!$OMP  END PARALLEL DO
c
      pi = 2.*asin(1.)
      r0 = 0.5
      th0 = 0.*pi/2.
      ph0 = 0. !pi/2.
      x0 = r0*cos(th0)*sin(ph0)
      y0 = r0*sin(th0)*sin(ph0)
      z0 = r0*cos(ph0)

      icon = 0
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1m,q1,q2,q3,rc,thetac,thetam,zz,zm,x0,y0,z0)
!$OMP$ PRIVATE(i,j,k,xt,yt,zt,rho,vx1,vx2,vx3)
      do 810 k=1,n3
        do 810 j=1,n2
          do 810 i=1,n1m
            xt = rc(j)*cos(thetac(i))*sin(zz(k))
            yt = rc(j)*sin(thetac(i))*sin(zz(k))
            zt = rc(j)*cos(zz(k))
            rho = sqrt((xt-x0)**2+(yt-y0)**2+(zt-z0)**2)
            if(rho.le.0.25) then
              vx1 =  -1.0!sqrt(2.)/2.! 0.0!*(1.-rm(j))*sin(zm(k))
              vx2 =  0.0!*(1.-rc(j))*sin(zm(k))
              vx3 =  0.0 !sqrt(2.)/2.! 1.0!*(1.-rm(j))*sin(zz(k))
              q1(i,j,k)=(vx2*cos(thetac(i))-vx1*sin(thetac(i)))
              q2(i,j,k)=((vx1*cos(thetam(i))+vx2*sin(thetam(i)))
     %                   *sin(zm(k))+vx3*cos(zm(k)))*rc(j)**2
              q3(i,j,k)=((vx1*cos(thetam(i))+vx2*sin(thetam(i)))
     %                   *cos(zz(k))-vx3*sin(zz(k)))*sin(zz(k))
            end if
c     if(rm(j).le.0.35) then
c     icon = icon+1
c     write(889,1888) icon,i,j,k,thetac(i),rm(j),zm(k),
c    %                q1(i,j,k),q2(i,j,k),q3(i,j,k)
c     end if
c1888 format(i5,3(1x,i3),8(1x,e12.5))
C             q3(i,j,k)=1.0*(
C    %       sin(2.*thetam(i))+sin(3.*thetam(i)))*sin(zz(k))*
C    %        1.414*(rm(j)*(1.-rm(j)/rc(n2))/rc(n2))**0.25
C             q2(i,j,k)=0.
C             q1(i,j,k)=1.0*sin(5.*zm(k))*
C    %        (sin(1.*thetac(i))+sin(4.*thetac(i)))*
C    %        1.414*(rm(j)*(1.-rm(j)/rc(n2))/rc(n2))**0.25
c    %        *rm(j)*(1.-rm(j)/rc(n2))/rc(n2)
c             q1(i,j,k)=0.1*cos(5.*thetac(i))*sin(2.*zm(k))
  810 continue
!$OMP  END PARALLEL DO
c      close(889)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      open(69,file='CXY0.dat',form='formatted',status='unknown')
      rewind(69)
      write(69,*)'Variables = X, Y, Z, Vx, Vy, Vz, Pr'
      write(69,*)'ZONE F=POINT, I= ',n1,', J= ',n2m, ', K= ',1
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
         isy = isym(ic)
         ksy = ksym(kc)
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
          vp = 0.5*(q3(ic,jc,kc)+q3(ic,jc,kp))/sin(zm(kc))
          ct = cos(thetam(ic))
          st = sin(thetam(ic))
          cp = cos(zm(kc))
          sp = sin(zm(kc))
                      write(69,*)xxx,yyy,zzz
     %           ,vr*ct*sp-vt*rm(jc)*st*sp+vp*rm(jc)*ct*cp
     %           ,vr*st*sp+vt*rm(jc)*ct*sp+vp*rm(jc)*st*cp
     %           ,vr*cp-vp*rm(jc)*sp
     %           ,pr(ic,jc,kc)
            end do
          end do
      close(69)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c      call divgck(qmax,qtot)
c
cc
      return
      end
c
c
c***********************************************************************
c                                                                      *
c                CONDIZIONI INIZIALI: Rayleigh–Bénard                  *
c                                                                      *
c***********************************************************************
      subroutine inqRB
      include'param.f'
      pi=2.*asin(1.)
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,q1,q2,q3,dens)
!$OMP$ PRIVATE(i,j,k)
      do 80 k=1,n3
         do 80 j=1,n2
            do 80 i=1,n1
              q1(i,j,k)=0.
              q2(i,j,k)=0.
              q3(i,j,k)=0.
              dens(i,j,k)=0.5 + 0.05*rand()
  80  continue
!$OMP  END PARALLEL DO
c
      return
      end
c
c
c***********************************************************************
c                                                                      *
c            CONDIZIONI INIZIALI: Homogeneous Rotation                 *
c                                                                      *
c***********************************************************************
      subroutine inqHR
      include'param.f'
      pi=2.*asin(1.)
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,q1,q2,q3,dens)
!$OMP$ PRIVATE(i,j,k)
      do 85 k=1,n3
        do 85 j=1,n2
          do 85 i=1,n1
            q1(i,j,k)=0.01
            q2(i,j,k)=0.
            q3(i,j,k)=0.
            dens(i,j,k)=0.5
  85  continue
!$OMP  END PARALLEL DO
      do 86 k=1,n3
        do 86 j=1,n2
          q1(1,j,k)=0.
          q1(n1,j,k)=0.
  86  continue  
c
      return
      end
c
c************************************************************************
c                                                                       *
c ****************************** subrout coetar  ********************** *
c                                                                       *
c    this subroutine calculates the coefficients for the                *
c    integration in the radial direction with non-uniform coor. trasf.  *
c                                                                       *
c************************************************************************
      subroutine coetar
      include 'param.f'
c
c  ***********  coefficients for q1
c
c    1    d  r^2  d  q_t
c   ---  ---      ----------
c   r^2   dr       dr
c
      do 151 jc=2,n2m-1
      jp=jc+1
      jm=jc-1
      a22=dx2q/(g2rm(jc)*rm(jc)**2)
      a22p= +a22*rc(jp)**2/g2rc(jp)
      a22m= +a22*rc(jc)**2/g2rc(jc)
      ap1j(jc)=a22p
      am1j(jc)=a22m
      ac1j(jc)=-(a22p+a22m)
  151 continue
c
c    internal sphere bound. conditions
c
      jc=1
      jp=jc+1
      a22=dx2q/(g2rm(jc)*rm(jc)**2)
      a22p= +a22*rc(jp)**2/g2rc(jp)
      a22m= +a22*rc(jc)**2/g2rc(jc)
       if(inslwri.eq.1) then        ! free-slip
        ap1j(jc)=a22p
        ac1j(jc)=-a22p
        am1j(jc)=0.
       else                         ! no-slip
        ap1j(jc)= a22p
        am1j(jc)=0.
        ac1j(jc)=-(a22p+2.*a22m)
       end if
c
c    external sphere bound. conditions
c
      jc=n2m
      jp=jc+1
      a22=dx2q/(g2rm(jc)*rm(jc)**2)
      a22p= +a22*rc(jp)**2/g2rc(jp)
      a22m= +a22*rc(jc)**2/g2rc(jc)
       if(inslwr.eq.1) then         ! free-slip
        ap1j(jc)=0.
        ac1j(jc)=-a22m
        am1j(jc)=a22m
       else                         ! no-slip
        ap1j(jc)= 0.
        am1j(jc)=a22m
        ac1j(jc)=-(2.*a22p+a22m)
       end if
c
c         1      d   sin phi  d  q_t            q_t
c   -----------  ---         ----------  - ------------
c   r^2 sin phi  dp           dp           r^2 sin^2 phi
c
      do 152 kc=1,n3m
      kp=kc+1
      do 152 jc=1,n2m
      a33=dx3q/(g3rm(kc)*rm(jc)**2*sinpm(kc))
      a33p= +a33*sinpc(kp)/g3rc(kp)
      a33m= +a33*sinpc(kc)/g3rc(kc)
      ap1k(jc,kc)=a33p
      am1k(jc,kc)=a33m
      ac1k(jc,kc)=-(a33p+a33m) -usrmq(jc)*usinpm(kc)**2
      aprk(jc,kc)=a33p
      amrk(jc,kc)=a33m
      acrk(jc,kc)=-(a33p+a33m) 
  152 continue
c
c  ***********  coefficients for q2
c
c
c     d  r^2  d  q_r/r^2       2 q_r
c    ---      ----------    - ------
c     dr       dr               r^2
c
c
c WARNING: only for jc=2 and only for the implicit term, in order to
c          avoid division by rc(1)=0 the above derivative is expanded
c          as:
c              d^2 q_r   -2 d q_r
c              ------     - -----
c              d r^2      r  d r
c
c
c
      am2j(1)=0.
      ap2j(1)=0.
      ac2j(1)=1.
      am2j(n2)=0.
      ap2j(n2)=0.
      ac2j(n2)=1.
      do 251 jc=2,n2m
      jp=jc+1
      jm=jc-1
      a22=dx2q/g2rc(jc)
      a22p= +a22*rm(jc)**2/(g2rm(jc)*rc(jp)**2)
      a22m= +a22*rm(jm)**2/(g2rm(jm)*rc(jm)**2)
      ap2j(jc)=a22p
      am2j(jc)=a22m
      ac2j(jc)=-(a22*(rm(jc)**2/g2rm(jc)+rm(jm)**2/g2rm(jm))+2.)
     %                   /rc(jc)**2
  251 continue
c
c         1      d   sin phi  d  q_r
c   -----------  ---         ----------
c   r^2 sin phi  dp           dp
c
      do 252 kc=1,n3m
      ap2k(1,kc)=0.
      am2k(1,kc)=0.
      ac2k(1,kc)=1.
      ap2k(n2,kc)=0.
      am2k(n2,kc)=0.
      ac2k(n2,kc)=1.
      kp=kc+1
      do 252 jc=2,n2m
      a33=dx3q/(g3rm(kc)*rc(jc)**2*sinpm(kc))
      a33p= +a33*sinpc(kp)/g3rc(kp)
      a33m= +a33*sinpc(kc)/g3rc(kc)
      ap2k(jc,kc)=a33p
      am2k(jc,kc)=a33m
      ac2k(jc,kc)=-(a33p+a33m)
  252 continue
c
c  ***********  coefficients for q3
c
c    1    d  r^2  d  q_p
c   ---  ---      ----------
c   r^2   dr       dr
c
      do 351 jc=2,n2m-1
      jp=jc+1
      a22=dx2q/(g2rm(jc)*rm(jc)**2)
      a22p= +a22*rc(jp)**2/g2rc(jp)
      a22m= +a22*rc(jc)**2/g2rc(jc)
      ap3j(jc)=a22p
      am3j(jc)=a22m
      ac3j(jc)=-(a22p+a22m)
  351 continue
c
c    internal sphere bound. conditions
c
      jc=1
      jp=jc+1
      a22=dx2q/(g2rm(jc)*rm(jc)**2)
      a22p= +a22*rc(jp)**2/g2rc(jp)
      a22m= +a22*rc(jc)**2/g2rc(jc)
       if(inslwri.eq.1) then
        ap3j(jc)=a22p
        ac3j(jc)=-a22p
        am3j(jc)=0.
       else
        ap3j(jc)= a22p
        am3j(jc)= 0.
        ac3j(jc)=-(a22p+2.*a22m)
       end if
c
c    external sphere bound. conditions
c
      jc=n2m
      jp=jc+1
      a22=dx2q/(g2rm(jc)*rm(jc)**2)
      a22p= +a22*rc(jp)**2/g2rc(jp)
      a22m= +a22*rc(jc)**2/g2rc(jc)
       if(inslwr.eq.1) then
        ap3j(jc)=0.
        ac3j(jc)=-a22m
        am3j(jc)=a22m
       else
        ap3j(jc)= 0.
        am3j(jc)=a22m
        ac3j(jc)=-(2.*a22p+a22m)
       end if
c
c    1    d  (sin phi  d  q_p/sin(p) )          q_p
c   ---  --- (       --------------- )    - ------------
c   r^2   dp (         dp            )      r^2 sin^2 phi
c
c
c
c WARNING: only for kc=2,n3m and only for the implicit term, in order to
c          avoid division by sin(kc)=0 the above derivative is expanded
c          as:
c        1  (  d^2 q_p    1   d q_p )
c      ---- (  ------  - ---- ----- )
c       r^2 (  d p^2     tg(p) d p  )
c
c
c
c
c
      do 352 kc=3,n3m-1
      kp=kc+1
      km=kc-1
      do 352 jc=1,n2m
      a33=dx3q/(g3rc(kc)*rm(jc)**2)
      a33p= +a33*sinpm(kc)/g3rm(kc)
      a33m= +a33*sinpm(km)/g3rm(km)
      ap3k(jc,kc)=a33p/sinpc(kp)
      am3k(jc,kc)=a33m/sinpc(km)
      ac3k(jc,kc)=-(a33p+a33m)/sinpc(kc)-usrmq(jc)*usinpc(kc)**2
c     a33=dx3/(g3rc(kc)*rm(jc)**2)
c     ap3k(jc,kc)=a33*(dx3/g3rm(kc)-utanpc(kc)/2.)
c     am3k(jc,kc)=a33*(dx3/g3rm(km)+utanpc(kc)/2.)
c     ac3k(jc,kc)=-a33*dx3*(1./g3rm(kc)+1./g3rm(km))
  352 continue
      kc=2
      kp=kc+1
      km=kc-1
      do jc=1,n2m
      a33=dx3/(g3rc(kc)*rm(jc)**2)
      ap3k(jc,kc)=a33*(dx3/g3rm(kc)-utanpc(kc)/2.)
      am3k(jc,kc)=a33*(dx3/g3rm(km)+utanpc(kc)/2.)
      ac3k(jc,kc)=-a33*dx3*(1./g3rm(kc)+1./g3rm(km))
      end do
      kc=n3m
      kp=kc+1
      km=kc-1
      do jc=1,n2m
      a33=dx3/(g3rc(kc)*rm(jc)**2)
      ap3k(jc,kc)=a33*(dx3/g3rm(kc)-utanpc(kc)/2.)
      am3k(jc,kc)=a33*(dx3/g3rm(km)+utanpc(kc)/2.)
      ac3k(jc,kc)=-a33*dx3*(1./g3rm(kc)+1./g3rm(km))
      end do
      do 353 jc=1,n2m
      ap3k(jc,1)=0.
      am3k(jc,1)=0.
      ac3k(jc,1)=1.
      ap3k(jc,n3)=0.
      am3k(jc,n3)=0.
      ac3k(jc,n3)=1.
  353 continue
c
C
C     write(6,*) ' centred coeff. K '
C     do k=1,n3
C     write(6,*) k, am3ck(k), ac3ck(k), ap3ck(k)
C     end do
C     pause
C     write(6,*) ' staggered coeff. K '
C     do k=1,n3m
C     write(6,*) k, am3sk(k), ac3sk(k), ap3sk(k)
C     end do
C     pause
C     write(6,*) ' staggered coeff. K (2)'
C     do k=1,n3m
C     write(6,*) k, am3ssk(k), ac3ssk(k), ap3ssk(k)
C     end do
C     pause
      return
      end
c
c
************************************************************************
c                       SUBROUTINE INVTRRO
c   This subroutine performs the computation of Q~~ for the Temperature
c   equation  by a factored implicit scheme.
c   Viscous terms are treated implicitly, nonlinear terms explicitly.
c
c        ~~     n
c   dQ = Q  -  Q
c
c         alp dt
c   bet = -------
c          2 Re
c
c   The first equation in the introduction of THSCHEM then becomes
c
c
c              2      [            n         n       n-1   alp       2  n ]
c  (1-bet*nabla ) dQ =[-alp*grad (P ) + gam*H + rho*H   + ----- nabla (Q )]*dt
c                     [                                    Re             ]
c
c                                                                      3
c   The left hand side of this equation is then factored at an order dt
c   as follows
c
c           2           2           2
c  (1 -bet*d_x)*(1-bet*d_r)*(1-bet*d_th) dQ = RHS
c
c           2     2       2
c  Where   d_x,  d_r and d_th  are the discrete differential operators of
c  the viscous terms in the axial, radial nad azimuthal directions
c  respectively. RHS is the righ hand side of the momentum equation as
c  is written above.
c
      subroutine invtrro
      include'param.f'
C     REAL etime_,t(2)
      alre=al/pec
c
c  compute the rhs of the factored equation everything at i,j+1/2,k+1/2
c
C     timhd0 = etime_(t)
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,kmv,kpv,n2m,jpv,jmv,n1m,ipv,imv,rc)
!$OMP$ SHARED(dens,ap1j,ac1j,am1j,aprk,acrk,amrk,pr)
!$OMP$ SHARED(ga,ro,hro,ruro,alre,dt,al,usrm,rhs)
!$OMP$ SHARED(dx1,dx1q,dx2,usinpm,fvisc1e,fvisc2e,fvisc3e)
!$OMP$ PRIVATE(km,kp,jp,jm,ip,im,d11q1,d22q1,d33q1,dcq1,dpx11,iadd)
!$OMP$ PRIVATE(h11,udx1,udx1nl,udx1q,app,acc,amm,ic,jc,kc)
        do 25 kc=1,n3m
          km=kmv(kc)
          kp=kpv(kc)
          do 25 jc=2,n2m-1
             udx1=dx1*usrm(jc)*usinpm(kc)*al
             udx1nl = dx1*usrm(jc)*0.25*al
             udx1q = dx1q*(usrm(jc)*usinpm(kc))**2
             app=ap1j(jc)
             acc=ac1j(jc)
             amm=am1j(jc)
             jp=jpv(jc)
             jm=jmv(jc)
             do 25 ic=1,n1m
                ip=ipv(ic)
                im=imv(ic)
c
c   11 second deriv. of dens(n)
c
            d11q1=(dens(ip,jc,kc)-dens(ic,jc,kc)*2.+dens(im,jc,kc))
     %            *udx1q
c
c   22 second deriv. of dens(n)
c
            d22q1=dens(ic,jp,kc)*app
     %           +dens(ic,jc,kc)*acc
     %           +dens(ic,jm,kc)*amm
c
c   33 second deriv. of dens(n)
c
            d33q1=dens(ic,jc,kp)*aprk(jc,kc)
     %           +dens(ic,jc,kc)*acrk(jc,kc)
     %           +dens(ic,jc,km)*amrk(jc,kc)
c
c   all viscous terms
c
c           dcq1 = d11q1 + d22q1 + d33q1
            dcq1 = d11q1*fvisc1e(jc,kc) + d22q1*fvisc2e(jc,kc) 
     %           + d33q1*fvisc3e(jc,kc)
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c   evaluation of RHS
c
            rhs(iadd)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc)
     %                +alre*dcq1)*dt
c
c   updating the nonlinear terms
c
                ruro(ic,jc,kc)=hro(ic,jc,kc)
   25 continue
!$OMP  END PARALLEL DO
      do 26 kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        jc=1
        udx1=dx1*usrm(jc)*usinpm(kc)*al
        udx1nl = dx1*usrm(jc)*0.25*al
        udx1q = dx1q*(usrm(jc)*usinpm(kc))**2
        app=ap1j(jc)
        acc=ac1j(jc)
        amm=am1j(jc)
        jp=jpv(jc)
        jm=jmv(jc)
        do 26 ic=1,n1m
          ip=ipv(ic)
          im=imv(ic)
c
c   11 second deriv. of dens(n)
c
          d11q1=(dens(ip,jc,kc)-dens(ic,jc,kc)*2.+dens(im,jc,kc))*
     %      udx1q
c
c   22 second deriv. of dens(n)
c
          d22q1=((dens(ic,jp,kc)-dens(ic,jc,kc))
     %          *dx2*rc(jp)**2/g2rc(jp)
     %          -(dens(ic,jc,kc)-denbs(ic,kc))
     %          *dx2*rc(jc)**2/g2rc(jc)*2.)
     %         *dx2/(g2rm(jc)*rm(jc)**2)
c
c   33 second deriv. of dens(n)
c
          d33q1=dens(ic,jc,kp)*aprk(jc,kc)
     %         +dens(ic,jc,kc)*acrk(jc,kc)
     %         +dens(ic,jc,km)*amrk(jc,kc)
c
c   all viscous terms
c
c         dcq1 = d11q1 + d22q1 + d33q1
          dcq1 = d11q1*fvisc1e(jc,kc) + d22q1*fvisc2e(jc,kc) 
     %         + d33q1*fvisc3e(jc,kc)
c
          iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c   evaluation of RHS
c
          rhs(iadd)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc)
     %              +alre*dcq1)*dt
c
c   updating the nonlinear terms
c
          ruro(ic,jc,kc)=hro(ic,jc,kc)
   26 continue
      do 27 kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        jc=n2m
        udx1=dx1*usrm(jc)*usinpm(kc)*al
        udx1nl = dx1*usrm(jc)*0.25*al
        udx1q = dx1q*(usrm(jc)*usinpm(kc))**2
        app=ap1j(jc)
        acc=ac1j(jc)
        amm=am1j(jc)
        jp=jc+1
        jm=jmv(jc)
        do 27 ic=1,n1m
          ip=ipv(ic)
          im=imv(ic)
c
c   11 second deriv. of dens(n)
c
          d11q1=(dens(ip,jc,kc)-dens(ic,jc,kc)*2.+dens(im,jc,kc))
     %          *udx1q
c
c   22 second deriv. of dens(n)
c
          d22q1=((denbn(ic,kc)-dens(ic,jc,kc))
     %           *dx2*rc(jp)**2/g2rc(jp)*2.
     %           -(dens(ic,jc,kc)-dens(ic,jm,kc))
     %           *dx2*rc(jc)**2/g2rc(jc))
     %          *dx2/(g2rm(jc)*rm(jc)**2)
c
c   33 second deriv. of dens(n)
c
          d33q1=dens(ic,jc,kp)*aprk(jc,kc)
     %         +dens(ic,jc,kc)*acrk(jc,kc)
     %         +dens(ic,jc,km)*amrk(jc,kc)
c
c   all viscous terms
c
c         dcq1 = d11q1 + d22q1 + d33q1
          dcq1 = d11q1*fvisc1e(jc,kc) + d22q1*fvisc2e(jc,kc) 
     %         + d33q1*fvisc3e(jc,kc)
c
          iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c   evaluation of RHS
c
          rhs(iadd)=(ga*hro(ic,jc,kc)+ro*ruro(ic,jc,kc)
     %                 +alre*dcq1)*dt
c
c   updating the nonlinear terms
c
          ruro(ic,jc,kc)=hro(ic,jc,kc)
   27 continue
c
c       timhd1 = etime_(t)
c
      call solroi
c       timhd2 = etime_(t)
      call solroj   !good
c       timhd3 = etime_(t)
      call solrok
c       timhd4 = etime_(t)
      return
      end
c
c
************************************************************************
c                       SUBROUTINE INVTR1
c   This subroutine performs the computation of Q~~ for the q1 momentum
c   equation (azimuthal direction) by a factored implicit scheme.
c   Viscous terms are treated implicitly, nonlinear terms explicitly.
c
c        ~~     n
c   dQ = Q  -  Q
c
c         alp dt
c   bet = -------
c          2 Re
c
c   The first equation in the introduction of THSCHEM then becomes
c
c
c              2      [            n         n       n-1   alp       2  n ]
c  (1-bet*nabla ) dQ =[-alp*grad (P ) + gam*H + rho*H   + ----- nabla (Q )]*dt
c                     [                                    Re             ]
c
c                                                                      3
c   The left hand side of this equation is then factored at an order dt
c   as follows
c
c           2           2           2
c  (1 -bet*d_x)*(1-bet*d_r)*(1-bet*d_th) dQ = RHS
c
c           2     2       2
c  Where   d_x,  d_r and d_th  are the discrete differential operators of
c  the viscous terms in the axial, radial nad azimuthal directions
c  respectively. RHS is the righ hand side of the momentum equation as
c  is written above.
c
      subroutine invtr1
      include'param.f'
C     REAL etime_,t(2)
      alre=al/ren
c
c  compute the rhs of the factored equation everything at i,j+1/2,k+1/2
c
C       timhd0 = etime_(t)
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,kmv,kpv,n2m,jpv,jmv,n1m,ipv,imv)
!$OMP$ SHARED(q1,ap1j,ac1j,am1j,ap1k,ac1k,am1k,pr)
!$OMP$ SHARED(ga,ro,dq,ru1,alre,dt,al,usrm,rhs)
!$OMP$ SHARED(dx1,dx1q,usinpm,fvisc1e,fvisc2e,fvisc3e)
!$OMP$ PRIVATE(km,kp,jp,jm,ip,im,d11q1,d22q1,d33q1,dcq1,dpx11,iadd)
!$OMP$ PRIVATE(h11,udx1,udx1nl,udx1q,app,acc,amm,ic,jc,kc)
      do 25 kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        do 25 jc=1,n2m
          udx1=dx1*usrm(jc)*usinpm(kc)*al
          udx1nl = dx1*usrm(jc)*0.25*al
          udx1q = dx1q*(usrm(jc)*usinpm(kc))**2
          app=ap1j(jc)
          acc=ac1j(jc)
          amm=am1j(jc)
          jp=jpv(jc)
          jm=jmv(jc)
          do 25 ic=1,n1m
            ip=ipv(ic)
            im=imv(ic)
c
c   11 second deriv. of q1(n)
c
            d11q1=(q1(ip,jc,kc)-q1(ic,jc,kc)*2.+q1(im,jc,kc))
     %            *udx1q
c
c   22 second deriv. of q1(n)
c
            d22q1=q1(ic,jp,kc)*app+q1(ic,jc,kc)*acc+q1(ic,jm,kc)*amm
c
c   33 second deriv. of q1(n)
c
            d33q1=q1(ic,jc,kp)*ap1k(jc,kc)
     %           +q1(ic,jc,kc)*ac1k(jc,kc)
     %           +q1(ic,jc,km)*am1k(jc,kc)
c
c   all viscous terms
c
c           dcq1 = d11q1 + d22q1 + d33q1
            dcq1 = d11q1*fvisc1e(jc,kc) + d22q1*fvisc2e(jc,kc) 
     %            + d33q1*fvisc3e(jc,kc)
c
c   grad(pr) along 1
c
            dpx11=(pr(ic,jc,kc)-pr(im,jc,kc))*udx1
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
c   evaluation of RHS
c
            rhs(iadd)=(ga*dq(ic,jc,kc)+ro*ru1(ic,jc,kc)
     %                 +alre*dcq1-dpx11)*dt
c
c   updating the nonlinear terms
c
            ru1(ic,jc,kc)=dq(ic,jc,kc)
   25 continue
!$OMP  END PARALLEL DO
c
c       timhd1 = etime_(t)
c
      call solq1i
c       timhd2 = etime_(t)
      call solq1j   !good
c       timhd3 = etime_(t)
      call solq1k
c       timhd4 = etime_(t)
c       do kc=1,n3m
c        do jc=1,n2m
c         do ic=1,n1m
c           iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c           q1(ic,jc,kc)=q1(ic,jc,kc) + rhs(iadd)
c         end do
c        end do
c       end do
c     write(6,*) ' TIME_I1 = ',timhd1-timhd0
c     write(6,*) ' TIME_S1I = ',timhd2-timhd1
c     write(6,*) ' TIME_S1J = ',timhd3-timhd2
c     write(6,*) ' TIME_S1K = ',timhd4-timhd3
c     write(6,*) ' TIME_TOT = ',timhd4-timhd0
      return
      end
************************************************************************
c                       SUBROUTINE INVTR2
c   This subroutine performs the computation of Q~~ for the q2 momentum
c   equation (radial direction) by a factored implicit scheme.
c   For details see the introduction of INVTR1
c
c
      subroutine invtr2
      include'param.f'
      alre=al/ren
c
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j,k+1/2
c
c    points inside the flowfield
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,kmv,kpv,n2m,n1m,ipv,imv,n2,rc)
!$OMP$ SHARED(dx1,dx1q,usinpm,fvisc1e,fvisc2e,fvisc3e)
!$OMP$ SHARED(q2,ap2j,ac2j,am2j,ap2k,ac2k,am2k,pr)
!$OMP$ SHARED(ga,ro,dph,ru2,alre,dt,udx2c,usrc,al,rhs)
!$OMP$ PRIVATE(km,kp,jp,jm,ip,im,d11q2,d22q2,d33q2)
!$OMP$ PRIVATE(aap,aac,aam,dcq2,dpx22,iadd)
!$OMP$ PRIVATE(h21,udx1,udx2,udx1q,ic,jc,kc)
      do 26 kc=1,n3m
        km=kmv(kc)
        kp=kpv(kc)
        do 26 jc=2,n2m
          jm=jc-1
          jp=jc+1
          udx1q =  dx1q*(usrc(jc)*usinpm(kc))**2
          udx2 = rc(jc)**2*al*udx2c(jc)
          udx1 = dx1*0.25*usrc(jc)*al*0.5
          do 26 ic=1,n1m
            im=imv(ic)
            ip=ipv(ic)
c
c   11 second derivative of q2
c
            d11q2=(q2(ip,jc,kc)-2.*q2(ic,jc,kc)+q2(im,jc,kc))
     %            *udx1q
c
c   22 second derivative of q2
c
            aap=ap2j(jc)
            aac=ac2j(jc)
            aam=am2j(jc)
            d22q2=q2(ic,jp,kc)*ap2j(jc)
     %           +q2(ic,jc,kc)*ac2j(jc) 
     %           +q2(ic,jm,kc)*am2j(jc)
c
c   33 second derivative of q2
c
            d33q2=q2(ic,jc,kp)*ap2k(jc,kc)
     %           +q2(ic,jc,kc)*ac2k(jc,kc)
     %           +q2(ic,jc,km)*am2k(jc,kc)
c
            dcq2 = d11q2*fvisc1e(jc,kc) 
     %            +d22q2*fvisc2e(jc,kc)
     %            +d33q2*fvisc3e(jc,kc)
c
c   component of grad(pr) along 2 direction
c
            dpx22=(pr(ic,jc,kc)-pr(ic,jm,kc))*udx2
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
c
            rhs(iadd)=(ga*dph(ic,jc,kc)+ro*ru2(ic,jc,kc)
     %                    +alre*dcq2-dpx22)*dt
c
            ru2(ic,jc,kc)=dph(ic,jc,kc)
   26 continue
!$OMP  END PARALLEL DO
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n1m,n2)
!$OMP$ SHARED(rhs)
!$OMP$ PRIVATE(ic,kc,iaddi,iaddf)
      do kc=1,n3m
        do ic=1,n1m
         iaddi=ic+(kc-1)*n1m*n2
         iaddf=ic+(n2-1)*n1m+(kc-1)*n1m*n2
         rhs(iaddi)=0.
         rhs(iaddf)=0.
        end do
      end do
!$OMP  END PARALLEL DO
      call solq2i
      call solq2k
      call solq2j
      return
      end
c
c************************************************************************
c                       SUBROUTINE INVTR3
c   This subroutine performs the computation of Q~~ for the q3 momentum
c   equation (axial direction) by a factored implicit scheme.
c   Viscous terms are treated implicitly, nonlinear terms explicitly.
c   For details see the introduction of INVTR1
c
      subroutine invtr3
      include'param.f'
      alre=al/ren
c
      small = 1.e-10
      icon = 0
c
c  compute the rhs of the factored equation
c  everything at i+1/2,j+1/2,k
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(kmv,jpv,jmv,n3m,n2m,n1m,ipv,imv)
!$OMP$ SHARED(q3,ap3j,ac3j,am3j,ap3k,ac3k,am3k,pr)
!$OMP$ SHARED(sinpc,usinpc,dx1,dx1q)
!$OMP$ SHARED(fvisc1e,fvisc2e,fvisc3e)
!$OMP$ SHARED(al,ga,ro,qcap,ru3,dq,alre,dt,udx3c,usrm,rhs)
!$OMP$ PRIVATE(km,kp,jp,jpp,jmm,ip,ipp,imm,dq31,dq32,dq33,dcq3,dpx33)
!$OMP$ PRIVATE(h31,iadd,aap,aac,aam,udx1q,udx1,udx3,ic,jc,kc)
      do 18 kc=3,n3m-1
c     do 18 kc=2,n3m
        km=kmv(kc)
        kp=kc+1
        do 18 jc=1,n2m
          udx3 = al*udx3c(kc)*sinpc(kc)*usrm(jc)
          jmm=jmv(jc)
          jpp=jpv(jc)
          jp=jc+1
          aap=ap3j(jc)
          aam=am3j(jc)
          aac=ac3j(jc)
          udx1=dx1*usrm(jc)*0.25*al*0.5
          udx1q = dx1q*(usrm(jc)*usinpc(kc))**2
          do 18 ic=1,n1m
            imm=imv(ic)
            ipp=ipv(ic)
            ip=ipp
c
c   11 second derivatives of q3
c
            dq31=(q3(ipp,jc,kc)-2.*q3(ic,jc,kc)+q3(imm,jc,kc))*udx1q
c
c   22 second derivatives of q3
c
            dq32=aam*q3(ic,jmm,kc)+aac*q3(ic,jc,kc)+aap*q3(ic,jpp,kc)
c
c   33 second derivatives of q3
c
            dq33=q3(ic,jc,kp)*ap3k(jc,kc)
     %          +q3(ic,jc,kc)*ac3k(jc,kc)
     %          +q3(ic,jc,km)*am3k(jc,kc)
C     if(rm(jc).le.0.35) then
C     icon = icon+1
C     fff = (1.-cos(2.*zz(kc)))/rm(jc)**2
C     if(abs(fff-dq33).gt.1.e-9) then 
C      ferr = (fff-dq33)/((abs(fff)+small))
C     else
C      ferr = small
C     end if
C     write(888,1888) icon,ic,jc,kc,thetac(ic),rm(jc),zm(kc),
C    % q3(ic,jc,kc),q3(ic,jc,kc),fff,dq33,ferr
C     end if
c     if(rm(jc).le.0.35) then
c     icon = icon+1
c     fff = -cos(thetam(ic))*(sin(2.*zz(kc))+cos(zz(kc))/sin(zz(kc)))
c    %       /rm(jc)**2
c     if(abs(fff-dq33).gt.1.e-9) then 
c      ferr = (fff-h11)/((abs(fff)+small))
c     else
c      ferr = small
c     end if
c     write(888,1888) icon,ic,jc,kc,thetac(ic),rm(jc),zm(kc),
c    % q1(ic,jc,kc),q1(ic,jc,kc),fff,h11,ferr
c     end if
c1888 format(i5,3(1x,i3),8(1x,e12.5))
c
c   viscous terms
c
            dcq3=dq31*fvisc1e(jc,kc)+dq32*fvisc2e(jc,kc)
     %          +dq33*fvisc3e(jc,kc)
c
c  component of grad(pr) along x3 direction
c
            dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3
c
c  right hand side of the momentum equation
c
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
            rhs(iadd)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc)
     %                    +alre*dcq3-dpx33)*dt
c
c  updating of the non-linear terms
c
            ru3(ic,jc,kc)=qcap(ic,jc,kc)
   18 continue
!$OMP  END PARALLEL DO

c Special treatment for north and south pole

c     go to 5566

      kc=2
      kp=kc+1
      km=kc-1
!$OMP  PARALLEL DO
!$OMP$ SHARED(kc,km,kp,jpv,jmv,n2m,n1m,ipv,imv)
!$OMP$ SHARED(ksyc,isym,qaxN,dx1,dx1q,dx3q)
!$OMP$ SHARED(fvisc1e,fvisc2e,fvisc3e)
!$OMP$ SHARED(q3,ap3j,ac3j,am3j,pr,rm,g3rc,g3rm)
!$OMP$ SHARED(sinpc,usinpc)
!$OMP$ SHARED(al,ga,ro,qcap,ru3,alre,dt,udx3c,usrm,rhs)
!$OMP$ PRIVATE(knsy,insy,insym)
!$OMP$ PRIVATE(jp,jpp,jmm,ip,ipp,imm)
!$OMP$ PRIVATE(q3sym1,dq31,dq32,dq33,dcq3,dpx33)
!$OMP$ PRIVATE(h31,iadd,aap,aac,aam,udx1q,udx1,udx3,ic,jc)
      do jc=1,n2m
        udx3 = al*udx3c(kc)*sinpc(kc)*usrm(jc)
        jmm=jmv(jc)
        jpp=jpv(jc)
        jp=jc+1
        aap=ap3j(jc)
        aam=am3j(jc)
        aac=ac3j(jc)
        udx1=dx1*usrm(jc)*0.25*al*0.5
        udx1q = dx1q*(usrm(jc)*usinpc(kc))**2
        do ic=1,n1m
          imm=imv(ic)
          ipp=ipv(ic)
          ip=ipp
          knsy = ksyc(kc)
          insy = isym(ic)
          insym = isym(imm)
c
c   33 second derivatives of q3 modified for poles
c

c         q3sym1=0.5*(usinpc(kc)*q3(ic,jc,kc)-
c    %                usinpc(kc)*q3(insy,jc,kc))
c         q3sym1=0.5*(q3(ic,jc,kc)-q3(insy,jc,kc))/sin(zz(kc))
          q3sym1 = qaxN(3,ic,jc)
          dq33=dx3q/(g3rc(kc)*rm(jc)**2)*(
     %    sin(zm(kc))*(q3(ic,jc,kp)/sin(zz(kp))
     %                -q3(ic,jc,kc)/sin(zz(kc)))/g3rm(kc)-
     %    sin(zm(km))*(q3(ic,jc,kc)/sin(zz(kc))-q3sym1)/g3rm(km)) 
     %    -q3(ic,jc,kc)/(rm(jc)*sin(zz(kc)))**2
C     if(rm(jc).le.0.35) then
C     icon = icon+1
C     fff = (1.-cos(2.*zz(kc)))/rm(jc)**2
C     if(abs(fff-dq33).gt.1.e-9) then 
C      ferr = (fff-dq33)/((abs(fff)+small))
C     else
C      ferr = small
C     end if
C     write(888,1888) icon,ic,jc,kc,thetac(ic),rm(jc),zm(kc),
C    % q3(ic,jc,kc),q3(ic,jc,kc),fff,dq33,ferr
C     end if
C1888 format(i5,3(1x,i3),8(1x,e12.5))
c           dq33=dx3q/(rm(jc)**2*g3rc(kc))*(
c    %      sinpm(kc)*(q3(ic,jc,kp)*usinpc(kp)-q3(ic,jc,kc)*usinpc(kc))
c    %               /g3rm(kc)-
c    %      sinpm(km)*(q3(ic,jc,kc)*usinpc(kc)-q3sym1)/g3rm(km)) -
c    %       q3(ic,jc,kc)/(rm(jc)*sinpc(kc))**2
c           dq33 = 0.
c           dq33=dx3/rc(jc)**2*(dx3/g3rm(kc)*sinpm(kc)*
c    %      (q3(ic, jc, kp)*usinpc(kp)-q3(ic, jc, kc)*usinpc(kc))-
c    %      dx3/g3rm(kc)*sinpm(km)*(q3(ic, jc, kc)*usinpc(kc)-q3sym1))
c           dq33=0.
c           write(*,*) ic,jc,q3sym1,dq33
c           write(*,*) usinpc(kc),g3rm(kc),g3rc(kc)
c           write(*,*) 
c    %      (q3(ic,jc,kp)*usinpc(kp)-q3(ic,jc,kc)*usinpc(kc))/g3rm(kc)
c           write(*,*) 
c    %      (q3(ic,jc,kc)*usinpc(kc)-q3sym1)/g3rm(km) 
c           write(*,*) 
c    %       q3(ic,jc,kc)/(rm(jc)*sinpc(kc))**2
c           pause

c!    a33=dx3q/(g3rc(kc)*rm(jc)**2)
c!    a33p= +a33*sinpm(kc)/g3rm(kc)
c!    a33m= +a33*sinpm(km)/g3rm(km)
c!    ap3k(jc,kc)=a33p/sinpc(kp)
c!    am3k(jc,kc)=a33m/sinpc(km)
c!    ac3k(jc,kc)=-(a33p+a33m)/sinpc(kc)-usrmq(jc)*usinpc(kc)**2
c
c   11 second derivatives of q3
c
          dq31=(q3(ipp,jc,kc)-2.*q3(ic,jc,kc)+q3(imm,jc,kc))*udx1q
c
c   22 second derivatives of q3
c
          dq32=aam*q3(ic,jmm,kc)+aac*q3(ic,jc,kc)+aap*q3(ic,jpp,kc)
c
c   viscous terms
c
          dcq3=dq31*fvisc1e(jc,kc)+dq32*fvisc2e(jc,kc)
     %        +dq33*fvisc3e(jc,kc)
c
c  component of grad(pr) along x3 direction
c
          dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3
c
c  right hand side of the momentum equation
c
          iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
          rhs(iadd)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc)
     %              +alre*dcq3-dpx33)*dt
c
c  updating of the non-linear terms
c
            ru3(ic,jc,kc)=qcap(ic,jc,kc)
        end do
      end do
!$OMP  END PARALLEL DO

      kc=n3m
      kp=kc+1
      km=kc-1
!$OMP  PARALLEL DO
!$OMP$ SHARED(kc,km,kp,jpv,jmv,n2m,n1m,ipv,imv)
!$OMP$ SHARED(ksym,isym,qaxN,dx1,dx1q,dx3q)
!$OMP$ SHARED(fvisc1e,fvisc2e,fvisc3e)
!$OMP$ SHARED(q3,ap3j,ac3j,am3j,pr,rm,g3rc,g3rm)
!$OMP$ SHARED(sinpc,sinpm,usinpc)
!$OMP$ SHARED(al,ga,ro,qcap,ru3,alre,dt,udx3c,usrm,rhs)
!$OMP$ PRIVATE(knsy,insy,insym)
!$OMP$ PRIVATE(jp,jpp,jmm,ip,ipp,imm)
!$OMP$ PRIVATE(q3symn3m,dq31,dq32,dq33,dcq3,dpx33)
!$OMP$ PRIVATE(h31,iadd,aap,aac,aam,udx1q,udx1,udx3,ic,jc)
      do jc=1,n2m
        udx3 = al*udx3c(kc)*sinpc(kc)*usrm(jc)
        jmm=jmv(jc)
        jpp=jpv(jc)
        jp=jc+1
        aap=ap3j(jc)
        aam=am3j(jc)
        aac=ac3j(jc)
        udx1=dx1*usrm(jc)*0.25*al*0.5
        udx1q = dx1q*(usrm(jc)*usinpc(kc))**2
        do ic=1,n1m
          imm=imv(ic)
          ipp=ipv(ic)
          ip=ipp
          knsy = ksym(kc)
          insy = isym(ic)
          insym = isym(imm)

c
c   11 second derivatives of q3
c
          dq31=(q3(ipp,jc,kc)-2.*q3(ic,jc,kc)+q3(imm,jc,kc))*udx1q
c
c   22 second derivatives of q3
c
          dq32=aam*q3(ic,jmm,kc)+aac*q3(ic,jc,kc)+aap*q3(ic,jpp,kc)
c
c   33 second derivatives of q3 modified for north pole
c

c         q3symn3m=0.5*usinpc(kc)*(q3(ic, jc, kc)-q3(insy, jc, kc))
c         dq33=dx3/rc(jc)**2*(dx3/g3rm(kc)*sinpm(kc)*
c    %    (q3symn3m-q3(ic, jc, kc)*usinpc(kc))- dx3/g3rm(kc)*sinpm(km)
c    %    *(q3(ic, jc, kc)*usinpc(kc)-q3(ic, jc, kp)*usinpc(kp)))
c         q3symn3m=0.5*(usinpc(kc)*q3(ic,jc,kc)-
c    %                    usinpc(kc)*q3(insy,jc,kc))
          q3symn3m = qaxN(3,ic,jc)
          dq33=dx3q/(rm(jc)**2*g3rc(kc))
     %      *(sinpm(kc)*(q3symn3m-q3(ic,jc,kc)*usinpc(kc))/g3rm(kc)
     %        -sinpm(km)*(q3(ic,jc,kc)*usinpc(kc)
     %                    -q3(ic,jc,km)*usinpc(km))/g3rm(km))
     %      -q3(ic,jc,kc)/(rm(jc)*sinpc(kc))**2
c     if(rm(jc).le.0.35) then
c     icon = icon+1
c     fff = (1.-cos(2.*zz(kc)))/rm(jc)**2
c     if(abs(fff-dq33).gt.1.e-9) then 
c      ferr = (fff-dq33)/((abs(fff)+small))
c     else
c      ferr = small
c     end if
c     write(888,1888) icon,ic,jc,kc,thetac(ic),rm(jc),zm(kc),
c    % q3(ic,jc,kc),q3(ic,jc,kc),fff,dq33,ferr
c     end if
c           dq33 = 0.
c
c   viscous terms
c
          dcq3=dq31*fvisc1e(jc,kc)+dq32*fvisc2e(jc,kc)
     %        +dq33*fvisc3e(jc,kc)
c
c  component of grad(pr) along x3 direction
c
          dpx33=(pr(ic,jc,kc)-pr(ic,jc,km))*udx3
c
c  right hand side of the momentum equation
c
          iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c
          rhs(iadd)=(ga*qcap(ic,jc,kc)+ro*ru3(ic,jc,kc)
     %              +alre*dcq3-dpx33)*dt
c
c  updating of the non-linear terms
c
          ru3(ic,jc,kc)=qcap(ic,jc,kc)
        end do
      end do
c 5566   continue
!$OMP  END PARALLEL DO


c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n2m,n1m,n3)
!$OMP$ SHARED(rhs)
!$OMP$ PRIVATE(ic,jc,iaddi,iaddf)
      do jc=1,n2m
        do ic=1,n1m
         iaddi=ic+(jc-1)*n1m
         iaddf=ic+(jc-1)*n1m+(n3-1)*n1m*n2m
         rhs(iaddi)=0.
         rhs(iaddf)=0.
        end do
      end do
!$OMP  END PARALLEL DO
      call solq3i
      call solq3j
      call solq3k
c
      return
      end
c
c
c
c  ****************************** subrout trvb  **********************
c
      subroutine trvb(am,ac,ap,f,mm,nn,n,m)
c
      dimension am(mm,nn),ac(mm,nn),ap(mm,nn),f(mm,nn)
c
c  ******** reduction of trid. matrix to an upper rigth matrix
c
      do 1 i=2,n
      do 2 k=1,m
      ac(k,i)=ac(k,i)*ac(k,i-1)-ap(k,i-1)*am(k,i)
      ap(k,i)=ap(k,i)*ac(k,i-1)
      f(k,i)=f(k,i)*ac(k,i-1)-f(k,i-1)*am(k,i)
    2 continue
    1 continue
c  ******** calculation of the unknown by backward elimination
c
      do 3 k=1,m
      f(k,n)=f(k,n)/ac(k,n)
    3 continue
      nm=n-1
      do 10 ii=1,nm
      i=n-ii
      do 11 k=1,m
      f(k,i)=(f(k,i)-ap(k,i)*f(k,i+1))/ac(k,i)
   11 continue
   10 continue
      return
      end

c
c
c
c     it inverts in 1 direction, vectors of n1m length
c     from 1, n2m
c
c
c
c***********************************************************************
      subroutine newtrip(ami,aci,api,fi,n1i,n1f,n2i,n2f,nn1,nn2)
      dimension ami(nn1,nn2),aci(nn1,nn2),api(nn1,nn2),fi(nn1,nn2)
     %         ,fei(nn1,nn2),q(nn1,nn2),s(nn1,nn2)
C      include'param.f'
c
c     vectorized for right hand side and coefficients
c
      dimension fn(nn2),p(nn2)
      ia = n1i + 1
      ii = n1i + n1f
c
c   COEFFICIENTS FOR TRIDIAGONAL INVERSION
c

c
c  THE INVERSION STARTS
c
      do 20 j=n2i,n2f
        q(n1i,j) = -api(n1i,j)/aci(n1i,j)
        s(n1i,j) = -ami(n1i,j)/aci(n1i,j)
        fn(j) = fi(n1f,j)
        fi(n1i,j) = fi(n1i,j)/aci(n1i,j)
   20 continue
c
c     forward elimination sweep
c
      do 10 i=ia,n1f
        do 21 j=n2i,n2f
          p(j) =1./( aci(i,j) + ami(i,j)*q(i-1,j))
          q(i,j) = - api(i,j)*p(j)
          s(i,j) = - ami(i,j)*s(i-1,j)*p(j)
          fi(i,j) = ( fi(i,j) - ami(i,j)*fi(i-1,j))*p(j)
   21   continue
   10 continue
c
c     backward pass
c
      do 22 j=n2i,n2f
        s(n1f,j) = 1.
        fei(n1f,j) = 0.
   22 continue
      do 11 l=ia,n1f
        i = ii - l
        do 23 j=n2i,n2f
          s(i,j) = s(i,j) + q(i,j)*s(i+1,j)
          fei(i,j) = fi(i,j) + q(i,j)*fei(i+1,j)
   23   continue
   11 continue
      do 24 j=n2i,n2f
        fi(n1f,j)=(fn(j)-api(i,j)*fei(n1i,j) -
     %      ami(i,j)*fei(n1f-1,j))/(api(i,j)*s(n1i,j) +
     %      ami(i,j)*s(n1f-1,j)+aci(i,j))
   24 continue
c
c     backward elimination pass
c
      do 12 l=ia,n1f
        i = ii -l
        do 25 j=n2i,n2f
          fi(i,j) = fi(n1f,j)*s(i,j) + fei(i,j)
   25   continue
   12 continue
      return
      end

c
c
************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c         direction x1
c
      subroutine solroi
      include'param.f'
      dimension amil(m2,m1),acil(m2,m1),apil(m2,m1),fil(m2,m1)
c
c  ********* compute dq1*  sweeping in the x1 direction
c            inflow outflow
c
c
c   Q_1
c
      betadx=0.5*al*dt/pec*dx1q
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,rhs,usrm,usinpm,fvisc1)
!$OMP$ PRIVATE(amil,apil,acil,fil,iadd,faci)
!$OMP$ PRIVATE(ugmmm,ugm2,ic,jc,kc)
      do kc=1,n3m
        do ic=1,n1m
          do jc=1,n2m
            ugmmm=betadx*(usrm(jc)*usinpm(kc))**2
            faci=1.+ugmmm*2.*fvisc1(jc,kc)
            apil(jc,ic)=-ugmmm*fvisc1(jc,kc)/faci
            amil(jc,ic)=-ugmmm*fvisc1(jc,kc)/faci
            acil(jc,ic)=1.
            iadd = ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fil(jc,ic) = rhs(iadd)/faci
          end do
        end do
c 
        call tripvmy(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)
c 
        do ic=1,n1m
          do jc=1,n2m
            iadd = ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            rhs(iadd) = fil(jc,ic)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
      return
      end
c
c
************************************************************************
c
c
c
c************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c
      subroutine solroj
      include'param.f'
      dimension amjl(m2),apjl(m2),acjl(m2),fjl(m1,m2)
c
c  ************ compute dq1  sweeping along the x2 direction
c               wall boundaries interior points
c
      betadx=0.5*al*dt/pec
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,rhs,ap1j,ac1j,fvisc2)
!$OMP$ PRIVATE(amjl,apjl,acjl,fjl,iadd,faci,nm,ic,ii,jc,kc)
      do kc=1,n3m
        do jc=1,n2m
          faci = 1.-ac1j(jc)*betadx*fvisc2(jc,kc)
          apjl(jc)= -ap1j(jc)*betadx*fvisc2(jc,kc)/faci
          acjl(jc)=1.
          amjl(jc)= -am1j(jc)*betadx*fvisc2(jc,kc)/faci
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fjl(ic,jc)=rhs(iadd)/faci
          end do
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        do jc=2,n2m
          acjl(jc)=acjl(jc)*acjl(jc-1)-apjl(jc-1)*amjl(jc)
          apjl(jc)=apjl(jc)*acjl(jc-1)
          do ic=1,n1m
            fjl(ic,jc)=fjl(ic,jc)*acjl(jc-1)-fjl(ic,jc-1)*amjl(jc)
          end do
        end do
c
        do ic=1,n1m
          iadd=ic+(n2m-1)*n1m+(kc-1)*n1m*n2m
          fjl(ic,n2m)=fjl(ic,n2m)/acjl(n2m)
          rhs(iadd)=fjl(ic,n2m)
        end do
c
        nm=n2m-1
        do ii=1,nm
          jc=n2m-ii
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fjl(ic,jc)=(fjl(ic,jc)-apjl(jc)*fjl(ic,jc+1))/acjl(jc)
            rhs(iadd)=fjl(ic,jc)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     call ntrvb(amjl,acjl,apjl,fjl,m1,m2,n2m,n1m)
c
c       do jc=1,n2m
c         do ic=1,n1m
c           iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c           rhs(iadd) = fjl(ic,jc)
c         end do
c       end do
c     end do
      return
      end
c
c
************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c       direction x3
c
      subroutine solrok
      include'param.f'
      dimension amkl(m3),apkl(m3),ackl(m3),fkl(m1,m3)
c
c   compute  from dq** sweeping along the x3 periodic direction
c
      betadx=0.5*al*dt/pec
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,amrk,acrk,aprk)
!$OMP$ SHARED(betadx,rhs,dens,fvisc3)
!$OMP$ PRIVATE(amkl,apkl,ackl,fkl,iadd,faci,nm,ic,ii,jc,kc)
      do jc=1,n2m
        do kc=1,n3m
          faci = 1.-acrk(jc,kc)*betadx*fvisc3(jc,kc)
          amkl(kc)=-amrk(jc,kc)*betadx*fvisc3(jc,kc)/faci
          ackl(kc)=1.
          apkl(kc)=-aprk(jc,kc)*betadx*fvisc3(jc,kc)/faci
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fkl(ic,kc)=rhs(iadd)/faci
          end do
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        do kc=2,n3m
          ackl(kc)=ackl(kc)*ackl(kc-1)-apkl(kc-1)*amkl(kc)
          apkl(kc)=apkl(kc)*ackl(kc-1)
          do ic=1,n1m
            fkl(ic,kc)=fkl(ic,kc)*ackl(kc-1)-fkl(ic,kc-1)*amkl(kc)
          end do
        end do
c
        do ic=1,n1m
          fkl(ic,n3m)=fkl(ic,n3m)/ackl(n3m)
          dens(ic,jc,n3m)=dens(ic,jc,n3m) + fkl(ic,n3m)
        end do
c
        nm=n3m-1
        do ii=1,nm
          kc=n3m-ii
          do ic=1,n1m
c           q1(ic,jc,kc)=q1(ic,jc,kc) +
            fkl(ic,kc)=(fkl(ic,kc)-apkl(kc)*fkl(ic,kc+1))/ackl(kc)
            dens(ic,jc,kc)=dens(ic,jc,kc) + fkl(ic,kc)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     call ntrvb(amkl,ackl,apkl,fkl,m1,m3,n3m,n1m)
c
c       do kc=1,n3m
c         do ic=1,n1m
c           q1(ic,jc,kc)=q1(ic,jc,kc) + fkl(ic,kc)
c         end do
c       end do
c     end do
c
c
      return
      end
c
c
c
c
************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c         direction x1
c
      subroutine solq1i
      include'param.f'
      dimension amil(m2,m1),acil(m2,m1),apil(m2,m1),fil(m2,m1)
c
c  ********* compute dq1*  sweeping in the x1 direction
c            inflow outflow
c
c
c   Q_1
c
      betadx=beta*dx1q*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,rhs,usrm,usinpm,fvisc1)
!$OMP$ PRIVATE(amil,apil,acil,fil,iadd,faci)
!$OMP$ PRIVATE(ugmmm,ugm2,ic,jc,kc)
      do kc=1,n3m
        do ic=1,n1m
          do jc=1,n2m
            ugmmm=betadx*(usrm(jc)*usinpm(kc))**2
            faci=1.+ugmmm*2.*fvisc1(jc,kc)
            apil(jc,ic)=-ugmmm*fvisc1(jc,kc)/faci
            amil(jc,ic)=-ugmmm*fvisc1(jc,kc)/faci
            acil(jc,ic)=1.
            iadd = ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fil(jc,ic) = rhs(iadd)/faci
          end do
        end do
c
        call tripvmy(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)
c
        do ic=1,n1m
          do jc=1,n2m
            iadd = ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            rhs(iadd) = fil(jc,ic)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
      return
      end
c
c
************************************************************************
c
c
c
c************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c
      subroutine solq1j
      include'param.f'
      dimension amjl(m2),apjl(m2),acjl(m2),fjl(m1,m2)
c
c  ************ compute dq1  sweeping along the x2 direction
c               wall boundaries interior points
c
      betadx=beta*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,rhs,ap1j,ac1j,am1j,fvisc2)
!$OMP$ PRIVATE(amjl,apjl,acjl,fjl,iadd,faci,nm,ic,ii,jc,kc)
      do kc=1,n3m
        do jc=1,n2m
          faci = 1.-ac1j(jc)*betadx*fvisc2(jc,kc)
          apjl(jc)= -ap1j(jc)*betadx*fvisc2(jc,kc)/faci
          acjl(jc)=1.
          amjl(jc)= -am1j(jc)*betadx*fvisc2(jc,kc)/faci
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fjl(ic,jc)=rhs(iadd)/faci
          end do
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        do jc=2,n2m
          acjl(jc)=acjl(jc)*acjl(jc-1)-apjl(jc-1)*amjl(jc)
          apjl(jc)=apjl(jc)*acjl(jc-1)
          do ic=1,n1m
            fjl(ic,jc)=fjl(ic,jc)*acjl(jc-1)-fjl(ic,jc-1)*amjl(jc)
          end do
        end do
c
        do ic=1,n1m
          iadd=ic+(n2m-1)*n1m+(kc-1)*n1m*n2m
          fjl(ic,n2m)=fjl(ic,n2m)/acjl(n2m)
          rhs(iadd)=fjl(ic,n2m)
        end do
c
        nm=n2m-1
        do ii=1,nm
          jc=n2m-ii
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fjl(ic,jc)=(fjl(ic,jc)-apjl(jc)*fjl(ic,jc+1))/acjl(jc)
            rhs(iadd)=fjl(ic,jc)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     call ntrvb(amjl,acjl,apjl,fjl,m1,m2,n2m,n1m)
c
c       do jc=1,n2m
c         do ic=1,n1m
c           iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c           rhs(iadd) = fjl(ic,jc)
c         end do
c       end do
c     end do
      return
      end
c
c
************************************************************************
c
c   this subroutine performs the inversion of the q1 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q1
c   are treated implicitly
c       direction x3
c
      subroutine solq1k
      include'param.f'
      dimension amkl(m3),apkl(m3),ackl(m3),fkl(m1,m3)
c
c   compute  from dq** sweeping along the x3 periodic direction
c
      betadx=beta*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,am1k,ac1k,ap1k,fvisc3)
!$OMP$ SHARED(q1,betadx,rhs)
!$OMP$ PRIVATE(amkl,apkl,ackl,fkl,iadd,faci,nm,ic,ii,jc,kc)
      do jc=1,n2m
        do kc=1,n3m
          faci = 1.-ac1k(jc,kc)*betadx*fvisc3(jc,kc)
          amkl(kc)=-am1k(jc,kc)*betadx*fvisc3(jc,kc)/faci
          ackl(kc)=1.
          apkl(kc)=-ap1k(jc,kc)*betadx*fvisc3(jc,kc)/faci
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fkl(ic,kc)=rhs(iadd)/faci
          end do
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        do kc=2,n3m
          ackl(kc)=ackl(kc)*ackl(kc-1)-apkl(kc-1)*amkl(kc)
          apkl(kc)=apkl(kc)*ackl(kc-1)
          do ic=1,n1m
            fkl(ic,kc)=fkl(ic,kc)*ackl(kc-1)-fkl(ic,kc-1)*amkl(kc)
          end do
        end do
c
        do ic=1,n1m
          fkl(ic,n3m)=fkl(ic,n3m)/ackl(n3m)
          q1(ic,jc,n3m)=q1(ic,jc,n3m) + fkl(ic,n3m)
        end do
c
        nm=n3m-1
        do ii=1,nm
          kc=n3m-ii
          do ic=1,n1m
c           q1(ic,jc,kc)=q1(ic,jc,kc) +
            fkl(ic,kc)=(fkl(ic,kc)-apkl(kc)*fkl(ic,kc+1))/ackl(kc)
            q1(ic,jc,kc)=q1(ic,jc,kc) + fkl(ic,kc)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c     call ntrvb(amkl,ackl,apkl,fkl,m1,m3,n3m,n1m)
c
c       do kc=1,n3m
c         do ic=1,n1m
c           q1(ic,jc,kc)=q1(ic,jc,kc) + fkl(ic,kc)
c         end do
c       end do
c     end do
c
c
      return
      end
c
c
************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c        direction x1
c
      subroutine solq2i
      include'param.f'
      dimension amil(m2,m1),acil(m2,m1),apil(m2,m1),fil(m2,m1)
c  ************ compute dq2** sweeping along the x1 direction
c
      betadx=beta*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n2,n1m,betadx,rhs,usrc,udx1qc)
!$OMP$ SHARED(usinpm,fvisc1)
!$OMP$ PRIVATE(amil,apil,acil,fil,iadd,faci)
!$OMP$ PRIVATE(ugmmm,ic,jc,kc)
      do kc=1,n3m
       do ic=1,n1m
        do jc=2,n2m
         ugmmm=betadx*udx1qc(jc)*usinpm(kc)**2
         faci=1.+ugmmm*2.*fvisc1(jc,kc)
         apil(jc,ic)=-ugmmm*fvisc1(jc,kc)/faci
         amil(jc,ic)=-ugmmm*fvisc1(jc,kc)/faci
         acil(jc,ic)=1.
         iadd = ic+(jc-1)*n1m+(kc-1)*n1m*n2
         fil(jc,ic)=rhs(iadd)/faci
        end do
       end do
       do ic=1,n1m
        apil(1,ic)=0.
        amil(1,ic)=0.
        acil(1,ic)=1.
        fil(1,ic)=0.
        apil(n2,ic)=0.
        amil(n2,ic)=0.
        acil(n2,ic)=1.
        fil(n2,ic)=0.
       end do
c
       call tripvmy(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)
c
       do ic=1,n1m
        do jc=1,n2m
         iadd = ic+(jc-1)*n1m+(kc-1)*n1m*n2
         rhs(iadd)=fil(jc,ic)
        end do
       end do
      end do
!$OMP  END PARALLEL DO

c
      return
      end
c
c
************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c         direction x2
c
      subroutine solq2j
      include'param.f'
      dimension amjl(m2),apjl(m2),acjl(m2),fjl(m1,m2)
c  ********* compute the dq2* sweeping in the x2 direction
c            wall boundaries direction
      betadx=beta*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,rhs,n2,ap2j,ac2j,am2j)
!$OMP$ SHARED(q2,fvisc2)
!$OMP$ PRIVATE(amjl,apjl,acjl,fjl,iadd,faci,ic,ii,jc,kc)
      do kc=1,n3m
        do jc=2,n2m
          faci =1.-ac2j(jc)*betadx*fvisc2(jc,kc)
          apjl(jc)=-ap2j(jc)*betadx*fvisc2(jc,kc)/faci
          acjl(jc)=1.
          amjl(jc)=-am2j(jc)*betadx*fvisc2(jc,kc)/faci
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
            fjl(ic,jc)=rhs(iadd)/faci
          end do
        end do
c
        amjl(1)=0.
        apjl(1)=0.
        acjl(1)=1.
        amjl(n2)=0.
        apjl(n2)=0.
        acjl(n2)=1.
        do ic=1,n1m
          fjl(ic,1)=0.
          fjl(ic,n2)=0.
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        do jc=2,n2
          acjl(jc)=acjl(jc)*acjl(jc-1)-apjl(jc-1)*amjl(jc)
          apjl(jc)=apjl(jc)*acjl(jc-1)
          do ic=1,n1m
            fjl(ic,jc)=fjl(ic,jc)*acjl(jc-1)-fjl(ic,jc-1)*amjl(jc)
          end do
        end do
c
        do ic=1,n1m
          fjl(ic,n2)=fjl(ic,n2)/acjl(n2)
          q2(ic,n2,kc)=q2(ic,n2,kc) + fjl(ic,n2)
        end do
c
        nm=n2-1
        do ii=1,nm
          jc=n2-ii
          do ic=1,n1m
            fjl(ic,jc)=(fjl(ic,jc)-apjl(jc)*fjl(ic,jc+1))/acjl(jc)
            q2(ic,jc,kc)=q2(ic,jc,kc) + fjl(ic,jc)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       call ntrvb(amjl,acjl,apjl,fjl,m1,m2,n2,n1m)
c
c
c       do jc=1,n2
c         do ic=1,n1m
c           q2(ic,jc,kc)=q2(ic,jc,kc) + fjl(ic,jc)
c         end do
c       end do
c     end do
c       do jc=1,n2
c         do ic=1,n1m
c           iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
c           rhs(iadd) = fjl(ic,jc)
c         end do
c       end do
c     end do
      return
      end
c
c
************************************************************************
c   this subroutine performs the inversion of the q2 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q2
c   are treated implicitly
c   in the first part the rhs is calculated
c        direction x3
c
      subroutine solq2k
      include'param.f'
      dimension amkl(m3),apkl(m3),ackl(m3),fkl(m1,m3)
c
c  ************ compute dq2 sweeping along the x3 direction
c               periodic
c
      betadx=beta*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,n2,am2k,ac2k,ap2k,fvisc3)
!$OMP$ SHARED(q2,betadx,rhs)
!$OMP$ PRIVATE(amkl,apkl,ackl,fkl,iadd,faci,ic,ii,jc,kc)
      do jc=2,n2m
        do kc=1,n3m
          faci=1.-ac2k(jc,kc)*betadx*fvisc3(jc,kc)
          amkl(kc)=-am2k(jc,kc)*betadx*fvisc3(jc,kc)/faci
          ackl(kc)=1.
          apkl(kc)=-ap2k(jc,kc)*betadx*fvisc3(jc,kc)/faci
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2
            fkl(ic,kc)=rhs(iadd)/faci
          end do
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        do kc=2,n3m
          ackl(kc)=ackl(kc)*ackl(kc-1)-apkl(kc-1)*amkl(kc)
          apkl(kc)=apkl(kc)*ackl(kc-1)
          do ic=1,n1m
            fkl(ic,kc)=fkl(ic,kc)*ackl(kc-1)-fkl(ic,kc-1)*amkl(kc)
          end do
        end do
c
        do ic=1,n1m
          fkl(ic,n3m)=fkl(ic,n3m)/ackl(n3m)
          iadd = ic+(jc-1)*n1m+(n3m-1)*n1m*n2
          rhs(iadd)= fkl(ic,n3m)
        end do
c
        nm=n3m-1
        do ii=1,nm
          kc=n3m-ii
          do ic=1,n1m
            iadd = ic+(jc-1)*n1m+(kc-1)*n1m*n2
            fkl(ic,kc)= (fkl(ic,kc)-apkl(kc)*fkl(ic,kc+1))/ackl(kc)
            rhs(iadd)= fkl(ic,kc)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     call ntrvb(amkl,ackl,apkl,fkl,m1,m3,n3m,n1m)
c
c       do kc=1,n3m
c         do ic=1,n1m
c           q2(ic,jc,kc)=q2(ic,jc,kc) + fkl(ic,kc)
c         end do
c       end do
c     end do
c
      return
      end
c
c
c
************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x1
c
      subroutine solq3i
      include'param.f'
      dimension amil(m2,m1),acil(m2,m1),apil(m2,m1),fil(m2,m1)
c
c  ************ compute dq3** sweeping along the x1 direction
c               inflow outflow
c
      betadx=beta*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,rhs)
!$OMP$ SHARED(usrm,udx1qm,imv,ipv,usinpc,fvisc1)
!$OMP$ PRIVATE(amil,apil,acil,fil,iadd,faci)
!$OMP$ PRIVATE(ugmmm,ugm2,im,ip,km,ic,jc,kc)
      do kc=2,n3m
          do ic=1,n1m
            im = imv(ic)
            ip = ipv(ic)
        do jc=1,n2m
          ugmmm=betadx*udx1qm(jc)*usinpc(kc)**2
            faci=1.+ugmmm*2.*fvisc1(jc,kc)
            apil(jc,ic)=-ugmmm*fvisc1(jc,kc)/faci
            amil(jc,ic)=-ugmmm*fvisc1(jc,kc)/faci
            acil(jc,ic)=1.
            iadd = ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fil(jc,ic)=rhs(iadd)/faci
          end do
        end do
c
        call tripvmy(amil,acil,apil,fil,1,n1m,1,n2m,m1,m2)
c
        do ic=1,n1m
          do jc=1,n2m
            iadd = ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            rhs(iadd)=fil(jc,ic)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
      return
      end
c
************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x2
c
      subroutine solq3j
      include'param.f'
      dimension amjl(m2),apjl(m2),acjl(m2),fjl(m1,m2)
c
c  ************ compute dq3 sweeping along the x2 direction
c
      betadx=beta*al
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,betadx,rhs,ap3j,am3j,ac3j,fvisc2)
!$OMP$ PRIVATE(amjl,apjl,acjl,fjl,iadd,faci,nm,ic,ii,jc,kc)
      do kc=2,n3m
        do jc=1,n2m
          faci=1.-ac3j(jc)*betadx*fvisc2(jc,kc)
          apjl(jc)=-ap3j(jc)*betadx*fvisc2(jc,kc)/faci
          amjl(jc)=-am3j(jc)*betadx*fvisc2(jc,kc)/faci
          acjl(jc)=1.
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fjl(ic,jc)=rhs(iadd)/faci
          end do
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        do jc=2,n2m
          acjl(jc)=acjl(jc)*acjl(jc-1)-apjl(jc-1)*amjl(jc)
          apjl(jc)=apjl(jc)*acjl(jc-1)
          do ic=1,n1m
            fjl(ic,jc)=fjl(ic,jc)*acjl(jc-1)-fjl(ic,jc-1)*amjl(jc)
          end do
        end do
c
        do ic=1,n1m
          fjl(ic,n2m)=fjl(ic,n2m)/acjl(n2m)
          iadd=ic+(n2m-1)*n1m+(kc-1)*n1m*n2m
          rhs(iadd)= fjl(ic,n2m)
        end do
c
        nm=n2m-1
        do ii=1,nm
          jc=n2m-ii
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fjl(ic,jc)=(fjl(ic,jc)-apjl(jc)*fjl(ic,jc+1))/acjl(jc)
            rhs(iadd)=fjl(ic,jc)
          end do
        end do
      end do
!$OMP  END PARALLEL DO
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       call ntrvb(amjl,acjl,apjl,fjl,m1,m2,n2m,n1m)
c
c       do jc=1,n2m
c         do ic=1,n1m
c           iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
c           rhs(iadd) = fjl(ic,jc)
c         end do
c       end do
c     end do
      return
      end
c
c
************************************************************************
c   this subroutine performs the inversion of the q3 momentum equation
c   by a factored implicit scheme, only the derivatives 11,22,33 of q3
c   are treated implicitly
c       direction x3
c
      subroutine solq3k
      include'param.f'
      dimension amkl(m3),apkl(m3),ackl(m3),fkl(m1,m3)
c  ********* compute the dq3* sweeping in the x3 direction
c
c
      betadx=beta*al
c
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3m,n2m,n1m,n3)
!$OMP$ SHARED(q3,betadx,rhs,am3k,ac3k,ap3k,fvisc3)
!$OMP$ PRIVATE(amkl,apkl,ackl,fkl,iadd,faci,nm,ic,ii,jc,kc)
      do jc=1,n2m
        do kc=2,n3m
          faci=1.-ac3k(jc,kc)*betadx*fvisc3(jc,kc)
          amkl(kc)=-am3k(jc,kc)*betadx*fvisc3(jc,kc)/faci
          ackl(kc)=1.
          apkl(kc)=-ap3k(jc,kc)*betadx*fvisc3(jc,kc)/faci
          do ic=1,n1m
            iadd=ic+(jc-1)*n1m+(kc-1)*n1m*n2m
            fkl(ic,kc)=rhs(iadd)/faci
          end do
        end do
c
        amkl(1)=0.
        apkl(1)=0.
        ackl(1)=1.
        amkl(n3)=0.
        apkl(n3)=0.
        ackl(n3)=1.
        do ic=1,n1m
          fkl(ic,1)= 0.
          fkl(ic,n3)= 0.
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        do kc=2,n3
          ackl(kc)=ackl(kc)*ackl(kc-1)-apkl(kc-1)*amkl(kc)
          apkl(kc)=apkl(kc)*ackl(kc-1)
          do ic=1,n1m
            fkl(ic,kc)=fkl(ic,kc)*ackl(kc-1)-fkl(ic,kc-1)*amkl(kc)
          end do
        end do
c
        do ic=1,n1m
          fkl(ic,n3)=fkl(ic,n3)/ackl(n3)
        end do
c
        nm=n3-1
        do ii=1,nm
          kc=n3-ii
          do ic=1,n1m
            fkl(ic,kc)=(fkl(ic,kc)-apkl(kc)*fkl(ic,kc+1))/ackl(kc)
          end do
        end do
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       call ntrvb(amkl,ackl,apkl,fkl,m1,m3,n3,n1m)
c
        do kc=1,n3
          do ic=1,n1m
            q3(ic,jc,kc)=q3(ic,jc,kc) + fkl(ic,kc)
          end do
        end do
      end do
c
!$OMP  END PARALLEL DO
      return
      end
c
c  ****************************** subrout ntrvb  **********************
c
c     call ntrvb(amkl,ackl,apkl,fkl,m1,m3,n3m,n1m)
      subroutine ntrvb(am,ac,ap,f,mm,nn,n,m)
c
      dimension am(nn),ac(nn),ap(nn),f(mm,nn)
c
c  ******** reduction of trid. matrix to an upper rigth matrix
c
      do 1 i=2,n
      ac(i)=ac(i)*ac(i-1)-ap(i-1)*am(i)
      ap(i)=ap(i)*ac(i-1)
      do 2 k=1,m
      f(k,i)=f(k,i)*ac(i-1)-f(k,i-1)*am(i)
    2 continue
    1 continue
c  ******** calculation of the unknown by backward elimination
c
      do 3 k=1,m
      f(k,n)=f(k,n)/ac(n)
    3 continue
      nm=n-1
      do 10 ii=1,nm
      i=n-ii
      do 11 k=1,m
      f(k,i)=(f(k,i)-ap(i)*f(k,i+1))/ac(i)
   11 continue
   10 continue
      return
      end

c
c
c
c     it inverts in 1 direction, vectors of n1m length
c     from 1, n2m
c
c
c
c
c***********************************************************************
      subroutine tripvmy(ami,aci,api,rrr,n1i,n1f,n2i,n2f,m1,m2)
c
c     vectorized for right hand side and coefficients
c
      dimension fn(m2),p(m2)
      dimension ami(m2,m1),aci(m2,m1),api(m2,m1),rrr(m2,m1)
      dimension q(m2,m1),s(m2,m1),fei(m2,m1)
      ia = n1i + 1
      ii = n1i + n1f
c
c   COEFFICIENTS FOR TRIDIAGONAL INVERSION
c

c
c  THE INVERSION STARTS
c
      do 20 j=n2i,n2f
        q(j,n1i) = -api(j,n1i)/aci(j,n1i)
        s(j,n1i) = -ami(j,n1i)/aci(j,n1i)
        fn(j) = rrr(j,n1f)
        rrr(j,n1i) = rrr(j,n1i)/aci(j,n1i)
   20 continue
c
c     forward elimination sweep
c
      do 10 i=ia,n1f
        do 21 j=n2i,n2f
          p(j) =1./( aci(j,i) + ami(j,i)*q(j,i-1))
          q(j,i) = - api(j,i)*p(j)
          s(j,i) = - ami(j,i)*s(j,i-1)*p(j)
          rrr(j,i) = ( rrr(j,i) - ami(j,i)*rrr(j,i-1))*p(j)
   21   continue
   10 continue
c
c     backward pass
c
      do 22 j=n2i,n2f
        s(j,n1f) = 1.
        fei(j,n1f) = 0.
   22 continue
      do 11 l=ia,n1f
        i = ii - l
        do 23 j=n2i,n2f
          s(j,i) = s(j,i) + q(j,i)*s(j,i+1)
          fei(j,i) = rrr(j,i) + q(j,i)*fei(j,i+1)
   23   continue
   11 continue
      do 24 j=n2i,n2f
        rrr(j,n1f)=(fn(j)-api(j,i)*fei(j,n1i) -
     %      ami(j,i)*fei(j,n1f-1))/(api(j,i)*s(j,n1i) +
     %      ami(j,i)*s(j,n1f-1)+aci(j,i))
   24 continue
c
c     backward elimination pass
c
      do 12 l=ia,n1f
        i = ii -l
        do 25 j=n2i,n2f
          rrr(j,i) = rrr(j,n1f)*s(j,i) + fei(j,i)
   25   continue
   12 continue
      return
      end



*************************************************************************
*                                                                       *
*                         Vortice di Hill                               *
*                                                                       *
*************************************************************************


      subroutine inqHill
      include'param.f'
      pi=2.*asin(1.)
!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2,n1,q1,q2,q3)
!$OMP$ PRIVATE(i,j,k)
      do 80 k=1,n3
         do 80 j=1,n2
            do 80 i=1,n1
              q1(i,j,k)=0.
              q2(i,j,k)=0.
              q3(i,j,k)=0.
  80  continue
!$OMP  END PARALLEL DO

      pi = 2.*asin(1.)
      r0 = 0.1 ! 1.00 per 6test   !0.1 per 5test
      th0 = 0.*pi/4.
      ph0 = 2.*pi/4. !1.*pi/6. per 6test !cosi solo asse, 2.*pi/4. passa per il centro  5test
      x0 = r0*cos(th0)*sin(ph0)
      y0 = r0*sin(th0)*sin(ph0)
      z0 = r0*cos(ph0)

C a1 = rotazione sistema riferimento attorno a versore x
C a2 = rotazione sistema riferimento attorno a versore y
C a3 = rotazione sistema riferimento attorno a versore z

      a1=-2.*pi/4.
      a2=-2.*pi/4.
      a3=pi/2.

C L const01 = u_0. const02=a^2

      const01 = 0.5
      radHi=0.5 !0.2
      const02 = radHi**2

!$OMP  PARALLEL DO
!$OMP$ SHARED(n3,n2m,n1m,q1,q2,q3,rc,thetac,zz)
!$OMP$ SHARED(x0,y0,z0,a1,a2,a3,radHi)
!$OMP$ PRIVATE(i,j,k,xt,yt,zt,xp,yp,zp,Ta,Tb,Tc)
!$OMP$ PRIVATE(xTilde,yTilde,zTilde,rho,raggio,h)
!$OMP$ PRIVATE(vh,vraggio,val01,psi,vx1,vy1,vz1)
!$OMP$ PRIVATE(VxTilde,VyTilde,VzTilde,TaV,TbV,TcV)
      do 810 k=1,n3
         do 810 j=1,n2m
            do 810 i=1,n1m
               xt = rc(j)*cos(thetac(i))*sin(zz(k))
               yt = rc(j)*sin(thetac(i))*sin(zz(k))
               zt = rc(j)*cos(zz(k))

               xp=xt-x0
               yp=yt-y0
               zp=zt-z0


              Ta=cos(a3)*xp-sin(a3)*yp
              Tb=sin(a3)*xp+cos(a3)*yp
              Tc=cos(a2)*zp-sin(a2)*Ta

              xTilde=cos(a2)*Ta+sin(a2)*zp
              yTilde=cos(a1)*Tb-sin(a1)*Tc
              zTilde=sin(a1)*Tb+cos(a1)*Tc

               rho = sqrt((xTilde)**2+(yTilde)**2+(zTilde)**2)
               raggio = sqrt((xTilde)**2+(yTilde)**2)
               h=zTilde

              if(rho.le.radHi) then
                vh=1.5*const01*(1-(2*raggio**2+h**2)/const02)+const01
                vraggio=1.5*const01*(h*raggio)/const02

                !vh=7
                !vraggio= 7
              else
                val01=((const02/(h**2+raggio**2)))**2.5
                vh=const01*(val01*(2*h**2-raggio**2)/(2*const02)) !qui sarebbe bello poter togliere il -1...
                vraggio=1.5*const01*(h*raggio)/const02*val01

                !vh=0
                !vraggio=0
              end if

              psi = atan2(yTilde,xTilde)
                            VxTilde = vraggio*cos(psi)
                            VyTilde = vraggio*sin(psi)
                            VzTilde = vh

              TaV=cos(a1)*VzTilde-sin(a1)*VyTilde
              TbV=sin(a1)*VzTilde+cos(a1)*VyTilde
              TcV=cos(a2)*VxTilde-sin(a2)*TaV

              vx1=sin(a3)*TbV+cos(a3)*TcV
              vy1=cos(a3)*TbV-sin(a3)*TcV
              vz1=cos(a2)*TaV+sin(a2)*VxTilde

              if(k.ne.n3) then
              q1(i,j,k)=-vx1*sin(thetac(i))
     %                  +vy1*cos(thetac(i))
              q2(i,j,k)=(vx1*cos(thetam(i))*sin(zm(k))
     %                  +vy1*sin(thetam(i))*sin(zm(k))
     %                  +vz1*cos(zm(k))
     %                  )*rc(j)**2
               end if
              q3(i,j,k)=(vx1*cos(thetam(i))*cos(zz(k))
     %                  +vy1*sin(thetam(i))*cos(zz(k))
     %                  -vz1*sin(zz(k))
     %                  )*sin(zz(k))

  810 continue
!$OMP  END PARALLEL DO
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c      call outcoordinate
      open(69,file='CXY0.dat',form='formatted',status='unknown')
      rewind(69)
      write(69,*)'Variables = X, Y, Z, Vx, Vy, Vz, Pr, r, theta, phi'
      write(69,*)'ZONE F=POINT, I= ',n1,', J= ',n2m, ', K= ',1
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
         isy = isym(ic)
         ksy = ksym(kc)
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
          vp = 0.5*(q3(ic,jc,kc)+q3(ic,jc,kp))/sin(zm(kc))
          ct = cos(thetam(ic))
          st = sin(thetam(ic))
          cp = cos(zm(kc))
          sp = sin(zm(kc))
                      write(69,*)xxx,yyy,zzz
     %           ,vr*ct*sp-vt*st*sp+vp*ct*cp
     %           ,vr*st*sp+vt*ct*sp+vp*st*cp
     %           ,vr*cp-vp*sp
     %           ,pr(ic,jc,kc)
            end do
          end do
      close(69)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c      call divgck(qmax,qtot)
c
cc
      return
      end
c
c***********************************************************************
c                                                                      *
c                       CONDIZIONI AL CONTORNO                         *
c                                                                      *
c***********************************************************************
      subroutine densbo
      include'param.f'
         do 810 k=1,n3
            do 810 i=1,n1
              denbn(i,k)=0.         ! temp BC on ext. sphere
              denbs(i,k)=1.         ! temp BC on int. sphere
  810 continue
      return
      end
C

