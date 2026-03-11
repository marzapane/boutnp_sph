c***********************************************************************
c                                                                      *
c            ************** SPECT  ************************            *
c                                                                      *
c***********************************************************************
      subroutine init_spec
      include 'param.f'
      integer n1mh

      call fftv
      n1mh=n1m/2+1

      do i=1,n1mh
        do k=1,n3m
          mode(k,i)=real(i-1)/sin(zm(k))
        end do
      end do
      nspec=0
      do i = 1, n1md
        do k = 1, n3m
          spec(k,i) = 0.0d0
        end do
      end do
      end subroutine init_spec
c
c***********************************************************************
c
      subroutine spect
      include 'param.f'
      dimension rwc(m2),rwm(m2),pwc(m3),pwm(m3)
      real      rtr,rsc,rsm,pi
      real   xxr(m1m+2,m3m),wwork(m1m+1,m3m)
      integer n1mh

      call fftv
      dvol = 1./(dx2 * dx3)
      pi = 2. * asin (1.)
      n1mh=n1m/2+1
      n1md=n1m+2
      rsc = 0.0d0
      rsm = 0.0d0
      do j = 1,n2
        rwc(j) = g2rc(j)/dx2
c r^2 factor simplifies between integration measure (r^2 dr)
c and field (v_r = q2/r^2)
        rsc = rsc + rc(j)*rc(j) * rwc(j)
        rwm(j) = rm(j)*rm(j) * g2rm(j)/dx2
        rsm = rsm + rwm(j)
      end do
      do k = 1,n3
        pwc(k) = g3rc(k)/dx3
c sin(phi) factor simplifies between integration measure (sin(phi) dphi)
c and field (v_phi = q3/sin(phi))
        pwm(k) = sinpm(k) * g3rm(k)/dx3
      end do
c average v_r in xxr
      do j = 1,n2
        do k = 1,n3m
          do i = 1,n1md
            xxr(i,k) = 0.0d0
          end do
          xxr(1,k) = q2(n1m,j,k)/rc(j)**2
          do i = 1,n1m
            is = i+1
            xxr(is,k) = q2(i,j,k)/rc(j)**2
          end do
          xxr(n1md,k) = q2(1,j,k)/rc(j)**2
        end do
c
c     fft applied to the (q_2) by fft99
c   from physical to wave number space
c
        call fft99(xxr,wwork,trigxv,ifxv,1,m1+1,n1m,n3m,-1)
        do 31 i=1,n1mh
          ip=2*i
          id=2*i-1
          do 41 k=1,n3m
            spec(k,i) = spec(k,i)
     1                + (xxr(id,k)**2 + xxr(ip,k)**2)
     2                * rc(j)**2 * rwc(j) * pwm(k)
  41      continue
  31    continue
      end do
c
c
      do 2 j=1,n2m
        do 12 k=1,n3m
          do i = 1,n1md
            xxr(i,k) = 0.0d0
          end do
          xxr(1,k) = q1(n1m,j,k)
          do 22 i=1,n1m
            is=i+1
            xxr(is,k) = q1(i,j,k)
  22      continue
          xxr(n1md,k) = q1(1,j,k)
  12    continue
c
c     fft applied to the (q_1) by fft99
c   from physical to wave number space
c
        call fft99(xxr,wwork,trigxv,ifxv,1,m1+1,n1m,n3m,-1)
        do 32 i=1,n1mh
          ip=2*i
          id=2*i-1
          do 42 k=1,n3m
            spec(k,i) = spec(k,i)
     1                + (xxr(id,k)**2 + xxr(ip,k)**2)
     2                * rwm(j) * pwm(k)
  42      continue
  32    continue
  2   continue
c
c
      do 3 j=1,n2m
        do 13 k=1,n3m
          do i = 1,n1md
            xxr(i,k) = 0.0d0
          end do
          xxr(1,k) = q3(n1m,j,k)*usinpc(k)
          do 23 i=1,n1m
            is=i+1
            xxr(is,k) = q3(i,j,k)*usinpc(k)
   23     continue
          xxr(n1md,k) = q3(1,j,k)*usinpc(k)
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
            spec(k,i) = spec(k,i)
     1                + (xxr(id,k)**2 + xxr(ip,k)**2)
     2                * rwm(j) * pwc(k) * sinpc(k)
  43      continue
  33    continue
   3  continue
      
      nspec = nspec + 1

      end subroutine spect
c
c***********************************************************************
c
      subroutine print_spec(time)
      include 'param.f'
      character*30 filsp
      integer n1mh

      n1mh=n1m/2+1
      write(filsp,'(A5,I6.6,A4)') 'spctr', int(10.*time), '.bin'
      open(39,file=filsp,form='unformatted',status='unknown')
      write(39) n3m, n1mh
      write(39) (((mode(k,i)),k=1,n3m),i=1,n1mh)
      write(39) (((spec(k,i)/nspec),k=1,n3m),i=1,n1mh)
      write(39) time
      close(39)
      nspec=0
      do i = 1, n1md
        do k = 1, n3m
          spec(k,i) = 0.0d0
        end do
      end do
      end subroutine print_spec
