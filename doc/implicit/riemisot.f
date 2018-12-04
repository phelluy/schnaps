      program tracriem
      implicit double precision (a-h,o-z)
      parameter(c=0.6,p0=0,rhoa=0.d0,rhow=0.d0)

      open(11,file='rho_u.dat')

      xmin = -0.5
      xmax = 0.5
      tmax = 0.4
      np = 10000
      dx = (xmax-xmin)/(np-1)

      rg = 2
      pg = c*c*rg
      rd = 1
      pd = c*c*rd
      ug = 0.d0
      ud = 0.d0


      do i=1,np
         x = (i-1)*dx + xmin
         xi = x/tmax
         call  riemisot(rg,ug,pg,rd,ud,pd,xi,r,u,p)
         write(11,*) x,r,u
      enddo
      close(11)

      end



c=======================================================================
      subroutine riemisot(rg,ug,pg,rd,ud,pd,xi,r,u,p)
c_______________________________________________________________________
c
c Résolution du problemen de Riemann Isotherme
c=======================================================================
      implicit double precision (a-h,o-z)

c c = Vitesse du son 
c coefficients de la loi d'etat du type isotherme
c p=p0+c**2*(rho- (phi*rhoa+(1-phi)*rhow))
c c: vitesse du son
c rhoa, rhow: masses volumiques de l'air et de l'eau
c p0: pression de reference
      parameter(c=0.6,p0=0,rhoa=0.d0,rhow=0)
c initialisations

      phig=phif(rg,pg)
      phid=phif(rd,pd)
      if(rd.le.0.d0)write(*,*)'rd',rd
      if(rg.le.0.d0)write(*,*)'rg',rg
      phig=0
      phid=0
      
c pression minimale
c (rho=0 dans l'air)
c      fmin=phig
c      if (phid.lt.fmin) fmin=phid
c      rmin=rhof(p0,fmin)

      pming=pf(0.d0,phig)
      pmind=pf(0.d0,phid)
      write(*,*) 'pmin',pming,pmind
      if(pming.le.pmind)then
        pmin=pmind
      else
        pmin=pming
      endif

      r0g=rhof(p0,phig)
      r0d=rhof(p0,phid)
      write(*,*) 'rho0 ',r0g,r0d
      drg=rg-r0g
      drd=rd-r0d

c initialisation par le solveur acoustique
      p=(ug*rg*rd+c*drg*rd-ud*rg*rd+c*drd*rg)/(rd+rg)*c+p0
      write(*,*) 'start p=',p

      dp=1
      ff=1
      iter=1
      itermax=100

      eps=1.d-10

      rgini=rg
      ugini=ug
      pgini=pg
      rdini=rd
      udini=ud
      pdini=pd

      do while(dabs(dp).gt.eps.and.iter.lt.itermax)

         if (p.lt.pmin) p=pmin+eps
         iter=iter+1

         r2=rhof(p,phid)
         r1=rhof(p,phig)
	 write(*,*) 'r1=',r1,' r2=',r2
       if(r1.le.0.d0)write(*,*)r1,p,phig
       if(r2.le.0.d0)write(*,*)r2,p,phid
         ff=ud-ug-h(rd,r2)-h(rg,r1)

         df=-dhdr(rd,r2)-dhdr(rg,r1)
         df=df/c**2

         dp=-ff/df

         p=p+dp

      end do

      if (iter.eq.itermax) then
         write(*,*) 'pas de convergence',dabs(dp/p0)
         write(*,*) 'rg,ug,pg',rg,ug,pg,phif(rg,pg)
         write(*,*) 'rd,ud,pd',rd,ud,pd,phif(rd,pd)
         write(*,*) 'rgi,ugi,pgi',rgini,ugini,pgini,phif(rgini,pgini)
         write(*,*) 'rdi,udi,pdi',rdini,udini,pdini,phif(rdini,pdini)
        stop
      endif

c intermediate densities and velocity
      r2=rhof(p,phid)
      r1=rhof(p,phig)
      um=ug+h(rg,r1)
c vitesse d'interface
      ui=um
      pm=p
c      u=ud-h(rd,r2)

c wave velocities

c 1-shock or 1-rarefaction
      if (r1.gt.rg) then
         s1m=ug-c*dsqrt(r1/rg)
         s1p=um-c*dsqrt(rg/r1)
      else
         s1m=ug-c
         s1p=um-c
      endif

c 2-shock or 2-rarefaction
      if (r2.gt.rd) then
         s2m=ud+c*dsqrt(r2/rd)
         s2p=um+c*dsqrt(rd/r2)
      else
         s2m=um+c
         s2p=ud+c
      endif

c computation of the solution at xi=x/t

      if (xi.lt.s1m) then
         r=rg
         u=ug
         phi=phig
      else if (xi.lt.s1p) then
         u=xi+c
         r=rg*exp((ug-u)/c)
         phi=phig
      else if (xi.lt.um) then
         r=r1
         u=um
         phi=phig
      else if (xi.lt.s2m) then
         r=r2
         u=um
         phi=phid
      else if (xi.lt.s2p) then
         u=xi-c
         r=rd*exp((u-ud)/c)
         phi=phid
      else
         r=rd
         u=ud
         phi=phid
      endif

      p=pf(r,phi)


      return
      end

c=======================================================================
      function dhdr(rL,r)
c=======================================================================
      implicit double precision (a-h,o-z)

c      include 'riemparam.f'
      parameter(c=0.6,p0=0,rhoa=0.d0,rhow=0)

       if(r*rl.le.0.d0)then
         write(*,*)'pb fonction dh'
	 stop
      endif

      if (r.gt.rL) then
         t2 = dsqrt(rL*r)
         t7 = t2**2
         dhdr = -c/t2-c*(rL-r)/t7/t2*rL/2.D0
      else
         dhdr=-c/r
      endif

      return
      end

c=======================================================================
      function h(rL,r)
c=======================================================================
      implicit double precision (a-h,o-z)

c      include 'riemparam.f'
      parameter(c=0.6,p0=0,rhoa=0.d0,rhow=0)


       if(r*rl.le.0.d0)then
         write(*,*)'pb fonction h',r,rl
	 stop
      endif
      if (r.gt.rL) then
         h=c*(rL-r)/dsqrt(rL*r)
      else
         h=c*log(rL/r)
      endif

      return
      end

c=======================================================================
      function rhof(p,phi)
c=======================================================================
      implicit double precision (a-h,o-z)

c      include 'riemparam.f'
      parameter(c=0.6,p0=0,rhoa=0.d0,rhow=0)



      rhof=(p-p0)/c**2+phi*rhoa+(1.d0-phi)*rhow

      return
      end

c=======================================================================
      function pf(rho,phi)
c=======================================================================
      implicit double precision (a-h,o-z)

c      include 'riemparam.f'
      parameter(c=0.6,p0=0,rhoa=0.d0,rhow=0)


      pf=p0+c**2*(rho-(phi*rhoa+(1.d0-phi)*rhow))

      return
      end

c=======================================================================
      function phif(rho,p)
c=======================================================================
      implicit double precision (a-h,o-z)

c      include 'riemparam.f'
      parameter(c=0.6,p0=0,rhoa=0.d0,rhow=0)


      phif=((p-p0)/c**2+rhow-rho)/(rhow-rhoa)

      return
      end
