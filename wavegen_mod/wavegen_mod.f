c compile with gfortran -o wavegen_mod wavegen_mod.f
      program atom
c
c Self-consistent solution of Kohn-Sham equation for an atom
c with a scalar-relativistic Hamilton operator
c
c   Based on a program written by
c   D. R. Hamann
c
c
c
      implicit none
      double precision al, alam, amesh, csc, csc1, dmax1, drl
      double precision e, emax, ep, et, f, fp, rc, rcmax, rct
      double precision rpk, sf, vsc, vsc1, z, zion,fs
      double precision bplus,bminus,btotal,es,diffa,rr,ltest
      integer i,iexc,iexct,iprint, it,l,l1,lmax,lt,max,mch
      integer mmax, n, nc, nin, np, nrl, nv
      character*3 chr
c
      common/par/z,amesh,al,mmax,nin,mch,iexc
      common/besetz/bplus,bminus,btotal,fs(30,2)
      dimension n(30),l(30),f(30),e(30),rpk(30)
      dimension np(5),fp(5),ep(5),rc(5),es(30,2)
c
c     Rydberg-Konstante
c
      rr=13.60569193
c
c********************************************************************
c basic input data
c  LDA   ( or GGA)   type of DFT functional
c  z       atomic number
c  n(i)    principal quantum number
c  l(i)    angular momentum
c  f(i)    occupancy (non-integral ok) - must sum to z or less
c********************************************************************
c
      do 111 i=1,30
      n(i)=0
      l(i)=0
      f(i)=0
      fs(i,1)=0
      fs(i,2)=0
 111  continue
      open(7,file='wavegen.dat',status='old')
      rewind 7
      write(6,*)
      write(6,*)
      read(7,51)chr
 51   format(a3)
      iexct=0
      if(chr.eq.'LDA') iexct=3
      if(chr.eq.'GGA') iexct=4
      if(iexct.eq.0) then
      write(6,*)  ' Mode not LDA or GGA'
      stop
      endif
      read (7,*) z
      nc=0
      if(iexct.eq.3) write(6,*)'LDA calculation of atom with z=',z
      if(iexct.eq.4) write(6,*)'GGA calculation of atom with z=',z
      zion=z
      sf=0.0d0
      lmax=0
      do 10 i=1,30
        read (7,*,end=300) n(i),l(i),fs(i,1),fs(i,2)
        f(i)=fs(i,1)+fs(i,2)
        ltest= 2 *l(i) +1+ 0.000001d0
        if(fs(i,1).gt.ltest.or.fs(i,2).gt.ltest) stop
     & 'Occupation is not allowed'
        if(i .eq. 1) go to 8
        if(n(i) .lt. n(i-1)) stop 'data: state order wrong'
    8   if(n(i) .le. l(i)) stop 'data: l .le. n'
        if(l(i) .gt. 3) stop 'data: l .gt. 3'
        lmax=max(l(i),lmax)
        sf=sf+f(i)
c       if(sf .gt. z) stop 'data: negative ions not allowed'
   10   continue
        write(6,*)'Error: more than 30 shells'
        stop
c
 300    nv=i-1
        close(7)
        do 1200 i=1,nc+nv
        bplus=bplus+fs(i,1)
        bminus=bminus+fs(i,2)
 1200   continue
        btotal=bplus+bminus
        Write(6,*)
        Write(6,*)'Occupation'
        write(6,*)'       up         down        total       z atom '
        write(6,501)bplus,bminus,btotal, z
        diffa=abs(bplus-bminus)
        if(diffa.gt.0.00001d0) write(6,*) 'The atom is spin-polarized'
        diffa=abs(btotal-z)
        if(diffa.gt.0.00001d0) write(6,*) 'The atom is not neutral '
  500 format(f12.8,3i6)
  501 format(4f12.3)
  510 format(2i6,f12.8)
c
      if(iexct .ne. 0) then
        iexc=iexct
      else
        iexc=3
      end if
c
c full potential atom solution
c
      call sratom(n,l,e,f,rpk,nc+nv,it,es)
c
c output printing
c     write(6,600)
      write(6,*)
      write(6,620)
      write(6,621)
      do 202 i=1,nc+nv
  202 write(6,630) n(i),l(i),fs(i,1),fs(i,2),es(i,1)*rr*2,es(i,2)*rr*2
      write(6,632) mmax,it
      if(it .ge. 100) stop 'potential not converged'
c
  600 format('1at. no. z   no. core    no. val     exc type    mesh')
  610 format(f6.2,6x,4(i6,6x))
  620 format(' n     l       population           eigen values [eV]')
  621 format('               up   down         up        down'/)
  630 format(2(i2,4x),2f6.2,2f13.3)
  632 format(//'mesh:', i4, ' iterations',i6)
      stop
c
      end
c
c self-consistent scalar-relativistic  all-electron atom
c calculation using log mesh
c
      subroutine sratom(n,l,e,f,rpk,ncv,it,es)
c
      implicit none
      double precision aei, aeo, aii, aio, al, amesh, bl,vis
      double precision dlog, dr, e, eeel, eexc, et, f, r, rho,esum
      double precision rl, rl1, rpk, sd, sf, sn, tfapot, thl, u, uld
      double precision up, vi, vi1, vn, vo, vo1, z, zion,es,vos,fs
      double precision bplus,bminus,btotal,vpseu, rmax
      double precision,allocatable,dimension(:,:,:) :: uall
      integer i, iexc, it, j, l, mch, mmax, n, ncv, nin,ispin
      integer moder
      logical convg
c
      common/arr/r(8192),vi(8192),u(8192),up(8192),vis(8192,2)
      common/par/z,amesh,al,mmax,nin,mch,iexc
      common/besetz/bplus,bminus,btotal,fs(30,2)
      dimension n(30),l(30),f(30),e(30),rpk(30),es(30,2)
      dimension rho(8192,2),vi1(8192,2),vo1(8192,2),vo(8192,2)
      dimension vpseu(3,8192)
	  allocate(uall(8192,30,2))
      bl=0.5d0
c
c     moder =1  :  with relativistic corrections
c     moder =2  :  without relativistic corrections
c
      it=0
      moder=2
c
c    Parameters of mesh
c
c define maximal meshpoints mmax, with r(mmax=4608)=45au
c original contribution for mmax from the program atom.f
c      mmax=dlog(7200.0d0*z)/al
	  mmax=2**13
	  rmax=45.
	  amesh=(rmax*160.*z)**(1./mmax)
c TODO set mmax appropriately (reduce it) after issue for mmax<2**13 has been resolved
c rmax gives the maximum r of the lattice
c mmax<8192 is still possible and reduces computation time greatly.
	  al=dlog(amesh)
	  write(6,*)'Zahl der Netzpunkte (mmax) =',mmax
c
c generate mesh and starting Thomas-Fermi potential
c
      r(1)=.00625d0/z
      do 12 i=2,mmax
        r(i)=amesh*r(i-1)
   12   continue
      do 20 i=1,mmax
        vi(i)=tfapot(r(i),z)
        vis(i,1)=vi(i)
 20     vis(i,2)=vi(i)
c
c starting approximation for energies
      sf=0.0d0
      do 21 i=1,ncv
        sf=sf+f(i)
        zion=z+1.0d0-sf
        e(i)=-0.5d0*(zion/n(i))**2
        if(e(i) .gt. vi(mmax)) e(i)=2.0d0*vi(mmax)
        es(i,1)=e(i)
        es(i,2)=e(i)
   21   continue
c
c return point for self-consistency loop
   24 it=it+1
      convg=.true.
c
c     loop for spin up and spin down
c
      do 1200 ispin=1,2
      do 22 j=1,mmax
   22   rho(j,ispin)=0.0d0
 1200 continue
c
c solve for states in turn
      do 1201 ispin=1,2
      do 1202 j=1,mmax
 1202 vi(j)=vis(j,ispin)
      do 1203 j=1,ncv
 1203 e(j)=es(j,ispin)
      do 40 i=1,ncv
        et=e(i)
c     lww=l(i)+1
        call lschps(moder,n(i),l(i),et,uld)
        if(e(i) .ne. et) convg=.false.
        e(i)=et
c accumulate charge
        do 30 j=1,mmax
        uall(j,i,ispin)=u(j)
   30     rho(j,ispin)=rho(j,ispin) + fs(i,ispin)*(u(j)/r(j))**2
c
c find outermost peak of wavefunction
        do 32 j=nin-1,1,-1
          if(up(j)*up(j+1) .lt. 0.0d0) then
            rpk(i)=r(j)
            go to 34
          end if
   32   continue
   34   continue
c
   40   continue
c
c
c output potential
      do 1204 j=1,ncv
 1204 es(j,ispin)=e(j)
 1201 continue
c
c  Calculation of potential
c    vo = Vion + Vcoul + Vxc
c
      call vout(rho,vo,sf-z,eeel,eexc)
c
c generate next iteration using  anderson's
c method for mixing of old and new potential
      do 1210 ispin=1,2
      thl=0.0d0
      if(it .le. 1) go to 130
      sn=0.0d0
      sd=0.0d0
      do 120 i=1,mmax
        rl=vo(i,ispin)-vis(i,ispin)
        rl1=vo1(i,ispin)-vi1(i,ispin)
        dr=rl-rl1
        sn=sn + rl*dr*r(i)**2
  120   sd=sd + dr*dr*r(i)**2
      thl=sn/sd
  130 do 140 i=1,mmax
        vn=(1.0d0-bl)*((1.0d0-thl)*vis(i,ispin) + thl*vi1(i,ispin))
     &   + bl*((1.0d0-thl)*vo(i,ispin) + thl*vo1(i,ispin))
        vi1(i,ispin)=vis(i,ispin)
        vo1(i,ispin)=vo(i,ispin)
        vis(i,ispin)=vn
  140   continue
 1210   continue
c
      if(.not. convg .and. it .lt. 100) go to 24
c
c     Calculation of coulomb energy and exchange-correlation energy
c
      call vout(rho,vo,sf,eeel,eexc)
      esum=0.0d0
      write(6,*)'Eigenvalues    up   down'
      do 2666 j=1,ncv
2666  write(6,*)es(j,1),es(j,2),es(j,1)*27.2116d0,es(j,2)*27.2116d0
      do 2605 j=1,ncv
      esum=esum+fs(j,1)*es(j,1)+fs(j,2)*es(j,2)
 2605 continue
      write(6,*)
      write(6,*)'Energies:'
      write(6,*)'========='
      write(6,*)'Coulomb energy :',eeel, ' [Hartree]'
      write(6,*)'Exchange-correlation energy:',eexc,  ' [Hartree]'
      write(6,*)'Sum of eigenvalues:',esum, ' [Hartree]'
      esum=esum-0.5d0*eeel+eexc
      write(6,*)'Total energy ',esum,' [Hartree]'
      esum=esum*2.D0
      write(6,*)'Total energy ',esum,' [Rydberg]'
      esum=esum*13.60569193D0
      write(6,*)'Total energy ',esum,' [eV]'
c
c     output of wave funtions
c
      open(13,file='waveup.dat')
c	  output of atomic number
	  write(13,*) 'atomic_number',z
c 	  output of energy in eV
	  write(13,*) 'eigenvalues_eV',(es(i,1)*27.21138386d0,i=1,ncv)
c     output radial wavefunction
	  write(13,*) 'unl(r)=Rnl(r)*r for spin=up in 1/sqrt(au)'
      do 200 j=1,mmax
c     output of u(r)  for spin up
      	write(13,*) r(j),vis(j,1),(uall(j,i,1),i=1,ncv)
 200  continue
      close(13)

      open(13,file='wavedown.dat')
c	  output of atomic number
	  write(13,*) 'atomic_number',z
	  write(13,*) 'eigenvalues_eV',(es(i,2)*27.21138386d0,i=1,ncv)
c     output radial wavefunction
	  write(13,*) 'unl(r)=Rnl(r)*r for spin=down in 1/sqrt(au)'
      do 201 j=1,mmax
c 		output of u(r)  for spin down
      	write(13,*) r(j),vis(j,2),(uall(j,i,2),i=1,ncv)
 201  continue
 203  format(f15.8,f30.8,31f15.8)
 204  format(2A11,31f20.8)
      close(13)

c output of r(j) and vis(j,ispin)
c	  open(13,file='r_v_value.dat')
c		write(13,206) 'number of meshpoints: mmax=',mmax
c		write(13,207) 'mesh spacing: amesh=',amesh
c		write(13,207) 'atomic number: Z=',z
c		write(13,208) 'r(j) in [au]','v(j) in [Hartree] spin=up'
c     &   ,'v(j) in [Hartree] spin=down'

c		do 210 j=1,mmax
c			write(13,205) r(j),vis(j,1),vis(j,2)
c 210	continue
c	  close(13)
 205  format(f15.8,2f30.8)
 206  format(A10,I10)
 207  format(A10,f15.8)
 208  format(A15,2A30)
      return
c
      end
c
c
c lschps
      subroutine lschps(mode,n,l,e,uld)
c integrates radial pauli-type scalar-relativistic equation
c on a logarithmic mesh
c
c mode = 1 is for full potential bound state
c mode = 2  without scalar relativistic term
c
      implicit none
      double precision aei, aeo, aii, aio, al, als, amesh, cf, cn,vis
      double precision dabs, de, dmax1, dmin1, dsqrt, dv, e, emax, emin
      double precision eps, exp, fr, frp, fss, gamma, r, ro, sc
      double precision sls, sn, tfapot, u, uld, uout, up, upin, upout
      double precision upp, v, xkap, z
      integer i, iexc, it, l, mch, mmax,mode,n,nin,nint,node
c
      common/arr/r(8192),v(8192),u(8192),up(8192),vis(8192,2)
      common/par/z,amesh,al,mmax,nin,mch,iexc
      dimension upp(8192),cf(8192),dv(8192),fr(8192),frp(8192)
c
c convergence factor for solution of schroedinger eq.  if calculated
c correction to eigenvalue is smaller in magnitude than eps times
c the magnitude of the current guess, the current guess is not changed.
      eps=1.0d-8
c
c relativistic - non-relativistic switch
c
	  if(mode .eq. 1 .or. mode .eq. 3) then
        fss=(1.0d0/137.036d0)**2
        if(l .eq. 0) gamma=dsqrt(1.0d0-fss*z**2)
        if(l .gt. 0) gamma=(l*dsqrt(l**2-fss*z**2) +
     &   (l+1)*dsqrt((l+1)**2-fss*z**2))/(2*l+1)
      else
c        fss=1.0d-20
		fss=0.0d0
c 	fss is the fine strukture constante fss=0-->
c	Kohn-Sham-equation(fss=0)=Schr\F6dinger-equation
        gamma=l+1
      end if
c
      sls=l*(l+1)
c
      if(mode .eq. 1 .or. mode .eq. 2) then
        emax=v(mmax)+0.5d0*sls/r(mmax)**2
        emin=0.0d0
        do 6 i=1,mmax
    6     emin=dmin1(emin,v(i)+0.5d0*sls/r(i)**2)
        if(e .gt. emax) e=1.25d0*emax
        if(e .lt. emin) e=0.75d0*emin
        if(e .gt. emax) e=0.5d0*(emax+emin)
      else if(mode .eq. 4) then
        emax=e +  10.0d0
        emin=e - 10.0d0
      end if
c
c null arrays to remove leftover garbage
      do 8 i=1,4
        u(i)=0.0d0
        up(i)=0.0d0
    8   upp(i)=0.0d0
      nint=0
      als=al**2
c
c return point for bound state convergence
   10 nint=nint+1
      if(nint .gt. 60) stop 'lschpp convergence error'
c
c coefficient array for u in differential eq.
      do 20 i=1,mmax
   20   cf(i)=als*sls + 2.0d0*als*(v(i)-e)*r(i)**2
c
c calculate dv/dr for darwin correction
      dv(1)=(-50.d0*v(1)+96.d0*v(2)-72.d0*v(3)+32.d0*v(4)
     &       -6.d0*v(5))/(24.d0*al*r(1))
      dv(2)=(-6.d0*v(1)-20.d0*v(2)+36.d0*v(3)-12.d0*v(4)
     &       +2.d0*v(5))/(24.d0*al*r(2))
c
      do 22 i=3,mmax
 22     dv(i)=(2.d0*v(i-2)-16.d0*v(i-1)+16.d0*v(i+1)
     &         -2.d0*v(i+2))/(24.d0*al*r(i))
c
c  relativistic coefficient arrays for u (fr) and up (frp).
      do 24 i=1,mmax
        fr(i)=als*(r(i)**2)*(-fss*(v(i)-e)**2 + 0.5d0*fss*dv(i)/
     &  (r(i)*(1.0d0+0.5d0*fss*(e-v(i)))))
 24     frp(i)=-al*r(i)*0.5d0*fss*dv(i)/(1.0d0+0.5d0*fss*(e-v(i)))
c
c find classical turning point for matching
      if(mode .eq. 1 .or. mode .eq. 2) then
        do 30 i=mmax,2,-1
          if(cf(i-1) .le. 0.d0 .and. cf(i) .gt. 0.d0) then
            mch=i
            go to 40
          end if
   30   continue
        stop 'no classical turning point'
   40   continue
      else
        nin=mch
      end if
c
c start wavefunction with series
c
      do 50 i=1,4
        u(i)=r(i)**gamma
        up(i)=al*gamma*r(i)**gamma
   50   upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
c
c outward integration using predictor once, corrector
c twice
      node=0
c
      do 70 i=4,mch-1
        u(i+1)=u(i)+aeo(up,i)
        up(i+1)=up(i)+aeo(upp,i)
        do 60 it=1,2
          upp(i+1)=(al+frp(i+1))*up(i+1)+(cf(i+1)+fr(i+1))*u(i+1)
          up(i+1)=up(i)+aio(upp,i)
   60     u(i+1)=u(i)+aio(up,i)
        if(u(i+1)*u(i) .le. 0.0d0) node=node+1
   70 continue
c
      uout=u(mch)
      upout=up(mch)
c
c
      if(node-n+l+1 .eq. 0 .or. mode .eq. 3 .or. mode .eq. 5) then
c
        if(mode .eq. 1 .or. mode .eq. 2) then
c
c start inward integration at 10*classical turning
c point with simple exponential
          nin=mch+2.3d0/al
          if(nin+4 .gt. mmax) nin=mmax-4
          xkap=dsqrt(sls/r(nin)**2 + 2.0d0*(v(nin)-e))
c
          do 110 i=nin,nin+4
            u(i)=exp(-xkap*(r(i)-r(nin)))
            up(i)=-r(i)*al*xkap*u(i)
  110       upp(i)=(al+frp(i))*up(i)+(cf(i)+fr(i))*u(i)
c
c integrate inward
c
          do 130 i=nin,mch+1,-1
            u(i-1)=u(i)+aei(up,i)
            up(i-1)=up(i)+aei(upp,i)
            do 130 it=1,2
              upp(i-1)=(al+frp(i-1))*up(i-1)+(cf(i-1)+fr(i-1))*u(i-1)
              up(i-1)=up(i)+aii(upp,i)
  130         u(i-1)=u(i)+aii(up,i)
c
c scale outside wf for continuity
          sc=uout/u(mch)
c
          do 150 i=mch,nin
            up(i)=sc*up(i)
  150       u(i)=sc*u(i)
c
          upin=up(mch)
c
        else
c
          upin=uld*uout
c
        end if
c
c perform normalization sum
c
        ro=r(1)/dsqrt(amesh)
        sn=ro**(2.0d0*gamma+1.0d0)/(2.0d0*gamma+1.0d0)
c
        do 160 i=1,nin-3
  160     sn=sn+al*r(i)*u(i)**2
c
        sn=sn + al*(23.0d0*r(nin-2)*u(nin-2)**2
     &            + 28.0d0*r(nin-1)*u(nin-1)**2
     &            +  9.0d0*r(nin  )*u(nin  )**2)/24.0d0
c
c normalize u
        cn=1.0d0/dsqrt(sn)
        uout=cn*uout
        upout=cn*upout
        upin=cn*upin
c
        do 180 i=1,nin
          up(i)=cn*up(i)
  180     u(i)=cn*u(i)
        do 190 i=nin+1,mmax
  190     u(i)=0.0d0
c
c exit for fixed-energy calculation
c
        if(mode .eq. 3 .or. mode .eq. 5) return

c perturbation theory for energy shift
        de=0.5d0*uout*(upout-upin)/(al*r(mch))
c
c convergence test and possible exit
c
        if(dabs(de) .lt. dmax1(dabs(e),0.2d0)*eps) return
c
        if(de .gt. 0.0d0) then
          emin=e
        else
          emax=e
        end if
        e=e+de
        if(e .gt. emax .or. e .lt. emin) e=0.5d0*(emax+emin)
c
c loop back to converge e
c
        go to 10
c
      else if(node-n+l+1 .lt. 0) then
c too few nodes
        emin=e
        e=0.5d0*(emin+emax)
        go to 10
c
      else
c too many nodes
        emax=e
        e=0.5d0*(emin+emax)
        go to 10
      end if
c
      end
c
c
c
c
c
c adams extrapolation and interpolation formulas for
c outward and inward integration, abramowitz and
c stegun, p. 896
      double precision function aeo(y,j)
c
      double precision y
      integer j
c
      dimension y(8192)
      aeo=(4.16666666667d-2)*(55.0d0*y(j)-59.0d0*y(j-1)
     & +37.0d0*y(j-2)-9.0d0*y(j-3))
      return
      end
c
      double precision function aio(y,j)
c
      double precision y
      integer j
c
      dimension y(8192)
      aio=(4.16666666667d-2)*(9.0d0*y(j+1)+19.0d0*y(j)
     & -5.0d0*y(j-1)+y(j-2))
      return
      end
c
      double precision function aei(y,j)
c
      double precision y
      integer j
c
      dimension y(8192)
      aei=-(4.16666666667d-2)*(55.0d0*y(j)-59.0d0*y(j+1)
     & +37.0d0*y(j+2)-9.0d0*y(j+3))
      return
      end
c
      double precision function aii(y,j)
c
      double precision y
      integer j
c
      dimension y(8192)
      aii=-(4.16666666667d-2)*(9.0d0*y(j-1)+19.0d0*y(j)
     & -5.0d0*y(j+1)+y(j+2))
      return
      end
c
c
c
c tfapot
      double precision function tfapot(r,z)
c generalized thomas fermi atomic potential
c
      double precision b, dsqrt, r, t, x, xs, z
c
      b=(0.69395656d0/z)**.33333333d0
      x=r/b
      xs=dsqrt(x)
c
      t=z/(1.0d0+xs*(0.02747d0 - x*(0.1486d0 - 0.007298d0*x))
     &   + x*(1.243d0 + x*(0.2302d0 + 0.006944d0*x)))
c
      if(t .lt. 1.0d0) t=1.0d0
      tfapot=-t/r
      return
      end

      subroutine vout(rhos,vos,zion,eeel,eexc)
c
c output electrostatic and exchange-correlation potential
c
      double precision aii, al, aln, amesh, datan, den, difxc,rhos, vis
      double precision dlog, dsqrt, ecp, eeel, eexc, exc, foutv, foutx
      double precision fxc, pi4, r, rh, rho, rs, rsl, rsm1, rv,pi,xalf
      double precision rvp, sqrs,tv,u,up,vi,vi1,vo,vo1,vos
      double precision x, z, zion,denp,excp,fxcp,foutxm,difxcm,fs
      double precision bplus,bminus,btotal,xxsi,fxsi,abfxsi,e43,e13
      double precision dritt,alpha,beta,gamma,delta,ep,ef,xf,xp
      integer i, iexc, mch, mmax, nin,contmmax
c
      common/arr/r(8192),vi(8192),u(8192),up(8192),vis(8192,2)
      common/par/z,amesh,al,mmax,nin,mch,iexc,contmmax
      common/besetz/bplus,bminus,btotal,fs(30,2)
      dimension rho(8192),rvp(8192),rv(8192),rhos(8192,2)
      dimension vo(8192),vi1(8192),vo1(8192),vos(8192,2)
      dimension difxc(8192),difxcm(8192)
c
      foutv(i)=rho(i)*vo(i)*r(i)**3
      foutx(i)=rhos(i,1)*difxc(i)*r(i)**3+rhos(i,2)*difxcm(i)*r(i)**3

c
      e43=4.d0/3.d0
      e13=1.d0/3.d0
      pi4=16.0d0*datan(1.0d0)
      pi=4.d0*datan(1.0d0)
c
c      total charge
c
      do 10 i=1,mmax
 10   rho(i)=rhos(i,1)+rhos(i,2)
c
c integration for electrostatic potential
      do 60 i=1,mmax
   60   rvp(i)=rho(i)*al*r(i)**3
c
      rv(mmax)=zion
      rv(mmax-1)=zion
      rv(mmax-2)=zion
c
      do 70 i=mmax-2,2,-1
   70   rv(i-1)=rv(i)+aii(rvp,i)
c
      do 80 i=1,mmax
   80   rvp(i)=rho(i)*al*r(i)**2
c
      tv=0.0d0
      do 90 i=mmax-2,2,-1
        tv=tv+aii(rvp,i)
   90   rv(i-1)=rv(i-1)-r(i-1)*tv
c
      do 100 i=1,mmax
 100    vo(i)=rv(i)/r(i)
c
c electron-electron interaction for total energy
c
      eeel=(9.0d0*foutv(1) + 23.0d0*foutv(2)
     &   + 28.0d0*foutv(3))/28.0d0
      do 102 i=4,mmax
  102   eeel=eeel + foutv(i)
      eeel=al*eeel + foutv(1)/3.0d0
c
      call gradk(rhos,rho,vos,mmax,r,al,eexc,iexc)
      do 888 i=1,mmax
      vos(i,1)=vos(i,1)+vo(i)
      vos(i,2)=vos(i,2)+vo(i)
 888   continue
      i=125
      return
      end
      subroutine  GRADK(RHOS,RHO,VOS,MMAX,R,al,EEXC,iexc)
      parameter(mm=2**13)
c TODO find mmax related issue
c changed from parameter (mm=8192)
c test what happens with mm!=mmax -> seems fine, output is according to mmax
c if mm<2**13 execution fails (for carbon) due to no convergence (independet of mmax, it seems)
c seems to be some sort of extended lattice...?
c mm>2**13 doesn't seem to work
      implicit real*8 (a-h,o-z)
      real*8 kf,ks
      dimension r(mm),rho(mm),rhos(mm,2),vos(mm,2),ea(mm),za(mm),
     &da(mm),x(mm),exx(mm),difxc(mm,2),dz(mm)
      common/perloc/iwahl
      common/gas/kf,ks,g,ec,ecrs,eczet
      grenze=1.d-10
      iwahl=iexc
      pi=4.d0*atan(1.d0)
      pi4=4.d0*pi
      ed=1.d0/3.d0
      zd=2.d0/3.d0
      fak=3.d0/pi4
      fak2=(9.d0/4.d0*pi)**ed
      fak3=4.d0/pi
c
c     Austauschanteil    1/2 E[2nup] + 1/2 E[2ndn]
c
      do 10 j=1,2
      do 1 i=1,mmax
  1   x(i)=2.d0*rhos(i,j)/pi4
      call ezabl(x,mmax,mm,ea,za,r,al)
      do 2 i=1,mmax
 2    x(i)=dabs(ea(i))
      call eabl(x,mmax,mm,da,r,al)
c     do 2 i=1,mmax
c     da(i)=za(i)
c2    if(ea(i).le.0.d0)da(i)=-da(i)
      do 3 i=1,mmax
      d=2.d0*rhos(i,j)/pi4
      if(d.lt.1.d-15) then
      vx=0.d0
      ex=0.0d0
      goto 30
      endif
      rs=(fak/D)**ed
      kf=fak2/rs
      s=x(i)/2.d0/kf/d
      u=ea(i)*DA(i)/D/d/(2.d0*kf)**3
      v=za(i)+2.d0/r(i)*ea(i)
      v=v/(D*(2.d0*kf)**2)
c     write(6,*)r(i),d,s,u,v
      call exch(d,s,u,v,ex,vx)
 30    difxc(i,j)=(ex-vx)*rhos(i,j)
c30    difxc(i,j)=(ex)*rhos(i,j)
      vos(i,j)=vx
c     write(6,*)i,vos(i,j)
 3    continue

 10   continue
      do 4 i=1,mmax
 4    exx(i)=(difxc(i,1)+difxc(i,2))*r(i)**3

      eexc=(9.d0*exx(1)+23.d0*exx(2)+28.d0*exx(3))/28.d0
      do 111 i=4,mmax
 111  eexc=eexc+exx(i)
      eexc=eexc*al+exx(1)/3.d0
c     write(6,*)'Austauschanteil',eexc
      i=125
c     write(6,*)r(i),vos(i,1),vos(i,2),eexc


c
c     Korrelationsanteil
c
      do 15 i=1,mmax
15    x(i)=rho(i)/pi4
      call ezabl(x,mmax,mm,ea,za,r,al)
      do 16 i=1,mmax
      if(rho(i).gt.1.d-18)then
       da(i)=(rhos(i,1)-rhos(i,2))/rho(i)
       else
       da(i)=0.d0
       endif
 16   x(i)=dabs(ea(i))
      call eabl(da,mmax,mm,dz,r,al)
c     call eabl(x,mmax,mm,da,r,al)
      do 17 i=1,mmax
      da(i)=za(i)
 17   if(ea(i).le.0.d0)da(i)=-da(i)
      do 20 i=1,mmax
      d=rho(i)/pi4
      if(d.lt.grenze) then
      ec=00.0d0
      vcup=0.0d0
      vcdn=0.0d0
      h=0.0d0
      dvcup=00.0d0
      dvcdn=0.0d0
      goto 40
      endif
      rs=(fak/d)**ed
      kf=fak2/rs
      ks=dsqrt(fak3*kf)
         if(rho(i).gt.1.d-15) then
         zet=(rhos(i,1)-rhos(i,2))/rho(i)
         else
         zet=0.d0
         endif
      g=((1.d0+zet)**zd+(1.d0-zet)**zd)/2.d0
      gks=2.d0*ks*g
      t=x(i)/(D*gks)
      uu=ea(i)*da(i)/(d**2*gks**3)
      vv=za(i)+2.d0/r(i)*ea(i)
      vv=vv/(D*gks**2)
      ww=ea(i)*dz(i)/(D*gks**2)
      call corlsd(rs,zet,ec,vcup,vcdn,ecrs,eczet,alfc)
      h=0.0d0
      dvcup=0.0d0
      dvcdn=0.0d0
      if(abs(ec).gt.1.d-12) call corgga(rs,zet,t,uu,vv,ww,h,dvcup,dvcdn)
c     write(6,*)'rs,zet,ec,vcup,vcdn',rs,zet,ec,vcup,vcdn
c     write(6,*)i,d,t,uu,vv,ww,h,dvcup,dvcdn,gks
 40   continue
      vos(i,1)=vos(i,1)+dvcup+vcup
      vos(i,2)=vos(i,2)+dvcdn+vcdn
c     if(i.eq.125) write(6,*)'gradp',dvcup, dvcdn
c     difxc(i,1)=(h+ec)*rhos(i,1)
c     difxc(i,2)=(h+ec)*rhos(i,2)
      difxc(i,1)=(h+ec-dvcup-vcup)*rhos(i,1)
      difxc(i,2)=(h+ec-dvcdn-vcdn)*rhos(i,2)
c20   exx(i)=(exx(i)+difxc(i,1)+difxc(i,2))*r(i)**3
 20   exx(i)=(difxc(i,1)+difxc(i,2))*r(i)**3

c
c     Aufsummation
c
      eec=(9.d0*exx(1)+23.d0*exx(2)+28.d0*exx(3))/28.d0
      do 11 i=4,mmax
 11   eec=eec+exx(i)
      eec=eec*al+exx(1)/3.d0
c     write(6,*)'Korrelationsanteil',eec
      eexc=eexc+eec
c
      return
      end

      SUBROUTINE EXCH(D,S,U,V,EX,VX)
C  GGA91 EXCHANGE FOR A SPIN-UNPOLARIZED ELECTRONIC SYSTEM
C  INPUT D : DENSITY
C  INPUT S:  ABS(GRAD D)/(2*KF*D)
C  INPUT U:  (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KF)**3)
C  INPUT V: (LAPLACIAN D)/(D*(2*KF)**2)
C  OUTPUT:  EXCHANGE ENERGY PER ELECTRON (EX) AND POTENTIAL (VX)
      IMPLICIT REAL*8 (A-H,O-Z)
c
c     parameter nach Perdew & Wang PW-GGA II
c
      DATA A1,A2,A3,A4/0.19645D0,0.27430D0,0.15084D0,100.D0/
      DATA AX,A,B1/-0.7385588D0,7.7956D0,0.004D0/
c
c     Parameter nach Becke
c
c     DATA A1,A2,A3,A4/0.1559208d0,0.21770d0,0.0,1.0d0/
c     DATA AX,A,B1/-0.7385588D0,6.1873354d0,0.0d0/
c
      DATA THRD,THRD4/0.333333333333D0,1.33333333333D0/
      common/perloc/iwahl
      FAC = AX*D**THRD
      S2 = S*S
      S3 = S2*S
      S4 = S2*S2
      P0 = 1.D0/DSQRT(1.D0+A*A*S2)
      P1 = DLOG(A*S+1.D0/P0)
      P2 = DEXP(-A4*S2)
      P3 = 1.D0/(1.D0+A1*S*P1+B1*S4)
      P4 = 1.D0+A1*S*P1+(A2-A3*P2)*S2
      F = P3*P4
      EX = FAC*F
C  LOCAL EXCHANGE OPTION
      if(iwahl.eq.3) ex=fac
c     EX = FAC
C  ENERGY DONE. NOW THE POTENTIAL:
      P5 = B1*S2-(A2-A3*P2)
      P6 = A1*S*(P1+A*S*P0)
      P7 = 2.D0*(A2-A3*P2)+2.D0*A3*A4*S2*P2-4.D0*B1*S2*F
      FS = P3*(P3*P5*P6+P7)
      P8 = 2.D0*S*(B1-A3*A4*P2)
      P9 = A1*P1+A*A1*S*P0*(3.D0-A*A*S2*P0*P0)
      P10 = 4.D0*A3*A4*S*P2*(2.D0-A4*S2)-8.D0*B1*S*F-4.D0*B1*S3*FS
      P11 = -P3*P3*(A1*P1+A*A1*S*P0+4.D0*B1*S3)
      FSS = P3*P3*(P5*P9+P6*P8)+2.D0*P3*P5*P6*P11+P3*P10+P7*P11
      VX = FAC*(THRD4*F-(U-THRD4*S3)*FSS-V*FS)
C  LOCAL EXCHANGE OPTION:
      if(iwahl.eq.3) vx=fac*THRD4
c     VX = FAC*THRD4
      RETURN
      END
      SUBROUTINE CORLSD(RS,ZET,EC,VCUP,VCDN,ECRS,ECZET,ALFC)
C  UNIFORM-GAS CORRELATION OF PERDEW AND WANG 1991
C  INPUT: SEITZ RADIUS (RS), RELATIVE SPIN POLARIZATION (ZET)
C  OUTPUT: CORRELATION ENERGY PER ELECTRON (EC), UP- AND DOWN-SPIN
C     POTENTIALS (VCUP,VCDN), DERIVATIVES OF EC WRT RS (ECRS) & ZET (ECZET)
C  OUTPUT: CORRELATION CONTRIBUTION (ALFC) TO THE SPIN STIFFNESS
      IMPLICIT REAL*8 (A-H,O-Z)
      DATA GAM,FZZ/0.5198421D0,1.709921D0/
      DATA THRD,THRD4/0.333333333333D0,1.333333333333D0/
      F = ((1.D0+ZET)**THRD4+(1.D0-ZET)**THRD4-2.D0)/GAM
      CALL GCOR(0.0310907D0,0.21370D0,7.5957D0,3.5876D0,1.6382D0,
     1    0.49294D0,1.00D0,RS,EU,EURS)
      CALL GCOR(0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     1    0.62517D0,1.00D0,RS,EP,EPRS)
      CALL GCOR(0.0168869D0,0.11125D0,10.357D0,3.6231D0,0.88026D0,
     1    0.49671D0,1.00D0,RS,ALFM,ALFRSM)
C  ALFM IS MINUS THE SPIN STIFFNESS ALFC
      ALFC = -ALFM
      Z4 = ZET**4
      EC = EU*(1.D0-F*Z4)+EP*F*Z4-ALFM*F*(1.D0-Z4)/FZZ
C  ENERGY DONE. NOW THE POTENTIAL:
      ECRS = EURS*(1.D0-F*Z4)+EPRS*F*Z4-ALFRSM*F*(1.D0-Z4)/FZZ
      FZ = THRD4*((1.D0+ZET)**THRD-(1.D0-ZET)**THRD)/GAM
      ECZET = 4.D0*(ZET**3)*F*(EP-EU+ALFM/FZZ)+FZ*(Z4*EP-Z4*EU
     1        -(1.D0-Z4)*ALFM/FZZ)
      COMM = EC -RS*ECRS/3.D0-ZET*ECZET
      VCUP = COMM + ECZET
      VCDN = COMM - ECZET
      RETURN
      END
      SUBROUTINE GCOR(A,A1,B1,B2,B3,B4,P,RS,GG,GGRS)
C  CALLED BY SUBROUTINE CORLSD
      IMPLICIT REAL*8 (A-H,O-Z)
      P1 = P + 1.D0
      Q0 = -2.D0*A*(1.D0+A1*RS)
      RS12 = DSQRT(RS)
      RS32 = RS12**3
      RSP = RS**P
      Q1 = 2.D0*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RSP)
      Q2 = DLOG(1.D0+1.D0/Q1)
      GG = Q0*Q2
      Q3 = A*(B1/RS12+2.D0*B2+3.D0*B3*RS12+2.D0*B4*P1*RSP)
      GGRS = -2.D0*A*A1*Q2-Q0*Q3/(Q1**2+Q1)
      RETURN
      END
      SUBROUTINE CORGGA(RS,ZET,T,UU,VV,WW,H,DVCUP,DVCDN)
C  GGA91 CORRELATION
C  INPUT RS: SEITZ RADIUS
C  INPUT ZET: RELATIVE SPIN POLARIZATION
C  INPUT T: ABS(GRAD D)/(D*2.*KS*G)
C  INPUT UU: (GRAD D)*GRAD(ABS(GRAD D))/(D**2 * (2*KS*G)**3)
C  INPUT VV: (LAPLACIAN D)/(D * (2*KS*G)**2)
C  INPUT WW:  (GRAD D)*(GRAD ZET)/(D * (2*KS*G)**2
C  OUTPUT H: NONLOCAL PART OF CORRELATION ENERGY PER ELECTRON
C  OUTPUT DVCUP,DVCDN:  NONLOCAL PARTS OF CORRELATION POTENTIALS
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/GAS/FK,SK,G,EC,ECRS,ECZET
      common/perloc/iwahl
      DATA XNU,CC0,CX,ALF/15.75592D0,0.004235D0,-0.001667212D0,0.09D0/
      DATA C1,C2,C3,C4/0.002568D0,0.023266D0,7.389D-6,8.723D0/
      DATA C5,C6,A4/0.472D0,7.389D-2,100.D0/
      DATA THRDM,THRD2/-0.333333333333D0,0.666666666667D0/
      BET = XNU*CC0
      DELT = 2.D0*ALF/BET
      G3 = G**3
      G4 = G3*G
      PON = -DELT*EC/(G3*BET)
      B = DELT/(DEXP(PON)-1.D0)
      B2 = B*B
      T2 = T*T
      T4 = T2*T2
      T6 = T4*T2
      RS2 = RS*RS
      RS3 = RS2*RS
      Q4 = 1.D0+B*T2
      Q5 = 1.D0+B*T2+B2*T4
      Q6 = C1+C2*RS+C3*RS2
      Q7 = 1.D0+C4*RS+C5*RS2+C6*RS3
      CC = -CX + Q6/Q7
      R0 = (SK/FK)**2
      R1 = A4*R0*G4
      COEFF = CC-CC0-3.D0*CX/7.D0
      R2 = XNU*COEFF*G3
      R3 = DEXP(-R1*T2)
      H0 = G3*(BET/DELT)*DLOG(1.D0+DELT*Q4*T2/Q5)
      H1 = R3*R2*T2
      H = H0+H1
c     write(6,*)'per2',h,h0,h1
c     write(6,*)'g3,ec,bet,xnu',g3,ec,bet,xnu
c     write(6,*)b,t2,pon
c     write(6,*)r2,r3,delt,q4,t2,q5,delt*q4*t2/q5
C  LOCAL CORRELATION OPTION:
      if(iwahl.eq.3)H=0.0d0
c     H = 0.0D0
C  ENERGY DONE. NOW THE POTENTIAL:
      CCRS = (C2+2.*C3*RS)/Q7 - Q6*(C4+2.*C5*RS+3.*C6*RS2)/Q7**2
      RSTHRD = RS/3.D0
      R4 = RSTHRD*CCRS/COEFF
      if(zet.gt.0.99999999d0) then
      GZ=2.d0**THRDM/3.d0
      else
         if(zet.lt.-0.99999999d0) then
         GZ=-2.d0**THRDM/3.d0
         else
            GZ = ((1.D0+ZET)**THRDM - (1.D0-ZET)**THRDM)/3.D0
         endif
      endif
      FAC = DELT/B+1.D0
      BG = -3.D0*B2*EC*FAC/(BET*G4)
      BEC = B2*FAC/(BET*G3)
      Q8 = Q5*Q5+DELT*Q4*Q5*T2
      Q9 = 1.D0+2.D0*B*T2
      H0B = -BET*G3*B*T6*(2.D0+B*T2)/Q8
      H0RS = -RSTHRD*H0B*BEC*ECRS
      FACT0 = 2.D0*DELT-6.D0*B
      FACT1 = Q5*Q9+Q4*Q9*Q9
      H0BT = 2.D0*BET*G3*T4*((Q4*Q5*FACT0-DELT*FACT1)/Q8)/Q8
      H0RST = RSTHRD*T2*H0BT*BEC*ECRS
      H0Z = 3.D0*GZ*H0/G + H0B*(BG*GZ+BEC*ECZET)
      H0T = 2.*BET*G3*Q9/Q8
      H0ZT = 3.D0*GZ*H0T/G+H0BT*(BG*GZ+BEC*ECZET)
      FACT2 = Q4*Q5+B*T2*(Q4*Q9+Q5)
      FACT3 = 2.D0*B*Q5*Q9+DELT*FACT2
      H0TT = 4.D0*BET*G3*T*(2.D0*B/Q8-(Q9*FACT3/Q8)/Q8)
      H1RS = R3*R2*T2*(-R4+R1*T2/3.D0)
      FACT4 = 2.D0-R1*T2
      H1RST = R3*R2*T2*(2.D0*R4*(1.D0-R1*T2)-THRD2*R1*T2*FACT4)
      H1Z = GZ*R3*R2*T2*(3.D0-4.D0*R1*T2)/G
      H1T = 2.D0*R3*R2*(1.D0-R1*T2)
      H1ZT = 2.D0*GZ*R3*R2*(3.D0-11.D0*R1*T2+4.D0*R1*R1*T4)/G
      H1TT = 4.D0*R3*R2*R1*T*(-2.D0+R1*T2)
      HRS = H0RS+H1RS
      HRST = H0RST+H1RST
      HT = H0T+H1T
      HTT = H0TT+H1TT
      HZ = H0Z+H1Z
      HZT = H0ZT+H1ZT
      COMM = H+HRS+HRST+T2*HT/6.D0+7.D0*T2*T*HTT/6.D0
      PREF = HZ-GZ*T2*HT/G
      FACT5 = GZ*(2.D0*HT+T*HTT)/G
      COMM = COMM-PREF*ZET-UU*HTT-VV*HT-WW*(HZT-FACT5)
      DVCUP = COMM + PREF
      DVCDN = COMM - PREF
C  LOCAL CORRELATION OPTION:
      if(iwahl.eq.3)DVCUP=0.0d0
      if(iwahl.eq.3)DVCDN=0.0d0
c     DVCUP = 0.0D0
c     DVCDN = 0.0D0
      RETURN
      END
      subroutine EABL(x,mmax,mm,ea,r,al)
      implicit real*8 (a-h,o-z)
      dimension x(mm),ea(mm),r(mm)
      mmax2=mmax-2
      mmax1=mmax-1
      ea(1)=(-50.d0*x(1)+96.d0*x(2)-72.d0*x(3)+32.d0*x(4)
     &       -6.d0*x(5))/(24.d0*al*r(1))
      ea(2)=(-6.d0*x(1)-20.d0*x(2)+36.d0*x(3)-12.d0*x(4)
     &       +2.d0*x(5))/(24.d0*al*r(2))
c
      do 22 i=3,mmax2
 22     ea(i)=(2.d0*x(i-2)-16.d0*x(i-1)+16.d0*x(i+1)
     &         -2.d0*x(i+2))/(24.d0*al*r(i))
      i=mmax2
      ea(mmax1)=(-2.d0*x(i-2)+12.d0*x(i-1)-36.d0*x(i)+20.d0*x(i+1)
     &         +6.d0*x(i+2))/(24.d0*al*r(i+1))
      ea(mmax)=(6.d0*x(i-2)-32.d0*x(i-1)+72.d0*x(i)-96.d0*x(i+1)
     &        +50.d0*x(i+2))/(24.d0*al*r(i+2))
      return
      end
      subroutine ZABL(x,mmax,mm,za,r,al)
      implicit real*8 (a-h,o-z)
      dimension x(mm),za(mm),r(mm)
      mmax2=mmax-2
      mmax1=mmax-1
      za(1)=(35.d0*x(1)-104.d0*x(2)+114.d0*x(3)-56.d0*x(4)
     &       +11.d0*x(5))/(12.d0*al*al*r(1)*r(1))
      za(2)=(11.d0*x(1)-20.d0*x(2)+6.d0*x(3)+4.d0*x(4)
     &       -1.d0*x(5))/(12.d0*al*al*r(2)*r(2))
c
      do 22 i=3,mmax2
 22     za(i)=(-1.d0*x(i-2)+16.d0*x(i-1)-30.d0*x(i)+16.d0*x(i+1)
     &         -1.d0*x(i+2))/(12.d0*al*al*r(i)*r(i))
      i=mmax2
      za(mmax1)=(-1.d0*x(i-2)+4.d0*x(i-1)+ 6.d0*x(i)-20.d0*x(i+1)
     &         +11.d0*x(i+2))/(12.d0*al*al*r(i+1)*r(i+1))
      za(mmax)=(11.d0*x(i-2)-56.d0*x(i-1)+114.d0*x(i)-104.d0*x(i+1)
     &        +35.d0*x(i+2))/(12.d0*al*al*r(i+2)*r(i+2))
      return
      end
      subroutine EZABL(x,mmax,mm,ea,za,r,al)
      implicit real*8 (a-h,o-z)
      dimension x(mm),ea(mm),r(mm),za(mm)
      mmax2=mmax-2
      mmax1=mmax-1
      ea(1)=(-50.d0*x(1)+96.d0*x(2)-72.d0*x(3)+32.d0*x(4)
     &       -6.d0*x(5))/(24.d0*al*r(1))
      ea(2)=(-6.d0*x(1)-20.d0*x(2)+36.d0*x(3)-12.d0*x(4)
     &       +2.d0*x(5))/(24.d0*al*r(2))
c
      do 22 i=3,mmax2
 22     ea(i)=(2.d0*x(i-2)-16.d0*x(i-1)+16.d0*x(i+1)
     &         -2.d0*x(i+2))/(24.d0*al*r(i))
      i=mmax2
      ea(mmax1)=(-2.d0*x(i-2)+12.d0*x(i-1)-36.d0*x(i)+20.d0*x(i+1)
     &         +6.d0*x(i+2))/(24.d0*al*r(i+1))
      ea(mmax)=(6.d0*x(i-2)-32.d0*x(i-1)+72.d0*x(i)-96.d0*x(i+1)
     &        +50.d0*x(i+2))/(24.d0*al*r(i+2))
C
      za(1)=(35.d0*x(1)-104.d0*x(2)+114.d0*x(3)-56.d0*x(4)
     &       +11.d0*x(5))/(12.d0*al*al*r(1)*r(1))
      za(2)=(11.d0*x(1)-20.d0*x(2)+6.d0*x(3)+4.d0*x(4)
     &       -1.d0*x(5))/(12.d0*al*al*r(2)*r(2))
c
      do 23 i=3,mmax2
 23     za(i)=(-1.d0*x(i-2)+16.d0*x(i-1)-30.d0*x(i)+16.d0*x(i+1)
     &         -1.d0*x(i+2))/(12.d0*al*al*r(i)*r(i))
      i=mmax2
      za(mmax1)=(-1.d0*x(i-2)+4.d0*x(i-1)+ 6.d0*x(i)-20.d0*x(i+1)
     &         +11.d0*x(i+2))/(12.d0*al*al*r(i+1)*r(i+1))
      za(mmax)=(11.d0*x(i-2)-56.d0*x(i-1)+114.d0*x(i)-104.d0*x(i+1)
     &        +35.d0*x(i+2))/(12.d0*al*al*r(i+2)*r(i+2))
      do 24 i=1,mmax
24    za(i)=za(i)-ea(i)/r(i)
      za(mmax-1)=0.0d0
      za(mmax-2)=0.0d0
      za(mmax-3)=0.0d0
      ea(mmax-1)=0.0d0
      ea(mmax-2)=0.0d0
      return
      end
