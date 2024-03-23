program Ex3
implicit none
integer, parameter :: natoms=125
real*8 :: r(natoms,3), rold(natoms,3), v(natoms,3),accel(natoms,3), rho , L, epot=0.0, rc, pot, ecin=0.0, dt, Temp,Tinst,kb,velmod
real*8 :: mass, epsilon, momentum(3), ppot, Pressure, rhoOld, g(100), ri(100), delg, nid,pi,vb
integer :: i,j,is, js, M, jj, iii, nhis,iind,s
character(len=1) :: aux
! Set parameters
rho = 0.7
L=(real(natoms,8)/rho)**(1./3.)
rc = L/2.0
mass = 1.
epsilon = 1.0

!open(5,file="thermo005.dat")
!open(4,file="thermo01.dat")
!open(3,file="thermo02.dat")
!open(2,file="thermo04.dat")
open(6,file="thermo06.dat")
!open(1,file="thermo08.dat")
!open(11,file="gdr005.dat")
!open(10,file="gdr01.dat")
!open(9,file="gdr02.dat")
!open(8,file="gdr04.dat")
open(12,file="gdr06.dat")
!open(7,file="gdr08.dat")
!open(17, file="Trajectories005.xyz")
!open(14, file="Trajectories04.xyz")
!open(13, file="Trajectories08.xyz")

! Read the disordered initial positions
open(20,file='InitConf.xyz') ! This is for rho' = 0.7
read(20,*) 
read(20,*)
do i=1,natoms
   read(20,*) aux, (r(i,j), j=1,3)
end do
close(20)

!rho = 0.8
rho = 0.6
rhoOld = 0.7
do jj=6,6
	
	L=(real(natoms,8)/rho)**(1./3.)
	rc = L/2.0
	
! 	inici gdr
	nhis = 100
	delg =L/(2*nhis) !bin size
	do iii = 1,nhis !nhis = total number of bins
	g(iii) = 0.0
	ri(iii)=delg*DFLOAT(iii -1) !distance r
	end do
	
	

	!velmod = 17.111!dsqrt(epsilon*Temp/mass)
	do i=1,natoms
	  do j=1,3
	    r(i,j) = (rhoOld/rho)**(1./3.) * r(i,j)
	    v(i,j) = 0.0 
	  end do
	end do


	dt=0.001
	Temp=1.2

	!open(2,file="thermodynamicsVerlet.dat")
	do i=1,500000
	  call velocityverlet(natoms,r,v,accel,L,rc,dt,ecin,epot,Temp,momentum,ppot,g,delg)  
	  Tinst = 2.*ecin/(3.*real(natoms,8)-3.)
	  Pressure = rho*Tinst + (1.d0/(3.0d0*(L**3.0d0)))*ppot
	  if (mod(i,10).eq.0) then
	    write(jj,*) i*dt, ecin, epot, ecin+epot, Pressure
	  end if
	  
	  ! Save part of the trajectory
	  if ((i>=490000) .and. (i<=500000) .and. (mod(i,100).eq.0) .and. ((jj.eq.1) .or. (jj.eq.2) .or. (jj.eq.5)) ) then
	    write(12+jj,*) natoms
            write(12+jj,*) i-1,'TIMESTEP:',jj 
            do iind=1,natoms
              write(12+jj,*) 'A',(r(iind,s),s=1,3)
            end do
	  end if 
	  
	end do
	close(jj)
	!close(2)
	!open(9,file="velocitiesVerlet.dat")
	!do is=1,natoms
	!  write(9,*)  (v(is,js), js=1,3)
	!end do
	! close(9)
	
!   Normalize gdr
	pi=3.14159265359
	do i=1,nhis
	 ri(i)=delg*DFLOAT(i -1) !distance r
	 vb = ((i+1)**3-i**3)*delg**3 !volume between bin i+1 and i
	 nid = (4.0/3.0)*pi*vb*rho !number of ideal gas part . in vb
	 g(i) =g(i)/(500000*natoms*nid) !normalize g(r)
	write(jj+6,*)  ri(i), g(i)
	enddo
	close(jj+6)
	rhoOld = rho
	rho = rho/2.0

end do
 
stop
end program Ex3


!------------------------------------
!---------SUBROUTINES----------------
!------------------------------------

!---------SC----------------

subroutine SCbuild(N,r,L, M)
implicit none
integer :: N
real*8, dimension(N,3) :: r  !, v(N,3)
real*8 ::  L, a
integer :: i, j, k, ind=1, M
M = N**(1./3.)
a=L/M
do i=0,M-1
   do j=0,M-1
      do k=0,M-1
        r(ind, 1) = real(i,8)*a
        r(ind, 2) = real(j,8)*a
        r(ind, 3) = real(k,8)*a
        ind=ind+1
      end do
   end do
end do
end subroutine SCbuild

!----------LJ--------------------

subroutine lj(N,is,js,r, accel, boxlength, rc, pot,press,g,delg)
implicit none
integer:: N
integer :: is, js, l, aux, ig
real*8, dimension(N,3) :: r, accel
real*8 :: rr, rr2, pot, rijl, ynvrr1, ynvrr2, ynvrr6, ynvrr12
real*8 :: rij(3), rc, boxlength, force, forcedist, press, g(100), delg

rr2 = 0.d0
pot = 0.d0
press = 0.d0
do l = 1,3
 rijl = r(js,l) - r(is,l)
 rij(l) = rijl - boxlength*dnint(rijl/boxlength) !with pbc's
 rr2 = rr2 + rij(l)*rij(l)
end do

rr = dsqrt(rr2)
if (rr.lt.rc) then
ynvrr1 = 1.d0/rr
ynvrr2 = 1.d0/rr2
ynvrr6 = ynvrr2*ynvrr2*ynvrr2
ynvrr12 = ynvrr6*ynvrr6
forcedist = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr2
force = 24.d0*(2.d0*ynvrr12-ynvrr6)*ynvrr1
pot = 4.d0*(ynvrr12-ynvrr6) - 4.d0*(1./(rc**12)-1./(rc**6))
press = force*rr
  do l = 1,3
    accel(is,l) = accel(is,l) - forcedist*rij(l)
    accel(js,l) = accel(js,l) + forcedist*rij(l)
 end do
end if

! 	Part de gdr
if (rr.lt.boxlength/2.) then !only within half the box length
 ig = int(rr/delg)
 g(ig) = g(ig) + 2 !contribution for particle i and j
 endif

return
end


!--------------FORCES----------------

subroutine forces(natoms,r,boxlength,accel,rc,epot,ppot,g,delg)
implicit none
integer :: natoms, is, js, l
real*8 :: r(natoms,3),accel(natoms,3), epot,rc,boxlength, pot,press, ppot, g(100), delg

do is = 1,natoms
 do l = 1,3
    accel(is,l) = 0.d0 !sets accelerations to 0
 end do 
end do 
epot = 0.d0 
ppot = 0.d0
do is = 1,natoms
 do js = is+1,natoms
    call lj(natoms,is,js,r, accel, boxlength, rc, pot,press,g,delg)
    epot = epot + pot
    ppot = ppot + press
 end do
end do

return
end


!-----------Velocity Verlet-------------

subroutine velocityverlet(natoms,r,vel,accel,boxlength, rc,dt,ecin,epot,Temp,momentum,ppot,g,delg)
implicit none
integer :: natoms, is, js, l
real*8 :: r(natoms,3), vel(natoms,3), rnew_old, ecin, accel(natoms,3), accelold(natoms,3)
real*8 :: boxlength, epot, rc, dt, v2, Temp, momentum(3), ppot, g(100), delg
ecin = 0.0
epot = 0.0 

call forces(natoms,r,boxlength,accel,rc,epot,ppot,g,delg)
do is = 1,natoms

 do l = 1,3
    r(is,l) = r(is,l) + vel(is,l)*dt + 0.5*accel(is,l)*dt*dt
    vel(is,l) = vel(is,l) + 0.5*(accel(is,l))*dt
    if (r(is,l) .lt. 0.0d0) then 
       r(is,l) = r(is,l) + boxlength
    else if (r(is,l) .gt. boxlength) then
       r(is,l) = r(is,l) - boxlength
    end if
  end do
end do
epot = 0.0 
ecin = 0.0

call forces(natoms,r,boxlength,accel,rc,epot,ppot,g,delg) 
momentum(1) = 0.0
momentum(2) = 0.0
momentum(3) = 0.0
do is = 1,natoms
 v2 = 0.d0
  do l=1,3
    vel(is,l) = vel(is,l) + 0.5*(accel(is,l))*dt
  end do
  call therm_Andersen(natoms,is,vel, 0.1d0, Temp)
  do l=1,3
    v2 = v2 + vel(is,l)*vel(is,l)
    momentum(l) = momentum(l) + vel(is,l)
  end do
 ecin = ecin + 0.5d0*v2
end do


end subroutine velocityverlet

!-----------Euler---------------

subroutine Euler(natoms,r,vel,accel,boxlength, rc,dt,ecin,epot,Temp,momentum,ppot,g,delg)
implicit none
integer :: natoms, is, js, l
real*8 :: r(natoms,3), rold(natoms,3), vel(natoms,3), rnew_old, ecin, accel(natoms,3), accelold(natoms,3)
real*8 :: boxlength, epot, rc, dt, v2, ppot, momentum(3), Temp, g(100), delg
ecin = 0.0
epot = 0.0 
momentum(1) = 0.0
momentum(2) = 0.0
momentum(3) = 0.0
call forces(natoms,r,boxlength,accel,rc,epot,ppot,g,delg)
do is = 1,natoms
 v2 = 0.d0
 do l = 1,3
    r(is,l) = r(is,l) + vel(is,l)*dt + 0.5*accel(is,l)*dt*dt
    vel(is,l) = vel(is,l) + 0.5*(accel(is,l))*dt
    v2 = v2 + vel(is,l)*vel(is,l)
    momentum(l) = momentum(l) + vel(is,l)
 end do
 ecin = ecin + 0.5d0*v2
 
end do

end subroutine Euler


!-----------Andersen thermostat----

subroutine therm_Andersen(natoms,is,vel, nu, Temp)
implicit none
integer :: is, js,l, natoms
real*8 :: vel(natoms,3), nu, sigma, Temp, PI, randnum1, randnum2, randvel
PI = 4.*atan(1.0d0)
sigma = sqrt(Temp)
call random_number(randnum1)
if (randnum1<nu) then
do l=1,3
  call random_number(randnum1)
  call random_number(randnum2)
  vel(is,l) = sigma*dsqrt(-2d0*(dlog(1d0-randnum1)))*dcos(2d0*PI*randnum2)
end do
end if
end subroutine therm_Andersen
