        program rocket
        implicit none

        character*80 version
        parameter (version='Program Rocket, version 3.0, 12/09/87')
*-------------------------------------------------------------------------
*       Program Overview
*-------------------------------------------------------------------------
*       The purpose of this program is to calculate the altitude and
*       velocity of a rocket at a given time t, taking into account
*       air resistance, changing air density, changing gravity,
*       the thrust of the rocket, and the changing mass of the rocket.
*       It would be impossible to determine a closed solution for
*       this problem, so using numerical methods is the only possible
*       way to find the solution. In this program, the Runge-Kutta
*       method is used for integration (4th order). 
*-------------------------------------------------------------------------
*       Variables
*-------------------------------------------------------------------------
*       Earthradius (m) Radius of Earth
*       Initgravity (m/s/s) Gravity at surface
*       Initdensity (kg/m/m/m) Air density at surface
*       R (Joules/kg/K) Gas constant for air
*       Temp (K) Mean temperature of atmosphere
*       Endtime (s) When to stop plotting
*       Rocketmass (kg) Mass of rocket without fuel
*       Fuelmass (kg) Mass of fuel
*       Impulse (s) Impulse of engine
*                   300 - Kerosene/Oxygen
*                   360 - Hydrogen/Oxygen
*                   490 - Hydrogen/Fluorine
*       Burnrate (kg/s) Amount of fuel burned per second
*       Endthrust (s) Total time of thrust
*       Velocity (m/s) Velocity of material ejected from rocket nozzle
*       Initthrust (N) Thrust of engine (constant)
*       Cd (dimensionless) Coefficient of Drag
*       SurfaceArea (m**2) Frontal surface area
*       Initheight (m) Initial altitude of rocket
*       Initvelocity (m/s) Initial velocity of rocket
*       time (s) Current time
*       dt (s) Step size
*       Vold,Vnew (m/s) Velocity
*       Hold, Hnew (m) Altitude
*       k1-4, L1-4 Used in the Runge-Kutta routine
*       storage provides storage for any variables you might wish to
*               save for later output.
*       snum is a counter to be used with storage.
*       function g computes the overall acceleration of the rocket at
*               a given time, altitude, and velocity.
*
*-------------------------------------------------------------------------

        real Earthradius,Initgravity,Initdensity,R,Temp,Endtime
        real Rocketmass,Fuelmass,Impulse,Burnrate,Endthrust
        real Velocity,Initthrust,Cd,SurfaceArea,Initheight,Initvelocity
        real time,dt,Vold,Vnew,Hold,Hnew,k1,k2,k3,k4,L1,L2,L3,L4
        real storage (0:1,1000),Vterminal
        integer i,snum
        real g,Thrust,Mass,Gravity,Drag,Density


        common /c1/Rocketmass,Fuelmass,Burnrate,Initthrust,Endthrust
        common /c2/Cd,SurfaceArea
        common /c3/Earthradius,Initdensity,Initgravity,R,Temp

*-------------------------------------------------------------------------
*       Initialization
*-------------------------------------------------------------------------

        Earthradius=6.356766e6
        Initgravity=9.80
        Initdensity=1.2250
        R=287.0
        Temp=260.0

        Rocketmass=100.0
        Fuelmass=30.0
        Impulse=300.0
        Burnrate=1.0
        Endthrust=Fuelmass/Burnrate
        Velocity=Initgravity*Impulse
        Initthrust=Burnrate*Velocity
        Cd=0.32
        SurfaceArea=0.1
        Initheight=0.0
        Initvelocity=0.0

        dt=0.2
        Endtime=1.e9
*       (by setting Endtime=1.e9, we can watch the rocket crash)

        open (unit=11,file='graph1.dat',status='new')
        write (11,*) 'LINE NOMARKER XFORMAT I3 YFORMAT I5'
        write (11,*) 'TITLE "ALTITUDE VS. TIME"'
        write (11,*) 'XLABEL "TIME (SEC)"'
        write (11,*) 'YLABEL "ALTITUDE (METERS)"'
        write (11,*) 'NEWPEN 2'

        open (unit=12,file='graph2.dat',status='new')
        write (12,*) 'LINE NOMARKER XFORMAT I3 YFORMAT I5'
        write (12,*) 'TITLE "VELOCITY VS. TIME"'
        write (12,*) 'XLABEL "TIME (SEC)"'
        write (12,*) 'YLABEL "VELOCITY (METERS/SEC)"'
        write (12,*) 'NEWPEN 4'
        
        open (unit=13,file='graph3.dat',status='new')
        write (13,*) 'LINE NOMARKER XFORMAT I3 YFORMAT I5'
        write (13,*) 'TITLE "ACCELERATION VS. TIME"'
        write (13,*) 'XLABEL "TIME (SEC)"'
        write (13,*) 'YLABEL "ACCELERATION (METERS/SEC/SEC)"'
        write (13,*) 'NEWPEN 3'

*-------------------------------------------------------------------------
*       Main Program
*-------------------------------------------------------------------------

        print *,version
        
        Hold=Initheight
        Vold=Initvelocity
        time=0
        snum=1

100     continue
        write (11,*) time,Hold
        write (12,*) time,Vold
        write (13,*) time,g(time,Hold,Vold)

        k1=Vold
        L1=g(time,Hold,k1)
        k2=Vold+dt*L1*(0.5)
        L2=g(time+dt*(0.5),Hold+dt*k1*(0.5),k2)
        k3=Vold+dt*L2*(0.5)
        L3=g(time+dt*(0.5),Hold+dt*k2*(0.5),k3)
        k4=Vold+dt*L3
        L4=g(time+dt,Hold+dt*k3,k4)

        Hnew=Hold+dt/6.*(k1+2*k2+2*k3+k4)
        Vnew=Vold+dt/6.*(L1+2*L2+2*L3+L4)

        if (Hnew .lt. 0) then
            print *,'The rocket has crashed.'
            go to 200
        endif

        Hold=Hnew
        Vold=Vnew
        time=time+dt

        if (time .lt. Endtime) then
            go to 100
        endif

        print *,'Long enough.'

200     continue

        write (11,*) 'NEWPEN 1'
        write (12,*) 'NEWPEN 1'
        write (13,*) 'NEWPEN 1'

        end

*-------------------------------------------------------------------------
*       Functions
*-------------------------------------------------------------------------

        function g(t,h,v)
        implicit none

*       This function computes the total acceleration of the 
*       rocket at a given time, height, and velocity.

        real t,h,v,g,Thrust,Drag,Mass,Gravity

        g=(Thrust(t)+Drag(v,h))/Mass(t)-Gravity(h)

        return
        end




*-------------------------------------------------------------------------

        function Thrust(t)
        implicit none

*       This function computes the thrust at a given time. Note that
*       a more complex thrust curve could be put here rather than just
*       a constant value.

        common /c1/Rocketmass,Fuelmass,Burnrate,Initthrust,Endthrust
        real Rocketmass,Fuelmass,Burnrate,Initthrust,Endthrust
        real Thrust,t

        Thrust=Initthrust
        if (t .ge. Endthrust) then
            Thrust=0
        endif

        return
        end
        
*-------------------------------------------------------------------------

        function Mass(t)
        implicit none

*       This function computes the mass of the rocket at a given 
*       time t.

        common /c1/Rocketmass,Fuelmass,Burnrate,Initthrust,Endthrust
        real Rocketmass,Fuelmass,Burnrate,Initthrust,Endthrust
        real Mass,t
        
        Mass=Rocketmass
        if (t .lt. Endthrust) then
            Mass=Rocketmass+Fuelmass-Burnrate*t
        endif

        return
        end

*-------------------------------------------------------------------------

        function Gravity(h)
        implicit none

*       This function computes the acceleration of gravity at an
*       altitude h, given that the gravity at the surface is
*       Initgravity.

        common /c3/Earthradius,Initdensity,Initgravity,R,Temp
        real Earthradius,Initdensity,Initgravity,R,Temp
        real Gravity,h

        Gravity=Initgravity*(Earthradius/(Earthradius+h))**2

        return
        end



*-------------------------------------------------------------------------

        function Drag(v,h)
        implicit none

*       Function Drag computes the air drag on the rocket for a certain
*       velocity and altitude. Note that drag acts in the direction 
*       opposite of the given velocity v. 

        common /c2/Cd,SurfaceArea
        real Cd,SurfaceArea
        real Drag,v,h
        real Density

        Drag=Sign(Cd*(0.5)*Density(h)*v*v*SurfaceArea,-v)

        return
        end

*-------------------------------------------------------------------------

        function Density(h)
        implicit none
        
*       Function Density finds the density at a given altitude h.
*       Note that a more complex formula could be used, but this
*       is the same model that NASA used in reentry studies in the
*       50's and 60's, so it's good enough for me.

        common /c3/Earthradius,Initdensity,Initgravity,R,Temp
        real Earthradius,Initdensity,Initgravity,R,Temp
        real Density,h

        Density=Initdensity*exp(-Initgravity*h/R/Temp)

        return
        end

