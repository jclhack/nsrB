* nsr_Bdd_20250226: version for UIUC cluster
* nsr_Bdd_20210301: geometry to better match proposed sample       
* nsr_Bdd_20210222: B-field components, classical dipole-dipole, source above sample in z. 
* ptb_Beff_20200701: effective B-field components, PTB, geometry from 6/20 emails
* ptb_Beff_20200630: effective B-field components, PTB, imp. samp., one source slab with exclusion zone
* ariadne_Beff_20200526: effective B-field components, imp. samp., one source slab with exclusion zone
* ariadne_Beff_20200525: effective B-field components, no imp. samp., one source slab
* ariadne_signal_20200523: basic ariadne geometry; source simple z-translation, geometry tracking, importance
* sampling flag  
* yukawa_simple_20200522
* based on Yukawa_simple_20200423: added variables for spin vectors in detector and source; gap in
* terms of closest approach; no "Nflag"
* based on yukawa_simple_20200422: essential elements only 
*	  
      program nsr_Bdd_20250226
*
      implicit none
*
* variables
*
      real*4 Results(4,100,100) ! results stored in Results(i,j,k)
				! k runs over games, j runs over phase points
				! Results(1,j,k)=Bx: computed effective field
				! Results(2,j,k)=By: computed effective field
				! Results(3,j,k)=Bz: computed effective field
				! Results(4,j,k)=points: number of inside points
*
      real*4 Stats(7,100) ! statistics stored in Stats(i,j)
	  		  ! j runs over phase points
	  		  ! Stats(1,j)=omega: source mass position
			  ! Stats(2,j)=Bxmean: mean effective x-field at this phase
			  ! Stats(3,j)=Bymean: mean effective y-field at this phase
			  ! Stats(4,j)=Bzmean: mean effective z-field at this phase
			  ! Stats(5,j)=Bxdev: stan dev of mean, x-field
			  ! Stats(6,j)=Bydev: stan dev of mean, y-field
			  ! Stats(7,j)=Bzdev: stan dev of mean, z-field
*	  						 
      real*4 range ! range of force
*
      integer*4 throws ! number of trial points per game
      integer*4 games  ! number of games per phase point
      integer*4 phpnts ! number of phase points
*
      integer*4 points ! number of points inside in a game
*
      real*4 Bx,By,Bz             ! effective field components
      real*4 Bxampl,Byampl,Bzampl ! Fourier amplitudes of components
*
      real*4 bounds1(6) ! limits for point in mass 1,
			! cartesian coordinates, order is
			! x1 upper, x1 lower, y1 upper, y1 lower
			! z1 upper, z1 lower
*
      real*4 bounds2(6) ! limits for point in mass 2
			! spherical coords about z axis, order is
			! R upper, R lower, theta upper, theta lower,
			! phi upper, phi lower
*
      real*4 shape1(36)	! shape of mass 1 (detector) in its cartesian
	  	        ! includes detector mode shape
      real*4 shape2(36)	! shape of mass 2 (source) in its cartesians
	  		! includes source mode shape
*
      real*4 xmax,ymax,zmax ! maximum separations of points in x, y, z for bounding box 
      real*4 omega	    ! source mass phase
      real*4 asinom         ! amp*sin(omega)
      real*4 D(3)	    ! source position relative to detector
	  	            ! in detector cartesians
*
      integer*4 i,j,k ! index
*
      real*4 volume    ! bounds box volume
      real*4 point1(3) ! point in box 1; x, y, z
      real*4 point2(3) ! point in box 2; scaled R, theta, phi
      real*4 modez     ! z mode shape function
*
      real*4 snth, csth   ! sin(theta), cos(theta)
      real*4 snphi, csphi ! sin(phi), cos(phi)
      real*4 R		  ! radial coord, mass 2
*
      logical*4 inside ! .true. if both points are inside masses
*	  
      real*4 dBx,dBy,dBz ! differential field components in integrand
*
* spin-dependent variables
*  
      real*4 theta1, phi1   ! polar, azimuthal angle of sample (e) spin in source coordinates
      real*4 snth1, csth1   ! sin(theta1), cos(theta1)
      real*4 snphi1, csphi1 ! sin(phi1), cos(phi1)
*
      real*4 theta2, phi2   ! polar, azimuthal angle of source (n) spin
      real*4 snth2, csth2   ! sin(theta2), cos(theta2)
      real*4 snphi2, csphi2 ! sin(phi2), cos(phi2)
*
      integer*4 sflag ! importance sampling flag
*
      real*4 Pi ! c/d of a circle
*
* initialize
*
      Pi=3.1415926535 ! c/d of a circle
*
      sflag=0 ! set to 1 for importance sampling (not used initially) 
*	
      throws=1000000
      games=80
      phpnts=5
*
      range=10000.0 ! not used in classical case
*
      bounds2(3)=Pi/2.0 ! theta upper
      bounds2(4)=0.0    ! theta lower
      bounds2(5)=2.0*Pi ! phi upper
      bounds2(6)=0.0    ! phi lower
*
      shape1(1)=16882.6 ! detector x size; area matches disk of radius 19 mm
      shape1(2)=16882.6	! detector y size; area matches disk of radius 19 mm
      shape1(3)=10000.0	! detector z size
*
      bounds1(1)= 1.0             ! x upper	
      bounds1(2)=-(shape1(1)+1)   ! x lower
      bounds1(3)= (shape1(2)/2+1) ! y upper
      bounds1(4)=-(shape1(2)/2+1) ! y lower
      bounds1(5)= (shape1(3)/2+1) ! z upper
      bounds1(6)=-(shape1(3)/2+1) ! z lower
*
      shape1(4)=0 ! detector surface curvature polynomial:
      shape1(5)=0 ! s1(4) + s1(5)x + s1(6)y
      shape1(6)=0 ! + s1(7)xy + s1(8)x**2 + s1(9)y**2
      shape1(7)=0 ! all curvature and mode polynomials
      shape1(8)=0 ! have this form and all are
      shape1(9)=0 ! in microns of displacement
*	  
      shape1(10)=1.0 ! detector z mode shape polynomial
      shape1(11)=0 
      shape1(12)=0  
      shape1(13)=0 
      shape1(14)=0 
      shape1(15)=0 
*	  
      shape1(16)=0 ! detector x virtual displacement
      shape1(17)=0 ! detector y virtual displacement
*	  
      shape2(1)=0.1*shape1(1)	 ! source x size
      shape2(2)=0.1*shape1(2)	 ! source y size
      shape2(3)=0.1*shape1(3)	 ! source z size
*	  
      shape2(4)=0 ! source surface curvature polynomial
      shape2(5)=0 
      shape2(6)=0 
      shape2(7)=0 
      shape2(8)=0 
      shape2(9)=0 
*	  
      shape2(10)=1.0 ! source z mode shape polynomial
      shape2(11)=0   ! 
      shape2(12)=0   ! conventional normalization is 
      shape2(13)=0   ! one micron at cartesian origin
      shape2(14)=0   ! 
      shape2(15)=0   ! 
*
      shape2(18)=25000.0   ! source exclusion zone x, not used
      shape2(19)=25000.0 ! source exclusion zone y, not used
      shape2(20)=25000.0  ! source exclusion zone z, not used
*	  
      D(1)=-shape1(1)/2.0-shape2(1)/2.0  ! initial source displacement x
      D(2)=0.0	                         ! y
      D(3)=shape1(3)/2.0+shape2(3)/2.0   ! z
*
      xmax=shape1(1)/2.0+shape2(1)/2.0 ! max distance between any 2 points in x
      ymax=shape1(2)/2.0+shape2(2)/2.0 ! max distance between any 2 points in y
      zmax=6.0*shape1(3)+shape2(3)         ! max distance between any 2 points in z
*
      bounds2(1)=sqrt(xmax**2+ymax**2+zmax**2) ! R upper, based on max between any points
      bounds2(2)=1.0    		       ! R lower, closest approach between source & detector
*
* replace R limits with q (scaled R) limits (importance sampling)
*
      if (sflag .gt. 0) then 
      	bounds2(1)=exp(-bounds2(1)/range)
      	bounds2(2)=exp(-bounds2(2)/range)
      end if
* 
* output files 
*
      open(9,file='W_out.txt',status='unknown')
      open(10,file='W_phase.txt',status='unknown')
      open(11,file='inside_points.txt',status='unknown')
*                
* compute volume
*
      volume=(bounds1(1)-bounds1(2))*(bounds1(3)-bounds1(4))*
     >(bounds1(5)-bounds1(6))*(bounds2(1)-bounds2(2))*
     >(bounds2(3)-bounds2(4))*(bounds2(5)-bounds2(6))
*
* loop over source mass phase values
*
      do 90 j=1,phpnts
      	print *,'start phase loop',j
      	omega=(5.0*shape1(3))*(1.0+1.0/phpnts)/phpnts
      	asinom=(j-1)*omega
      	Stats(1,j)=D(3)+asinom
*
* loop over games
*
      	do 80 k=1,games
		points=0
      		Bx=0.0
      		By=0.0
      		Bz=0.0
*
* loop over throws
*
      		do 50 i=1,throws	
      			call dice(point1,point2,bounds1,bounds2) ! get a point
*      
      			snth=sin(point2(2))
      			csth=cos(point2(2))
      			snphi=sin(point2(3))
      			csphi=cos(point2(3))
                        if (sflag .gt. 0) then
      				R=-range*log(point2(1))
                        else
                        	R=point2(1)
                        end if
*
      			call Plates(shape1,shape2,D,point1,point2,   ! is the point inside?
     >asinom,R,snth,csth,snphi,csphi,inside)
*
      			if (inside) then
      				points=points+1                      ! bump points count
*
                                if (sflag .gt. 0) then
					dBx=0.0     ! compute field components 
					dBy=0.0
					dBz=0.0
				else
      				write(11,45) point1(1),point1(2),
     >point1(3),point2(1),point2(2),point2(3),j
45    format(e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,i4)
					dBx=(3.0/R)
     >*snth**2*csth*csphi				! compute field components 
					dBy=(3.0/R)
     >*snth**2*csth*snphi	  		
					dBz=(1.0/R)
     >*snth*(3.0*csth**2-1.0)	  		
				end if
      				Bx=Bx+dBx
      				By=By+dBy
      				Bz=Bz+dBz
      			end if	  
50    		continue
*
* normalize and record results
*
      		Bx=Bx*volume/throws
      		By=By*volume/throws
      		Bz=Bz*volume/throws
     		Results(1,j,k)=Bx
     		Results(2,j,k)=By
     		Results(3,j,k)=Bz
		Results(4,j,k)=points
*
* report progress
*
*     		write(9,65) j,k
*65    		format (i4,i4)
*	  
80    	continue ! close loop over games		
*
90    continue ! close loop over phase points
*
* compute means
*
      do 85 j=1,phpnts
      	Stats(2,j)=0
      	Stats(3,j)=0
      	Stats(4,j)=0
      	do 82 k=1,games
      		Stats(2,j)=Stats(2,j)+Results(1,j,k)	
      		Stats(3,j)=Stats(3,j)+Results(2,j,k)	
      		Stats(4,j)=Stats(4,j)+Results(3,j,k)	
82    	continue
      	Stats(2,j)=Stats(2,j)/games	
      	Stats(3,j)=Stats(3,j)/games	
      	Stats(4,j)=Stats(4,j)/games	
85    continue
*
* standard deviation of mean
*
      do 95 j=1,phpnts
      	Stats(5,j)=0
      	Stats(6,j)=0
      	Stats(7,j)=0
      	do 92 k=1,games
      		Stats(5,j)=Stats(5,j)+(Results(1,j,k)-Stats(2,j))**2  
      		Stats(6,j)=Stats(6,j)+(Results(2,j,k)-Stats(3,j))**2  
      		Stats(7,j)=Stats(7,j)+(Results(3,j,k)-Stats(4,j))**2  
92    	continue
      	Stats(5,j)=(Stats(5,j)**0.5)/games	
      	Stats(6,j)=(Stats(6,j)**0.5)/games	
      	Stats(7,j)=(Stats(7,j)**0.5)/games	
95    continue
*
* write Results example
*
*      write(9,125)
*125   format(/,'phase  game   Bx   points',/)
*
*      do 135 j=1,phpnts
*      	do 136 k=1,games
*      	write(9,100) j,k,Results(1,j,k),Results(2,j,k)
*100    format(i5,i5,e12.5,e12.5)	  
*136 	continue		
*135   continue
*
* write Stats
*
*      write(9,225)
      write(10,225)
225   format('  omega        Bxmean       Bymean       Bzmean
     >Bxdev       Bydev       Bzdev')
      do 235 j=1,phpnts
*     	write(9,200) Stats(1,j),Stats(2,j),Stats(3,j)
      	write(10,200) Stats(1,j),Stats(2,j),Stats(3,j),Stats(4,j),
     >Stats(5,j),Stats(6,j),Stats(7,j)
200   	format(e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5)	  		
235   continue
*
* compute amplitude of effective field and write
*
      Bxampl=0.0
      Byampl=0.0
      Bzampl=0.0
      do 245 j=1,phpnts
      	omega=6.2832*(j-1)/phpnts
      	Bxampl=Bxampl+sin(omega)*Stats(2,j) ! Fourier amplitude		
      	Byampl=Byampl+sin(omega)*Stats(3,j)		
      	Bzampl=Bzampl+sin(omega)*Stats(4,j)		
245   continue
      Bxampl=Bxampl*2.0/phpnts ! scale Fourier amplitude
      Byampl=Byampl*2.0/phpnts
      Bzampl=Bzampl*2.0/phpnts
      write(9,255)
255   format('Bxampl       Byampl       Bzampl')
      write(9,260) Bxampl,Byampl,Bzampl
260   format(e13.5,e13.5,e13.5)
*
      print *, "done"
*
      end
*
* functions and subroutines
*
*** modez ***
*
*  mode shape function
*
      function modez(point,shape)
      real*4 modez,point(3),shape(36)
      modez = shape(10) + shape(11)*point(1) + shape(12)*point(2)
     > + shape(13)*point(1)*point(2)
     > + shape(14)*point(1)**2 + shape(15)*point(2)**2
      end
*	  
*** curve ***
*
*  curvature function
*
      function curve(point,shape)
      real*4 curve,point(3),shape(36)
      curve = shape(4) + shape(5)*point(1) + shape(6)*point(2)
     > + shape(7)*point(1)*point(2)
     > + shape(8)*point(1)**2 + shape(9)*point(2)**2
      end
*
*** Plates.f ***
*
*  returns inside=.true. iff point1 is inside mass 1 according to shape1
*  and point2 is inside mass 2 according to shape2.  This version is for
*  our planar torsional oscillator.
*
      subroutine Plates(shape1,shape2,D,point1,point2,asinom,
     > R,snth,csth,snphi,csphi,inside)
      implicit none
*
* variables passed
*
      real*4 shape1(36)	  ! shape of mass 1 (detector) in its cartesian
	  		  ! includes detector mode shape
      real*4 shape2(36)	  ! shape of mass 2 (source) in its cartesians
	  		  ! includes source mode shape
      real*4 D(3)	  ! source position relative to detector
	  		  ! in detector cartesians
      real*4 point1(3)    ! point in box 1; x, y, z
      real*4 point2(3)	  ! point in box 2; scaled R, theta, phi
      real*4 asinom	  ! amp*sin(omega)
      real*4 R		  ! radial coord, mass 2
      real*4 snth, csth	  ! sin(theta), cos(theta)
      real*4 snphi, csphi ! sin(phi), cos(phi)

      logical*4 inside	  ! .true. if both points are inside masses
*
* variables
*	  
      real*4  modez  ! mode shape function
      real*4  curve  ! curvature function
      real*4  sdelta ! source deflection
      real*4  ddelta ! detector deflection
      real*4  pn2(3) ! cartesian coords of point 2 in mass 2 cartesians
*
* compute point2 in source mass cartesians
*
      pn2(1)=point1(1)+R*snth*csphi-D(1)
      pn2(2)=point1(2)+R*snth*snphi-D(2)
      pn2(3)=point1(3)+R*csth-D(3)
*
*  deflections
*
      sdelta=asinom*modez(pn2,shape2)+curve(pn2,shape2)
      ddelta=curve(point1,shape1)
*
* return if outside
*
      inside=.false.
*     print *,'plates called; start plates',sdelta,ddelta
	  
      if (pn2(3) .gt. shape2(3)/2.0+sdelta) then	   ! source z max
      	return
      else if (pn2(3) .lt. -(shape2(3)/2.0)+sdelta) then   ! source z min
      	return
      else if (pn2(1) .gt. shape2(1)) then		   ! source x max
      	return
      else if (pn2(1) .lt. 0.0) then			   ! source x min
      	return
      else if (pn2(2) .gt. shape2(2)/2.0) then		   ! source y max
      	return
      else if (pn2(2) .lt. -shape2(2)/2.0) then		   ! source y min
   	return
      else if (point1(3) .gt. (shape1(3)/2.0)+ddelta) then  ! detector z max
      	return
      else if (point1(3) .lt. -(shape1(3)/2.0)+ddelta) then ! detector z min
      	return
      else if (point1(1) .gt. 0.0) then			    ! detector x max
      	return
      else if (point1(1) .lt. -shape1(1)) then		    ! detector x min
      	return
      else if (point1(2) .gt. shape1(2)/2.0) then	    ! detector y max
      	return
      else if (point1(2) .lt. -shape1(2)/2.0) then	    ! detector y min
      	return
*
      else
      	inside=.true.					    ! both points inside
      end if
*	
      return
*	
      end
*	
*** dice.f ***
*
* returns point1 inside the box defined by bounds1 and point2 inside
* the box defined by bounds2.
*
      subroutine dice(point1,point2,bounds1,bounds2)
      implicit none
*
* variables passed
*
      real*4 point1(3)	! point in box 1: x, y, z
      real*4 point2(3)	! point in box 2: R, theta, phi
      real*4 bounds1(6)	! box around mass 1, order is:
			!    x1 upper, x1 lower, y1 upper, y1 lower
			!    z1 upper, z1 lower
      real*4 bounds2(6)	! box around polar coordinates, order is:
			!    R upper, R lower, theta upper, theta lower
			!    phi upper, phi lower
*
* variables
*
      integer*4 randi       ! function for getting random integers
      integer*4 rand1,rand2 ! scratch random numbers	
      integer*4 i	    ! index
*
* throw the dice
*
      do 20 i=1,3
*      print *,'dice called; start dice loop' 
      	rand1=randi()
      	rand2=randi()
*     	print *, rand1,rand2
      	point1(i)=((bounds1(2*i-1)+bounds1(2*i))+
     > (bounds1(2*i-1)-bounds1(2*i))*(rand1/32767.0))/2.0
      	point2(i)=((bounds2(2*i-1)+bounds2(2*i))+
     > (bounds2(2*i-1)-bounds2(2*i))*(rand2/32767.0))/2.0
20    continue
      return
      end
*	  
***  rand.f  ***
*
*  function returns a random integer [-32768,+32767]
*  derived from Absoft example random.f
*
      FUNCTION randi()  
      implicit none
      REAL*4 rando, rand
      integer*4 randi
*
      rando=rand(0)
      randi = -32768 + INT(65536 * rando)
*
      END
*
