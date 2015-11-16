!-------------------
!
!20150820
!
!runge_kutta_method
!van der Pol equation
!
!funcdy:dy/dx=F1(x,y,z)速度にあたる
!             F1(T,X,V)
!funcdz:dz/dx=F2(x,y,z)解きたい2次の微分方程式
!             F2(T,X,V)
!-------------------

program main
implicit none

!set parameter
double precision,parameter::tini=0.0d0		!t初期値(時間)
double precision,parameter::xini=0.0d0		!x初期値(位置)
double precision,parameter::vini=10.d0		!v初期値(初速度)
double precision,parameter::lamda=0.1d0		!定数lamdaの値

double precision,parameter::tfin=100.0d0		!xをどこまで計算するか
integer,parameter::tNbin=100000				!刻み数

integer::i
double precision::t,x,v
double precision::dt
double precision::k1,k2,k3,k4,q1,q2,q3,q4

character*5,parameter::INPX='X0000'
character*8,parameter::INPVini='Vini0010'
character*9,parameter::INPLambda='Lamda00k1'

open(10,file='RungeKutta_van_'//INPX//INPLambda//INPVini)

dt=(tfin-tini)/dble(tNbin)		!刻む幅を計算
t=tini 							!初期値を代入
x=xini							!初期値を代入
v=vini 							!初期値を代入

!ルンゲクッタ法の計算をする
do i=1,tNbin
	k1=dt*funcdy(t,x,v)
	q1=dt*funcdz(t,x,v)

	k2=dt*funcdy(t+0.50d0*dt,x+0.50d0*k1,v+0.50d0*q1)
	q2=dt*funcdz(t+0.50d0*dt,x+0.50d0*k1,v+0.50d0*q1)
	
	k3=dt*funcdy(t+0.50d0*dt,x+0.50d0*k2,v+0.50d0*q2)
	q3=dt*funcdz(t+0.50d0*dt,x+0.50d0*k2,v+0.50d0*q2)
	
	k4=dt*funcdy(t+dt,x+k3,v+q3)
	q4=dt*funcdz(t+dt,x+k3,v+q3)

	x=x+(k1+2d0*k2+2d0*k3+k4)/6d0
	v=v+(q1+2d0*q2+2d0*q3+q4)/6d0
    t=tini+dt*dble(i)
    write(10,*) t,x,v
end do

close(10)

!-----------------------------------------

stop
contains


!dy/dx
function funcdy(T,X,V) result(dydx)		
	double precision,intent(in)::T,X,V
	double precision dydx

	dydx=V  		

end function funcdy

!dz/dx
function funcdz(T,X,V) result(dzdx)
	double precision,intent(in)::T,X,V
	double precision dzdx

	dzdx=-lamda*(X*X-1)*V-X  		

end function funcdz

end program main