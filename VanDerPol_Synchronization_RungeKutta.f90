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
double precision,parameter::tfin=100.0d0	!tをどこまで計算するか

double precision,parameter::x1ini=2.0d0		!振動子1のx初期値(位置)
double precision,parameter::v1ini=0.0d0		!振動子1のv初期値(初速度)

double precision,parameter::x2ini=0.0d0		!振動子2のx初期値(位置)
double precision,parameter::v2ini=0.0d0		!振動子2のv初期値(初速度)

double precision,parameter::lamda=1.0d0		!定数lamdaの値(摩擦の大きさ)
double precision,parameter::alpha=0.0d0		!相互作用の強さを表すパラメータ
double precision,parameter::omega1=1.0d0	!振動子1の周波数
double precision,parameter::omega2=0.99d0	!振動子2の周波数

integer,parameter::tNbin=1000				!t時間の刻み数

integer::i
double precision::t,x1,v1,x2,v2
double precision::dt
double precision::k1,k2,k3,k4,q1,q2,q3,q4
double precision::l1,l2,l3,l4,r1,r2,r3,r4

!character*5,parameter::INPX='X0000'
!character*8,parameter::INPVini='Vini0010'
!character*9,parameter::INPLambda='Lamda0001'
character*9,parameter::INPALPHA='Alpha0000'

open(10,file='RungeKutta_VanderPol_Synchronization'//INPALPHA)

dt=(tfin-tini)/dble(tNbin)		!刻む幅を計算
t=tini 							!初期値を代入
x1=x1ini							!初期値を代入
v1=v1ini 							!初期値を代入
x2=x2ini
v2=v2ini


!ルンゲクッタ法の計算をする
do i=1,tNbin
	k1=dt*funcdy1(t,x1,v1,x2,v2)
	q1=dt*funcdz1(t,x1,v1,x2,v2)
	l1=dt*funcdy2(t,x1,v1,x2,v2)
	r1=dt*funcdz2(t,x1,v1,x2,v2)

	k2=dt*funcdy1(t+0.50d0*dt,x1+0.50d0*k1,v1+0.50d0*q1,x2+0.50d0*l1,v2+0.50d0*r1)
	q2=dt*funcdz1(t+0.50d0*dt,x1+0.50d0*k1,v1+0.50d0*q1,x2+0.50d0*l1,v2+0.50d0*r1)
	l2=dt*funcdy2(t+0.50d0*dt,x1+0.50d0*k1,v1+0.50d0*q1,x2+0.50d0*l1,v2+0.50d0*r1)
	r2=dt*funcdz2(t+0.50d0*dt,x1+0.50d0*k1,v1+0.50d0*q1,x2+0.50d0*l1,v2+0.50d0*r1)
	
	k3=dt*funcdy1(t+0.50d0*dt,x1+0.50d0*k2,v1+0.50d0*q2,x2+0.50d0*l2,v2+0.50d0*r2)
	q3=dt*funcdz1(t+0.50d0*dt,x1+0.50d0*k2,v1+0.50d0*q2,x2+0.50d0*l2,v2+0.50d0*r2)
	l3=dt*funcdy2(t+0.50d0*dt,x1+0.50d0*k2,v1+0.50d0*q2,x2+0.50d0*l2,v2+0.50d0*r2)
	r3=dt*funcdz2(t+0.50d0*dt,x1+0.50d0*k2,v1+0.50d0*q2,x2+0.50d0*l2,v2+0.50d0*r2)
	
	k4=dt*funcdy1(t+dt,x1+k3,v1+q3,x2+l3,v2+r3)
	q4=dt*funcdz1(t+dt,x1+k3,v1+q3,x2+l3,v2+r3)
	l4=dt*funcdy2(t+dt,x1+k3,v1+q3,x2+l3,v2+r3)
	r4=dt*funcdz2(t+dt,x1+k3,v1+q3,x2+l3,v2+r3)

	x1=x1+(k1+2d0*k2+2d0*k3+k4)/6d0
	v1=v1+(q1+2d0*q2+2d0*q3+q4)/6d0
	x2=x2+(l1+2d0*l2+2d0*l3+l4)/6d0
	v2=v2+(r1+2d0*r2+2d0*r3+r4)/6d0    

    t=tini+dt*dble(i)

    write(10,*) t,x1,v1,x2,v2
end do

close(10)

!-----------------------------------------

stop
contains


!振動子1についての式のdy/dx一階微分
function funcdy1(T,X1,V1,X2,V2) result(dydx)		
	double precision,intent(in)::T,X1,V1,X2,V2
	double precision dydx
	dydx=V1  		
end function funcdy1

!!振動子1についての式のdz/dx二階微分
function funcdz1(T,X1,V1,X2,V2) result(dzdx)
	double precision,intent(in)::T,X1,V1,X2,V2
	double precision dzdx
	dzdx=-lamda*(X1*X1-1)*V1-omega1*omega1*X1+alpha*(X2-X1)  		
end function funcdz1

!!振動子2についての式のdy/dx一階微分
function funcdy2(T,X1,V1,X2,V2) result(dydx)		
	double precision,intent(in)::T,X1,V1,X2,V2
	double precision dydx
	dydx=V2  		
end function funcdy2

!!振動子2についての式のdz/dx二階微分
function funcdz2(T,X1,V1,X2,V2) result(dzdx)
	double precision,intent(in)::T,X1,V1,X2,V2
	double precision dzdx
	dzdx=-lamda*(X2*X2-1)*V2-omega2*omega2*X2+alpha*(X1-X2)  		
end function funcdz2

end program main