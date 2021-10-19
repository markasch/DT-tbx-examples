function xdot=duffing(t,x)

global gamma epsilon GAM OMEG

xdot(2)=-gamma*x(2)-epsilon*x(1)^3+GAM*cos(OMEG*t);
xdot(1)=x(2);

xdot=xdot';

end