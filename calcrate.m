function dpdt = calcrate(t,y,const)
dpdt = [y(1)*const.r-const.c*y(2)*y(1);y(2)*(-const.d)+const.e*y(2)*y(1)];
end