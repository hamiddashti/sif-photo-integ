clear
clc
x0 = 1.5
x2=5;
x3=6;
func = @(x) f(x,x2,x3)
[ci_final, err2, fcounter] = fixedp_brent_ari(func,x0)

function [err, y,y2] = f(x,x2,x3)
y=x^2-3*x+4
y2=x2+x3;
err = y-x;
end
