
function [] = tmp ()
# test
pkg load symbolic
syms a b c
varname = genvarname("g")
# syms ["d"; "e"; "f"]
disp(varname);
syms varname
x = [2*varname b; b c] * [1 0 ; 1 1]
disp(x)

y = [ sym("x_1") 0; 0 1]
disp(y);
end