

function [d_vec] = derivative_vector(filename)

[compounds, lhs, rhs] = human_parser(filename);
Gamma = rhs - lhs;

v_dot = calculate_vdot(lhs);

d_vec = Gamma * v_dot;


end