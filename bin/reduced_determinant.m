
% Run the main function called 'reduced_determinant'
% with the filename as input
% This method will call the parser and the appropriate methods to covert
% the data and then calculate the reduced determinant.

% This file creates appropriate tags for HTML formatting.


function [reduced_determinant_val] = reduced_determinant(txtName)

% Setup the font.
fprintf(stdout(), "<p style='font-family: monospace;'>\n");

% Parse the file
[r, c, Gamma, V_num] = human_parser(txtName);
V = general_kinetics(r, c, V_num);
MA = mass_action(r, c, V_num);

fprintf(stdout(), "<br/><br/>    Reduced Determinant Using General Kinetics<br/><br/><br/>\n");
disp(reduced_determinant_calc(Gamma, V));

fprintf(stdout(), "<br/><br/>    Reduced Determinant Using Mass Action<br/><br/><br/>\n");
disp(reduced_determinant_calc(Gamma, MA));

% End the formatting.
fprintf(stdout(), "\n</p><br/><br/><br/>\n");
end




%% Completed
%   This function takes two arguments, the V and Gamma,
%   and computes the reduceddeterminant
function [reduced_determinant_val] = reduced_determinant_calc(Gamma, V)
r = rank(Gamma);
n = size(Gamma, 1);
m = size(Gamma, 2);
product = Gamma * V;
reduced_determinant_val = sym("0");
combos = nchoosek(1:n, r);

% take all different arrangements of prinicpal minors
for ii=1:length(combos)
  reduced_determinant_val = det(product(combos(ii, :),
          combos(ii, :))) + reduced_determinant_val;
end

end
