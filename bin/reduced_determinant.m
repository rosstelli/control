function [reduced_determinant_val] = reduced_determinant(txtName)
fprintf(stdout(), "<p style='font-family: monospace;'>");
[Gamma, V_num, V, MA] = human_parser(txtName);
%r = rank(Gamma)
%n = size(Gamma, 1);
%m = size(Gamma, 2);
%product = Gamma * V_num;
%product2 = Gamma * V
%reduced_determinant_val = sym("0");
%combos = nchoosek(1:n, r);
%for ii=1:length(combos)
%  reduced_determinant_val = det(product2(combos(ii, :), combos(ii, :))) + reduced_determinant_val;
%end
%fprintf(stdout(), "The Gamma is:\n");
%disp(Gamma);
%%fprintf(stdout(), "The V numeric is:\n");
%%disp(V_num);
%fprintf(stdout(), "The V transpose is:\n");
%disp(V);
%fprintf(stdout(), "The reduced_determinant is:");
%disp(reduced_determinant_val);

fprintf(stdout(), "<br/><br/>    Reduced Determinant Using General Kinetics<br/><br/><br/>");


disp(reduced_determinant_calc(Gamma, V));

fprintf(stdout(), "<br/><br/>    Reduced Determinant Using Mass Action<br/><br/><br/>");

disp(reduced_determinant_calc(Gamma, MA));
fprintf(stdout(), "</p>");
fprintf(stdout(), "<br/><br/><br/>");
end


function [reduced_determinant_val] = reduced_determinant_calc(Gamma, V)

r = rank(Gamma);
n = size(Gamma, 1);
m = size(Gamma, 2);
product = Gamma * V;
reduced_determinant_val = sym("0");
combos = nchoosek(1:n, r);
for ii=1:length(combos)
  reduced_determinant_val = det(product(combos(ii, :), combos(ii, :))) + reduced_determinant_val;
end
%fprintf(stdout(), "The Gamma is:\n");
%disp(Gamma);
%fprintf(stdout(), "The V transpose is:\n");
%disp(V);
%fprintf(stdout(), "The reduced_determinant is:");
%disp(reduced_determinant_val);



end
