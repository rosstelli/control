function [V_gen] = general_kinetics (reactions, compounds, V_numeric)

reactions(1,:)
num_of_reactions = length(reactions(:,1))
num_of_compounds = compounds.size();

counter = 1;
#V = zeros(num_of_reactions, num_of_compounds);
V_gen = repmat(sym("0"), num_of_reactions, num_of_compounds);
for ii=1:num_of_reactions
    for jj=1:num_of_compounds
        if V_numeric(ii,jj)
            V_gen(ii,jj) = sym(["x_", num2str(counter), "^1"]);
            counter = counter + 1;
        else 
            V_gen(ii,jj) = sym("0");
        end
    end
end

end