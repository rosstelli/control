


% need the reactions, and compounds

function [Mass_Action] = mass_action(reactions, compounds, V_numeric)


num_of_reactions = length(reactions(:,1));
num_of_compounds = compounds.size();


# Calculate the MASS ACTION
original_reaction_concentrations = repmat(sym("1"), num_of_reactions,1);
for ii=1:num_of_reactions
   reaction = reactions(ii,:);
   #curr_reaction_concentrations = ""
   
   original_reaction_concentrations(ii) = original_reaction_concentrations(ii) * sym(["k_", num2str(ii), "^1"]);
   
   if isempty(findstr('-->', reaction)) 
     # Need the reactants from the left hand side
     rhs = reaction((findstr('<--', reaction)+3):length(reaction));
     spaces = [1 strfind(rhs, ' ') length(rhs)];
     # TODO: avoid + and other signs
     for jj=1:(length(spaces)-1) 
       substring = strrep(rhs((spaces(jj)):(spaces(jj+1))), " ", ""); # remove all the spaces
       if strcmp('',substring) == 0 && strcmp('+', substring) == 0 && !isempty(substring) #&& strcmp(' ',substring) == 0 && strcmp("",substring) == 0
         # get the coefficient
         count = 1; #find the coefficient
         while substring(count) <= '9' && substring(count) >= '0' && count <= length(substring)
             count = count + 1;
         end
         
         coefficient = (substring(1:(count-1)));
         if strcmp(coefficient, "") || isempty(coefficient)
             coefficient = '1';
         end
         substring = substring(count:length(substring));
         original_reaction_concentrations(ii) = original_reaction_concentrations(ii) * sym([substring "^" coefficient]);
       end # End if empty
     end # end for loop
   elseif isempty(findstr('<--', reaction))
     # Need the reactants from the left hand side
     lhs = reaction(1:(findstr('-->', reaction)-1));
     spaces = [1 strfind(lhs, ' ') length(lhs)];
     # TODO: avoid + and other signs
     for jj=1:(length(spaces)-1) 
       substring = strrep(lhs((spaces(jj)):(spaces(jj+1))), " ", ""); # remove all the spaces
       if strcmp('',substring) == 0 && strcmp('+', substring) == 0 #&& strcmp(' ',substring) == 0 && strcmp("",substring) == 0
         # get the coefficient
         count = 1; #find the coefficient
         while substring(count) <= '9' && substring(count) >= '0' && count <= length(substring)
             count = count + 1;
         end
         
         coefficient = (substring(1:(count-1)));
         if strcmp(coefficient, "") || isempty(coefficient)
             coefficient = '1';
         end
         substring = substring(count:length(substring));
         original_reaction_concentrations(ii) = original_reaction_concentrations(ii) * sym([substring "^" coefficient]);
       end # End if empty
     end # end for loop
   end# end checking sides
   
end # end creating original_reaction_concentrations vector
fprintf(stdout(), "<br/>");
original_reaction_concentrations
fprintf(stdout(), "<br/>");

Mass_Action = repmat(sym("0"), num_of_reactions, num_of_compounds);
for ii = 1:num_of_reactions
  for jj = 1:num_of_compounds
    Mass_Action(ii,jj) = sym(num2str(V_numeric(ii,jj))) * original_reaction_concentrations(ii);
  end
end

# Get the key set of compounds, i.e. the compound names
setOfCompounds = compounds.keySet();
iterator = setOfCompounds.iterator();

while iterator.hasNext()
  nextCompound = iterator.next();
  Mass_Action(:, compounds.get(nextCompound)) = diff(Mass_Action(:, compounds.get(nextCompound)), nextCompound);
end
fprintf(stdout(), "<br/>");
Mass_Action
fprintf(stdout(), "<br/>");
# END CALCULATE MASS ACTION



end