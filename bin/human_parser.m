function [Gamma, V_numeric, V, Mass_Action] = human_parser(txtName)
pkg load symbolic;
# [ K, S, reactants, targets, gamma]

num_of_reactions = 0;
num_of_compounds = 0;
reactions = [""];
compounds = javaObject("java.util.Hashtable");

fid = fopen(txtName);

# Go through the file, line by line
tline = fgets(fid);
while ischar(tline)
    disp("line:")
    fprintf(stdout(), "<br/>");
    disp(tline)
    fprintf(stdout(), "<br/>");
    if isempty(tline) tline = fgets(fid); continue end
    num_of_reactions = num_of_reactions + 1;
    single_reaction = true;
    
    
    tline = strrep(tline, "\n", "")
    spaces = [1 strfind(tline, ' ') length(tline)];
    for ii=1:(length(spaces)-1) 
        
        substring = strtrim(tline((spaces(ii)):(spaces(ii+1))));
    
        if strcmp(substring, "") continue end
        #remove coefficients
        count = 1; #find the coefficient
        while substring(count) <= '9' && substring(count) >= '0' && count <= length(substring)
            count = count + 1;
        end
        #Remove the coefficient you found.
        substring = substring(count:length(substring));
        
        if ~(strcmp(substring, "-->") || strcmp(substring, "<-->") || strcmp(substring, "<--") || strcmp(substring, "+")) && isempty(compounds.get(substring))
            num_of_compounds = num_of_compounds + 1;
            compounds.put(substring, num_of_compounds);
        elseif (strcmp(substring, "<-->"))
            single_reaction = false;
            num_of_reactions = num_of_reactions + 1;
            reactions = [reactions; regexprep(tline, '<', ''); regexprep(tline, '>', '')];
        end 
 
    end

    if single_reaction
        reactions = [reactions; tline];
    end
    tline = fgets(fid); # Get the next line.
end

fclose(fid);

#reactions = cellstr(reactions);

# Calculate Matrices
Gamma = zeros(num_of_compounds, num_of_reactions);
V_numeric = zeros(num_of_reactions, num_of_compounds);

for rxn=1:num_of_reactions
    curr = reactions(rxn,:);
    spaces = [1 strfind(curr, ' ') length(curr)];
    side_G=1;
    side_V=1;
    for ii=1:(length(spaces)-1) 
        
        substring = strtrim(curr((spaces(ii)):(spaces(ii+1))));

        if strcmp(substring, "") continue end
        #remove coefficients
        count = 1; #find the coefficient
        while substring(count) <= '9' && substring(count) >= '0' && count <= length(substring)
            count = count + 1;
        end
        
        #Get the coefficient
        coefficient = (substring(1:(count-1)));
        if strcmp(coefficient, "") || isempty(coefficient)
            coefficient = '1';
        end
        coefficient = str2num(coefficient);
        
        #Remove the coefficient you found.
        substring = substring(count:length(substring));
        
        compound = compounds.get(substring);
        if ~isempty(compound)
            #compound = compounds.get(substring);
            Gamma(compound,rxn) = (side_G * coefficient);
            V_numeric(rxn,compound) = side_V;
        elseif strcmp(substring, "-->")
            side_V=0;
            Gamma(:,rxn)=-1*Gamma(:,rxn) + 0;
        elseif strcmp(substring, "<--")
            side_G=-1;
            V_numeric(rxn,:)=0;
        elseif strcmp(substring, "<-->")
            fprintf(stdout(), "<br/>Warning: Attempting to parse reversible reactions.<br/>");
        end
    end
end

counter = 1;
#V = zeros(num_of_reactions, num_of_compounds);
V = repmat(sym("0"), num_of_reactions, num_of_compounds);
for ii=1:num_of_reactions
    for jj=1:num_of_compounds
        if V_numeric(ii,jj)
            V(ii,jj) = sym(["x_", num2str(counter), "^1"]);
            counter = counter + 1;
        else 
            V(ii,jj) = sym("0");
        end
    end
end


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
  #disp(nextCompound)
  #disp(compounds.get(nextCompound))
  #disp(diff(Mass_Action(:,compounds.get(nextCompound)), nextCompound))
  Mass_Action(:, compounds.get(nextCompound)) = diff(Mass_Action(:, compounds.get(nextCompound)), nextCompound);
end
fprintf(stdout(), "<br/>");
Mass_Action
fprintf(stdout(), "<br/>");
# END CALCULATE MASS ACTION

end






