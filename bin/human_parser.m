# input is a filename

# output is left stoichiometric matrix 
# right stoichiometric matrix 

function [reactions, compounds, Gamma, V_numeric] = human_parser(txtName)
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



end
