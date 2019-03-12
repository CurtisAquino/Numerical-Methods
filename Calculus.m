classdef Calculus
    methods(Static)
        %%
function D = Hessian(f,vars,Order,x0)
            
        %f = @(x,y) x.^2+y;     
        %Order = 2;    
        %x0 = [1,2];
        
        % Initialization
        g = f;
        if isrow(x0); x0 = x0'; end
        numVar      = nargin(g);
        if length(vars) ~= numVar; sprintf('Error'); return; end
        
        AtHandVar = '@('; AtHandDif = '@(a1,'; 
        HandVar   = '('; HandDif = '(a1,'; 
        for i = 1:numVar
            if i == numVar
                AtHandVar   = [AtHandVar,vars{i},')'];  
                AtHandDif   = [AtHandDif,vars{i},')'];
                HandVar     = [HandVar,vars{i},')'];
                HandDif     = [HandDif,vars{i},')'];
            else 
                AtHandVar   = [AtHandVar,vars{i},','];  
                AtHandDif   = [AtHandDif,vars{i},','];
                HandVar     = [HandVar,vars{i},','];
                HandDif     = [HandDif,vars{i},','];
            end
        end
        
        syms a1
        for i = 1:numVar
            dy{i} = matlabFunction(subs(f,vars{i},{vars{i}-a1}));
            eval(sprintf('d{%i} = %s (g %s - dy{%i} %s)./a1;',i,AtHandDif,HandVar,i,HandDif));
        end
        
        h           = sqrt((x0~=0).*x0.*10^(-16))+(x0==0).*10^(-16);
        
        for i = 1:numVar
            temp = num2cell([h(i),x0']);
            D(i) = d{i}(temp{:});
        end
        
        
        
        end
    
        
    end
end