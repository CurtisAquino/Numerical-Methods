function J = derivative(f,x0,method)
        
        % ********************************
        % Corrects the dimensions of input
        % ********************************
        
%         f = {@(x,y,z) x*y*z^2;@(x,y,z) x*y*z^2};
%         x0 = [1,2,3];
        
        if ~isrow(x0)
            x0 = x0'; 
        end
        
        % ********************************************************
        % Determines the number of arguments per inputted function
        % ********************************************************
        
        if length(f) == 1
            V = nargin(f);
        else
            for i = 1:length(f)
                V(i) = nargin(f{i});
            end
        end
        
        % **************************************************
        % Terminates if variables are not properly specified
        % **************************************************
        
        if range(V) ~= 0 || length(x0) ~= V(1)
            fprintf('Error: incorrect number of variables')
            return
        end
        
        % ************************
        % Initializes the Jacobian 
        % ************************
        
        h       = (x0 ~= 0).*x0*10^(-8) + (x0 == 0).*10^(-16);
        x0      = num2cell(x0);  
        J       = zeros(length(f),V(1));
        
        % ***********************************
        % Default: 2-point forward difference
        % ***********************************
        
        if nargin < 3
            xh = repmat([x0{:}],V(1),1)+eye(V(1)).*h; 
            for i = 1:length(f)
                if length(f) == 1
                    F = f;
                else
                    F = f{i};
                end
                for j = 1:V(1)
                    x1 = num2cell(xh(j,:));
                    J(i,j) = F(x1{:})/h(j);
                end
                J(i,:) = J(i,:)-repmat(F(x0{:}),1,V(1))./h;
            end
        
        % ***************************
        % User has specified a method
        % ***************************
        
        else
            switch method
                case 'f2'
                    xh = repmat([x0{:}],V(1),1)+eye(V(1)).*h; 
                    for i = 1:length(f)
                        if length(f) == 1
                            F = f;
                        else
                            F = f{i};
                        end
                        for j = 1:V(1)
                            x1 = num2cell(xh(j,:));
                            J(i,j) = F(x1{:})/h(j);
                        end
                        J(i,:) = J(i,:)-repmat(F(x0{:}),1,V(1))./h;
                    end
                case 'b2'
                    xh1 = repmat([x0{:}],V(1),1)+eye(V(1)).*h; 
                    xh2 = repmat([x0{:}],V(1),1)-eye(V(1)).*h; 
                    for i = 1:length(f)
                        if length(f) == 1
                            F   = f;
                        else
                            F   = f{i};
                        end
                        for j = 1:V(1)
                            x11     = num2cell(xh1(j,:));
                            x12     = num2cell(xh2(j,:));
                            J(i,j)  = (F(x11{:})-F(x12{:}))/(2*h(j));
                        end
                    end
                case 'f3'
                    xh1 = repmat([x0{:}],V(1),1)+2*eye(V(1)).*h; 
                    xh2 = repmat([x0{:}],V(1),1)+eye(V(1)).*h; 
                    for i = 1:length(f)
                        if length(f) == 1
                            F   = f;
                        else
                            F   = f{i};
                        end
                        for j = 1:V(1)
                            x11     = num2cell(xh1(j,:));
                            x12     = num2cell(xh2(j,:));
                            J(i,j)  = (-F(x11{:})+4*F(x12{:}))/(2*h(j));
                        end
                        J(i,:)      = J(i,:)-3*repmat(F(x0{:}),1,V(1))./(2*h);
                    end
                case 'b3'
                    xh1 = repmat([x0{:}],V(1),1)+2*eye(V(1)).*h; 
                    xh2 = repmat([x0{:}],V(1),1)+eye(V(1)).*h; 
                    for i = 1:length(f)
                        if length(f) == 1
                            F   = f;
                        else
                            F   = f{i};
                        end
                        for j = 1:V(1)
                            x11     = num2cell(xh1(j,:));
                            x12     = num2cell(xh2(j,:));
                            J(i,j)  = (3*F(x11{:})-4*F(x12{:}))/(2*h(j));
                        end
                        J(i,:)      = J(i,:)+repmat(F(x0{:}),1,V(1))./(2*h);
                    end
            end
        end
end