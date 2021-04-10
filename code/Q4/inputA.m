function val=inputA(t)
% piecewise linear signal
%----------------------------------------------------------
% Definitions:
[t1] = deal(0.03);

if t<t1   %Seg#1
    val=0;
else
    val=1;
end

end 
