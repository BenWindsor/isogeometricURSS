function val = nonUniformCircleYParam( t )

if t>=0 && t<0.5
    val=sin(3*pi*t); %first three quarders of circle 
else
    val=sin(pi*(t-0.5)+(3/2)*pi); %last quarter of circle
end
    

end

