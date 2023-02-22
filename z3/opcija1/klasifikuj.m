function klasa = klasifikuj(x, V, v0)
if ((V'*x + v0) > 0) 
    klasa = 1;
else
    klasa = -1;
end
end