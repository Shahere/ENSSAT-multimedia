 function [Zout,y] = filtrage(alpha11,alpha12,x,Zin)
    [y,Zout] = filter(alpha11,alpha12,x,Zin);