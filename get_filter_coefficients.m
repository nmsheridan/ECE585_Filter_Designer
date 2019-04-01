function [G, err] = get_filter_coefficients(N,type)

%This functions serves the sole purpose of returning the filter
%coefficiencts for the LPF prototypes

G = NaN;
err = 'ERROR!';

if((N>10)||(isempty(N))||(isnan(N)))
   err = 'Filter Order Out of Bounds! (must be between 1 & 10)';
   return
end

if(strcmp(type,'butterworth'))
    err = '';
    T =[2     ,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        1.4142,1.4142,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        1     ,2     ,1     ,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        0.7654,1.8478,1.8478,0.7654,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        0.618 ,1.618 ,2     ,1.681 ,0.618 ,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        0.5176,1.4142,1.9318,1.9318,1.4142,0.5176,1     ,NaN   ,NaN   ,NaN   ,NaN;...
        0.445 ,1.247 ,1.8019,2     ,1.8019,1.247 ,0.445 ,1     ,NaN   ,NaN   ,NaN;...
        0.3902,1.1111,1.6629,1.9615,1.9615,1.6629,1.1111,0.3902,1     ,NaN   ,NaN;...
        0.3473,1     ,1.5321,1.8794,2     ,1.8794,1.5321,1     ,0.3473,1     ,NaN;...
        0.3129,0.908 ,1.4142,1.782 ,1.9754,1.9754,1.782 ,1.4142,0.908 ,0.3129,1];   
end

if(strcmp(type,'chebyshev05'))
    err = '';
    T =[0.6986,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        1.4029,0.7071,1.9841,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        1.5963,1.0967,1.5963,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        1.6703,1.1926,2.3661,0.8419,1.9841,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        1.7058,1.2296,2.5408,1.2296,1.7058,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        1.7254,1.2479,2.6064,1.3137,2.4758,0.8696,1.9841,NaN   ,NaN   ,NaN   ,NaN;...
        1.7372,1.2583,2.6381,1.3444,2.6381,1.2583,1.7372,1     ,NaN   ,NaN   ,NaN;...
        1.7451,1.2647,2.6564,1.359 ,2.6964,1.3389,2.5093,0.8796,1.9841,NaN   ,NaN;...
        1.7504,1.269 ,2.6678,1.3673,2.7239,1.3673,2.6678,1.269 ,1.7504,1     ,NaN;...
        1.7543,1.2721,2.6754,1.3725,2.7392,1.3806,2.7231,1.3485,2.5238,0.8842,1.9841];
end

if(strcmp(type,'chebyshev3'))
    err = '';
    T =[1.9953,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        3.1013,0.5339,5.8095,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        3.3487,0.7117,3.3487,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        3.4389,0.7483,4.3471,0.592 ,5.8095,NaN   ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        3.4817,0.7618,4.5381,0.7618,3.4817,1     ,NaN   ,NaN   ,NaN   ,NaN   ,NaN;...
        3.5045,0.7685,4.6061,0.7929,4.4641,0.6033,5.8095,NaN   ,NaN   ,NaN   ,NaN;...
        3.5182,0.7723,4.6386,0.8039,4.6386,0.7723,3.5182,1     ,NaN   ,NaN   ,NaN;...
        3.5277,0.7745,4.6575,0.8089,4.699 ,0.8018,4.499 ,0.6073,5.8095,NaN   ,NaN;...
        3.534 ,0.776 ,4.6692,0.8118,4.7272,0.8118,4.6692,0.776 ,3.534 ,1     ,NaN;...
        3.5384,0.7771,4.6768,0.8136,4.7425,0.8164,4.726 ,0.8051,4.5142,0.6091,5.8095];
end

G = T(N,1:(N+1));