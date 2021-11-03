%--------------------------------------------------------------------------
%                           Ejercicio 1 :
%--------------------------------------------------------------------------
f0=2;
n_points=1001 % Nuestro número de muestras.
t = linspace(-1,1,n_points); % generamos n_points numero de muestras entre -1 y 1.
xc0=cos(t*f0*2*pi);
figure(1);
hold on;
plot(t,xc0,'b');
[xd0,ts0]=conv_cd(xc0,t,0.1);
xlabel('time'), ylabel('y(t)') ;
stem(ts0,xd0,'r') ;
hold off

%--------------------------------------------------------------------------
%                           Ejercicio 2 :
%--------------------------------------------------------------------------

%Meta : calcular las transformadas de fouried de una misma señal en forma
%continua y discreta .

hold on;
figure(2);  
%  1º Generamos xc(t) : t € [-2,2] con una variacion sobre t de 0,001 .
f0=2;
muestras=1001;
Ts1=0.1;Ts2=0.01;
vt=0.001;
t=-2:vt:2;
xc = cos(2*pi*f0*t); %Nuestra funcion senoidal con frecuencia 2 .

% 2º Evaluamos Xc(jw) : w € [-6pi,6pi].

w = linspace (-6*pi,6*pi,muestras)

for n = 1:length(w);
   Xc(n) = sum(xc.*exp(-j*w(n)*t)*vt);
end 

[xd,ts0]=conv_cd(Xc,w,0.1); %Transformada de foulier de xc.
subplot(3,1,1);
plot(ts0,xd); 
%3º Muestrear la señal Xc(t) con T1=0.1 obteniendo Xd [n]=xc(nTs) entre w €
% [-pi,pi] con 1001 puntos.
wd= linspace(-pi,pi,muestras); 
%Calculamos Xd [t] con xc (t) segundo con un periodo de 0.1.
[xd1,tsampl]=conv_cd(xc,t,Ts1);


%Calculamos Xd [t] con xc (t) segundo con un periodo de 0.01.
[xd2,tsamp3]=conv_cd(xc,t,Ts2);
%Normalizamos el numero de muestras   
n1=-floor(size(tsampl,2)/2):floor(size(tsampl,2)/2);
n2=-floor(size(tsamp3,2)/2):floor(size(tsamp3,2)/2);

%Calculamos la transformada de fulier de la señal con periodo 0,1.
for k =1:length(wd)
    Xd1(k)=sum(xd1.*exp(-j*n1*wd(k)));
end 

subplot(3,1,2)
plot(wd,real(Xd1)) % la ponemos en segundo lugar de la figura

%Calculamos la transformada de fulier de la señal con periodo 0,01.

for k =1:length(wd)
    Xd2(k)=sum(xd2.*exp(-j*n2*wd(k)))
end   

subplot(3,1,3)
plot(wd,real(Xd2))%La ponemos en el ultimo lugar de la figura.
hold off;
%--------------------------------------------------------------------------
%                           Ejercicio 3 :
%--------------------------------------------------------------------------
% 1º Calculamos las señales reconstruidas para tener todos los
% ingredientes.

for k =1:length(t)
    hil=filtro_ideal(t(k)-tsampl,Ts1); %Para la convolucion [xd1,tsamp1] ; T=0,1
    xr1(k)=sum(xd1.*hil);
end

for k =1:length(t)
    hil=filtro_ideal(t(k)-tsamp3,Ts2);%Para la convolucion [xd1,tsamp1] ; T=0,01
    xr2(k)=sum(xd2.*hil);
end
%Comenzamos a plotear :
%Nos pide dos figuras , estas son:


%   (2,1,1) --> T = 0,1
figure(3);

subplot(2,1,1); 
plot(t,xc,'b'); % Representamos la señal senoidal en el tiempo.
hold on;
stem(tsampl,xd1,'r') ; % Le añadimos las muestras discretas de la señal con T=0,1
plot(t,xr1,'g'); % Añadimos en verde la reconstruccion de la señal. 
hold off;

%   (2,1,2) --> T = 0,01

subplot(2,1,2);
plot(t,xc,'b');  % Representamos la señal senoidal en el tiempo.
hold on;
stem(tsamp3,xd2,'r') ; % Le añadimos las muestras discretas de la señal con T=0,1
plot(t,xr2,'g'); % Añadimos en verde la reconstruccion de la señal. 
hold off;
%--------------------------------------------------------------------------
%                           Ejercicio 4 :
%--------------------------------------------------------------------------
%
%
%
%--------------------------------------------------------------------------
%                            APARTADO A 
%--------------------------------------------------------------------------
for k =1:length(t)
    hil=filtro_orden0(t(k)-tsampl,Ts1); %Para la convolucion [xd1,tsamp1] ; T=0,1
    xr1(k)=sum(xd1.*hil);
end
figure(3);

for k =1:length(t)
    hil=filtro_orden0(t(k)-tsamp3,Ts2);%Para la convolucion [xd1,tsamp1] ; T=0,01
    xr2(k)=sum(xd2.*hil);
end

% 2º Comenzamos a plotear :
%Nos pide dos figuras , estas son:


%   (2,1,1) --> T = 0,1
figure(4);

subplot(2,1,1); 
plot(t,xc,'b'); % Representamos la señal senoidal en el tiempo.
hold on;
stem(tsampl,xd1,'r') ; % Le añadimos las muestras discretas de la señal con T=0,1
plot(t,xr1,'g'); % Añadimos en verde la reconstruccion de la señal. 
hold off;

%   (2,1,2) --> T = 0,01

subplot(2,1,2);
plot(t,xc,'b');  % Representamos la señal senoidal en el tiempo.
hold on;
stem(tsamp3,xd2,'r') ; % Le añadimos las muestras discretas de la señal con T=0,1
plot(t,xr2,'g'); % Añadimos en verde la reconstruccion de la señal. 
hold off;

%--------------------------------------------------------------------------
%                            APARTADO B 
%--------------------------------------------------------------------------

% 1º Obtenemos las señales reconstruidas utilizando el interpolador lineal. 

for k =1:length(t)
    hil=filtro_lineal(t(k)-tsampl,Ts1); %Para la convolucion [xd1,tsamp1] ; T=0,1
    xr1(k)=sum(xd1.*hil);
end
figure(3);

for k =1:length(t)
    hil=filtro_lineal(t(k)-tsamp3,Ts2);%Para la convolucion [xd1,tsamp1] ; T=0,01
    xr2(k)=sum(xd2.*hil);
end

% 2º Comenzamos a plotear :
%Nos pide dos figuras , estas son:


%   (2,1,1) --> T = 0,1
figure(5);

subplot(2,1,1); 
plot(t,xc,'b'); % Representamos la señal senoidal en el tiempo.
hold on;
stem(tsampl,xd1,'r') ; % Le añadimos las muestras discretas de la señal con T=0,1
plot(t,xr1,'g'); % Añadimos en verde la reconstruccion de la señal. 
hold off;

%   (2,1,2) --> T = 0,01

subplot(2,1,2);
plot(t,xc,'b');  % Representamos la señal senoidal en el tiempo.
hold on;
stem(tsamp3,xd2,'r') ; % Le añadimos las muestras discretas de la señal con T=0,1
plot(t,xr2,'g'); % Añadimos en verde la reconstruccion de la señal. 
hold off;

%--------------------------------------------------------------------------