
%%*****************************************************************************
 %  hybrid plasmonic ring resonator in 2.5D by siavash mirzaei
%%*****************************************************************************
% in this code we analysis hybrid plasmonics dielectric losded ring resonator in  it's cross
% section. This method is named as Body of revolution-FDTD (BOR-FDTD). we consider axisymmetric strustures.the core is silicon, the width of silicon is 200 nm and the
% height of it, is  150 nm. The Ag with height of 40 nm is considered as substrate.
% a mat=rial with low index as sio2  is placed between silicon and Ag.the hybrid
% plasmonic modes propagate in low index.
% for dispersive model of silicon we use three resonance of lorentz model
% for permitivity model of Ag we use drude-lorentz model.one term of drude model and two term of
% lorentz model is used. the range of frequency is
% visible(300nm -900 nm).the boundary condition is PEC.the spatial
% resolution is 1 nm.the meshgrid is uniform using fine meshes.

  
clear all
close all
clc   
%***********************************************************************  
%     Fundamental constants  
%***********************************************************************     
c0=2.99792458e8;            %speed of light in free space  
mu0=4.0*pi*1.0e-7;          %permeability of free space  
eps0=1.0/(c0*c0*mu0);       %permittivity of free space  
%***********************************************************************  
%     Grid parameters       
%***********************************************************************  
N_r = 200;            %number of grid cells in x-direction  
N_z = 150;            %number of grid cells in z-direction    
N_t = 3000000;          %total number of time steps  

d_r = 2e-9; d_z = d_r;    % space resolution
R_ring = 1e-6;
n_rc= N_r/2;              % the center of meshgrid in r direction
n_zc = N_z/2;             % the center of meshgrid in z direction
%****************************************************
% metal parameters
%****************************************************
width_Ag = N_r*d_r; 
Height_Ag = 60e-9;
n_z1= round(Height_Ag/d_z);

%************************************************
% silicon-stipe
%************************************************
width_si_stripe= 200e-9; Height_si_stipe= 10e-9;
                        
n_r2=round(n_rc-(width_si_stripe)/(2*d_r));
n_r3=round(n_rc+(width_si_stripe)/(2*d_r));
n_z2= round((Height_si_stipe)/d_z)+ n_z1;

%***************************************************
 % dielectric (si,sio2) parameters in geometrical
%*************************************************** 
width_si_core= 200e-9; Height_si_core= 150e-9;
gap= 30e-9;
n_r1= round(n_rc-(width_si_core)/(2*d_r));
n_r4= n_rc+(width_si_core)/(2*d_r);
n_z3= (gap/d_z)+n_z2;
n_z4= round((Height_si_core/d_z)+ n_z3);  

%*********************************************
% observation point
%*********************************************
n_r_point_1= n_rc;          %observation point
n_z_point_1=round((n_z2+n_z3)/2);   %observation point

n_r_point_2= n_rc;
n_z_point_2=n_zc;

n_r_point_3= n_r3+5;
n_z_point_3= n_z2;

%*************************************************
% azimuthal mode number
%**************************************************
m_azimuthal= 23;                        
%*******************************************************
%stability condition
%*******************************************************
% d_t= 0.99*(1/((sqrt((1/d_r)^2+(m/(d_r))^2+(1/d_z)^2))*c0)); 
d_t= 0.99*d_r/((1+m_azimuthal)*c0);

%******************************************************************
% drude-lorentz model parameters for Ag
%******************************************************************
eps_inf= 1;
omega_d= 1.3877e16;
gama_d= 1.87e13;
omega_L1= 3.254e15;
omega_L2= 7.758e15;
gama_L1= 1.165e15;
gama_L2= 3.46e14;
delta_epsL1= 0.089;
delta_epsL2= 2.066;
%******************************************************************
% lorentz model parameters for silicon
%******************************************************************
cc0= 3e14*2*pi;
eps_infinite= 1;
omega_Lorentz1= 3.64*cc0;
omega_Lorentz2= 2.76*cc0;
omega_Lorentz3= 1.73*cc0;
gama_Lorentz1= 0;
gama_Lorentz2= 0.063*cc0;
gama_Lorentz3= 2.5*cc0;
delta_epsiLon1= 8;
delta_epsiLon2= 2.85;
delta_epsiLon3= -0.107;
%******************************************************
% CPML parameters
%******************************************************
gradingorder= 5; ma= 1; 
k= 1.2;
epsr= 3;
sigma_r_opt = (0.8*(gradingorder+1)./(d_r*((mu0/eps0)*epsr)^0.5)); 
sigma_z_opt = (0.8*(gradingorder+1)./(d_z*((mu0/eps0)*epsr)^0.5)); 
sigma_r_max= k.*sigma_r_opt;
sigma_z_max= k.*sigma_z_opt;
alpha_r_max = 0.24; 
alpha_z_max = alpha_r_max;    
kappa_r_max = 2.0;
kappa_r_max_prim = 5.0;
kappa_z_max = kappa_r_max;  

%****************************************************************************
 % Material parameters 
%*********************************************************************** 

n=ones(N_r,N_z)*1.4519; % refractive index of air
epsilonr= n.*n;
c= c0./n; 
%************************************************
% frequency range
%************************************************
% source parameters
% fmin_1= 100e12; fmax_1= 750e12;
% f= linspace(fmin_1,fmax_1,500);
% f0=(fmin_1+fmax_1)/2;
% lambda0= c0/f0;
%***************************************************
%  filter diagonalization parmeters
%***************************************************
fmin= 300e12;       % minimum frequency
fmin_1=fmin;
fmax= 360e12;       % maximum frequency
fmax_1=fmax;
Fs= N_t/(N_t*d_t);   % sapling frequency rate        
eps= 1e-2;
L= N_t-1;              % length of signal
f= linspace(fmin_1,fmax_1,500);
f0=(fmin_1+fmax_1)/2;
lambda0= c0/f0;

 %******************************************************
 %  coeificients
 %******************************************************
 % coeficient of magnetic field update
 ca= d_t/d_r;
 cb= d_t/d_z;
 cc= (d_t*m_azimuthal);
 cz= d_t/(mu0*d_z);
 cr= d_t/(mu0*d_r);
 ct=(d_t*m_azimuthal)/mu0;
 ck= d_t/2;
%***************************************************************
 % coeficients in core (si) region
%***************************************************************
alpha1=(1-gama_Lorentz1*d_t)/(1+gama_Lorentz1*d_t);
alpha2=(1-gama_Lorentz2*d_t)/(1+gama_Lorentz2*d_t);
alpha3=(1-gama_Lorentz3*d_t)/(1+gama_Lorentz3*d_t);
beta1=(((omega_Lorentz1)^2)*d_t)/(1+gama_Lorentz1*d_t);
beta2=(((omega_Lorentz2)^2)*d_t)/(1+gama_Lorentz2*d_t);
beta3=(((omega_Lorentz3)^2)*d_t)/(1+gama_Lorentz3*d_t);
ceta1=(eps0*delta_epsiLon1*d_t*((omega_Lorentz1)^2))/(1+gama_Lorentz1*d_t);
ceta2=(eps0*delta_epsiLon2*d_t*((omega_Lorentz2)^2))/(1+gama_Lorentz2*d_t);
ceta3=(eps0*delta_epsiLon3*d_t*((omega_Lorentz3)^2))/(1+gama_Lorentz3*d_t);
%*****************************************************************
% coificients in metal meduim
%*****************************************************************
a1=(2+(gama_L1*d_t)-(((omega_L1*d_t)^2)/2))/(1+(gama_L1*d_t)+(((omega_L1*d_t)^2)/2));
a2=(2+(gama_L2*d_t)-(((omega_L2*d_t)^2)/2))/(1+(gama_L2*d_t)+(((omega_L2*d_t)^2)/2));
b1= 1/(1+(gama_L1*d_t)+(((omega_L1*d_t)^2)/2));
b2= 1/(1+(gama_L2*d_t)+(((omega_L2*d_t)^2)/2));
c1=((eps0*delta_epsL1*(omega_L1*d_t)^2)/2)/(1+(gama_L1*d_t)+(((omega_L1*d_t)^2)/2));
c2=((eps0*delta_epsL2*(omega_L2*d_t)^2)/2)/(1+(gama_L2*d_t)+(((omega_L2*d_t)^2)/2));
d1=(2+(gama_d*d_t))/(1+(gama_d*d_t));
x1= 1/(1+(gama_d*d_t));
G1=((eps0*((omega_d*d_t)^2))/2)/(1+(gama_d*d_t));

k1=(1/(eps0*eps_inf))/(1+((c1+c2+G1)/(eps0*eps_inf)));
k2=(d1/(eps0*eps_inf))/(1+((c1+c2+G1)/(eps0*eps_inf)));
k3=(x1/(eps0*eps_inf))/(1+((c1+c2+G1)/(eps0*eps_inf)));
k4=(a1/(eps0*eps_inf))/(1+((c1+c2+G1)/(eps0*eps_inf)));
k5=(b1/(eps0*eps_inf))/(1+((c1+c2+G1)/(eps0*eps_inf)));
k6=(a2/(eps0*eps_inf))/(1+((c1+c2+G1)/(eps0*eps_inf)));
k7=(b2/(eps0*eps_inf))/(1+((c1+c2+G1)/(eps0*eps_inf)));
k8= ((c1+c2+G1)/(eps0*eps_inf))/(1+((c1+c2+G1)/(eps0*eps_inf)));

%*************************************************************
% initialization
%*************************************************************
er1= zeros(N_r-1,N_z);
er2= zeros(N_r-1,N_z);

JD_r0= zeros(N_r-1,N_z);
JD_r1= zeros(N_r-1,N_z);
JD_r2= zeros(N_r-1,N_z);

QL_r0= zeros(N_r-1,N_z);
QL_r1= zeros(N_r-1,N_z);
QL_r2= zeros(N_r-1,N_z);

et1= zeros(N_r,N_z);
et2= zeros(N_r,N_z);

JD_t0= zeros(N_r,N_z);
JD_t1= zeros(N_r,N_z);
JD_t2= zeros(N_r,N_z);

QL_t0= zeros(N_r,N_z);
QL_t1= zeros(N_r,N_z);
QL_t2= zeros(N_r,N_z);

ez1= zeros(N_r,N_z-1);
ez2= zeros(N_r,N_z-1);

JD_z0= zeros(N_r,N_z-1);
JD_z1= zeros(N_r,N_z-1);
JD_z2= zeros(N_r,N_z-1);

QL_z0= zeros(N_r,N_z-1);
QL_z1= zeros(N_r,N_z-1);
QL_z2= zeros(N_r,N_z-1);

Q_r0= zeros(N_r-1,N_z);
Q_r1= zeros(N_r-1,N_z);
Q_r2= zeros(N_r-1,N_z);

Q_t0= zeros(N_r,N_z);
Q_t1= zeros(N_r,N_z);
Q_t2= zeros(N_r,N_z);

Q_z0= zeros(N_r,N_z-1);
Q_z1= zeros(N_r,N_z-1);
Q_z2= zeros(N_r,N_z-1);

hr1= zeros(N_r,N_z-1); 
hr2= zeros(N_r,N_z-1);  

ht1= zeros(N_r-1,N_z-1);
ht2= zeros(N_r-1,N_z-1);
 
hz1= zeros(N_r-1,N_z);
hz2= zeros(N_r-1,N_z);

Pr= zeros(N_r-1,N_z);
Pr1= zeros(N_r-1,N_z);
Pr1_1= zeros(N_r-1,N_z);
Pr1_2= zeros(N_r-1,N_z);
Pr2= zeros(N_r-1,N_z);
Pr2_1= zeros(N_r-1,N_z);
Pr2_2= zeros(N_r-1,N_z);
Pr3= zeros(N_r-1,N_z);
Pr3_1= zeros(N_r-1,N_z);
Pr3_2= zeros(N_r-1,N_z);

Pt= zeros(N_r,N_z);
Pt1= zeros(N_r,N_z);
Pt1_1= zeros(N_r,N_z);
Pt1_2= zeros(N_r,N_z);
Pt2= zeros(N_r,N_z);
Pt2_1= zeros(N_r,N_z);
Pt2_2= zeros(N_r,N_z);
Pt3= zeros(N_r,N_z);
Pt3_1= zeros(N_r,N_z);
Pt3_2= zeros(N_r,N_z);

Pz= zeros(N_r,N_z-1);
Pz1= zeros(N_r,N_z-1);
Pz1_1= zeros(N_r,N_z-1);
Pz1_2= zeros(N_r,N_z-1);
Pz2= zeros(N_r,N_z-1);
Pz2_1= zeros(N_r,N_z-1);
Pz2_2= zeros(N_r,N_z-1);
Pz3= zeros(N_r,N_z-1);
Pz3_1= zeros(N_r,N_z-1);
Pz3_2= zeros(N_r,N_z-1);

r= zeros(N_r,N_z); 
r_1= zeros(N_r,N_z); 
r_2= zeros(N_r,N_z);
r_3= zeros(N_r,N_z);
r_4= zeros(N_r,N_z); 
r_5= zeros(N_r,N_z);
r_a= zeros(N_r,N_z);
r_b= zeros(N_r,N_z);
r_c= zeros(N_r,N_z);
r_d= zeros(N_r,N_z);

Jr1= zeros(N_r-1,N_z);
Jr1_1= zeros(N_r-1,N_z);
Jr1_2= zeros(N_r-1,N_z);
Jr2= zeros(N_r-1,N_z);
Jr2_1= zeros(N_r-1,N_z);
Jr2_2= zeros(N_r-1,N_z);
Jr3= zeros(N_r-1,N_z);
Jr3_1= zeros(N_r-1,N_z);
Jr3_2= zeros(N_r-1,N_z);

Jt1= zeros(N_r,N_z);
Jt1_1= zeros(N_r,N_z);
Jt1_2= zeros(N_r,N_z);
Jt2= zeros(N_r,N_z);
Jt2_1= zeros(N_r,N_z);
Jt2_2= zeros(N_r,N_z);
Jt3= zeros(N_r,N_z);
Jt3_1= zeros(N_r,N_z);
Jt3_2= zeros(N_r,N_z);

Jz1= zeros(N_r,N_z-1);
Jz1_1= zeros(N_r,N_z-1);
Jz1_2= zeros(N_r,N_z-1);
Jz2= zeros(N_r,N_z-1);
Jz2_1= zeros(N_r,N_z-1);
Jz2_2= zeros(N_r,N_z-1);
Jz3= zeros(N_r,N_z-1);
Jz3_1= zeros(N_r,N_z-1);
Jz3_2= zeros(N_r,N_z-1);

Dr1= zeros(N_r-1,N_z);
Dr2= zeros(N_r-1,N_z);

Dt1= zeros(N_r,N_z);
Dt2= zeros(N_r,N_z);

Dz1= zeros(N_r,N_z-1);
Dz2= zeros(N_r,N_z-1);

 
 F_ez2= zeros(1,length(f));
 FF_ez2= zeros(1,length(f));
 F1_ez2= zeros(1,length(f));
 F_Jz= zeros(1,length(f));
 
%**************************************************************************
% the raduis of regions
%**************************************************************************
 r_min= R_ring-(N_r*d_r)/2;

for i=1:N_r
    r(i,1:N_z)= r_min+ d_r*(i-1);
end


%     for i=1:ib 
%         ez0(i,:,:)=besselj(0,t*(N_r+i-1)*d_r)*bessely(0,t*d_r*2*N_r)-bessely(0,t*(N_r+i-1)*d_r)*besselj(0,t*2*d_r*N_r);  
%     end 


sigM1 = 0.0;  


%*************************************************************
% source excitation
%*************************************************************
n_r_source_1= n_rc;       %location of z-directed current source 
n_z_source_1= round((n_z3 + n_z4)/2);       %location of z-directed current source  

% n_r_source_2=  n_rc;
% n_z_source_2= (n_z1 + n_z2)/2;



%***********************************************************************  
%     BEGIN TIME-STEPPING LOOP  
%***********************************************************************  
   tic

for n_t=1:N_t;  
     


%********************************************************************
% source parameters
%********************************************************************
tau= 0.5/fmax_1;
t_delay= 3*tau;

% J_z(1,n_t)=sin(2*pi*f0*n_t*d_t)* exp(-((n_t*d_t-t_delay)^2/tau^2));
  J_z(1,n_t)= sin(2*pi*(c0/lambda0)*n_t*d_t)* exp(-((n_t*d_t-t_delay)^2/tau^2)); 
% if n_t<4163
%     J_z(1,n_t)=sin(2*pi*f1*n_t*d_t);
% else 
%     J_z(1,n_t)=0;
% end
%     J_z(1,n_t)= exp(-((n_t*d_t-t_delay)^2/tau^2));  

%  t= n_t*d_t;
%***********************************************************************  
%     Update magnetic fields  
%*********************************************************************** 
hr2(2:N_r-1,1:N_z-1)=hr1(2:N_r-1,1:N_z-1)+... 
                    (ct./r(2:N_r-1,1:N_z-1)).*(ez1(2:N_r-1,1:N_z-1))+cz*(et1(2:N_r-1,2:N_z)-et1(2:N_r-1,1:N_z-1)); 
 
ht2(1:N_r-1,1:N_z-1)=ht1(1:N_r-1,1:N_z-1)-...  
                   cz*(er1(1:N_r-1,2:N_z)-er1(1:N_r-1,1:N_z-1))+cr*(ez1(2:N_r,1:N_z-1)-ez1(1:N_r-1,1:N_z-1)); 

               
hz2(1:N_r-1,2:N_z-1)=hz1(1:N_r-1,2:N_z-1)-(d_t./(2*mu0.*(r(1:N_r-1,2:N_z-1)+d_r/2))).*(et1(1:N_r-1,2:N_z-1)+et1(2:N_r,2:N_z-1))-... 
                     cr*(et1(2:N_r,2:N_z-1)-et1(1:N_r-1,2:N_z-1))-(ct./(r(1:N_r-1,2:N_z-1)+d_r/2)).*(er1(1:N_r-1,2:N_z-1));                                                                
     
%***************************************************************************
%  update  first term polarization curent density in region 3 (stripe region)
%***************************************************************************                
 Jr1_2(n_r2:n_r3,n_z1:n_z2)= alpha1* Jr1_1(n_r2:n_r3,n_z1:n_z2)-beta1* Pr1_1(n_r2:n_r3,n_z1:n_z2)+...
       ceta1* er1(n_r2:n_r3,n_z1:n_z2);                 
                  
 Jt1_2(n_r2:n_r3,n_z1:n_z2)= alpha1* Jt1_1(n_r2:n_r3,n_z1:n_z2)-beta1* Pt1_1(n_r2:n_r3,n_z1:n_z2)+...
       ceta1* et1(n_r2:n_r3,n_z1:n_z2); 
   
 Jz1_2(n_r2:n_r3,n_z1:n_z2)= alpha1* Jz1_1(n_r2:n_r3,n_z1:n_z2)-beta1* Pz1_1(n_r2:n_r3,n_z1:n_z2)+...
       ceta1* ez1(n_r2:n_r3,n_z1:n_z2);            
%*********************************************************************
%  update second term  polarization curent density in region 3 (stripe region)
%*********************************************************************                
  Jr2_2(n_r2:n_r3,n_z1:n_z2)=alpha2*Jr2_1(n_r2:n_r3,n_z1:n_z2)-beta2*Pr2_1(n_r2:n_r3,n_z1:n_z2)+...
       ceta2*er1(n_r2:n_r3,n_z1:n_z2);                 
                 
 Jt2_2(n_r2:n_r3,n_z1:n_z2)=alpha2*Jt2_1(n_r2:n_r3,n_z1:n_z2)-beta2*Pt2_1(n_r2:n_r3,n_z1:n_z2)+...
      ceta2*et1(n_r2:n_r3,n_z1:n_z2); 
   
 Jz2_2(n_r2:n_r3,n_z1:n_z2)=alpha2*Jz2_1(n_r2:n_r3,n_z1:n_z2)-beta2*Pz2_1(n_r2:n_r3,n_z1:n_z2)+...
       ceta2*ez1(n_r2:n_r3,n_z1:n_z2);                                                                      
%*********************************************************************
%  update third term  polarization curent density in region 3 (stipe region)
%*********************************************************************                
  Jr3_2(n_r2:n_r3,n_z1:n_z2)=alpha3*Jr3_1(n_r2:n_r3,n_z1:n_z2)-beta3*Pr3_1(n_r2:n_r3,n_z1:n_z2)+...
       ceta3*er1(n_r2:n_r3,n_z1:n_z2);                 
                  
  Jt3_2(n_r2:n_r3,n_z1:n_z2)=alpha3*Jt3_1(n_r2:n_r3,n_z1:n_z2)-beta3*Pt3_1(n_r2:n_r3,n_z1:n_z2)+...
       ceta3*et1(n_r2:n_r3,n_z1:n_z2); 
   
  Jz3_2(n_r2:n_r3,n_z1:n_z2)=alpha3*Jz3_1(n_r2:n_r3,n_z1:n_z2)-beta3*Pz3_1(n_r2:n_r3,n_z1:n_z2)+...
       ceta3*ez1(n_r2:n_r3,n_z1:n_z2);           
 %*********************************************************************************
 % update first term Polarization vector in stripe region- region 3
 %*********************************************************************************
    Pr1_2(n_r2:n_r3,n_z1:n_z2)= Pr1_1(n_r2:n_r3,n_z1:n_z2)+d_t*Jr1_2(n_r2:n_r3,n_z1:n_z2);
    Pt1_2(n_r2:n_r3,n_z1:n_z2)= Pt1_1(n_r2:n_r3,n_z1:n_z2)+d_t*Jt1_2(n_r2:n_r3,n_z1:n_z2);  
    Pz1_2(n_r2:n_r3,n_z1:n_z2)= Pz1_1(n_r2:n_r3,n_z1:n_z2)+d_t*Jz1_2(n_r2:n_r3,n_z1:n_z2);  
 %*********************************************************************************
 % update second term Polarization vector in stripe region- region 3
 %*********************************************************************************
    Pr2_2(n_r2:n_r3,n_z1:n_z2)= Pr2_1(n_r2:n_r3,n_z1:n_z2)+d_t*Jr2_2(n_r2:n_r3,n_z1:n_z2);
    Pt2_2(n_r2:n_r3,n_z1:n_z2)= Pt2_1(n_r2:n_r3,n_z1:n_z2)+d_t*Jt2_2(n_r2:n_r3,n_z1:n_z2);  
    Pz2_2(n_r2:n_r3,n_z1:n_z2)= Pz2_1(n_r2:n_r3,n_z1:n_z2)+d_t*Jz2_2(n_r2:n_r3,n_z1:n_z2);  
 %*********************************************************************************
 % update third term Polarization vector in stripe region-region 3
 %*********************************************************************************
    Pr3_2(n_r2:n_r3,n_z1:n_z2)= Pr3_1(n_r2:n_r3,n_z1:n_z2)+d_t*Jr3_2(n_r2:n_r3,n_z1:n_z2);
    Pt3_2(n_r2:n_r3,n_z1:n_z2)= Pt3_1(n_r2:n_r3,n_z1:n_z2)+d_t*Jt3_2(n_r2:n_r3,n_z1:n_z2);  
    Pz3_2(n_r2:n_r3,n_z1:n_z2)= Pz3_1(n_r2:n_r3,n_z1:n_z2)+d_t*Jz3_2(n_r2:n_r3,n_z1:n_z2);    
 %*********************************************************************************
 % update orginal  Polarization vector in stripe region -region 3
 %*********************************************************************************
    Pr(n_r2:n_r3,n_z1:n_z2)= Pr1_2(n_r2:n_r3,n_z1:n_z2)+ Pr2_2(n_r2:n_r3,n_z1:n_z2)+ Pr3_2(n_r2:n_r3,n_z1:n_z2);
    Pt(n_r2:n_r3,n_z1:n_z2)= Pt1_2(n_r2:n_r3,n_z1:n_z2)+ Pt2_2(n_r2:n_r3,n_z1:n_z2)+ Pt3_2(n_r2:n_r3,n_z1:n_z2);  
    Pz(n_r2:n_r3,n_z1:n_z2)= Pz1_2(n_r2:n_r3,n_z1:n_z2)+ Pz2_2(n_r2:n_r3,n_z1:n_z2)+ Pz3_2(n_r2:n_r3,n_z1:n_z2);  
    
%***************************************************************************
%  update  first term polarization curent density in region 4 (core region)
%***************************************************************************                
 Jr1_2(n_r1:n_r4,n_z3:n_z4)= alpha1* Jr1_1(n_r1:n_r4,n_z3:n_z4)-beta1* Pr1_1(n_r1:n_r4,n_z3:n_z4)+...
       ceta1* er1(n_r1:n_r4,n_z3:n_z4);                 
                  
 Jt1_2(n_r1:n_r4,n_z3:n_z4)= alpha1* Jt1_1(n_r1:n_r4,n_z3:n_z4)-beta1* Pt1_1(n_r1:n_r4,n_z3:n_z4)+...
       ceta1* et1(n_r1:n_r4,n_z3:n_z4); 
   
 Jz1_2(n_r1:n_r4,n_z3:n_z4)= alpha1* Jz1_1(n_r1:n_r4,n_z3:n_z4)-beta1* Pz1_1(n_r1:n_r4,n_z3:n_z4)+...
       ceta1* ez1(n_r1:n_r4,n_z3:n_z4);            
%*********************************************************************
%  update second term  polarization curent density in region 4(core region)
%*********************************************************************                
 Jr2_2(n_r1:n_r4,n_z3:n_z4)=alpha2*Jr2_1(n_r1:n_r4,n_z3:n_z4)-beta2*Pr2_1(n_r1:n_r4,n_z3:n_z4)+...
       ceta2*er1(n_r1:n_r4,n_z3:n_z4);                 
                 
 Jt2_2(n_r1:n_r4,n_z3:n_z4)=alpha2*Jt2_1(n_r1:n_r4,n_z3:n_z4)-beta2*Pt2_1(n_r1:n_r4,n_z3:n_z4)+...
      ceta2*et1(n_r1:n_r4,n_z3:n_z4); 
   
 Jz2_2(n_r1:n_r4,n_z3:n_z4)=alpha2*Jz2_1(n_r1:n_r4,n_z3:n_z4)-beta2*Pz2_1(n_r1:n_r4,n_z3:n_z4)+...
       ceta2*ez1(n_r1:n_r4,n_z3:n_z4);                                                                      
%*********************************************************************
%  update third term  polarization curent density in region 4 (core region)
%*********************************************************************                
 Jr3_2(n_r1:n_r4,n_z3:n_z4)=alpha3*Jr3_1(n_r1:n_r4,n_z3:n_z4)-beta3*Pr3_1(n_r1:n_r4,n_z3:n_z4)+...
       ceta3*er1(n_r1:n_r4,n_z3:n_z4);                 
                  
 Jt3_2(n_r1:n_r4,n_z3:n_z4)=alpha3*Jt3_1(n_r1:n_r4,n_z3:n_z4)-beta3*Pt3_1(n_r1:n_r4,n_z3:n_z4)+...
       ceta3*et1(n_r1:n_r4,n_z3:n_z4); 
   
 Jz3_2(n_r1:n_r4,n_z3:n_z4)=alpha3*Jz3_1(n_r1:n_r4,n_z3:n_z4)-beta3*Pz3_1(n_r1:n_r4,n_z3:n_z4)+...
       ceta3*ez1(n_r1:n_r4,n_z3:n_z4);           
 %*********************************************************************************
 % update first term Polarization vector in core region- region8
 %*********************************************************************************
    Pr1_2(n_r1:n_r4,n_z3:n_z4)= Pr1_1(n_r1:n_r4,n_z3:n_z4)+d_t*Jr1_2(n_r1:n_r4,n_z3:n_z4);
    Pt1_2(n_r1:n_r4,n_z3:n_z4)= Pt1_1(n_r1:n_r4,n_z3:n_z4)+d_t*Jt1_2(n_r1:n_r4,n_z3:n_z4);  
    Pz1_2(n_r1:n_r4,n_z3:n_z4)= Pz1_1(n_r1:n_r4,n_z3:n_z4)+d_t*Jz1_2(n_r1:n_r4,n_z3:n_z4);  
 %*********************************************************************************
 % update second term Polarization vector in core region- region 8
 %*********************************************************************************
    Pr2_2(n_r1:n_r4,n_z3:n_z4)= Pr2_1(n_r1:n_r4,n_z3:n_z4)+d_t*Jr2_2(n_r1:n_r4,n_z3:n_z4);
    Pt2_2(n_r1:n_r4,n_z3:n_z4)= Pt2_1(n_r1:n_r4,n_z3:n_z4)+d_t*Jt2_2(n_r1:n_r4,n_z3:n_z4);  
    Pz2_2(n_r1:n_r4,n_z3:n_z4)= Pz2_1(n_r1:n_r4,n_z3:n_z4)+d_t*Jz2_2(n_r1:n_r4,n_z3:n_z4);  
 %*********************************************************************************
 % update third term Polarization vector in core region-region 8
 %*********************************************************************************
    Pr3_2(n_r1:n_r4,n_z3:n_z4)= Pr3_1(n_r1:n_r4,n_z3:n_z4)+d_t*Jr3_2(n_r1:n_r4,n_z3:n_z4);
    Pt3_2(n_r1:n_r4,n_z3:n_z4)= Pt3_1(n_r1:n_r4,n_z3:n_z4)+d_t*Jt3_2(n_r1:n_r4,n_z3:n_z4);  
    Pz3_2(n_r1:n_r4,n_z3:n_z4)= Pz3_1(n_r1:n_r4,n_z3:n_z4)+d_t*Jz3_2(n_r1:n_r4,n_z3:n_z4);    
 %*********************************************************************************
 % update orginal  Polarization vector in core region -region 8
 %*********************************************************************************
    Pr(n_r1:n_r4,n_z3:n_z4)= Pr1_2(n_r1:n_r4,n_z3:n_z4)+ Pr2_2(n_r1:n_r4,n_z3:n_z4)+ Pr3_2(n_r1:n_r4,n_z3:n_z4);
    Pt(n_r1:n_r4,n_z3:n_z4)= Pt1_2(n_r1:n_r4,n_z3:n_z4)+ Pt2_2(n_r1:n_r4,n_z3:n_z4)+ Pt3_2(n_r1:n_r4,n_z3:n_z4);  
    Pz(n_r1:n_r4,n_z3:n_z4)= Pz1_2(n_r1:n_r4,n_z3:n_z4)+ Pz2_2(n_r1:n_r4,n_z3:n_z4)+ Pz3_2(n_r1:n_r4,n_z3:n_z4);      
    
    
%**********************************************************************
% Update electric flux density  fields in region 1(metal) 
%*********************************************************************** 
Dr2(1:N_r-1,2:n_z1-1)= Dr1(1:N_r-1,2:n_z1-1)-...
     cb*(ht2(1:N_r-1,2:n_z1-1)-ht2(1:N_r-1,1:n_z1-2))+... 
                  (cc./(r(1:N_r-1,2:n_z1-1)+d_r/2)).*(hz2(1:N_r-1,2:n_z1-1)); 
                                
Dt2(2:N_r-1,2:n_z1-1)= Dt1(2:N_r-1,2:n_z1-1)+...
    cb*(hr2(2:N_r-1,2:n_z1-1)-hr2(2:N_r-1,1:n_z1-2))- ca.*(hz2(2:N_r-1,2:n_z1-1)- hz2(1:N_r-2,2:n_z1-1)); 
 
Dz2(2:N_r-1,1:n_z1-1)= Dz1(2:N_r-1,1:n_z1-1)+...  
                   ca*(ht2(2:N_r-1,1:n_z1-1)-ht2(1:N_r-2,1:n_z1-1))+(ck./r(2:N_r-1,1:n_z1-1)).*(ht2(2:N_r-1,1:n_z1-1)+ht2(1:N_r-2,1:n_z1-1))-...  
                   (cc./r(2:N_r-1,1:n_z1-1)).*(hr2(2:N_r-1,1:n_z1-1));      
               
%***********************************************************************  
%     Update electric flux density fields in region 2 ( silica meduim)
%***********************************************************************
 Dr2(1:N_r-1,n_z1:N_z-1)= Dr1(1:N_r-1,n_z1:N_z-1)-...
     cb*(ht2(1:N_r-1,n_z1:N_z-1)-ht2(1:N_r-1,n_z1-1:N_z-2))+... 
                  (cc./(r(1:N_r-1,n_z1:N_z-1)+d_r/2)).*(hz2(1:N_r-1,n_z1:N_z-1)); 
                                
Dt2(2:N_r-1,n_z1:N_z-1)= Dt1(2:N_r-1,n_z1:N_z-1)+...
    cb*(hr2(2:N_r-1,n_z1:N_z-1)-hr2(2:N_r-1,n_z1-1:N_z-2))-ca*(hz2(2:N_r-1,n_z1:N_z-1)-hz2(1:N_r-2,n_z1:N_z-1)); 
 
Dz2(2:N_r-1,n_z1:N_z-1)= Dz1(2:N_r-1,n_z1:N_z-1)+...  
                   ca*(ht2(2:N_r-1,n_z1:N_z-1)-ht2(1:N_r-2,n_z1:N_z-1))+(ck*(ht2(2:N_r-1,n_z1:N_z-1)+ht2(1:N_r-2,n_z1:N_z-1)))./r(2:N_r-1,n_z1:N_z-1)-...  
                   (cc*(hr2(2:N_r-1,n_z1:N_z-1)))./r(2:N_r-1,n_z1:N_z-1);
                
 %*****************************************************************************
   % update electric displacement vector in region 3 -(stripe region)
 %******************************************************************************       
   Dr2(n_r2:n_r3,n_z1:n_z2)= Dr1(n_r2:n_r3,n_z1:n_z2)-...
     cb*(ht2(n_r2:n_r3,n_z1:n_z2)-ht2(n_r2:n_r3,n_z1-1:n_z2-1))+... 
                  (cc./(r(n_r2:n_r3,n_z1:n_z2)+d_r/2)).*(hz2(n_r2:n_r3,n_z1:n_z2)); 
                                
   Dt2(n_r2:n_r3,n_z1:n_z2)= Dt1(n_r2:n_r3,n_z1:n_z2)+...
    cb*(hr2(n_r2:n_r3,n_z1:n_z2)-hr2(n_r2:n_r3,n_z1-1:n_z2-1))- ca .*(hz2(n_r2:n_r3,n_z1:n_z2)- hz2(n_r2-1:n_r3-1,n_z1:n_z2)); 
 
   Dz2(n_r2:n_r3,n_z1:n_z2)= Dz1(n_r2:n_r3,n_z1:n_z2)+...  
                   ca*(ht2(n_r2:n_r3,n_z1:n_z2)-ht2(n_r2-1:n_r3-1,n_z1:n_z2))+(ck*(ht2(n_r2:n_r3,n_z1:n_z2)+ht2(n_r2-1:n_r3-1,n_z1:n_z2)))./r(n_r2:n_r3,n_z1:n_z2)-...  
                   (cc./r(n_r2:n_r3,n_z1:n_z2)).*(hr2(n_r2:n_r3,n_z1:n_z2));   
          
%*****************************************************************************
   % update electric displacement vector in region 4 -silicon-(core region)
 %******************************************************************************       
   Dr2(n_r1:n_r4,n_z3:n_z4)= Dr1(n_r1:n_r4,n_z3:n_z4)-...
     cb*(ht2(n_r1:n_r4,n_z3:n_z4)-ht2(n_r1:n_r4,n_z3-1:n_z4-1))+... 
                  (cc./(r(n_r1:n_r4,n_z3:n_z4)+d_r/2)).*(hz2(n_r1:n_r4,n_z3:n_z4)); 
                                
   Dt2(n_r1:n_r4,n_z3:n_z4)= Dt1(n_r1:n_r4,n_z3:n_z4)+...
    cb*(hr2(n_r1:n_r4,n_z3:n_z4)-hr2(n_r1:n_r4,n_z3-1:n_z4-1))- ca .*(hz2(n_r1:n_r4,n_z3:n_z4)- hz2(n_r1-1:n_r4-1,n_z3:n_z4)); 
 
   Dz2(n_r1:n_r4,n_z3:n_z4)= Dz1(n_r1:n_r4,n_z3:n_z4)+...  
                   ca*(ht2(n_r1:n_r4,n_z3:n_z4)-ht2(n_r1-1:n_r4-1,n_z3:n_z4))+(ck*(ht2(n_r1:n_r4,n_z3:n_z4)+ht2(n_r1-1:n_r4-1,n_z3:n_z4)))./r(n_r1:n_r4,n_z3:n_z4)-...  
                   (cc./r(n_r1:n_r4,n_z3:n_z4)).*(hr2(n_r1:n_r4,n_z3:n_z4));    
               

%***********************************************************************  
%     Update electric fields in region 1 (metal) 
%***********************************************************************  
  er2(1:N_r-1,2:n_z1-1)= k1* Dr2(1:N_r-1,2:n_z1-1)-k2* JD_r1(1:N_r-1,2:n_z1-1)+k3* JD_r0(1:N_r-1,2:n_z1-1)-...
      k4* QL_r1(1:N_r-1,2:n_z1-1)+k5* QL_r0(1:N_r-1,2:n_z1-1)-k6* Q_r1(1:N_r-1,2:n_z1-1)+...
      k7* Q_r0(1:N_r-1,2:n_z1-1)-k8* er1(1:N_r-1,2:n_z1-1);
  
  et2(2:N_r-1,2:n_z1-1)= k1* Dt2(2:N_r-1,2:n_z1-1)-k2* JD_t1(2:N_r-1,2:n_z1-1)+k3* JD_t0(2:N_r-1,2:n_z1-1)-...
      k4* QL_t1(2:N_r-1,2:n_z1-1)+k5* QL_t0(2:N_r-1,2:n_z1-1)-k6* Q_t1(2:N_r-1,2:n_z1-1)+...
      k7* Q_t0(2:N_r-1,2:n_z1-1)-k8* et1(2:N_r-1,2:n_z1-1);
  
  ez2(2:N_r-1,1:n_z1-1)= k1* Dz2(2:N_r-1,1:n_z1-1)-k2* JD_z1(2:N_r-1,1:n_z1-1)+k3* JD_z0(2:N_r-1,1:n_z1-1)-...
      k4* QL_z1(2:N_r-1,1:n_z1-1)+k5* QL_z0(2:N_r-1,1:n_z1-1)-k6* Q_z1(2:N_r-1,1:n_z1-1)+...
      k7* Q_z0(2:N_r-1,1:n_z1-1)-k8* ez1(2:N_r-1,1:n_z1-1); 
                 
%***********************************************************************  
%     Update electric fields in region 1 (silica-meduim)
%***********************************************************************                
   er2(1:N_r-1,n_z1:N_z-1)=(1./(eps0.*epsilonr(1:N_r-1,n_z1:N_z-1))).*Dr2(1:N_r-1,n_z1:N_z-1);
   et2(2:N_r-1,n_z1:N_z-1)=(1./(eps0.*epsilonr(2:N_r-1,n_z1:N_z-1))).*Dt2(2:N_r-1,n_z1:N_z-1);
   ez2(2:N_r-1,n_z1:N_z-1)=(1./(eps0.*epsilonr(2:N_r-1,n_z1:N_z-1))).*Dz2(2:N_r-1,n_z1:N_z-1);
   
                        
%***********************************************************************  
%     Update electric fields in region 3 (stripe- silicon) 
%***********************************************************************               
  er2(n_r2:n_r3,n_z1:n_z2)=(1/(eps0*eps_infinite))*(Dr2(n_r2:n_r3,n_z1:n_z2)-Pr(n_r2:n_r3,n_z1:n_z2));
  et2(n_r2:n_r3,n_z1:n_z2)=(1/(eps0*eps_infinite))*(Dt2(n_r2:n_r3,n_z1:n_z2)-Pt(n_r2:n_r3,n_z1:n_z2));
  ez2(n_r2:n_r3,n_z1:n_z2)=(1/(eps0*eps_infinite))*(Dz2(n_r2:n_r3,n_z1:n_z2)-Pz(n_r2:n_r3,n_z1:n_z2));
   
 %*********************************************************************** 
%     Update electric fields in region 4 (core- silicon) 
%***********************************************************************               
  er2(n_r1:n_r4,n_z3:n_z4)=(1/(eps0*eps_infinite))*(Dr2(n_r1:n_r4,n_z3:n_z4)-Pr(n_r1:n_r4,n_z3:n_z4));
  et2(n_r1:n_r4,n_z3:n_z4)=(1/(eps0*eps_infinite))*(Dt2(n_r1:n_r4,n_z3:n_z4)-Pt(n_r1:n_r4,n_z3:n_z4));
  ez2(n_r1:n_r4,n_z3:n_z4)=(1/(eps0*eps_infinite))*(Dz2(n_r1:n_r4,n_z3:n_z4)-Pz(n_r1:n_r4,n_z3:n_z4));
     
%***********************************************************
% update JD in region 1(metal)
%************************************************************
JD_r2(1:N_r-1,2:n_z1-1)= d1* JD_r1(1:N_r-1,2:n_z1-1)-x1* JD_r0(1:N_r-1,2:n_z1-1)+...
    G1*(er2(1:N_r-1,2:n_z1-1)+er1(1:N_r-1,2:n_z1-1));

JD_t2(2:N_r-1,2:n_z1-1)= d1* JD_t1(2:N_r-1,2:n_z1-1)-x1* JD_t0(2:N_r-1,2:n_z1-1)+...
    G1*(et2(2:N_r-1,2:n_z1-1)+et1(2:N_r-1,2:n_z1-1));

JD_z2(2:N_r-1,1:n_z1-1)= d1* JD_z1(2:N_r-1,1:n_z1-1)-x1* JD_z0(2:N_r-1,1:n_z1-1)+...
    G1*(ez2(2:N_r-1,1:n_z1-1)+ez1(2:N_r-1,1:n_z1-1));

%*******************************************************************
% update QL( first term of lorentz model) in  region 1(metal)
%******************************************************************* 
QL_r2(1:N_r-1,2:n_z1-1)= a1* QL_r1(1:N_r-1,2:n_z1-1)-b1* QL_r0(1:N_r-1,2:n_z1-1)+...
    c1*(er2(1:N_r-1,2:n_z1-1)+er1(1:N_r-1,2:n_z1-1));

QL_t2(2:N_r-1,2:n_z1-1)= a1* QL_t1(2:N_r-1,2:n_z1-1)-b1* QL_t0(2:N_r-1,2:n_z1-1)+...
    c1*(et2(2:N_r-1,2:n_z1-1)+et1(2:N_r-1,2:n_z1-1));

QL_z2(2:N_r-1,1:n_z1-1)= a1* QL_z1(2:N_r-1,1:n_z1-1)-b1* QL_z0(2:N_r-1,1:n_z1-1)+...
    c1*(ez2(2:N_r-1,1:n_z1-1)+ez1(2:N_r-1,1:n_z1-1));

%*******************************************************************
% update Q (second term of lorentz model) in  region 1(metal)
%*******************************************************************
Q_r2(1:N_r-1,2:n_z1-1)= a2*Q_r1(1:N_r-1,2:n_z1-1)-b2* Q_r0(1:N_r-1,2:n_z1-1)+...
    c2*(er2(1:N_r-1,2:n_z1-1)+er1(1:N_r-1,2:n_z1-1));

Q_t2(2:N_r-1,2:n_z1-1)= a2*Q_t1(2:N_r-1,2:n_z1-1)-b2* Q_t0(2:N_r-1,2:n_z1-1)+...
    c2*(et2(2:N_r-1,2:n_z1-1)+et1(2:N_r-1,2:n_z1-1));

Q_z2(2:N_r-1,1:n_z1-1)= a2*Q_z1(2:N_r-1,1:n_z1-1)-b2* Q_z0(2:N_r-1,1:n_z1-1)+...
    c2*(ez2(2:N_r-1,1:n_z1-1)+ez1(2:N_r-1,1:n_z1-1));

%*****************************************************************************
% source: modulated sinusoidal gaussian pulse
%**************************************************************************
  ez2(n_r_source_1,n_z_source_1)= ez2(n_r_source_1,n_z_source_1)+ J_z(1,n_t);    
%   ez2(n_r_source_2,n_z_source_2)= ez2(n_r_source_2,n_z_source_2)+ J_z(1,n_t);

%*********************************************************
% update martixes
%*********************************************************

 et1= et2;
 er1= er2;
 ez1= ez2;

 hz1= hz2;
 hr1= hr2;
 ht1= ht2;
 
QL_r0= QL_r1;
QL_r1= QL_r2;
QL_t0= QL_t1;
QL_t1= QL_t2;
QL_z0= QL_z1;
QL_z1= QL_z2;

Q_r0= Q_r1;
Q_r1= Q_r2;
Q_t0= Q_t1;
Q_t1= Q_t2;
Q_z0= Q_z1;
Q_z1= Q_z2;

Dr1= Dr2;
Dt1= Dt2;
Dz1= Dz2;

JD_r0= JD_r1;
JD_r1= JD_r2;
JD_t0= JD_t1;
JD_t1= JD_t2;
JD_z0= JD_z1;
JD_z1= JD_z2;

Jr1_1=Jr1_2;
Jt1_1=Jt1_2;
Jz1_1=Jz1_2;

Jr2_1=Jr2_2;
Jt2_1=Jt2_2;
Jz2_1=Jz2_2;


Jr3_1=Jr3_2;
Jt3_1=Jt3_2;
Jz3_1=Jz3_2;

Pr1_1=Pr1_2;
Pt1_1=Pt1_2;
Pz1_1=Pz1_2;

Pr2_1=Pr2_2;
Pt2_1=Pt2_2;
Pz2_1=Pz2_2;

Pr3_1=Pr3_2;
Pt3_1=Pt3_2;
Pz3_1=Pz3_2;


%**************************
% time signals
%**************************
sum(1,n_t)=(ez2(n_r_point_1,n_z_point_1));
sum1(1,n_t)= ez2(n_r_point_2,n_z_point_2);
sum2(1,n_t)= ez2(n_r_point_3,n_z_point_3);   

       
end
toc
%*******************************
% fourier transform
%*******************************
% %%
%      for ii=1:length(f)
%        F_ez2(1,ii)= F_ez2(1,ii)+ exp(-1j*2*pi*f(ii)*n_t*d_t)* ez2(n_r_point_1,n_z_point_1)*d_t;  
% %        FF_ez2(1,ii)= FF_ez2(1,ii)+ exp(-1j*2*pi*f(ii)*n_t*d_t)* ez2(n_r_point_2,n_z_point_2)*d_t;
% %        F1_ez2(1,ii)= F1_ez2(1,ii)+ exp(-1j*2*pi*f(ii)*n_t*d_t)* ez2(n_r_point_3,n_z_point_3)*d_t;
% %        F_Jz(1,ii)= F_Jz(1,ii)+ exp(-1j*2*pi*f(ii)*n_t*d_t)* J_z(1,n_t)*d_t;
%       
%      end
% 
% 
% %  fft
% end 
%  toc
%  
% %*********************
% % results
% %*********************
% % source signal
% figure(1) 
% plot(J_z,'r--')
% legend('Ez(t)')
% xlabel('Nt');
% ylabel('Jz(t)');
% legend('source')
% title('signal time Jz(t)');
% % hold on
% 
% % fourier transform of source
% figure(2)
% plot(f/1e12,(abs(F_Jz).^2));
% xlabel('f(THz)');
% ylabel('fourier transform of source');
% %*********************
% %time signal in core region
% 
% figure(3)
% t=0:d_t:(N_t-1)*d_t;
% plot(sum,'b')
% xlabel('N_{t}');
% ylabel('e_{z}(t)');
% title('signal time ez(t)');
% hold on  
% 
% figure(4)
% % fourier transform signal in core region
% plot(f/1e12,real(abs(((F_ez2)).^2)));
% xlabel('Frequency(THz)');
% ylabel('Fourier transform');
% legend('M=20')
% title('resonance frequency of ring resonator')
% 
% % normalized fourier transform in core region
% figure(5)
% plot(f/1e12,(1/12*1000*((abs(F_ez2)).^2)./((abs(F_Jz)).^2)),'b');
% xlabel('Frequency(THz)');
% ylabel('Fourier transform');
% legend('m=20')
% title('resonance frequency of ring resonator')
% 
% 
% 
% 
% plot(real(ez2(:,((n_z1+n_z1)/2)+0.5)));
% xlabel('x')
% ylabel('Energy density')
% %********************
% figure(7)
% t=0:d_t:(N_t-1)*d_t;
% plot(t*1e15,sum1,'k','linewidth',1);
% xlabel('t (Femtosecond)');
% ylabel('e_{z}(t)');
% legend('m=20')
% 
% plot(f/1e12,(1/12*1000*((abs(FF_ez2)).^2)./((abs(F_Jz)).^2)),'b');
% xlabel('f(THz)');
% ylabel('|F(ez)|^2/|F(Jz)|^2');
% 
% 
% 
% 
% %***********************
% figure(9)
% plot(sum2,'b')
% xlabel('nt');
% ylabel('Ez(t)');
% legend('Ez(t)')
% 
% plot(f/1e12,(1/40*10e6*real(((abs(F1_ez2)).^2)./((abs(F_Jz)).^2))),'b');
% xlabel('f(THz)');
% ylabel('|F(ez)|^2/|F(Jz)|^2');
% legend('m=20')
% 
% %%
% 
% % [freq,growth,amp,phase]=FDM(sum,Fs,fmin,fmax,eps);
% % %%
% % [freq1,growth1,amp1,phase1]=FDM2(sum1,Fs,fmin,fmax,eps);
% % 
% % [freq2,growth2,amp2,phase2]=FDM3(sum2,Fs,fmin,fmax,eps);
% 
% %%
% x=linspace(-1,1,400);
% y=linspace(-1,1,329); 
% [X,Y]=meshgrid(x,y);
% surf(x,y,real(ez2.'));shading interp;
% colorbar 
% view(0,90)