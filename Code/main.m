%Or Bahari

function main

clear;
clc;
global Bragg n_eff NG1 NG2;
Bragg=1550e-9;  % Bragg wavelength                                     
span=40e-9;  % Set the wavelength span for the simultion
resolution=0.0125e-9;  % Set the wavelength resolution
N=span/resolution;
disp([ 'Number of points: ' num2str(N) ])
Lambda=zeros(N+1,1);
R=zeros(N+1,1);
T=zeros(N+1,1);

        for i=1:N+1
            wavelength=Bragg+(i-1-N/2)*resolution;    % Wavelength sweep
            Lambda(i)=wavelength*1e9;% in nm
            Grating_Parameters(Lambda(i))
            [r,t]=Grating_RT(wavelength); % Calculate the R and T  
            R(i)=r;
            T(i)=t;
        end
        figure(1);
        T_db=10*log10(T);
        R_db=10*log10(R);
        %subplot(2,5,j);
        
        plot(Lambda,[R_db T_db],'LineWidth',2);
        set(gca,'FontSize',14);
        xlabel('Wavelength (nm)','FontSize',14);
        ylabel('Response','FontSize',14);
        title(['cavity with NG1=',num2str(NG1),' NG2=',num2str(NG2)]);
        legend('Reflection','Transmission');
        grid on;
        [T_peak,Q]=FP_transmission(Lambda*1e-9,R);
        
        insertion_Loss=10*log(T_peak);
   
a=3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Grating_Parameters(Lambda)
%Set the parameters
global Period NG1 NG2 delta_n n1 n2 n_eff loss C_length;

%arm_length=2.5e-6;
Period=318e-9;  % Bragg period
n_eff=2.44 - 1.13 * ((Lambda/1000.0)-1.550)-0.043*((Lambda/1000.0)-1.550)^2;
NG1=100;    % Number of grating periods #1
NG2=210;    % Number of grating periods #2


L1=NG1*Period;    % Grating length #1
L2=NG2*Period;    % length #2
C_length=L1+L2; % total cavity length

delta_n=0.05;    % Index contrast between n1 and n2, matches to dw=40nm
n1=n_eff-delta_n/2;
n2=n_eff+delta_n/2;

loss=590; % (what was said in class we should use...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [R,T]=Grating_RT(wavelength)
%Calculate the R and T for a certain wavelength

M=Grating_Matrix(wavelength);

T=abs(1/M(1,1))^2;
R=abs(M(2,1)/M(1,1))^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Rpeak]=Center_peak_extraction(R,width)
Length=length(R);
Rpeak=min(R(Length*0.5:width+Length*0.5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TP,Q]= FP_transmission(wavelengths,R)
    %wavelengths Units: m
    global  n_eff C_length loss
    delta=-4*pi* 1.94*1e-6*n_eff./wavelengths;
    % 1.94*1e-6
    
    %finding the peak of R:
    Rpeak=Center_peak_extraction(R,120);
    % modify the following line to output the result of your calculation, for the auto-grader
    output = ((1-Rpeak)^2)./(((1-Rpeak)^2)+4*Rpeak*(sin(delta/2)).^2);
    output10=10*log10(output);
    figure()
    %subplot(1,2,1);
    plot(wavelengths,output10);
    xlabel('wavelengths [m]');
    ylabel('It/Ii [db]');
    title('power transmission of Fabry-Perot');
    TP=max(output);
    %subplot(1,2,2);
    bw=powerbw(output,wavelengths); %Provides Bandwidth 3db
    Q=TP/bw;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T]=Grating_Matrix(wavelength)
% Calculate the total transfer matrix of the gratings

global Period NG1 NG2;
global n1 n2 n_eff loss arm_length;

l=Period/2;

T_hw1=HomoWG_Matrix(wavelength,l,n1,loss);
T_is12=IndexStep_Matrix(n1,n2);
T_hw2=HomoWG_Matrix(wavelength,l,n2,loss);
T_is21=IndexStep_Matrix(n2,n1);
T1=T_hw1*T_is12*T_hw2*T_is21;

%code for back-to-back BGs
T_is22=IndexStep_Matrix(n2,n2);
T2=T_hw1*T_is12*T_hw2*T_is22*T_hw2*T_is21;

%code to silicon separtor between 2 back-to-back BGs
%T_is2n= IndexStep_Matrix(n2,n_eff);
%T_hwn=HomoWG_Matrix(wavelength,arm_length,n_eff,loss);
%T_isn1=IndexStep_Matrix(n_eff,n1);
%T2=T_hw1*T_is12*T_hw2*T_is2n*T_hwn*T_isn1;
%T=T1^(NG1)*T2*T1^(NG2);
T=T1^(NG1)*T_hw1*T1^(NG2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T_hw=HomoWG_Matrix(wavelength,l,neff,loss)
% Calculate the transfer matrix of a homogeneous waveguide.
%Complex propagation constant
beta=2*pi.*neff/wavelength-1i*loss/2;

v=[exp(1i*beta*l) exp(-1i*beta*l)];
T_hw=diag(v);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T_is=IndexStep_Matrix(n1,n2)
% Calculate the transfer matrix for a index step from n1 to n2.

a=(n1+n2)/(2*sqrt(n1*n2));
b=(n1-n2)/(2*sqrt(n1*n2));
T_is=[a b; b a];
        