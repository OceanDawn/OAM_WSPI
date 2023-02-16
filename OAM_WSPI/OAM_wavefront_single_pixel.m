clear;
close all;
clc;
addpath(['TVAL3' filesep 'Utilities']) % TVAL3
addpath(['TVAL3' filesep 'Solver']) % TVAL3
addpath(['TVAL3' filesep 'Fast_Walsh_Hadamard_Transform']) % TVAL3
addpath('pattern') % TVAL3

method = input('Please select a reconstruction algorithm:(0)SOC (1)TVAL3' );
%% Set projection and imaging resolution
P_pixels=768; % Projection resolution
Tpixels=128;
pixels=64; % Imaging resolution
%% The generation of the target of complex field
I_abs = imread('heart.png');
I_abs=imresize(I_abs,[P_pixels,P_pixels],'nearest');
I_abs = double(I_abs(:,:,1));
I_abs = I_abs / max(max(I_abs));
% I_abs=ones(P_pixels,P_pixels);
% I_abs = LBMask(P_pixels,P_pixels,1,'circle',385,385);
% 
phase = imread('house.png');
phase=imresize(phase,[P_pixels,P_pixels],'nearest');
phase = double(phase(:,:,1)); 
phase = phase/max(max(phase));
phase = (2*pi*phase)*1i;

XX = I_abs.*exp(phase);
Amp=abs(XX);
phase=W_angle(XX);

figure(1);
subplot(121);
imshow(Amp,[]);
subplot(122);
imshow(phase,[]);

%%  MASK  (checkerboard)
ref = zeros(Tpixels,Tpixels);
for i =1:Tpixels
    for j = 1:Tpixels
        if(mod(floor((i-1)/1)+floor((j-1)/1),2)==0)       
            ref(i,j)=1;
        else
            ref(i,j)=0;
        end
    end
end

sig = ~ref;
%% Patterns and sample
percent = input('Sampling rate：');
M=round(pixels^2*percent);
load('R1.mat');
load('R2.mat');
load('R3.mat');
load('R4.mat');
load('R5.mat');
load('R6.mat');
load('R7.mat');
load('R8.mat');
R=[R1;R2;R3;R4;R5;R6;R7;R8];
R=R(1:M,:);

% P=15;L=65;    %P=34:64   9:32
% Poam=zeros(P*L,pixels^2);  
% v=1;
% for p=0:1:P-1
%      for m=-((L-1)/2):1:((L-1)/2)       
% %         v=L*p+m+(L-1)/2+1;   % (L-1)/2+1
%         fprintf('The %d pattern is generteting \n',v);
%         p1 = LGpatterns(p,m,0.16e-3,4.096e-3,pixels);
%         Poam(v,:)=reshape(p1,1,pixels^2); 
%         v=v+1;
%       end
% end
% R=Poam(1:M,:);

%%  ----------------------------------------LG----------------------------------------
D=zeros(1,4);
D_N=zeros(M,1);
One=ones(2,2);
if method == 0
    I_R = zeros(pixels,pixels);
end
         
for v=1:1:M
    fprintf('The %d pattern is working \n',v);
    
        p1=reshape(R(v,:),pixels,pixels);

        for x=1:pixels
            for y=1:pixels
            
                Patt(x*2-1:x*2,y*2-1:y*2)=One*p1(x,y);
                
            end
        end
        
        for i=0:3
            
            pat=Patt.*sig+ref.*exp(1i*pi/2*i);  % phase shifting
            pat=imresize(pat,[P_pixels,P_pixels],'nearest');

            %sample----------------
            FT_sample=pat.*XX;
            FT_sample = fftshift(fft2(fftshift(FT_sample)));
             
            D(1,i+1) = (abs(FT_sample(385,385)))^2;
         
        end
         
         D_N(v,1)=((D(1,1)-D(1,3))+1i*(D(1,2)-D(1,4)));
         
end

%%  reconstruction
tic;
if method==0
    I_R =zeros(pixels);
    for j = 1:1:M
        p1 = reshape(R(j,:),[pixels pixels]);
        I_R=I_R+(D_N(j,1).*p1)/(pixels*pixels);
    end
    
elseif method==1
    clear opts
    opts.mu = 2^8;%（2^4）~（2^13）
    opts.beta = 2^5;%（2^4）~（2^13）
    opts.mu0 = 2^4;      % trigger continuation shceme
    opts.beta0 = 2^-2;    % trigger continuation scheme
    opts.tol = 1.e-6; %1.e-2    determine the solution accuracy. Their smaller values result in a longer elapsed time and usually a better solution quality.
    opts.tol_inn = 1.e-3;   %determine the solution accuracy. Their smaller values result in a longer elapsed time and usually a better solution quality.
    opts.maxit = 300;  
    opts.maxcnt = 10;  
    opts.TVnorm = 1;     
    % opts.nonneg = true; 
    opts.TVL2 = true;   %switch for TV/L2 models
    % opts.isreal = false;    % switch for real signals/images 
    opts.init = 1;  %1: A*b,,  0: 0,, U0: U0

    [I_R, out1] = TVAL3(R,D_N,pixels,pixels,opts);   
end
toc;

% x=I_R(:);
% xt=XX(:);
% alpha = (x'*xt)/(x'*x);
% x = alpha * x;
% I_R=reshape(x,pixels,pixels);

%% imshow
% I_R=conj(I_R);
if method == 0
    I_R_Z_abs = flipud(abs(I_R));
    I_R_Z_ang = flipud(angle(I_R));
else
    I_R_Z_abs = abs(I_R);
    I_R_Z_ang = angle(I_R);
end

figure(2);
subplot(121);
imshow(I_R_Z_abs,[]);
subplot(122);
imshow(I_R_Z_ang,[]);