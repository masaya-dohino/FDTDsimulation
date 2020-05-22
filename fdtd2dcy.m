%%%%%_______________2DTM simulation of a plane wave source impinging on a dielectric cylinder __________%%%%%%%
%_________________________SOUGATA_CHATTERJEE________________________________
%%%%%%%%%_________________SAMEER_KOLKATA_CENTER___________________%%%%%%%%%
ie = 50;
je = 50;
nfreqs = 3;
ez_inc_low_m1 = 0.;
ez_inc_low_m2 = 0.;
ez_inc_high_m1 = 0.;
ez_inc_high_m2 = 0.;
ic=fix(ie/2);
jc=fix(je/2);
ia = 8;
ib = ie - ia;
ja = 7;
jb = je - ja;
epsz=8.8e-12;
ddx=0.1;
dt=ddx/6e8;
%%%%%%-----Time_interval-----------%%%%%%%%%%%%
nsteps=100;
%%%%%%-----Time_interval-----------%%%%%%%%%%%%
freq = zeros(nfreqs);
arg = zeros(nfreqs);
real_pt = zeros(nfreqs,ie,je);
imag_pt = zeros(nfreqs,ie,je);
amp = zeros(ie,je);
phase = zeros(ie,je);
real_in = zeros(5);
imag_in = zeros(5);
amp_in = zeros(5);
phase_in = zeros(5);
ez_inc = zeros(je);
hx_inc = zeros(je);
ez=zeros(ie,je);
dz=zeros(ie,je);
hy=zeros(ie,je);
ihx=zeros(ie,je);
ihy=zeros(ie,je);
hx=zeros(ie,je);
ga=ones(ie,je);
gb=ones(ie,je);

%Parameters for the Fourier Transform
freq(1) = 50.e6;
freq(2) = 300.e6;
freq(3) = 700.e6;

for n = 1:nfreqs
    arg(n) = 2*pi*freq(n)*dt;
end
%%%%%%-----wave_specification-----------%%%%%%%%%%%%
t0 = 20.0;
spread = 8.0;
%------------calculated the PML parameter---------------%

    gi2=ones(ie);
    gi3=ones(ie);
    fi1=zeros(ie);
    fi2=ones(ie);
    fi3=ones(ie);

    gj2=ones(ie);
    gj3=ones(ie);
    fj1=zeros(ie);
    fj2=ones(ie);
    fj3=ones(ie);
    
npml=20;
for i=1:npml
    xnum=npml-i;
    xd=npml;
    xxn=xnum/xd;
    xn=0.33*((xxn)^3);
    gj2(i)=1.0/(1.0+xn);
    gi2(ie-1-i)=1.0/(1.0+xn);
    gi3(i)=(1.0-xn)/(1.0+xn);
    gi3(ie-i-1)=(1.0-xn)/(1.0+xn);
    xxn=(xnum-0.5)/xd;
    xn=0.25*((xxn)^3);
    fi1(i)=xn;
    fi1(ie-2-i)=xn;
    fi2(i)=1.0/(1.0+xn);
    fi2(ie-2-i)=1.0/(1.0+xn);
    fi3(i)=(1.0-xn)/(1.0+xn);
    fi3(ie-2-i)=(1.0-xn)/(1.0+xn);
end
for j=1:npml
    xnum=npml-j;
    xd=npml;
    xxn=xnum/xd;
    xn=0.33*((xxn)^3);
    gj2(j)=1.0/(1.0+xn);
    gi2(je-1-j)=1.0/(1.0+xn);
    gj3(j)=(1.0-xn)/(1.0+xn);
    gj3(je-j-1)=(1.0-xn)/(1.0+xn);
    xxn=(xnum-0.5)/xd;
    xn=0.25*((xxn)^3);
    fj1(j)=xn;
    fj1(je-2-j)=xn;
    fj2(j)=1.0/(1.0+xn);
    fj2(je-2-j)=1.0/(1.0+xn);
    fj3(j)=(1.0-xn)/(1.0+xn);
    fj3(je-2-j)=(1.0-xn)/(1.0+xn);
end
T=0;
radius = 10;
epsilon = 30;
sigma = 0.3;

for j = ja:jb
    for i = ia:ib
        xdist = ic -i;
        ydist = jc-j;
        dist = sqrt(xdist^2 + ydist^2);% Specfy the dielectric cylinder %

        if dist <= radius
            ga(i,j) = 1.0/ (epsilon + (sigma*dt/epsz));
            gb(i,j) = sigma*dt/epsz;
        end
    end
end

%%%%%%------------Main Loop Begins----------------%%%%%%%%%%%%for n = 1:nsteps 
for n = 1:nsteps
T = T + 1;
    
    %Calculate the Incident Ez
    for j = 2 :je
        ez_inc(j) = ez_inc(j)+0.5*(hx_inc(j-1)-hx_inc(j));
    end
    
    %Fourier Transform of the incident field
    for m = 1: nfreqs
        real_in(m) = real_in(m) + cos(arg(m)*T)*ez_inc(ja-1);
        imag_in(m) = imag_in(m) - sin(arg(m)*T)*ez_inc(ja-1);
    end
    %ABC for the incident bufer
    ez_inc(1) = ez_inc_low_m2;
    ez_inc_low_m2 = ez_inc_low_m1;
    ez_inc_low_m1 = ez_inc(2);
    
    ez_inc(je) = ez_inc_high_m2;
    ez_inc_high_m2 = ez_inc_high_m1;
    ez_inc_high_m1 = ez_inc(je-1);
    
    %Calculate Dz field
    for j=2:je
        for i=2:ie
            dz(i,j)=gi3(i)*gj3(j)*dz(i,j)+gi2(i)*gj2(j)*0.5*(hy(i,j)-hy(i-1,j)-hx(i,j)+hx(i,j-1));
        end
    end
   
    %pulse=-2.0*((t0-T)./spread).*exp(-1.*((t0-T)./spread)^2);
    pulse = exp(-0.5*((t0-T)./spread)^2);
   % pulse = sin(2*pi*1500*1e6*dt*T);

   % dz(ic-5,ic-5)= pulse;
    ez_inc(ja) = pulse;
    
    %incident the Dz values
    for i= ia:ie
        dz(i,ja) = dz(i,ja) + 0.5*hx_inc(ja-1);
        dz(i,jb) = dz(i,jb) - 0.5*hx_inc(jb);
    end
    
    %calculate Ez field
    for j=2:je
        for i=2:ie
            ez(i,j)=ga(i,j).*dz(i,j);
        end
    end
  %calcukate the fourier transform ex
  for j = 1:je
      for i = 1:ie
          for m = 1:nfreqs
              real_pt(m,i,j) = real_pt(m,i,j) + cos(arg(m)*T)*ez_inc(i,j);
              imag_pt(m,i,j) = imag_pt(m,i,j) + sin(arg(m)*T)*ez_inc(i,j);
          end
      end
  end
%%%%%%-----ABC-----------%%%%%%%%%%%%
  for j=1:je-1
     ez(1,j)=0;
     ez(ie-1,j)=0;
  end
 for i=1:ie-1
     ez(i,1)=0;
     ez(i,je-1)=0;
 end     
 %%%%%%-----ABC-----------%%%%%%%%%%%%
 for j = 1:je
     hx_inc(j) = hx_inc(j) + .5*(ez_inc(j) - ez_inc(j));
 end
 
 %calculate Hx field
    for j=1:je-1
        for i=1:ie-1
            curl_e=ez(i,j)-ez(i,j+1);
            ihx(i,j)=ihx(i,j)+fi1(i)*curl_e;
            hx(i,j)=fj3(j)*hx(i,j)+fj2(j)*0.5*(curl_e+ihx(i,j));
        end
    end
   %Incident Hx values
   for i = ia:ib
       hx(i,ja-1) = hx(i,ja-1) + .5*ez_inc(ja);
       hx(i,jb) = hx(i,jb) - .5*ez_inc(jb);
   end
   
   %calculate Hy field
     for j=1:je-1
        for i=1:ie-1
            curl_e=ez(i+1,j)-ez(i,j);
            ihy(i,j)=ihy(i,j)+fj1(j)*curl_e;
            hy(i,j)=fi3(i)*hy(i,j)+fi2(i)*0.5*(curl_e+ihy(i,j));
        end
     end
    %Incident Hy values
    for j = ja:jb
        hy(ia-1,j) = hy(ia-1,j) - 0.5*ez_inc(j);
        hy(ib,j) = hy(ib,j) + 0.5*ez_inc(j);
    end
   
    %calculate fourier amplitude and phase of the incident pulse
    
    for m = 1:nfreqs
        amp_in(m) = sqrt(real_in(m)^2) + imag_in(m)^2;
        phase_in(m) = atan2(imag_in(m),real_in(m));
    end
    
     for m = 1:nfreqs 
            for j = ja:jb
                  if ga(ic,j) < 1.00
                            amp(ic,j) = (1./amp_in(m))*sqrt(real_pt(m,ic,j)^2.0 + imag_pt(m,ic,j)^2.);
                            phase(ic,j) = atan2(imag_pt(m,ic,j),real_pt(m,ic,j));
                  end
            end
     end
     
                        
                    
    
 %%imagesc(ez);
 %title(['Time = ',num2str(n)]);
 %colorbar
 %pause(0.02);%%
 timestep=int2str(T);
    z = ez;
    zlim = ([-2.0 2.0]);
    surf(z(1:ie,1:je));
    title(['Ez at time step = ',timestep]);
    colorbar;
    pause(0.01)
 end