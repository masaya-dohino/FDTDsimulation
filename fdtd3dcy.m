
format long g
clear all;
close all;
IE = 50;
JE = 50;
KE = 50;
ia = 8;
ja = 8;
ka = 8;
ib = IE - ia;
jb = JE - ja;
kb = KE - ka;
NFREQS = 3;
ez_low_m1 = 0.;
ez_low_m2 = 0.;
ez_high_m1 = 0.;
ez_high_m2 = 0.;

freq = zeros(NFREQS);
arg = zeros(NFREQS);
real_pt = zeros(NFREQS,IE,JE);
imag_pt = zeros(NFREQS,IE,JE);
amp = zeros(IE,JE);
phase = zeros(IE,JE);
real_in = zeros(5);
imag_in = zeros(5);
amp_in = zeros(5);
phase_in = zeros(5);
ez_inc = zeros(JE);
hx_inc = zeros(JE);

radius = zeros(10);
epsilon = zeros(10);
sigma = zeros(10);


freq(1) = 10.e6;
freq(2) = 100.e6;
freq(3) = 433.e6;

cc = 2.99792458e8;                %speed of light 
mu_0 = 4.0*pi*1.0e-7;             %permeability of free space
eps_0 = 1.0/(cc*cc*mu_0);         %permitivity of free space
ic = floor(IE/2);
jc = floor(JE/2);
kc = floor(KE/2);
ddx = 0.01;
ddy = ddx;
dt = ddx/6e8;
ex = zeros(IE,JE,KE);
ey = zeros(IE,JE,KE);
ez = zeros(IE,JE,KE);

ix = zeros(IE,JE,KE);
iy = zeros(IE,JE,KE);
iz = zeros(IE,JE,KE);

dx = zeros(IE,JE,KE);
dy = zeros(IE,JE,KE);
dz = zeros(IE,JE,KE);

hx = zeros(IE,JE,KE);
hy = zeros(IE,JE,KE);
hz = zeros(IE,JE,KE);

idx = zeros(IE,JE,KE);
ihx = zeros(IE,JE,KE);
idy = zeros(IE,JE,KE);
ihy = zeros(IE,JE,KE);
idz = zeros(IE,JE,KE);
ihz = zeros(IE,JE,KE);

gax = ones(IE,JE,KE);
gay = ones(IE,JE,KE);
gaz = ones(IE,JE,KE);
gbx = zeros(IE,JE,KE);
gby = zeros(IE,JE,KE);
gbz = zeros(IE,JE,KE);
%--------------------specify the Dipole------------------------------------
for k = 11:kc+10
    gaz(:,k) = 0.;
end
gaz(ic,jc,kc) = 1.0;
%------------calculated the PML parameter-----------------------------
    gi1 = zeros(IE);
    gi2 = ones(IE);
    gi3 = ones(IE);
    fi1 = zeros(IE);
    fi2 = ones(IE);
    fi3 = ones(IE);
    
    gj1 = zeros(JE);
    gj2 = ones(JE);
    gj3 = ones(JE);
    fj1 = zeros(JE);
    fj2 = ones(JE);
    fj3 = ones(JE);
    
    gk1 = zeros(KE);
    fk1 = zeros(KE);
    gk2 = ones(KE);
    fk2 = ones(KE);
    gk3 = ones(KE);
    fk3 = ones(KE);
    
    
npml = 8;
for i = 1:npml
    xnum = npml-i;
    xd = npml;
    xxn = xnum/xd;
    xn = 0.33*((xxn)^3);
    fi1(i) = xn;
    fi1(IE-i-1) = xn;
    gi2(i) = 1.0/(1.0+xn);
    gi2(IE-1-i) = 1.0/(1.0+xn);
    gi3(i) = (1.0-xn)/(1.0+xn);
    gi3(IE-i-1) = (1.0-xn)/(1.0+xn);
    xxn = (npml-i-0.5)/npml;
    xn = 0.33*(xxn)^3.0;
    gi1(i) = xn;
    gi1(IE-2-i) = xn;
    fi2(i) = 1.0/(1.0+xn);
    fi2(IE-2-i) = 1.0/(1.0+xn);
    fi3(i) = (1.0-xn)/(1.0-xn);
    fi3(IE-2-i) = (1.0-xn)/(1.0+xn);
end
    
   for j = 1:npml
       xnum = npml-j;
    xd = npml;
    xxn = xnum/xd;
    xn = 0.33*((xxn)^3);
    fj1(j) = xn;
    fj1(JE-j-1) = xn;
    gj2(j) = 1.0/(1.0+xn);
    gj2(JE-1-j) = 1.0/(1.0+xn);
    gj3(j) = (1.0-xn)/(1.0+xn);
    gj3(JE-j-1) = (1.0-xn)/(1.0+xn);
    xxn = (npml-j-0.5)/npml;
    xn = 0.33*(xxn)^3.0;
    gj1(j) = xn;
    gj1(JE-2-j) = xn;
    fj2(j) = 1.0/(1.0+xn);
    fj2(JE-2-j) = 1.0/(1.0+xn);
    fj3(j) = (1.0-xn)/(1.0-xn);
    fj3(JE-2-j) = (1.0-xn)/(1.0+xn);
        
   end
    
for k = 1:npml
        xnum = npml-k;
    xd = npml;
    xxn = xnum/xd;
    xn = 0.33*((xxn)^3);
    fk1(k) = xn;
    fk1(KE-k-1) = xn;
    gk2(k) = 1.0/(1.0+xn);
    gk2(KE-1-k) = 1.0/(1.0+xn);
    gk3(k) = (1.0-xn)/(1.0+xn);
    gk3(KE-k-1) = (1.0-xn)/(1.0+xn);
    xxn = (npml-i-0.5)/npml;
    xn = 0.33*(xxn)^3.0;
    gk1(k) = xn;
    gk1(KE-2-k) = xn;
    fk2(k) = 1.0/(1.0+xn);
    fk2(KE-2-k) = 1.0/(1.0+xn);
    fk3(k) = (1.0-xn)/(1.0-xn);
    fk3(KE-2-k) = (1.0-xn)/(1.0+xn);
end

%----Specify dielectric sphere---%
epsilon(1) = 1.0;
sigma(1) = 0.;
numsph = 20;

RAD = 20;
SIGM = 0.3;
EPS = 30;

for n = 1:1:numsph
    
end

%Calculate gax,gbx
for i = ia:1:ib
    for j = ja:1:jb
        for k = ka:1:kb
            eps = epsilon(1);
            cond = sigma(1);
            ydist = jc-j-0.5;
            xdist = ic-i;
            zdist = kc-k;
           dist = sqrt((xdist)^2 + (ydist)^2 + (zdist)^2);
           for n = 2:1:numsph
               if dist <= RAD
                   epsilon(n) = EPS;
                   sigma(n) = SIGM;
                   eps = epsilon(n);
                   cond = sigma(n);
               end
           end
           gax(i,j,k) = 1./(eps + cond*dt/eps_0);
           gbx(i,j,k) = cond*dt/eps_0;
        end
    end
end

%Calculate gay,gby
for i = ia:1:ib
    for j = ja:1:jb
        for k = ka:1:kb
            eps = epsilon(1);
            cond = sigma(1);
            ydist = jc-j-0.5;
            xdist = ic-i;
            zdist = kc-k;
           dist = sqrt((xdist)^2 + (ydist)^2 + (zdist)^2);
           for n = 2:1:numsph
               if dist <= RAD
                   epsilon(n) = EPS;
                   sigma(n) = SIGM;
                   eps = epsilon(n);
                   cond = sigma(n);
               end
           end
           gay(i,j,k) = 1./(eps + cond*dt/eps_0);
           gby(i,j,k) = cond*dt/eps_0;
        end
    end
end

%Calculate gaz,gbz
for i = ia:1:ib
    for j = ja:1:jb
        for k = ka:1:kb
            eps = epsilon(1);
            cond = sigma(1);
            ydist = jc-j-0.5;
            xdist = ic-i;
            zdist = kc-k;
           dist = sqrt((xdist)^2 + (ydist)^2 + (zdist)^2);
           for n = 2:1:numsph
               if dist <= RAD
                   epsilon(n) = EPS;
                   sigma(n) = SIGM;
                   eps = epsilon(n);
                   cond = sigma(n);
               end
           end
           gaz(i,j,k) = 1./(eps + cond*dt/eps_0);
           gbz(i,j,k) = cond*dt/eps_0;
        end
    end
end

%---------------------Time instance specification-------------------------
t0 = 40.0;
spread = 10.0;
T = 0;
nsteps = 200;

for n = 1:1:nsteps
    T = T+1;
    %----------------------start of main FDTD loop-----------------------------
    %calculate incident buffer
    for j = 2:1:JE
        ez_inc(j) = ez_inc(j) + 0.5*(hx_inc(j-1) - hx_inc(j));
    end
    %Fourier Transform of the incident buffer
    for m = 1:1:NFREQS
        real_in(m) = real_in(m) + cos(arg(m)*T)*ez_inc(ja);
        imag_in(m) = imag_in(m) - sin(arg(m)*T)*ez_inc(ja);
    end
    
    %--Source
    pulse = exp(-0.5*((t0-T)/spread)^2.0);
    ez_inc(4) = pulse;
    
    %-- Boundary condition for the incident buffer
    ez_inc(1) = ez_low_m2;
    ez_low_m2 = ez_low_m1;
    ez_low_m1 = ez_inc(2);
    
    ez_inc(JE) = ez_high_m2;
    ez_high_m2 = ez_high_m1;
    ez_high_m1 = ez_inc(JE-1);
    
    %----------------------Calculate the Dx field------------------------------
    for k = 2:1:KE
        for j = 2:1:JE
            for i = 1:1:IE 
                curl_h = (hz(i,j,k) - hz(i,j-1,k) - hy(i,j,k) + hy(i,j,k-1));
                idx(i,j,k) = idx(i,j,k) + curl_h;
                dx(i,j,k) = gj3(j)*gk3(k)*dx(i,j,k) + gj2(j)*gk2(k)*0.5*(curl_h + gi1(i)*idx(i,j,k));
            end
        end
    end
    %------------------------Calculate of the Dy field--------------------------
    for k = 2:1:KE
        for j = 1:1:JE
            for i = 2:1:IE
                curl_h = (hx(i,j,k) - hx(i,j,k-1) - hz(i,j,k) + hz(i-1,j,k));
                idy(i,j,k) = idy(i,j,k) + curl_h;
                dy(i,j,k) = gi3(i)*gk3(k)*dy(i,j,k) + gi2(i)*gk2(k)*0.5*(curl_h + gj1(j)*idy(i,j,k));
            end 
        end
    end
    
    %--- Incident Dy --%
    for i = ia:1:ib
        for j = ja:1:jb
            dy(i,j,ka) = dy(i,j,ka) - 0.5*hx_inc(j);
            dy(i,j,kb+1) = dy(i,j,kb+1) + 0.5*hx_inc(j);
        end
    end
    %---------------------Calculate the Dz field--------------------------------
    for k=1:1:KE
        for j=2:1:JE
            for i=2:1:IE
                curl_h = (hy(i,j,k) - hy(i-1,j,k) - hx(i,j,k) + hx(i,j-1,k));
                idz(i,j,k) = idz(i,j,k) + curl_h;
                dz(i,j,k) = gi3(i)*gj3(j)*dz(i,j,k) + gi2(i)*gj2(j)*0.5*(curl_h + gk1(k)*idz(i,j,k));
            end
        end
    end
    
    %--Incident Dz--%
    for i = ia:1:ib
        for k = ka:1:kb
            dz(i,ja,k) = dz(i,ja,k) + 0.5*hx_inc(ja-1);
            dz(i,jb,k) = dz(i,jb,k) - 0.5*hx_inc(jb);
        end
    end
    %--source ?
    
    %--------------Calculate the E from D field---------------------------------
    for k=1:1:KE-1
        for j=1:1:JE-1
            for i=1:1:IE-1
                ex(i,j,k) = gax(i,j,k).*(dx(i,j,k)- ix(i,j,k));
                ix(i,j,k) = ix(i,j,k) + gbx(i,j,k)*ex(i,j,k);
                ey(i,j,k) = gay(i,j,k).*(dy(i,j,k)-iy(i,j,k));
                iy(i,j,k) = iy(i,j,k) + gby(i,j,k)*ey(i,j,k);
                ez(i,j,k) = gaz(i,j,k).*(dz(i,j,k)-iz(i,j,k));
                iz(i,j,k) = iz(i,j,k) + gbz(i,j,k)*ez(i,j,k);
            end
        end
    end     
    %-----Calculate the Fourier Transform of Ex---%
    for j = 1:1:JE
        for i = 1:1:IE
            for m = 1:1:NFREQS
                real_pt(m,i,j) = real_pt(m,i,j) + cos(arg(m)*T)*ez(i,j,kc);
                imag_pt(m,i,j) = imag_pt(m,i,j) + sin(arg(m)*T)*ez(i,j,kc);
            end
        end
    end
    
    %---Calculate the Incident field---%
    for j = 1:1:JE-1
        hx_inc(j) = hx_inc(j) + 0.5*(ez_inc(j) - ez_inc(j+1));
    end
    
    %----------------Calculate the Hx field------------------------------------
    for k=1:1:KE-1
        for j=1:1:JE-1
            for i=1:1:IE
                curl_e = (ey(i,j,k+1)-ey(i,j,k)-ez(i,j+1,k)+ez(i,j,k));
                ihx(i,j,k) = ihx(i,j,k) + curl_e;
                hx(i,j,k) = fj3(j)*fk3(k)*hx(i,j,k) + fj2(j)*fk2(k)*0.5*(curl_e + fi1(i)*ihx(i,j,k));
            end
        end
    end
    
    
    %--Incident Hx--%
    for i = ia:1:ib
        for k =ka:1:kb
            hx(i,ja-1,k) = hx(i,ja-1,k) + 0.5*ez_inc(ja);
            hx(i,jb,k) = hx(i,jb,k) - 0.5*ez_inc(jb);
        end
    end
    %--------------Calculate the Hy field--------------------------------------
    for k=1:1:KE-1
        for j=1:1:JE
            for i=1:1:IE-1
                curl_e = (ez(i+1,j,k)-ez(i,j,k)-ex(i,j,k+1)+ex(i,j,k));
                ihy(i,j,k) = ihy(i,j,k) + curl_e;
                hy(i,j,k) = fi3(i)*fk3(k)*hy(i,j,k) + fi2(i)*fk3(k)*0.5*(curl_e + fj1(j)*ihy(i,j,k));
            end
        end
    end
    
    %---Incident Hy field---%
    for j = ja:1:jb
        for k =ka:1:kb
            hy(ia-1,j,k) = hy(ia-1,j,k) + 0.5*ez_inc(j);
            hy(ib,j,k) = hy(ib,j,k) - 0.5*ez_inc(j);
        end
    end
    %----------------Calculate the Hz field------------------------------------
    for k=1:1:KE
        for j=1:1:JE-1
            for i=1:1:IE-1
                curl_e = (ex(i,j+1,k)-ex(i,j,k)-ey(i+1,j,k)+ey(i,j,k));
                ihz(i,j,k) = ihz(i,j,k) + curl_e;
                hz(i,j,k) = fi3(i)*fj3(j)*hz(i,j,k) + fi2(i)*fj2(j)*0.5*(curl_e + fk1(k)*ihz(i,j,k));
            end
        end
    end
    %----------------End of the main FDTD_loop--------------------------------
    
    %Calculate Fourier Amplitude and the phase of the total field
    for m = 1:1:NFREQS
        for j = ja:1:jb
         if gaz(ic,j,kc) < 1.00
            amp(ic,j) = (1.0/amp_in(m))*sqrt(real_pt(m,ic,j)^2 + imag_pt(m,ic,j)^2);
         end
        end
    end
    timestep=int2str(T);
    surf(ez(1:IE,1:JE,kc));
    title(['Ez at time step = ',timestep]);
    colorbar;
    pause(0.01)
end