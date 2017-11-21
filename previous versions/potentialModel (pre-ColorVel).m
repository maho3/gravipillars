clear all; close all

% Models the position, velocity, and acceleration of any N number of
% randomly placed 1 kg particles in a gravitational field of strength mag
% (G*mass) and planetary size of radius a. The source of the gravitational
% field moves diagonally to model a dynamic system.

%%
graphTitle='graviPillars';  %set title of graph

T=50000;                    %set number of frames
step=0.01;                  %set step time length (inv. proportional to precision)
catLength=10;               %set length of caterpillars
N=8;                        %set number of projectiles


o=[0,0,0];                  %set start position of origin
wscale=2;                   %set scale of plot in relation to spawn window

xmax=10;                    %set spawn window
xmin=-10;
ymax=10;
ymin=-10;
zmax=10;
zmin=-10;


toggColors=1;               %toggle colors
cc=lines(N+1);              %define type of colors


toggField=1;                %gravity=0, ideal spring=1, real spring=2
magO=5000;                  %set strength (G*mass) of gravity field
a=5;                        %set origin planet width
    
springBreak=3*a;            %set breakage radius of real spring

toggBounce=0;               %toggle bounce at edges of plot

toggFrictionField=1;        %toggle friction field (confines particles to sphere)
frictionRadius=3*a;         %set radius of friction field
frictionStrength=.1;        %set strength of friction field

toggParticleField=1;        %turn particle field on (will produce same type as toggField)
magP=100;                   %set magnitude of particle field
pA=1;                       %set radius of particle "planets"

toggParticleDeflect=0;      %turn on particle deflection (bounce particles off each other)
pDRad=1;                    %set deflection radius
dBuffer=3;                  %set deflection frame buffer

toggParticleRepulse=0;      %turn on particle repulsion (repel as a function of 1/distance^3)
magR=100;                   %set magnitude of repulsion


toggSphere=0;               %toggle planet sphere plot (will slow comp)
toggFricSphere=0;           %toggle friction sphere plot (will slow comp)

toggMoveOrigin=0;           %toggle origin movement in diagonal fashion
oStep=0.01;                 %set movement of origin

%%

graphXMin=wscale*xmin;       %derive graph window
graphXMax=wscale*xmax;
graphYMin=wscale*ymin;
graphYMax=wscale*ymax;
graphZMin=wscale*zmin;
graphZMax=wscale*zmax;


figure,hold on

g(1)=plot3(o(1),o(2),o(3),'g+');    %plot initial origin
grid on
xlabel('x')
ylabel('y')
zlabel('z')
view([1,1,1])
title(graphTitle)

if toggSphere==1                    %plot planet sphere
    syms u v t;
    u=a*sin(t);
    v=a*cos(t);
    i(1)=ezplot3(u+o(1),v+o(2),o(3));
    i(2)=ezplot3(o(1),u+o(2),v+o(3));
    i(3)=ezplot3(u+o(1),o(2),v+o(3));
    title(graphTitle)
end

if toggFricSphere==1                %plot friction sphere
    syms c d e;
    c=frictionRadius*sin(e);
    d=frictionRadius*cos(e);
    j(1)=ezplot3(c+o(1),d+o(2),o(3));
    j(2)=ezplot3(o(1),c+o(2),d+o(3));
    j(3)=ezplot3(c+o(1),o(2),d+o(3));
    title(graphTitle)
end

axis([graphXMin,graphXMax,graphYMin,graphYMax,graphZMin,graphZMax])

for n=1:N                           %set random inital values
    x{n}(1)=rand*(xmax-xmin)+xmin;
    y{n}(1)=rand*(ymax-ymin)+ymin;
    z{n}(1)=rand*(zmax-zmin)+zmin;
    
    posVecInit=[o(1)-x{n}(1),o(2)-y{n}(1),o(3)-z{n}(1)];
    escVel=sqrt(magO/norm(posVecInit));
    
    posVecInitDirec=posVecInit/norm(posVecInit);
    
    vx{n}(1)=(rand*2*escVel-escVel)*posVecInitDirec(1); %avoid initial escape velocity
    vy{n}(1)=(rand*2*escVel-escVel)*posVecInitDirec(2);
    vz{n}(1)=(rand*2*escVel-escVel)*posVecInitDirec(3);
    
    for m=1:N   
        prevDeflect{n}(m)=-dBuffer;
    end
    
    if toggColors==1
        for m=1:N
            color{m}=cc(m,:);
        end
    end
end


for k = 1:T
    for n=1:N                       %math and physics to calculate position
        posVec=[o(1)-x{n}(k),o(2)-y{n}(k),o(3)-z{n}(k)];
        d=norm(posVec);             %determine radius
        normPosVec=posVec/d;
        
        if toggField==0             %calculate origin field force
            if (d>=a)
                f=magO/d^2;
            else
                f=magO*d/a^3;
            end
        elseif toggField==1
            f=magO*d/a^3;
        elseif toggField==2
            c=4/a;
            aa=a/2;                 %don't bother deriving this equation
            q=1/2*exp(-2*c*aa)*(a*c*exp(c*aa+aa)-sqrt(a*c*(a*c-4)*exp(2*c*aa+2*aa))-2*exp(c*aa+aa));
            r=q*magO/a^2;
            if d<=a/2
                f=magO/a^3*(d-aa)+r/(2*q);
            elseif d>=springBreak
                f=0;
            else
                f=r/(q+exp(-d+a/2));
            end
        end
        
        ax{n}(k)=f*normPosVec(1);   %determine acceleration in 3D
        ay{n}(k)=f*normPosVec(2);
        az{n}(k)=f*normPosVec(3);
        
        if toggParticleField==1;        %calculate particle fields
            for b=1:N                   %between all particles and not themselves
                if b~=n
                    partPosVec=[x{b}(k)-x{n}(k),y{b}(k)-y{n}(k),x{b}(k)-z{n}(k)];
                    pD=norm(partPosVec);             %determine distance between particles
                    normPartPosVec=partPosVec/pD;
                    
                    if toggField==0             %calculate particle field force
                        if (pD>=pA)
                            pF=magP/pD^2;
                        else
                            pF=magP*pD/pA^3;
                        end
                    elseif toggField==1
                        pF=magP*pD/pA^3;
                    elseif toggField==2
                        pC=4/pA;
                        pAA=pA/2;               %don't bother deriving this equation
                        pQ=1/2*exp(-2*pC*pAA)*(pA*pC*exp(pC*pAA+pAA)-sqrt(pA*pC*(pA*pC-4)*exp(2*pC*pAA+2*pAA))-2*exp(pA*pAA+pAA));
                        pR=pQ*magP/pA^2;
                        if pD<=pA/2
                            pF=magP/pA^3*(pD-pAA)+pR/(2*pQ);
                        elseif pD>=springBreak
                            pF=0;
                        else
                            pF=pR/(pQ+exp(-pD+pA/2));
                        end
                    end
                                    %add particle acceleration vector to
                                    %total acceleration vector
                    ax{n}(k)=ax{n}(k)+pF*normPartPosVec(1);
                    ay{n}(k)=ay{n}(k)+pF*normPartPosVec(2);
                    az{n}(k)=az{n}(k)+pF*normPartPosVec(3);
                end
            end
        end
        
        if toggParticleRepulse==1;
            for rB=1:N
                if rB~=n
                    rPosVec=[x{rB}(k)-x{n}(k),y{rB}(k)-y{n}(k),x{rB}(k)-z{n}(k)];
                    rD=norm(rPosVec);             %determine distance between particles
                    normRPosVec=rPosVec/rD;
                    
                    rF=magR/rD^3;
                                    %add acceleration vector
                    ax{n}(k)=ax{n}(k)-rF*normRPosVec(1);
                    ay{n}(k)=ay{n}(k)-rF*normRPosVec(2);
                    az{n}(k)=az{n}(k)-rF*normRPosVec(3);
                end
            end
        end
        
        if toggFrictionField==1     %calculate and add friction accel
            if d>=frictionRadius
                velMag=norm([vx{n}(k),vy{n}(k),vz{n}(k)]);
                velDirect=[vx{n}(k),vy{n}(k),vz{n}(k)]/velMag;
                
                ff=frictionStrength*(velMag)^2;
                
                ax{n}(k)=ax{n}(k)-ff*velDirect(1);
                ay{n}(k)=ay{n}(k)-ff*velDirect(2);
                az{n}(k)=az{n}(k)-ff*velDirect(3);
                
                
            end
        end
        
        vx{n}(k+1)=vx{n}(k)+step*ax{n}(k);  %calculate velocity
        vy{n}(k+1)=vy{n}(k)+step*ay{n}(k);
        vz{n}(k+1)=vz{n}(k)+step*az{n}(k);
        
        if toggBounce==1            %velocity reverses at edges
            if x{n}(k)>=graphXMax
                vx{n}(k+1)=-vx{n}(k+1);
            elseif x{n}(k)<=graphXMin
                vx{n}(k+1)=-vx{n}(k+1);
            end
            if y{n}(k)>=graphYMax
                vy{n}(k+1)=-vy{n}(k+1);
            elseif y{n}(k)<=graphYMin
                vy{n}(k+1)=-vy{n}(k+1);
            end
            if z{n}(k)>=graphZMax
                vz{n}(k+1)=-vz{n}(k+1);
            elseif z{n}(k)<=graphZMin
                vz{n}(k+1)=-vz{n}(k+1);
            end
        end
        
        if toggParticleDeflect==1
            for mm=1:(n-1)
                if (k-prevDeflect{n}(mm))>=dBuffer  %ensure particles don't trap themselves in alternating velocity
                    dPosVecNM=[x{n}(k)-x{mm}(k),y{n}(k)-y{mm}(k),z{n}(k)-z{mm}(k)];
                    magDPosVec=norm(dPosVecNM);     %determine distance between particles
                    dPosVecDirecNM=dPosVecNM/magDPosVec;    %determine direction vectors of contact
                    dPosVecDirecMN=-dPosVecDirecNM;
                    
                    if magDPosVec<=pDRad            %if particles are within a each other's radii
                        vn=[vx{n}(k+1),vy{n}(k+1),vz{n}(k+1)];  %determine velocity vectors
                        vm=[vx{mm}(k+1),vy{mm}(k+1),vz{mm}(k+1)];
                        
                        vnCollide=dot(vn,dPosVecDirecNM)*dPosVecDirecNM;    %determine velocity vectors in direction of contact
                        vmCollide=dot(vm,dPosVecDirecMN)*dPosVecDirecMN;
                        
                        vx{n}(k+1)=vx{n}(k+1)-2*vnCollide(1);   %reverse contact velocity vectors from total velocities
                        vy{n}(k+1)=vy{n}(k+1)-2*vnCollide(2);
                        vz{n}(k+1)=vz{n}(k+1)-2*vnCollide(3);
                        
                        vx{mm}(k+1)=vx{mm}(k+1)-2*vmCollide(1);
                        vy{mm}(k+1)=vy{mm}(k+1)-2*vmCollide(2);
                        vz{mm}(k+1)=vz{mm}(k+1)-2*vmCollide(3);
                        
                        prevDeflect{n}(mm)=k;               %mark that a collision has occurred
                        prevDeflect{mm}(n)=k;
                    end
                end
            end
        end
        
        x{n}(k+1)=x{n}(k)+step*vx{n}(k+1);  %calculate position
        y{n}(k+1)=y{n}(k)+step*vy{n}(k+1);
        z{n}(k+1)=z{n}(k)+step*vz{n}(k+1);
    end
    
    
    for n=1:N                       %plot
        if toggColors==1
            h{n}(k)=plot3(x{n}(k),y{n}(k),z{n}(k),'*','color',color{n});
        else
            h{n}(k)=plot3(x{n}(k),y{n}(k),z{n}(k),'*');
        end
        axis([graphXMin,graphXMax,graphYMin,graphYMax,graphZMin,graphZMax])
        
        if k>catLength
            delete(h{n}(k-catLength));
        end
    end
    if toggMoveOrigin==1
        if rem(k,200)<=100
            o(1)=o(1)+oStep;        %move origin
            o(2)=o(2)+oStep;
            o(3)=o(3)+oStep;
        elseif rem(k,200)>100
            o(1)=o(1)-oStep;        %move origin
            o(2)=o(2)-oStep;
            o(3)=o(3)+oStep;
        end
        
        g(k+1)=plot3(o(1),o(2),o(3),'g+');  %plot origin
        delete(g(k));
        
        if toggSphere==1            %plot planet sphere
            delete(i(3*(k-1)+1));
            delete(i(3*(k-1)+2));
            delete(i(3*(k-1)+3));
            i(3*k+1)=ezplot3(u+o(1),v+o(2),o(3));
            i(3*k+2)=ezplot3(o(1),u+o(2),v+o(3));
            i(3*k+3)=ezplot3(u+o(1),o(2),v+o(3));
            title(graphTitle)
        end
        if toggFricSphere==1        %plot friction sphere
            delete(j(3*(k-1)+1));
            delete(j(3*(k-1)+2));
            delete(j(3*(k-1)+3));
            j(3*k+1)=ezplot3(c+o(1),d+o(2),o(3));
            j(3*k+2)=ezplot3(o(1),c+o(2),d+o(3));
            j(3*k+3)=ezplot3(c+o(1),o(2),d+o(3));
            title(graphTitle)
        end
    end
    
    disp(k)
    pause(0.000000001)               %pause to allow for animation
end

close all;