clear all; close all;

% Models the position, velocity, and acceleration of any N number of
% randomly placed 1 kg particles in a gravitational field of strength mag
% (G*mass) and planetary size of radius a. The source of the gravitational
% field moves diagonally to model a dynamic system.
%%

s = osc_new_server(3334);

init{1} = 'go';

aD = osc_new_address('127.0.0.1', 3335);


%%
graphTitle='graviPillars';  %set title of graph

T=50000;                    %set number of frames
step=0.01;                  %set step time length (inv. proportional to precision)
catLength=10;               %set length of caterpillars
N=16;                        %set number of projectiles

numOrigin=2;
o{1}=[0,0,0];                  %set start position of origin
o{2}=[0,0,0];
o{3}=[0,0,0];                 %set start position of origin
wscale=2;                   %set scale of plot in relation to spawn window

toggColors=1;               %toggle colors
cc=lines(N+1);              %define type of colors


toggField=0;                %gravity=0, ideal spring=1, real spring=2
magO(1)=1000000;                  %set strength (G*mass) of gravity field
magO(2)=1000000;
magO(3)=1000000;
a=5;                        %set origin planet width

xmax=a;                    %set spawn window
xmin=-a;
ymax=a;
ymin=-a;
zmax=a;
zmin=-a;

springBreak=3*a;            %set breakage radius of real spring

toggBounce=0;               %toggle bounce at edges of plot

toggFrictionField=1;        %toggle friction field (confines particles to sphere)
frictionRadius=6*a;         %set radius of friction field
frictionStrength=.1;        %set strength of friction field

toggParticleField=0;        %turn particle field on (will produce same type as toggField)
magP=100;                   %set magnitude of particle field
pA=1;                       %set radius of particle "planets"

toggParticleDeflect=0;      %turn on particle deflection (bounce particles off each other)
pDRad=1;                    %set deflection radius

toggParticleRepulse=0;      %turn on particle repulsion (repel as a function of 1/distance^3)
magR=1000;                   %set magnitude of repulsion


toggSphere=0;               %toggle planet sphere plot (will slow comp)
toggFricSphere=0;           %toggle friction sphere plot (will slow comp)

toggMoveOrigin=0;           %toggle origin movement in diagonal fashion
oStep=0.01;                 %set movement of origin
toggMagOAdjust=1;

%%

graphXMin=wscale*xmin;       %derive graph window
graphXMax=wscale*xmax;
graphYMin=wscale*ymin;
graphYMax=wscale*ymax;
graphZMin=wscale*zmin;
graphZMax=wscale*zmax;

k=0;
while strcmp(init{1},'go')
    k=k+1;
    
    oSC=osc_recv(s);
    holdVar=oSC{1,1}.data;
    
    if oSC{1}.path=='stop'
        init{1}='stop';
    else
        toggBounce=holdVar{1,1};
        toggParticleField=holdVar{1,2};
        toggParticleDeflect=holdVar{1,3};
        toggParticleRepulse=holdVar{1,4};
        toggMagOAdjust=holdVar{1,5};
        step=holdVar{1,6};
        magR=holdVar{1,7};
        magP=holdVar{1,8};
        frictionStrength=holdVar{1,9};
        frictionRadius=holdVar{1,10};
        a=holdVar{1,11};
        
        for m=1:numOrigin
            magO(m)=holdVar{1,8+4*m};
            
            prevO{m}(1)=o{m}(1);
            prevO{m}(2)=o{m}(2);
            prevO{m}(3)=o{m}(3);
            
            o{m}(1)=holdVar{1,9+4*m}*20;
            o{m}(2)=holdVar{1,10+4*m}*20;
            o{m}(3)=holdVar{1,11+4*m}*20;
        end
    end
    
    if k==1
        for n=1:N                           %set random inital values
            x{n}(1)=rand*(xmax-xmin)+xmin;
            y{n}(1)=rand*(ymax-ymin)+ymin;
            z{n}(1)=rand*(zmax-zmin)+zmin;
            
            posVecInit=[-x{n}(1),-y{n}(1),-z{n}(1)];
            escVel=sqrt((magO(1)+magO(2))/norm(posVecInit));
            
            posVecInitDirec=posVecInit/norm(posVecInit);
            
            v=rand*2*escVel-escVel;
            
            vx{n}(1)=v*posVecInitDirec(1); %avoid initial escape velocity
            vy{n}(1)=v*posVecInitDirec(2);
            vz{n}(1)=v*posVecInitDirec(3);
            
            for m=1:N
                color{m}=cc(m,:);
            end
        end
    end
    
    if toggMagOAdjust==1
        if k~=1
            for m=1:numOrigin
                diffO{m}(1)=abs(o{m}(1)-prevO{m}(1));
                diffO{m}(2)=abs(o{m}(2)-prevO{m}(2));
                diffO{m}(3)=abs(o{m}(3)-prevO{m}(3));
                
                dO(m)=norm(diffO{m});
                
                magO(m)=magO(m)*(1+dO(m)/a);
            end
        end
    end
    
            
            
            
    for n=1:N                       %math and physics to calculate position
        
        ax{n}(k)=0;   %determine acceleration in 3D
        ay{n}(k)=0;
        az{n}(k)=0;
        
        for m=1:numOrigin
            posVec=[o{m}(1)-x{n}(k),o{m}(2)-y{n}(k),o{m}(3)-z{n}(k)];
            d=norm(posVec);             %determine radius
            normPosVec=posVec/d;
            
            if toggField==0             %calculate origin field force
                if (d>=a)
                    f=magO(m)/d^2;
                else
                    f=magO(m)*d/a^3;
                end
            elseif toggField==1
                f=magO(m)*d/a^3;
            elseif toggField==2
                c=4/a;
                aa=a/2;                 %don't bother deriving this equation
                q=1/2*exp(-2*c*aa)*(a*c*exp(c*aa+aa)-sqrt(a*c*(a*c-4)*exp(2*c*aa+2*aa))-2*exp(c*aa+aa));
                r=q*magO(m)/a^2;
                if d<=a/2
                    f=magO(m)/a^3*(d-aa)+r/(2*q);
                elseif d>=springBreak
                    f=0;
                else
                    f=r/(q+exp(-d+a/2));
                end
            end
            
            ax{n}(k)=ax{n}(k)+f*normPosVec(1);   %determine acceleration in 3D
            ay{n}(k)=ay{n}(k)+f*normPosVec(2);
            az{n}(k)=az{n}(k)+f*normPosVec(3);
        end
        
        if toggParticleField==1;        %calculate particle fields
            for b=1:N                   %between all particles and not themselves
                if b~=n
                    partPosVec=[x{b}(k)-x{n}(k),y{b}(k)-y{n}(k),z{b}(k)-z{n}(k)];
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
                    rPosVec=[x{rB}(k)-x{n}(k),y{rB}(k)-y{n}(k),z{rB}(k)-z{n}(k)];
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
            friction=1;
            for m=1:numOrigin
                posVecF=[o{m}(1)-x{n}(k),o{m}(2)-y{n}(k),o{m}(3)-z{n}(k)];
                dF(m)=norm(posVecF);
                
                if dF(m)<=frictionRadius
                    friction=0;
                end
            end
            
            if friction==1;
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
            if (x{n}(k)>=graphXMax)&&(vx{n}(k+1)>0)
                vx{n}(k+1)=-vx{n}(k+1);
            elseif (x{n}(k)<=graphXMin)&&(vx{n}(k+1)<0)
                vx{n}(k+1)=-vx{n}(k+1);
            end
            if (y{n}(k)>=graphYMax)&&(vy{n}(k+1)>0)
                vy{n}(k+1)=-vy{n}(k+1);
            elseif (y{n}(k)<=graphYMin)&&(vy{n}(k+1)<0)
                vy{n}(k+1)=-vy{n}(k+1);
            end
            if (z{n}(k)>=graphZMax)&&(vz{n}(k+1)>0)
                vz{n}(k+1)=-vz{n}(k+1);
            elseif (z{n}(k)<=graphZMin)&&(vz{n}(k+1)<0)
                vz{n}(k+1)=-vz{n}(k+1);
            end
        end
        
        if toggParticleDeflect==1
            for mm=1:(n-1)
                dPosVecNM=[x{n}(k)-x{mm}(k),y{n}(k)-y{mm}(k),z{n}(k)-z{mm}(k)];
                magDPosVec=norm(dPosVecNM);     %determine distance between particles
                dPosVecDirecNM=dPosVecNM/magDPosVec;    %determine direction vectors of contact
                dPosVecDirecMN=-dPosVecDirecNM;
                
                vn=[vx{n}(k+1),vy{n}(k+1),vz{n}(k+1)];  %determine velocity vectors
                vm=[vx{mm}(k+1),vy{mm}(k+1),vz{mm}(k+1)];
                
                vnCollide=dot(vn,dPosVecDirecNM)*dPosVecDirecNM;    %determine velocity vectors in direction of contact
                vmCollide=dot(vm,dPosVecDirecMN)*dPosVecDirecMN;
                
                if magDPosVec<=pDRad           %if particles are within a each other's radii
                    if dot(vn,dPosVecDirecNM)<0 && dot(vm,dPosVecDirecMN)<0
                        vx{n}(k+1)=vx{n}(k+1)-2*vnCollide(1);   %reverse contact velocity vectors from total velocities
                        vy{n}(k+1)=vy{n}(k+1)-2*vnCollide(2);
                        vz{n}(k+1)=vz{n}(k+1)-2*vnCollide(3);
                        
                        vx{mm}(k+1)=vx{mm}(k+1)-2*vmCollide(1);
                        vy{mm}(k+1)=vy{mm}(k+1)-2*vmCollide(2);
                        vz{mm}(k+1)=vz{mm}(k+1)-2*vmCollide(3);
                    end
                end
            end
        end
        
        x{n}(k+1)=x{n}(k)+step*vx{n}(k+1);  %calculate position
        y{n}(k+1)=y{n}(k)+step*vy{n}(k+1);
        z{n}(k+1)=z{n}(k)+step*vz{n}(k+1);
    end
    
    
    for n=1:N                       %plot
        pos{n}=struct('path',sprintf('/cat%d',n),'tt','fff','data',{{x{n}(k),y{n}(k),z{n}(k)}});
        disp(step);
    end
    
    err=osc_send(aD,pos);
    
end

osc_free_server(s);
close all;