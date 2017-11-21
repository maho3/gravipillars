clear all; close all

% Models the position, velocity, and acceleration of any N number of 
% randomly placed 1 kg particles in a gravitational field of strength mag 
% (G*mass) and planetary size of radius a. The source of the gravitational
% field moves diagonally to model a dynamic system. 

%% 
T=1000;                 %set number of frames
o=[-2.5,-2.5];          %set start position of origin
oStep=0.05;             %set movement of origin
mag=7000;               %set strength (G*mass) of gravity field
step=0.01;              %set step length (inv. proportional to precision)
a=3;                    %set planet width
N=10;                   %set number of projectiles

xmax=5;                 %set spawn window
xmin=-5;
ymax=5;
ymin=-5;
%% 

graphXMax=5*xmax;       %derive graph window
graphXMin=5*xmin;
graphYMax=5*ymax;
graphYMin=5*ymin;


%  syms v w;
% w = sqrt(a^2-(v-o(1))^2)+o(2);

figure,hold on

g(1)=plot(o(1),o(2),'g+');                          %plot origin

%  i(1)=ezplot(w);
% axis([graphXMin,graphXMax,graphYMin,graphYMax])

for n=1:N                                           %set inital values
    x{n}(1)=rand*(xmax-xmin)+xmin;
    y{n}(1)=rand*(ymax-ymin)+ymin;
    
    posVecInit=[o(1)-x{n}(1),o(2)-y{n}(1)];
    escVel=sqrt(mag/norm(posVecInit));
    vx{n}(1)=rand*2*escVel-escVel;                           %avoid initial escape velocity
    vy{n}(1)=rand*2*escVel-escVel;
end

                                                      %math and physics
for k = 1:T
    for n=1:N
        posVec=[o(1)-x{n}(k),o(2)-y{n}(k)];
        d=norm(posVec);
        normPosVec=posVec/d;
        
        if (d>=a)
            f=mag/d^2;
        else
            f=mag*d/a^3;
        end
        
        if x{n}(k)>=graphXMax
            vx{n}(k)=-vx{n}(k);
        elseif x{n}(k)<=graphXMin
            vx{n}(k)=-vx{n}(k);
        end
        
        if y{n}(k)>=graphYMax
            vy{n}(k)=-vy{n}(k);
        elseif y{n}(k)<=graphYMin
            vy{n}(k)=-vy{n}(k);
        end
        
        ax{n}(k)=f*normPosVec(1);
        ay{n}(k)=f*normPosVec(2);
        
        vx{n}(k+1)=vx{n}(k)+step*ax{n}(k);
        vy{n}(k+1)=vy{n}(k)+step*ay{n}(k);
        
        x{n}(k+1)=x{n}(k)+step*vx{n}(k);
        y{n}(k+1)=y{n}(k)+step*vy{n}(k);

    end
                                                      %plot
    for n=1:N 
        h{n}(k)=plot(x{n}(k),y{n}(k),'*');
        axis([graphXMin,graphXMax,graphYMin,graphYMax])
        
        if k>10
            delete(h{n}(k-10));
        end
    end
    
    if rem(k,200)<=100
        o(1)=o(1)+oStep;                             %move origin
        o(2)=o(2)+oStep;
    elseif rem(k,200)>100
        o(1)=o(1)-oStep;                             %move origin
        o(2)=o(2)-oStep;
    end
    
    g(k+1)=plot(o(1),o(2),'g+');
    delete(g(k));
    
%      w = sqrt(a^2-(v-o(1))^2)+o(2);
%      i(k+1)=ezplot(w);
%      axis([graphXMin,graphXMax,graphYMin,graphYMax])
%      delete(i(k));
   
    
    pause(0.00001)
end

close all;