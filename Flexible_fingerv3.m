clc
clear all
%%defines r values
r1=[8*10^-3 8*10^-3 8*10^-3 8*10^-3 8*10^-3] ;
r2=[10*10^-3 12*10^-3 12*10^-3 12*10^-3 12*10^-3] ;
r3=[10*10^-3 8*10^-3 8*10^-3 8*10^-3 8*10^-3] ;
r4=[10*10^-3 13*10^-3 13*10^-3 13*10^-3 13*10^-3];
%%hand dimensions
HL=0.175;
HB=0.075;
%%dimensinos with respect to base
thumb=[r1(1)+3*10^-3 0.251*HL-0.015 0.196*HL 0.158*HL];
index=[sqrt((0.374*HL)^2+(0.126*HB)^2)-0.055 0.265*HL 0.143*HL 0.097*HL];
middle=[0.373*HL-0.055 0.277*HL 0.170*HL 0.108*HL];
little=[sqrt((0.336*HL)^2+(0.077*HB)^2)-0.03 0.259*HL 0.165*HL 0.107*HL];
ring=[sqrt((0.295*HL)^2+(0.179*HB)^2)-0.03 0.206*HL 0.117*HL 0.093*HL];
k=14*10^-3 ;
l=14*10^-3 ;
E=3500000000 ; %modulus of elasticity
FP=[2 2 2 2 2]; %deflection forces
F=FP';
step_size=1000;

x=zeros(5,step_size+1); %position vector

L=[thumb; index; middle; little; ring] %dimension matrix


for m=1:5 %dimension
    l1(m)=L(m,1);
    l2(m)=L(m,2);
    l3(m)=L(m,3);
    l4(m)=L(m,4);
    
end
for m=1:1 %total length
    TL(m)=sum(L(m,1:4));
end
for m=2:5 %total length
    TL(m)=sum(L(m,:))+r4(m);
end
for m=1:5 %position vector
   for i=1:step_size
    x(m,i+1)=x(m,i)+TL(1,m)/step_size;
   end
end
for m=1:1 %inertia matrix for thumb
for i=1:step_size+1 %inertia matrix
    if x(m,i)<=L(m,1)-r1(m)
        I(m,i)=1/12*k^3*l;
    end
    
    if x(m,i)<L(m,1) && x(m,i)>=L(m,1)-r1(m)
        I(m,i)=1/12*l*(k-sqrt(r1(m)^2-(L(m,1)-x(m,i))^2))^3;
    end
     if x(m,i)>=L(m,1) && x(m,i)<L(m,1)+r1(m)
        I(m,i)=1/12*l*(k-sqrt(r1(m)^2-(x(m,i)-L(m,1))^2))^3;
     end
    
     if x(m,i)>=L(m,1)+r1(m) && x(m,i)<L(m,1)+L(m,2)-r3(m)
        I(m,i)=1/12*k^3*l;
     end
   
    if x(m,i)<L(m,1)+L(m,2) && x(m,i)>=L(m,1)+L(m,2)-r3(m)
        I(m,i)=1/12*l*(k-sqrt(r3(m)^2-((L(m,1)+L(m,2)-x(m,i)))^2))^3;
    end
    if x(m,i)>=L(m,1)+L(m,2) && x(m,i)<L(m,1)+L(m,2)+r3(m)
        I(m,i)=1/12*l*(k-sqrt(r3(m)^2-((x(m,i)-L(m,1)-L(m,2)))^2))^3;
    end 
      if x(m,i)>=L(m,1)+L(m,2)+r3(m) && x(m,i)<L(m,1)+L(m,2)+L(m,3)
        I(m,i)=1/12*k^3*l;
     end
end
end
for m=2:5%inertia matrix for other fingers
    for i=1:step_size+1 
    if x(m,i)<=L(m,1)-r1(m)
        I(m,i)=1/12*k^3*l;
    end
    
    if x(m,i)<L(m,1) && x(m,i)>=L(m,1)-r1(m)
        I(m,i)=1/12*l*(k-sqrt(r1(m)^2-(L(m,1)-x(m,i))^2))^3;
    end
     if x(m,i)>=L(m,1) && x(m,i)<L(m,1)+r1(m)
        I(m,i)=1/12*l*(k-sqrt(r1(m)^2-(x(m,i)-L(m,1))^2))^3;
     end
    
     if x(m,i)>=L(m,1)+r1(m) && x(m,i)<L(m,1)+L(m,2)-r3(m)
        I(m,i)=1/12*k^3*l;
     end
   
    if x(m,i)<L(m,1)+L(m,2) && x(m,i)>=L(m,1)+L(m,2)-r3(m)
        I(m,i)=1/12*l*(k-sqrt(r3(m)^2-((L(m,1)+L(m,2)-x(m,i)))^2))^3;
    end
    if x(m,i)>=L(m,1)+L(m,2) && x(m,i)<L(m,1)+L(m,2)+r3(m)
        I(m,i)=1/12*l*(k-sqrt(r3(m)^2-((x(m,i)-L(m,1)-L(m,2)))^2))^3;
    end
    
    
     if x(m,i)>=L(m,1)+L(m,2)+r3(m) && x(m,i)<L(m,1)+L(m,2)+L(m,3)-r2(m)
        I(m,i)=1/12*k^3*l;
     end
    
    if x(m,i)>=L(m,1)+L(m,2)+L(m,3)-r2(m) && x(m,i)<L(m,1)+L(m,2)+L(m,3)
         I(m,i)=1/12*l*(k-sqrt(r2(m)^2-((L(m,1)+L(m,2)+L(m,3)-x(m,i)))^2))^3;
    end
    if x(m,i)>=L(m,1)+L(m,2)+L(m,3) && x(m,i)<L(m,1)+L(m,2)+L(m,3)+r2(m)
        I(m,i)=1/12*l*(k-sqrt(r2(m)^2-((x(m,i)-(L(m,1)+L(m,2)+L(m,3))))^2))^3;
    end
    
    if x(m,i)>=L(m,1)+L(m,2)+L(m,3)+r2(m) && x(m,i)<=L(m,1)+L(m,2)+L(m,3)+L(m,4)-r4(m)
        I(m,i)=1/12*k^3*l;
    end
     if x(m,i)>=L(m,1)+L(m,2)+L(m,3)+L(m,4)-r4(m) && x(m,i)<L(m,1)+L(m,2)+L(m,3)+L(m,4)
        I(m,i)=1/12*l*(k-sqrt(r4(m)^2-((L(m,1)+L(m,2)+L(m,3)+L(m,4)-x(m,i)))^2))^3;
     end
     if x(m,i)>=L(m,1)+L(m,2)+L(m,3)+L(m,4) && x(m,i)<=L(m,1)+L(m,2)+L(m,3)+L(m,4)+r4(m)
        I(m,i)=1/12*l*(k-sqrt(r2(m)^2-((x(m,i)-(L(m,1)+L(m,2)+L(m,3)+L(m,4))))^2))^3;
    end
end
end
for m=1:5 %inertia to deflection
I(m,step_size+1)=I(m,step_size);
for i=1:step_size+1 %shear matrix
    
    V(m,i)=-F(m);
end

M(m,1)=F(m)*TL(1,m);
for i=1:step_size
    M(m,i+1)=M(i)-(x(m,i+1)-x(m,i))*F(m);
end
for i=1:step_size+1
    C(m,i)=M(i)/(E*I(i));
end
Mc(1)=0;
for i=2:step_size 
    Mc(m,i)= Mc(i-1)+trapz(x(1:i),C(1:i));
end
for i=2:step_size
    deflection(m,i)=trapz(x(1:i),C(1:i))*180/pi;
end
end
%plot(x(1,1:1000),radtodeg(Mc(1,:)))
%%
for m=1:5 %slopeeqn
    slopeeqn(m,:)=polyfit(x(m,1:step_size+1),C(m,:),2);
end
figure
hold
legend('Show')
for z=1:5
y_1 = @(x) (slopeeqn(z,3)+slopeeqn(z,2)*x+slopeeqn(z,1)*x.^2)/3;
x_1 = 0:0.001:TL(z);
%plot(x(z,2:step_size+1),Mc(z,:))
legend('Show')
plot(x_1,y_1(x_1))
end
% 
% figure
% plot(x_1(1,:),y_1(x_1(1,:)),'DisplayName','Thumb')
% grid
% hold
% plot(x(2,2:step_size+1),Mc(2,:),'DisplayName','Index')
% plot(x(3,2:step_size+1),Mc(3,:),'DisplayName','Middle')
% plot(x(4,2:step_size+1),Mc(4,:),'DisplayName','Ring')
% plot(x(5,2:step_size+1),Mc(5,:),'DisplayName','Little')
% Create xlabel
% xlabel({'Length(m)'},'FontSize',20);
% Create ylabel
% ylabel({'Deflection(m)'},'FontSize',20);
% 
