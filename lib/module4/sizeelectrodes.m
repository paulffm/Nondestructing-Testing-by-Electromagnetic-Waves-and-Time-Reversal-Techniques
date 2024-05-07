function [XA,XB]=sizeelectrodes(b,nx,ny,nz)
np=nx*ny*nz;

 XA=zeros(3*np,1);
 XB=zeros(3*np,1);
 
 %b=0 Standard
 %b=1 beide klein
 %b=2 B klein
 %b=3 A klein
 
switch b
    case 0
        % Kanten Elektrode A
ia=3;
  for ja=3:6
        na=ia+(ja-1)*nx;
 
       
        XA(na)=1;
        XA(na+1)=1;
        XA(na+nx)=1;
        XA(na+nx+1)=1;

        XA(na+np)=1;
        XA(na+1+np)=1;
         XA(na+2+np)=1;
 

   
  end


% Kanten Elektrode B
ib=7;
  for jb=3:7
 
nb=ib+(jb-1)*nx;

    if jb<7
        XB(nb)=1;
        XB(nb+1)=1;
 
        XB(nb+np)=1;
        XB(nb+1+np)=1;
        XB(nb+2+np)=1;
        
    else 
        XB(nb)=1;
        XB(nb+1)=1;
    end 
    
  end

    case 1
        
% Kanten Elektrode A klein B klein
    ia=3;
    
    for ja=4:5 
        na=ia+(ja-1)*nx;
        %x
        XA(na)=1;
        XA(na+nx)=1;
        %y
        XA(na+np)=1;
        XA(na+1+np)=1;

    end
    
    
% Kanten Elektrode B
    ib=8;
   for jb=4:5
        nb=ib+(jb-1)*nx;
         %x
         XB(nb)=1;
         XB(nb+nx)=1;
         %y
         XB(nb+np)=1;
         XB(nb+1+np)=1;
   end
    case 2 % A groß B klein x
        
% Kanten Elektrode A
   ia=3;
  for ja=3:7
        na=ia+(ja-1)*nx;
   if ja<7
       
        XA(na)=1;
        XA(na+1)=1;

        XA(na+np)=1;
        XA(na+1+np)=1;
         XA(na+2+np)=1;
 
   else
       
        XA(na)=1;
        XA(na+1)=1;
       
   end
   
  end

% Kanten Elektrode B
    ib=8;
   for jb=4:5
        nb=ib+(jb-1)*nx;
         %x
         XB(nb)=1;
         XB(nb+nx)=1;
         %y
         XB(nb+np)=1;
         XB(nb+1+np)=1;
   
   end
   
    case 3 %A klein B groß x
        
 % Kanten Elektrode A
    ia=3;
    
    for ja=4:5 
        na=ia+(ja-1)*nx;
        %x
        XA(na)=1;
        XA(na+nx)=1;
        %y
        XA(na+np)=1;
        XA(na+1+np)=1;

    end   
% Kanten Elektrode B
ib=7;
  for jb=3:7
 
nb=ib+(jb-1)*nx;

    if jb<7
        XB(nb)=1;
        XB(nb+1)=1;
 
        XB(nb+np)=1;
        XB(nb+1+np)=1;
        XB(nb+2+np)=1;
        
    else 
        XB(nb)=1;
        XB(nb+1)=1;
    end 
    
  end

end
